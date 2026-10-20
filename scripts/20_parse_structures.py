#!/usr/bin/env python3
"""20_parse_structures.py — 从 mmCIF 提取链级记录。

关键设计（2026-10-20 两次修正）：
  1. 不能依赖 gemmi 的 entity→subchain 映射：实测 10BE 里 entity 声明 A 是 RNA，
     而 chain A 实际是 792 个残基的蛋白，说明 label/auth 链名在这种结构里是错位的。
     因此按化学性质判链：有 P 与 C4'/C3' 视为核酸，全链存在 O2' 视为 RNA。
  2. 完整度用残基的 label_seq 跨度作分母（= 声明编号长度），不用实体序列长度：
     completeness = 已解析核酸残基数 / (max label_seq - min label_seq + 1)。
  序列本身仍由残基重建，已知局限：边界上 1–2 个未解析残基无法体现。

输入：raw/pdb/*.cif.gz + meta/rcsb_entries.jsonl
输出：interim/chains.jsonl
"""
import glob
import gzip
import json
import os
import sys
import tempfile

import gemmi

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
META = os.path.join(ROOT, "meta")
RAW = os.path.join(ROOT, "raw", "pdb")
INTERIM = os.path.join(ROOT, "interim")

STD_MAP = {"A": "A", "C": "C", "G": "G", "U": "U",
           "ADE": "A", "CYT": "C", "GUA": "G", "URI": "U"}
BACKBONE = ["P", "O5'", "C5'", "C4'", "C3'", "O3'"]


def open_structure(path):
    try:
        return gemmi.read_structure(path)
    except Exception:
        tmp = tempfile.NamedTemporaryFile(suffix=".cif", delete=False)
        with gzip.open(path, "rb") as src:
            tmp.write(src.read())
        tmp.close()
        try:
            return gemmi.read_structure(tmp.name)
        finally:
            os.unlink(tmp.name)


def is_nucleotide(res):
    names = {a.name.upper() for a in res}
    return "P" in names and ("C4'" in names or "C3'" in names)


def _cif_block(path):
    """读源文件的 CIF block（不能用 st.make_mmcif_document：它会重生成，丢掉 label 字段）。"""
    try:
        return gemmi.cif.read(path).sole_block()
    except Exception:
        tmp = tempfile.NamedTemporaryFile(suffix=".cif", delete=False)
        with gzip.open(path, "rb") as src:
            tmp.write(src.read())
        tmp.close()
        try:
            return gemmi.cif.read(tmp.name).sole_block()
        finally:
            os.unlink(tmp.name)


def label_to_auth(st, path):
    """label_asym_id → auth_asym_id（= gemmi 的 chain.name）。

    ⚠ 实测坑：10BE 里 label_asym 是 A/B/C，auth_asym 是 B/C/A；
    `entity.subchains` 登记的是 **label** 名，而 `chain.name` 是 **auth** 名。
    两者不直接相等，必须建映射，否则序列会被贴到错误的链上。
    """
    m = {}
    try:
        t = _cif_block(path).find_mmcif_category("_atom_site.")
        # gemmi 的 Table 是可索引对象，按 "._atom_site.label_asym_id" 全名取列
        for row in t:
            la, aa = row["_atom_site.label_asym_id"], row["_atom_site.auth_asym_id"]
            if la:
                m.setdefault(la, aa or la)
    except Exception:
        pass
    return m


def rna_entity_sequences(st, path):
    """chain 名 → 实体声明序列（只收 RNA 实体）。

    ⚠ 实测坑：gemmi 在 auth 链名重复/缺省时会自己给链加后缀，
    导致 `chain.name` 与 `entity.subchains` 的名字错位。
    例：10BE 文件里实际链是 B/C/A，gemmi 读成 A/B/C，于是 A 链拿到了 B 链的序列。
    因此这里同时用 label_asym 与 auth_asym 两个名字登记，由调用方按相似度二选一。
    """
    out = {}
    l2a = label_to_auth(st, path)
    for ent in st.entities:
        if ent.entity_type != gemmi.EntityType.Polymer:
            continue
        if ent.polymer_type != gemmi.PolymerType.Rna:
            continue
        chars = []
        for tok in ent.full_sequence:
            t = str(tok).strip().upper()
            chars.append(t if len(t) == 1 and t in "ACGU" else "N")
        seq = "".join(chars)
        for sub in ent.subchains:
            names = []
            if isinstance(sub, str):
                names = [sub]
            elif isinstance(sub, (tuple, list)):
                names = [x for x in sub if isinstance(x, str)]
            for name in names:
                if name:
                    out.setdefault(l2a.get(name, name), seq)
                    out.setdefault(name, seq)
    return out


def parse_one(path, entry_meta):
    st = open_structure(path)
    entry_id = os.path.basename(path).split(".")[0]
    declared_by_chain = rna_entity_sequences(st, path)

    # label_seq→符号 的全局映射：直接读源 mmCIF 的 label_seq_id + label_comp_id
    label_sym = {}
    try:
        t = _cif_block(path).find_mmcif_category("_atom_site.")
        for row in t:
            ls = row["_atom_site.label_seq_id"].strip()
            ci = row["_atom_site.label_comp_id"].strip().upper()
            if ls.isdigit():
                label_sym.setdefault(int(ls), STD_MAP.get(ci, "N"))
    except Exception:
        label_sym = {}

    n_protein = n_dna = 0
    chains = []
    for chain in st[0]:
        residues = [r for r in chain if is_nucleotide(r)]
        if not residues:
            if any("CA" in {a.name for a in r} for r in chain):
                n_protein += 1
            continue

        # RNA / DNA 判定：全链任一核酸残基带 O2' 即为 RNA（只看前 20 个残基会误判）
        if not any("O2'" in {a.name.upper() for a in r} for r in residues):
            n_dna += 1
            continue

        chars, mods, complete = [], [], 0
        for r in residues:
            names = {a.name.upper() for a in r}
            ch = STD_MAP.get(r.name.upper(), "N")
            if ch == "N":
                mods.append(r.name.upper())
            chars.append(ch)
            if all(b in names for b in BACKBONE):
                complete += 1

        seq = "".join(chars)
        labels = [r.label_seq for r in residues if getattr(r, "label_seq", 0) and r.label_seq > 0]
        declared_len = (max(labels) - min(labels) + 1) if labels else len(residues)
        declared_len = max(declared_len, len(residues))

        # 序列口径一：由已解析残基重建（保留，作对照）
        # 序列口径二：按 label_seq 跨度铺满，缺口填 N（未解析残基 → N）
        # 序列口径三：有实体声明就用实体声明（最权威），否则退回口径二
        seq_spanned = ""
        if labels:
            lo, hi = min(labels), max(labels)
            seen = {r.label_seq: STD_MAP.get(r.name.upper(), "N") for r in residues
                    if getattr(r, "label_seq", 0) and r.label_seq > 0}
            seq_spanned = "".join(seen.get(i) or label_sym.get(i, "N")
                                  for i in range(lo, hi + 1))
        declared = declared_by_chain.get(chain.name)
        if declared and seq_spanned:
            from difflib import SequenceMatcher
            same = SequenceMatcher(None, seq_spanned, declared).ratio()
            # 映射错位时（相似度低）宁可不用声明序列，也不要把别人的序列贴上来
            if same < 0.9:
                declared = None
        else:
            same = None
        seq_final = declared if declared else (seq_spanned or seq)

        chains.append({
            "entry_id": entry_id,
            "chain_id": chain.name,
            "seq": seq_final,
            "seq_rebuilt": seq,
            "seq_source": "entity" if declared else ("label_seq" if seq_spanned else "rebuilt"),
            "span_vs_entity_similarity": None if same is None else round(same, 4),
            "seq_len": len(seq_final),
            "declared_len": declared_len,
            "n_resolved_residues": len(residues),
            "n_backbone_complete": complete,
            "backbone_completeness": round(complete / declared_len, 4),
            "residue_completeness": round(len(residues) / declared_len, 4),
            "n_modified": len(mods),
            "modified_types": sorted(set(mods)),
            "canonical_fraction": round(
                sum(1 for c in seq_final if c in "ACGU") / max(len(seq_final), 1), 4),
            "unknown_fraction": round(
                sum(1 for c in seq_final if c == "N") / max(len(seq_final), 1), 4),
            "n_rna_chains": None,
            "n_protein_chains": None,
            "n_dna_chains": None,
            "method": entry_meta.get("method"),
            "resolution": entry_meta.get("resolution"),
            "release_date": entry_meta.get("release_date"),
        })

    n_rna = len(chains)
    for c in chains:
        c["n_rna_chains"] = n_rna
        c["n_protein_chains"] = n_protein
        c["n_dna_chains"] = n_dna
    return chains


def _work(job):
    path, entry_meta = job
    entry_id = os.path.basename(path).split(".")[0]
    try:
        return entry_id, parse_one(path, entry_meta), None
    except Exception as exc:
        return entry_id, [], f"{type(exc).__name__}: {exc}"


def main():
    import argparse
    from concurrent.futures import ProcessPoolExecutor

    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    os.makedirs(INTERIM, exist_ok=True)
    meta = {}
    with open(os.path.join(META, "rcsb_entries.jsonl"), encoding="utf-8") as fh:
        for line in fh:
            r = json.loads(line)
            meta[r["entry_id"]] = r

    files = sorted(glob.glob(os.path.join(RAW, "*.cif.gz")))
    print(f"待解析 {len(files)} 个文件，{args.workers} 进程", file=sys.stderr)

    out_path = os.path.join(INTERIM, "chains.jsonl")
    n_chains = n_err = 0
    jobs = [(p, meta.get(os.path.basename(p).split(".")[0], {})) for p in files]
    with open(out_path, "w", encoding="utf-8") as out, \
            ProcessPoolExecutor(max_workers=args.workers) as pool:
        for i, (entry_id, chains, err) in enumerate(pool.map(_work, jobs, chunksize=20), 1):
            if err:
                n_err += 1
                print(f"  ERR {entry_id}: {err}", file=sys.stderr)
            for c in chains:
                out.write(json.dumps(c, ensure_ascii=False) + "\n")
                n_chains += 1
            if i % 500 == 0 or i == len(files):
                print(f"  {i}/{len(files)}  chains={n_chains} err={n_err}", file=sys.stderr)

    print(json.dumps({"files": len(files), "chains": n_chains, "parse_errors": n_err,
                      "out": os.path.relpath(out_path, ROOT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
