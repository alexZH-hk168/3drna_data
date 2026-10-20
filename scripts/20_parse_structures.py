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


def rna_entity_sequences(st):
    """chain 名 → 实体声明序列（只收 RNA 实体）。

    gemmi 的 `Entity.subchains` 是 chain 名字符串列表；
    `Entity.full_sequence` 是**列表**，每个元素是一个残基的单字母/三字母代号（如 UMS），
    不是字符串。非规范残基统一记为 N，以便与规范序列比例统计对齐。
    """
    out = {}
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
            name = sub if isinstance(sub, str) else (
                sub[0] if isinstance(sub, (tuple, list)) else getattr(sub, "name", None))
            if name:
                out[name] = seq
    return out


def parse_one(path, entry_meta):
    st = open_structure(path)
    entry_id = os.path.basename(path).split(".")[0]

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

        chains.append({
            "entry_id": entry_id,
            "chain_id": chain.name,
            "seq": seq,
            "seq_len": len(seq),
            "declared_len": declared_len,
            "n_resolved_residues": len(residues),
            "n_backbone_complete": complete,
            "backbone_completeness": round(complete / declared_len, 4),
            "residue_completeness": round(len(residues) / declared_len, 4),
            "n_modified": len(mods),
            "modified_types": sorted(set(mods)),
            "canonical_fraction": round(
                sum(1 for c in seq if c in "ACGU") / max(len(seq), 1), 4),
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


def main():
    os.makedirs(INTERIM, exist_ok=True)
    meta = {}
    with open(os.path.join(META, "rcsb_entries.jsonl"), encoding="utf-8") as fh:
        for line in fh:
            r = json.loads(line)
            meta[r["entry_id"]] = r

    files = sorted(glob.glob(os.path.join(RAW, "*.cif.gz")))
    print(f"待解析 {len(files)} 个文件", file=sys.stderr)

    out_path = os.path.join(INTERIM, "chains.jsonl")
    n_chains = n_err = 0
    with open(out_path, "w", encoding="utf-8") as out:
        for i, path in enumerate(files, 1):
            entry_id = os.path.basename(path).split(".")[0]
            try:
                chains = parse_one(path, meta.get(entry_id, {}))
            except Exception as exc:
                n_err += 1
                print(f"  ERR {entry_id}: {type(exc).__name__}: {exc}", file=sys.stderr)
                continue
            for c in chains:
                out.write(json.dumps(c, ensure_ascii=False) + "\n")
                n_chains += 1
            if i % 100 == 0 or i == len(files):
                print(f"  {i}/{len(files)}  chains={n_chains} err={n_err}", file=sys.stderr)

    print(json.dumps({"files": len(files), "chains": n_chains, "parse_errors": n_err,
                      "out": os.path.relpath(out_path, ROOT)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
