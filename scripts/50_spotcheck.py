#!/usr/bin/env python3
"""50_spotcheck.py — 独立复核：重建序列是否落在 mmCIF 声明的 RNA 序列里。

前两版失败原因：
  - 用正则读 `_entity_poly` 循环，遇到多行引号序列会截断（于是出现
    declared_lens=[50,25] 这种明显不合理的值）。
  - 用 gemmi 的 entity→subchain 映射，实测 10BE 里 label/auth 错位。

本版改用 Bio.PDB.MMCIF2Dict（独立的 CIF 词典实现，正确处理多行值），
比对方式为「重建序列是否是某条声明 RNA 序列的子序列」——因为未解析残基
会让重建序列短于声明序列，不要求逐位相等。

输出：meta/spotcheck.json
"""
import argparse
import csv
import gzip
import json
import os
import random
import tempfile
from collections import Counter

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PROCESSED = os.path.join(ROOT, "processed")
RAW = os.path.join(ROOT, "raw", "pdb")
META = os.path.join(ROOT, "meta")


def declared_rna_seqs(path):
    # MMCIF2Dict 不会自己解 gz，先解到临时文件
    if path.endswith(".gz"):
        tmp = tempfile.NamedTemporaryFile(suffix=".cif", delete=False)
        with gzip.open(path, "rb") as src:
            tmp.write(src.read())
        tmp.close()
        try:
            d = MMCIF2Dict(tmp.name)
        finally:
            os.unlink(tmp.name)
    else:
        d = MMCIF2Dict(path)
    types = d.get("_entity_poly.type", [])
    seqs = d.get("_entity_poly.pdbx_seq_one_letter_code_can", [])
    if isinstance(types, str):
        types = [types]
    if isinstance(seqs, str):
        seqs = [seqs]
    out = []
    for t, s in zip(types, seqs):
        # mmCIF 的取值是 polyribonucleotide（不是 ribonucleic）
        tl = t.lower()
        if "ribonucleotide" in tl or "ribonucleic" in tl:
            out.append(s.replace("\n", "").replace(" ", "").replace("T", "U").upper())
    return out


def is_subseq(a, b):
    it = iter(b)
    return all(ch in it for ch in a)


def matches(ours, s):
    """重建序列与声明序列是否相容。

    重建序列里修饰残基记为 N，而声明序列用的是母体碱基（如甲基化 A 写作 A），
    所以长度相等时允许 ours 的 N 匹配任意碱基；长度不等则要求 ours 是 s 的子序列
    （未解析残基会让重建序列变短）。
    """
    if len(s) == len(ours):
        return all(a == b or a == "N" for a, b in zip(ours, s))
    # 长度不等时：先去掉 ours 里的 N（修饰残基）再做子序列判定，
    # 否则长链 rRNA 的修饰位点会直接把匹配打断。
    return is_subseq(ours.replace("N", ""), s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=20)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(os.path.join(PROCESSED, "chains_filtered.tsv"),
                                    newline="", encoding="utf-8"), delimiter="\t"))
    random.Random(args.seed).shuffle(rows)
    sample = rows[: args.n]

    outcomes, details, cache = Counter(), [], {}
    for r in sample:
        eid = r["entry_id"]
        path = os.path.join(RAW, f"{eid}.cif.gz")
        if not os.path.exists(path):
            outcomes["file_missing"] += 1
            continue
        if eid not in cache:
            try:
                cache[eid] = declared_rna_seqs(path)
            except Exception as exc:
                cache[eid] = None
                details.append({"chain": f"{eid}_{r['chain_id']}", "status": f"cif_error:{type(exc).__name__}"})
                outcomes["cif_error"] += 1
                continue
        seqs = cache[eid]
        ours = r["seq"].replace("T", "U").upper()
        if seqs is None:
            continue
        if not seqs:
            outcomes["no_rna_entity"] += 1
            details.append({"chain": f"{eid}_{r['chain_id']}", "status": "no_rna_entity"})
            continue
        if ours in seqs:
            outcomes["exact_match"] += 1
            details.append({"chain": f"{eid}_{r['chain_id']}", "status": "ok", "len": len(ours)})
        elif any(matches(ours, s) for s in seqs):
            outcomes["subsequence"] += 1
            details.append({"chain": f"{eid}_{r['chain_id']}", "status": "ok_subseq", "len": len(ours)})
        else:
            outcomes["mismatch"] += 1
            details.append({"chain": f"{eid}_{r['chain_id']}", "status": "mismatch",
                            "ours_len": len(ours),
                            "declared_lens": [len(s) for s in seqs],
                            "longest_overlap_hint": max((len(os.path.commonprefix([ours, s])) for s in seqs), default=0)})

    passed = outcomes["exact_match"] + outcomes["subsequence"]
    report = {"sampled": len(sample), "seed": args.seed, "outcomes": dict(outcomes),
              "pass": passed == len(sample), "details": details}
    with open(os.path.join(META, "spotcheck.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1)
    print(json.dumps({k: v for k, v in report.items() if k != "details"}, ensure_ascii=False, indent=1))
    for d in details:
        if d.get("status") not in ("ok", "ok_subseq"):
            print("  !", d)


if __name__ == "__main__":
    main()
