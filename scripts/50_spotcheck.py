#!/usr/bin/env python3
"""50_spotcheck.py — 独立复核序列是否与 mmCIF 声明一致。

做法：不走 gemmi，也不走 Biopython，直接对 mmCIF 文本解析 `_entity_poly` 循环，
取 `pdbx_seq_one_letter_code_can`，与 processed/chains_filtered.tsv 里的序列比较。

这一步防的是「解析器自己把序列拼错了还以为对了」。

输出：meta/spotcheck.json
"""
import argparse
import csv
import gzip
import json
import os
import random
import re
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PROCESSED = os.path.join(ROOT, "processed")
RAW = os.path.join(ROOT, "raw", "pdb")
META = os.path.join(ROOT, "meta")


def entity_poly_seqs(path):
    """从 mmCIF 文本里抽 _entity_poly 的 (entity_id, type, one_letter_seq)。"""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    m = re.search(r"loop_\s*\n(.*?)(?=\n#)", text, re.S)
    # 只处理包含 _entity_poly. 的 loop
    for block in re.finditer(r"loop_\s*\n((?:\s*_entity_poly\.[^\n]*\n)+)(.*?)(?=\n#)", text, re.S):
        tags = re.findall(r"_entity_poly\.(\S+)", block.group(1))
        body = block.group(2)
        if "pdbx_seq_one_letter_code_can" not in tags or "entity_id" not in tags:
            continue
        i_id = tags.index("entity_id")
        i_seq = tags.index("pdbx_seq_one_letter_code_can")
        i_type = tags.index("type") if "type" in tags else None
        rows = []
        for line in body.split("\n"):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) <= max(i_id, i_seq):
                continue
            rows.append((parts[i_id], parts[i_type] if i_type is not None else "", parts[i_seq]))
        return rows
    return []


def is_subseq(a, b):
    """a 是否是 b 的子序列（未解析残基会让重建序列短于声明序列）。"""
    it = iter(b)
    return all(ch in it for ch in a)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=20)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(os.path.join(PROCESSED, "chains_filtered.tsv"), newline="", encoding="utf-8"), delimiter="\t"))
    random.Random(args.seed).shuffle(rows)
    sample = rows[: args.n]

    outcomes, details = Counter(), []
    cache = {}
    for r in sample:
        eid = r["entry_id"]
        path = os.path.join(RAW, f"{eid}.cif.gz")
        if not os.path.exists(path):
            outcomes["file_missing"] += 1
            continue
        if eid not in cache:
            cache[eid] = entity_poly_seqs(path)
        seqs = [s.replace("T", "U") for (_iid, _ty, s) in cache[eid]]
        ours = r["seq"].replace("T", "U")
        if not seqs:
            outcomes["no_entity_poly"] += 1
            details.append({"chain": f"{eid}_{r['chain_id']}", "status": "no_entity_poly"})
            continue
        if any(is_subseq(ours, s) for s in seqs):
            outcomes["subsequence_of_declared"] += 1
            details.append({"chain": f"{eid}_{r['chain_id']}", "status": "ok", "len": len(ours)})
        elif ours in seqs:
            outcomes["exact_match"] += 1
            details.append({"chain": f"{eid}_{r['chain_id']}", "status": "ok", "len": len(ours)})
        else:
            outcomes["mismatch"] += 1
            details.append({"chain": f"{eid}_{r['chain_id']}", "status": "mismatch",
                            "ours_len": len(ours), "declared_lens": [len(s) for s in seqs]})

    report = {"sampled": len(sample), "seed": args.seed, "outcomes": dict(outcomes),
              "pass": (outcomes["exact_match"] + outcomes["subsequence_of_declared"]) == len(sample),
              "details": details}
    with open(os.path.join(META, "spotcheck.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1)
    print(json.dumps({k: v for k, v in report.items() if k != "details"}, ensure_ascii=False, indent=1))
    for d in details:
        if d.get("status") != "ok":
            print("  !", d)


if __name__ == "__main__":
    main()
