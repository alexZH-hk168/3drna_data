#!/usr/bin/env python3
"""24_bgsu_crosscheck.py — 把过滤后的链对到 BGSU 等价类，验证 C4 口径。

BGSU 成员写法有两种：单链 `PDB|model|chain`（如 8D28|1|B），
多链复合体用 `+` 连成一个单元（如 11AA|1|1+11AA|1|2）。
本流水线键：`entry_chain`（如 8D28_B）。
输出：meta/bgsu_crosscheck.json
"""
import csv
import json
import os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PROCESSED = os.path.join(ROOT, "processed")
META = os.path.join(ROOT, "meta")


def bgsu_keys(token):
    """把 BGSU 成员 token 拆成单个 entry_chain 键；多链用 + 连接的逐个拆开。"""
    out = []
    for part in token.split("+"):
        p = part.strip().split("|")
        if len(p) >= 3:
            out.append(f"{p[0]}_{p[2]}")
    return out


def main():
    classes = defaultdict(dict)  # cutoff -> key -> class_id
    with open(os.path.join(META, "bgsu_classes.tsv"), newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cutoff = row["cutoff"]
            for tok in row["members"].split(";"):
                for k in bgsu_keys(tok):
                    classes[cutoff][k] = row["class_id"]

    ours = []
    with open(os.path.join(PROCESSED, "chains_filtered.tsv"), newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ours.append(f"{row['entry_id']}_{row['chain_id']}")

    report = {"chains_checked": len(ours), "by_cutoff": {}}
    for cutoff in ["1.5A", "2.0A", "2.5A", "3.0A", "3.5A", "4.0A", "20.0A", "all"]:
        m = classes.get(cutoff, {})
        hit = [c for c in ours if c in m]
        report["by_cutoff"][cutoff] = {
            "matched": len(hit),
            "coverage": round(len(hit) / len(ours), 3) if ours else 0,
            "distinct_classes": len({m[c] for c in hit}),
        }

    matched_all = [c for c in ours if c in classes.get("all", {})]
    unmatched = [c for c in ours if c not in classes.get("all", {})]
    report["unmatched_examples"] = unmatched[:10]
    report["note"] = (
        "对不上 BGSU 的链可能原因：该结构尚未进入 BGSU 4.61、链级过滤口径不同、"
        "或 BGSU 只收 resolution ≤ 其档位要求的条目"
    )
    with open(os.path.join(META, "bgsu_crosscheck.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1)
    print(json.dumps(report, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
