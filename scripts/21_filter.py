#!/usr/bin/env python3
"""21_filter.py — F1–F6 过滤漏斗。

输入：interim/chains.jsonl
输出：processed/chains_filtered.tsv、meta/funnel.json、meta/rejection_reasons.json

阈值默认取数据合同 v0.1 的建议值，可用命令行覆盖。
"""
import argparse
import json
import os
from collections import Counter, OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
INTERIM = os.path.join(ROOT, "interim")
PROCESSED = os.path.join(ROOT, "processed")
META = os.path.join(ROOT, "meta")

FIELDS = ["entry_id", "chain_id", "seq_len", "backbone_completeness", "canonical_fraction",
          "n_modified", "method", "resolution", "release_date",
          "n_rna_chains", "n_protein_chains", "n_dna_chains", "seq"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res-max", type=float, default=4.0, help="X-ray/cryo-EM 分辨率上限")
    ap.add_argument("--res-ext", type=float, default=6.0, help="扩展池上限")
    ap.add_argument("--len-min", type=int, default=20)
    ap.add_argument("--len-max", type=int, default=10000,
                    help="链长上限；导师要求覆盖 >1000 nt 的长链 rRNA，默认放到 10000")
    ap.add_argument("--backbone-min", type=float, default=0.90)
    ap.add_argument("--canonical-min", type=float, default=0.95)
    args = ap.parse_args()

    os.makedirs(PROCESSED, exist_ok=True)
    rows = []
    with open(os.path.join(INTERIM, "chains.jsonl"), encoding="utf-8") as fh:
        for line in fh:
            rows.append(json.loads(line))

    funnel = OrderedDict()
    reasons = Counter()
    funnel["F0_all_chains"] = len(rows)

    # F1 含 RNA 链（解析阶段已只保留 RNA 链，这里保留占位以便对齐口径）
    stage = rows
    funnel["F1_has_rna_chain"] = len(stage)

    # F2 方法 + 分辨率
    keep, ext = [], []
    for r in stage:
        method = " ".join(r.get("method") or []).upper()
        res = r.get("resolution")
        if "NMR" in method:
            keep.append(r)          # NMR 无分辨率，单独标记
            r["_pool"] = "nmr"
            continue
        if res is None:
            reasons["F2_no_resolution"] += 1
            continue
        if res <= args.res_max:
            r["_pool"] = "main"
            keep.append(r)
        elif res <= args.res_ext:
            r["_pool"] = "extended"
            ext.append(r)
        else:
            reasons["F2_resolution_too_low"] += 1
    funnel["F2_method_resolution"] = len(keep)
    funnel["F2_extended_pool"] = len(ext)

    # F3 链长
    s3 = []
    for r in keep:
        if r["seq_len"] < args.len_min:
            reasons["F3_too_short"] += 1
        elif r["seq_len"] > args.len_max:
            reasons["F3_too_long"] += 1
        else:
            s3.append(r)
    funnel["F3_length"] = len(s3)

    # F4 主链完整性
    s4 = []
    for r in s3:
        if r["backbone_completeness"] < args.backbone_min:
            reasons["F4_incomplete_backbone"] += 1
        else:
            s4.append(r)
    funnel["F4_backbone_completeness"] = len(s4)

    # F5 规范序列比例
    s5 = []
    for r in s4:
        if r["canonical_fraction"] < args.canonical_min:
            reasons["F5_too_many_modified"] += 1
        else:
            s5.append(r)
    funnel["F5_canonical_sequence"] = len(s5)

    # F6 完全重复合并（同一序列只留一条，记录被合并的条目）
    seen = {}
    s6 = []
    dup = 0
    for r in sorted(s5, key=lambda x: (-x["seq_len"], x["entry_id"], x["chain_id"])):
        key = r["seq"]
        if key in seen:
            dup += 1
            seen[key]["merged_from"] = seen[key].get("merged_from", 0) + 1
            continue
        r["merged_from"] = 0
        seen[key] = r
        s6.append(r)
    funnel["F6_exact_dedup"] = len(s6)
    funnel["F6_merged_duplicates"] = dup

    out = os.path.join(PROCESSED, "chains_filtered.tsv")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(FIELDS + ["merged_from", "pool"]) + "\n")
        for r in s6:
            fh.write("\t".join(str(r.get(k, "")) for k in FIELDS + ["merged_from", "_pool"]) + "\n")

    with open(os.path.join(META, "funnel.json"), "w", encoding="utf-8") as fh:
        json.dump({"thresholds": vars(args), "funnel": funnel,
                   "rejection_reasons": dict(reasons)}, fh, ensure_ascii=False, indent=1)

    print(json.dumps({"thresholds": vars(args), "funnel": funnel,
                      "rejection_reasons": dict(reasons)}, ensure_ascii=False, indent=1))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
