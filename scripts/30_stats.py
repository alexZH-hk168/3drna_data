#!/usr/bin/env python3
"""30_stats.py — 过滤后数据集的分布统计。

输入：processed/chains_filtered.tsv、meta/funnel.json、meta/cluster_summary.json
输出：meta/distribution_stats.json
"""
import csv
import json
import os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PROCESSED = os.path.join(ROOT, "processed")
META = os.path.join(ROOT, "meta")

LEN_BINS = [(20, 50), (50, 100), (100, 200), (200, 400), (400, 700), (700, 1001)]
RES_BINS = [(0, 2.0), (2.0, 2.5), (2.5, 3.0), (3.0, 3.5), (3.5, 4.01)]


def label(v, bins, fmt="{:g}-{:g}"):
    for lo, hi in bins:
        if lo <= v < hi:
            return fmt.format(lo, hi)
    return "other"


def main():
    rows = list(csv.DictReader(open(os.path.join(PROCESSED, "chains_filtered.tsv"), newline="", encoding="utf-8"), delimiter="\t"))
    n = len(rows)

    lens = [int(r["seq_len"]) for r in rows]
    res = [float(r["resolution"]) for r in rows if r["resolution"] not in ("", "None")]

    def method_bucket(r):
        m = (r["method"] or "").upper()
        if "ELECTRON" in m:
            return "cryo-EM"
        if "X-RAY" in m:
            return "X-ray"
        if "NMR" in m:
            return "NMR"
        return "other"

    chains = Counter()
    for r in rows:
        if int(r["n_protein_chains"] or 0) == 0 and int(r["n_dna_chains"] or 0) == 0:
            chains["RNA-only"] += 1
        elif int(r["n_protein_chains"] or 0) > 0:
            chains["with protein"] += 1
        else:
            chains["with DNA"] += 1

    stats = {
        "n_chains": n,
        "n_entries": len({r["entry_id"] for r in rows}),
        "length": {
            "min": min(lens), "max": max(lens),
            "median": sorted(lens)[n // 2],
            "mean": round(sum(lens) / n, 1),
            "histogram": dict(sorted(Counter(label(v, LEN_BINS) for v in lens).items())),
        },
        "method": dict(Counter(method_bucket(r) for r in rows).most_common()),
        "resolution": {
            "n_with_resolution": len(res),
            "median": sorted(res)[len(res) // 2] if res else None,
            "histogram": dict(sorted(Counter(label(v, RES_BINS, "{:.1f}-{:.1f}") for v in res).items())),
        },
        "chain_composition": dict(chains),
        "merged_duplicates": sum(int(r["merged_from"]) for r in rows),
        "pool": dict(Counter(r["pool"] for r in rows)),
    }

    for name in ("funnel.json", "cluster_summary.json"):
        p = os.path.join(META, name)
        if os.path.exists(p):
            stats[name.replace(".json", "")] = json.load(open(p, encoding="utf-8"))

    with open(os.path.join(META, "distribution_stats.json"), "w", encoding="utf-8") as fh:
        json.dump(stats, fh, ensure_ascii=False, indent=1)
    print(json.dumps({k: v for k, v in stats.items() if k != "funnel"}, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
