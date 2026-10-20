#!/usr/bin/env python3
"""22b_cluster_summary.py — 汇总五档聚类结果。

输入：processed/clusters/c{100,90,70,50,30}_cluster.tsv（MMseqs easy-cluster 输出，两列：代表元 / 成员）
输出：meta/cluster_summary.json、processed/clusters/cluster_members.tsv
"""
import glob
import json
import os
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CLU = os.path.join(ROOT, "processed", "clusters")
META = os.path.join(ROOT, "meta")

TIERS = [100, 90, 70, 50, 30]


def main():
    summary = {}
    rows = []
    for tier in TIERS:
        path = os.path.join(CLU, f"c{tier}_cluster.tsv")
        if not os.path.exists(path):
            summary[tier] = {"status": "missing"}
            continue
        clusters = defaultdict(list)
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                rep, member = parts[0], parts[1]
                clusters[rep].append(member)
        sizes = [len(v) for v in clusters.values()]
        summary[tier] = {
            "n_sequences": sum(sizes),
            "n_clusters": len(clusters),
            "singleton_clusters": sum(1 for s in sizes if s == 1),
            "largest_cluster": max(sizes) if sizes else 0,
            "mean_cluster_size": round(sum(sizes) / len(sizes), 2) if sizes else 0,
        }
        for rep, members in clusters.items():
            for m in members:
                rows.append((tier, rep, m))

    with open(os.path.join(CLU, "cluster_members.tsv"), "w", encoding="utf-8") as fh:
        fh.write("tier\trepresentative\tmember\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    with open(os.path.join(META, "cluster_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)
    print(json.dumps(summary, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
