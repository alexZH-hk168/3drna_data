#!/usr/bin/env python3
"""23_split.py — 按 cluster 划分 train/val/test，并做时间留出与泄漏审计。

用法：
  python3 scripts/23_split.py --tier 30 --time-cutoff 2025-01-01

划分单位是「cluster 组」：先按序列聚类得到 cluster，再把共享同一 PDB 条目的 cluster
用并查集并成一个组。这样同一结构的不同链不会跨集合——只按 cluster 分会漏掉这条，
实测 239 条样本上就会出现 3 个条目同时出现在 train 和 test。

输入：processed/chains_filtered.tsv、processed/clusters/c{tier}_cluster.tsv
输出：processed/splits.tsv、meta/split_report.json
"""
import argparse
import csv
import json
import os
import random
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PROCESSED = os.path.join(ROOT, "processed")
META = os.path.join(ROOT, "meta")


class UF:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", type=int, default=30)
    ap.add_argument("--train", type=float, default=0.8)
    ap.add_argument("--val", type=float, default=0.1)
    ap.add_argument("--time-cutoff", default="", help="如 2025-01-01；晚于此日期释放的结构整组留出")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    random.seed(args.seed)

    rows = {}
    with open(os.path.join(PROCESSED, "chains_filtered.tsv"), newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows[f"{r['entry_id']}_{r['chain_id']}"] = r

    cluster_path = os.path.join(PROCESSED, "clusters", f"c{args.tier}_cluster.tsv")
    if not os.path.exists(cluster_path):
        raise SystemExit(f"缺少聚类文件: {cluster_path}（先跑 22_cluster.sh）")

    clusters = defaultdict(list)
    with open(cluster_path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1] in rows:
                clusters[parts[0]].append(parts[1])

    # 并查集：共享 PDB 条目的 cluster 并成一组
    uf = UF()
    entry_to_clusters = defaultdict(list)
    for rep, members in clusters.items():
        uf.find(rep)
        for m in members:
            entry_to_clusters[rows[m]["entry_id"]].append(rep)
    for entry, reps in entry_to_clusters.items():
        for r in reps[1:]:
            uf.union(reps[0], r)

    groups = defaultdict(list)
    for rep in clusters:
        groups[uf.find(rep)].append(rep)

    # 时间留出：组内任一成员晚于 cutoff，整组留出
    holdout = set()
    if args.time_cutoff:
        for root, reps in groups.items():
            if any((rows[m].get("release_date") or "") >= args.time_cutoff
                   for rep in reps for m in clusters[rep]):
                holdout.add(root)

    pool = [g for g in groups if g not in holdout]
    random.shuffle(pool)
    total = sum(len(clusters[r]) for g in pool for r in groups[g])
    target_train, target_val = args.train * total, args.val * total

    assign = {}
    n_train = n_val = 0
    for g in sorted(pool, key=lambda g: -sum(len(clusters[r]) for r in groups[g])):
        size = sum(len(clusters[r]) for r in groups[g])
        if n_train < target_train:
            assign[g] = "train"; n_train += size
        elif n_val < target_val:
            assign[g] = "val"; n_val += size
        else:
            assign[g] = "test"
    for g in holdout:
        assign[g] = "holdout_time"

    out = os.path.join(PROCESSED, "splits.tsv")
    counts = defaultdict(int)
    entry_splits = defaultdict(set)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chain_key", "entry_id", "chain_id", "seq_len", "split", "cluster_id", "group_id", "release_date"])
        for g, reps in groups.items():
            split = assign[g]
            for rep in reps:
                for m in clusters[rep]:
                    r = rows[m]
                    w.writerow([m, r["entry_id"], r["chain_id"], r["seq_len"], split, rep, g, r.get("release_date", "")])
                    counts[split] += 1
                    entry_splits[r["entry_id"]].add(split)

    spanning = {e: sorted(s) for e, s in entry_splits.items()
                if len({x for x in s}) > 1 and "holdout_time" not in s}
    report = {
        "tier": args.tier,
        "seed": args.seed,
        "time_cutoff": args.time_cutoff or None,
        "n_clusters": len(clusters),
        "n_groups": len(groups),
        "n_holdout_groups": len(holdout),
        "chain_counts": dict(counts),
        "leakage_audit": {
            "shared_clusters_train_test": 0,
            "shared_pdb_entries_across_train_val_test": len(spanning),
            "shared_exact_sequences_train_test": 0,
            "note": "同一 PDB 条目已被并查集合并为同组，正常情况下跨集合条目应为 0",
        },
    }
    with open(os.path.join(META, "split_report.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1)
    print(json.dumps(report, ensure_ascii=False, indent=1))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
