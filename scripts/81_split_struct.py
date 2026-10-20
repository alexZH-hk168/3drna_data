#!/usr/bin/env python3
"""81_split_struct.py — 以「结构等价类 + PDB 条目 + 完全同序列」为分组单元的划分。

背景：序列级去冗余只保证「序列不同簇不跨集合」，防不住「序列差异、结构趋同」。
本脚本把 BGSU 4.0 Å 等价类直接并进分组约束（用户要求的第 1 步），
使同一结构类的成员必定落在同一集合。

分组单元 = 并查集合并以下关系：
  1. 共享同一 BGSU 4.0A 等价类
  2. 共享同一 PDB 条目（同一结构的不同链）
  3. 序列完全相同

划分按长度分层（<=200 / 200-1000 / >1000），比例 80/10/10。

用法：在 from_server 目录下 python3 scripts/81_split_struct.py
输出：processed/splits.tsv、meta/split_report_struct.json
"""
import csv
import json
import os
import random
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROC = os.path.join(ROOT, "processed")
META = os.path.join(ROOT, "meta")

TARGETS = {"train": 0.8, "val": 0.1, "test": 0.1}
SEED = 42


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


def load_chains():
    rows = {}
    with open(os.path.join(PROC, "chains_filtered.tsv"), newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows[f"{r['entry_id']}_{r['chain_id']}"] = r
    return rows


def load_bgsu(cutoff="4.0A"):
    m = {}
    for path in (os.path.join(META, "bgsu_classes.tsv"), os.path.join(ROOT, "meta", "bgsu_classes.tsv")):
        if not os.path.exists(path):
            continue
        with open(path, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["cutoff"] != cutoff:
                    continue
                for tok in r["members"].split(";"):
                    for part in tok.split("+"):
                        q = part.strip().split("|")
                        if len(q) >= 3:
                            m[f"{q[0]}_{q[2]}"] = r["class_id"]
        break
    return m


def main():
    random.seed(SEED)
    rows = load_chains()
    bgsu = load_bgsu("4.0A")

    uf = UF()
    for k in rows:
        uf.find(k)

    # 关系 1：同一结构等价类
    by_class = defaultdict(list)
    for k in rows:
        c = bgsu.get(k)
        if c:
            by_class[c].append(k)
    n_class_edges = 0
    for c, ks in by_class.items():
        for k in ks[1:]:
            uf.union(ks[0], k)
            n_class_edges += 1

    # 关系 2：同一 PDB 条目
    by_entry = defaultdict(list)
    for k, r in rows.items():
        by_entry[r["entry_id"]].append(k)
    n_entry_edges = 0
    for e, ks in by_entry.items():
        for k in ks[1:]:
            uf.union(ks[0], k)
            n_entry_edges += 1

    # 关系 3：完全同序列
    by_seq = defaultdict(list)
    for k, r in rows.items():
        by_seq[r["seq"]].append(k)
    for s, ks in by_seq.items():
        for k in ks[1:]:
            uf.union(ks[0], k)

    groups = defaultdict(list)
    for k in rows:
        groups[uf.find(k)].append(k)

    def bucket(k):
        L = int(rows[k]["seq_len"])
        return "long" if L > 1000 else ("mid" if L > 200 else "short")

    assign = {}
    for b in ("short", "mid", "long"):
        pool = [g for g in groups if bucket(groups[g][0]) == b]
        random.shuffle(pool)
        total_b = sum(len(groups[g]) for g in pool)
        cap = {s: TARGETS[s] * total_b for s in TARGETS}
        used = {s: 0 for s in TARGETS}
        # 大簇优先，投给剩余容量最大的集合，减少大簇导致的过冲
        for g in sorted(pool, key=lambda g: -len(groups[g])):
            s = max(TARGETS, key=lambda x: cap[x] - used[x])
            assign[g] = s
            used[s] += len(groups[g])

    counts = defaultdict(int)
    with open(os.path.join(PROC, "splits.tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chain_key", "entry_id", "chain_id", "seq_len", "split", "group_id", "release_date"])
        for g, ks in groups.items():
            for k in ks:
                r = rows[k]
                w.writerow([k, r["entry_id"], r["chain_id"], r["seq_len"], assign[g], g, r.get("release_date", "")])
                counts[assign[g]] += 1

    report = {
        "grouping": ["BGSU_4.0A_class", "PDB_entry", "identical_sequence"],
        "n_chains": len(rows),
        "n_groups": len(groups),
        "class_edges_merged": n_class_edges,
        "entry_edges_merged": n_entry_edges,
        "chains_by_split": dict(counts),
        "chains_mapped_to_bgsu": sum(1 for k in rows if k in bgsu),
    }
    with open(os.path.join(META, "split_report_struct.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1)
    print(json.dumps(report, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
