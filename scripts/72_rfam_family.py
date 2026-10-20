#!/usr/bin/env python3
"""72_rfam_family.py — 从 cmscan 的 tblout 生成簇代表元的家族归属表。

输出：meta/rfam_family.tsv（query, rfam_acc, rfam_name, evalue, coverage, desc）
      每个 query 只保留 E-value 最小的那条命中。
"""
import sys, csv
from collections import defaultdict

tbl = sys.argv[1] if len(sys.argv) > 1 else "meta/rfam_cmscan.tsv"
fasta = sys.argv[2] if len(sys.argv) > 2 else "processed/clusters/c30_rep_sorted.fasta"

best = {}
for line in open(tbl, encoding="utf-8", errors="replace"):
    if line.startswith("#") or not line.strip():
        continue
    f = line.split()
    if len(f) < 18:
        continue
    rfam, acc, query = f[0], f[1], f[2]
    try:
        ev = float(f[15])
    except ValueError:
        continue
    try:
        cov = (int(f[8]) - int(f[7]) + 1)
    except ValueError:
        cov = 0
    if query not in best or ev < best[query][2]:
        best[query] = (rfam, acc, ev, cov, " ".join(f[17:]))

total = 0
name = None
for line in open(fasta):
    if line.startswith(">"):
        total += 1

with open("meta/rfam_family.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["query", "rfam_acc", "rfam_name", "evalue", "hit_len", "description"])
    for q, (nm, acc, ev, cov, desc) in sorted(best.items()):
        w.writerow([q, acc, nm, f"{ev:.2e}", cov, desc[:120]])

fams = defaultdict(int)
for q, v in best.items():
    fams[v[0]] += 1
print(f"代表元总数 {total}")
print(f"有 Rfam 命中 {len(best)} ({len(best)/total*100:.1f}%)")
print(f"命中不同家族 {len(fams)}")
print("Top10 家族:")
for nm, n in sorted(fams.items(), key=lambda x: -x[1])[:10]:
    print(f"  {nm:24s} {n}")
