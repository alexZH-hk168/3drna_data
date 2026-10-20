#!/usr/bin/env python3
"""71_sort_reps.py — 按序列长度升序重排簇代表元 FASTA。

cmscan 在全长 rRNA（>1500 nt）上会卡住很久；把短序列排在前面，
即使整轮没跑完也能先覆盖大部分语料。

用法：python3 scripts/71_sort_reps.py [in.fasta] [out.fasta]
"""
import sys

src = sys.argv[1] if len(sys.argv) > 1 else "processed/clusters/c30_rep_seq.fasta"
dst = sys.argv[2] if len(sys.argv) > 2 else "processed/clusters/c30_rep_sorted.fasta"

recs, name, buf = [], None, []
for line in open(src):
    line = line.rstrip("\n")
    if line.startswith(">"):
        if name:
            recs.append((name, "".join(buf)))
        name, buf = line, []
    else:
        buf.append(line.strip())
if name:
    recs.append((name, "".join(buf)))

recs.sort(key=lambda x: len(x[1]))
with open(dst, "w") as fh:
    for n, s in recs:
        fh.write(n + "\n" + s + "\n")

lens = [len(s) for _, s in recs]
print(f"{len(recs)} reps -> {dst}")
print(f"shortest={lens[0]} median={lens[len(lens)//2]} longest={lens[-1]}")
print(f"<=500nt: {sum(1 for x in lens if x <= 500)} | >1500nt: {sum(1 for x in lens if x > 1500)}")
