#!/usr/bin/env python3
"""22a_clstr_to_tsv.py — 把 CD-HIT 的 .clstr 转成两列 cluster.tsv（代表元 / 成员）。

用法：python3 scripts/22a_clstr_to_tsv.py in.clstr out_prefix
产物：<out_prefix>_cluster.tsv
"""
import re
import sys

TOKEN = re.compile(r">([^\s.]+)")


def main():
    src, prefix = sys.argv[1], sys.argv[2]
    rep = None
    members = []
    rows = []
    with open(src, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">Cluster"):
                if rep:
                    rows.append((rep, rep))
                    rows.extend((rep, m) for m in members)
                rep, members = None, []
                continue
            m = TOKEN.search(line)
            if not m:
                continue
            name = m.group(1)
            if line.rstrip().endswith("*"):
                rep = name
            else:
                members.append(name)
    if rep:
        rows.append((rep, rep))
        rows.extend((rep, m) for m in members)

    with open(f"{prefix}_cluster.tsv", "w", encoding="utf-8") as fh:
        for r, m in rows:
            fh.write(f"{r}\t{m}\n")
    n_clusters = len({r for r, _ in rows})
    print(f"{prefix}: {len(rows)} sequences -> {n_clusters} clusters")


if __name__ == "__main__":
    main()
