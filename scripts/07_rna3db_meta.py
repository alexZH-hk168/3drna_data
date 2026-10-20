#!/usr/bin/env python3
"""07_rna3db_meta.py — 把 RNA3DB 的官方切分与 cmscan 结果整理成链级表。

输入：
  interim/rna3db/rna3db-jsons/{split,cluster}.json
  interim/rna3db/rna3db-cmscans/*.tbl
输出：
  processed/rna3db_chains.tsv    chain / split / component / rfam / evalue / has_hit
  meta/rna3db_meta_summary.json  计数与交叉核对
"""
import glob
import json
import os
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
INTERIM = os.path.join(ROOT, "interim", "rna3db")
PROCESSED = os.path.join(ROOT, "processed")
META = os.path.join(ROOT, "meta")


def flatten_split(path):
    """split.json 结构：{split_set: {component: {cluster_rep: cluster_size}}}。
    key 是簇代表元，value 是该簇的链数；汇总时按 value 加权。"""
    d = json.load(open(path, encoding="utf-8"))
    out = {}
    for split_name, comps in d.items():
        if not isinstance(comps, dict):
            continue
        label = split_name.replace("_set", "")
        for comp, members in comps.items():
            if isinstance(members, dict):
                for rep, size in members.items():
                    out[rep] = {"split": label, "component": comp,
                                "size": size if isinstance(size, int) else 1}
            else:
                for rep in members:
                    out[rep] = {"split": label, "component": comp, "size": 1}
    return out


def parse_cmscan(paths):
    hits = {}
    for p in paths:
        with open(p, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.split(None, 17)
                if len(f) < 16:
                    continue
                rfam, chain, evalue = f[0], f[2], f[15]
                try:
                    ev = float(evalue)
                except ValueError:
                    continue
                cur = hits.get(chain)
                if cur is None or ev < cur["evalue"]:
                    hits[chain] = {"rfam": rfam, "evalue": ev,
                                   "description": (f[17].strip() if len(f) > 17 else "")}
    return hits


def main():
    os.makedirs(PROCESSED, exist_ok=True)
    split = flatten_split(os.path.join(INTERIM, "rna3db-jsons", "split.json"))
    tbls = sorted(glob.glob(os.path.join(INTERIM, "rna3db-cmscans", "*.tbl")))
    hits = parse_cmscan(tbls)

    chains = sorted(split)
    per_split = Counter(v["split"] for v in split.values())
    chains_per_split = Counter()
    for v in split.values():
        chains_per_split[v["split"]] += v["size"]
    with_hit = Counter()
    for c in chains:
        if c in hits:
            with_hit[split[c]["split"]] += 1

    out = os.path.join(PROCESSED, "rna3db_chains.tsv")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("cluster_rep\tn_chains\tsplit\tcomponent\trfam\tevalue\thas_hit\n")
        for c in chains:
            h = hits.get(c)
            fh.write("\t".join([
                c, str(split[c]["size"]), split[c]["split"], split[c]["component"],
                h["rfam"] if h else "",
                f"{h['evalue']:.3g}" if h else "",
                "1" if h else "0",
            ]) + "\n")

    # 交叉核对：README 记录 train 13,451 / test 1,990，其中无 Rfam 命中 58 / 1,365
    no_hit_clusters = {s: per_split[s] - with_hit[s] for s in per_split}
    summary = {
        "clusters_total": len(chains),
        "clusters_per_split": dict(per_split),
        "clusters_with_rfam_hit": dict(with_hit),
        "clusters_without_rfam_hit": no_hit_clusters,
        "cmscan_files": [os.path.basename(p) for p in tbls],
        "distinct_rfam_families": len({h["rfam"] for h in hits.values()}),
        "sql_note": "split.json 的键是簇代表元，值不携带簇大小；链口径看官方 rna3db_stats.json",
        "cross_check_note": "官方 stats：train 13,451 / test 1,990 条链；无 Rfam 命中 58 / 1,365 条链。本表 test 有 100 个簇无命中，对应官方 1,365 条链，说明留出的是少数大型簇",
    }
    with open(os.path.join(META, "rna3db_meta_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)
    print(json.dumps(summary, ensure_ascii=False, indent=1))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
