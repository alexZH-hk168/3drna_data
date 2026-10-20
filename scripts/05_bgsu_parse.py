#!/usr/bin/env python3
"""05_bgsu_parse.py — 把 BGSU nrlist 各分辨率档的 CSV 整理成标准等价类表。

输入：raw/bgsu/nrlist_4.59_{cutoff}.csv（无表头，列序：class_id, representative, members）
输出：
  meta/bgsu_classes.tsv          每条等价类一行
  meta/bgsu_summary.json         各档的类数与覆盖链数

这些等价类是后面结构级聚类（C4）的外部对照标尺。
"""
import csv
import glob
import json
import os
import re
import sys
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RAW = os.path.join(ROOT, "raw", "bgsu")
META = os.path.join(ROOT, "meta")

ORDER = ["1.5A", "2.0A", "2.5A", "3.0A", "3.5A", "4.0A", "20.0A", "all"]
RELEASE = sys.argv[1] if len(sys.argv) > 1 else "4.61"


def main():
    os.makedirs(META, exist_ok=True)
    rows = []
    summary = OrderedDict()
    for cutoff in ORDER:
        path = os.path.join(RAW, f"nrlist_{RELEASE}_{cutoff}.csv")
        if not os.path.exists(path):
            summary[cutoff] = {"status": "missing"}
            continue
        n_classes = 0
        n_members = 0
        sizes = []
        with open(path, newline="", encoding="utf-8", errors="replace") as fh:
            for rec in csv.reader(fh):
                if len(rec) < 3:
                    continue
                class_id = rec[0].strip()
                members = [m.strip() for m in rec[2].split(",") if m.strip()]
                if not class_id or not members:
                    continue
                n_classes += 1
                n_members += len(members)
                sizes.append(len(members))
                rows.append(
                    {
                        "cutoff": cutoff,
                        "class_id": class_id,
                        "representative": rec[1].strip(),
                        "n_members": len(members),
                        "members": ";".join(members),
                    }
                )
        summary[cutoff] = {
            "n_classes": n_classes,
            "n_chains": n_members,
            "largest_class": max(sizes) if sizes else 0,
            "singleton_classes": sum(1 for s in sizes if s == 1),
            "mean_class_size": round(n_members / n_classes, 2) if n_classes else 0,
        }

    out = os.path.join(META, "bgsu_classes.tsv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["cutoff", "class_id", "representative", "n_members", "members"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    with open(os.path.join(META, "bgsu_summary.json"), "w", encoding="utf-8") as fh:
        json.dump({"release": RELEASE, "by_cutoff": summary}, fh, ensure_ascii=False, indent=1)

    print(json.dumps({"release": RELEASE, "by_cutoff": summary}, ensure_ascii=False, indent=1))
    print(f"\nwrote {out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
