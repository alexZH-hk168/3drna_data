#!/usr/bin/env python3
"""80_struct_leakage.py — 用 BGSU 结构等价类检查划分是否存在「结构跨集合」。

序列级去冗余只保证「序列不同簇不跨集合」。RNA 的结构保守性强，
序列相似度低于 30% 的两条链仍可能折叠成同一拓扑。这里用 BGSU 的
结构等价类（官方按几何相似性归并）作为结构分组的代理，量化这一层泄漏。

输入：from_server/processed/{chains_filtered.tsv, splits.tsv}
      from_server/meta/bgsu_classes.tsv
输出：meta/struct_leakage.json（屏幕打印摘要）
"""
import csv
import json
import os
from collections import Counter, defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "from_server")
META = os.path.join(ROOT, "meta")


def load_splits():
    out = {}
    with open(os.path.join(SRC, "processed", "splits.tsv"), newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["chain_key"]] = r["split"]
    return out


def load_bgsu(cutoff):
    m = {}
    p = os.path.join(SRC, "meta", "bgsu_classes.tsv")
    with open(p, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["cutoff"] != cutoff:
                continue
            for tok in r["members"].split(";"):
                for part in tok.split("+"):
                    q = part.strip().split("|")
                    if len(q) >= 3:
                        m[f"{q[0]}_{q[2]}"] = r["class_id"]
    return m


def main():
    splits = load_splits()
    report = {"note": "结构分组代理 = BGSU 等价类；holdout_time 不参与冲突判定",
              "chains_in_splits": len(splits),
              "by_cutoff": {}}

    for cutoff in ("4.0A", "all"):
        bgsu = load_bgsu(cutoff)
        cls2splits = defaultdict(set)
        cls2chains = defaultdict(list)
        unmapped = 0
        for key, split in splits.items():
            c = bgsu.get(key)
            if c is None:
                unmapped += 1
                continue
            cls2splits[c].add(split)
            cls2chains[c].append(key)

        core = {"train", "val", "test"}
        cross = {c: s for c, s in cls2splits.items()
                 if len(s & core) > 1 and "holdout_time" not in s}
        involved = sum(len(cls2chains[c]) for c in cross)
        pairs = Counter()
        for c, s in cross.items():
            pairs["+".join(sorted(s & core))] += 1

        # 每个 test 链是否与某条 train 链同属一个结构类
        train_cls = {c for c, s in cls2splits.items() if "train" in s}
        test_keys = [k for k, s in splits.items() if s == "test"]
        test_hit = sum(1 for k in test_keys if bgsu.get(k) in train_cls)

        report["by_cutoff"][cutoff] = {
            "chains_mapped_to_bgsu": len(splits) - unmapped,
            "unmapped": unmapped,
            "n_bgsu_classes_used": len(cls2splits),
            "classes_spanning_splits": len(cross),
            "chains_in_spanning_classes": involved,
            "spanning_patterns": dict(pairs.most_common()),
            "test_chains_sharing_structure_with_train": test_hit,
            "test_chains_total": len(test_keys),
            "test_structural_contamination_rate": round(test_hit / max(len(test_keys), 1), 3),
        }

    os.makedirs(META, exist_ok=True)
    with open(os.path.join(META, "struct_leakage.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, ensure_ascii=False, indent=1)
    print(json.dumps(report, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
