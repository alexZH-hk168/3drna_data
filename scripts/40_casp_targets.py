#!/usr/bin/env python3
"""40_casp_targets.py — 从 CASP 靶标清单里抽出核酸相关靶标与其 PDB 编号。

输入：raw/isolated/casp/casp{15,16}_targetlist.html
输出：
  meta/casp_targets.tsv         target / type / length / pdb_id / description
  meta/casp_holdout_pdb.json    用于训练隔离的 PDB 条目清单

CASP 靶标页每行形如：
  276.  R1296 NucA  72  A1  2024-07-29 ... - G34U mutant of M1209 complex PDB code 9c75
所以可以用「目标号 + 类型 + 长度」定位，再向后找 PDB code。
"""
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RAW = os.path.join(ROOT, "raw", "isolated", "casp")
META = os.path.join(ROOT, "meta")

ENTRY = re.compile(r"(\d+)\.\s+([A-Z]\d{3,4}\w*)\s+([A-Za-z][A-Za-z\-/ ]{0,20}?)\s+(\d{1,6})\s+([A-Z0-9a-z]+)")
PDB = re.compile(r"PDB code\s+([0-9][A-Za-z0-9]{3})", re.I)
# CASP16 用 NucA / Prot-NucA；CASP15 用 RNA / RNA /Ligand，两套写法都要收进去
NUA = ("NucA", "Prot-NucA", "NucA/Ligand", "RNA", "Prot-RNA", "DNA", "Prot-DNA")


def clean(html):
    t = re.sub(r"<script.*?</script>", " ", html, flags=re.S | re.I)
    t = re.sub(r"<style.*?</style>", " ", t, flags=re.S | re.I)
    t = re.sub(r"<[^>]+>", " ", t)
    return re.sub(r"[ \t]+", " ", t)


def parse(path, casp):
    txt = clean(open(path, encoding="utf-8", errors="replace").read())
    ms = list(ENTRY.finditer(txt))
    rows = []
    for i, m in enumerate(ms):
        end = ms[i + 1].start() if i + 1 < len(ms) else len(txt)
        tail = txt[m.end():end]
        pm = PDB.search(tail)
        desc = re.sub(r"\s+", " ", tail[:160]).strip()
        rows.append({
            "casp": casp,
            "target": m.group(2),
            "type": m.group(3).strip(),
            "length": int(m.group(4)),
            "stoichiometry": m.group(5),
            "pdb_id": pm.group(1).upper() if pm else "",
            "description": desc,
        })
    return rows


def main():
    all_rows = []
    for casp in ("15", "16"):
        p = os.path.join(RAW, f"casp{casp}_targetlist.html")
        if not os.path.exists(p):
            print(f"missing {p}")
            continue
        rows = parse(p, casp)
        print(f"CASP{casp}: {len(rows)} targets parsed")
        all_rows.extend(rows)

    seen, uniq = set(), []
    for r in all_rows:
        k = (r["casp"], r["target"])
        if k in seen:
            continue
        seen.add(k)
        uniq.append(r)

    with open(os.path.join(META, "casp_targets.tsv"), "w", encoding="utf-8") as fh:
        fh.write("casp\ttarget\ttype\tlength\tstoichiometry\tpdb_id\tdescription\n")
        for r in uniq:
            fh.write("\t".join(str(r[k]) for k in
                               ("casp", "target", "type", "length", "stoichiometry", "pdb_id", "description")) + "\n")

    nuc = [r for r in uniq if any(r["type"].startswith(x) for x in NUA)]
    pdb_ids = sorted({r["pdb_id"] for r in nuc if r["pdb_id"]})
    summary = {
        "targets_parsed": len(uniq),
        "nucleic_acid_targets": len(nuc),
        "with_pdb_id": len(pdb_ids),
        "by_type": {t: sum(1 for r in nuc if r["type"].startswith(t)) for t in NUA},
        "by_casp": {c: sum(1 for r in nuc if r["casp"] == c) for c in ("15", "16")},
        "holdout_pdb_ids": pdb_ids,
        "note": "用于训练隔离：这些 PDB 条目及其同 cluster 成员不得进入训练集",
    }
    with open(os.path.join(META, "casp_holdout_pdb.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)
    print(json.dumps({k: v for k, v in summary.items() if k != "holdout_pdb_ids"}, ensure_ascii=False, indent=1))
    print("pdb ids:", len(pdb_ids))


if __name__ == "__main__":
    main()
