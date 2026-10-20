#!/usr/bin/env python3
"""60_download_pdbredo.py — 按过滤后的条目清单定向下载 PDB-REDO 精修坐标。

入口：https://pdb-redo.eu/db/{id}/{id}_final.cif（实测 1ehz 返回 200，323 KB）
整库 tar.bz2 太大，按需取子集即可。

输入：processed/chains_filtered.tsv
输出：raw/pdbredo/{ID}.cif、meta/pdbredo_manifest.tsv、meta/pdbredo_summary.json
"""
import argparse
import csv
import http.client
import json
import os
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PROCESSED = os.path.join(ROOT, "processed")
META = os.path.join(ROOT, "meta")
OUT = os.path.join(ROOT, "raw", "pdbredo")
URL = "https://pdb-redo.eu/db/{lo}/{lo}_final.cif"


def fetch(entry_id):
    lo = entry_id.lower()
    path = os.path.join(OUT, f"{entry_id}.cif")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return entry_id, path, "cached"
    last = None
    for attempt in range(3):
        try:
            req = urllib.request.Request(URL.format(lo=lo), headers={"User-Agent": "rna3d-prep/0.1"})
            with urllib.request.urlopen(req, timeout=180) as resp:
                data = resp.read()
            if not data:
                raise ValueError("empty")
            with open(path, "wb") as fh:
                fh.write(data)
            return entry_id, path, "downloaded"
        except urllib.error.HTTPError as exc:
            # 实测：不在 PDB-REDO 里的条目返回 500（不是 404），直接记为不可用，不重试
            if exc.code in (404, 500):
                return entry_id, None, "not_in_pdbredo"
            last = exc
            time.sleep(1.5 * (attempt + 1))
        except (urllib.error.URLError, TimeoutError, http.client.HTTPException, OSError, ValueError) as exc:
            last = exc
            time.sleep(1.5 * (attempt + 1))
    return entry_id, None, f"failed:{type(last).__name__}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    entries = sorted({r["entry_id"] for r in csv.DictReader(
        open(os.path.join(PROCESSED, "chains_filtered.tsv"), newline="", encoding="utf-8"), delimiter="\t")})
    print(f"待取 {len(entries)} 个条目", file=sys.stderr)

    manifest = os.path.join(META, "pdbredo_manifest.tsv")
    fh = open(manifest, "w", encoding="utf-8")
    fh.write("entry_id\tstatus\tbytes\n")

    counts, total = {}, 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(fetch, e): e for e in entries}
        for i, fut in enumerate(as_completed(futures), 1):
            entry_id, path, status = fut.result()
            key = status.split(":")[0]
            counts[key] = counts.get(key, 0) + 1
            size = os.path.getsize(path) if path else 0
            total += size
            fh.write(f"{entry_id}\t{status}\t{size}\n")
            if i % 200 == 0 or i == len(entries):
                fh.flush()
                print(f"  {i}/{len(entries)} {counts}", file=sys.stderr)
    fh.close()

    summary = {"entries": len(entries), "outcomes": counts,
               "bytes": total, "gib": round(total / 1024 ** 3, 3)}
    json.dump(summary, open(os.path.join(META, "pdbredo_summary.json"), "w"), ensure_ascii=False, indent=1)
    print(json.dumps(summary, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
