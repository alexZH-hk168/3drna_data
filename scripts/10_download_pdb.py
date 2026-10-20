#!/usr/bin/env python3
"""10_download_pdb.py — 按索引清单下载含 RNA 的 PDB 条目（mmCIF.gz）。

用法：
  python3 scripts/10_download_pdb.py --limit 50          # 小样试跑
  python3 scripts/10_download_pdb.py                     # 全量
  python3 scripts/10_download_pdb.py --max-resolution 4.0 --methods X-RAY,ELECTRON

特性：断点续传（已存在且非空则跳过）、并发、逐条 sha256 记账。
"""
import argparse
import gzip
import hashlib
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
META = os.path.join(ROOT, "meta")
RAW = os.path.join(ROOT, "raw", "pdb")
URL = "https://files.rcsb.org/download/{}.cif.gz"


def load_index(max_resolution=None, methods=None):
    rows = []
    with open(os.path.join(META, "rcsb_entries.jsonl"), encoding="utf-8") as fh:
        for line in fh:
            r = json.loads(line)
            if max_resolution is not None:
                res = r.get("resolution")
                if res is not None and res > max_resolution:
                    continue
            if methods:
                m = " ".join(r.get("method") or []).upper()
                if not any(k in m for k in methods):
                    continue
            rows.append(r)
    return rows


def fetch(entry_id, retries=3):
    path = os.path.join(RAW, f"{entry_id}.cif.gz")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        with open(path, "rb") as fh:
            data = fh.read()
        return entry_id, path, data, "cached"
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(URL.format(entry_id), headers={"User-Agent": "rna3d-prep/0.1"})
            with urllib.request.urlopen(req, timeout=90) as resp:
                data = resp.read()
            if not data:
                raise ValueError("empty body")
            # 落盘前校验 gzip 头
            if data[:2] != b"\x1f\x8b":
                raise ValueError("not gzip")
            with open(path, "wb") as fh:
                fh.write(data)
            return entry_id, path, data, "downloaded"
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError,
                ValueError, http.client.HTTPException, OSError) as exc:
            last = exc
            time.sleep(1.5 * (attempt + 1))
        except Exception as exc:  # 不吞异常，但别让单条杀死整个任务
            last = exc
            time.sleep(1.0)
    return entry_id, path, None, f"failed:{last}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--max-resolution", type=float, default=None)
    ap.add_argument("--methods", default=None, help="逗号分隔，如 X-RAY,ELECTRON")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--shuffle", action="store_true",
                    help="小样试跑前先按固定种子打乱，避免只取到最新那批（索引默认按发布日期倒序）")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    os.makedirs(RAW, exist_ok=True)
    methods = [m.strip().upper() for m in args.methods.split(",")] if args.methods else None
    rows = load_index(args.max_resolution, methods)
    if args.shuffle:
        import random
        random.Random(args.seed).shuffle(rows)
    if args.limit:
        rows = rows[: args.limit]
    print(f"待下载 {len(rows)} 条", file=sys.stderr)

    manifest_path = os.path.join(META, "download_manifest.tsv")
    new_file = not os.path.exists(manifest_path)
    fh = open(manifest_path, "a", encoding="utf-8")
    if new_file:
        fh.write("entry_id\tstatus\tbytes\tsha256\tpath\n")
        fh.flush()

    ok = fail = 0
    started = time.time()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(fetch, r["entry_id"]): r["entry_id"] for r in rows}
        for i, fut in enumerate(as_completed(futures), 1):
            entry_id, path, data, status = fut.result()
            if data is None:
                fail += 1
                fh.write(f"{entry_id}\t{status}\t0\t-\t-\n")
            else:
                ok += 1
                digest = hashlib.sha256(data).hexdigest()
                rel = os.path.relpath(path, ROOT)
                fh.write(f"{entry_id}\t{status}\t{len(data)}\t{digest}\t{rel}\n")
            if i % 50 == 0 or i == len(rows):
                fh.flush()
                rate = i / max(time.time() - started, 1e-6)
                print(f"  {i}/{len(rows)}  ok={ok} fail={fail}  {rate:.1f}/s", file=sys.stderr)
    fh.close()
    print(json.dumps({"requested": len(rows), "ok": ok, "fail": fail,
                      "seconds": round(time.time() - started, 1)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
