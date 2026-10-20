#!/usr/bin/env python3
"""00_rcsb_index.py — 拉取 RCSB 中含 RNA 链的全部条目索引。

产出：
  meta/rcsb_entries.jsonl        每条 entry 一行（id、方法、分辨率、链组成、日期）
  meta/rcsb_index_summary.json   计数汇总与抓取时间

只读网络请求，不下载结构文件。
"""
import json
import os
import sys
import time
import urllib.error
import urllib.request
from collections import Counter

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"
PAGE = 10000
BATCH = 500

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
META = os.path.join(ROOT, "meta")


def _post(url, payload, timeout=120, retries=3):
    body = json.dumps(payload).encode()
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, data=body, headers={"Content-Type": "application/json"}
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read())
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
            last = exc
            if isinstance(exc, urllib.error.HTTPError) and exc.code == 400:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError(f"request failed after {retries} attempts: {last}")


def list_rna_entries():
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.polymer_entity_count_RNA",
                "operator": "greater",
                "value": 0,
            },
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": PAGE},
            "results_verbosity": "compact",
        },
    }
    ids = []
    start = 0
    total = None
    while True:
        query["request_options"]["paginate"]["start"] = start
        r = _post(SEARCH_URL, query)
        total = r.get("total_count", total)
        batch = r.get("result_set", []) or []
        if not batch:
            break
        ids.extend(batch)
        start += len(batch)
        if len(batch) < PAGE:
            break
    return ids, total


GQL = """
query($ids:[String!]!) {
  entries(entry_ids:$ids) {
    rcsb_id
    struct { title }
    exptl { method }
    rcsb_entry_info {
      resolution_combined
      polymer_entity_count_RNA
      polymer_entity_count_protein
      polymer_entity_count_DNA
      polymer_entity_count_nucleic_acid
      structure_determination_methodology
    }
    rcsb_accession_info { initial_release_date deposit_date }
  }
}
"""


def fetch_meta(ids):
    out = []
    for i in range(0, len(ids), BATCH):
        chunk = ids[i : i + BATCH]
        r = _post(GRAPHQL_URL, {"query": GQL, "variables": {"ids": chunk}})
        if "errors" in r:
            raise RuntimeError(json.dumps(r["errors"])[:800])
        for e in r["data"]["entries"] or []:
            info = e.get("rcsb_entry_info") or {}
            acc = e.get("rcsb_accession_info") or {}
            res = info.get("resolution_combined") or [None]
            out.append(
                {
                    "entry_id": e["rcsb_id"],
                    "title": (e.get("struct") or {}).get("title"),
                    "method": [m.get("method") for m in (e.get("exptl") or [])],
                    "resolution": min([x for x in res if x is not None], default=None),
                    "methodology": info.get("structure_determination_methodology"),
                    "n_rna": info.get("polymer_entity_count_RNA"),
                    "n_protein": info.get("polymer_entity_count_protein"),
                    "n_dna": info.get("polymer_entity_count_DNA"),
                    "n_nucleic": info.get("polymer_entity_count_nucleic_acid"),
                    "release_date": acc.get("initial_release_date"),
                    "deposit_date": acc.get("deposit_date"),
                }
            )
        sys.stderr.write(f"  fetched {min(i+BATCH, len(ids))}/{len(ids)}\n")
    return out


def main():
    os.makedirs(META, exist_ok=True)
    ids, total = list_rna_entries()
    sys.stderr.write(f"search total_count={total}, collected={len(ids)}\n")
    if total is not None and len(ids) != total:
        sys.stderr.write("WARN: collected count != total_count\n")

    rows = fetch_meta(ids)
    path = os.path.join(META, "rcsb_entries.jsonl")
    with open(path, "w", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row, ensure_ascii=False) + "\n")

    def bucket(method_list):
        m = (method_list or ["?"])[0] or "?"
        for key in ("X-RAY", "ELECTRON MICROSCOPY", "NMR", "HYBRID", "NEUTRON"):
            if key in m.upper():
                return key
        return "OTHER"

    summary = {
        "fetched_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "search_total_count": total,
        "collected": len(ids),
        "written": len(rows),
        "by_methodology": dict(
            Counter(bucket(r["method"]) for r in rows).most_common()
        ),
        "rna_only_entries": sum(
            1 for r in rows if r["n_rna"] and not r["n_protein"] and not r["n_dna"]
        ),
        "with_resolution": sum(1 for r in rows if r["resolution"] is not None),
        "release_year_min": min(
            (r["release_date"][:4] for r in rows if r.get("release_date")), default=None
        ),
        "release_year_max": max(
            (r["release_date"][:4] for r in rows if r.get("release_date")), default=None
        ),
    }
    with open(os.path.join(META, "rcsb_index_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, ensure_ascii=False, indent=1)
    print(json.dumps(summary, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
