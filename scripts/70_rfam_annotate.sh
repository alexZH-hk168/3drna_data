#!/usr/bin/env bash
# 70_rfam_annotate.sh — 用 Infernal cmscan 给 30% 档 cluster 代表元打 Rfam 家族
#
# 为什么只注释代表元：每个簇内的成员序列高度相似，同一个簇共用一个家族标签即可；
# 4,180 条全量跑比 1,657 条代表元慢约 2.5 倍，收益很小。
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
export PATH="$ROOT/tools/mmseqs/bin:$PATH"

CM_DIR="$ROOT/tools/rfam"
CLUSTER_TSV="$ROOT/processed/clusters/cluster_30_cluster.tsv"
OUT_TSV="$ROOT/meta/rfam_hits.tsv"
REP_FA="$ROOT/interim/reps30.fasta"

mkdir -p "$CM_DIR" logs

echo "== 1) 解压并索引 Rfam.cm"
if [ ! -s "$CM_DIR/Rfam.cm" ]; then gunzip -c "$CM_DIR/Rfam.cm.gz" > "$CM_DIR/Rfam.cm"; fi
if [ ! -s "$CM_DIR/Rfam.cm.i1f" ]; then cmpress "$CM_DIR/Rfam.cm"; fi
ls -lh "$CM_DIR"/Rfam.cm*

echo "== 2) 抽代表元序列"
python3 - "$CLUSTER_TSV" "$ROOT/processed/chains.fasta" "$REP_FA" <<'PY'
import sys
ctsv, fasta, out = sys.argv[1:4]
reps = {line.split('\t')[0] for line in open(ctsv, encoding='utf-8') if line.strip()}
seqs, name, buf = {}, None, []
for line in open(fasta, encoding='utf-8'):
    if line.startswith('>'):
        if name: seqs[name] = ''.join(buf)
        name, buf = line[1:].strip(), []
    else:
        buf.append(line.strip())
if name: seqs[name] = ''.join(buf)
n = 0
with open(out, 'w', encoding='utf-8') as fh:
    for r in reps:
        if r in seqs:
            fh.write(f">{r}\n{seqs[r]}\n"); n += 1
print(f"代表元 {len(reps)}，写入序列 {n}", file=sys.stderr)
PY

echo "== 3) cmscan（RNA 序列 → Rfam 模型）"
cmscan --cpu 8 --noali --tblout "$OUT_TSV" "$CM_DIR/Rfam.cm" "$REP_FA" > logs/cmscan.log 2>&1 || true
echo "命中行数: $(grep -vc '^#' "$OUT_TSV" 2>/dev/null || echo 0)"
head -3 "$OUT_TSV" 2>/dev/null || true
