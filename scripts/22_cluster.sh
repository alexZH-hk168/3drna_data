#!/usr/bin/env bash
# 22_cluster.sh — 多阈值序列聚类（CD-HIT-EST 为主，MMseqs2 备选）
#
# 用法：
#   bash scripts/22_cluster.sh              # 默认 CD-HIT-EST
#   bash scripts/22_cluster.sh mmseqs       # 改用 MMseqs2
#
# 说明：部分容器里 MMseqs2 会误判可用内存并拒绝运行
#       （报 "Cannot fit databases into …"），CD-HIT-EST 不依赖该判断，因此设为主路径。
#
# 产物：processed/clusters/c{tier}_cluster.tsv（两列：代表元 / 成员），供 22b、23 使用
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
cd "$ROOT"

TOOL="${1:-cdhit}"
TSV="processed/chains_filtered.tsv"
FASTA="processed/chains.fasta"
OUTDIR="processed/clusters"
TMP="$ROOT/interim/cluster_tmp"
mkdir -p "$OUTDIR" "$TMP" "$ROOT/logs" "$ROOT/tools"

# 1) 生成 FASTA
python3 - <<'PY'
import csv
src, dst = "processed/chains_filtered.tsv", "processed/chains.fasta"
n = 0
with open(src, newline="", encoding="utf-8") as fh, open(dst, "w", encoding="utf-8") as out:
    for row in csv.DictReader(fh, delimiter="\t"):
        out.write(f">{row['entry_id']}_{row['chain_id']}\n{row['seq']}\n")
        n += 1
print(f"FASTA: {n} sequences -> {dst}")
PY

TIERS="100 90 70 50 30"

case "$TOOL" in
  cdhit)
    CDHIT="$(command -v cd-hit-est || find "$ROOT/tools" -maxdepth 4 -name cd-hit-est -type f 2>/dev/null | head -1)"
    if [ -z "${CDHIT}" ]; then echo "ERROR: 找不到 cd-hit-est，先跑 scripts/setup_env.sh" >&2; exit 127; fi
    echo "==> 使用 CD-HIT-EST: ${CDHIT}"
    # 实测：V4.8.1 的 cd-hit-est 强制要求 -c >= 0.8，低于 0.8 直接 Fatal Error。
    # 因此 CD-HIT 只负责 100/90 两档；70/50/30 必须用 MMseqs2。
    for ID in 100 90; do
      case "$ID" in
        100) N=10 ;; 90) N=8 ;;
      esac
      PREFIX="$OUTDIR/c${ID}"
      [ -s "${PREFIX}_cluster.tsv" ] && { echo "==> c${ID} 已存在，跳过"; continue; }
      echo "==> 聚类 ${ID}% (c=$([ "$ID" = 100 ] && echo 1.0 || echo "0.${ID}") n=${N})"
      timeout 3600 "$CDHIT" -i "$FASTA" -o "$TMP/c${ID}" \
        -c "$([ "$ID" = 100 ] && echo 1.0 || echo "0.${ID}")" -n "$N" \
        -M 0 -T "$(nproc)" -d 0 -l 10 > "logs/cluster_c${ID}.log" 2>&1
      RC=$?
      if [ $RC -ne 0 ]; then echo "    WARN: c${ID} 退出码 ${RC}，见 logs/cluster_c${ID}.log" >&2; continue; fi
      python3 scripts/22a_clstr_to_tsv.py "$TMP/c${ID}.clstr" "$PREFIX"
    done
    echo "==> 注意：70/50/30 档 CD-HIT-EST 不支持（-c 下限 0.8），须用 mmseqs 跑"
    ;;
  mmseqs)
    if ! command -v mmseqs >/dev/null 2>&1; then echo "ERROR: mmseqs 不在 PATH" >&2; exit 127; fi
    echo "==> 使用 MMseqs2"
    for ID in $TIERS; do
      MINID=$(python3 -c "print(${ID}/100)")
      PREFIX="$OUTDIR/c${ID}"
      [ -s "${PREFIX}_cluster.tsv" ] && { echo "==> c${ID} 已存在，跳过"; continue; }
      mkdir -p "$TMP/c${ID}"
      echo "==> 聚类 ${ID}%"
      timeout 3600 mmseqs easy-cluster "$FASTA" "$PREFIX" "$TMP/c${ID}" \
        --min-seq-id "$MINID" -c 0.8 --cov-mode 0 --threads "$(nproc)" -v 1 \
        > "logs/cluster_c${ID}.log" 2>&1
      RC=$?
      [ $RC -ne 0 ] && echo "    WARN: c${ID} 退出码 ${RC}" >&2
    done
    ;;
  *) echo "ERROR: 未知工具 $TOOL（可选 cdhit / mmseqs）" >&2; exit 2 ;;
esac

echo "==> 汇总"
python3 scripts/22b_cluster_summary.py
