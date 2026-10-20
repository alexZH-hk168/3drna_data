#!/usr/bin/env bash
# setup_env.sh — 在数据机器上准备目录与工具链。
# 用法：bash scripts/setup_env.sh [DATA_ROOT]
# 默认 DATA_ROOT 为本仓库目录（脚本所在目录的上一级）。
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${1:-$(dirname "$HERE")}"
cd "$ROOT"

echo "==> 项目根目录: $ROOT"
mkdir -p raw/{pdb,rna3db,rnasolo,bgsu,rnaglib,pdbredo,isolated} \
         interim processed meta scripts logs tools
echo "==> 目录已建"

# ---------- 1) 静态二进制：MMseqs2 ----------
if ! command -v mmseqs >/dev/null 2>&1; then
  echo "==> 安装 MMseqs2 (static)"
  ( cd tools && \
    curl -fsSL -o mmseqs.tar.gz https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz && \
    tar xzf mmseqs.tar.gz && rm -f mmseqs.tar.gz )
  MMSEQS_BIN="$(find "$ROOT/tools" -maxdepth 3 -type f -name mmseqs | head -1 || true)"
  if [ -n "${MMSEQS_BIN}" ]; then
    ln -sf "$MMSEQS_BIN" "$ROOT/tools/mmseqs-link"
    export PATH="$(dirname "$MMSEQS_BIN"):$PATH"
    echo "    mmseqs -> $MMSEQS_BIN"
  else
    echo "    WARN: mmseqs 解包后未找到可执行文件" >&2
  fi
else
  echo "==> MMseqs2 已存在: $(command -v mmseqs)"
fi

# ---------- 2) CD-HIT（对照用） ----------
if ! command -v cd-hit-est >/dev/null 2>&1; then
  echo "==> 安装 CD-HIT"
  ( cd tools && \
    curl -fsSL -o cdhit.tar.gz https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz && \
    tar xzf cdhit.tar.gz && rm -f cdhit.tar.gz && \
    cd cd-hit-v4.8.1-2019-0228 && make -j"$(nproc)" >/dev/null 2>&1 ) || \
    echo "    WARN: CD-HIT 编译失败，聚类仍可用 MMseqs2" >&2
else
  echo "==> CD-HIT 已存在: $(command -v cd-hit-est)"
fi

# ---------- 3) Python 依赖 ----------
echo "==> 安装 Python 依赖 (gemmi / biopython / pandas / pyarrow)"
python3 -m pip install --quiet --upgrade pip >/dev/null 2>&1 || true
python3 -m pip install --quiet gemmi biopython pandas pyarrow || \
  echo "    WARN: pip 安装失败，请手动安装 gemmi biopython pandas pyarrow" >&2

# ---------- 4) Infernal（家族注释，可选） ----------
if ! command -v cmscan >/dev/null 2>&1; then
  echo "==> 安装 Infernal (conda 优先)"
  if command -v conda >/dev/null 2>&1; then
    conda install -y -q -c bioconda -c conda-forge infernal >/dev/null 2>&1 && \
      echo "    infernal 已装（conda）" || echo "    WARN: conda 安装 infernal 失败" >&2
  else
    echo "    SKIP: 无 conda，Infernal 稍后手动安装（可复用 RNA3DB 自带 cmscan 结果）"
  fi
else
  echo "==> Infernal 已存在: $(command -v cmscan)"
fi

# ---------- 5) 版本锁 ----------
LOCK="$ROOT/meta/env_versions.txt"
{
  echo "# 生成时间: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "# host: $(hostname)"
  echo "# nproc: $(nproc)"
  echo "python3: $(python3 -V 2>&1)"
  command -v mmseqs >/dev/null 2>&1 && echo "mmseqs: $(mmseqs version 2>/dev/null || echo unknown)" || echo "mmseqs: MISSING"
  command -v cd-hit-est >/dev/null 2>&1 && echo "cd-hit-est: $(cd-hit-est -h 2>&1 | head -1)" || echo "cd-hit-est: MISSING"
  command -v cmscan >/dev/null 2>&1 && echo "cmscan: $(cmscan -h 2>&1 | head -1)" || echo "cmscan: MISSING"
  python3 - <<'PY' 2>/dev/null || true
import importlib
for m in ("gemmi", "Bio", "pandas", "pyarrow", "numpy"):
    try:
        mod = importlib.import_module(m)
        print(f"{m}: {getattr(mod, '__version__', 'unknown')}")
    except Exception:
        print(f"{m}: MISSING")
PY
} > "$LOCK"
echo "==> 版本锁写入 $LOCK"
cat "$LOCK"
