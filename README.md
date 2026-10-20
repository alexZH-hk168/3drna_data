# RNA 3D 数据准备流水线

把台账里 12 条三维数据源接入，做成可复跑的「下载 → 解析 → 过滤 → 聚类 → 划分 → 泄漏审计」链路，产出可用于「序列 → 3D 坐标」训练的语料。

数据合同（阈值与字段定义）见 [数据合同.md](./数据合同.md)，逐源状态见 [meta/source_registry.json](./meta/source_registry.json)。

## 最终规模（2026-10-20）

| 项目 | 数值 |
|---|---|
| RCSB 含 RNA 链的条目 | 10,502 |
| 下载体积 | 15.0 GiB |
| 解析出的 RNA 链 | 28,949 |
| 过滤后可用语料 | **4,180 条链 / 3,072 个结构** |
| 其中 >1000 nt 长链 | 1,251 条（26.5%） |
| 30% 档序列簇 | 1,661（最大 418） |
| 划分分组 | 1,282 组 |
| train / val / test | 3,344 / 418 / 418（80 / 10 / 10） |
| 结构泄漏（4.0 Å 结构类跨集合） | 0 |

## 目录

```
raw/        原样下载（pdb/ rna3db/ pdbredo/ bgsu/ isolated/）
interim/    解析中间件（chains.jsonl、聚类临时库）
processed/  最终数据集（chains_filtered.tsv、splits.tsv、clusters/）
meta/       清单、台账基线、统计与验收记录
scripts/    全部脚本
logs/       每次运行的日志
```

## 跑法

```bash
bash scripts/setup_env.sh                 # 装 mmseqs / cd-hit / gemmi / biopython，写版本锁

python3 scripts/00_rcsb_index.py          # RCSB 索引 → meta/rcsb_entries.jsonl
python3 scripts/05_bgsu_parse.py 4.61     # BGSU 等价类 → meta/bgsu_classes.tsv
python3 scripts/07_rna3db_meta.py         # RNA3DB 切分与 cmscan → processed/rna3db_chains.tsv
python3 scripts/40_casp_targets.py        # CASP 靶标 → meta/casp_holdout_pdb.json

python3 scripts/10_download_pdb.py --limit 50 --shuffle   # 先小样，且必须打散
python3 scripts/20_parse_structures.py --workers 10
python3 scripts/21_filter.py

python3 scripts/10_download_pdb.py --workers 16           # 全量
python3 scripts/20_parse_structures.py --workers 10
python3 scripts/21_filter.py

bash scripts/22_cluster.sh mmseqs                         # 五档聚类（100/90/70/50/30）
python3 scripts/81_split_struct.py --tier 90              # 三层约束划分
python3 scripts/24_bgsu_crosscheck.py                     # 与 BGSU 等价类对照
```

## 过滤漏斗（全量）

```
28,949 条 RNA 链
 → 26,165  方法与分辨率（X 射线 / 冷冻电镜 ≤ 4.0 Å，NMR 单列；1,599 条进扩展池）
 → 16,885  链长 20–10,000 nt（9,280 条短于 20 nt；上限按覆盖长链 rRNA 的要求放开）
 → 15,150  主链原子完整度 ≥ 90%（1,735 条剔除）
 → 13,993  规范碱基比例 ≥ 95%（1,157 条修饰过多）
 →  4,180  完全重复合并（合并 9,277 条，占 66%）
```

## 划分方式

分组单元 = 并查集合并四种关系：同一序列簇 ∪ 同一 BGSU 结构等价类 ∪ 同一 PDB 条目 ∪ 完全同序列。任一组完整落入一个集合。

用 **90% 序列档**而不是 30%：30% 档下长链 rRNA 会并成少数超大簇，val 只剩 1 条长链、test 29 条，长链性能无法监控；90% 档 test 长链升到 53 条。

## 已知限制（均为实测踩到，非推测）

1. CD-HIT-EST V4.8.1 的 `-c` 下限是 0.8，低于 0.8 直接 Fatal Error。70/50/30 档必须用 MMseqs2。
2. MMseqs2 在部分容器里误判可用内存（23 条序列也报 `Cannot fit databases into 1G`）。跑前先做 20 条 smoke test。
3. 不能用索引前 N 条做小样：RCSB 索引按发布日期倒序，前 N 条全是新结构，取样严重偏置。必须 `--shuffle`。
4. mmCIF 有两套链名（label_asym_id / auth_asym_id），实体声明序列登记的是 label，解析库 `chain.name` 返回 auth。直接按名字取序列会把序列贴错链，且不报错。必须建 label→auth 映射。
5. BGSU 成员写法里多链复合体用 `+` 连成一个单元（`11AA|1|1+11AA|1|2`），只按 `;` 拆会大面积假性不匹配。
6. 只按序列簇分集合会漏结构泄漏：同一结构的不同链可能跨集合。必须同时锁 PDB 条目与结构等价类。
7. 长期偏置：长度分布双峰，200–1000 nt 几乎空白；56% 的簇无 Rfam 家族归属（多为蛋白结合态短片段）。

## 未完成

- 长链评测集：val/test 长链占比 4.3% / 12.7%，语料整体 27%。rRNA 天然聚成少数大簇，需单独按释放日期构造长链 holdout。
- 30% 序列档作为敏感性口径（当前定稿为 90%）。
- RNAsolo 数据、rnaglib 基准包、RNA-Puzzles 的 PDB 映射未取得，原因见 source_registry.json。
