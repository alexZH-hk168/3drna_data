# RNA 3D 数据准备流水线

把台账里 12 条三维数据源接入，做成可复跑的「下载 → 解析 → 过滤 → 聚类 → 划分 → 审计」链路。

数据合同（阈值与字段定义）见 [数据合同.md](../数据合同.md)，源登记见 [meta/source_registry.json](../meta/source_registry.json)。

## 目录

```
raw/        原样下载（pdb/ rna3db/ bgsu/ isolated/ …）
interim/    解析中间件（chains.jsonl、聚类临时库）
processed/  最终数据集（chains_filtered.tsv、splits.tsv、clusters/）
meta/       清单、台账基线、统计与验收记录
scripts/    全部脚本
logs/       每次运行的日志
tools/      本地二进制（mmseqs、cd-hit）
```

## 建议机器

- CPU ≥ 16 核，内存 ≥ 64 GB（32 GB 为下限），数据盘 ≥ 500 GB SSD
- 实测：10,502 条 mmCIF 约 12 GB（50 条样本 58 MB，均值 1.16 MB/条）
- 本阶段不需要 GPU

## 跑法

```bash
# 0) 装环境（建目录、装 mmseqs / cd-hit / gemmi / biopython，写版本锁）
bash scripts/setup_env.sh

# 1) 索引（只读网络，不占磁盘）
python3 scripts/00_rcsb_index.py          # RCSB 含 RNA 条目 → meta/rcsb_entries.jsonl
python3 scripts/05_bgsu_parse.py          # BGSU 等价类 → meta/bgsu_classes.tsv

# 1.5) RNA3DB 元数据（不需要 mmCIF 也能先做：切分 + 家族注释）
#      先下 raw/rna3db 里的 rna3db_LICENSE / _releases.json / -jsons.tar.gz / -cmscans.tar.gz
python3 scripts/07_rna3db_meta.py         # → processed/rna3db_chains.tsv

# 2) 小样试通（先 50 条，别一上来全量）
python3 scripts/10_download_pdb.py --limit 50 --shuffle --workers 8
python3 scripts/20_parse_structures.py
python3 scripts/21_filter.py

# 3) 全量
python3 scripts/10_download_pdb.py --workers 16
python3 scripts/20_parse_structures.py
python3 scripts/21_filter.py

# 4) 聚类（高阈值用 CD-HIT，低阈值必须用 MMseqs）
bash scripts/22_cluster.sh cdhit          # 100% / 90%
bash scripts/22_cluster.sh mmseqs         # 70 / 50 / 30%

# 5) 划分与泄漏审计
python3 scripts/23_split.py --tier 30 --time-cutoff 2025-01-01

# 6) 结构级对照（自建聚类 vs BGSU 等价类）
python3 scripts/24_bgsu_crosscheck.py      # → meta/bgsu_crosscheck.json
```

## 各步实测（2026-10-20，云电脑 4 核环境小样）

| 步骤            | 输入                          | 输出                          | 实测                  |
| ------------- | --------------------------- | --------------------------- | ------------------- |
| 00 RCSB 索引    | —                           | 10,502 条                    | 22 批 GraphQL，约 2 分钟 |
| 05 BGSU       | 8 个 CSV                     | 5,212 类 / 22,887 链          | 秒级                  |
| 07 RNA3DB 元数据 | cmscan 88 MB + jsons 6.9 MB | 2,199 簇 / 890 个 Rfam 家族     | 秒级                  |
| 10 下载         | 300 条（打散抽样）                 | 462 MB（累计 349 个文件）          | 58.5 秒，5.1 条/秒      |
| 20 解析         | 349 个文件                     | 941 条 RNA 链                 | 0 解析错误              |
| 21 过滤         | 941 链                       | 239 链                       | 见下                  |
| 22 CD-HIT     | 239 序列                      | 208 簇(100%) / 174 簇(90%)    | 秒级                  |
| 24 BGSU 对照    | 239 链                       | 4.0 Å 档覆盖 89.5%、all 档 96.7% | 秒级                  |

过滤漏斗（349 结构 / 941 链，阈值 4.0 Å / 20–1000 nt / 主链 90%）：

```
941 条 RNA 链
 → 869  方法与分辨率（52 条进扩展池、20 条分辨率不达标）
 → 459  链长（271 条太短、139 条太长）
 → 459  主链完整性（无一被剔）
 → 418  规范序列比例（41 条修饰过多）
 → 239  完全重复合并（合并掉 179 条，占 43%）
```

过滤后 239 条的构成（`meta/distribution_stats.json`）：

- 长度：20–50 nt 85 条、50–100 nt 84、100–200 nt 55、200–400 nt 7、400–700 nt 5、700–1000 nt 3。**94% 在 200 nt 以下。**
- 方法：cryo-EM 149、X-ray 75、NMR 15。
- 分辨率：中位数 3.03 Å；2.5–3.0 Å 65 条、3.0–3.5 Å 91 条、3.5–4.0 Å 36 条。
- 链组成：含蛋白复合物 187 条、纯 RNA 51 条、含 DNA 1 条。**78% 是蛋白结合的 RNA。**

## 已知限制（都是实测踩到的，不是推测）

1. **CD-HIT-EST V4.8.1 的 `-c` 下限是 0.8**，低于 0.8 直接 Fatal Error。所以 70/50/30 三档必须用 MMseqs2，CD-HIT 只能做 100/90 的交叉验证。
2. **MMseqs2 在部分容器里会误判可用内存**：23 条序列也报 `Cannot fit databases into 1G`。目标机的 64 GB 内存下应正常，跑之前先用 20 条序列做 smoke test 确认。
3. **不要用索引前 50 条做小样**：RCSB 索引默认按发布日期倒序，前 50 条全是 2026 年新结构，取样严重偏置（会导致时间切分把全部样本都划进 holdout）。用 `--shuffle`。
4. **RNA3DB 的数据许可未解决**：仓库随附 MIT LICENSE，但正文通篇写 Software，属代码许可；底层 PDB 数据的再分发仍无声明。
5. **RNA3DB 的家族留出集中在一类结构上**：测试集 100 个无 Rfam 命中的簇全部落在同一个结构连通分量 component_0。照抄它的切分方式会得到一个覆盖很窄的“未见过家族”测试集，我们自己的划分要把留出集合的多样性单独报出来。
6. 小样只有 23 条时，`23_split.py` 的 test 集可能为空——贪心分配在小样本上会先把 train/val 填满。全量不受影响。
7. **BGSU 的成员会写成多链复合体**：格式是 `PDB|model|chain`，但多链时用 `+` 连成一个单元，如 `11AA|1|1+11AA|1|2`。只按 `;` 拆会把这整串当一个 token，导致大批链"对不上"。必须 `+` 和 `;` 都拆。
8. **只按 cluster 分集合会漏一种泄露**：同一 PDB 结构的不同链可能落进不同 cluster，于是同一个结构同时出现在 train 和 test（239 条样本上实测 3 例）。`23_split.py` 已用并查集把共享条目的 cluster 并成一组（174 簇 → 139 组），审计里这一项必须恒为 0。
9. 以下源尚未打通批量入口：RNAsolo（前端壳页）、EMDB 批量检索（403）、PDB-REDO 条目清单（500 / 需鉴权）、ModelArchive（壳页）、RNANet（官网死）。详见源登记表。

## 已完成

- 已验证：索引、下载、解析、过滤漏斗、CD-HIT 高阈值聚类、划分与泄漏审计脚本。
- 未验证：全量下载、70/50/30 三档聚类（等 64 GB 机器）、RNA3DB 与 rnaglib 的下载与整合、Infernal 家族注释。
