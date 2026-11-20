# auto_RNA_Seq_analysis.v4 转录组自动化分析流程使用说明

> **文档版本**: 2026-10-15（基于迁移后的当前部署状态编写）
> **原作者**: 岳尧（yueyao，2019-2020，见 `2020-3-25 RNA-Seq自动化流程使用说明.pdf`）
> **迁移维护**: yuzhij（2026-07-16 从 liudy 账户迁移部署至本人账户）
> **本文档取代 2020 年 PDF 说明**：PDF 中的路径均为原作者旧服务器路径（`/media/sdc/yueyao/...`），已全部失效；本文档反映迁移后的真实部署。

---

## 1. 流程概述

本流程是基于 **Snakemake** 的 Bulk RNA-seq 全自动分析流程，从 clean FASTQ 开始，经过质控、比对、定量、差异分析，到富集分析（GO/KEGG/Reactome/GSEA/GSVA）、PPI 网络和聚类可视化，最终可一键生成整理好的结果目录（`Report/`）和 HTML 交互报告。

**核心设计**：
- 主控脚本（Python + click）按 YAML 配置动态生成 Snakemake 主文件
- 各功能模块拆分为独立规则文件（`rule/*.snakemake.py`），通过 YAML 中的 step 开关自由组合
- 支持任务并行调度、断点续跑（Snakemake 自带）

---

## 2. 部署位置与运行环境

| 项目 | 路径 |
|---|---|
| 流程根目录 | `/media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/` |
| 主控脚本 | `tools/auto_RNA-Seq_analysis.v4.py`、`tools/make_yaml.py` |
| Snakemake 规则 | `rule/`（20 个规则文件） |
| 各物种依赖 | `depends_<物种>_<参考基因组版本>/`（9 个目录） |
| 通用依赖脚本 | `depends_script/`（R/Perl/Python 脚本 + 各物种 KEGG .ko 文件） |
| 报告模板 | `report/part/`、`report/pipline/`、`report/image/` |
| Python 运行环境 | `/media/nfs/nfs02/yuzhij/biosoft/mambaforge/envs/RNA-seq/bin/python` |
| Snakemake 运行环境 | 同上（RNA-seq 环境内含 snakemake） |
| R 运行环境 | `/media/nfs/nfs02/yuzhij/biosoft/mambaforge/envs/RNA-seq-R/bin/R`（R 4.2.3） |

> **迁移说明**：`biosoft/mambaforge` 整个目录是从 liudy 账户复制过来的独立 Miniforge 安装，
> 其中 `RNA-seq`（Python 侧：snakemake/hisat2/stringtie/samtools/fastqc）和
> `RNA-seq-R`（R 侧：DESeq2/clusterProfiler/GSVA 等）两个环境已在迁移期间修复为可在本账户直接使用。

### 软件版本（实测）

| 软件 | 版本 | 用途 |
|---|---|---|
| Snakemake | 3.13.3 | 流程调度 |
| FastQC | 0.11.5 | 原始数据质控 |
| HISAT2 | 2.2.1 | 序列比对（`--dta --dta-cufflinks`） |
| StringTie | 1.3.3b | 转录本组装与定量 |
| samtools | 1.16.1 | BAM 排序/处理 |
| SOAP | 2.21release | Reads 在基因上的分布统计（ReadsRandom） |
| GSEA | 3.0（gsea-3.0.jar） | 基因集富集分析 |
| R | 4.2.3（RNA-seq-R 环境） | DESeq2 差异分析、GO/KEGG/Reactome 富集、GSVA、PPI、聚类 |

---

## 3. 快速开始（标准用法）

流程以"项目目录"为单位运行。每个项目目录（如 `<日期>/<项目名>/02.raw_FC2/`）内需要准备
`all.sh`（入口脚本）和 `sample.txt`（样本信息），然后依次执行三步。

### 第 1 步：准备项目目录

在项目目录中放置两个文件：

**`sample.txt`**（制表符分隔，3 列）：

```
<原始数据目录名>   <分组名>   <重命名后的样本名>
SRR4415752    Control    Control-1
SRR4415753    Control    Control-2
SRR4415758    Test       Test-1
comparison    Control    Test          ← 以 comparison 开头的行定义比较（第2列=对照，第3列=处理）
```

- 第 1 列：原始 FASTQ 所在子目录名（流程会到 `RAW_DATA_DIR/<第1列>/` 下找 `.fq/.fq.gz`；PE 数据找 `*_1.fq*` 和 `*_2.fq*`）
- 第 2 列：分组名（比较时用）
- 第 3 列：输出使用的样本名（**建议带组名前缀 + 数字后缀**，如 `Control-1`；后续 Cluster 步骤按 `-` 后的数字还原重复编号）
- comparison 行：生成 `对照_vs_处理` 的比较组合，可写多行

**`all.sh`**（从模板修改 4 个参数即可）：

```bash
SPECIES="fca"                      # 物种：hsa/mmu/rno/ssc/mfa/fca/cfa
SEQ_TYPE="SE"                      # SE（单端）或 PE（双端）
BATCH="20260813"                   # 测序批次（数据传输站目录名）
PROJECT="HTP_ER_PRF_2607_0130"     # 项目名称
THREADS=48                         # 总线程数（单样本比对固定用 12 线程）
RAW_DATA_DIR="/media/nfs/nfs03/seq/Data_transfer_station/${BATCH}/${PROJECT}/2.clean/"
```

> 模板位置：`/media/nfs/nfs02/yuzhij/Project/HTP/all.sh`（内含 7 个物种的 YAML 路径映射表）。
> 注意数据传输站有 nfs01 和 nfs03 两个入口，按实际填写。

### 第 2 步：生成配置并启动

```bash
cd <项目目录>
sh all.sh          # 依次完成：链接工具 → 生成 fca_complemented.yaml → 生成 run.sh/cp.sh
sh run.sh          # 运行 Snakemake 主流程（日志：logo.snakemake 标准输出 / loge.snakemake 错误）
```

### 第 3 步：整理结果与生成报告

```bash
sh cp.sh           # 生成 Report/ 整理目录
# 可选：HTML 报告
ln -s /media/nfs/nfs02/yuzhij/pipeline/auto_RNA_Seq_analysis.v4/report/part/RNAseq_report_html.sh .
sh RNAseq_report_html.sh    # 生成 RNA_seq_report.html（需要先有 Report/ 和 group.txt）
```

---

## 4. 配置文件机制（make_yaml.py）

`all.sh` 调用 `make_yaml.py`，将**物种模板 YAML + sample.txt** 合成为当前项目专用的
`<物种>_complemented.yaml`：

```
python make_yaml.py yaml -s <物种> -e <SE|PE> -y <物种模板yaml> -r <raw数据路径> -f <foldchange> -sa sample.txt
```

| 参数 | 说明 |
|---|---|
| `-s` | 物种（决定使用哪个 depends 目录） |
| `-e` | SE/PE |
| `-y` | 物种模板 YAML 路径（all.sh 中已按物种映射好） |
| `-r` | 原始数据路径 |
| `-f` | 差异分析 foldchange 阈值（默认 2） |
| `-sa` | 样本信息文件（默认 sample.txt） |

生成的 `<物种>_complemented.yaml` 结构（5 个部分）：

1. **database**：HISAT2 索引（genome）、SOAP 基因组索引（gene）、GTF、基因描述文件（des）、depends 目录、PPI 文件路径
2. **software**：所有软件的绝对路径（指向本账户 mambaforge 环境）
3. **step**：12 个步骤的开关（1=运行，0=跳过）
4. **comparison / group / DEseq2**：比较组合、分组样本、差异阈值（foldchange、padj 0.05）
5. **samples**：样本名 → FASTQ 绝对路径映射
6. **outdir**：固定为运行时的工作目录（make_yaml.py 自动改写）

> 如需重跑部分步骤或调整参数，直接编辑 `<物种>_complemented.yaml` 后重新执行
> `auto_RNA-Seq_analysis.v4.py makeshell -c <物种>_complemented.yaml -t <线程数>` 再 `sh run.sh`。

---

## 5. 分析步骤详解（12 步）

| # | 目录 | 步骤 | 主要软件 | 输入 → 输出 |
|---|---|---|---|---|
| 0 | `0.FastQC/` | 质控 | FastQC | FASTQ → html/zip 质控报告 |
| 1 | `1.ReadsRandom/` | Reads 基因分布 | SOAP + ReadsRandomInGene.pl | FASTQ → soap 比对列表 + Reads_distribution 分布图 |
| 2 | `2.Alignment/` | 基因组比对 | HISAT2 + samtools | FASTQ → `{sample}.bam` + `{sample}.log`（含 alignment rate） |
| 3 | `3.Quantify/` | 定量 | StringTie + prepDE.py | BAM → 每样本 abundance/gtf；汇总 FPKM/TPM 矩阵、count 矩阵、ensembl2symbol 映射 |
| 4 | `4.Comparison/` | 差异分析 | DESeq2（R） | count 矩阵 → 每个比较的 `*.differentially_expressed_genes.xls`、`total_genes.xls`、DEG 统计图；去重复/去低表达过滤矩阵 |
| 5 | `5.Cluster/` | 聚类 | R（PCA/层次聚类） | normalized counts → `PCA.<比较>.pdf/png`、`hclust.<比较>.*`、`<比较>.fpkm.txt`（供 GSEA） |
| 6 | `6.GSEA/` | GSEA | gsea-3.0.jar（Java） | fpkm + cls 分组文件 → GSEA 富集结果 |
| 7 | `7.KEGG/` | KEGG 富集 | clusterProfiler（R） | DEG → `*.kegg_enrichment.xls`、heatmap、pathview 通路图 |
| 8 | `8.GO/` | GO 富集 | R | DEG → MF/BP/CC 三大类 `*_enrichment.xls`、heatmap |
| 9 | `9.Reactome/` | Reactome 富集 | R | DEG → `*.reactome_enrichment.xls`、heatmap |
| 10 | `10.PPI/` | 蛋白互作 | STRING links + PPI.R | DEG → 互作网络 `*.network.relation.xls`、top100、网络图 |
| 11 | `11.GSVA/` | 基因集变异 | GSVA + limma（R） | normalized counts → `all.gsva_score.xls`、每比较的差异通路 |

**规则文件与步骤的对应**（`rule/` 目录）：

- `FastQC.IMA.{paired,single}.snakemake.py` — 步骤 0（按 SE/PE 选择）
- `SOAP.IMA.{paired,single}.snakemake.py` — 步骤 1
- `Align.IMA.{paired,single}.snakemake.py` — 步骤 2（每样本固定 12 线程）
- `Quantify.IMA.snakemake.py` — 步骤 3（stringtie / mergefpkm / all_deseq_normal 三个规则）
- `Comparison.IMA.v2.snakemake.py` — 步骤 4（pre_comparison / deg / add_des 三个规则）
- `Cluster.IMA.v2.snakemake.py` — 步骤 5
- `GSEA.IMA.snakemake.py` — 步骤 6（pregsea 生成 cls / gsea）
- `GO.IMA.snakemake.py`、`KEGG.IMA.snakemake.py`、`REACTOME.IMA.snakemake.py` — 步骤 7-9（各自 pre 规则复制脚本 + 主规则）
- `PPI.IMA.snakemake.py` — 步骤 10
- `GSVA.IMA.snakemake.py` — 步骤 11（gsva_all / gsva_deg）

> 注意：`Comparison.IMA.snakemake.py`（非 v2）和 `Cluster.IMA.snakemake.py`（非 v2）为旧版备份，
> 当前主控生成的主 Snakefile 引用的是 **v2 版本**。

### 运行日志与计时

- `log/`：各步骤的 `.e/.o` 错误/标准输出
- `time/`：各任务的开始/结束时间戳
- `logo.snakemake` / `loge.snakemake`：Snakemake 全局标准输出/错误输出（**排查报错先看 loge.snakemake**）

---

## 6. 输出目录结构

```
<项目目录>/
├── 0.FastQC/<sample>/            # 每样本 FastQC 报告（可用 delete_dirs.sh 清理）
├── 1.ReadsRandom/                # soap 比对 + Reads_distribution/
├── 2.Alignment/                  # bam + 比对日志（可用 delete_dirs.sh 清理）
├── 3.Quantify/
│   ├── DESeq2/<sample>/          # StringTie 每样本 abundance + gtf
│   ├── all.genes.fpkm.xls / all.genes.tpm.xls / all.ensembl.*.xls
│   ├── gene_symbol_count_matrix.csv / gene_count_matrix.csv / transcript_count_matrix.csv
│   └── ensembl2symbol.txt
├── 4.Comparison/
│   ├── comparisons.txt           # 分组还原后的比较列表
│   ├── <比较>.differentially_expressed_genes.xls      # DEG 列表（含描述）
│   ├── <比较>.total_genes.xls                         # 所有基因差异结果
│   ├── <比较>.DEG_stat.png / All_comparison.DEG_stat.png / DEG_stat.xls
│   └── all.genes.no_dup_exp.fpkm.xls 等过滤后矩阵
├── 5.Cluster/                    # PCA/hclust（每比较 + allsamples）
├── 6.GSEA/  7.KEGG/  8.GO/  9.Reactome/  10.PPI/  11.GSVA/
├── Report/                       # cp.sh 生成的整理目录
│   ├── 1.QC/（FastQC、Reads_distribution）
│   ├── 2.Cluster/（PCA、Hierarchical_clustering）
│   ├── 3.DEG/
│   ├── 4.Statistics/（Alignment_rate.xls、DEG_stat）
│   ├── 5.Expression/（fpkm/tpm/normalized counts）
│   └── 6.Function/（GSEA、KEGG、GO、Reactome、PPI、GSVA）
├── RNA_seq_report.html           # 可选 HTML 报告（RNAseq_report_html.sh 生成）
├── run.sh / cp.sh / RNAseq.snakemake.py / <物种>_complemented.yaml
└── logo.snakemake / loge.snakemake / log/ / time/
```

---

## 7. 支持的物种与 depends 目录

| 物种代码 | depends 目录 | 参考基因组 | 注释来源 |
|---|---|---|---|
| hsa | `depends_hsa_GRCh38.p13_Ensembl_release_89` | GRCh38.p13 | Ensembl R89 |
| mmu | `depends_mmu_GRCm38.p6_Ensembl_release_89` | GRCm38.p6 | Ensembl R89 |
| mmu | `depends_mmu_GRCm38.p6_Ensembl_release_96` | GRCm38.p6 | Ensembl R96（备用） |
| rno | `depends_rno_Rnor_6.0_Ensembl_release_90` | Rnor_6.0 | Ensembl R90 |
| ssc | `depends_Sscrofa11.1_NCBI_release_106` | Sscrofa11.1 | NCBI R106 |
| mfa | `depends_Mfa_MFA1912RKSv2_NCBI_release_102` | MFA1912RKSv2 | NCBI R102 |
| fca | `depends_F.catus_Fca126_mat1.0_NCBI_release_105` | Felis_catus 9.0 | NCBI R105 |
| cfa | `depends_ROS_Cfam_1.0_NCBI_release_106` | ROS_Cfam_1.0 | NCBI R106 |
| mcc | `depends_Mmul_10_NCBI_release_103`（目录存在，all.sh 未映射） | Mmul_10 | NCBI R103 |

每个 depends 目录内含：该物种的 `*.yaml` 模板、`gene_description.txt`、KEGG GMT（`*.kegg_optimal.symbol.gmt`）、
STRING PPI 文件（`rename.uniq.<tax_id>.protein.links.*.txt`）、`mRNA_length.txt` 以及从 `depends_script/`
链接过来的分析脚本。

### 新增物种

参考 `new_depends.sh`（原作者的建目录脚本）：新建 `depends_<物种>` 目录，从 `depends_script/`
软链各分析脚本，放入该物种的基因描述/KEGG GMT/STRING/PPI 文件和 yaml 模板，再在 `all.sh` 的
`YAML_PATHS` 中加一行映射。

---

## 8. 辅助脚本（`/media/nfs/nfs02/yuzhij/Project/HTP/`）

| 脚本 | 用途 |
|---|---|
| `all.sh` | 项目运行入口模板（含物种→YAML 映射） |
| `samples_process.py` | 从《高通量项目执行表-*.xlsx》自动生成 sample.txt |
| `png_make.py` / `auto_png.py` | 生成 HTML 报告用的样本信息图/分析策略图 |
| `delete_dirs.sh` | 分析完成后批量删除大体积中间目录（0.FastQC/1.ReadsRandom/2.Alignment）；默认 dry-run，`bash delete_dirs.sh false` 才真删 |
| `HTP_0030_step.sh` | 分步执行示例 |

---

## 9. 迁移历史与已知问题（重要）

### 9.1 迁移时间线（2026-07-16 从 liudy 账户迁移）

1. 复制 `pipeline/auto_RNA_Seq_analysis.v4`（流程本体）和 `biosoft/mambaforge`（软件环境）
2. 修复 mambaforge 环境的 R 路径问题（`.conda/envs/RNA-seq-R/bin/R` 报"没有那个文件或目录"，
   原因是复制的 conda 环境内部路径前缀失效，已修复为 `/media/nfs/nfs02/yuzhij/biosoft/mambaforge/...`）
3. 修改各物种 yaml 中的软件/数据库路径指向本账户
4. 2026-07-24 在 `Project/test/HTP_ER_PRF_2512_0233/02.raw_FC2` 验证通过
5. 2026-08-20~22 在 `Project/20260813/HTP_ER_PRF_2607_0130` 实战运行，期间修复：
   - `sample.txt` 样本名/分组问题
   - `fca_complemented.yaml` 中缺失的 STRING PPI 文件：liudy 目录缺
     `rename.uniq.9685.protein.links.v10.5.txt`，已从 STRING 官网下载 **v12.0** 版替代
     （`depends_F.catus.../rename.uniq.9685.protein.links.v12.0.txt`）
   - 基因 symbol 矩阵修复（`fix_gene_symbol_matrices.R`，用后已不在项目目录中）
   - `Quantify.IMA.snakemake.py`、`Comparison.IMA.v2.snakemake.py`、`Cluster.IMA.snakemake.py`
     等规则文件在 7-8 月迁移期间有修改（保留有 `.backup/.bak` 旧版可对比）

### 9.2 ✅ 已修复（2026-10-20）：liudy/旧服务器依赖彻底清除

**原问题**：迁移时 `depends_*` 各物种目录是按目录复制的，每个目录有 21-27 个符号链接指向
liudy 账户（或已废弃的旧服务器 /media/sda、/media/sdc）；另外 hsa/mfa/cfa/ssc 四个物种的
yaml 模板整个还是 liudy 路径，mmul/mmu_96 两个备用 yaml 还是原作者旧服务器路径，
ssc.yaml 还存在文件末尾缺换行符的原版 bug（会导致猪物种配置生成失败）。

**2026-10-20 已全部修复**，详见 `迁移修复报告_20261020.md`：
- 197 个符号链接改指本地 `depends_script/`（其中 147 个内容经 md5 校验一致，44 个死链复活）
- 6 个数据文件替换为本地库真副本（md5 校验一致），猕猴 mcc gmt 从本地 KEGG 库恢复
- 狗的 STRING PPI 文件 `rename.uniq.9615.protein.links.v11.5.txt` 从本地 STRING v11.5
  raw+info 重新生成（格式与 fca 2026-08 实战跑通文件一致）
- 9 个 yaml 全部修复为本地路径（原文件备份为 `*.bak_liudy_20261020` / `*.bak_oldserver_20261020`）
- 修复后验证：0 个 liudy/旧服务器链接、0 断链；9 物种假数据端到端 Snakemake 干跑全部通过；
  15 个 R + 3 个 Perl + 1 个 Python 脚本语法全部通过

**回滚方式**：符号链接按 `修复备份_symlink清单_20261020.tsv` 重建；yaml 用各目录下的 .bak 备份。

> mmu_96 备用目录说明：原 R96 注释索引无任何本地副本，其 yaml 的 genome/gtf 已改用本地
> release-89 等价路径（yaml 内有注释），如确需 R96 需另行下载重建。

### 9.3 其他注意点

- **0/1/2 目录清理**：跑完并确认结果无误后可用 `delete_dirs.sh` 删除中间大文件（默认 dry-run 安全模式）
- **比对线程**：Align 规则每样本固定 12 线程，`THREADS` 控制的是 Snakemake 并行任务数上限
- **`-k` 参数**：run.sh 中 snakemake 带 `-k`（某任务失败仍继续其他任务），排错时看 `loge.snakemake` 末尾
- **padj 阈值**：固定 0.05；foldchange 由 `all.sh` 的 `-f` 参数控制（默认 2，即 |log2FC|≥1）
- **GSEA 需要 Java**：使用系统 `/usr/bin/java` + `biosoft/gsea-3.0.jar`
- **gsva 脚本路径**：所有物种 yaml 的 `gsva` 均指向 `depends_F.catus...` 目录（gsva 脚本与物种无关，历史遗留写法，勿动）

---

## 10. 常见报错排查（迁移期间实际遇到）

| 报错 | 原因 | 解决 |
|---|---|---|
| `.conda/envs/RNA-seq-R/bin/R: 没有那个文件或目录` | 复制的 conda 环境路径前缀失效 | 已修复（迁移期间）；若再现检查 mambaforge/envs/RNA-seq-R/bin/ 实际文件 |
| PPI 步骤报缺 `rename.uniq.*.protein.links.v10.5.txt` | liudy 目录本身缺该物种 STRING 文件 | 从 STRING 数据库下载对应物种 links 文件，改名放入对应 depends 目录（fca 已用 v12.0 解决） |
| 样本找不到 FASTQ | sample.txt 第 1 列与 RAW_DATA_DIR 下子目录名不一致 | 核对数据目录名与扩展名（`.fq.gz`），PE 数据需有 `_1/_2` 后缀 |
| Cluster 图样本编号错乱 | 样本名不符合 `分组-数字` 命名规则 | 统一用 `<组名>-<数字>` 命名（pre_comparison 按首个 `-` 后数字还原重复） |
| 基因 symbol/ensembl 混乱 | GTF 版本与基因描述文件不匹配 | 用修复脚本重新生成映射（见 9.1 迁移历史第 5 条） |
| Snakemake 中断后重跑卡住 | 旧版本 Snakemake 3.13.3 的锁文件 | 删除 `.snakemake/` 目录后重跑 |

---

## 11. 与 2020 年 PDF 旧说明的主要差异

| 旧 PDF（2020） | 当前状态（本文档） |
|---|---|
| 直接手动调用 make_yaml.py + 主脚本 | 提供 `all.sh` 一键封装（含物种映射） |
| 路径 `/media/sdc/yueyao/...`、`module load snakemake` | 全部部署在 `/media/nfs/nfs02/yuzhij/` 自有账户，环境自带 |
| 支持 mmu/hsa/rno/ssc 4 物种 | 7 物种已映射（hsa/mmu/rno/ssc/mfa/fca/cfa），另 2 个目录备用 |
| run.sh 需手写 snakemake 命令 | 主控脚本自动生成 run.sh/cp.sh |
| 无报告整理说明 | cp.sh → Report/ 结构 + RNAseq_report_html.sh HTML 报告 |
| 无迁移/故障排查内容 | 第 9-10 节完整记录迁移修复与排查经验 |

---

## 附：参考项目实例

- 验证项目：`/media/nfs/nfs02/yuzhij/Project/test/HTP_ER_PRF_2512_0233/02.raw_FC2/`
- 实战项目（猫，SE，fca）：`/media/nfs/nfs02/yuzhij/Project/20260813/HTP_ER_PRF_2607_0130/02.raw_FC2/`
  （含完整输出 0-11 步骤、Report/、RNA_seq_report.html、ThreeGroup_Comparison 三组比较扩展分析）
- 模板目录：`/media/nfs/nfs02/yuzhij/Project/HTP/`
