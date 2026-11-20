# 转录组分析流程（Bulk RNA-seq）

服务器部署代码的分类归档，快照日期 **2026-11-20**。clean FASTQ → FastQC → HISAT2 → StringTie → DESeq2 → GO/KEGG/Reactome/GSEA/GSVA/PPI → HTML 报告。

## 使用入口

```bash
# 使用与服务器部署兼容的 RNA-seq Python / Snakemake 环境。
python tools/make_yaml.py --help
python tools/auto_RNA-Seq_analysis.v4.py --help
```

实际项目按 [服务器使用手册](README_流程使用说明.md) 的三步执行。
`config/species/` 保存物种 YAML 配置模板；必须修改其中的 `ruledir`、软件与数据库路径，再通过 `make_yaml.py yaml -y /绝对路径/物种.yaml ...` 指定。
规则文件是 Snakemake DSL，不是可直接执行的普通 Python 程序。
原部署采用旧版 Snakemake 3.13.3；不能假定新版语法兼容。

## 目录

核心目录：`tools/、rule/、depends_script/、report/、config/species/`。

- [README_流程使用说明.md](README_流程使用说明.md)：原部署使用说明。
- [部署与资源](docs/DEPLOYMENT.md)：迁移前需配置的路径、软件和外部资源。
- [验证状态](docs/VALIDATION.md)：本次实际执行的检查及未执行事项。
- [来源与署名](NOTICE.md)：来源、修改范围及第三方材料说明。
- [源文件校验清单](docs/source_manifest.json)：已归档文件与原服务器快照的SHA-256。

## 三类组学导航

- [转录组分析流程（Bulk RNA-seq）](https://github.com/black662189/bulk-rnaseq-pipeline)
- [定量蛋白质组下游分析流程](https://github.com/black662189/quantitative-proteomics-pipeline)
- [非靶向代谢组分析流程](https://github.com/black662189/untargeted-metabolomics-pipeline)

## 版本范围

本次发布为部署源码快照，不承诺任意机器零配置运行。大型数据库、参考基因组、索引、Conda/venv环境、真实样本矩阵及项目结果未上传。
服务器原目录未修改；未为上传重新运行真实项目。使用现有说明时，请将服务器绝对路径替换为自己的部署路径。
