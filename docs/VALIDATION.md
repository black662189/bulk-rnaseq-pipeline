# 发布验证记录

验证日期：2026-11-20。验证对象为从服务器复制出的代码快照，在独立临时目录运行；未改动生产流程或真实项目结果。

- 发布前检查未发现私钥、GitHub令牌或云访问密钥。
- 已归档源文件SHA-256与服务器快照一致；源码统计逻辑未修改。
- shell脚本通过 `bash -n`；R/r脚本通过R语法解析。
- 未重新跑真实数据全流程，未验证全物种、全模块或新机器部署。
- 不把历史VALIDATION文档中的结果冒充本次测试。

## 本次结果

- 两个Python CLI入口的 `--help` 通过。
- 5个Python 3源码文件语法解析通过。
- `depends_script/PPI/fd_relaion.py`保留Python 2语法；规则明确用`params.python2`调用，已用服务器Python 2编译检查通过。这不是Python 3兼容脚本。
- 19个规则文件属于Snakemake DSL，没有用Python AST误判；本次未跑Snakemake全流程。
- 参考基因组不同版本的YAML分目录保存，避免覆盖同名物种配置。

Git发布时统一文本换行为LF，源清单保留原始服务器字节校验值；publication_sha256.json记录Git索引中发布文件的SHA-256（不含清单自身）。换行标准化不改变统计实现。
