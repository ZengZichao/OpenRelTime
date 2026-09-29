# OpenRelTime documentation / 文档目录

Every document listed below exists in **both** languages, with the English file
as the canonical one and the Chinese counterpart carrying the `-zh` suffix; each
file links to its counterpart at the top.  This index itself is bilingual.  The
project `README.md` has its counterpart one level up (`README-zh.md`).

下列每一篇文档都有**中、英两个版本**：以英文文件为准，中文版以 `-zh` 结尾，每篇
文档开头都有一条指向另一语言版本的链接。本目录本身中英合排。项目根目录的
`README.md` 对应上一层目录的 `README-zh.md`。

## Manual / 使用手册

| | English | 中文 | 内容 / Content |
|---|---|---|---|
| User guide | [usage-en.md](usage-en.md) | [usage-zh.md](usage-zh.md) | Installation, inputs, every subcommand, calibration format, output inventory, desktop-GUI pointer, FAQ, troubleshooting 安装、输入准备、全部子命令、校正格式、输出清单、图形界面指引、FAQ 与故障排查 |
| Parameter handbook | [parameters.md](parameters.md) | [parameters-zh.md](parameters-zh.md) | Semantics, default, range and recommendation for **every** option, plus a per-subcommand option matrix 全部可选参数的语义、默认值、取值范围与建议，含各子命令选项对照表 |
| Methods | [methods.md](methods.md) | [methods-zh.md](methods-zh.md) | What the RRF computes, relation to MEGA RelTime, calibrations, CIs, assumptions and limitations 方法原理、与 MEGA 的关系、校正、置信区间、假设与局限 |
| Validation | [validation.md](validation.md) | [validation-zh.md](validation-zh.md) | What the test suite actually proves, gate-by-gate; shipped-but-unvalidated modules; tested vs claimed platforms 测试套件实际证明了什么、逐门槛地图、未验证模块、已测与仅声明的平台 |

## Tutorials / 教程

| English | 中文 |
|---|---|
| [Tutorial 1 — five minutes with OpenRelTime](tutorials/tutorial1-quickstart.md) | [教程 1 —— 五分钟上手](tutorials/tutorial1-quickstart-zh.md) |
| [Tutorial 2 — end-to-end dating with calibrations and CIs](tutorials/tutorial2-end-to-end.md) | [教程 2 —— 带校正与置信区间的端到端定年](tutorials/tutorial2-end-to-end-zh.md) |
| [Tutorial 3 — MCMCTree priors with ddBD and the MEGA-CC bridge](tutorials/tutorial3-priors-megacc.md) | [教程 3 —— ddBD 树先验与 MEGA-CC 桥接](tutorials/tutorial3-priors-megacc-zh.md) |

## Architecture records / 架构决策记录

| English | 中文 | 内容 / Content |
|---|---|---|
| [adr/README.md](adr/README.md) | [adr/README-zh.md](adr/README-zh.md) | ADR index (ADR-001 … ADR-005) 架构决策记录索引 |
| [adr/ADR-004-tree-operations.md](adr/ADR-004-tree-operations.md) | [adr/ADR-004-tree-operations-zh.md](adr/ADR-004-tree-operations-zh.md) | Outgroup rooting + pruning semantics 外群置根与剔除的语义 |

## Release and provenance records / 版本与来源记录（仓库根目录）

| English | 中文 | 内容 / Content |
|---|---|---|
| `CHANGELOG.md` | `CHANGELOG-zh.md` | Release notes 版本记录 |
| `THIRD-PARTY-NOTICES.md` | `THIRD-PARTY-NOTICES-zh.md` | Reference material and adapted constants, per source 参考材料与被改造的常数，逐来源说明 |
| `data/PROVENANCE.md` | `data/PROVENANCE-zh.md` | Origin, licence and citation of every bundled third-party file 随包第三方文件的来源、许可与引用 |
| `data/examples/README.md` | `data/examples/README-zh.md` | What the bundled example trees, calibration and classification tables are, and which of their numbers are illustrative only 随包示例数据说明，以及其中哪些数值只是示意值 |

## Desktop GUI lives elsewhere / 图形界面文档不在此处

The point-and-click front end, **OpenRelTime Studio**, is an independently
installed product that depends on this engine's public API.  It ships as the
separate distribution `OpenRelTime-Studio`
(<https://github.com/ZengZichao/OpenRelTime-Studio>), and its user manual lives
in that project.  Nothing below documents its windows; `usage-en.md` §7 /
`usage-zh.md` §7 only point at it.

图形界面程序 **OpenRelTime Studio** 是一款独立安装的产品，把本引擎的公开 API 当作
普通依赖使用，它的使用手册也在该项目内。下面这份索引里的文档不描述界面细节，只有
§7 一节给出指引。

## Distribution / 分发方式

> **Distribution.** This release is distributed as source. Install with
> `pip install --editable ".[dev,plot]"` from a checkout; an index release will
> be announced here.

> **分发方式。** 本版本以源码形式分发，请在检出目录中用
> `pip install --editable ".[dev,plot]"` 安装；发布到包索引后会在此处公告。

## Which document to read first / 阅读顺序建议

New user 新用户 → user guide 使用说明 → tutorial 1 教程 1.
Choosing options 需要挑参数 → parameter handbook 参数手册.
Assessing a result 复核结果 → methods 方法说明 + validation 验证说明.
Regenerating the reference data 重新生成参考标准数据 → `reproduce/golden_r3f.R`.
