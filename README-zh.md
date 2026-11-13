# OpenRelTime

**语言：中文** · [English](README.md)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23053145.svg)](https://doi.org/10.5281/zenodo.23053145)
[![PyPI](https://img.shields.io/pypi/v/openreltime)](https://pypi.org/project/openreltime/)
[![Python](https://img.shields.io/badge/python-3.10--3.13-blue)](https://pypi.org/project/openreltime/)
[![CI](https://github.com/ZengZichao/OpenRelTime/actions/workflows/ci.yml/badge.svg)](https://github.com/ZengZichao/OpenRelTime/actions/workflows/ci.yml)
[![License: GPL-3.0](https://img.shields.io/badge/License-GPL--3.0-blue.svg)](LICENSE)
[![ORCID](https://img.shields.io/badge/ORCID-0000--0001--6553--970X-a6ce39?logo=orcid&logoColor=white)](https://orcid.org/0000-0001-6553-970X)

**OpenRelTime** 是用于**非贝叶斯**分子定年的**相对速率框架**（Relative Rate
Framework, RRF；Tamura et al. 2018, *MBE* 35:1770-1782）的开源 Python 实现；
该框架即 MEGA 中 RelTime 所依据的方法。

软件从带枝长的系统树出发，估计谱系相对速率与节点分化时间。用户给出的校正点
（硬边界或概率密度）会转换为绝对分化时间，并给出解析置信区间。软件还内置三项
扩展功能：CorrTest 速率自相关检验、ddBD 物种分化树先验估计、BLB（bag of
little bootstraps，轻量自助重采样）管线。所有功能都提供 Python API 与命令行
两种入口。

- 中文文档（完整手册均有中文版）：[文档目录](docs/README.md) ·
  [使用说明](docs/usage-zh.md) · [参数手册](docs/parameters-zh.md)
- English documentation: [README.md](README.md) ·
  [docs/README.md](docs/README.md)（使用说明、参数手册、方法、验证、教程）·
  [Usage guide](docs/usage-en.md)
- Python API 与命令行完全等价；金标准（golden）参照数据可复现；采用 GPL-3.0 许可

```
@software{openreltime,
  author    = {曾, 子超},
  title     = {OpenRelTime: relative rate framework based molecular dating in Python},
  year      = {2026},
  version   = {v0.1.1},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.23053145},
  url       = {https://github.com/ZengZichao/OpenRelTime}
}
```

## 为什么需要 OpenRelTime

| 工具 | 语言 | 开源 | 校正 | 解析置信区间 | CorrTest | ddBD | API 与 CLI |
|---|---|---|---|---|---|---|---|
| MEGA 11 与 12（RelTime） | Pascal | GPL-3 | ✓ | ✓ | ✓ | – | 弱 |
| R3F | R | GPL-3 | – | – | ✓ | ✓ | – |
| **OpenRelTime** | **Python** | **GPL-3** | **✓（边界+密度）** | **✓** | **✓** | **✓** | **✓** |

OpenRelTime 在 274 分类群哺乳动物基准树上与 R3F 数值一致：回归斜率 1.0000，速率与
时间的实测中位相对偏差分别为 6.8e-16 和 8.7e-16。自动验收门槛 G2a 与 G2b 的判据为
斜率 ≥ 0.999、中位相对偏差 ≤ 1e-6。实测最大相对偏差属于观测值，不在该门槛的判定
范围内，其中速率约 4.8e-15、时间 3.6e-14，详见 `docs/validation-zh.md`。在共享
核心之上，OpenRelTime 补齐了 R3F 缺失的校正转换与解析置信区间（msz236），并
提供端到端的 BLB 管线。

## 安装

推荐从 PyPI 安装（`[plot]` 会带上 matplotlib，供 `--plot` 使用）：

```bash
pip install openreltime             # 引擎（numpy、scipy、pandas、click）
pip install "openreltime[plot]"     # + matplotlib 可视化
```

从检出安装（`dev` 为开发者加装 pytest/ruff/mypy）：

```bash
git clone https://github.com/ZengZichao/OpenRelTime
cd OpenRelTime
pip install --editable ".[dev,plot]"   # 核心（numpy、scipy、pandas、click）
                                       # + matplotlib 可视化（[plot]）
                                       # + pytest/ruff/mypy（[dev]，开发者）
```

桌面图形界面是独立产品，另行安装（见下文）。

`plot` 与 `dev` 是本包声明的全部可选依赖组，需要多组时合并书写，例如
`"openreltime[plot,dev]"`。

各版本的 sdist 与 wheel 归档也继续挂在对应的 GitHub Release 上（如
[v0.1.1](https://github.com/ZengZichao/OpenRelTime/releases/tag/v0.1.1)），
供溯源使用。

环境要求：Python ≥ 3.10，支持 Linux、macOS 和 Windows。项目不含编译扩展，
安装只需数秒。CI 在 Linux（`ubuntu-latest`）上跨 Python 3.10、3.11、3.12
与 3.13 运行测试套件，在 macOS（`macos-latest`）上**只在 Python 3.13** 上运行
——`.github/workflows/ci.yml` 的矩阵排除了 macOS 上较旧的三个解释器，以缩
短排队时间，所以 macOS 只是单解释器检查，而不是 3.10～3.13 的覆盖声明。
Windows（`windows-latest`）以一个独立作业在 Python 3.13 上同时运行引擎套件与端到端
套件；该作业自 0.1.1 加入，因此其首个证据要等下一次 Actions 运行（已测与仅声明的对照、
以及套件覆盖不到的一条 Windows 局限，见 `docs/validation-zh.md`）。

## OpenRelTime Studio（桌面图形界面）——独立产品

**OpenRelTime Studio** 是一款原生桌面程序，面向偏好图形化点选操作的研究者：加载
Newick 或 NEXUS 树，运行 RRF，在画布上点选内部节点设置校正点，然后依次完成校
正、解析置信区间、CorrTest 与 ddBD。全部计算都经由本引擎的公开 Python API 执
行，结果可导出为 CSV、NEXUS、JSON 和 PNG 四种格式，并提供可复现的等效命令行面板。
界面完整支持中英双语并即时切换，内置浅色/深色两套主题、手绘矢量 SVG 图标，
以及一棵自带演示校正点的内置示例树，初次使用无需自备数据。

它**不属于**本仓库，也**不是**从本仓库安装的：Studio 是一款独立项目，把本引擎当作
普通第三方依赖来使用。引擎从 PyPI 安装；Studio 本身暂未上包索引，其 wheel 从它的
Release 安装：

```bash
pip install "openreltime[plot]" \
  "OpenRelTime-Studio @ https://github.com/ZengZichao/OpenRelTime-Studio/releases/download/v0.1.1/openreltime_studio-0.1.1-py3-none-any.whl"
openreltime-studio
```

Apple Silicon 上还可以直接取 Studio 的 Release 附件
`OpenRelTimeStudio-v0.1.1-macOS-arm64.zip`，完全不需要 Python。

主页：<https://github.com/ZengZichao/OpenRelTime-Studio> · 图形界面手册见该项目
自己的文档。只安装本引擎时，上文的 Python API 与全部命令行子命令都不受影响
——Studio 提供的是入口界面，不新增任何分析代码。

## 30 秒上手

```python
import openreltime as ort

tree = ort.read_tree("example.nwk", outgroup=["Out1", "Out2"])  # 置根并剔除外群
times = ort.rrf_rates_times(tree)                # 相对速率 + 相对节点时间
times.write("quickstart", with_rate=True, nexus=True)  # CSV 表 + NEXUS 时间树
```

等价命令行：

```bash
openreltime rates-times -i example.nwk --outgroup "Out1,Out2" -o quickstart
```

从比对出发的 BLB 管线（需要 IQ-TREE）：

```bash
openreltime blb -a alignment.fasta --iqtree iqtree --seed 42 -o blb_out
```

## 功能总览

- **RRF 核心**：提供两套平均约定——几何平均（msy044 式 28～42，R3F 与本项目
  默认）与算术平均（msy044 式 1～27，PNAS 2012 的原始语义）。同时支持
  极端速率保护、多分叉处理，以及 Newick 与 NEXUS（含 `[&rate=...]` 注释）读写。
- **校正转换**：支持硬边界的下限与上限，以及 `uniform`、`exponential`、
  `normal`、`lognormal` 四类概率密度，经全局时间因子 *f* 转换为绝对时间。
  同时支持 msz236 的有效边界（effective bounds）重采样法。
- **解析置信区间**：实现 msz236 的 delta 方法（微分法，式 7～14）。抽样方差
  按显式提供的枝长方差 → Poisson 近似 → 零 的顺序取第一个已提供的来源，结果按硬
  边界截断。有效边界法的重采样复制还给出经验的节点年龄 CI。
- **CorrTest**：基于已发表的固定系数逻辑回归模型的速率自相关检验（Tao et al. 2019）。
- **ddBD**：为贝叶斯定年估计出生-死亡物种分化树先验（出生率、死亡率、
  抽样比例）。
- **BLB 管线**：对位点执行两级重采样（无放回子采样，再用 bootstrap 上采样至原始
  长度），经 IQ-TREE 推断树后按各节点所代表单系群（clade）聚合时间。
- **MEGA-CC 桥接**（可选）：自动生成 `.mao`、调用 `megacc`、解析
  `megacc` 的 RelTime 输出并回写为 OpenRelTime 表格，便于 MEGA 用户迁移与第三
  方对照。

## 命令行一览

```
openreltime rates        -i tree.nwk [--outgroup og.txt] [--mean geometric] [-o prefix]
openreltime times        -i tree.nwk [--normalize] [--no-guard] [--r3f-compat] [-o prefix]
openreltime rates-times  -i tree.nwk [--plot timetree.png] [--r3f-compat] [-o prefix]
openreltime calibrate    -i tree.nwk -c calibrations.tsv [--method effective] [-o prefix]
openreltime ci           -c <calibrated prefix> [--n-sites 1000] [-o prefix]
openreltime corrtest     -i tree.nwk [--sister-resample 100] [-o prefix]
openreltime ddbd         -i tree.nwk [--anchor-time 1.85] [-o prefix]
openreltime tree2table   -i tree.nwk [--time] [-o prefix]
openreltime monophyly   -i tree.nwk [--taxon-table taxa.tsv] [--groups A,B] [-o prefix]
openreltime blb          -a aln.fasta --iqtree iqtree [--gamma 0.7 ...] [-o dir]
openreltime megacc       -i tree.nwk -c calibrations.tsv --megacc /path/to/megacc
```

参数会改变数值结果的分析步骤有七个：`rates`、`times`、`rates-times`、
`calibrate`、`ci`、`corrtest` 和 `ddbd`。它们都把运行参数、警告与随机种子合并
写入 `<prefix>_report.json`。共用同一前缀的各次运行
在该文件中各占一个顶层键。其余四个子命令（`tree2table`、`monophyly`、`blb`、
`megacc`）只输出自己的结果表。

## 校正文件格式（`calibrations.tsv`，制表符分隔）

```
node_id	taxon_set	min_bound	max_bound	density	density_params
	Homo_sapiens|Pan_troglodytes	6.0	8.5	.
	Elephas_maximus|Loxodonta_africana	.	.	exponential	offset=60;mean=20
	Dasypus_novemcinctus|Choloepus_didactylus	70	105	.
```

节点既可以用分类群集合的 MRCA 指定，也可以用内部节点 id 指定。校正可以只给
硬边界、只给密度，或者两者同时给出。密度支持 `uniform`、`exponential`、
`normal` 和 `lognormal` 四种，参数写成 `key=value;...` 形式。

上面三行只演示**语法**。作为数值，它们是占位符，任何一棵树无论怎样缩放速率都
无法同时满足；`calibrate` 现在会在这里停下来并点名有问题的节点，而不是把算坏的
年龄返回给你。只给密度、不给硬边界的行，在默认的 `--method bounds` 下同样会被
拒绝。想拿到开箱即用的校正表，请用
`data/examples/example_calibrations_feasible.tsv`；这两条规则的完整讲解见
[`docs/tutorials/tutorial2-end-to-end-zh.md`](docs/tutorials/tutorial2-end-to-end-zh.md)。

## 验证

本项目的验收门槛内置于测试套件。每个门槛实际断言什么，以及哪些随包交付的模块
不在覆盖范围内，都记录在 [`docs/validation-zh.md`](docs/validation-zh.md)。一键复现：

```bash
python -m pytest tests/
```

- **G1**：3 与 4 分类群 RRF 手算金标准值（两套均值约定，误差 < 1e-12）；
- **G2a 与 G2b**：274 分类群哺乳动物基准树与 R3F 回归（自动验收门槛的判据为
  斜率 ≥ 0.999、中位相对偏差 ≤ 1e-6；实测中位为速率 6.8e-16、时间 8.7e-16）；
- **G2c**：NEXUS `[&rate=...]` 往返互操作；
- **G6 与 G7**：断言 CorrTest 与 R3F 已存储的 CorrScore（0.9996）一致；覆盖
  校正边界与冲突检测；
- 金标准（golden）文件、示例输入与生成金标准的 R 脚本随包发布（`data/golden/r3f/`、
  `data/examples/`（见该目录的 `README-zh.md`）、`data/PROVENANCE-zh.md`、
  `reproduce/golden_r3f.R`）。

## 文档

完整手册同时提供**中英两个版本**；每篇文档开头都有指向另一语言版本的链接。
索引见 [docs/README.md](docs/README.md)。

| 文档 | 内容 | English |
|---|---|---|
| [docs/usage-zh.md](docs/usage-zh.md) | 详细使用说明（全部子命令、输出文件、FAQ 与报错信息） | [docs/usage-en.md](docs/usage-en.md) |
| [docs/parameters-zh.md](docs/parameters-zh.md) | 参数手册（全部可选参数，含各子命令选项对照表） | [docs/parameters.md](docs/parameters.md) |
| [docs/methods-zh.md](docs/methods-zh.md) | 方法、假设与局限 | [docs/methods.md](docs/methods.md) |
| [docs/validation-zh.md](docs/validation-zh.md) | 验证门槛、覆盖缺口与已测平台 | [docs/validation.md](docs/validation.md) |
| [docs/tutorials/](docs/tutorials/) | 三篇逐步教程（`*_zh.md`） | [docs/tutorials/](docs/tutorials/) |
| [docs/adr/README-zh.md](docs/adr/README-zh.md) | 架构决策记录（ADR）索引 | [docs/adr/README.md](docs/adr/README.md) |
| [CHANGELOG-zh.md](CHANGELOG-zh.md) | 版本记录 | [CHANGELOG.md](CHANGELOG.md) |
| [THIRD-PARTY-NOTICES-zh.md](THIRD-PARTY-NOTICES-zh.md) | 参考材料与被改造的常数，逐来源说明 | [THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md) |
| [data/PROVENANCE-zh.md](data/PROVENANCE-zh.md) | 每一个随包第三方文件的来源、许可与引用 | [data/PROVENANCE.md](data/PROVENANCE.md) |

## 许可与复用

本项目采用 GPL-3.0-or-later（见 `LICENSE`）。对 R3F、MEGA 源码与 ape 的行为
参考，记录在 `THIRD-PARTY-NOTICES-zh.md`；随包第三方文件的逐文件来源说明见
`data/PROVENANCE-zh.md`。CorrTest 仓库没有提供许可证，本项目未复制它的任何代码，
只引用 Tao et al. (2019) 已发表的常数。

## 作者

曾子超（Zichao Zeng）· [zengzichao@sjtu.edu.cn](mailto:zengzichao@sjtu.edu.cn) ·
[ORCID 0000-0001-6553-970X](https://orcid.org/0000-0001-6553-970X)

首选引用方式记录在 [`CITATION.cff`](CITATION.cff) 中。跨版本引用本软件整体请引用
引擎的 Zenodo *concept* DOI
[10.5281/zenodo.23053145](https://doi.org/10.5281/zenodo.23053145)；
[10.5281/zenodo.23053146](https://doi.org/10.5281/zenodo.23053146) 是已存档的
**v0.1.0** 版本记录，0.1.1 对应的版本记录在切出 0.1.1 发布时生成。
