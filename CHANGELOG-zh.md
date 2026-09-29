# 版本记录

**语言：中文** · [English](CHANGELOG.md)

OpenRelTime 遵循语义化版本规范。本文件记录各版本的发布内容，最新版本在最前。

## 0.1.0 — 2026-09-29

OpenRelTime 引擎的首个版本：按 Tamura、Tao 与 Kumar（2018, *MBE* 35:1770-1782）
提出的**相对速率框架（RRF）**实现的 GPL-3 Python 程序，也就是 MEGA 中 RelTime
的解析核心。软件从带枝长的系统树出发估计谱系相对速率与节点时间，把校正点转换为
绝对时间并给出解析置信区间，同时提供 CorrTest、ddBD 与轻量自助重采样（bag of
little bootstraps）管线。全部分析既可在 Python 中调用，也可在命令行中运行。

运行环境为 Python ≥ 3.10，依赖 numpy、scipy、pandas 与 click；不含编译扩展。

### 分析引擎

- **RRF 速率与时间估计**（`openreltime.rrf`、`rates`、`times`；
  `rrf_rates`、`rrf_times`、`rrf_rates_times`）：后序遍历中以闭式求解局部的 3
  与 4 谱系问题，再以前序遍历做一次调整；子树折叠为组合枝长；内群根锚定在平均速
  率 1。
- **两套平均约定**：几何平均（msy044 式 28～42，即 R3F 与 RelTime 使用的约定，
  默认）与算术平均（msy044 式 1～27，即 2012 年论文的语义）。折叠算子随所选平均
  约定改变。
- **极端速率保护**：调整后速率越过阈值带的节点改取最近的未越界祖先的时间（R3F
  语义）。默认阈值 `20.0`，用 `--no-guard` 关闭，受影响节点数以
  `n_rate_guarded_nodes` 记录。
- **与 R3F 共享的数值处理**：在 R3F 精确相同的位置代入 epsilon，零枝长占比超过
  10% 时发出警告。全部默认值登记在 `openreltime/_constants.py`。

### 校正与置信区间

- **校正为绝对时间**（`openreltime.calibrate`）：`bounds` 方法（全局时间因子 *f*
  取可行区间交集的中点；对有责任的谱系按比例缩放速率并向后代传播，上限 100
  轮；冲突约束集连同逐节点诊断一并报告）与 `effective` 方法（msz236 的密度重采
  样，默认 10,000 次复制）。
- **概率密度**：`uniform`、`exponential`、`normal`、`lognormal`，按 `key=value`
  严格解析（不做表达式求值），可与硬边界 `min_bound`/`max_bound` 一同写在制表符分
  隔的 `calibrations.tsv` 中。
- **节点定位**：可用内部节点 id，也可用某一分类群集合的 MRCA；集合中的成员既可
  以是末端名，也可以是类群名。
- **可取消的校正**：`calibrate(..., progress=...)` 的回调返回 `False` 时抛出
  `CalibrationCancelled`。
- **解析置信区间**（`openreltime.ci`）：msz236 的 delta method（式 7～14）。枝方
  差拆分为抽样分量与速率异质性分量；偏导数由**整条估计管线**的中心差分求得；区间
  在硬边界处截断并在 0 处收口。抽样方差可来自逐枝的 `var` 表、Poisson 近似
  `vS(b) = b/L`，或者完全不给定。
- **一套估计器贯穿始终**：求导过程复现产生点估计时所采用的平均约定、保护阈值与校
  正边界（包括 `effective` 结果的有效边界）。
- **可靠性标注**：年龄恰为自己约束不动点的被校正节点，改按其校正本身断定的不确定
  度报告，并标记 `se_reliable = False`，`notes` 列说明原因。有效边界法另给出经验
  2.5 与 97.5 节点年龄分位数。

### 速率检验、树先验与重采样

- **CorrTest**（`openreltime.corrtest`）：Tao et al. (2019) 的速率自相关统计量，
  采用已发表的固定系数逻辑回归模型，输出 P 值档位，小树可选 sister 对重抽样，并可
  选绝对速率锚定。
- **ddBD**（`openreltime.ddbd`）：估计出生-死亡物种分化树先验（出生率、死亡率、抽
  样比例），以 SSE/KL 打分的初值网格配合 L-BFGS-B 拟合，可选 anchor 节点与锚定时
  间，并可绘制拟合密度。
- **轻量自助重采样**（`openreltime.bootstrap`）：对比对做两级位点重采样，复制树由
  外部 IQ-TREE 可执行文件推断，节点年龄按 clade 聚合为中位数与 2.5/97.5 分位数并
  附上包含该 clade 的复制比例；运行种子会传给每一次 IQ-TREE 调用。

### 树的读写与分类学工具

- **自建 Newick 与 NEXUS 读写**（`openreltime.treeio`、`openreltime.tree`）：
  `[&rate=...]` 节点注释读入时保留、写出时重新生成（R3F 与 FigTree 方言）；注释与
  带引号的标签不会破坏名称。
- **外群置根与剔除**（`root_outgroup`）复现 `ape::root(resolve.root = TRUE)` 加
  `ape::drop.tip` 的可观测行为，包括一元节点收缩与根边丢弃；多末端外群需通过完整
  簇检查，内群至少保留 3 个末端。
- **ape cladewise 节点编号**（`PhyloNode.assign_ids`：末端 1..n，内部节点按前序，
  根 = n+1），因此节点表与 `--anchor-node` 取值可与 R3F、ape、MEGA 的输出直接对照。
- **多分叉处理**：默认报错，也可用零枝长随机二叉化。
- **单系性与分类表工具**（`openreltime.taxonomy`）：`check_monophyly`、
  `parse_taxon_table`（TSV/CSV 的末端到类群对照表）与 `detect_groups_from_labels`
  （按 `Genus_species` 惯例自动识别属名）。类群名同样可用作外群与校正的
  `taxon_set` 取值。
- **节点表**（`openreltime.table.tree2table`）：按 R3F 的版式输出枝长表或节点年龄
  表。

### 输出与可复现性

- **结果对象**（`openreltime.report`）提供 `to_pandas()`、`to_json()` 与
  `write(prefix)`。CSV、Newick、NEXUS 与 JSON 输出与 CLI 走同一条代码路径；
  `--r3f-compat` 切换为 R3F 的文件名与列名。
- **可复现性报告**：参数会影响数值的七个分析步骤（`rates`、`times`、
  `rates-times`、`calibrate`、`ci`、`corrtest`、`ddbd`）把自己的参数、警告与随机
  种子合并写入 `<prefix>_report.json`，共用同一前缀的每次运行各占一个顶层键；
  `ci` 会重放其 `calibrate` 运行记录下的设置。
- **全部随机性**（密度抽样、BLB、sister 重抽样）都经由同一个带种子的
  `numpy.random.Generator`。

### 使用入口

- **Python API**：`import openreltime as ort`，导出上述函数以及 `read_tree`、
  `parse_calibrations` 与各结果类型。
- **命令行**：`openreltime`（别名 `ort`），11 个子命令与 API 一一对应：`rates`、
  `times`、`rates-times`、`calibrate`、`ci`、`corrtest`、`ddbd`、`tree2table`、
  `monophyly`、`blb`、`megacc`；报错为可读消息而非 traceback。
- **MEGA-CC 桥接**（`openreltime.megacc`）：生成 `.mao` 分析描述文件、调用外部
  `megacc`、解析其输出并返回与 R3F 同构的节点时间表，其中 `source` 列逐行标注来
  源。
- **可选绘图**（`openreltime.viz`，`[plot]` extra）：对数色阶的速率着色时间树、节
  点年龄区间图，以及按所绘坐标轴归一化的 ddBD 密度。
- **外部进程不经 shell 启动**，转传的每个参数都先经严格白名单校验
  （`openreltime._external`）。

### 验证与数据

- **参考工具黄金标准**：R3F 在 274-tip 哺乳动物基准树上的输出
  （`data/golden/r3f/`）、生成它们的 R 脚本（`reproduce/golden_r3f.R`）、示例输入
  （`data/examples/`）以及逐文件来源记录（`data/PROVENANCE.md`）。
- **自动化门槛**（`tests/`）：G1 两套约定的 3 与 4 分类群手算值；G2a 与 G2b 274-tip
  树上的 R3F 速率与时间回归；G2c NEXUS `[&rate=...]` 往返；G3 CLI 与 API 等
  价；G4 输入校验；G5 与 ape 一致的置根及分类学；G6 相对 R3F 未舍入得分的
  CorrTest；G7 校正与 CI 行为；G8 ddBD 参数；G9 运行时间测量（记录而不做断言）。
  数值回归测试位于 `tests/test_regressions_engine.py`、
  `tests/test_regressions_calibration_ci.py`、
  `tests/test_regressions_rrf_engine.py` 与
  `tests/test_regressions_treeio_peripherals.py`。
- **端到端套件**（`e2e_tests/`）：用随包示例数据跑完整 CLI 与 API 流程，并与同一
  套黄金值比对；随源码包一起发布。

### 文档

- 手册全部提供**中英两个版本**：`README.md`、`docs/usage-en.md`、
  `docs/parameters.md`、`docs/methods.md`、`docs/validation.md`、
  `docs/tutorials/` 中的三篇教程、`docs/adr/` 的架构决策记录索引、`CHANGELOG.md`、
  `THIRD-PARTY-NOTICES.md` 与 `data/PROVENANCE.md`。

### 本版本的适用范围

- **只做解析 RRF。** MEGA 的 RelTime 在此之上还有一步基于似然比检验的速率约束精
  修，本引擎不实现该步骤，因此 OpenRelTime 与 MEGA 是高度相关（在已发表基准上约
  0.98～0.99）而非逐节点相同。与 MEGA 比较时请报告相关性。
- **与参考工具的比对范围限于随包基准树。** CorrTest 与 ddBD 的门槛是相对 R3F 在
  274-tip 树上已存储输出的回归检查，不能证明在其他拓扑上也一致；速率与时间回归的
  门槛断言的是斜率与*中位*相对偏差，最大值属于测量并报告的观测值，不作断言。
- **`bootstrap`、`megacc` 与 `viz` 随包交付但没有数值参考门槛。** BLB 与 MEGA-CC
  路径需要外部可执行文件，测试套件并不安装它们；这三个模块的测试覆盖参数传递、输出
  解析与结构性质。MEGA-CC 的一致性依据已发表的 R3F/MEGA 结果与随包黄金文件确立，而
  不是依据新做的头对头运行。
- **平台。** CI 在 Linux 与 macOS 上跨 Python 3.10～3.13 运行测试套件；Windows 是
  声明的支持目标，CI 不在该平台运行。所报告的计时来自单台机器，测试套件止于 274-tip
  黄金树。
- **外群可选。** `--outgroup` 可以省略；此时树按读入原样使用，其根位置决定所有年龄
  且不发出警告，因此调用方须确认输入树已在预期的单系群上置根（R3F 强制要求外群）。
