# 版本记录

**语言：中文** · [English](CHANGELOG.md)

OpenRelTime 遵循语义化版本规范。本文件记录各版本的发布内容，最新版本在最前。

## 0.1.1 — 2026-10-08

投稿前审阅驱动的缺陷修复与加固版本。下列每一项都是先把稿件陈述逐条对照代码、
再复现出问题后修正的。其中前三项原本会**静默给出错误结果**，因此按缺陷修复处理，
而非按改进处理。

### 修复

- **不可满足的嵌套校正现在会报错。** 要求后代节点比其祖先更老的界值，此前只追加
  一条警告并返回一个**违反所请求 `min_bound`** 的解，同时还放弃修复其余校正行。
  现在抛出 `ValueError` 并点名节点与上下界，与既有的「根无法变轻」路径一致。
- **`corrtest(anchor_node=...)` 不再产出 `inf`/`NaN`。** 其 `anchor_time` 默认
  为 `0.0`，锚点缩放因此除以零，并把 `{"mean.rate": Infinity, "sd.rate": NaN}`
  写进 JSON 报告——这不是合法的严格 JSON，符合规范的解析器会直接拒绝。现在给定
  锚点节点而 `anchor_time` 非正时前置拒绝，且报告以 `allow_nan=False` 序列化。
- **未知的密度参数键现在报错而非被忽略。** `density='normal'` 配
  `density_params='mu=5;sigma=1'` 此前被接受并按默认的 `N(0, 1)` 抽样，静默套用
  了错误的先验。现在每个密度族都有键白名单（`uniform`: min,max；
  `exponential`: offset,mean；`normal`: mean,sd；`lognormal`:
  offset,meanlog,sdlog），未知键或缺键一律报错。
- **退化输入不再被判定为显著。** 未定义的 Spearman 相关系数此前被替换为 `0.0`
  再送入固定系数的 logistic 模型，于是全零枝长的树会得到 `CorrScore ≈ 1` 且落在
  `P < 0.001` 档——用零信息给出最强结论。此类树现在被拒绝。
- **拒绝 `n_effective=1`**，不再「跑 2 次却报告 1 次」。
- **`openreltime.<模块>` 名字遮蔽。** 包内重导出的函数与其同名模块互相遮蔽，
  导致 `openreltime.calibrate.Calibration` 抛 `AttributeError`。现在这些属性是
  可调用的模块对象：`openreltime.calibrate(...)` 与
  `openreltime.calibrate.Calibration` 都可用。
- **`--seed` 现在能抵达随机多分叉解析。** `read_tree()` 一直接受 `seed`，
  但命令行的读树辅助函数从不传它，于是除 `blb` 外每个子命令的
  `--resolve random` 都不可复现、也不进报告。现在凡是有 `--resolve` 的命令
  都同时有 `--seed` 并向下传递；种子会记进报告（只写 CSV、没有报告块的
  `monophyly` 是文档化的例外），`ci` 则用它所依据的那次标定记录下的种子重读树。
- `ddbd` 与 `corrtest` 中归一化极大值处的除零保护。

### 新增

- `blb` 与 `megacc` 现在写出 JSON 报告，并记录所调用外部程序的版本横幅，
  使运行结果可归因到具体的 IQ-TREE / MEGA-CC 构建。仍然不锁定任何版本。
- 面向 FigTree 可读性的 NEXUS 结构测试（有根标记、`TAXA` 块无注释、每个节点
  一条 `[&rate=...]`、枝长为有限数值、可被独立重解析），并显式测试普通 `.nwk`
  输出**不**携带速率注释。
- `data/examples/example_calibrations_feasible.tsv`：一份真正可运行、彼此自洽的
  示例。原有 `example_calibrations.tsv` 只用于演示语法——它的纯密度行使其在默认
  `method='bounds'` 下直接抛错，而其互相矛盾的界值会把 4 条谱系压到非正时长。
- `tools/check_gates.py`：把 `docs/validation.md` 的门禁表与实测收集到的测试 id
  及 `@pytest.mark.gate` 标记互相校验；`tools/count_tests.py` 报告各套件计数。
- 单元测试套件范围内的 `@pytest.mark.gate(...)` 标记，并已在 `pyproject.toml` 注册。
- `figures/`：手稿四张图的绘图脚本、S2/S3 所读取的两份 `-data.json`，以及
  `figures/README.md`。这些文件与产出已投稿图面的源文件逐字节一致——为每个脚本
  加上 `--outdir` 重跑，四张已发布 PNG 均可逐字节复现。`MANIFEST.in` 随源码分发
  一并打包它们，使 Data availability 中点名 `figures/` 的表述对存档本身成立，
  而不仅对作者的工作副本成立。渲染产物不进入分发包。

### 变更

- `docs/validation.md`、`docs/validation-zh.md` 与 `README.md` 原称 CI 在 Ubuntu
  与 macOS 上都跑 Python 3.10–3.13。`ci.yml` 实际排除了 macOS 的
  3.10/3.11/3.12，文档已改为真实矩阵：Linux 3.10–3.13，macOS 仅 3.13。
- `ZERO_BRLEN_WARN_FRACTION` 与 `MIN_INGROUP_TIPS` 现在真正被使用它们的代码引用，
  不再是与字面量并存的第二份真相。
- `ci.py` 模块文档串不再声称区间使用「标定后的绝对速率」；实现一直使用的是尺度
  不变的相对 RRF 速率。
- 门禁 id 现在可机器校验，且 `tests/test_taxonomy_monophyly.py` 的标注与文档
  分配给它的门禁一致。
- 文档已与代码和测试对齐：G5 一行不再声称与 `ape` 置根「字节一致」（断言比的是
  规范化后的 `(末端, 深度)` 二元组）；G2a/G2b 的说明拆成两条，因为 G2b 确实约束了
  逐节点最大值；FigTree 可读性契约改引用真实存在的测试 id；`seed` 一节点名多分叉
  解析用的是标准库 `random.Random`，不再承诺「一律经 `numpy.random.Generator`」；
  `data/golden/r3f/example_RRF_rates.nwk` 也补上了 `PROVENANCE.md` 开头对每个文件
  都作出承诺的那一行。

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
