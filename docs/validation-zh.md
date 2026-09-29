# 验证门槛、覆盖缺口与已测平台

**语言：中文** · [English](validation.md)

本文件集中说明三件事：OpenRelTime 的自动化测试套件*究竟证明了什么*，*没有证明
什么*，以及 CI 究竟在哪些平台上跑过。此处不淡化任何一点，
每一条“只是声明、并未实证”的说法都会明确标出。

复现引擎测试套件的命令如下（需要 `dev` 可选依赖）：

```bash
python -m pytest            # 收集范围遵循 [tool.pytest.ini_options]：tests/
```

桌面图形界面不在本文范围内：它作为独立产品 `OpenRelTime-Studio` 发行，其验证
由该项目自己的测试套件与 CI 负责，不由本套件覆盖。

## 1. 门槛地图（G1～G9）——每个门槛断言了什么

表中 `where` 指明承载该断言的测试文件。**断言（Asserted）**指存在 `pytest` 断言。
**实测（measured）**指软件计算并打印了某个数值，但没有对它施加通过或失败的界限。

| 门槛 | 主题 | 断言的条件 | 位置 |
|---|---|---|---|
| G1 | 手算 RRF（msy044） | 4 分类群几何平均加算术平均的节点值与手算常数一致到 1e-12；节点 id 与 ape cladewise 契约 | `tests/test_g1_hand_computed.py` |
| G2a | R3F 速率回归（274-tip 树） | 斜率 ≥ 0.999 **且** 中位相对偏差 ≤ 1e-6 | `tests/test_regression_golden.py::test_g2a_rate_regression` |
| G2b | R3F 时间回归（274-tip 树） | 斜率 ≥ 0.999 **且** 所有 clade 时间按 clade 键匹配 | `tests/test_regression_golden.py::test_g2b_time_regression` |
| G2c | NEXUS `[&rate=...]` 往返 | 读**与**写两个方向都保留速率注释（R3F 黄金 NEXUS 保住它自己的速率） | `tests/test_calibration_ci_edges.py::test_nexus_roundtrip_with_rate_annotations`、`tests/test_regressions_treeio_peripherals.py::test_r3f_annotated_golden_nexus_keeps_its_rates` |
| G3 | CLI 与 API 等价 | API 与 CLI 共用一条代码路径（`report.write()`）；CLI 选项确实抵达引擎；CLI 报错是干净消息而非 traceback。**跨进程字节一致只在运行测试的那台机器上、固定种子下成立**——它*不是*跨平台声明 | `tests/test_regressions_engine.py`（`test_*_cli_*`、`test_report_json_merges_across_analyses`）、`tests/test_regressions_treeio_peripherals.py`（`test_cli_*`、`test_branch_var_help_matches_the_header_the_parser_requires`、`test_rates_help_and_report_admit_the_guard_is_ignored`） |
| G4 | 输入校验 | 负枝长、重复末端名、多分叉（拒绝，或用 `random` 二叉化）、内群少于 3 个末端，以及“内部节点退化后残留为末端”的守卫都会报错 | `tests/test_calibration_ci_edges.py`（`test_negative_*`、`test_duplicate_*`、`test_polytomy_*`、`test_two_tip_ingroup_rejected`）、`tests/test_regressions_treeio_peripherals.py::test_validation_rejects_an_internals_degenerated_into_a_tip` |
| G5 | 外群置根与分类学 | 置根加剔除与 `ape::root(resolve.root=TRUE)` 加 `drop.tip` 字节一致；完整簇检查（error 与 warn 两档）；群名外群与群名校正；分类信息表解析 | `tests/test_regression_golden.py::test_rooted_topology_matches_ape`、`tests/test_taxonomy_monophyly.py` |
| G6 | CorrTest | CorrScore 与 R3F**未舍入**的参考值 0.99959505940763882 一致（随包的 `example_sr0_corrtest.txt` 只保留 5 位，即 0.9996）；lag-2 与 lag-3 衰减特征非退化；锚点缩放；CorrTest 使用的每个速率都等于引擎自己的逐节点速率 | `tests/test_regression_golden.py::test_corrttest_matches_r3f_score`、`::test_corrttest_lag_decay_features_are_not_degenerate`、`tests/test_regressions_engine.py`（`test_corrtest_anchor_*`） |
| G7 | 校正加解析 CI | 边界精确命中；联合约束；冲突约束报错；有效边界与按 taxon_set 校正；零方差时的 CI 恒等、随抽样方差变宽、硬边界截断；**无边界校正行被拒绝；`RV(R)` 为正且可饱和；被钉住节点拿到校正自身的不确定度而非零宽；CI 求导跟随被校正的 `mean` 与 `method`** | `tests/test_calibration_ci_edges.py`（`test_bounds_*`、`test_two_calibrations_*`、`test_conflicting_*`、`test_calibration_by_taxon_set_and_effective`、`test_ci_*`）、`tests/test_regressions_calibration_ci.py`（`test_calibration_without_any_bound_is_rejected`、`test_rate_heterogeneity_*`、`test_pinned_nodes_*`、`test_arithmetic_point_estimate_*`、`test_rrf_and_implied_rates_*`、`test_clamp_warning_*`、`test_divergent_trial_ids_are_detected`） |
| G8 | ddBD | 拟合出的 (birth, death, sampling) 在优化器容差内与 R3F 的运行结果一致；抽样比例语义；网格顺序；L-BFGS-B 异常终止的处理 | `tests/test_regression_golden.py::test_ddbd_close_to_r3f`、`tests/test_regressions_engine.py::test_ddbd_*`、`tests/test_regressions_treeio_peripherals.py::test_ddbd_*` |
| G9 | 运行时间与规模 | **实测，未断言。** 逐节点计时与“1,000 tips < 1 s”这一数字来自单台 Apple M 系列笔记本；没有任何 CI 作业断言运行时间上限，测试套件也从不跑 5,000 个以上分类群的树 | —（无自动化门槛） |

**关于门槛强度的说明（不做隐瞒）**

* **G2a 与 G2b**：*门槛*断言的是中位相对偏差 ≤ 1e-6，从不检查最大值。在 274-tip
  基准树上实测的中位约为 6.8e-16（速率）与 8.7e-16（时间），最大值分别约为
  4.8e-15 与 3.6e-14；这些最大值是观测值，不是断言。
* **G1**：算术平均约定只在浅层（2 层）树上检验，此时几何与算术折叠恰好重合。
  更深层的算术行为未在此设手算门槛。
* **G3**：“字节一致”是同机、同种子的性质（CSV 区域设置、换行符、JSON 键序与
  zip 元数据都可能跨平台不同）。本文件不把它表述为跨平台保证。
* **G6 与 G8**：这些是针对随包基准树上 R3F 已存储输出的回归检查，不能证明在
  其他拓扑上也一致。

## 2. 随包交付模块及其验证状态

`openreltime/` 交付下列 Python 模块。状态列说明测试套件*实际证明了*什么，不转述
模块自己的宣称。

| 模块 | 状态 | 依据 |
|---|---|---|
| `rrf`、`times`、`rates` | 相对 R3F 与手算**已验证** | G1、G2a、G2b |
| `tree`、`treeio` | **已验证**（拓扑、置根、节点 id 契约、NEXUS） | G2c、G4、G5 |
| `calibrate` | **已验证**（边界、冲突、effective） | G7 |
| `ci` | **已验证**（结构：恒等、变宽、截断） | G7 |
| `corrtest` | 在基准树上相对 R3F 得分**已验证** | G6 |
| `ddbd` | 在基准树上相对 R3F 参数**已验证** | G8 |
| `report` | **已验证**（写出开关、JSON 合并） | G3 |
| `taxonomy` | **已验证**（单系性、群名解析、表解析） | G5（`test_taxonomy_monophyly.py`） |
| `table` | 相对 R3F `tree2table()` 黄金值**已验证** | `test_tree2table_matches_golden` |
| `_external`、`_constants`、`cli` | 基础设施；行为经 CLI 与 G4 测试冒烟检验，没有独立的参考门槛 | — |
| **`bootstrap`（BLB）** | **已交付但数值上未验证。** CI 从不安装 IQ-TREE；唯一的一个测试用桩程序替换可执行文件，只检查 argv 与 `-seed` 的传递，以及退化复制的守卫。重采样与聚合得到的*分布*没有同任何参考比较过 | 无参考门槛 |
| **`megacc`** | **已交付但未验证。** 桥接需要真实的 `megacc` 二进制，CI 既不安装也不运行；解析测试的输入是*合成*输出。MEGA-CC 的一致性依据已发表的 R3F/MEGA 结果与随包黄金文件确立，而不是依据头对头运行 | 无参考门槛 |
| **`viz`** | **已交付但未与基线比较。** 测试断言的是结构性质（对数色阶；ddBD 密度在所绘坐标轴上积分为 1）；没有黄金图像比对 | 无参考门槛 |

因此本手册其他处把 `bootstrap`、`megacc` 和 `viz` 描述为**已交付但本文未加验证**；
`taxonomy` 与 `table` 已在上表覆盖。

## 3. 平台与 Python 版本——实测与声明

依据 `.github/workflows/ci.yml`：

| 维度 | CI 实测 | 仅在别处声明 |
|---|---|---|
| 操作系统 | `ubuntu-latest`、`macos-latest`（`test` 作业矩阵） | **Windows 属于声明**：出现在 `README` 与 `pyproject` 分类器中，但**没有 CI 作业**——是声明，不是实证 |
| Python | Ubuntu 与 macOS 上的 3.10、3.11、3.12 和 3.13 | — |
| 桌面图形界面 | **不在此处检验**——图形界面是独立产品 `OpenRelTime-Studio` | 其无头与交互式覆盖见该项目自己的测试套件与 CI |
| 外部工具 | **不安装**——`blb`（IQ-TREE）与 `megacc`（MEGA-CC）的路径从不在 CI 中执行 | — |
| 运行规模 | 测试套件止于 274-tip 黄金树 | “1,000 tips < 1 s”来自一台笔记本，**不是** CI 断言（G9） |

上表就是 `.github/workflows/ci.yml` 声明并在每次 push 与 pull request 时执行的产
品矩阵；运行日志发布在仓库的 Actions 页面。矩阵之外的平台属于声明的支持目标，不是
实证过的支持目标。

## 4. 外部工具——引用与假定的版本

**IQ-TREE**：`openreltime/bootstrap.py` 为 BLB 管线调用 IQ-TREE，代码与 CI 都没有
固定它的版本。

* 引用：Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Wood MR, von Haeseler A,
  Lanfear R. *IQ-TREE 2: New Models and Efficient Methods for Phylogenetic
  Inference in the Genomic Era.* Mol Biol Evol 37:1534-1538 (2020),
  doi:10.1093/molbev/msaa156。
* 实际假定的版本：**CI 既不固定也不安装任何版本。** 所用的旗标（`-s`、`-pre`、
  `-nt AUTO`、可选 `-seed`、`-m`）与 `<prefix>.treefile` 输出属于 IQ-TREE **2.x**
  的命令行约定（IQ-TREE 3.x 予以保留）。开发与画图阶段用的是 IQ-TREE 2.x 与 3.x
  系列。请把具体小版本视为操作者的选择，并在 BLB 运行记录中一并保存
  `iqtree --version` 的输出。

**MEGA-CC 与 MEGA 的等价性边界（一句话）**：OpenRelTime 只实现解析 RRF 层。MEGA
的 RelTime 额外施加了基于似然比检验的速率约束精修，OpenRelTime 没有这一步，因此
二者预期高度相关（约 0.98～0.99）但不相等；MEGA-CC 的一致性依据已发表的 R3F/MEGA
结果与随包黄金文件确立，而不是依据新做的头对头运行。

## 5. 门槛未能捕获的已知行为局限

以下条目在此重述，以免这些门槛被过度解读：

* CLI 中 `--outgroup` 是可选的。软件照单接受未置根的输入，输入树的根位置因此静
  默决定所有年龄。R3F 要求提供外群。这一条不在门槛覆盖之内：不带外群运行
  `openreltime times -i tree.nwk`，软件会顺利完成且不给任何警告。调用方必须自行
  确认输入树已在预期的单系群上置根。

## 6. 已纳入回归门槛的校正与 CI 行为

下列每条行为都在 `tests/test_regressions_calibration_ci.py` 中拥有会在回归时失败
的断言，因此它们不属于覆盖缺口；面向用户的完整描述以 `docs/parameters-zh.md` 与
`docs/methods-zh.md` 为准。

| 已断言的行为 | 测试 |
|---|---|
| 一行既不给 `min_bound` 也不给 `max_bound` 的校正在解析校正表时即按行号拒绝；非有限的可行因子会报错 | `test_calibration_without_any_bound_is_rejected`、`test_parse_calibrations_names_the_offending_line`、`test_feasible_factor_raises_instead_of_returning_inf` |
| `RV(R) = max(Vobs(R) - SV(R)/N, 0)` 按式 (9) 的逐谱系尺度计算，因此 msz236 的异质性分量在少位点区间为正，并随位点数增长饱和到 `Vobs(R)`；当夹取仍然生效时发出警告 | `test_rate_heterogeneity_*`、`test_clamp_warning_*` |
| 钉在硬边界上的节点改按校正自身断定的不确定度报告，而不是宽度为 0 的区间，并标记 `se_reliable = False`，`notes` 列说明原因 | `test_pinned_nodes_*`、`test_midpoint_fixed_calibration_*`、`test_unsatisfied_constraint_*` |
| 求导管线复现被校正记录的 `mean`、保护阈值，以及 `method="effective"` 时的有效边界 | `test_arithmetic_point_estimate_*`、`test_effective_result_is_differentiated_along_its_effective_bounds`、`test_divergent_trial_ids_are_detected` |
| 校正后的谱系速率与 RRF 速率分别命名，方差分解改用与 RRF 一致的定义 | `test_rrf_and_implied_rates_are_named_apart`、`test_variance_uses_the_rrf_rates_*` |

面向用户的相应说明见 `docs/parameters-zh.md`（*校正*、*置信区间*两节）与
`docs/methods-zh.md`（§4）。

## 7. CorrTest 参考值的精度

`data/golden/r3f/example_sr0_corrtest.txt` 只存了 5 位有效数字，因为 R3F 打印的是
`format(score, digits = 5)`。未舍入的参考 double 存于
`data/golden/r3f/example_corrtest_full_precision.txt`（如何在不触碰数值路径的前提
下捕获它，见 `reproduce/golden_r3f_fullprec.R`）。门槛 G6 与该未舍入值断言相等，
也就是逐位一致。G6 还单独断言：CorrTest 使用的每个速率都等于引擎自己的逐节点速
率。黄金目录本身在记录的 R 环境中重新生成，12 个文件中有 11 个逐字节复现，NEXUS
时间树的差异仅来自 ape 嵌入的写出时间戳。
