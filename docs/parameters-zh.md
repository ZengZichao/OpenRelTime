# 参数手册

**语言：中文** · [English](parameters.md)

OpenRelTime 的每一个可选参数：语义、默认值、取值范围与使用建议。

图例——**API**：关键字参数；**CLI**：命令行选项；**where**：所在模块。

---

## 各子命令支持哪些选项

本表等价于 `openreltime <子命令> --help` 的结果，同时充当下文各条目的快速索引。
除 `blb` 使用 `-o/--output-dir` 之外，所有子命令都有 `-i/--input`、
`-o/--output` 与 `-h/--help`，表中不再重复。

| 选项 | rates | times | rates-times | calibrate | ci | corrtest | ddbd | tree2table | monophyly | blb | megacc |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `--outgroup` | ✓ | ✓ | ✓ | ✓ | –¹ | ✓ | ✓ | – | – | ✓ | ✓ |
| `--fmt` | ✓ | ✓ | ✓ | ✓ | –¹ | ✓ | ✓ | ✓ | ✓ | – | – |
| `--resolve` | ✓ | ✓ | ✓ | ✓ | –¹ | ✓ | ✓ | – | ✓ | – | – |
| `--taxon-table` | ✓ | ✓ | ✓ | ✓ | –¹ | ✓ | ✓ | – | ✓ | – | – |
| `--outgroup-check` | ✓ | ✓ | ✓ | ✓ | –¹ | ✓ | ✓ | – | – | – | – |
| `--mean` | ✓ | ✓ | ✓ | ✓ | ✓² | – | – | – | – | – | – |
| `--rate-ratio-threshold` | ✓³ | ✓ | ✓ | ✓ | ✓² | – | – | – | – | – | – |
| `--no-guard` | – | ✓ | – | ✓ | ✓² | – | – | – | – | – | – |
| `--normalize` | – | ✓ | ✓ | – | – | – | – | – | – | – | – |
| `--r3f-compat` | – | ✓ | ✓ | – | – | – | – | – | – | – | – |
| `--plot` | – | – | ✓ | – | – | – | ✓ | – | – | – | – |
| `--calibrations` | – | – | – | ✓ | – | – | – | – | – | – | ✓⁴ |
| `--method` | – | – | – | ✓ | – | – | – | – | – | – | – |
| `--n-effective` | – | – | – | ✓ | – | – | – | – | – | – | – |
| `--seed` | – | – | – | ✓ | – | ✓ | – | – | – | ✓ | – |
| `--branch-var` | – | – | – | – | ✓ | – | – | – | – | – | – |
| `--n-sites` | – | – | – | – | ✓ | – | – | – | – | – | – |
| `--level` | – | – | – | – | ✓ | – | – | – | – | – | – |
| `--sister-resample` | – | – | – | – | – | ✓ | – | – | – | – | – |
| `--anchor-node`、`--anchor-time` | – | – | – | – | – | ✓ | ✓ | – | – | – | – |
| `--sampling-frac` | – | – | – | – | – | – | ✓ | – | – | – | – |
| `--measure` | – | – | – | – | – | – | ✓ | – | – | – | – |
| `--time` | – | – | – | – | – | – | – | ✓ | – | – | – |
| `--groups`、`--tips`、`--group-sep` | – | – | – | – | – | – | – | – | ✓ | – | – |
| `--alignment`、`--iqtree`、`--gamma`、`--n-subsample`、`--n-replicate`、`--iqtree-arg` | – | – | – | – | – | – | – | – | – | ✓ | – |
| `--megacc`、`--mao` | – | – | – | – | – | – | – | – | – | – | ✓ |
| `--workdir` | – | – | – | – | – | – | – | – | – | ✓ | ✓ |

¹ `ci` 不在命令行接受树输入，而是从 `calibrate` 那次运行留下的
`<prefix>_report.json` 中读回树、外群、输入格式、多分叉处理与外群检查方式（见
`calibrated`）。
² 覆盖 `calibrate` 记录的值（见 `mean`、`rate_ratio_threshold`）；不指定即“沿用
记录值”。
³ 速率计算过程会**忽略**这个选项——见 `rate_ratio_threshold`。
⁴ 对 `megacc` 是可选项：不给时，桥接只按枝长运行 RelTime。

---

## 输入处理

### `fmt` — 输入格式
* 取值：`"newick"` 或 `"nexus"`。
* 默认：`"newick"`（`openreltime.read_tree(path, fmt="newick")`；CLI 的
  `--fmt` 同样默认 `newick`，且只接受这两个值）。
* NEXUS 文件可携带 `[&rate=...]` 注释（R3F 与 FigTree 的方言）；读入时保留这些
  注释，写出时重新生成。

### `outgroup` — 置根外群
* 类型：末端名列表（API），逗号分隔字符串或每行一个名称的文本文件
  （CLI：`--outgroup`）。
* 默认：`None`——不置根也不剔除；树按读入原样使用，且必须已经是有根二叉树。
  CLI 选项同样可选，不传即 `None`。
* 软件在外群的 MRCA 处置根，随后删除外群末端，输出仅含内群末端。行为与
  `ape::root(resolve.root=TRUE)` 加 `ape::drop.tip` 一致（ADR-004）。
* 建议：使用 1～5 个演化较慢、构成单系群的分类群。
* 两种写法：**单个末端**，或**多个末端**且必须构成完整簇（见
  `outgroup_check`）。不是末端名的记号会按**类群名称**解析（经 `taxon_table`
  或末端标签前缀自动识别），并展开为该名称的全部成员。
* 省略该选项不会触发任何警告，而所有 RelTime 量都锚定在根上：输入文件的根位置
  直接决定全部节点年龄。请务必自行确认拓扑（`docs/validation-zh.md` §5）。

### `taxon_table` — 末端分类信息
* 类型：`None` 或 `{类群: [末端, ...]}` 映射（API：
  `openreltime.parse_taxon_table(path)`；CLI：`--taxon-table PATH`）。
* TSV 或 CSV，含一列末端名（`tip`、`label`、`name` 等）与一列类群名（`group`、
  `clade`、`taxa`、`classification` 等）；两个表头都识别不出时，第一列为末端、
  第二列为类群。不给表时，软件按末端标签前缀自动识别类群（`Genus_species` →
  `Genus`）。
* 供群名外群、群名校正 `taxon_set` 取值以及 `monophyly` 子命令使用。

### `outgroup_check` — 多末端外群的完整性检查
* 取值：`"error"`（默认）或 `"warn"`。CLI：`--outgroup-check`。
* 多末端外群必须是完整簇：外群 MRCA 之下的所有末端都要属于该外群。
  `"error"` 抛出 `ValueError` 并列出多余末端；`"warn"` 记录警告，并仍以 MRCA
  置根。单末端外群天然完整，不做检查。

### `monophyly` 子命令
* `openreltime monophyly -i tree.nwk [--taxon-table taxa.tsv] [--groups A,B]
  [--tips A,B,C] [--group-sep REGEX] [--fmt ...] [--resolve ...] -o PREFIX`。
* 判断每个类群是否单系——即该类群末端的 MRCA 之下不再包含其他末端（与
  `outgroup_check` 同一判据）。
* `--groups` 指定要检验的类群（默认检验分类信息中的全部类群）；`--tips` 改为
  检验一个临时给定的末端集合。两者可各自独立地收窄检验范围。
* `--group-sep REGEX`——自动识别标签前缀时使用的分隔模式，默认取 `_`、`|`、
  `@`、空白中的第一个。仅在未提供 `--taxon-table` 时有意义。
* 写出 `<prefix>_monophyly.csv`（列 `group`、`n_tips`、`is_monophyletic`、
  `n_mrca_tips`、`extra_tips`、`missing_names`、`tips`）；对非单系类群，另外在
  stderr 记录警告。`missing_names` 列出请求了但树中不存在的名称。

### `resolve_polytomy` — 非二叉树
* 取值：`"error"`（默认）或 `"random"`。CLI：`--resolve error|random`。
* `"random"` 以零枝长把多分叉随机二叉化（记录警告）。零枝长会让受影响节点的
  RRF 速率与时间不可靠（R3F 在零枝占比超过 10% 时警告；OpenRelTime 采用同一
  阈值）。

## RRF 估计

### `mean` — 平均约定
* 取值：`"geometric"`（默认）或 `"arithmetic"`。
* 几何平均 = msy044 式 28～42，与 R3F 兼容。算术平均 = PNAS 2012 原始语义
  （式 1～27），保留用于交叉验证；两者输出会有差异。

### `rate_ratio_threshold` — 极端速率保护
* 类型：float 或 `None`。默认：`rrf_rates` 为 `None`，`rrf_times` 与
  `rrf_rates_times` 为 `20.0`。CLI：`--rate-ratio-threshold` 或 `--no-guard`
  （关于 `rates` 的特别说明见下）。
* 前序调整之后，若某节点的调整后速率超过阈值（或低于阈值的倒数），软件会把该
  节点的时间替换为最近的未超限祖先节点的时间（R3F 语义）。`None` 关闭保护。
* **保护只属于时间计算过程。** `rrf_rates` 保留该关键字只为与 `rrf_times` 签名
  对称，从不应用这个关键字：R3F 参考实现只在重写节点时间的那一步设保护
  （`rrf_times.R` 353～373 行），而 `rrf_rates.R` 既无保护也无时间可重写。因此
  给 `rrf_rates` 传数值会记录一条警告，该警告同时进入 `RateResult.warnings` 与
  `<prefix>_report.json`。`openreltime rates --rate-ratio-threshold` 走的正是这
  条代码路径。
* 敏感性分析：分别用 `None`、20、10 重复分析，可量化受影响节点数（报告 JSON 中
  的 `n_rate_guarded_nodes`）。
* `openreltime ci` 会重跑区间背后的估计器，因此它也暴露 `--mean`、
  `--rate-ratio-threshold` 与 `--no-guard` 来覆盖 `calibrate` 记录的值；不指定即
  完全复现记录设置。`--no-guard` 与 `--rate-ratio-threshold` 互斥。

### `normalize` — 节点时间重缩放
* 默认：`False`。CLI：`--normalize`。
* 把节点时间除以其最大值（根时间 = 1.0），即 MEGA 与 PNAS 2012 的报告惯例。
  R3F 报告的是调整后的原始时间（本项目默认亦如此）。

## 校正（`calibrate`）

### `method`
* `"bounds"`（默认）——直接使用 `min` 与 `max` 约束。全局时间因子 `f` 取可行区间
  `[min/t, max/t]` 交集的中点。残余的违反会触发按比例缩放肇事谱系的速率，并向
  后代传播，最多 100 轮（`openreltime._constants.CALIBRATION_MAX_ITER`）。
  不可行的约束集会抛出 `ValueError` 并给出逐节点诊断；若算出的可行因子是非有限
  值，同样报错，而不会把 `f = inf` 传播进每一个绝对年龄。只给密度、不给硬边界
  的行在 `method="bounds"` 下不构成任何约束，解析校正表时会拒绝该行（改用
  `--method effective`，或补上 `min_bound` 与 `max_bound`）。
* `"effective"`——msz236：对每个密度校正，每次复制抽取两个日期作为
  (min, max)；分析重复 `n_effective` 次；每个被校正节点的 2.5 与 97.5 分位数成
  为最终估计所使用的*有效边界*。同时报告各节点的经验 2.5 与 97.5 年龄分位数。

### `n_effective`
* 默认 `10_000`（msz236 协议）；大树的探索性运行可调低（≥ 200）。复制之间复用
  已缓存的组合枝长，因此每次复制的代价是 O(#校正数 + n)。

### `seed`
* 一切随机性（密度抽样、BLB、sister 重抽样）都经由
  `numpy.random.Generator`；请把种子（写入 `<prefix>_report.json`）一并记录，
  以便复现。

### 校正文件格式（`calibrations.tsv`，TSV）
列：`node_id`（可选）、`taxon_set`（`|` 分隔，可选；两者恰好给一个）、
`min_bound`、`max_bound`（`.` 表示缺失）、
`density`（可选：`uniform|exponential|normal|lognormal`）、
`density_params`（以 `;` 分隔的 `key=value` 对）：
* `uniform`：`min,max`
* `exponential`：`offset,mean`
* `normal`：`mean,sd`
* `lognormal`：`offset,meanlog,sdlog`

`taxon_set` 的条目可以是末端标签，也可以是单个**类群名**（经 `taxon_table` 或末
端标签前缀自动识别解析）。一旦某个 `taxon_set` 在内群树上判定为非单系，软件就记
录一条警告，并用该 `taxon_set` 末端的 MRCA 作为校正节点。

参数一律按严格的键值对解析，绝不做表达式求值。

**每一行都必须构成约束。** 解析阶段会拒绝既无 `min_bound` 又无 `max_bound`（也没
有密度）的行，并报出出错的行号——否则这一行会把全局因子 `f` 推向 `inf`，让每一
个绝对年龄都变成 `inf`。只给密度的行在 `method="bounds"` 下同样会遭到拒绝；请改
用 `method="effective"`，或补上硬边界。缺列、无法解析的数值、重复的节点 id 也都
会报错。

## 置信区间（`confidence_interval`）

### `calibrated`、CLI `-c` 与 `--calibrated`
* 输入是 `CalibratedResult`（API），或**上一次 `openreltime calibrate` 运行的
  输出前缀**（CLI）。CLI 会从 `<prefix>_report.json` 读回树、外群、输入格式、多
  分叉处理与外群检查方式，使区间与它所定性的那次估计保持一致；找不到报告或原始
  树文件时，CLI 给出可读的错误。
* 一套估计器贯穿始终：求导管线复现点估计**实际使用**的平均约定（`mean`）、速率
  比保护以及校正边界。`method="effective"` 的结果沿其有效边界求导（msz236 正是
  把有效边界当作约束），`mean="arithmetic"` 的结果沿算术平均的导数求导。`method`
  取 `bounds` 与 `effective` 之外的值会报错；没有有效边界的 effective 结果（例如
  `n_effective < 2`）同样报错。
* 输出 `<prefix>_ci.csv` 每个内部节点一行：`node_id`、`label`、`time`、`se`、
  `lower`、`upper`、`width`、`se_reliable`、`notes`。凡区间**不是**来自 msz236
  eq (7) 的枝方差分解的行，`se_reliable` 为 `False`，并由 `notes` 说明原因（校正
  边界钉住、求解器未满足约束、下界截到 0）。这类节点的数量在报告 JSON 中记为
  `n_nodes_se_not_reliable`。

### `branch_var`、`n_sites`、`seq_length` — 抽样方差 vS(b)
* 默认：三者皆为 `None`（`ci.py:confidence_interval`）。`None` 不是一个方差取值，
  而是表示“未提供该来源”。什么都不提供时，回退产生的才是零方差，它不是默认值。
* 行为——按优先级取第一个被提供的来源：
  1. 显式 `branch_var` 字典 `{子节点 id: 方差}`；
  2. Poisson 近似 `vS(b) = b/L`，`L = seq_length`，否则 `L = n_sites`（先检查
     `seq_length`，同时给定时以 `seq_length` 为准）；
  3. 三者都为 `None`（默认调用）时 `vS(b) = 0`，区间因此只含速率异质性方差
     （msz236 的模拟协议）。
* CLI：`--branch-var TSV`（`node_id`、`var` 两列）与 `--n-sites L`；
  `seq_length` 没有对应的 CLI 选项（仅 API）。
* 实际使用的来源记录在 `<prefix>_report.json` 的
  `confidence_interval.parameters.v_s_source`。Poisson 近似按其原文属于
  近似方法（msz236 “sampling variance” 一节）。

### `level`
* 默认 `0.95`。区间在硬校正边界处截断；若对称区间的下界小于 0，则截到 0（分化
  时间不可能为负；`notes` 会记录这次截断）。

### 速率异质性分量 `RV(R)`（不由用户设定）
* `RV(R) = max(Vobs(R) - SV(R)/N, 0)`，取在 msz236 eq (9) 的**逐谱系**尺度上：
  `SV(R)` 是跨 `N` 个谱系的求和，因此要先除以 `N` 再与逐谱系平均值 `Vobs(R)`
  相减。直接拿原始和去减一个平均值，会让软件在任何真实树上都把 `RV(R)` 夹到 0，
  整个异质性分量随之消失。
* `Vobs(R)`、`SV(R)`、`RV(R)` 在报告 JSON 中记为 `Vobs_R`、`SV_R`、`RV_R`；
  method 字符串还会回显所用的抽样方差来源、被复现的 `method` 与 `mean`。
* 当 `SV(R)/N > Vobs(R)`（抽样主导的树：位点数相对于观测速率离散太少）时
  `RV(R)` 夹到 0，`v(b) = vS(b)`，并发出警告说明——此时区间只反映抽样噪声。
* 年龄恰好等于该节点自身约束不动点的节点（硬边界上，或可行区间中点）通过
  eq (7) 的导数恒为 0。这类节点不再报告“零不确定度”，而是给出该节点校正本身断
  定的不确定度（有效边界、密度分位数，或用户给定的 min 与 max），并标记
  `se_reliable = False`。

## CorrTest（`corrtest`）

### `sister_resample`
* 默认 `0`。末端数少于约 50 的树应使用 ≥ 50 次 sister 对分配重抽样；取抽样间
  Spearman rho 的均值（R3F 行为）。

### `anchor_node`、`anchor_time`
* 默认：`anchor_node=None`、`anchor_time=0.0`（`corrtest.py:corrtest`）。
* `anchor_node=None`（默认）完全跳过绝对速率一步，此时软件不使用 `anchor_time`。
  节点 id 从 1 开始，所以 `0` 在这里**不是**“不设锚点”的哨兵（与 R3F 不同，R3F
  的 `anchor.node = 0` 表示不设锚点）：`anchor_node=0` 抛 `ValueError`，不是内群
  树内部节点的 id 抛 `KeyError`。
* CLI：`--anchor-node INT`（不传 = `None`）与 `--anchor-time FLOAT`
  （默认 `0.0`）。
* 用于围绕一个已知节点年龄把相对速率换算为绝对速率；绝对速率的均值与标准差追加
  进报告。

### 参考值与精度（不由用户设定）
* 门 G6 与 R3F 在随包基准树上的 CorrScore 比较。R3F 只打印
  `format(score, digits = 5)`，所以
  `data/golden/r3f/example_sr0_corrtest.txt` 里是 `0.9996`；未舍入的参考值
  `0.99959505940763882` 存于
  `data/golden/r3f/example_corrtest_full_precision.txt`（由
  `reproduce/golden_r3f_fullprec.R` 捕获）。G6 断言与未舍入值逐位相等，并额外
  断言 CorrTest 所用的每个速率都等于引擎自己的逐节点速率。
* 使用 `--sister-resample > 0` 时只能在 P 值档位上比较：R3F 的重抽样随机数来自
  R 的 RNG，即使种子相同，OpenRelTime 也无法逐位复现。

## ddBD（`ddbd`）

### `sampling_frac`
* `None`（默认）报告拟合得到的 rho。给定某个值时，该值是**报告**用的抽样比例：
  与 R3F 一致（`ddbd.R` 491～505 行），似然优化始终以 rho 为自由参数（边界
  `(0, 1)`），给定值只是替换报告中的抽样比例，出生率与死亡率仍是自由 rho 拟合的
  最优点。也就是说 `sampling_frac` **不**约束拟合，只是替换一个报告值；自由 rho
  的最优值在报告参数 `sampling_frac_fitted` 中仍可见。
* 取值范围 `(0, 1]`，否则抛 `ValueError`。注意本项目不接受 R3F 表示“待估计”的
  哨兵值 `0`，请用 `None`。

### `anchor_node`、`anchor_time`
* 默认：`anchor_node=None`、`anchor_time=1.0`（`ddbd.py:ddbd`）。注意与
  `corrtest` 不对称，`corrtest` 的 `anchor_time` 默认 `0.0`。
* CLI：`--anchor-node INT`（不传 = `None`）、`--anchor-time FLOAT`
  （默认 `1.0`）。
* `anchor_node=None`（默认）保持相对时间不缩放（缩放因子 1），因此报告的速率是
  “每单位相对时间”。与 `corrtest` 相同，`anchor_node=0` 抛 `ValueError`（节点 id
  从 1 开始；`0` 是 R3F 的“不设锚点”哨兵，不是本项目的合法 id），非内部节点的
  id 抛 `KeyError`。
* 给定锚点时，软件缩放相对时间，使该节点具有指定年龄。为数值稳定，请保持最大锚
  点年龄 ≤ 10 个时间单位（R3F 建议）。

### `measure`
* `"SSE"`（默认）或 `"KL"`：初值网格的打分准则。网格：出生率 1.1～10.1 步长 1，
  死亡率 1～10 步长 1，抽样比例 {0.001, 0.01, 0.1, 0.5, 0.9}，并过滤为
  birth ≥ death；用 L-BFGS-B 在 birth-death 似然上优化，按打分顺序尝试初值
  （最多 50 个）。

## BLB（`blb`）

### `gamma`
* 默认 `0.7`，允许范围 [0.6, 0.9]：子采样规模 `ceil(L**gamma)`。
* 默认值落在已发表 RelTime-JA 抽样器的许可范围*之内*（`lb_sampler.R` 第 4
  行）；但它**不**是那次已发表运行所用的取值。该运行使用 `g = 0.78`、`20` 个
  子采样乘 `20` 个复制（`lbs_codelines.r` 第 14 行），而 OpenRelTime 的默认是
  `0.7` 与 `10 x 10`。因此默认值复现的是*方法*，不是那个具体的已发表协议（见
  `docs/validation-zh.md`）。

### `n_subsample`、`n_replicate`
* 默认 10 与 10（即 100 次 IQ-TREE 运行）。`seed` 控制全部重采样，包括每次
  IQ-TREE 的 `-seed`。这些是便利默认值，不是 RelTime-JA 的 `20 x 20` 协议；要
  匹配已发表的 BLB 预算请调高它们（并把 `gamma` 设为 0.78）。

### `iqtree_args`
* 类型：字符串列表，例如 `["-m", "GTR+G", "-B", "1000"]`。
* 默认：`None`（`bootstrap.py:blb`），按“无额外参数”处理
  （`iqtree_args or []`）；空列表等价。
* CLI：`--iqtree-arg` 可重复
  （`--iqtree-arg=-m --iqtree-arg=GTR+G`），不传得到空列表。
* 原样转发给 IQ-TREE。记号经过严格白名单校验，且进程不经 shell 启动。

### `workdir` 与 `output_dir`（BLB）
* CLI `--workdir PATH` 与 `-o/--output-dir PATH`（默认 `ORT_blb`）。
  中间产出的重采样比对与树写入 `workdir`（默认是事后删除的临时目录）；
  `blb_summary.csv` 与 `blb_replicates.csv` 写入 `output_dir`，目录不存在时创建。
  需要逐复制审计时请指定 `--workdir`。

## MEGA-CC 桥接（`megacc`）

### `megacc`
* CLI `--megacc PATH`，必填：`megacc` 可执行文件的路径。CI 既不安装也不下载
  它——见 `docs/validation-zh.md`。

### `calibrations`
* CLI `-c/--calibrations PATH`，可选：与 `calibrate` 相同的
  `calibrations.tsv` 格式。不给出时桥接只按枝长运行 RelTime
  （`ppRelTimeBLens = true`）。

### `mao`
* CLI `--mao PATH`，可选：用手写的 `.mao` 分析描述文件替代生成的最小控制文件。
  软件把生成的文件写在 `workdir` 下，并在运行输出中给出该文件的路径。

### `workdir` 与输出
* CLI `--workdir PATH`：`.mao` 与 megacc 的 stdout、stderr 落在这里。
* 结果写入 `<prefix>_megacc_times.csv`，是与 R3F 同构的节点时间表，该表的
  `source` 列记录每一行的来源。所有参数都经白名单校验，且 megacc 不经 shell 启
  动；若 megacc 没有产生可解析的输出，运行会报错，而不会悄悄改用 OpenRelTime 自
  己的估计。

## 数值常数（进阶）

本项目不把这些常数暴露为选项。列在这里是因为它们决定可复现行为，并且与 R3F 对
应。全部位于 `openreltime/_constants.py`，并附有来源注释。

| 常数 | 取值 | 含义 | 来源 |
|---|---|---|---|
| `EPS_R3F` | 1e-19 | 在 R3F 写 `10e-20` 的位置精确使用的零枝替代值 | `rrf_times.R` |
| `EPS_R3F_SMALL` | 1e-20 | 在 R3F 写 `1e-20` 的位置使用的零枝替代值 | `rrf_times.R` |
| `ZERO_BRLEN_WARN_FRACTION` | 0.10 | 触发警告的零枝占比 | `rrf_times.R:67` |
| `RATE_RATIO_THRESHOLD` | 20.0 | `rrf_times` 与 `rates-times` 使用的速率比保护阈值；`rate_ratio_threshold` 唯一注册的默认值 | `rrf_times.R:354-356` |
| `RATE_RATIO_GUARD_DISABLED` | `None` | 关闭保护的哨兵值；纯速率计算过程不运行保护 | `rrf_rates.R:250-336` 对比 `rrf_times.R:353-373` |
| `CALIBRATION_MAX_ITER` | 100 | 上报冲突校正之前，违反传播的轮数上限 | msz236 与 R3F |
| `EFFECTIVE_BOUNDS_REPLICATES` | 10,000 | 有效边界过程的复制数 | msz236 |
