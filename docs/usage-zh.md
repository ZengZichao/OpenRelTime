# OpenRelTime 使用说明

**语言：中文** · [English](usage-en.md)

版本：v0.1.0 · 环境：Python ≥ 3.10 · 许可：GPL-3.0-or-later

本说明面向日常使用，覆盖安装、输入准备、全部子命令、校正文件格式、输出
文件清单、独立发行的桌面图形界面指引、常见问题与故障排查。全部可选参数的语义见
[参数手册](parameters-zh.md)；算法原理与局限见 [方法说明](methods-zh.md)；
验证门槛见 [验证说明](validation-zh.md)。

---

## 1. 安装

### 1.1 从源码目录安装

> **分发方式。** 本版本以源码形式分发，请在检出目录中用
> `pip install --editable ".[dev,plot]"` 安装；发布到包索引后会在此处公告。

```bash
git clone https://github.com/ZengZichao/OpenRelTime
cd OpenRelTime
python -m pip install --editable ".[dev,plot]"   # 核心 + 绘图 + 开发工具
```

`plot`（matplotlib 绘图）与 `dev`（pytest、pytest-cov、ruff、mypy，面向开发者）
即 `pyproject.toml` 中 `[project.optional-dependencies]` 声明的全部 extras；解析
置信区间只用 NumPy 与 SciPy，不需要额外依赖。需要多个可选依赖时合并书写，例如
`pip install --editable ".[plot,dev]"`。桌面图形界面**不是**本包的 extra——它是独
立的发行包 `OpenRelTime-Studio`，反过来依赖本引擎（见 §7）。

安装后可用 `openreltime --version`（或短别名 `ort --version`）验证。在 conda
或 micromamba 环境中，先激活目标环境再安装，包才会落进该环境：

```bash
micromamba activate <your-env>      # 或 conda activate <your-env>
pip install --editable ".[dev,plot]"
```

本包无编译扩展，在 Linux、macOS 与 Windows 上安装都只需数秒。

### 1.2 在检出目录之外使用源码

对检出的源码目录（checkout）做**可编辑安装**，是在仓库之外运行分析、基准测试或绘
图脚本时的推荐做法：可编辑安装完成后，检出的源码目录（即 `openreltime` 这个包，
以及 `openreltime`、`ort` 两个命令行入口）在任意工作目录下均可导入，脚本直接写

```python
import openreltime as ort
```

即可。**不要**手工把源码目录路径插进 `sys.path`（例如
`sys.path.insert(0, "/path/to/OpenRelTime-source")`）：这种做法绕开了声明的
依赖集合，命令行入口不可用，还会遮蔽已安装的包。由于是可编辑安装，修改源码后
立即生效，无需重装。确认脚本实际导入的是哪份代码：

```bash
python -c "import openreltime; print(openreltime.__version__, openreltime.__file__)"
```

## 2. 输入准备

### 2.1 系统树（必需）

带枝长的**有根二叉** Newick（或 NEXUS）树，枝长单位为每位点替换数（IQ-TREE、
RAxML 与 FastTree 的输出可直接使用）。要求：

- 至少 3 个内群分类群；
- 枝长非负且有限；缺失枝长按 0 处理（会给出警告）；
- 多分叉节点默认报错；`--resolve random` 可随机二叉化（新增枝长置 0）。

### 2.2 外群（推荐）

外群用于置根，置根后软件剔除外群，所有输出仅含内群。支持两种形式：

- **单末端外群**：只写一个末端名称；
- **多末端外群**：写多个名称（逗号分隔或每行一个的文本文件）。

```bash
# 单末端
--outgroup Ornithorhynchus_anatinus
# 多末端
--outgroup "Ornithorhynchus_anatinus,Zaglossus_bruijni,Tachyglossus_aculeatus"
--outgroup outgroup.txt
```

**不提供外群时**：该选项可以省略，此时既不置根也不剔除，树按读入原样使用。
RelTime 的一切量都锚定在根上，输入文件的根位置因此直接决定所有节点年龄，而
OpenRelTime 对此**不做任何警告**。请务必确认输入树已按预期的单系群置根（R3F 强制
要求外群；见 `docs/validation-zh.md` §5）。

**完整簇检查**：多末端外群必须构成一个完整的簇（单系群）——回溯到这些分支的最
近共同祖先（MRCA）后，该祖先下的**所有**末端都必须属于所指定的外群。若 MRCA 下
还混入了未指定的末端，软件默认报错并列出这些多余末端（`--outgroup-check warn`
可降级为警告并继续，置根位置仍为 MRCA）。单末端外群天然完整，不做此检查。

**按单系群名称指定外群**：`--outgroup` 的每个名称先按末端名匹配；若不是末端名，
则按类群名称解析（见 §2.4 分类信息），展开为该群全部末端并做同样的完整簇检查：

```bash
--outgroup Monotremata --taxon-table taxa.tsv
```

### 2.3 校正文件（calibrate 子命令需要）

制表符分隔的 TSV（见 §4），或 Python 端构造
`openreltime.Calibration` 对象列表。`taxon_set` 中除末端名外也可写类群名（先查
`--taxon-table`，再尝试末端标签前缀自动识别）；若该类群经判断不是单系群，软件会
发出警告，并仍以该 MRCA 为校正节点。

### 2.4 末端分类信息（单系性判断）

`monophyly` 子命令、群名外群与群名校正都需要“末端 → 类群”的对应关系，来源有两
种：

1. **用户提供分类信息表**（`--taxon-table`）：TSV 或 CSV，两列，表头可用 `tip`、
   `label`、`name` 之一与 `group`、`clade`、`taxa`、`classification` 之一（示例
   见 `data/examples/example_taxa.tsv`）；
2. **自动从树文件末端标签识别**（默认）：按 ``Genus_species`` 惯例，取每个末端
   标签第一个分隔符（`_`、空格、`|`、`@`）之前的词作为类群名（属名），例如
   `Mus_musculus` → 类群 `Mus`。

## 3. 命令行用法

### 3.1 相对速率与相对时间

```bash
# 仅相对速率（谱系速率，根 = 1.0）
openreltime rates -i tree.nwk --outgroup og.txt -o mam

# 仅相对节点时间（tip = 0）
openreltime times -i tree.nwk --outgroup og.txt -o mam --r3f-compat

# 一次算齐（速率 + 时间 + NEXUS 时间树 + 速率着色图）
openreltime rates-times -i tree.nwk --outgroup og.txt -o mam --plot mam_tree.png
```

- `--normalize`：将最大节点时间归一为 1.0（MEGA 与 PNAS 2012 语义）；
- `--mean geometric|arithmetic`：选择平均约定；几何平均与 R3F 一致，也是各处
  的默认值；
- `--no-guard`：关闭速率比阈值保护（默认 20，与 R3F 相同）；
  `--rate-ratio-threshold FLOAT`：设定该阈值。该保护只作用于**时间**计算过程：
  `rates` 出于签名对称而接受 `--rate-ratio-threshold`，但会忽略这个参数，并记录
  一条警告说明此事；
- `--r3f-compat`：输出的文件名与列名同 R3F 完全一致（`_RRF_times.csv` 等），
  便于流水线替换。

`--fmt`、`--resolve`、`--taxon-table`、`--outgroup-check`（见 §2）在所有基于树
的子命令上可用。`docs/parameters-zh.md` 开头给出了各子命令支持哪些选项的对照
表。

### 3.2 校正与绝对时间

```bash
openreltime calibrate -i tree.nwk -c calibrations.tsv --outgroup og.txt \
    --method effective --n-effective 10000 --seed 42 -o mam_cal
```

- `--method bounds`（默认）：`min` 与 `max` 硬边界；全局因子 f 取可行域中点，
  无法满足时按谱系速率缩放迭代调整（上限 100 次），仍不可行则报错并列出冲突
  节点；
- `--method effective`：密度校正按 msz236 抽样转换，`--n-effective` 控制复制数
  （默认 10,000），输出有效边界表 `<prefix>_effective_bounds.csv` 与各节点经验
  分位区间。只给密度而不给硬边界的行在 `bounds` 下不构成任何约束，解析阶段会
  拒绝这类行，因此这类行必须改用 `effective`；
- `--mean`、`--rate-ratio-threshold`、`--no-guard`：决定被校正的那次 RRF 估计。
  传入的取值都会写入报告，`openreltime ci` 会重放这些取值，从而保证置信区间与
  点估计出自同一套估计器。

### 3.3 解析置信区间

```bash
# 基于上一步的 calibrate 前缀
openreltime ci -c mam_cal --n-sites 1000 -o mam_ci
```

- `-c` 传入的是 `calibrate` 那一步的**输出前缀**而非树文件：树、外群与输入
  处理方式都从 `<prefix>_report.json` 读回；
- `--branch-var var.tsv`：用户提供每条边的抽样方差（只含 `node_id`、`var`
  两列）；
- `--n-sites L`：Poisson 近似 `vS(b)=b/L`；
- 两者都不给时仅含速率异质性方差（msz236 模拟协议）；
- `--level 0.95`：置信水平；CI 会在硬边界处截断，并把可能为负的下界截到 0；
- `--mean`、`--rate-ratio-threshold`、`--no-guard`：覆盖 `calibrate` 记录的
  设置。不指定时，区间沿产生该点估计的同一估计器求导——`--method effective`
  的结果也按其有效边界求导；
- `<prefix>_ci.csv` 的列为 `node_id, label, time, se, lower, upper, width,
  se_reliable, notes`。节点年龄若恰好钉在自身校正上，该节点的 eq-(7) 导数恰为
  0，因此软件改按该校正本身断定的不确定度报告，并标记
  `se_reliable=False`；`notes` 说明原因，报告 JSON 中以 `n_nodes_se_not_reliable`
  计数。

### 3.4 CorrTest 速率自相关检验

```bash
openreltime corrtest -i tree.nwk --outgroup og.txt -o mam
# 小树（<50 tips）建议重抽样：
openreltime corrtest -i tree.nwk --outgroup og.txt --sister-resample 100 --seed 1 -o mam
```

输出 `<prefix>_corrtest.txt`：CorrScore ∈ [0,1] 与 P 值档位
（≥0.92 → P<0.001；[0.83,0.92) → P<0.01；[0.5,0.83) → P<0.05；<0.5 →
P>0.05）。CorrScore ≥ 0.5 即提示存在显著的速率自相关，RelTime 类方法的
对数速率假设可能不成立。

### 3.5 ddBD 物种分化树先验

```bash
openreltime ddbd -i tree.nwk --outgroup og.txt --anchor-time 1.85 -o mam_bd
openreltime ddbd -i tree.nwk --outgroup og.txt --sampling-frac 0.5 -o mam_bd
```

- `--anchor-node N --anchor-time T`：把某内部节点的相对时间标定为 T
  （时间单位建议使最大节点年龄 ≤ 10）；
- `--measure SSE|KL`：初值网格选择准则；
- `--sampling-frac 0.5`：设定**报告**的抽样比例。与 R3F 一致，似然优化始终把
  rho 留在自由参数中（边界 (0,1)），给定值只替换报告中写出的抽样比例，并不约束
  优化；birth 与 death 仍是自由 rho 拟合的最优点，自由拟合值可在报告的
  `sampling_frac_fitted` 中查看；
- `--plot curve.png`：输出节点时间分布直方图与拟合密度曲线。

输出表：`birth.rate / death.rate / sampling.frac`，可直接作为 MCMCTree 与
BEAST 的树先验参数。

### 3.6 tree2table

```bash
openreltime tree2table -i tree.nwk -o mam_tbl          # 枝长表
openreltime tree2table -i tree.nwk --time -o mam_tbl   # 节点年龄表
```

### 3.7 单系性判断（monophyly）

判断类群是否为单系群（完整簇）。分类信息来自 `--taxon-table` 或末端标签自动识
别；`--groups` 指定要检验的类群名（默认检验全部），`--tips` 可检验任意末端集合：

```bash
# 自动按属名识别类群，检验全部
openreltime monophyly -i tree.nwk -o chk

# 检验指定类群（分类信息表）
openreltime monophyly -i tree.nwk --taxon-table taxa.tsv --groups Monotremata,Primates -o chk

# 检验任意末端集合
openreltime monophyly -i tree.nwk --tips A,B,C -o chk
```

输出 `<prefix>_monophyly.csv`：每个类群一行，含末端数、MRCA 末端数、多余末端
（`extra_tips`）与判定 `is_monophyletic`。对非单系类群，同时在 stderr 发出警告。
非单系判定同样适用于群名外群（§2.2）与群名校正（§4）。不提供 `--taxon-table`
时按末端标签前缀自动识别类群，`--group-sep REGEX` 可改变切分模式（默认
`[_|@\s]`，即取 `_`、`|`、`@` 或空白中的第一个）。

### 3.8 BLB 管线（需要外部 IQ-TREE）

```bash
openreltime blb -a alignment.fasta --iqtree /path/to/iqtree \
    --outgroup og.txt --gamma 0.7 --n-subsample 10 --n-replicate 10 \
    --iqtree-arg=-m --iqtree-arg=GTR+G --seed 42 -o blb_out
```

两级重采样：① 无放回抽取 ⌈L^0.7⌉ 个位点；② 每个子采样内有放回重抽 L 个位点。
其后是 IQ-TREE 建树、OpenRelTime 定时，最后按节点 clade 聚合，给出中位时间以及
2.5% 与 97.5% 分位区间（`blb_out/blb_summary.csv`）。

`topology_ok` 列给出包含该 clade 的复制比例，`blb_replicates.csv` 保留全部复制
的逐棵树结果。`--workdir` 用于保留中间的重采样比对与树（默认写入临时目录并在结
束后删除）；`--seed` 同时以 `-seed` 传给每一次 IQ-TREE 运行。默认
`--gamma 0.7`、`10 x 10` 落在 BLB 许可范围内，但**不**复现已发表的 RelTime-JA
运行（`g = 0.78`、`20 x 20`），详见 `docs/parameters-zh.md` 与
`docs/validation-zh.md`。

### 3.9 MEGA-CC 桥接（可选）

```bash
openreltime megacc -i tree.nwk -c calibrations.tsv --megacc /path/to/megacc \
    --outgroup og.txt --workdir megacc_run -o mam_mega
```

软件自动生成 `reltimeFromBranchLengths.mao`，调用 megacc，再把 megacc 的输出解析
为与 OpenRelTime 和 R3F 同构的节点时间表 `<prefix>_megacc_times.csv`，其中
`source` 列逐行标注来源。`-c/--calibrations` 可以省略（此时只按枝长跑
RelTime），`--mao custom.mao` 可用手写的控制文件替代生成文件。若 megacc 没有给出
可解析的输出，运行会直接报错——绝不会拿 OpenRelTime 自己的估计冒充 MEGA 结果。

桥接有两个限制。第一，它需要真实的 `megacc` 可执行文件，CI 既不安装也不运行该程
序，因此 MEGA-CC 的一致性依据已发表的 R3F/MEGA 结果与随包黄金文件确立，而不是依据
此处新做的头对头实测。第二，OpenRelTime 的解析 RRF 层缺少 MEGA 的
似然比精修，与 MEGA 的逐节点差异属预期之内。把它当作对照工具，不要当作等价性声
明（见 `docs/validation-zh.md`）。

## 4. 校正文件格式

制表符分隔（示例见 `data/examples/example_calibrations.tsv`——**仅演示文件格
式**，其中的数值是随意填写的示意值，见 `data/examples/README-zh.md`）。示例与黄金
数据随源码分发（sdist）发布，也随可编辑安装的源码目录（`data/`、`MANIFEST.in`）发
布；纯 wheel 不含它们，请在源码目录内运行示例。首行为表头；每行都含六个制表符分隔
字段（`.` 表示空单元格）。

| 列 | 必需 | 说明 |
|---|---|---|
| `node_id` | 二选一 | 内部节点 id（同 `<prefix>_times.csv` 的 NodeId） |
| `taxon_set` | 二选一 | `\|` 分隔的末端名，或单个类群名（`--taxon-table` 或自动识别，见 §2.4），取其 MRCA；类群经判断非单系时发出警告 |
| `min_bound` | 至少给一个边界或密度 | 最小年龄；`.` 表示缺失 |
| `max_bound` | | 最大年龄；`.` 表示缺失 |
| `density` | 可选 | `uniform`、`exponential`、`normal` 和 `lognormal` |
| `density_params` | 密度时必需 | `key=value;key=value`，如 `offset=60;mean=20` |

各密度参数键：uniform `min,max`；exponential `offset,mean`；normal
`mean,sd`；lognormal `offset,meanlog,sdlog`。参数仅做严格的键值解析，不做
表达式求值。

## 5. 输出文件清单

| 后缀 | 由谁写出 | 内容 |
|---|---|---|
| `_rates.csv` | `rates` | NodeLabel、NodeId、Des1、Des2、Rate |
| `_rates.nwk` | `rates` | 以枝注释携带速率的同一棵树 |
| `_times.csv` | `times` | 同上，列为 Time |
| `_timetree.nwk` | `times` | 相对时间树（枝长 = 父子时间差） |
| `_rates_times.csv` | `rates-times` | Rate 与 Time 两列 |
| `_timetree.nexus` | `rates-times` | 带 `[&rate=...]` 注释，FigTree 可直接打开 |
| `_RRF_*` 系列 | `--r3f-compat`（`times`、`rates-times`）或 API `write(r3f_compat=True)` | `_RRF_rates.csv`、`_RRF_times.csv`、`_RRF_timetree.nwk`、`_RRF_table.csv`、`_RRF_timetree.nexus` |
| `<prefix>.csv` | `tree2table` | 节点表（枝长；加 `--time` 则为节点年龄） |
| `_monophyly.csv` | `monophyly` | group、n_tips、is_monophyletic、n_mrca_tips、extra_tips、missing_names、tips |
| `_calibrated.csv/.nwk/.nexus` | `calibrate` | 绝对时间树与表 |
| `_effective_bounds.csv` | `calibrate --method effective` | 有效边界（msz236）与经验分位区间 |
| `_ci.csv` | `ci` | node_id、label、time、se、lower、upper、width、se_reliable、notes |
| `_corrtest.txt` | `corrtest` | CorrScore 与 P 值档位 |
| `_ddbd.txt` | `ddbd` | BD 参数 |
| `blb_summary.csv`、`blb_replicates.csv` | `blb` | 逐 clade 汇总；全部复制表 |
| `_megacc_times.csv` | `megacc` | MEGA-CC 节点时间，含 `source` 列 |
| `_report.json` | 见下 | 全部运行参数、警告与随机种子 |

`_report.json` 只由参数会影响数值的分析步骤写出：`rates`、`times`、
`rates-times`、`calibrate`、`ci`、`corrtest`、`ddbd`。共用同一前缀的多次运行会
**合并**进同一个报告文件，各占一个顶层键，因此整条流水线仍可集中检查与复现。
`tree2table`、`monophyly`、`blb`、`megacc` 只写出自己的结果表。

## 6. Python API 速查

```python
import openreltime as ort

tree  = ort.read_tree(path, fmt="newick", outgroup=[...])
rates = ort.rrf_rates(tree, mean="geometric", rate_ratio_threshold=None)
times = ort.rrf_times(tree, rate_ratio_threshold=20.0, normalize=False)
rt    = ort.rrf_rates_times(tree)
cal   = ort.calibrate(times, cals, method="effective", n_effective=10_000, seed=42)
ci    = ort.confidence_interval(cal, seq_length=1000, level=0.95)
corr  = ort.corrtest(tree, sister_resample=100, seed=1)
bd    = ort.ddbd(tree, anchor_time=1.85, measure="SSE")
table = ort.tree2table(tree, time=False)

# 单系性判断（分类表或末端标签前缀自动识别）
table_map  = ort.parse_taxon_table("taxa.tsv")
verdict    = ort.check_monophyly(tree, ["Mus_musculus", "Mus_spretus"], name="Mus")
verdict.is_monophyletic
groups     = ort.detect_groups_from_labels(tree.tip_labels())   # {属名: [末端...]}
# 群名外群 / 群名校正：outgroup=、taxon_set= 中写类群名即可
tree2 = ort.read_tree(path, outgroup=["Monotremata"], taxon_table=table_map)
cal2  = ort.calibrate(times, cals, taxon_table=table_map)
```

所有结果对象都有 `to_pandas()`、`to_json()` 与 `write(prefix)`，与 CLI 输出同源。

## 7. 图形界面——OpenRelTime Studio（独立产品）

本引擎的点击式图形前端 **OpenRelTime Studio** 不属于本仓库，也不从本仓库安装。
它是一个独立发行的产品，依赖本引擎的公开 API——与任何第三方脚本的用法完全一
致。因此在界面里得到的数值，就是用相同输入跑 CLI 或 API 得到的数值：分析代码始
终在本引擎内，装或不装图形界面都不会改变它。

```bash
pip install OpenRelTime-Studio     # 会自动装入本引擎作为依赖
openreltime-studio
```

* 主页：<https://github.com/ZengZichao/OpenRelTime-Studio>
* 界面能力：完整的中英双语界面（「视图 ▸ 语言」即时切换）、浅色/深色两套主题
  （「视图 ▸ 主题」）、根节点四向可调（「视图 ▸ 根节点方向」）、手绘矢量 SVG 图标，
  以及一棵自带演示校正点的内置示例树（「文件 ▸ 打开内置示例…」）。
* 功能范围：加载 Newick 或 NEXUS 树，运行 RRF，在画布上点选内部节点设置校正
  点，再依次完成校正、解析置信区间、CorrTest 与 ddBD；导出 CSV、NEXUS、JSON、
  PNG 四种结果与等效命令行。
* 窗口布局、各运行按钮的前置条件链、导出具体写了哪些文件，以及其余的界面操作
  手册：**见该项目自己的文档**。本说明只描述引擎。
* 它读入的校正文件就是 §4 的格式，它导出的文件就是 §5 清单里的文件。

## 8. 常见问题（FAQ）

**Q1 与 MEGA 的 RelTime 结果完全一致吗？**
不完全。OpenRelTime 实现解析 RRF 层（与 R3F 同层）；MEGA 在此之上还有
LRT 速率约束精修。二者应高度相关（相关系数报告为 0.98～0.99），但不会逐
位一致。详见 `docs/methods-zh.md` §2。

**Q2 与 R3F 一致到什么程度？**
同一算法、同一 epsilon 位置。G2a 与 G2b 门在 274-tip 哺乳动物树上断言回归斜率
≥ 0.999、中位相对偏差 ≤ 1e-6。实测中位约为 6.8e-16（速率）与 8.7e-16（时间），
最大偏差约为 4.8e-15 与 3.6e-14，这些最大值属于观测值，不是该门断言的内容（见
`docs/validation-zh.md`）。

**Q3 树有 10% 以上零枝长怎么办？**
软件会发出警告。零枝长使相连节点的速率与时间不可靠，建议改进树推断或删除质量差
的序列；OpenRelTime 的 epsilon 处理与 R3F 一致，可做敏感性对比（`--no-guard` 或
不同阈值）。

**Q4 多分叉树报错？**
这是有意设计（解析解只对二叉树成立）。确认拓扑后再用
`--resolve random`，报告中会记录受影响节点。

**Q5 BLB 管线很慢？**
BLB 需要跑 `n_subsample × n_replicate` 次 IQ-TREE。先估计单次建树耗时，
按预算选择网格；OpenRelTime 定时部分单树亚秒级（1,000 tips < 1 s）。

**Q6 如何引用？**
见 `CITATION.cff`，并请同时引用 msy044（RRF）、msz236（CI 与 effective
bounds）、msz014（CorrTest）与 btab307（ddBD）。

## 9. 故障排查

| 报错与现象 | 原因与处理 |
|---|---|
| `tree contains polytomies at ...` | 用 `--resolve random` 或先二叉化 |
| `negative branch length at node ...` | 检查树文件；负枝长非法 |
| `duplicate tip names` | 去重后再运行 |
| `ingroup has N tips after outgroup removal` | 内群不足 3 个分类群，调整外群 |
| `outgroup [...] is not a complete clade` | 多末端外群的 MRCA 下混有未指定末端（非完整簇）；补全外群成员，或改用 `--outgroup-check warn` 降级为警告 |
| `outgroup taxon [...] not found in tree (as tip or group)` | 外群名称既不是末端名也不是类群名；检查拼写或提供 `--taxon-table` |
| `conflicting calibration constraints` | 校正节点互相冲突（如子节点边界早于父节点），放宽边界 |
| `calibration constraints could not be satisfied within 100 iterations` | 同上，见报错中的节点列表 |
| `a calibration needs at least one finite bound … but neither was given` | 该行的 `min_bound` 与 `max_bound` 都是 `.`；补上硬边界，或改用 `--method effective` 处理只给密度的行 |
| `Under method='bounds' a density does not count as a bound` | `bounds` 下给了纯密度校正；改用 `--method effective` |
| `no 'calibrate' block in <prefix>_report.json` | `ci -c` 拿到的是一个从未跑过 `calibrate` 的前缀；用同样的 `-o` 先跑 `calibrate` |
| `cannot locate the original tree` | 树文件在 `calibrate` 之后移动了；用相同的 `-i` 与 `-o` 重跑 `calibrate` 以记录 `tree_file` |
| `group '…' not found in …; available: …` | 该名称既不在 `--taxon-table` 中，也不是自动识别出的前缀；报错会列出可用名称 |
| `sampling_frac must lie in (0, 1]` | `0` 是 R3F 表示“待估计”的哨兵值；请改为省略 `--sampling-frac` |
| `anchor node id must be a positive internal-node id` | 本项目的节点 id 从 1 开始，`0` 不是“不设锚点”的哨兵值 |
| `executable not found: iqtree/megacc` | 用完整路径传入可执行文件 |
| `argument ... contains characters outside the allowed set` | 外部工具参数含 shell 元字符；请改用不含空格与分号的路径 |
| CorrScore = 0 或 NaN | 树太小或极端速率过多；先 `--sister-resample`，仍异常则检查枝长 |

## 10. 可复现性

- 黄金标准与生成脚本：`data/golden/r3f/`、`reproduce/golden_r3f.R`
  （R 4.5.3、ape 5.8.1、phangorn 2.11.1、R.utils 2.13.0、FNN 1.1.4.1、
  RColorBrewer 1.1-3，见 `data/PROVENANCE.md`）；
- 测试：在仓库根目录运行 `python -m pytest`，收集范围是引擎测试
  （`tests/`）。桌面前端的自带测试随独立产品 `OpenRelTime-Studio` 一起维护（见
  §7）；
- 持续集成：`.github/workflows/ci.yml` 在 Ubuntu 与 macOS 上以 Python
  3.10、3.11、3.12 和 3.13 运行上述测试套件。Windows 是声明的支持目标，CI 不在该平台
  运行；外部工具链（`blb` 的 IQ-TREE、`megacc` 的 MEGA-CC）不在工作流中安装——
  已测与仅声明的对照见 `docs/validation-zh.md`；
- 随机性（有效边界、BLB、sister 重抽样）一律经 `--seed` 注入并写入
  `_report.json`。
