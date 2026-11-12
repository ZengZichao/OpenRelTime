# 插图源代码

**语言：中文** · [English](README.md)

OpenRelTime 论文四张插图的源代码与数据：**Figure 1**（分析流程）、**Figure S1**
（软件架构）、**Figure S2**（与 R3F 的一致性及墙钟时间）与 **Figure S3**（有效边界、
校正敏感度，以及随包示例校正文件中的不自洽）。文件名沿用稿件自身的命名、未作改动，
以便与投稿所用的成品图保持可追溯的对应关系。

| 文件 | 产出 |
|---|---|
| `make_figures.py` | 正文图1-分析流程图 (Figure 1) —— 同时提供 Figure S1 所调用的共享绘图引擎 |
| `附件图S1-软件架构图.py` | 附件图S1-软件架构图 (Figure S1) |
| `附件图S2-性能对比图.py` + `附件图S2-性能对比图-data.json` | 附件图S2-性能对比图 (Figure S2) |
| `附件图S3-置信区间图.py` + `附件图S3-置信区间图-data.json` | 附件图S3-置信区间图 (Figure S3) |

渲染出的 `.pdf`/`.png`/`.drawio` 成品图随论文投稿提交，刻意**不**纳入本仓库；
每个脚本都会重新测量自己交付的文件（页面宽度、栅格 dpi、bbox 裁剪）并打印报告。

## 运行

渲染只需要 matplotlib 和 numpy——即本项目的 `[plot]` extra；Figure S2 另用
pandas 与 numpy 做金标准文件比对。S2 与 S3 的渲染路径从不 import
`openreltime`，因此即使一台机器根本没有安装本包，也能复现已发表的图。只有
`--regen-data` 需要安装本包（见下文）。

```bash
python figures/make_figures.py            --outdir /tmp/figs
python figures/附件图S1-软件架构图.py       --outdir /tmp/figs
python figures/附件图S2-性能对比图.py       --outdir /tmp/figs
python figures/附件图S3-置信区间图.py       --outdir /tmp/figs
```

**务必传入 `--outdir`。** 每个脚本的默认输出目录就是脚本自身所在的目录，因此
不加参数直接运行会把 PDF/PNG（Figure 1 与 S1 还有 `.drawio`）写进 `figures/`。

各脚本支持的命令行参数：

| 脚本 | `--outdir` | `--engine-root` | `--regen-data` | `--regen-out` | `--force` | 其他 |
|---|---|---|---|---|---|---|
| `make_figures.py` | 支持 | — | — | — | — | — |
| `附件图S1-…​.py` | 支持 | — | — | — | — | `--figure-engine` |
| `附件图S2-…​.py` | 支持 | 支持 | 支持 | 支持 | 支持 | `--golden-dir` |
| `附件图S3-…​.py` | 支持 | 支持 | 支持 | 支持 | 支持 | `--golden-dir`（可接受，但被忽略） |

* **Figure S1 与 Figure 1 共用绘图引擎。** 它 import `make_figures.py`，查找
  顺序为：`--figure-engine PATH`、环境变量 `$OPENRELTIME_FIGURE_ENGINE`、脚本
  自身目录，然后是同级目录与工作目录。本目录同时存放这两个文件，因此无需任何
  参数。
* **`--engine-root`** 指定 OpenRelTime 仓库（即包含 `openreltime/` 与 `data/`
  的目录）；也可经 `$OPENRELTIME_ENGINE_ROOT` 读取，探测时只认目录内容、绝不认
  目录名。脚本从本目录运行时，其父目录*就是*仓库，故该参数可省；仅当仓库位于
  别处时才需要。`--golden-dir` 同理指定 `data/golden/r3f`
  （`$OPENRELTIME_GOLDEN_DIR`），供 Figure S2 读取。

## 已发表插图的渲染环境

已发表的 PDF/PNG 渲染自 **Python 3.11.15 与 matplotlib 3.11.0**（numpy 2.4.6、
pandas 2.3.3、Pillow 12.2.0），平台为 macOS/arm64。这是经过验证而非口头断言的：
从本目录重新运行全部四个脚本，可**逐字节**复现四张已发表 PNG，PDF 之间也仅在
内嵌的创建时间戳上有差异。

不要把它与 `附件图S2-性能对比图-data.json` 内的 `environment` 块混为一谈，其中
记录的是 Python 3.14.6 / matplotlib 3.11.1 / numpy 2.4.4 / pandas 3.0.3：那是
*执行*墙钟时间基准测试的版本，而非*绘制*该图的版本。

## 纳入版本库的数据文件即发表记录

Figure S2 与 S3 的每一个数值都取自其旁边纳入版本库的 `-data.json` 文件；渲染
路径不读取任何其他输入。`--regen-data` 则改用已安装的 `openreltime` **在本机
重新测量**；由于纳入版本库的文件就是发表记录，刷新结果默认只会在其旁边写入
一份*副本*（`…-data.regen.json`），除非用 `--regen-out PATH` 指定文件名，或
有意识地传入 `--force`（这将直接替换原文件）。

```bash
python figures/附件图S2-性能对比图.py --regen-data --regen-out /tmp/s2.json
python figures/附件图S3-置信区间图.py --regen-data --regen-out /tmp/s3.json
```

请把这一保护机制视为数据可用性声明的一部分。`--regen-data` 是重新核验
Figure S2 分图 (a) 的正当手段——与 R3F 的双精度一致性不依赖机器——但绝不能用它
复现分图 (b)/(c)，其中的计时数据属于文件中记录的那台机器。

## 已知限制：`附件图S3-置信区间图.py --regen-data` 在 0.1.1 上会中止

用 0.1.1 引擎重新测量 Figure S3 会**失败**。其分图 (a) 与分图 (c) 的输入是
`data/examples/example_calibrations.tsv`，其中三行占位数据彼此矛盾；0.1.0 会
通过缩放各谱系速率修复这种矛盾并带警告返回结果，而 0.1.1 会拿收敛后的节点年龄
与树核对并拒绝它们（`openreltime/calibrate.py`，`_check_time_order`）。在当前
仓库上实测：

```
$ python 附件图S3-置信区间图.py --engine-root <repo> --regen-data --regen-out /tmp/x.json
...
ValueError: conflicting calibration constraints: the repair loop converged on
node ages that are not nested in the tree -- 1 node(s) older than the present
(negative age): [307]; 2 parent-child inversion(s): node 306 at 3225.03 sits
below its ancestor node 305 at 103.04; node 372 at 143.919 sits below its
ancestor node 371 at 6.23598.  This means the calibration set is unsatisfiable
by a rate rescaling of this tree; relax the bounds or move them toward the
root.
```

（出问题的节点 id 与年龄取决于随机种子；中止这一行为本身则不依赖种子。）

`附件图S3-置信区间图-data.json` **早于该保护机制**——其中记录
`"openreltime_version": "0.1.0"`——而它*就是*发表记录：已录用插图中的数值与该
文件完全一致。插图本身无需任何改动，因为分图 (c) 画的恰恰就是这项检查如今会
拒绝的那份不自洽；受影响的只有刷新路径。
[`data/examples/README-zh.md`](../data/examples/README-zh.md)
解释了为什么该文件是一份预期会失败的格式示例，而
`data/examples/example_calibrations_feasible.tsv` 是其可运行的对应版本。

## 分发方式

`MANIFEST.in` 把 `figures/*.py`、`figures/*.json` 与本 README 一并纳入**源码
分发**，沿用它对 `data/`、`docs/`、`reproduce/` 已有的同一模式。这些目录不是
Python 包，setuptools 既不能也不会把它们放进 wheel：安装后的布局不变
（`openreltime` 包加 `openreltime`/`ort` 两个控制台脚本），插图源码经 sdist 或
仓库检出到达读者手中。渲染成品被排除在外，因此在本地用 `--outdir figures`
跑出的文件不会泄漏进由工作树构建的分发包。

## 代码风格

四个脚本与两个数据文件是产出已发表渲染图的源文件的**逐字节相同副本**，因此
不按本仓库的 ruff 配置重新格式化；`[tool.ruff].extend-exclude` 中列有
`figures`，lint 门槛因此维持在原有数量。它们也位于包的 import 图之外、测试
套件之外。
