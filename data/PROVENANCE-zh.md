# 数据来源记录

**语言：中文** · [English](PROVENANCE.md)

`data/` 下每一个第三方文件都在此登记其来源、许可与应引用的文献。黄金输出可用
`reproduce/golden_r3f.R` 重新生成（R 4.5.3、ape 5.8.1、phangorn 2.11.1、
R.utils 2.13.0、FNN 1.1.4.1、RColorBrewer 1.1-3）。

| 文件 | 来源 | 许可 | 引用 |
|---|---|---|---|
| `data/examples/example.nwk` | R3F（Tao, Sharma, Tamura & Kumar 2025）的 `data/example.nwk`，原始出处为 dos Reis et al. 2012（由 274 条线粒体基因组序列构建的 274-tip 哺乳动物树，7,370 个位点） | GPL-3（可再分发） | dos Reis M, Inoue J, Hasegawa M, Asher RJ, Donoghue PCJ, Yang Z. *Proc R Soc B* 279:3491-3500 (2012), doi:10.1098/rspb.2012.0683 |
| `data/golden/r3f/example_RRF_rates.csv` | 由 R3F `rrf_rates()` 生成（本仓库的 `reproduce/golden_r3f.R`） | GPL-3 | Tao Q, Sharma S, Tamura K, Kumar S. (2025) R3F, arXiv:2502.05004；软件包版本 0.1.0 |
| `data/golden/r3f/example_RRF_rates.nwk` | 由同一个 R3F `rrf_rates()` 调用生成 —— 它是该函数的第二个产物，一条分支长度即相对谱系速率的树 | GPL-3 | 同上 |
| `data/golden/r3f/example_RRF_times.csv` | 由 R3F `rrf_times()` 生成 | GPL-3 | 同上 |
| `data/golden/r3f/example_RRF_table.csv` | 由 R3F `rrf_rates_times()` 生成 | GPL-3 | 同上 |
| `data/golden/r3f/example_RRF_timetree.nwk` 与 `.nexus` | 由 R3F 生成 | GPL-3 | 同上 |
| `data/golden/r3f/example_sr0_corrtest.txt` | R3F `corrtest()`，`sister.resample = 0` | GPL-3 | Tao Q. et al. *Mol Biol Evol* 36:811-824 (2019) |
| `data/golden/r3f/example_sr50_corrtest.txt` | R3F `corrtest()`，`sister.resample = 50`，`set.seed(42)` | GPL-3 | 同上 |
| `data/golden/r3f/example_corrtest_full_precision.txt` | R3F `corrtest()`（sr0 与 sr50），重新 source 时只把得分的 `format(score, digits = 5)` 放宽到 17 位，从而打印出未舍入的 double；由 `reproduce/golden_r3f_fullprec.R` 生成，**不是** `reproduce/golden_r3f.R` 的产物——这也是该目录有第十三个文件、且重生成流程不写它的原因 | GPL-3 | 同上 |
| `data/golden/r3f/example_anchor_ddbd.txt` | R3F `ddbd()`，`anchor.node = 272, anchor.time = 1.85` | GPL-3 | Tao Q. et al. *Bioinformatics* 37:i102-i110 (2021) |
| `data/golden/r3f/example_noanchor_ddbd.txt` | R3F `ddbd()`，不设锚点 | GPL-3 | 同上 |
| `data/golden/r3f/example_tree2table.csv` | R3F `tree2table()` | GPL-3 | 同上 |
| `data/golden/r3f/example_rooted_dropped.nwk` | 对 `example.nwk` 执行 `ape::root(resolve.root = TRUE)` 加 `ape::drop.tip` | 开放（程序生成） | E. Paradis & K. Schliep, *Bioinformatics* 35:526-528 (2019) |
| 3 与 4 分类群 RRF 手算值（msy044 式 28-42 与式 1-27） | 以 pytest fixture 的形式写在 `tests/test_g1_hand_computed.py` 中（不另设黄金目录） | 本项目 | Tamura K., Tao Q., Kumar S. *Mol Biol Evol* 35:1770-1782 (2018) |

生成黄金值时使用的软件版本为：**R 4.5.3**、**ape 5.8.1**、**phangorn 2.11.1**、
**R.utils 2.13.0**、**FNN 1.1.4.1**、**RColorBrewer 1.1-3**。这些版本是记录值而非
强制要求（见 `reproduce/golden_r3f.R` 的文件头）。

## `example.nwk` 的序列登录号

该树由 274 条线粒体基因组序列构建（dos Reis et al. 2012,
doi:10.1098/rspb.2012.0683）。逐分类群的登录号清单**不属于**本仓库：
`data/examples/example.nwk` 只携带末端标签与枝长，而它随附的 R3F 来源信息指向的是
论文本身，不是机器可读的登录号表，因此原始序列通过该论文引用，本文不在源数据之外
断言任何登录号标识。

## CorrScore 精度说明

`corrtest()` 用 `format(score, digits = 5)` 写出得分，因此
`r3f/example_sr0_corrtest.txt` 记录的是 `score = 0.9996`，
`r3f/example_sr50_corrtest.txt` 记录的是 `0.99959`。这两者是参考结果的 5 位有效数
字，而非结果本身。未舍入的 double 存于
`r3f/example_corrtest_full_precision.txt`：做法是重新 source 同一批 R3F 函数，只把
该打印宽度放宽到 17 位（`reproduce/golden_r3f_fullprec.R`），数值路径未作改动。

## 本目录的可复现性

`Rscript reproduce/golden_r3f.R <project_root> <R3F/R>` 在 micromamba 环境
`r-4.5.3` 中运行（R 4.5.3、ape 5.8.1、phangorn 2.11.1、R.utils 2.13.0、
FNN 1.1.4.1、RColorBrewer 1.1-3）。12 个文件中有 11 个返回**逐字节一致**的结果。
唯一的例外是 `example_RRF_timetree.nexus`，它唯一的差异来自 ape 写进文件头那一行的
时间戳（`[R-package APE, <day date time year>]`）。该文件应视为**数值可复现而逐字
节不可复现**；与它比对时必须忽略文件头的时间戳。第十三个文件
`example_corrtest_full_precision.txt` 完全不由这条命令写出，见上文的
`reproduce/golden_r3f_fullprec.R`。

### 如何取得脚本所需的 R3F 源码

脚本不再默认指向作者私人工作目录里的路径。它按以下顺序解析 R3F，两者都找不到时
就停下来打印克隆命令：

```bash
git clone https://github.com/cathyqqtao/R3F
git -C R3F checkout 6c30d09b650c83f9b34e10387ac57c4ea60a1624   # 文中引用的提交
Rscript reproduce/golden_r3f.R . "$PWD/R3F/R"                  # 或用 R3F_SOURCE_DIR=…
Rscript reproduce/golden_r3f.R --dry-run                       # 只做解析，不写任何文件
```

已安装的 `R3F` 包作为第二选择；随包黄金值是由 source `R/` 下的源文件得到的，因此
那是第一选择。

**R 环境没有被锁定。** `reproduce/` 不附带 `renv.lock`，而 R3F 在上述提交里也只有
`renv/activate.R` 与 `renv/settings.dcf`——同样没有 lock 文件——所以本页开头那份版本
清单就是这些文件的完整记录。在与该环境不同的栈上重新生成，末位数字可能不同，那就
需要为它另写一条溯源记录。
