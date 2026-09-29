# 第三方声明

**语言：中文** · [English](THIRD-PARTY-NOTICES.md)

OpenRelTime 采用 **GNU General Public License v3.0 或更高版本**（见 `LICENSE`）。
本文件记录参考材料的来源，以及从 GPL 兼容来源改造而来的代码或常数，以满足
GPL-3 §4/§5 的要求，同时落实项目的许可决策 ADR-003（索引见
`docs/adr/README-zh.md`）。

## 1. R3F（GPL-3）——行为参考

* 仓库：<https://github.com/cathyqqtao/R3F>（Qiqing Tao、Shreya Sharma、Koichiro
  Tamura 与 Sudhir Kumar）
* 查阅的文件：`R/rrf_rates.R`、`R/rrf_times.R`、`R/rrf_rates_times.R`、
  `R/corrtest.R`、`R/ddbd.R`、`R/tree2table.R`。
* 复用性质：
  * RRF 遍历的后序/前序结构、零枝长 epsilon 守卫的放置位置、末端孙辈速率的回填、
    前序调整中祖辈速率的特例处理，以及速率比保护，都是按**语义**复刻的（ADR-003 中
    的 Level L2/L3），以便 OpenRelTime 在数值上复现 R3F（门槛 G2a 与 G2b，见
    `docs/validation-zh.md`）；相应代码内嵌了形如
    `R3F/R/rrf_times.R lines 114-131` 的行内引用；
  * CorrTest 的归一化常数、逻辑回归系数与 P 值档位取自 `corrtest.R` 第 555-584 行
    （事实与参数，引用到 Tao et al. 2019）；
  * ddBD 的初值网格、出生-死亡密度与拟合流程遵循 `ddbd.R` 第 399-524 行。
* 改动：所有复用的逻辑都以不同的树数据结构在 Python 中重新实现；不逐字分发任何 R
  代码。
* 版权：(c) Qiqing Tao, Sudhir Kumar。GPL-3.0。

## 2. MEGA 源码（GPL-3）——设计参考

* 仓库：<https://github.com/KumarMEGALab/MEGA-source-code>
* 查阅的文件：`MEGA12.1-source/reltime/mreltimecomputer.pas`
  （`TimeFactor`/`MinTimeFactor`/`MaxTimeFactor` 可行性设计与
  `PropagateConstraints`）、`mcalibrationsampler.pas`（有效边界复制循环）。
* 复用性质：仅为设计层面（L2）。`openreltime/calibrate.py` 中的校正求解器围绕同一
  个全局因子概念，是对已发表描述（msz236）的独立实现。未转录任何 Pascal 代码。
* 版权：(c) Sudhir Kumar, Koichiro Tamura, Glen Stecher 及 MEGA 开发团队。GPL-3.0。

## 3. ape（GPL-2+）——树操作的行为参考

* 外群置根加剔除的流程（`openreltime/treeio.py::root_outgroup`、
  `openreltime/tree.py::contract_unary_nodes`）复现 `ape::root(resolve.root = TRUE)`
  加 `ape::drop.tip`（ape 5.8.1）的可观测行为，使输出的树可互换。相关记录见
  `docs/adr/ADR-004-tree-operations.md`。
* 版权：(c) Emmanuel Paradis, Klaus Schliep。GPL-2+。

## 4. 第三方数据

随包数据的逐文件来源与许可见 `data/PROVENANCE-zh.md`（其中 example.nwk © dos Reis
et al. 2012，随 R3F 以 GPL-3 再分发）。

## 5. 未复用的部分

CorrTest 的参考仓库（`cathyqqtao/CorrTest`）未提供许可证；**本项目未复制它的任何
代码**，只使用 Tao et al. (2019) 已发表的常数。

## 6. 学术致谢

若使用 OpenRelTime，请引用第 1～3 节所列论文，并引用 OpenRelTime 本身（见
`CITATION.cff`），以回应原作者的善意请求。
