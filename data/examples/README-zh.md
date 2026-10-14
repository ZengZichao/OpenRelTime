# 示例数据

**语言：中文** · [English](README.md)

本目录中的文件是**示例演算**，供教程与测试套件使用。它们是输入用的 fixture，不是
分析结果，也不是推荐给真实研究使用的设置。

## `example.nwk`

一棵 274-tip 哺乳动物时间树（7,370 个位点），由 R3F（Tao, Sharma, Tamura & Kumar
2025）再分发，原始出处为 dos Reis et al. 2012。它是黄金标准回归测试的回归输入。来
源与许可见 `data/PROVENANCE-zh.md`。测试与教程使用的外群是三个单孔目末端
`Ornithorhynchus_anatinus`、`Zaglossus_bruijni`、`Tachyglossus_aculeatus`。

## `example_calibrations.tsv` 与 `example_taxa.tsv`

**这两个文件只用于演示文件格式。** 它们说明 `key=value` 校正表与「末端 → 类群」分
类信息表如何排布；其中的*数值*不可当作建议的化石约束，也不可当作经过整理的分类系
统。

**`example_calibrations_feasible.tsv` 是配套的可运行版本**：四条彼此自洽的硬边界，
能在文档默认的 `method='bounds'` 下求解（上面的文件做不到，因为它的第二行只有密度）。

具体说，`example_calibrations.tsv` 中的三条边界

| 分类群集合 | 此处使用的边界 |
|---|---|
| `Homo_sapiens|Pan_troglodytes` | 下限 6.0 / 上限 8.5 |
| `Elephas_maximus|Loxodonta_africana` | 指数密度 `offset=60;mean=20` |
| `Dasypus_novemcinctus|Choloepus_didactylus` | 下限 70.0 / 上限 105.0 |

是随意填写的示意值，目的在于把解析器的每一条分支都跑一遍（硬边界、密度、一条很宽
的硬边界）。若照字面理解，它们彼此矛盾：`Elephas`/`Loxodonta` 的密度所暗示的
`Dasypus`/`Choloepus` 节点远比 70-105 这条边界允许的更老，因此求解器必须把谱系速率
放大约 23 倍，才能让这个文件变得可行。这种不可行性是*示意值*的性质，不是文件格式的
性质；任何准备报告的分析，都应当给出真实、各自有据的边界。

**运行这个文件实际会发生什么。** 在文档默认的 `method='bounds'` 下，它会在进入求解器
之前就被拒绝，因为第二行只有密度，解析器会点名这一行。在 `method='effective'` 下，它
是在求解器收敛*之后*被拒绝的。在 274-tip 基准树上按 msz236 协议实测
（`n_effective=10000`、`seed=20260906`），有效边界求解落在：

| 量 | 取值 |
|---|---|
| `time_factor` | 4432.82 |
| 根年龄（节点 272） | 1162.86 Ma |
| 最年轻的节点（节点 307） | **-3998.56 Ma** |
| 父子次序倒置 | 4 处（节点 306 位于 1746.23 Ma，却在其祖先节点 305 的 105 Ma 之下；末端 32 与 33 位于节点 307 之下；节点 372 位于 75.49 Ma，却在其祖先节点 371 的 8.5 Ma 之下） |
| 为 `NaN` 的 `ImpliedRate` 条目 | 4 个 |

负年龄、以及后代比祖先更老，都不是“另一种答案”，而是“没有答案”。
`openreltime.calibrate` 现在会把收敛后的向量对照树的拓扑重新检查，并抛出
`ValueError` 点名相关节点 id，而不再像 0.1.0 那样带着四条日志警告把这些年龄返回给
用户（`calibrate.py` 的违反修复循环每轮只检查当轮最严重的违反者，因此在树的别处刚刚
造成的倒置永远不会被复查）。上表这些数字正是被该检查拒绝的内容，而且它们随种子变化
——这本身就说明成因是结构性的，而非某个种子的偶发问题。请把这个文件当作
**预期会失败的格式示意**，凡是要真正运行的东西，使用下文的
`example_calibrations_feasible.tsv`。

> 这条提示写在本 README 里而不是写在 TSV 的文件头里，原因是：读入器
> （`openreltime.calibrate.parse_calibrations`）用严格的制表符分隔 `key=value` 解析
> 器读取该文件，**没有注释语法**。以 `#` 开头的行会被当作数据行解析，并让整个文件
> 读取失败（它会被读成表头的第一个字段，随后必需的 `min_bound`/`max_bound` 列就缺失
> 了）。因此注释式文件头不是这条提示的安全位置，改用 `data/examples/README-zh.md`
> 承载。

## `example_calibrations_feasible.tsv` —— 随包即可直接运行的那一份

与上面那份只演示格式的文件不同，这个文件是**可用的示例**：四条硬边界、不含密度，
且四条边界可以被同一个全局时间因子同时满足。因此它能在文档默认的
`method='bounds'` 下运行，而 `example_calibrations.tsv` 不能（它的第二行只有密度，
会在进入求解器之前就被拒绝）：

```bash
openreltime calibrate -i data/examples/example.nwk \
    -c data/examples/example_calibrations_feasible.tsv \
    --outgroup "Ornithorhynchus_anatinus,Zaglossus_bruijni,Tachyglossus_aculeatus"
```

```python
cal = openreltime.calibrate(times, openreltime.parse_calibrations(
    "data/examples/example_calibrations_feasible.tsv"))          # method='bounds'
```

两种用法都给出 `time_factor = 371.60`，以及下面这些节点年龄：每条都严格落在自己的
边界内，结果里没有 warning，也没有缺失的隐含速率：

| 分类群集合 | 边界 | 定标年龄 |
|---|---|---|
| `Homo_sapiens\|Pan_troglodytes` | 5.0 – 8.0 | 5.59 |
| `Loxodonta_africana\|Elephas_maximus` | 2.0 – 4.0 | 2.70 |
| `Mus_musculus\|Rattus_norvegicus` | 12.0 – 18.0 | 14.76 |
| `Dasypus_novemcinctus\|Choloepus_didactylus` | 43.0 – 52.0 | 47.50 |

这些边界依然是示意性的软边界，不是整理过的化石定标集；它们是按「彼此自洽」选出的，
而不是按古生物学依据选出的。四个分类群集合里有三个是两尖集，其 MRCA 之下还挂着近缘
类群（`Pan_paniscus`；其余 `Mus`/`Rattus` 物种；其余异关节类），所以读入时会打印
只告警的「不是完整单系群」提示，同时仍然指向正确的 MRCA——这是 MRCA 定标的正常行为，
不是错误。
