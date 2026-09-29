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

> 这条提示写在本 README 里而不是写在 TSV 的文件头里，原因是：读入器
> （`openreltime.calibrate.parse_calibrations`）用严格的制表符分隔 `key=value` 解析
> 器读取该文件，**没有注释语法**。以 `#` 开头的行会被当作数据行解析，并让整个文件
> 读取失败（它会被读成表头的第一个字段，随后必需的 `min_bound`/`max_bound` 列就缺失
> 了）。因此注释式文件头不是这条提示的安全位置，改用 `data/examples/README-zh.md`
> 承载。
