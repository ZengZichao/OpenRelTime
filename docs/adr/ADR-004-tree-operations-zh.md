# ADR-004：外群置根与剔除的语义

**语言：中文** · [English](ADR-004-tree-operations.md)

**状态**：accepted（已采纳）

**背景**：R3F（以及 MEGA 中的 RelTime）消费的是**仅含内群**、有根且二叉的树。用
户提供的是一棵未置根或以外群置根的 Newick 树，外加一组外群分类群，因此 OpenRelTime
必须以可与参考工具互换的方式重新置根并剔除末端。

**决策**：`openreltime.treeio.root_outgroup` 复刻
`ape::root(..., resolve.root = TRUE)` 紧接着 `ape::drop.tip` 的可观测行为（已对照
ape 5.8.1 验证；见 `data/PROVENANCE-zh.md` 与黄金文件
`data/golden/r3f/example_rooted_dropped.nwk`）：

* 插入一个新的根，以外群的 MRCA 作为它的一个子节点；外群一侧带一条零长干枝，除非
  该 MRCA 恰为单个末端，此时保留它自己的末端枝；
* 反向路径上剩余的各边保留原有枝长，该路径上的一元节点被收缩（枝长累加）；
* 删除外群之后，残留的一元根被收缩，单子根之上的最后一条边被**丢弃**；
* 多末端外群必须构成一个完整的簇；违反时默认抛出 `ValueError`，取
  `outgroup_check="warn"` 时降级为警告；
* 内群至少保留 3 个末端。

`openreltime.tree.contract_unary_nodes` 实现共用的「节点收缩＋根边丢弃」原语；
`PhyloNode.assign_ids` 按 ape cladewise 顺序编号（末端 1..n，内部节点自 n+1 起按前
序，根 = n+1），因此导出的节点表与 `--anchor-node` 取值可与 R3F、ape、MEGA 直接对
照。

**后果**：与 ape/R3F 交换的树在置根加剔除之后拓扑与枝长完全一致；节点标识符在各工
具间含义相同。自建实现（ADR-001/ADR-002）让整个操作保持缓存一致且不使用递归。
