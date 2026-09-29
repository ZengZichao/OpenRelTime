# ADR 索引

**语言：中文** · [English](README.md)

| ADR | 标题 | 状态 |
|---|---|---|
| ADR-001 | 自研的二叉树结构（`openreltime.tree.PhyloNode`），不用 Bio.Phylo/DendroPy | accepted |
| ADR-002 | 自研 Newick/NEXUS 读写；DendroPy 仅作为可选的互操作 extra | accepted |
| ADR-003 | 许可：GPL-3.0-or-later（与 R3F 和 MEGA 源码一致） | accepted |
| ADR-004 | 外群置根与剔除的语义复刻 `ape::root(resolve.root=TRUE)` 加 `ape::drop.tip`，包括丢弃一元根上累加而成的根边（全文见 [ADR-004-tree-operations-zh.md](ADR-004-tree-operations-zh.md)） | accepted |
| ADR-005 | CI 偏导数用**整条估计管线**的中心差分求得（在协方差为零的假设下与 msz236 的闭式解加递推代数等价） | accepted |

ADR-004 与 ADR-005 的完整正文写在 `openreltime/treeio.py` 与 `openreltime/ci.py`
的 docstring 中，并见 `docs/methods-zh.md`。
