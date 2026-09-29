# ADR index

**Language: English** · [中文](README-zh.md)

| ADR | Title | Status |
|---|---|---|
| ADR-001 | Self-contained binary tree structure (`openreltime.tree.PhyloNode`) instead of Bio.Phylo/DendroPy | accepted |
| ADR-002 | In-house Newick/NEXUS I/O; DendroPy only as optional interoperability extra | accepted |
| ADR-003 | Licence: GPL-3.0-or-later (same as R3F and MEGA source) | accepted |
| ADR-004 | Outgroup rooting + pruning semantics replicate `ape::root(resolve.root=TRUE)` + `ape::drop.tip`, including discarding the accumulated root edge of a unary root (full text: [ADR-004-tree-operations.md](ADR-004-tree-operations.md)) | accepted |
| ADR-005 | CI partial derivatives by central finite differences of the whole estimation pipeline (algebraically equivalent to the msz236 closed forms + recursion under the no-covariance assumption) | accepted |

Full texts for ADR-004 and ADR-005 are embedded in the docstrings of
`openreltime/treeio.py` / `openreltime/ci.py` and in `docs/methods.md`.
