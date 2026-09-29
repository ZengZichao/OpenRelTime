# ADR-004: Outgroup rooting + pruning semantics

**Language: English** · [中文](ADR-004-tree-operations-zh.md)

**Status**: accepted

**Context**: R3F (and RelTime in MEGA) consume an *ingroup-only*, rooted,
binary tree.  Users supply an unrooted or outgroup-rooted Newick tree plus a
set of outgroup taxa, so OpenRelTime must re-root and prune in a way that is
interchangeable with the reference tooling.

**Decision**: `openreltime.treeio.root_outgroup` replicates the observed
behaviour of `ape::root(..., resolve.root = TRUE)` followed by
`ape::drop.tip` (verified against ape 5.8.1; see `data/PROVENANCE.md` and
the golden file `data/golden/r3f/example_rooted_dropped.nwk`):

* a new root is inserted with the outgroup MRCA as one child; the outgroup
  side carries a zero-length stem, except when the MRCA is a single tip,
  which keeps its own terminal branch;
* edges along the reversed remaining path keep their lengths and unary nodes
  on that path are contracted (accumulating lengths);
* after outgroup removal, leftover unary roots are contracted and the final
  edge above a single-child root is **discarded**;
* a multi-tip outgroup must form a complete clade; violations raise
  `ValueError` (default) or downgrade to a warning with
  `outgroup_check="warn"`;
* the ingroup must retain at least 3 tips.

`openreltime.tree.contract_unary_nodes` implements the shared
contraction/root-edge-discard primitive; `PhyloNode.assign_ids` numbers
nodes in ape cladewise order (tips 1..n, internal n+1.. in pre-order, root
= n+1) so exported node tables and `--anchor-node` values are directly
comparable with R3F/ape/MEGA.

**Consequences**: trees exchanged with ape/R3F produce identical topology
and branch lengths after rooting + pruning; node identifiers mean the same
thing across tools.  The in-house implementation (ADR-001/ADR-002) keeps the
whole operation cache-consistent and recursion-free.
