# Methods and limitations

**Language: English** · [中文](methods-zh.md)

## 1. What OpenRelTime computes

OpenRelTime implements the analytical solution of the **relative rate
framework (RRF)** of Tamura, Tao & Kumar (2018, MBE 35:1770-1782).  Given a
rooted binary tree whose branch lengths are expected substitutions per site,
it estimates (i) a relative rate for every lineage and (ii) a relative time
for every node, using the *minimum rate change principle*: the rate of an
ancestral lineage is the (geometric, by default) average of the rates of its
two daughter lineages.  Internal nodes are visited postorder (deepest
first); at each node the local 3- or 4-lineage problem is solved in closed
form, subtrees are folded into two composite lengths
(`l = sqrt(l_a * l_b) + l_stem`), and a preorder pass then multiplies local
rates by ancestral rates (dividing times), anchoring the ingroup root at
average rate 1.

Two averaging conventions are available: the **geometric mean** (msy044
eqs 28-42; used by R3F and the default here) and the **arithmetic mean**
(msy044 eqs 1-27; the 2012 paper's semantics, kept for cross-validation).

## 2. Relation to MEGA's RelTime

MEGA's RelTime applies the RRF analytical solution *and* an additional
likelihood-ratio-test (LRT) refinement in which daughter branches whose
rates do not differ significantly are constrained to the ancestral rate.
That step requires the sequence alignment and maximum-likelihood branch
re-optimisation.  **OpenRelTime 0.1.0 implements the analytical RRF layer
only** — the same layer that R3F implements — so results are expected to be
*highly correlated* with MEGA RelTime but not identical.  Report
correlations, not equality, when comparing against MEGA.

## 3. Calibrations

With calibration constraints, a global time factor `f` is chosen so that
`f * t_node` lies within every `[min, max]` interval (feasible-intersection
midpoint).  When constraints cannot be jointly satisfied, the lineage
responsible is adjusted: its rate is scaled so the node age reaches the
violated bound exactly, the scaling is propagated to all descendants, and
the procedure repeats (hard-bound semantics of Tamura et al. 2013 as
described by msz236), with a cap of 100 iterations after which conflicting
constraints are reported.

The **effective bounds** method (msz236) converts probability densities into
data-derived bounds: per replicate, two dates are sampled from each density
and used as a (min, max) pair; the point estimate is kept for calibrated
nodes; 10,000 replicates give the 2.5/97.5 percentiles that become the
effective bounds of the final analysis.

## 4. Confidence intervals

CIs follow the delta method of msz236: total branch variance
`v(b) = vS(b) + vR(b)`; the rate-heterogeneity variance `vR` is recovered by
subtracting the sampling contribution `SV(R)` from the observed rate
variance `Vobs(R)` (eqs 9-14) and distributed over branches proportional to
`b^2`; node-time variance is `sum_j (dt_i/db_j)^2 v(b_j)` (eq 7).  Partial
derivatives are evaluated by central finite differences of the whole
estimation pipeline (equivalent to the paper's closed forms + recursion
under the no-covariance assumption; ADR-005, indexed in
`docs/adr/README.md`).  Intervals use final
(calibrated) rates and are truncated at hard calibration bounds.

Two consistency rules keep the published equations from being mixed:

* *Per-lineage scale.*  eq (9) averages over lineages while eqs (11)-(13) sum
  `SV(R)` over them, so `RV(R) = Vobs(R) - SV(R)/N` is evaluated on the
  per-lineage scale of eq (9); subtracting the raw sum clamps `RV(R)` to zero
  on any real tree and drops the heterogeneity component entirely.
* *One estimator throughout.*  The differentiation pipeline reproduces the
  averaging convention (`mean`), the rate-ratio guard and the calibration
  bounds **actually used** for the point estimate, so an arithmetic-mean or
  `method="effective"` result is qualified by derivatives of that same
  estimator (for effective results the msz236 effective bounds act as the
  constraints).  Analytical intervals are therefore available for either
  calibration method.  A calibrated node whose age is a fixed point of its own
  constraint has an exactly zero eq-(7) derivative, and is instead reported
  with the uncertainty its calibration asserts, flagged `se_reliable = False`
  in `<prefix>_ci.csv`.

The effective-bounds method additionally reports its own uncertainty as
empirical 2.5/97.5 node-age quantiles in `<prefix>_effective_bounds.csv`; both
are available and answer slightly different questions (see
`docs/validation.md` and tutorial 2).

## 5. Assumptions and known limitations

* **Binary, rooted trees.**  Polytomies must be resolved first (optionally
  at random with zero-length branches, which weakens affected estimates).
* **Zero-length branches.**  Rates and times touching many zero branches are
  unreliable; OpenRelTime warns when >= 10% of branches are zero (R3F
  threshold) and substitutes tiny epsilons exactly where R3F does, so both
  programs degrade identically.
* **No LRT refinement** (see section 2) and no likelihood computation;
  branch lengths are taken as given.
* **Outgroup rates are not comparable.**  The outgroup is used for rooting
  and removed; the ingroup/outgroup equal-rate assumption (Kumar et al.
  2016) is not testable within RRF.
* **Internal-rate collapse.**  RRF's minimum-change approximation behaves
  close to a strict clock on internal branches; see Lozano-Fernandez et al.
  (2017, GBE) for an empirical discussion of when this matters.
* **CorrTest** is a fixed-coefficient logistic model trained by Tao et al.
  (2019) on simulated mammal-like trees; for small trees (< 50 tips) use
  sister-pair resampling and interpret the score cautiously.
* **ddBD** assumes species-level trees (no population samples) and a
  birth-death process with constant rates over the tree.
* **No warning when the outgroup is omitted.**  `--outgroup` is optional in
  the CLI, and an unrooted input is accepted as read; its root position then
  silently determines every age (R3F requires an outgroup).  See
  `docs/validation.md` §5.

## 6. Validation (test suite)

| Gate | Test | Result |
|---|---|---|
| G1 | hand-computed 3-/4-taxon values, both conventions, tol 1e-12 | pass |
| G2a | 274-tip mammal tree vs R3F rates: slope 1.0000, median rel. dev. 6.8e-16 | pass |
| G2b | same tree vs R3F times: slope 1.0000, max rel. dev. 3.6e-14 | pass |
| G2c | NEXUS `[&rate=]` round-trip | pass |
| G6 | CorrTest fixed-vector behaviour, bit-for-bit agreement with R3F's unrounded score 0.99959505940763882 (R3F's own golden file stores 5 digits, 0.9996), both P < 0.001 | pass |
| G7 | calibration hitting bounds, joint constraints, conflict detection, effective bounds | pass |
| G7+ | regression gates: bound-less rows rejected, `RV(R)` positive/saturating, pinned-node intervals, CI follows the calibrated estimator | pass |

Run `python -m pytest tests/` to reproduce.
