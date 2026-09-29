# Third-party notices

**Language: English** · [中文](THIRD-PARTY-NOTICES-zh.md)

OpenRelTime is licensed under the **GNU General Public License v3.0 or
later** (`LICENSE`).  This file records the provenance of reference material
and any code or constants adapted from GPL-compatible sources, as required by
GPL-3 §4/§5 and by the project licence decision ADR-003 (indexed in
`docs/adr/README.md`).

## 1. R3F (GPL-3) — behavioural reference

* Repository: <https://github.com/cathyqqtao/R3F> (Qiqing Tao, Shreya Sharma, Koichiro Tamura & Sudhir Kumar)
* Files consulted: `R/rrf_rates.R`, `R/rrf_times.R`, `R/rrf_rates_times.R`,
  `R/corrtest.R`, `R/ddbd.R`, `R/tree2table.R`.
* Nature of reuse:
  * the postorder/preorder structure of the RRF traversal, the placement of
    zero-length epsilon guards, the tip-grandchild rate backfill, the
    grandparent-rate special case in the preorder adjustment and the
    rate-ratio guard were replicated **semantically** (Level L2/L3 in
    ADR-003) so that OpenRelTime reproduces R3F numerically (Gate G2a/G2b,
    see `docs/validation.md`);
    corresponding code carries inline citations such as
    `R3F/R/rrf_times.R lines 114-131`;
  * CorrTest normalisation constants, logistic coefficients and P-value
    bands were extracted from `corrtest.R` lines 555-584 (facts/parameters,
    cited to Tao et al. 2019);
  * the ddBD initial-value grid, birth-death density and fitting procedure
    follow `ddbd.R` lines 399-524.
* Modifications: all reused logic was re-implemented in Python with a
  different tree data structure; no R code is distributed verbatim.
* Copyright: (c) Qiqing Tao, Sudhir Kumar.  GPL-3.0.

## 2. MEGA source code (GPL-3) — design reference

* Repository: <https://github.com/KumarMEGALab/MEGA-source-code>
* Files consulted: `MEGA12.1-source/reltime/mreltimecomputer.pas`
  (`TimeFactor`/`MinTimeFactor`/`MaxTimeFactor` feasibility design,
  `PropagateConstraints`), `mcalibrationsampler.pas` (effective-bound
  replicate loop).
* Nature of reuse: design-level (L2) only; the calibration solver in
  `openreltime/calibrate.py` is an independent implementation of the
  published description (msz236) organised around the same global-factor
  concept.  No Pascal code was transcribed.
* Copyright: (c) Sudhir Kumar, Koichiro Tamura, Glen Stecher, and the MEGA
  development team.  GPL-3.0.

## 3. ape (GPL-2+) — behavioural reference for tree operations

* The outgroup rooting + pruning procedure (`openreltime/treeio.py::root_outgroup`,
  `openreltime/tree.py::contract_unary_nodes`) reproduces the observed
  behaviour of `ape::root(resolve.root = TRUE)` + `ape::drop.tip` (ape 5.8.1)
  so that output trees are interchangeable.  Documented in
  `docs/adr/ADR-004-tree-operations.md`.
* Copyright: (c) Emmanuel Paradis, Klaus Schliep.  GPL-2+.

## 4. Third-party data

See `data/PROVENANCE.md` for per-file origin and licence of bundled data
(example.nwk © dos Reis et al. 2012, redistributed with R3F under GPL-3).

## 5. Not reused

The CorrTest reference repository (`cathyqqtao/CorrTest`) ships without a
licence; **no code from it was copied**.  Only the published constants of
Tao et al. (2019) are used.

## 6. Academic disclosure

If you use OpenRelTime, please cite the papers in section 1-3 as well as the
OpenRelTime release (see `CITATION.cff`), consistent with the goodwill
requests of the original authors.
