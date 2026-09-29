"""BLB pipeline: bag-of-little-bootstraps for time-tree confidence intervals.

Implements the two-level resampling scheme popularised by RelTime-JA
(Barba-Montoya et al. 2023, Front Bioinform), verified against its
``lb_sampler.R``:

* level 1: ``s`` without-replacement subsamples of ``ceil(L^gamma)`` sites
  (``gamma`` in [0.6, 0.9], default 0.7);
* level 2: ``r`` with-replacement resamples of the original length ``L``
  from each subsample;
* each replicate tree is inferred by an external tool (IQ-TREE), timed by
  :func:`openreltime.rrf_times` and aggregated per node clade (median and
  2.5/97.5 percentiles).  Nodes absent from some replicates are reported
  with their ``topology_ok`` fraction.
* the random component is two-fold: the resampling scheme is seeded by
  ``seed``, and IQ-TREE's own heuristics (ML search, subtree pruning and
  the bootstrap trees it draws internally) are seeded per replicate through
  ``-seed`` derived from the *same* generator.  Both are therefore
  reproducible from ``seed`` alone, as ``docs/usage-en.md`` promises.

External-process policy: the executable is resolved and validated up front;
every token of the command is checked against a strict allow-list (no shell
metacharacters) and the process is started from an argument *list* with
``shell=False``.  No command string is ever assembled.
"""

from __future__ import annotations

import logging
import math
import tempfile
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from . import _constants as C
from ._external import resolve_executable, run_tool
from .report import TimeResult
from .times import rrf_times
from .treeio import read_tree

logger = logging.getLogger("openreltime")

__all__ = ["blb"]

#: Upper bound (exclusive) of the 32-bit IQ-TREE seed window, so that the
#: derived values are accepted by every IQ-TREE version.  ``1`` keeps the
#: documented "seeds are positive integers" convention.
#: Source: RelTime-JA lb_sampler.R seeds are 32-bit integers; the range is
#: IQ-TREE's own ``-seed`` convention (a positive integer < 2**31).
_IQTREE_SEED_MAX: int = 2**31 - 1


def _read_fasta(path: Path) -> tuple[list[str], list[str]]:
    names: list[str] = []
    seqs: list[str] = []
    chunks: list[str] = []
    current: Optional[str] = None
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    seqs.append("".join(chunks))
                current = line[1:].split()[0]
                names.append(current)
                chunks = []
            else:
                chunks.append(line)
    if current is not None:
        seqs.append("".join(chunks))
    return names, seqs


def _write_resampled_alignment(
    names: list[str], seqs: list[str], weights: np.ndarray, destination: Path
) -> None:
    """Write an alignment whose column j is source column ``weights[j]``."""
    with destination.open("w") as out:
        for name, seq in zip(names, seqs):
            resampled = "".join(seq[int(i)] for i in weights)
            out.write(f">{name}\n")
            for start in range(0, len(resampled), 60):
                out.write(resampled[start : start + 60] + "\n")


def blb(
    alignment: str | Path,
    iqtree: str,
    *,
    outgroup: Optional[list[str]] = None,
    gamma: float = C.BLB_DEFAULT_GAMMA,
    n_subsample: int = 10,
    n_replicate: int = 10,
    iqtree_args: Optional[list[str]] = None,
    seed: Optional[int] = None,
    workdir: Optional[str | Path] = None,
    resolve_polytomy: str = "random",
    **params: Any,
) -> pd.DataFrame:
    """Run the BLB pipeline and return the per-node aggregate table.

    Parameters
    ----------
    alignment:
        FASTA alignment.
    iqtree:
        Path to (or name on PATH of) the IQ-TREE executable.
    outgroup:
        Outgroup tip names used for every OpenRelTime timing call.
    gamma:
        Subsampling exponent; subsample size is ``ceil(L**gamma)``.  The
        default 0.7 sits inside the admissible range [0.6, 0.9] of the
        RelTime-JA sampler (``lb_sampler.R`` line 4); note that the published
        RelTime-JA protocol itself used ``g = 0.78`` with 20 subsamples x 20
        replicates (``lbs_codelines.r`` line 14), so this default reproduces
        the *method*, not that run.
    n_subsample, n_replicate:
        Level-1 subsample count and level-2 bootstrap count per subsample.
        Both must be >= 1: with no replicate tree there is nothing to
        aggregate, and with a single replicate per clade no quantile interval
        exists (``ci_lower``/``ci_upper`` are then reported as NaN rather
        than as a zero-width interval).
    iqtree_args:
        Extra IQ-TREE arguments forwarded verbatim (e.g. ``["-m", "GTR+G"]``).
        A user-supplied ``-seed`` wins over the derived one.
    seed:
        Seed for the whole resampling scheme, including the per-replicate
        IQ-TREE ``-seed`` values.  ``None`` draws a fresh seed, which is then
        logged so the run stays reproducible.
    workdir:
        Directory keeping replicate files (default: a fresh temp directory).
    resolve_polytomy:
        How to handle non-binary replicate trees: ``"error"`` (fail) or
        ``"random"`` (resolve with zero-length branches). The BLB pipeline
        defaults to ``"random"`` because bootstrap trees frequently contain
        polytomies.

    Returns
    -------
    pandas.DataFrame
        Columns ``node_id, clade_size, n_replicates, topology_ok,
        median_time, ci_lower, ci_upper``.  Also writes
        ``blb_summary.csv`` and ``blb_replicates.csv`` into *workdir*.
    """
    alignment = Path(alignment)
    if not alignment.exists():
        raise FileNotFoundError(f"alignment not found: {alignment}")
    executable = resolve_executable(iqtree)
    if not (C.BLB_GAMMA_RANGE[0] - 1e-9 <= gamma <= C.BLB_GAMMA_RANGE[1] + 1e-9):
        raise ValueError(f"gamma must be within {C.BLB_GAMMA_RANGE}, got {gamma}")
    if int(n_subsample) < 1 or int(n_replicate) < 1:
        raise ValueError(
            "n_subsample and n_replicate must be >= 1 "
            f"(got n_subsample={n_subsample}, n_replicate={n_replicate})"
        )
    extra_args = [str(a) for a in (iqtree_args or [])]
    user_seed = any(a in ("-seed", "--seed") or a.startswith("-seed=") for a in extra_args)
    if seed is None:
        # every source of randomness must be traceable to a printed number
        seed = int(np.random.default_rng().integers(1, _IQTREE_SEED_MAX))
        logger.info(
            "BLB master seed was not supplied; using %d -- rerun with seed=%d "
            "to reproduce these replicate trees",
            seed,
            seed,
        )
    rng = np.random.default_rng(seed)

    names, seqs = _read_fasta(alignment)
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError("alignment sequences have unequal lengths")
    length_sites = lengths.pop()
    sub_size = int(math.ceil(length_sites**gamma))
    out_dir = Path(workdir) if workdir else Path(tempfile.mkdtemp(prefix="ort_blb_"))
    out_dir.mkdir(parents=True, exist_ok=True)

    replicate_results: list[TimeResult] = []
    total = n_subsample * n_replicate
    done = 0
    for s in range(n_subsample):
        subsample_idx = np.sort(rng.choice(length_sites, size=sub_size, replace=False))
        for r in range(n_replicate):
            weights = subsample_idx[rng.integers(0, sub_size, size=length_sites)]
            rep_dir = out_dir / f"rep{s * n_replicate + r:04d}"
            rep_dir.mkdir(exist_ok=True)
            rep_aln = rep_dir / "replicate.fasta"
            _write_resampled_alignment(names, seqs, weights, rep_aln)
            cmd = [
                executable,
                "-s",
                str(rep_aln),
                "-pre",
                str(rep_dir / "iqtree"),
                "-nt",
                "AUTO",
                *extra_args,
            ]
            if not user_seed:
                # ML search and the internal bootstrap trees are stochastic;
                # without -seed the same master seed would still give
                # different replicate trees
                cmd += ["-seed", str(int(rng.integers(1, _IQTREE_SEED_MAX)))]
            run_tool(cmd)
            tree = read_tree(
                rep_dir / "iqtree.treefile",
                outgroup=outgroup,
                resolve_polytomy=resolve_polytomy,
                seed=seed,
            )
            replicate_results.append(rrf_times(tree))
            done += 1
            logger.info("BLB replicate %d/%d done", done, total)

    if not replicate_results:  # pragma: no cover - guarded above, defensive
        raise RuntimeError(
            "no BLB replicate tree was inferred; check the IQ-TREE output and "
            "the n_subsample/n_replicate settings"
        )

    # -- aggregate by clade (robust to replicate topology differences) ----------
    reference = replicate_results[0]
    ref_clades = {
        n.node_id: n.clade for n in reference.tree.walk() if not n.is_tip()
    }
    agg_rows = []
    thin_clades: list[int] = []
    for nid, clade in ref_clades.items():
        values: list[float] = []
        present = 0
        for res in replicate_results:
            match = next(
                (n for n in res.tree.walk() if not n.is_tip() and n.clade == clade),
                None,
            )
            if match is None:
                continue
            present += 1
            values.append(res.times[match.node_id])
        if not values:
            continue
        if len(values) < 2:
            # a 2.5/97.5 percentile of one observation is undefined; reporting
            # the point itself on both sides would fake a zero-width interval
            thin_clades.append(nid)
            ci_lower = ci_upper = float("nan")
        else:
            ci_lower = float(np.quantile(values, 0.025))
            ci_upper = float(np.quantile(values, 0.975))
        agg_rows.append(
            {
                "node_id": nid,
                "clade_size": len(clade),
                "n_replicates": present,
                "topology_ok": present / len(replicate_results),
                "median_time": float(np.median(values)),
                "ci_lower": ci_lower,
                "ci_upper": ci_upper,
            }
        )
    if thin_clades:
        logger.warning(
            "no confidence interval can be formed for %d clade(s) recovered in "
            "fewer than two replicates (node ids %s): ci_lower/ci_upper are NaN",
            len(thin_clades),
            ", ".join(str(nid) for nid in thin_clades[:10])
            + (" ..." if len(thin_clades) > 10 else ""),
        )
    table = pd.DataFrame(agg_rows)
    table.to_csv(out_dir / "blb_summary.csv", index=False)
    rep_frames = []
    for i, res in enumerate(replicate_results):
        frame = res.to_pandas()
        frame.insert(0, "replicate", i)
        rep_frames.append(frame)
    pd.concat(rep_frames, ignore_index=True).to_csv(
        out_dir / "blb_replicates.csv", index=False
    )
    logger.info("BLB summary written to %s", out_dir / "blb_summary.csv")
    return table
