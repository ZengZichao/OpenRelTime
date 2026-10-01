"""Command line interface for OpenRelTime (``openreltime`` / ``ort``).

Every subcommand mirrors a public API function and writes the same files
(``<prefix>_*.csv/nwk/nexus/txt/json``); stderr carries structured warnings.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional

import click

from . import _constants as C

logger = logging.getLogger("openreltime")

_VERBOSE = False


def _setup_logging(progress: bool = False) -> None:
    """Configure stderr logging.

    ``-v/--verbose`` (top-level flag) enables DEBUG; subcommands that stream
    progress from long-running pipelines enable INFO via *progress*; the
    default level is WARNING.
    """
    if _VERBOSE:
        level = logging.DEBUG
    elif progress:
        level = logging.INFO
    else:
        level = logging.WARNING
    logging.basicConfig(
        level=level,
        format="%(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )


def _read_outgroup(outgroup: str | None) -> list[str] | None:
    if not outgroup:
        return None
    path = Path(outgroup)
    if path.exists():
        return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
    return [tok.strip() for tok in outgroup.split(",") if tok.strip()]


def _read_taxon_table(taxon_table: str | None) -> dict[str, list[str]] | None:
    if not taxon_table:
        return None
    from .taxonomy import parse_taxon_table

    return parse_taxon_table(taxon_table)


def _tree_input(
    tree: str,
    outgroup: str | None,
    fmt: str,
    resolve: str,
    taxon_table: str | None = None,
    outgroup_check: str = "error",
):
    from .treeio import read_tree

    return read_tree(
        tree,
        fmt=fmt,
        outgroup=_read_outgroup(outgroup),
        resolve_polytomy=resolve,
        taxon_table=_read_taxon_table(taxon_table),
        outgroup_check=outgroup_check,
    )


def _outgroup_options(command):
    """Attach the shared taxonomy/outgroup options to *command*."""
    command = click.option(
        "--taxon-table",
        "taxon_table",
        default=None,
        type=click.Path(),
        help="TSV/CSV tip->group table; enables group-name outgroups/calibrations.",
    )(command)
    command = click.option(
        "--outgroup-check",
        "outgroup_check",
        default="error",
        type=click.Choice(["error", "warn"]),
        show_default=True,
        help="Multi-tip outgroup that is not a complete clade: error or warn.",
    )(command)
    return command


#: Options that decide *which* RRF estimate a downstream stage reproduces
#: (``--mean`` and the rate-ratio guard); ``times`` documents the same pair.
def _rrf_estimate_options(command):
    """Attach the shared RRF estimator options to *command*."""
    command = click.option(
        "--no-guard",
        is_flag=True,
        help=f"Disable the rate-ratio guard (default threshold "
        f"{C.RATE_RATIO_THRESHOLD}).",
    )(command)
    command = click.option(
        "--rate-ratio-threshold",
        "rrt",
        default=None,
        type=float,
        help=f"Extreme-rate guard threshold (default: "
        f"{C.RATE_RATIO_THRESHOLD}; ignored with --no-guard).",
    )(command)
    command = click.option(
        "--mean",
        default="geometric",
        type=click.Choice(["geometric", "arithmetic"]),
        help="Averaging convention (geometric = R3F-compatible).",
    )(command)
    return command


def _guard_threshold(no_guard: bool, rrt: Optional[float]) -> Optional[float]:
    """Resolve ``--no-guard`` / ``--rate-ratio-threshold`` into a threshold."""
    if no_guard:
        if rrt is not None:
            raise click.ClickException(
                "--no-guard and --rate-ratio-threshold are mutually exclusive"
            )
        return None
    return C.RATE_RATIO_THRESHOLD if rrt is None else rrt


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(package_name="openreltime", prog_name="openreltime")
@click.option("-v", "--verbose", is_flag=True, help="Log progress to stderr.")
def cli(verbose: bool) -> None:
    """OpenRelTime: relative rate framework molecular dating in Python."""
    global _VERBOSE
    _VERBOSE = verbose


@cli.command()
@click.option("-i", "--input", "tree", required=True, help="Branch-length tree (Newick/NEXUS).")
@click.option("--outgroup", default=None, help="Comma-separated tips/group names or a one-per-line file.")
@click.option("--fmt", "fmt", default="newick", type=click.Choice(["newick", "nexus"]))
@click.option("--resolve", "resolve", default="error", type=click.Choice(["error", "random"]))
@click.option("--mean", default="geometric", type=click.Choice(["geometric", "arithmetic"]))
@click.option("--rate-ratio-threshold", "rrt", default=None, type=float,
              help="No effect here: the guard rewrites node times only, so the "
              "rate pass ignores this value and records a warning. Set it on "
              "'times' or 'rates-times' instead.")
@click.option("-o", "--output", "prefix", default="ORT", show_default=True)
@_outgroup_options
def rates(tree, outgroup, fmt, resolve, mean, rrt, prefix, taxon_table, outgroup_check):
    """Estimate relative lineage rates (rrf_rates)."""
    _setup_logging()
    from .rates import rrf_rates

    t = _tree_input(tree, outgroup, fmt, resolve, taxon_table, outgroup_check)
    result = rrf_rates(t, mean=mean, rate_ratio_threshold=rrt)
    for p in result.write(prefix):
        click.echo(str(p))


@cli.command()
@click.option("-i", "--input", "tree", required=True)
@click.option("--outgroup", default=None)
@click.option("--fmt", "fmt", default="newick", type=click.Choice(["newick", "nexus"]))
@click.option("--resolve", "resolve", default="error", type=click.Choice(["error", "random"]))
@click.option("--normalize", is_flag=True, help="Scale the root time to 1.0.")
@click.option("--r3f-compat", is_flag=True, help="Write R3F file/column names.")
@click.option("-o", "--output", "prefix", default="ORT", show_default=True)
@_outgroup_options
@_rrf_estimate_options
def times(tree, outgroup, fmt, resolve, mean, normalize, no_guard, rrt, r3f_compat, prefix, taxon_table, outgroup_check):
    """Estimate relative node times (rrf_times)."""
    _setup_logging()
    from .times import rrf_times

    t = _tree_input(tree, outgroup, fmt, resolve, taxon_table, outgroup_check)
    result = rrf_times(
        t,
        mean=mean,
        rate_ratio_threshold=_guard_threshold(no_guard, rrt),
        normalize=normalize,
    )
    for p in result.write(prefix, with_rate=False, r3f_compat=r3f_compat):
        click.echo(str(p))


@cli.command(name="rates-times")
@click.option("-i", "--input", "tree", required=True)
@click.option("--outgroup", default=None)
@click.option("--fmt", "fmt", default="newick", type=click.Choice(["newick", "nexus"]))
@click.option("--resolve", "resolve", default="error", type=click.Choice(["error", "random"]))
@click.option("--mean", default="geometric", type=click.Choice(["geometric", "arithmetic"]),
              help="Averaging convention (default: R3F-compatible geometric).")
@click.option("--normalize", is_flag=True)
@click.option("--plot", "plot", default=None, type=click.Path(), help="Optional rate-coloured timetree image path.")
@click.option("--r3f-compat", is_flag=True)
@click.option("-o", "--output", "prefix", default="ORT", show_default=True)
@_outgroup_options
def rates_times(tree, outgroup, fmt, resolve, mean, normalize, plot, r3f_compat, prefix, taxon_table, outgroup_check):
    """Estimate relative rates and node times (rrf_rates_times)."""
    _setup_logging()
    from .times import rrf_rates_times

    t = _tree_input(tree, outgroup, fmt, resolve, taxon_table, outgroup_check)
    result = rrf_rates_times(t, mean=mean, normalize=normalize)
    paths = result.write(prefix, with_rate=True, nexus=True, r3f_compat=r3f_compat)
    if plot:
        from .viz import plot_timetree

        paths.append(plot_timetree(result, plot))
    for p in paths:
        click.echo(str(p))


@cli.command(name="calibrate")
@click.option("-i", "--input", "tree", required=True)
@click.option("-c", "--calibrations", required=True, help="calibrations.tsv (see docs).")
@click.option("--outgroup", default=None)
@click.option("--fmt", "fmt", default="newick", type=click.Choice(["newick", "nexus"]))
@click.option("--resolve", "resolve", default="error", type=click.Choice(["error", "random"]))
@click.option("--method", default="bounds", type=click.Choice(["bounds", "effective"]))
@click.option("--n-effective", "n_effective", default=10000, show_default=True)
@click.option("--seed", default=None, type=int)
@click.option("-o", "--output", "prefix", default="ORT", show_default=True)
@_outgroup_options
@_rrf_estimate_options
def calibrate_cmd(tree, calibrations, outgroup, fmt, resolve, method, n_effective, seed, prefix, taxon_table, outgroup_check, mean, no_guard, rrt):
    """Convert relative times to absolute times with calibrations."""
    _setup_logging()
    from .calibrate import calibrate as run_calibrate
    from .times import rrf_times

    table = _read_taxon_table(taxon_table)
    t = _tree_input(tree, outgroup, fmt, resolve, taxon_table, outgroup_check)
    # the estimator settings are passed to rrf_times, so the mean/threshold
    # that calibrate records from times.params are the ones actually used
    times = rrf_times(
        t, mean=mean, rate_ratio_threshold=_guard_threshold(no_guard, rrt)
    )
    cals = _parse_calibrations_or_fail(calibrations)
    result = run_calibrate(
        times,
        cals,
        method=method,
        n_effective=n_effective,
        seed=seed,
        tree_file=str(Path(tree).resolve()),
        calibrations_file=str(Path(calibrations).resolve()),
        outgroup=_read_outgroup(outgroup),
        taxon_table_file=str(Path(taxon_table).resolve()) if taxon_table else None,
        taxon_table=table,
        rate_ratio_threshold=times.params.get("rate_ratio_threshold"),
        input_fmt=fmt,
        resolve_polytomy=resolve,
        outgroup_check=outgroup_check,
    )
    for p in result.write(prefix):
        click.echo(str(p))


def _parse_calibrations_or_fail(path: str):
    """``parse_calibrations`` with user-facing errors instead of tracebacks."""
    from .calibrate import parse_calibrations

    if not path or not Path(path).exists():
        raise click.ClickException(
            f"calibration file not found: {path!r} (expected the TSV described "
            "in docs/usage-en.md, section 4)"
        )
    try:
        cals = parse_calibrations(path)
    except (OSError, ValueError, KeyError) as exc:
        raise click.ClickException(f"cannot read calibrations from {path}: {exc}")
    if not cals:
        raise click.ClickException(f"no calibrations found in {path}")
    return cals


@cli.command(name="ci")
@click.option("-c", "--calibrated", required=True, help="Prefix of a previous calibrate run.")
@click.option(
    "--branch-var",
    default=None,
    type=click.Path(),
    help="TSV with a header row and exactly the columns 'node_id' and 'var' "
    "(the per-edge sampling variance); see docs.",
)
@click.option("--n-sites", "n_sites", default=None, type=int)
@click.option("--level", default=0.95, show_default=True)
@click.option("-o", "--output", "prefix", default="ORT", show_default=True)
@click.option(
    "--mean",
    default=None,
    type=click.Choice(["geometric", "arithmetic"]),
    help="Override the averaging convention recorded by calibrate.",
)
@click.option(
    "--rate-ratio-threshold",
    "rrt",
    default=None,
    type=float,
    help="Override the extreme-rate guard threshold recorded by calibrate.",
)
@click.option(
    "--no-guard",
    is_flag=True,
    help="Disable the rate-ratio guard (default threshold "
    f"{C.RATE_RATIO_THRESHOLD}).",
)
def ci(calibrated, branch_var, n_sites, level, prefix, mean, no_guard, rrt):
    """Analytical confidence intervals for a calibrated result.

    The tree is re-read and re-timed with the settings recorded in
    ``<calibrated>_report.json`` (outgroup, format, polytomy handling, mean,
    rate-ratio guard), so the intervals belong to the estimate that is being
    qualified.  ``--mean``/``--rate-ratio-threshold``/``--no-guard`` override
    the recorded values when given explicitly.
    """
    _setup_logging()
    import pandas as pd

    from .calibrate import calibrate as run_calibrate
    from .ci import confidence_interval
    from .report import read_report_block
    from .times import rrf_times
    from .treeio import read_tree

    report = read_report_block(calibrated, "calibrate")
    if not report:
        raise click.ClickException(
            f"no 'calibrate' block in {calibrated}_report.json: pass the prefix "
            "of a previous 'openreltime calibrate' run"
        )
    params = report.get("parameters", {})
    calibs_path = params.get("calibrations_file")
    tree_path = params.get("tree_file")
    outgroup = params.get("outgroup")
    taxon_table_file = params.get("taxon_table_file")
    fmt = params.get("input_fmt", "newick")
    resolve = params.get("resolve_polytomy", "error")
    # 'warn' must not degrade into 'error' at the CI stage
    outgroup_check = params.get("outgroup_check", "error")
    if tree_path is None or not Path(tree_path).exists():
        raise click.ClickException(
            "cannot locate the original tree; re-run the calibrate step with "
            "the same -i/--output so tree_file is recorded in the report"
        )
    if not calibs_path:
        raise click.ClickException(
            "the report does not record a calibrations_file; re-run 'openreltime "
            "calibrate -c <calibrations.tsv> ...' so the calibration table is "
            "reproducible"
        )
    table = _read_taxon_table(taxon_table_file)
    t = read_tree(
        tree_path,
        fmt=fmt,
        outgroup=outgroup,
        resolve_polytomy=resolve,
        taxon_table=table,
        outgroup_check=outgroup_check,
    )
    recorded_mean = params.get("mean") or "geometric"
    recorded_threshold = params.get("rate_ratio_threshold", C.RATE_RATIO_THRESHOLD)
    if no_guard and rrt is not None:
        raise click.ClickException(
            "--no-guard and --rate-ratio-threshold are mutually exclusive"
        )
    threshold = None if no_guard else (rrt if rrt is not None else recorded_threshold)
    times = rrf_times(t, mean=mean or recorded_mean, rate_ratio_threshold=threshold)
    cals = _parse_calibrations_or_fail(calibs_path)
    method = params.get("method", "bounds")
    cal = run_calibrate(
        times,
        cals,
        method=method,
        n_effective=params.get("n_effective") or 10000,
        seed=params.get("seed"),
        taxon_table=table,
    )
    branch_var_map = None
    if branch_var:
        try:
            frame = pd.read_csv(branch_var, sep="\t")
            branch_var_map = {
                int(r["node_id"]): float(r["var"]) for _, r in frame.iterrows()
            }
        except (OSError, ValueError, KeyError) as exc:
            raise click.ClickException(
                f"cannot read --branch-var {branch_var!r}: {exc}; expected a TSV "
                "with a header row and the columns 'node_id' and 'var'"
            )
    result = confidence_interval(
        cal,
        branch_var=branch_var_map,
        n_sites=n_sites,
        level=level,
    )
    for p in result.write(prefix):
        click.echo(str(p))


@cli.command()
@click.option("-i", "--input", "tree", required=True)
@click.option("--outgroup", default=None)
@click.option("--fmt", "fmt", default="newick", type=click.Choice(["newick", "nexus"]))
@click.option("--resolve", "resolve", default="error", type=click.Choice(["error", "random"]))
@click.option("--sister-resample", "sister_resample", default=0, type=int, help="Use >50 when the tree has <50 tips.")
@click.option("--seed", default=None, type=int)
@click.option("--anchor-node", "anchor_node", default=None, type=int)
@click.option("--anchor-time", "anchor_time", default=0.0, type=float)
@click.option("-o", "--output", "prefix", default="ORT", show_default=True)
@_outgroup_options
def corrtest(tree, outgroup, fmt, resolve, sister_resample, seed, anchor_node, anchor_time, prefix, taxon_table, outgroup_check):
    """CorrTest: test for autocorrelated evolutionary rates."""
    _setup_logging()
    from .corrtest import corrtest as run_corrtest

    t = _tree_input(tree, outgroup, fmt, resolve, taxon_table, outgroup_check)
    result = run_corrtest(
        t,
        sister_resample=sister_resample,
        seed=seed,
        anchor_node=anchor_node,
        anchor_time=anchor_time,
    )
    for p in result.write(prefix):
        click.echo(str(p))


@cli.command()
@click.option("-i", "--input", "tree", required=True)
@click.option("--outgroup", default=None)
@click.option("--fmt", "fmt", default="newick", type=click.Choice(["newick", "nexus"]))
@click.option("--resolve", "resolve", default="error", type=click.Choice(["error", "random"]))
@click.option("--sampling-frac", "sampling_frac", default=None, type=float)
@click.option("--anchor-node", "anchor_node", default=None, type=int)
@click.option("--anchor-time", "anchor_time", default=1.0, type=float)
@click.option("--measure", default="SSE", type=click.Choice(["SSE", "KL"]))
@click.option("--plot", "plot", default=None, type=click.Path())
@click.option("-o", "--output", "prefix", default="ORT", show_default=True)
@_outgroup_options
def ddbd(tree, outgroup, fmt, resolve, sampling_frac, anchor_node, anchor_time, measure, plot, prefix, taxon_table, outgroup_check):
    """ddBD: birth-death speciation tree prior from RRF times."""
    _setup_logging()
    from .ddbd import ddbd as run_ddbd
    from .times import rrf_times

    t = _tree_input(tree, outgroup, fmt, resolve, taxon_table, outgroup_check)
    result = run_ddbd(
        t,
        sampling_frac=sampling_frac,
        anchor_node=anchor_node,
        anchor_time=anchor_time,
        measure=measure,
    )
    paths = result.write(prefix)
    if plot:
        from .viz import plot_ddbd

        paths.append(plot_ddbd(result, rrf_times(t), plot))
    for p in paths:
        click.echo(str(p))


@cli.command(name="tree2table")
@click.option("-i", "--input", "tree", required=True)
@click.option("--fmt", "fmt", default="newick", type=click.Choice(["newick", "nexus"]))
@click.option("--time", "with_time", is_flag=True, help="Report node ages instead of branch lengths.")
@click.option("-o", "--output", "prefix", default="ORT", show_default=True)
def tree2table_cmd(tree, fmt, with_time, prefix):
    """Convert a tree to a node table (tree2table)."""
    _setup_logging()
    from .table import tree2table

    t = _tree_input(tree, None, fmt, "error")
    frame = tree2table(t, time=with_time)
    out = Path(f"{prefix}.csv")
    frame.to_csv(out, index=False)
    click.echo(str(out))


@cli.command()
@click.option("-i", "--input", "tree", required=True, help="Tree file whose tips are tested.")
@click.option("--fmt", "fmt", default="newick", type=click.Choice(["newick", "nexus"]))
@click.option("--taxon-table", "taxon_table", default=None, type=click.Path(),
              help="TSV/CSV tip->group classification table (default: auto-detect groups from tip-label prefixes).")
@click.option("--groups", default=None, help="Comma-separated group names to test (default: all groups).")
@click.option("--tips", default=None, help="Comma-separated tip names forming one ad-hoc set to test.")
@click.option("--group-sep", "group_sep", default=None,
              help="Regex separator for auto-detecting groups from tip labels (default: first of '_', '|', '@' or whitespace).")
@click.option("--resolve", "resolve", default="error", type=click.Choice(["error", "random"]))
@click.option("-o", "--output", "prefix", default="ORT", show_default=True)
def monophyly(tree, fmt, taxon_table, groups, tips, group_sep, resolve, prefix):
    """Test whether taxon groups are monophyletic (complete clades).

    Groups come from --taxon-table or are auto-detected from the tip labels
    (``Genus_species`` -> genus ``Genus``).  Each tested group is reported
    with its MRCA, any extra tips under the MRCA and a monophyly verdict;
    non-monophyletic groups additionally raise a warning on stderr.
    """
    import pandas as pd

    from .taxonomy import (
        detect_groups_from_labels,
        parse_taxon_table,
        warn_if_not_monophyletic,
    )
    from .treeio import read_tree

    _setup_logging()
    t = read_tree(tree, fmt=fmt, resolve_polytomy=resolve)
    if taxon_table:
        registry = parse_taxon_table(taxon_table)
        source = f"taxon table {taxon_table}"
    else:
        registry = detect_groups_from_labels(t.tip_labels(), sep=group_sep or r"[_|@\s]")
        source = "auto-detected tip-label prefixes"
    logger.info("monophyly groups sourced from %s (%d groups)", source, len(registry))

    targets: list[tuple[str, list[str]]] = []
    if groups:
        for name in [g.strip() for g in groups.split(",") if g.strip()]:
            if name not in registry:
                raise click.ClickException(
                    f"group {name!r} not found in {source}; available: "
                    f"{', '.join(sorted(registry))}"
                )
            targets.append((name, list(registry[name])))
    if tips:
        targets.append(("custom-tip-set", [s.strip() for s in tips.split(",") if s.strip()]))
    if not targets:
        targets = [(name, list(members)) for name, members in sorted(registry.items())]

    rows = []
    for name, members in targets:
        res = warn_if_not_monophyletic(t, members, name=name)
        rows.append(
            {
                "group": res.name,
                "n_tips": len(res.tips),
                "is_monophyletic": res.is_monophyletic,
                "n_mrca_tips": len(res.mrca_tips),
                "extra_tips": ";".join(res.extra_tips),
                "missing_names": ";".join(res.missing),
                "tips": ";".join(sorted(res.tips)),
            }
        )
    frame = pd.DataFrame(rows)
    out = Path(f"{prefix}_monophyly.csv")
    frame.to_csv(out, index=False)
    click.echo(str(out))


@cli.command()
@click.option("-a", "--alignment", required=True, help="FASTA alignment.")
@click.option("--iqtree", required=True, help="Path to the IQ-TREE executable.")
@click.option("--outgroup", default=None)
@click.option("--gamma", default=0.7, show_default=True)
@click.option("--n-subsample", "n_subsample", default=10, show_default=True)
@click.option("--n-replicate", "n_replicate", default=10, show_default=True)
@click.option("--iqtree-arg", "iqtree_args", multiple=True, help="Extra IQ-TREE args, e.g. --iqtree-arg=-m --iqtree-arg=GTR+G")
@click.option("--seed", default=None, type=int)
@click.option("--resolve", "resolve_polytomy", default="random", show_default=True, type=click.Choice(["error", "random"]))
@click.option("--workdir", default=None, type=click.Path())
@click.option("-o", "--output-dir", "output_dir", default="ORT_blb", show_default=True, type=click.Path())
def blb(alignment, iqtree, outgroup, gamma, n_subsample, n_replicate, iqtree_args, seed, resolve_polytomy, workdir, output_dir):
    """BLB pipeline: subsample + bootstrap alignment replicates through IQ-TREE."""
    _setup_logging(progress=True)
    from .bootstrap import blb as run_blb

    extra = []
    for item in iqtree_args:
        extra.extend(item.split())
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    # By default keep replicate files in the output directory so both the
    # summary and the per-replicate table are retained and reproducible.
    wd = Path(workdir) if workdir else out
    table = run_blb(
        alignment,
        iqtree,
        outgroup=_read_outgroup(outgroup),
        gamma=gamma,
        n_subsample=n_subsample,
        n_replicate=n_replicate,
        iqtree_args=extra,
        seed=seed,
        resolve_polytomy=resolve_polytomy,
        workdir=str(wd),
    )
    table.to_csv(out / "blb_summary.csv", index=False)
    # run_blb writes the per-replicate table into the workdir; copy it to the
    # user-facing output directory alongside the summary when the two differ.
    rep_src = wd / "blb_replicates.csv"
    rep_dst = out / "blb_replicates.csv"
    if rep_src.exists() and rep_src.resolve() != rep_dst.resolve():
        import shutil
        shutil.copy2(rep_src, rep_dst)
    click.echo(str(out / "blb_summary.csv"))


@cli.command()
@click.option("-i", "--input", "tree", required=True)
@click.option("-c", "--calibrations", default=None)
@click.option("--megacc", required=True, help="Path to the megacc executable.")
@click.option("--outgroup", default=None)
@click.option("--mao", default=None, help="Custom .mao file (default: generated).")
@click.option("--workdir", default=None, type=click.Path())
@click.option("-o", "--output", "prefix", default="ORT", show_default=True)
def megacc(tree, calibrations, megacc, outgroup, mao, workdir, prefix):
    """MEGA-CC bridge (optional extra): run megacc RelTime and parse its output."""
    _setup_logging(progress=True)
    from .megacc import parse_megacc_output, run_megacc

    og = _read_outgroup(outgroup)
    _mao, outputs = run_megacc(
        tree,
        megacc=megacc,
        calibrations=calibrations,
        outgroup=og,
        mao=mao,
        workdir=workdir,
    )
    frame = parse_megacc_output(outputs, outgroup=og)
    out = Path(f"{prefix}_megacc_times.csv")
    frame.to_csv(out, index=False)
    click.echo(str(out))


def main() -> None:  # pragma: no cover
    cli()


if __name__ == "__main__":  # pragma: no cover
    main()
