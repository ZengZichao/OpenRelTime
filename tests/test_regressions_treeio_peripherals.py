"""Regression gates for tree I/O and the peripheral modules.

Each test pins one behaviour that must hold:

* outgroup removal leaves no "ghost tip" behind: a labelled internal node that
  loses all of its descendants is dropped with them;
* the MEGA-CC bridge returns parsed MEGA results only, and every row carries a
  ``source`` recording where its numbers came from;
* the ddBD fitted density is plotted with the Jacobian of the time rescaling,
  so it is correct in absolute time as well as in shape;
* ``rrf_rates`` states that it applies no rate-ratio guard, matching what the
  rate pass really does;
* ``[&rate=...]`` annotations survive a NEXUS read instead of being stripped,
  and a bracket inside a quoted label is not eaten by comment stripping;
* ``sampling_frac`` follows R3F: rho stays free during the fit and is rewritten
  in the report;
* ddBD start points are accepted under R3F's rule, with its grid order and tie
  breaking;
* BLB forwards the run seed to every IQ-TREE call, an empty grid raises instead
  of indexing out of range, and a single replicate yields no interval rather
  than a fake one;
* ``calibrate`` and ``ci`` restore the averaging convention, guard threshold and
  outgroup check recorded by the run they follow, rather than re-running with
  defaults;
* ``--branch-var`` help matches the header the parser actually requires;
* ``--seed`` reaches ``read_tree`` on every subcommand that can binarise a
  polytomy, so a seeded ``--resolve random`` run reproduces its topology and
  names its seed in the report, while ``--resolve error`` still refuses it;
* the rate colour ramp in ``viz`` is the logarithmic one.
"""

from __future__ import annotations

import importlib

import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from openreltime import _constants as C
from openreltime.cli import cli
from openreltime.megacc import MegaccParseError, parse_megacc_output
from openreltime.tree import invalidate_caches
from openreltime.treeio import (
    _validate,
    node_rates,
    parse_newick,
    read_tree,
    root_outgroup,
    to_newick,
    write_nexus,
)

ROOT = Path(__file__).resolve().parents[1]
GOLDEN = ROOT / "data" / "golden" / "r3f"

SIMPLE = "((A:0.10,B:0.12):0.20,(C:0.15,(D:0.18,E:0.22):0.25):0.30):0.0;"

#: A labelled internal node whose every descendant falls inside the outgroup.
GHOST = (
    "(((A:0.10,X:0.12)OutClade:0.05,(B:0.20,C:0.30)Ingr:0.06):0.04,"
    "(D:0.15,E:0.17)O2:0.03);"
)

#: A megacc-style *time tree*: branch lengths are times, not substitutions.
MEGA_TIMETREE = (
    "(((A:0.0,B:0.0):30.2,(C:0.0,D:0.0):22.4):40.0,(E:0.0,F:0.0):55.0);"
)


def _simple_tree():
    tree = parse_newick(SIMPLE)
    tree.assign_ids()
    return tree


# ---------------------------------------------------------------------------
#  -- ghost tips after outgroup removal
# ---------------------------------------------------------------------------

@pytest.mark.gate("G5")
def test_labelled_outgroup_mrca_does_not_become_a_ghost_tip(tmp_path: Path) -> None:
    path = tmp_path / "ghost.nwk"
    path.write_text(GHOST)
    before = parse_newick(GHOST)
    assert sorted(before.tip_labels()) == ["A", "B", "C", "D", "E", "X"]

    pruned = root_outgroup(before, ["A", "X"])
    # reviewed behaviour: ['B','C','D','E','OutClade'] -- a clade name as a taxon
    assert sorted(pruned.tip_labels()) == ["B", "C", "D", "E"]
    assert pruned.n_tips() == 4
    assert all(not n.internal or n.children for n in pruned.walk())

    whole = read_tree(path, outgroup=["A", "X"])
    assert whole.tip_labels() == ["B", "C", "D", "E"]


@pytest.mark.gate("G4")
def test_validation_rejects_an_internals_degenerated_into_a_tip() -> None:
    """The guard behind : no internal node may survive as a tip."""
    tree = parse_newick(GHOST)
    tree.assign_ids()
    ghost = tree.find_internal_by_label("OutClade")
    assert ghost is not None
    for child in list(ghost.children):  # simulate a pruner that let it survive
        ghost.children.remove(child)
        child.parent = None
    invalidate_caches(tree)
    with pytest.raises(ValueError, match="degenerated into tips"):
        _validate(tree, "error", None)


# ---------------------------------------------------------------------------
#  -- MEGA-CC bridge provenance
# ---------------------------------------------------------------------------

def test_megacc_time_tree_is_read_as_megacc_numbers(tmp_path, monkeypatch) -> None:
    tree_file = tmp_path / "runMEGA-TT.nwk"
    tree_file.write_text(MEGA_TIMETREE + "\n")
    # the bridge must not touch OpenRelTime's estimator at all
    import openreltime.times as times_module

    def forbidden(*_args, **_kwargs):
        raise AssertionError("OpenRelTime's rrf_times must not stand in for MEGA")

    monkeypatch.setattr(times_module, "rrf_times", forbidden)

    frame = parse_megacc_output([tree_file])
    assert list(frame.columns) == ["NodeLabel", "NodeId", "Des1", "Des2", "Time", "source"]
    assert (frame["source"] == "megacc").all()
    by_id = {int(r.NodeId): float(r.Time) for r in frame.itertuples()}
    # node ages are the root-to-present path lengths of MEGA's own tree
    assert by_id[7] == pytest.approx(70.2)  # root: 40.0 + 30.2
    assert by_id[8] == pytest.approx(30.2)  # MRCA(A,B)
    assert by_id[1] == pytest.approx(0.0) and by_id[2] == pytest.approx(0.0)
    assert by_id[3] == pytest.approx(7.8)  # non-ultrametric tip keeps MEGA's depth

    ingroup = parse_megacc_output([tree_file], outgroup=["E", "F"])
    assert sorted(v for v in ingroup["NodeLabel"] if v != "-") == ["A", "B", "C", "D"]


def test_megacc_refuses_unparsable_and_ambiguous_output(tmp_path) -> None:
    log = tmp_path / "megacc.log"
    log.write_text("MEGA-CC (c) 2022; see (the) manual; done;\n")
    with pytest.raises(MegaccParseError, match="no usable megacc output"):
        parse_megacc_output([log])

    other = tmp_path / "secondTT.nwk"
    first = tmp_path / "firstTT.nwk"
    first.write_text(MEGA_TIMETREE)
    other.write_text(MEGA_TIMETREE)
    with pytest.raises(MegaccParseError, match="several megacc outputs"):
        parse_megacc_output([first, other])

    # a tree-like log line is not a tree: it must not be parsed as one, and
    # must not be mistaken for the (unique) table either
    with pytest.raises(MegaccParseError, match="no usable megacc output"):
        parse_megacc_output([log, tmp_path / "missing.nwk"])

    # a table with a parenthesised unit header is a table, not a tree
    table = tmp_path / "run_times.tsv"
    table.write_text("NodeId\tTime (Myr)\n6\t100.0\n7\t55.5\n")
    frame = parse_megacc_output([table])
    assert (frame["source"] == "megacc").all()
    assert [float(v) for v in frame["Time"]] == [100.0, 55.5]


# ---------------------------------------------------------------------------
#  -- ddBD density units on the plotted axis
# ---------------------------------------------------------------------------

def _plotted_ddbd_curve(ddbd_result, times_result, monkeypatch, target: Path):
    """Return ``(x, y)`` of the fitted curve drawn by ``plot_ddbd``."""
    import matplotlib.pyplot as plt

    from openreltime.viz import plot_ddbd

    captured: dict[str, object] = {}
    monkeypatch.setattr(plt, "close", lambda fig=None: captured.__setitem__("fig", fig))
    plot_ddbd(ddbd_result, times_result, target)
    ax = captured["fig"].axes[0]  # type: ignore[union-attr]
    line = ax.get_lines()[-1]
    return np.asarray(line.get_xdata(), dtype=float), np.asarray(
        line.get_ydata(), dtype=float
    )


@pytest.mark.gate("G8")
def test_ddbd_fitted_curve_integrates_to_one_on_the_plotted_axis(
    tmp_path, monkeypatch
) -> None:
    ddbd_module = importlib.import_module("openreltime.ddbd")
    tree = _simple_tree()
    times = importlib.import_module("openreltime.times").rrf_times(tree)
    result = ddbd_module.ddbd(tree, anchor_node=tree.node_id, anchor_time=8.0)
    assert result.scale_factor == pytest.approx(8.0)

    x, y = _plotted_ddbd_curve(result, times, monkeypatch, tmp_path / "bd.png")
    inside = (x >= 0.0) & (x <= result.scale_factor) & np.isfinite(y)
    area = float(np.trapezoid(y[inside], x[inside]))
    # the histogram drawn on the same axis is a density (mass 1), so the
    # fitted curve must have mass 1 too; without the 1/sf Jacobian it had
    # mass == scale_factor (8 here, 500 for a 500 Myr anchor)
    assert area == pytest.approx(1.0, abs=2e-2)


# ---------------------------------------------------------------------------
#  (viz part) -- rate colour ramp
# ---------------------------------------------------------------------------

def test_rate_colour_ramp_is_logarithmic_by_default() -> None:
    viz = importlib.import_module("openreltime.viz")
    rates = {1: 0.03, 2: 1.0, 3: 22.05, 4: 0.0}
    log_ramp = viz._rate_ramp(rates, "log")
    linear_ramp = viz._rate_ramp(rates, "linear")
    assert log_ramp[2] == pytest.approx(5, abs=1)  # the middle rate gets a middle colour
    assert linear_ramp[2] == 0  # ... and collapses onto the extreme under linear
    assert log_ramp[1] == 0 and log_ramp[3] == len(viz._RATE_PALETTE) - 1
    assert log_ramp[4] == 0  # a zero rate must not index the palette from the back


# ---------------------------------------------------------------------------
#  -- [&rate=...] survives a read
# ---------------------------------------------------------------------------

@pytest.mark.gate("G2c")
def test_nexus_rate_annotations_survive_a_read(tmp_path: Path) -> None:
    tree = _simple_tree()
    rates = {n.node_id: round(0.25 * n.node_id, 4) for n in tree.walk()}
    horizon = max(n.node_id for n in tree.walk())
    times = {
        n.node_id: 0.0 if n.is_tip() else (horizon - n.node_id + 1) * 0.05
        for n in tree.walk()
    }
    path = tmp_path / "rates.nexus"
    path.write_text(write_nexus(tree, rate_map=rates, time_map=times))

    back = read_tree(path, fmt="nexus")
    assert node_rates(back) == pytest.approx(rates, rel=1e-9)

    # and the round trip is stable in the write direction, too
    again = write_nexus(back, rate_map=node_rates(back), time_map=times)
    assert again.count("[&rate=") == len(rates)


@pytest.mark.gate("G2c")
def test_r3f_annotated_golden_nexus_keeps_its_rates() -> None:
    """Read direction of gate G2c: ape/R3F ``[&rate=]`` annotations are kept.

    The reviewed reader stripped every comment before parsing, so a written
    NEXUS came back with no rates at all (the gate only proved the write
    direction, and no test read this file).
    """
    tree = read_tree(GOLDEN / "example_RRF_timetree.nexus", fmt="nexus")
    read_rates = node_rates(tree)
    golden = pd.read_csv(GOLDEN / "example_RRF_table.csv")
    assert len(read_rates) == len(golden)
    assert max(
        abs(read_rates[int(nid)] - float(rate))
        for nid, rate in zip(golden["NodeId"], golden["Rate"])
    ) == 0.0


@pytest.mark.gate("G4")
def test_comments_do_not_cut_through_labels(tmp_path: Path) -> None:
    # a bracket inside a quoted label must survive comment stripping
    tree = parse_newick("((A:0.1,'B[X], sp.':0.2)N:0.3,(C:0.1,D:0.2)M:0.15):0.0;")
    assert sorted(tree.tip_labels()) == ["A", "B[X], sp.", "C", "D"]
    assert [t.blen for t in tree.tips()] == [0.1, 0.2, 0.1, 0.2]

    # NHOL allows nested comments; the parser must not stop at the inner ']'
    text = "[&R] ((A:0.1[a [nested] comment],B:0.2[&rate=0.5])N:0.3):0.0;"
    tree = parse_newick(text)
    assert tree.tip_labels() == ["A", "B"]
    assert tree.find_tip("A").rate is None
    assert tree.find_tip("B").rate == pytest.approx(0.5)
    assert tree.internal is True

    path = tmp_path / "plain.nwk"
    path.write_text("((A:0.1,B:0.2)N:0.3[&rate=1.25],(C:0.1,D:0.2)M:0.15):0.0;\n")
    read = read_tree(path)
    assert read.find_internal_by_label("N").rate == pytest.approx(1.25)


# ---------------------------------------------------------------------------
# issues 22 and 23 -- ddBD fit semantics vs R3F
# ---------------------------------------------------------------------------

@pytest.mark.gate("G8")
def test_ddbd_sampling_frac_leaves_rho_free_like_r3f(monkeypatch) -> None:
    """R3F keeps rho free and only rewrites the *reported* sampling fraction."""
    ddbd_module = importlib.import_module("openreltime.ddbd")
    tree = _simple_tree()

    seen_bounds: list = []
    real_minimize = ddbd_module.minimize

    def spy(func, x0, *args, **kwargs):
        seen_bounds.append(list(kwargs["bounds"]))
        return real_minimize(func, x0, *args, **kwargs)

    monkeypatch.setattr(ddbd_module, "minimize", spy)

    free = ddbd_module.ddbd(tree)
    fixed_high = ddbd_module.ddbd(tree, sampling_frac=0.9)
    fixed_low = ddbd_module.ddbd(tree, sampling_frac=0.1)

    # rho is optimised inside (0, 1) in every branch, never pinned to rho0
    assert seen_bounds, "no optimisation was attempted"
    assert all(bounds[2] == (0.0, 1.0) for bounds in seen_bounds), seen_bounds[0]
    assert fixed_high.sampling_frac == pytest.approx(0.9)
    assert fixed_low.sampling_frac == pytest.approx(0.1)
    # ... and the reported birth/death rates are the *same* free-rho optimum
    assert fixed_high.birth_rate == pytest.approx(free.birth_rate, rel=1e-12)
    assert fixed_low.birth_rate == pytest.approx(free.birth_rate, rel=1e-12)
    assert fixed_high.death_rate == pytest.approx(free.death_rate, rel=1e-12)
    # the actually-fitted rho stays inspectable
    assert fixed_high.params["sampling_frac_fitted"] == pytest.approx(
        free.sampling_frac, rel=1e-9
    )


@pytest.mark.gate("G8")
def test_ddbd_grid_follows_expand_grid_order(monkeypatch) -> None:
    """Birth rate varies fastest, and tied scores keep the first grid row."""
    ddbd_module = importlib.import_module("openreltime.ddbd")
    starts: list[list[float]] = []

    def always_fails(func, x0, *args, **kwargs):
        starts.append([float(v) for v in x0])
        raise RuntimeError("boom")

    monkeypatch.setattr(ddbd_module, "minimize", always_fails)
    # a parameter-independent density ties every start score
    monkeypatch.setattr(
        ddbd_module, "_bd_density", lambda t, *a, **k: np.ones_like(np.asarray(t))
    )
    with pytest.raises(RuntimeError, match="best parameter setting"):
        ddbd_module.ddbd(_simple_tree())

    assert len(starts) == C.DDBD_MAX_START_ATTEMPTS
    # expand.grid(Var1=birth, Var2=death, Var3=sampling) varies Var1 fastest,
    # and rows with birth < death are dropped (ddbd.R lines 399-403)
    assert starts[0] == pytest.approx([1.1, 1.0, 0.001])
    assert starts[1] == pytest.approx([2.1, 1.0, 0.001])
    assert starts[9] == pytest.approx([10.1, 1.0, 0.001])
    assert starts[10] == pytest.approx([2.1, 2.0, 0.001])


@pytest.mark.gate("G8")
def test_ddbd_accepts_abnormal_l_bfgs_b_termination(monkeypatch) -> None:
    """ABNORMAL_TERMINATION_IN_LNSRCH coefficients are usable; R3F keeps them."""
    ddbd_module = importlib.import_module("openreltime.ddbd")

    class Result:
        x = np.array([7.0, 0.9, 0.5])
        fun = 1000.0
        success = False
        status = 1
        message = "ABNORMAL_TERMINATION_IN_LNSRCH"

    monkeypatch.setattr(ddbd_module, "minimize", lambda *a, **k: Result())
    result = ddbd_module.ddbd(_simple_tree())
    assert result.birth_rate == pytest.approx(7.0)
    assert result.death_rate == pytest.approx(0.9)
    assert result.sampling_frac == pytest.approx(0.5)

    class Garbage(Result):
        x = np.array([np.nan, 0.9, 0.5])
        fun = np.inf
        success = True
        status = 0

    monkeypatch.setattr(ddbd_module, "minimize", lambda *a, **k: Garbage())
    with pytest.raises(RuntimeError, match="best parameter setting"):
        ddbd_module.ddbd(_simple_tree())


# ---------------------------------------------------------------------------
#  -- BLB reproducibility and degenerate inputs
# ---------------------------------------------------------------------------

BLB_TREES = [
    "(((A:0.10,B:0.12)X:0.20,(C:0.15,(D:0.18,E:0.22)Y:0.25)Z:0.30)W:0.04,F:0.30)R;",
    "(((A:0.11,B:0.13)X:0.21,(C:0.14,(D:0.19,E:0.21)Y:0.24)Z:0.29)W:0.05,F:0.31)R;",
]

#: Stand-in for the IQ-TREE banner ``blb`` records in its JSON report.
BLB_VERSION_BANNER = "iqtree 2.3.1 macOS-13.6.1-arm64 " \
    "(3919f8ae5c 2024-05-05T07:39:31-04:00) [5.0.0]"


@pytest.fixture()
def blb_env(tmp_path, monkeypatch):
    """Stub out IQ-TREE: record argv, write a replicate treefile.

    ``blb`` probes ``--version`` before it runs any replicate so the invoked
    release lands in the JSON report; that probe is recorded separately in
    ``version_calls`` so the argv/``-seed`` assertions below keep describing the
    *replicate* command lines, exactly as they did before the probe existed.
    """
    bootstrap = importlib.import_module("openreltime.bootstrap")
    alignment = tmp_path / "aln.fasta"
    alignment.write_text(
        "".join(f">{name}\n{'ACGT' * 15}\n" for name in "ABCDEF")
    )
    calls: list[list[str]] = []
    version_calls: list[list[str]] = []

    def fake_run_tool(argv, check=True):
        if "--version" in argv:
            version_calls.append(list(argv))

            class VersionDone:
                stdout = BLB_VERSION_BANNER.encode()
                stderr = b""
                returncode = 0

            return VersionDone()
        calls.append(list(argv))
        pre = argv[argv.index("-pre") + 1]
        Path(f"{pre}.treefile").write_text(
            BLB_TREES[len(calls) % len(BLB_TREES)] + "\n"
        )

        class Done:
            stdout = b""
            stderr = b""
            returncode = 0

        return Done()

    monkeypatch.setattr(bootstrap, "run_tool", fake_run_tool)
    blb_env.version_calls = version_calls
    return bootstrap, alignment, calls, tmp_path


def test_blb_seeds_iqtree_and_guards_degenerate_replicates(blb_env, caplog) -> None:
    import logging

    bootstrap, alignment, calls, tmp_path = blb_env

    with caplog.at_level(logging.INFO, logger="openreltime"):
        first = bootstrap.blb(
            alignment, "/bin/echo", n_subsample=1, n_replicate=2, seed=7,
            workdir=tmp_path / "run-a",
        )
    seeds = [c[c.index("-seed") + 1] for c in calls]
    assert len(seeds) == 2 and "-seed" in calls[0]
    assert all(a != b for a, b in zip(seeds, seeds[1:]))

    calls.clear()
    second = bootstrap.blb(
        alignment, "/bin/echo", n_subsample=1, n_replicate=2, seed=7,
        workdir=tmp_path / "run-b",
    )
    assert [c[c.index("-seed") + 1] for c in calls] == seeds  # reproducible
    assert second.equals(first)

    calls.clear()
    bootstrap.blb(
        alignment, "/bin/echo", n_subsample=1, n_replicate=1, seed=7,
        iqtree_args=["-m", "GTR+G", "-seed", "99"], workdir=tmp_path / "run-c",
    )
    assert sum(c.count("-seed") for c in calls) == 1  # the user's -seed wins
    assert calls[0][calls[0].index("-seed") + 1] == "99"

    # an unsupplied seed is still traceable
    calls.clear()
    with caplog.at_level(logging.INFO, logger="openreltime"):
        bootstrap.blb(
            alignment, "/bin/echo", n_subsample=1, n_replicate=1,
            workdir=tmp_path / "run-d",
        )
    assert "master seed" in caplog.text

    # zero replicates: a clean error, not an IndexError on replicate_results[0]
    with pytest.raises(ValueError, match="must be >= 1"):
        bootstrap.blb(
            alignment, "/bin/echo", n_subsample=0, n_replicate=10, seed=1,
            workdir=tmp_path / "run-e",
        )

    # a single replicate: no interval exists -> NaN, and the user is told
    with caplog.at_level(logging.WARNING, logger="openreltime"):
        lone = bootstrap.blb(
            alignment, "/bin/echo", n_subsample=1, n_replicate=1, seed=3,
            workdir=tmp_path / "run-f",
        )
    assert lone["ci_lower"].isna().all() and lone["ci_upper"].isna().all()
    assert (lone["ci_upper"] - lone["ci_lower"]).notna().sum() == 0
    assert "fewer than two replicates" in caplog.text
    assert (lone["median_time"] > 0).all()

    # two replicates do produce a (non-degenerate) interval again
    pair = bootstrap.blb(
        alignment, "/bin/echo", n_subsample=1, n_replicate=2, seed=3,
        workdir=tmp_path / "run-g",
    )
    assert pair["ci_lower"].notna().all()
    assert (pair["ci_upper"] > pair["ci_lower"]).all()


# ---------------------------------------------------------------------------
# issues 25 and 26 -- CLI fidelity
# ---------------------------------------------------------------------------

CALIBRATIONS = "taxon_set\tmin_bound\tmax_bound\nA|B\t1.0\t2.0\n"


@pytest.fixture()
def cli_inputs(tmp_path: Path) -> tuple[Path, Path, list[str]]:
    tree_file = tmp_path / "tree.nwk"
    tree_file.write_text(SIMPLE)
    cal_file = tmp_path / "cals.tsv"
    cal_file.write_text(CALIBRATIONS)
    return tree_file, cal_file, []


@pytest.mark.gate("G3")
def test_cli_calibrate_records_and_ci_restores_estimator_settings(
    cli_inputs, monkeypatch
) -> None:
    tree_file, cal_file, _ = cli_inputs
    times_module = importlib.import_module("openreltime.times")
    prefix = str(tree_file.parent / "cal")
    real_rrf = times_module.rrf_times

    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "calibrate", "-i", str(tree_file), "-c", str(cal_file),
            "--mean", "arithmetic", "--no-guard", "--outgroup-check", "warn",
            "-o", prefix,
        ],
    )
    assert result.exit_code == 0, result.output
    recorded = json.loads(Path(f"{prefix}_report.json").read_text())["calibrate"]
    params = recorded["parameters"]
    assert params["mean"] == "arithmetic"
    assert params["rate_ratio_threshold"] is None
    assert params["outgroup_check"] == "warn"
    assert params["input_fmt"] == "newick"
    assert params["resolve_polytomy"] == "error"

    seen: dict[str, object] = {}

    def spy(tree, **kwargs):
        seen.update(kwargs)
        seen["n_tips"] = tree.n_tips()
        return real_rrf(tree, **kwargs)

    monkeypatch.setattr(times_module, "rrf_times", spy)
    ci_result = runner.invoke(
        cli, ["ci", "-c", prefix, "--n-sites", "1000", "-o", f"{prefix}_ci"]
    )
    assert ci_result.exit_code == 0, ci_result.output
    # the CI stage must re-time the tree with the *recorded* settings
    assert seen["mean"] == "arithmetic"
    assert seen["rate_ratio_threshold"] is None

    # ... and an explicit override wins over the recorded value
    seen.clear()
    override = runner.invoke(
        cli,
        ["ci", "-c", prefix, "--rate-ratio-threshold", "5", "-o", f"{prefix}_ci2"],
    )
    assert override.exit_code == 0, override.output
    assert seen["rate_ratio_threshold"] == 5.0
    assert seen["mean"] == "arithmetic"


@pytest.mark.gate("G3")
def test_cli_errors_are_clean_not_tracebacks(cli_inputs) -> None:
    tree_file, cal_file, _ = cli_inputs
    runner = CliRunner()

    missing = runner.invoke(cli, ["calibrate", "-i", str(tree_file), "-c", "nope.tsv"])
    assert missing.exit_code == 1
    assert "calibration file not found" in missing.output
    assert "Traceback" not in missing.output

    bogus = tree_file.parent / "bogus.tsv"
    bogus.write_text("taxon_set\tmin_bound\nA|B\tnot-a-number\n")
    bad = runner.invoke(
        cli, ["calibrate", "-i", str(tree_file), "-c", str(bogus)]
    )
    assert bad.exit_code == 1 and "cannot read calibrations" in bad.output

    empty = runner.invoke(cli, ["ci", "-c", str(tree_file.parent / "nothing")])
    assert empty.exit_code == 1
    assert "no 'calibrate' block" in empty.output

    both = runner.invoke(
        cli,
        ["calibrate", "-i", str(tree_file), "-c", str(cal_file),
         "--no-guard", "--rate-ratio-threshold", "5"],
    )
    assert both.exit_code == 1 and "mutually exclusive" in both.output


@pytest.mark.gate("G3")
def test_branch_var_help_matches_the_header_the_parser_requires(cli_inputs) -> None:
    tree_file, cal_file, _ = cli_inputs
    option = next(
        p for p in cli.commands["ci"].params if p.name == "branch_var"
    )
    help_text = " ".join(str(option.help or "").split())
    assert "'node_id'" in help_text and "'var'" in help_text
    # the help must not advertise "TSV: node_id<TAB>vS(b)": the parser
    # rejected with a KeyError for exactly such a file
    assert "vS(b)" not in help_text

    prefix = str(tree_file.parent / "cal")
    runner = CliRunner()
    assert runner.invoke(
        cli,
        ["calibrate", "-i", str(tree_file), "-c", str(cal_file), "-o", prefix],
    ).exit_code == 0

    good = tree_file.parent / "v.tsv"
    good.write_text("node_id\tvar\n6\t0.004\n7\t0.002\n")
    ok = runner.invoke(
        cli, ["ci", "-c", prefix, "--branch-var", str(good), "-o", f"{prefix}_ci"]
    )
    assert ok.exit_code == 0, ok.output

    legacy = tree_file.parent / "legacy.tsv"
    legacy.write_text("node_id\tvS(b)\n6\t0.004\n")
    bad = runner.invoke(
        cli, ["ci", "-c", prefix, "--branch-var", str(legacy), "-o", f"{prefix}_x"]
    )
    assert bad.exit_code == 1
    assert "'node_id' and 'var'" in bad.output


@pytest.mark.gate("G3")
def test_rates_help_and_report_admit_the_guard_is_ignored(
    cli_inputs,
) -> None:
    """`rates --rate-ratio-threshold` must not promise a guard it never applies."""
    tree_file, _, _ = cli_inputs
    option = next(p for p in cli.commands["rates"].params if p.name == "rrt")
    help_text = " ".join(str(option.help or "").split())
    # the help must not say "Enable the guard (default off)", since rrf_rates
    # cannot do: the guard rewrites node times, and a rates-only pass has none
    assert "Enable the guard" not in help_text
    assert "node times only" in help_text

    runner = CliRunner()
    bare = str(tree_file.parent / "bare")
    guarded = str(tree_file.parent / "guarded")
    for prefix, extra in ((bare, []), (guarded, ["--rate-ratio-threshold", "10"])):
        res = runner.invoke(
            cli, ["rates", "-i", str(tree_file), "-o", prefix, *extra]
        )
        assert res.exit_code == 0, res.output

    block = json.loads(Path(f"{guarded}_report.json").read_text())["rrf_rates"]
    assert block["parameters"]["rate_ratio_threshold"] == 10
    assert any("ignored by rrf_rates" in w for w in block["warnings"])

    # ... and the threshold really changes nothing in the numbers
    left = pd.read_csv(f"{bare}_rates.csv")
    right = pd.read_csv(f"{guarded}_rates.csv")
    pd.testing.assert_frame_equal(left, right)


# ---------------------------------------------------------------------------
#  -- --seed reaches the polytomy resolver
# ---------------------------------------------------------------------------

#: The trichotomy gate G4 already refuses
#: (``test_calibration_ci_edges.py::test_polytomy_rejected_and_random_resolution``),
#: reused here so the CLI is tested against the same polytomy.
POLYTY = "((A:1,B:2,C:3):4,(D:5,E:6):7);"

#: A seven-way polytomy: it has ~5.7e4 binary resolutions, so a run that ignored
#: its seed could not land on the same topology as the seeded run by chance.
POLYTY_FAN = (
    "((A:1.0,B:2.0,C:3.0,D:4.0,E:5.0,F:6.0,G:7.0):8.0,(H:1.0,I:2.0):9.0);"
)

#: Every subcommand that can binarise a polytomy; ``tree2table`` reads trees but
#: has no ``--resolve``, and ``ci`` inherits both settings from the report.
RESOLVING_SUBCOMMANDS = [
    "rates",
    "times",
    "rates-times",
    "calibrate",
    "corrtest",
    "ddbd",
    "monophyly",
]

#: Subcommand -> the report block that must name the seed it read the tree
#: with.  ``monophyly`` is absent: it writes only ``<prefix>_monophyly.csv``.
REPORTED_SEED_CASES = [
    ("rates", "rrf_rates"),
    ("times", "rrf_times"),
    ("rates-times", "rrf_times"),
    ("calibrate", "calibrate"),
    ("corrtest", "corrtest"),
    ("ddbd", "ddbd"),
]

#: Subcommand -> the artefact whose bytes encode the resolved topology.
REPLAYED_TOPOLOGY_CASES = [
    ("rates", "rates.nwk"),
    ("times", "timetree.nwk"),
    ("rates-times", "timetree.nexus"),
]

#: Seed a trichotomy has to honour: with it the resolution is fixed, without it
#: the same command draws a different pairing on every run.
SEEDED = ["--resolve", "random", "--seed", "7"]


@pytest.fixture()
def poly_file(tmp_path: Path) -> Path:
    path = tmp_path / "poly.nwk"
    path.write_text(POLYTY)
    return path


@pytest.fixture()
def fan_file(tmp_path: Path) -> Path:
    path = tmp_path / "fan.nwk"
    path.write_text(POLYTY_FAN)
    return path


@pytest.fixture()
def read_tree_calls(monkeypatch) -> list[dict]:
    """Record the keyword arguments of every ``read_tree`` call the CLI makes."""
    treeio = importlib.import_module("openreltime.treeio")
    calls: list[dict] = []
    real = treeio.read_tree

    def spy(*args, **kwargs):
        calls.append(kwargs)
        return real(*args, **kwargs)

    monkeypatch.setattr(treeio, "read_tree", spy)
    return calls


def _argv(sub: str, poly: Path, prefix: Path, extra: list[str]) -> list[str]:
    argv = [sub, "-i", str(poly), *extra, "-o", str(prefix)]
    if sub == "calibrate":
        cal = poly.parent / "calibrations.tsv"
        cal.write_text(CALIBRATIONS)
        argv[3:3] = ["-c", str(cal)]
    return argv


def _artefacts(sub: str, poly: Path, prefix: Path, seed: int) -> dict[str, bytes]:
    res = CliRunner().invoke(
        cli, _argv(sub, poly, prefix, ["--resolve", "random", "--seed", str(seed)])
    )
    assert res.exit_code == 0, res.output
    return {
        p.name.removeprefix(f"{prefix.name}_"): p.read_bytes()
        for p in poly.parent.glob(f"{prefix.name}_*")
    }


@pytest.mark.gate("G3")
@pytest.mark.parametrize("sub", RESOLVING_SUBCOMMANDS)
def test_cli_seed_reaches_read_tree(poly_file, read_tree_calls, sub) -> None:
    prefix = poly_file.parent / f"seeded_{sub}"
    res = CliRunner().invoke(cli, _argv(sub, poly_file, prefix, SEEDED))
    assert res.exit_code == 0, res.output
    assert len(read_tree_calls) == 1, read_tree_calls
    assert read_tree_calls[0]["resolve_polytomy"] == "random"
    assert read_tree_calls[0]["seed"] == 7


@pytest.mark.gate("G3")
@pytest.mark.parametrize("sub,analysis", REPORTED_SEED_CASES)
def test_cli_records_the_tree_seed_in_the_report(poly_file, sub, analysis) -> None:
    prefix = poly_file.parent / f"reported_{sub}"
    res = CliRunner().invoke(cli, _argv(sub, poly_file, prefix, SEEDED))
    assert res.exit_code == 0, res.output
    report = json.loads(Path(f"{prefix}_report.json").read_text())
    assert report[analysis]["parameters"]["seed"] == 7
    assert report[analysis]["parameters"]["resolve_polytomy"] == "random"


@pytest.mark.gate("G3")
@pytest.mark.parametrize("sub,topology_artefact", REPLAYED_TOPOLOGY_CASES)
def test_cli_random_resolution_is_seed_reproducible(
    fan_file, sub, topology_artefact
) -> None:
    """Same seed -> byte-identical artefacts; another seed -> another topology."""
    first = _artefacts(sub, fan_file, fan_file.parent / "seed7a", 7)
    second = _artefacts(sub, fan_file, fan_file.parent / "seed7b", 7)
    assert first == second
    other = _artefacts(sub, fan_file, fan_file.parent / "seed5", 5)
    assert other[topology_artefact] != first[topology_artefact]


@pytest.mark.gate("G3")
def test_cli_ci_restores_the_recorded_polytomy_seed(poly_file, read_tree_calls) -> None:
    """``ci`` re-reads with calibrate's seed, so it qualifies the same topology."""
    prefix = poly_file.parent / "cal"
    res = CliRunner().invoke(cli, _argv("calibrate", poly_file, prefix, SEEDED))
    assert res.exit_code == 0, res.output
    read_tree_calls.clear()
    ci_prefix = poly_file.parent / "ci"
    ci = CliRunner().invoke(cli, ["ci", "-c", str(prefix), "-o", str(ci_prefix)])
    assert ci.exit_code == 0, ci.output
    assert len(read_tree_calls) == 1, read_tree_calls
    assert read_tree_calls[0]["resolve_polytomy"] == "random"
    assert read_tree_calls[0]["seed"] == 7
    report = json.loads(Path(f"{ci_prefix}_report.json").read_text())
    assert report["confidence_interval"]["parameters"]["seed"] == 7
    assert report["confidence_interval"]["parameters"]["resolve_polytomy"] == "random"


@pytest.mark.gate("G4")
@pytest.mark.parametrize("sub", RESOLVING_SUBCOMMANDS)
def test_cli_resolve_error_still_refuses_the_polytomy(poly_file, sub) -> None:
    """The seed plumbing must not have loosened the binary-tree requirement."""
    prefix = poly_file.parent / f"strict_{sub}"
    res = CliRunner().invoke(
        cli, _argv(sub, poly_file, prefix, ["--resolve", "error", "--seed", "7"])
    )
    assert res.exit_code != 0
    assert "polytom" in str(res.exception)
    assert not list(poly_file.parent.glob(f"{prefix.name}_*"))


# ---------------------------------------------------------------------------
#  -- FigTree readability of the emitted timetree (gate G2c)
# ---------------------------------------------------------------------------

@pytest.mark.gate("G2c")
def test_emitted_nexus_satisfies_the_figtree_readability_contract(
    tmp_path: Path,
) -> None:
    """Structural contract behind "FigTree opens OpenRelTime timetrees".

    FigTree is not importable here, so this asserts the properties a FigTree
    (and ape's ``read.nexus``) parser depends on rather than launching the GUI:

    * a ``#NEXUS`` header and balanced ``BEGIN``/``END`` blocks;
    * the rooted-tree token ``[&R]`` on the ``TREE`` line, without which
      FigTree treats the topology as unrooted;
    * a TAXA block whose labels are exactly the tip labels and which carries no
      node annotations, because annotations there are read as taxon comment
      fields and shift the tip mapping;
    * one ``[&rate=...]`` comment per node — keyed by label on the tips and
      attached positionally after the closing parenthesis on internal clades,
      which is the Newick-Hudson convention ape and R3F emit and FigTree's
      comment parser consumes;
    * numeric branch lengths throughout.

    The engine's own bracket-aware reader is used as the independent parser: a
    round trip that recovers every rate shows the comments are well formed and
    correctly placed, not merely present as text.
    """
    tree = _simple_tree()
    rates = {n.node_id: round(0.25 * n.node_id, 4) for n in tree.walk()}
    horizon = max(n.node_id for n in tree.walk())
    times = {
        n.node_id: 0.0 if n.is_tip() else (horizon - n.node_id + 1) * 0.05
        for n in tree.walk()
    }
    nexus = write_nexus(tree, rate_map=rates, time_map=times)

    # -- header and blocks
    assert nexus.lstrip().startswith("#NEXUS")
    assert nexus.count("BEGIN ") == nexus.count("END;") == 2
    taxa_block = nexus.split("BEGIN TAXA;")[1].split("END;")[0]
    trees_block = nexus.split("BEGIN TREES;")[1]

    # -- the rooted token FigTree needs
    tree_line = [ln for ln in trees_block.splitlines() if ln.strip().startswith("TREE")][0]
    assert "[&R]" in tree_line, tree_line

    # -- TAXA block: exactly the tips, and annotation-free
    declared = taxa_block.split("TAXLABELS")[1].replace(";", "").split()
    assert sorted(declared) == sorted(tree.tip_labels())
    assert len(declared) == tree.n_tips()
    assert "[&" not in taxa_block, "annotations in TAXA would shift FigTree's mapping"

    # -- one rate comment per node: tips by label, internal clades positionally
    assert trees_block.count("[&rate=") == len(rates)
    body = tree_line.split("[&R]", 1)[1]
    tip_annots = re.findall(r"([A-Za-z0-9_.|-]+)\[&rate=([-\d.eE+]+)\]", body)
    internal_annots = re.findall(r"\)\[&rate=([-\d.eE+]+)\]", body)
    assert sorted(lbl for lbl, _ in tip_annots) == sorted(tree.tip_labels())
    assert len(internal_annots) == tree.n_internal()
    assert len(tip_annots) + len(internal_annots) == len(rates)
    for _, value in tip_annots:
        assert float(value) > 0
    for value in internal_annots:
        assert float(value) > 0

    # -- the tip annotations really carry the rates we asked for
    tip_rates = {n.node_id: rates[n.node_id] for n in tree.walk() if n.is_tip()}
    assert {float(v) for _, v in tip_annots} == {float(v) for v in tip_rates.values()}

    # -- branch lengths are numeric everywhere, including the root edge
    for length in re.findall(r":([-\d.eE+]+)", body):
        assert math.isfinite(float(length))

    # -- an independent bracket-aware parse recovers every rate
    out = tmp_path / "figtree.nexus"
    out.write_text(nexus)
    back = read_tree(out, fmt="nexus")
    assert node_rates(back) == pytest.approx(rates, rel=1e-9)


@pytest.mark.gate("G2c")
def test_plain_newick_output_drops_rate_annotations(tmp_path: Path) -> None:
    """The documented limit: only ``.nexus`` keeps ``[&rate=...]``.

    ``to_newick()`` has no ``rate_map`` parameter, so a ``.nwk`` artefact is not
    FigTree-annotatable.  Asserting the limitation keeps it from silently
    becoming a claim, and points users at the ``.nexus`` artefact instead.
    """
    tree = _simple_tree()
    rates = {n.node_id: 1.0 for n in tree.walk()}
    assert "[&rate=" in write_nexus(tree, rate_map=rates)

    path = tmp_path / "plain.nwk"
    path.write_text(to_newick(tree))
    assert "[&rate=" not in path.read_text()
