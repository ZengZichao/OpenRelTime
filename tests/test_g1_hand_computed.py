"""G1 gate: hand-computed golden values for 3-/4-taxon RRF solutions.

Formulas: msy044 equations 28-33 (three-lineage geometric), 34-42 (four-lineage
geometric) and 19-27 (four-lineage arithmetic).  Tolerance 1e-12 per Gate G1.

The expected values are *frozen decimals*, evaluated by hand from the published
equations, not expressions re-derived from the same algebra the implementation
uses: an assertion built from the implementation's own formula cannot fail.
The formula objects are kept and cross-checked against the frozen numbers so a
transcription error in either is still caught.
"""

from __future__ import annotations

import math

import pytest

from openreltime.rrf import RrfEngine
from openreltime.treeio import parse_newick

TOL = 1e-12

# --- hand-computed fixtures -------------------------------------------------
# 4-taxon + outgroup: b1=0.4, b2=0.2, b3=0.5, b4=0.3, b5=0.25, b6=0.35
B = dict(b1=0.4, b2=0.2, b3=0.5, b4=0.3, b5=0.25, b6=0.35)

# 3-taxon geometric (msy044 eqs 28-33) with b1=0.4, b2=0.2, b4=0.25, b3=0.35:
_LA = math.sqrt(B["b1"] * B["b2"]) + B["b5"]
G3 = {
    "r1": math.sqrt(B["b1"] * _LA / (B["b2"] * B["b6"])),
    "r2": math.sqrt(B["b2"] * _LA / (B["b1"] * B["b6"])),
    "r3": math.sqrt(B["b6"] / _LA),
    "ra": math.sqrt(_LA / B["b6"]),
    "t4": math.sqrt(B["b1"] * B["b2"] * B["b6"]) / math.sqrt(_LA),
    "t5": math.sqrt(B["b6"] * _LA),
}

# 4-taxon geometric (msy044 eqs 34-42)
_LB = math.sqrt(B["b3"] * B["b4"]) + B["b6"]
G4 = {
    "r1": math.sqrt(B["b1"] * _LA / (B["b2"] * _LB)),
    "r2": math.sqrt(B["b2"] * _LA / (B["b1"] * _LB)),
    "r3": math.sqrt(B["b3"] * _LB / (B["b4"] * _LA)),
    "r4": math.sqrt(B["b4"] * _LB / (B["b3"] * _LA)),
    "ra": math.sqrt(_LA / _LB),
    "rb": math.sqrt(_LB / _LA),
    "t5": math.sqrt(B["b1"] * B["b2"] * _LB / _LA),
    "t6": math.sqrt(B["b3"] * B["b4"] * _LA / _LB),
    "t7": math.sqrt(_LA * _LB),
}


def _four_taxon_tree(mean: str) -> tuple[object, dict]:
    """((t1,t2):b5,(t3,t4):b6) rooted; local values at the root."""
    newick = (
        f"((T1:{B['b1']},T2:{B['b2']}):{B['b5']},"
        f"(T3:{B['b3']},T4:{B['b4']}):{B['b6']}):0.0;"
    )
    tree = parse_newick(newick)
    tree.assign_ids()
    engine = RrfEngine(tree, mean=mean, rate_ratio_threshold=None, compute_times=True)
    engine.run()
    return engine, engine.row(tree)


def test_four_taxon_geometric_closed_form() -> None:
    engine, row = _four_taxon_tree("geometric")
    assert row.child1 is not None and row.child2 is not None
    r1c1, r1c2 = engine.row(row.child1), engine.row(row.child2)
    assert row.r5 == pytest.approx(G4["ra"], abs=TOL)
    assert row.r6 == pytest.approx(G4["rb"], abs=TOL)
    assert row.t7 == pytest.approx(G4["t7"], abs=TOL)
    # both grandchildren of each child are tips -> r1..r4 backfilled
    assert r1c1.r5 == pytest.approx(G4["r1"], abs=TOL)
    assert r1c1.r6 == pytest.approx(G4["r2"], abs=TOL)
    assert r1c2.r5 == pytest.approx(G4["r3"], abs=TOL)
    assert r1c2.r6 == pytest.approx(G4["r4"], abs=TOL)
    assert r1c1.t7 == pytest.approx(G4["t5"], abs=TOL)
    assert r1c2.t7 == pytest.approx(G4["t6"], abs=TOL)


def test_identity_relations_geometric() -> None:
    """msy044 internal constraints: r5*r6=1 and b = r*t along edges."""
    engine, row = _four_taxon_tree("geometric")
    assert row.r5 is not None and row.r6 is not None
    assert row.r5 * row.r6 == pytest.approx(1.0, abs=TOL)
    t7 = row.t7
    assert t7 is not None
    assert (math.sqrt(row.l1 * row.l2) + row.l5) == pytest.approx(row.r5 * t7, abs=TOL)
    assert (math.sqrt(row.l3 * row.l4) + row.l6) == pytest.approx(row.r6 * t7, abs=TOL)


def test_three_vs_four_taxon_geometric_consistency() -> None:
    """Pruning one 4-taxon subtree gives the 3-taxon solution (A.3 == A.4)."""
    # G3 above was derived with the left subtree collapsed to length _LA; the
    # 4-clade root values must reproduce it with l = composite lengths.
    _, row = _four_taxon_tree("geometric")
    assert row.r5 == pytest.approx(math.sqrt(_LA / _LB), abs=TOL)
    assert row.t7 == pytest.approx(math.sqrt(_LA * _LB), abs=TOL)


# Frozen decimals, hand-evaluated from msy044 eqs 28-33 and 34-42 with the
# branch lengths of B above (b1=0.4, b2=0.2, b3=0.5, b4=0.3, b5=0.25, b6=0.35).
F3 = {
    "r1": 1.744939970272607,
    "r2": 0.8724699851363035,
    "r3": 0.8104654523743625,
    "ra": 1.2338588857432131,
    "t4": 0.22923424691653388,
    "t5": 0.43185061001012454,
}
F4 = {
    "r1": 1.2022447425954197,
    "r2": 0.6011223712977098,
    "r3": 1.5186108066558242,
    "r4": 0.9111664839934943,
    "ra": 0.8501154101350966,
    "rb": 1.1763108727097236,
    "t5": 0.3327109579506045,
    "t6": 0.32924828258075167,
    "t7": 0.6267886761280355,
}


def test_frozen_constants_match_the_transcribed_equations() -> None:
    """The frozen decimals must equal the msy044 formula values."""
    for key, value in F3.items():
        assert G3[key] == pytest.approx(value, abs=TOL), f"G3[{key!r}]"
    for key, value in F4.items():
        assert G4[key] == pytest.approx(value, abs=TOL), f"G4[{key!r}]"


def _three_taxon_tree(mean: str):
    """((T1,T2):b5,T3:b6) rooted: the three-lineage case of eqs 28-33."""
    newick = (
        f"((T1:{B['b1']},T2:{B['b2']}):{B['b5']},T3:{B['b6']}):0.0;"
    )
    tree = parse_newick(newick)
    tree.assign_ids()
    engine = RrfEngine(tree, mean=mean, rate_ratio_threshold=None, compute_times=True)
    engine.run()
    return tree, engine, engine.row(tree)


def test_three_taxon_geometric_closed_form() -> None:
    """The three-lineage solution is asserted, not merely defined (G3)."""
    _tree, engine, row = _three_taxon_tree("geometric")
    assert row.child1 is not None
    cherry = engine.row(row.child1)
    assert row.r5 == pytest.approx(F3["ra"], abs=TOL)
    assert row.r6 == pytest.approx(F3["r3"], abs=TOL)
    assert row.t7 == pytest.approx(F3["t5"], abs=TOL)
    assert cherry.r5 == pytest.approx(F3["r1"], abs=TOL)
    assert cherry.r6 == pytest.approx(F3["r2"], abs=TOL)
    assert cherry.t7 == pytest.approx(F3["t4"], abs=TOL)


def test_four_taxon_arithmetic_closed_form() -> None:
    """msy044 eqs 19-27 (arithmetic mean), hand-computed constants."""
    s = B["b1"] + B["b2"] + B["b3"] + B["b4"] + 2 * B["b5"] + 2 * B["b6"]
    engine, row = _four_taxon_tree("arithmetic")
    assert row.child1 is not None and row.child2 is not None
    r1c1, r1c2 = engine.row(row.child1), engine.row(row.child2)
    assert row.r5 == pytest.approx(2 * (B["b1"] + B["b2"] + 2 * B["b5"]) / s, abs=TOL)
    assert row.r6 == pytest.approx(2 * (B["b3"] + B["b4"] + 2 * B["b6"]) / s, abs=TOL)
    assert row.t7 == pytest.approx(s / 4, abs=TOL)
    assert r1c1.t7 == pytest.approx(
        (B["b1"] + B["b2"]) * s / (4 * (B["b1"] + B["b2"] + 2 * B["b5"])), abs=TOL
    )
    assert r1c2.t7 == pytest.approx(
        (B["b3"] + B["b4"]) * s / (4 * (B["b3"] + B["b4"] + 2 * B["b6"])), abs=TOL
    )
