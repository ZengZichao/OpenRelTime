"""OpenRelTime: relative rate framework (RRF) molecular dating in Python.

Public API::

    import openreltime

    tree = openreltime.read_tree("tree.nwk", outgroup=["Out1", "Out2"])
    rates = openreltime.rrf_rates(tree)
    times = openreltime.rrf_times(tree)
    cal = openreltime.calibrate(times, calibrations)
    ci = openreltime.confidence_interval(cal)
    ct = openreltime.corrtest(tree)
    bd = openreltime.ddbd(tree, anchor_time=1.0)

References
----------
* Tamura K, Tao Q, Kumar S. Mol Biol Evol 35:1770-1782 (2018) -- RRF.
* Tamura K et al. PNAS 109:19334-19339 (2012) -- RelTime.
* Tao Q et al. Mol Biol Evol 37:280-290 (2020) -- CIs and effective bounds.
* Tao Q et al. Mol Biol Evol 36:811-824 (2019) -- CorrTest.
* Tao Q et al. Bioinformatics 37:i102-i110 (2021) -- ddBD.
"""

from __future__ import annotations

from . import _constants
from ._constants import VERSION as __version__
from ._constants import VERSION
from .calibrate import (
    Calibration,
    CalibrationCancelled,
    calibrate,
    parse_calibrations,
)
from .ci import confidence_interval
from .corrtest import corrtest
from .ddbd import ddbd
from .rates import rrf_rates
from .report import (
    CalibratedResult,
    CIResult,
    CorrTestResult,
    DDBDResult,
    RateResult,
    TimeResult,
)
from .table import tree2table
from .taxonomy import (
    MonophylyResult,
    check_monophyly,
    detect_groups_from_labels,
    parse_taxon_table,
)
from .times import rrf_rates_times, rrf_times
from .tree import PhyloNode
from .treeio import parse_newick, read_tree, to_newick, write_nexus

__all__ = [
    "__version__",
    "VERSION",
    "read_tree",
    "parse_newick",
    "to_newick",
    "write_nexus",
    "rrf_rates",
    "rrf_times",
    "rrf_rates_times",
    "calibrate",
    "Calibration",
    "CalibrationCancelled",
    "parse_calibrations",
    "confidence_interval",
    "corrtest",
    "ddbd",
    "tree2table",
    "check_monophyly",
    "parse_taxon_table",
    "detect_groups_from_labels",
    "MonophylyResult",
    "PhyloNode",
    "RateResult",
    "TimeResult",
    "CalibratedResult",
    "CIResult",
    "CorrTestResult",
    "DDBDResult",
    "_constants",
]
