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

Namespace contract
------------------
Three analysis modules are named after the function they expose
(``calibrate``, ``corrtest``, ``ddbd``).  A package attribute can hold only one
object, so the *module* wins the binding -- ``openreltime.calibrate`` is the
submodule and ``openreltime.calibrate.Calibration`` therefore resolves -- and the
module object is made callable as its own eponymous entry point
(:class:`_EntryPointModule`).  Every one of these forms works::

    import openreltime
    openreltime.calibrate(times, cals)              # call the analysis
    openreltime.calibrate.Calibration(node_id=9)    # read the module's symbols
    from openreltime.calibrate import calibrate     # the plain function
    from openreltime import calibrate as calibrate_module   # the module

The callable-module binding is what keeps the documented ``openreltime.calibrate(
...)`` call form and the OpenRelTime-Studio adapter working unchanged; nothing
here silently rebinds a module to a function.  ``__all__`` lists the
function-level API; the submodules are reachable as attributes and through the
dotted import path.

References
----------
* Tamura K, Tao Q, Kumar S. Mol Biol Evol 35:1770-1782 (2018) -- RRF.
* Tamura K et al. PNAS 109:19334-19339 (2012) -- RelTime.
* Tao Q et al. Mol Biol Evol 37:280-290 (2020) -- CIs and effective bounds.
* Tao Q et al. Mol Biol Evol 36:811-824 (2019) -- CorrTest.
* Tao Q et al. Bioinformatics 37:i102-i110 (2021) -- ddBD.
"""

from __future__ import annotations

import sys
from types import ModuleType
from typing import Any

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


class _EntryPointModule(ModuleType):
    """Submodule that also *calls* the entry-point function it defines.

    ``openreltime.calibrate`` has to name the module for
    ``openreltime.calibrate.Calibration`` to resolve, while the documented
    ``openreltime.calibrate(times, calibrations)`` call form has to keep working
    for existing callers.  A module whose ``__call__`` forwards to the attribute
    named by ``__entry_point__`` satisfies both at once.  The forwarded call is
    resolved lazily, so nothing is duplicated and the function object stays the
    single implementation.
    """

    #: name of the attribute ``__call__`` forwards to (set per submodule)
    __entry_point__: str = ""

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        target = getattr(self, self.__entry_point__)
        return target(*args, **kwargs)


def _expose(module_path: str, entry_point: str) -> ModuleType:
    """Bind submodule *module_path* as a callable alias of *entry_point*.

    The package attribute keeps pointing at the module object (so all of its
    symbols stay reachable through it) and gains callability that forwards to
    the module's own eponymous function.
    """
    module = sys.modules[module_path]
    module.__class__ = _EntryPointModule  # type: ignore[assignment]
    module.__entry_point__ = entry_point  # type: ignore[attr-defined]
    globals()[module.__name__.rpartition(".")[2]] = module
    return module


for _path, _entry in (
    ("openreltime.calibrate", "calibrate"),
    ("openreltime.corrtest", "corrtest"),
    ("openreltime.ddbd", "ddbd"),
):
    _expose(_path, _entry)
del _path, _entry

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
