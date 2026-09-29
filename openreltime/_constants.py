"""Central constants for OpenRelTime.

Every constant documents its provenance. Constants extracted from the GPL-3
reference implementations are cited to the exact source location, as required
by the project licence workflow (see THIRD-PARTY-NOTICES.md).
"""

from __future__ import annotations

from typing import Optional

# ---------------------------------------------------------------------------
# Numerical guards
# ---------------------------------------------------------------------------

#: Zero branch-length substitutes used by R3F.  R3F uses two distinct
#: literals: ``10e-20`` (= 1e-19) in most spots and ``1e-20`` in two spots of
#: the 4-clade case.  Both are kept verbatim so that OpenRelTime reproduces
#: R3F numerically (validation gates G2a/G2b of docs/methods.md, enforced by
#: tests/test_regression_golden.py).
#: Source: R3F/R/rrf_times.R lines 114-131 (GPL-3).
EPS_R3F: float = 10e-20  # 1e-19, R3F literal ``10e-20``
EPS_R3F_SMALL: float = 1e-20  # R3F literal ``1e-20``

#: Fraction of zero-length branches above which a warning is raised.
#: Source: R3F/R/rrf_times.R lines 64-69.
ZERO_BRLEN_WARN_FRACTION: float = 0.10

#: Default rate-ratio protection threshold for relative node times.
#: Nodes whose adjusted rate exceeds the threshold (or its reciprocal) have
#: their node time replaced by the closest non-exceeding ancestor's time.
#: Source: R3F/R/rrf_times.R lines 353-373.
RATE_RATIO_THRESHOLD: float = 20.0

#: Sentinel default that switches the rate-ratio guard off, which is what R3F's
#: rates-only pass does: ``rrf_rates.R`` contains no ``rate.ratio`` block at all
#: and the guard in ``rrf_times.R`` only rewrites node times, never rates.  It
#: is the default of :func:`openreltime.rrf_rates` and of
#: ``RrfEngine(rate_ratio_threshold=...)``.
#: Source: R3F/R/rrf_rates.R lines 250-336 (no guard) vs rrf_times.R 353-373.
RATE_RATIO_GUARD_DISABLED: Optional[float] = None

#: Minimum number of ingroup tips required after outgroup removal.
MIN_INGROUP_TIPS: int = 3

#: Maximum iterations for the calibration violation-adjustment loop.
CALIBRATION_MAX_ITER: int = 100

#: Effective-bounds default number of replicates (msz236 protocol).
EFFECTIVE_BOUNDS_REPLICATES: int = 10_000

# ---------------------------------------------------------------------------
# CorrTest fixed-coefficient logistic regression (Tao et al. 2019)
# ---------------------------------------------------------------------------
# Extracted from R3F/R/corrtest.R lines 555-583 (GPL-3); original reference:
# Q. Tao et al. Mol Biol Evol 36:811-824 (2019), doi:10.1093/molbev/msz014.

#: Normalisation means for [rho_s, rho_ad, rho_ad_1_decay, rho_ad_2_decay].
CORRTEST_NORM_MEAN: tuple[float, float, float, float] = (
    0.436462708,
    0.828259994,
    -0.169515205,
    -0.292940377,
)

#: Normalisation standard deviations, same order as :data:`CORRTEST_NORM_MEAN`.
CORRTEST_NORM_SD: tuple[float, float, float, float] = (
    0.268015804,
    0.087506404,
    0.103192689,
    0.163616884,
)

#: Logistic regression intercept and coefficients (b0..b4).
CORRTEST_COEFFS: tuple[float, float, float, float, float] = (
    -0.07500227,
    6.02875922,
    -0.29265746,
    2.30731322,
    3.2276037,
)


def corrtest_p_band(score: float) -> str:
    """Map a CorrScore to its P-value band (R3F/R/corrtest.R lines 569-584)."""
    if score >= 0.92:
        return "P-value < 0.001"
    if score >= 0.83:
        return "P-value < 0.01"
    if score >= 0.5:
        return "P-value < 0.05"
    return "P-value > 0.05"


# ---------------------------------------------------------------------------
# ddBD grid-search constants
# ---------------------------------------------------------------------------
# Source: R3F/R/ddbd.R lines 399-403 (GPL-3); reference:
# Q. Tao et al. Bioinformatics 37:i102-i110 (2021), doi:10.1093/bioinformatics/btab307

DDBD_BIRTH_GRID: tuple[float, ...] = tuple(float(v) + 0.1 for v in range(1, 11))
DDBD_DEATH_GRID: tuple[float, ...] = tuple(float(v) for v in range(1, 11))
DDBD_SAMPLING_GRID: tuple[float, ...] = (0.001, 0.01, 0.1, 0.5, 0.9)
DDBD_MAX_START_ATTEMPTS: int = 50

# ---------------------------------------------------------------------------
# BLB (bag of little bootstraps) defaults
# ---------------------------------------------------------------------------

BLB_DEFAULT_GAMMA: float = 0.7
BLB_GAMMA_RANGE: tuple[float, float] = (0.6, 0.9)

# ---------------------------------------------------------------------------
# Software metadata
# ---------------------------------------------------------------------------

VERSION: str = "0.1.0"
SOFTWARE_NAME: str = "OpenRelTime"
SOFTWARE_URL: str = "https://github.com/ZengZichao/OpenRelTime"
