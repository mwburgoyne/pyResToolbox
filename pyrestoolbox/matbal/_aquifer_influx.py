"""
Cumulative aquifer influx by the Carter-Tracy (1960) method, using the Van
Everdingen-Hurst constant-terminal-rate dimensionless pressure from
simtools.influence_tables.

Carter & Tracy (1960), An Improved Method for Calculating Water Influx,
Trans. AIME 219, 415 (SPE-1626-G):
    We_n = We_{n-1} + (tD_n - tD_{n-1}) (B dp_n - We_{n-1} pD'_n) / (pD_n - tD_{n-1} pD'_n)
with dp_n = p_i - p_n the total pressure drop at the aquifer boundary. This is the
influx model behind ECLIPSE AQUCT; it is an approximation to Van Everdingen-Hurst
superposition that improves as the time steps shorten.
"""

import math
import warnings

import numpy as np

# Dimensionless-time grid for the pD table: points per decade, and margins below the
# first and above the last tD so the spline is not evaluated at its ends.
_PD_POINTS_PER_DECADE = 12
_PD_GRID_MARGIN = 2.0
# An infinite aquifer is represented by a radius whose boundary is not felt over the
# requested times: pD departs from the infinite-acting line near tD ~ 0.25 reD^2.
_INFINITE_RED_MIN = 1.0e4
_INFINITE_RED_FACTOR = 10.0


def aquifer_influx(t, p, B, td_per_t, red=None, M=7):
    """Cumulative aquifer influx history by the Carter-Tracy method.

    Units are consistent rather than fixed: We has the units of B times pressure
    (rb/psi x psi = rb, or rm3/bar x bar = rm3), and t is in whatever unit td_per_t
    was built for. For a radial aquifer in field units, B = 1.119 phi ct h ro^2
    (theta / 360) rb/psi and tD per day = 0.0063283 k / (phi mu ct ro^2), with k in
    mD, mu in cP, ct in 1/psi and ro, h in ft. The first entry of t and p is the
    initial state (influx zero). Shorter time steps track Van Everdingen-Hurst more
    closely; steps of 30 days or less are usually adequate.

    t: Times, increasing, starting at the initial state (any unit, matching td_per_t)
    p: Pressure at the reservoir-aquifer boundary at each time (psia | barsa, matching B)
    B: Aquifer influx constant (rb/psi | rm3/bar)
    td_per_t: Dimensionless time per unit of t
    red: Aquifer to reservoir radius ratio re/ro (> 1). None (default) for an
        infinite-acting aquifer
    M: Laplace inversion accuracy passed to simtools.influence_tables (default 7)

    Returns dict with 'We' (cumulative influx, B x pressure units), 'tD', 'pD' and
    'dpD' (pD derivative with respect to tD), each an array aligned with t.
    """
    t = np.asarray(t, dtype=float)
    p = np.asarray(p, dtype=float)
    if t.ndim != 1 or len(t) != len(p) or len(t) < 2:
        raise ValueError("t and p must be 1-D arrays of the same length (at least 2)")
    if not (np.all(np.isfinite(t)) and np.all(np.isfinite(p))):
        raise ValueError("t and p must be finite")
    if np.any(np.diff(t) <= 0):
        raise ValueError("t must be strictly increasing")
    if B <= 0 or td_per_t <= 0:
        raise ValueError(f"B and td_per_t must be positive, got {B} and {td_per_t}")
    if red is not None and red <= 1:
        raise ValueError(f"red must be greater than 1, got {red}")

    td = (t - t[0]) * td_per_t
    pd, dpd = _pd_and_derivative(td[1:], red, M)
    dp = p[0] - p

    we = np.zeros(len(t))
    for n in range(1, len(t)):
        pd_n, dpd_n = pd[n - 1], dpd[n - 1]
        we[n] = we[n - 1] + (td[n] - td[n - 1]) * (
            (B * dp[n] - we[n - 1] * dpd_n) / (pd_n - td[n - 1] * dpd_n))
    return {'We': we, 'tD': td, 'pD': np.concatenate([[0.0], pd]),
            'dpD': np.concatenate([[np.nan], dpd])}


def _pd_and_derivative(td, red, M):
    """pD and dpD/dtD at td (all > 0) from a log-spaced influence table spline."""
    from scipy.interpolate import CubicSpline
    from pyrestoolbox.simtools import influence_tables

    td_min, td_max = td.min() / _PD_GRID_MARGIN, td.max() * _PD_GRID_MARGIN
    if red is None:
        red = max(_INFINITE_RED_MIN, _INFINITE_RED_FACTOR * math.sqrt(td_max))
    elif td.max() > red ** 2:
        warnings.warn(f"tD reaches {td.max():.3g}, beyond reD^2 = {red ** 2:.3g}: the "
                      f"aquifer is depleting and the pD table is far into pseudo-steady "
                      f"state")
    n_incr = max(8, int(math.ceil(_PD_POINTS_PER_DECADE * math.log10(td_max / td_min))))
    td_grid, pds = influence_tables([red], min_td=td_min, max_td=td_max,
                                    n_incr=n_incr, M=M)
    spline = CubicSpline(np.log(np.asarray(td_grid, dtype=float)),
                         np.asarray(pds[0], dtype=float))
    x = np.log(td)
    return spline(x), spline(x, 1) / td
