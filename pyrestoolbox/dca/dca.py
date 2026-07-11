#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
    pyResToolbox - A collection of Reservoir Engineering Utilities
              Copyright (C) 2022, Mark Burgoyne

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    The GNU General Public License can be found in the LICENSE directory,
    and at  <https://www.gnu.org/licenses/>.

          Contact author at mark.w.burgoyne@gmail.com

Decline Curve Analysis (DCA) functions.

Functions
---------
arps_rate          Arps decline rate (exponential, hyperbolic, harmonic)
arps_cum           Arps cumulative production
mh_rate            Modified hyperbolic (hyperbolic-to-exponential) decline rate
mh_cum             Modified hyperbolic cumulative production
mh_eur             Modified hyperbolic estimated ultimate recovery
hyp2_rate          Two-segment hyperbolic (transition at telf) decline rate
hyp2_cum           Two-segment hyperbolic cumulative production
hyp2_eur           Two-segment hyperbolic estimated ultimate recovery
hyp2_from_eur      Solve two-segment hyperbolic parameters for a target EUR
duong_rate         Duong decline rate for unconventionals
duong_cum          Duong cumulative production
eur                Estimated ultimate recovery
fit_decline        Fit decline model to production data (time domain)
fit_decline_cum    Fit decline model to rate-vs-cumulative data
fit_ratio          Fit ratio model (e.g. GOR, WOR) to data
ratio_forecast     Evaluate fitted ratio model
forecast           Generate rate and cumulative forecast from fitted model

Classes
-------
DeclineResult   Result from fit_decline() or fit_decline_cum()
ForecastResult  Result from forecast()
RatioResult     Result from fit_ratio()
"""

__all__ = [
    'arps_rate', 'arps_cum', 'mh_rate', 'mh_cum', 'mh_eur',
    'hyp2_rate', 'hyp2_cum', 'hyp2_eur', 'hyp2_from_eur',
    'duong_rate', 'duong_cum',
    'eur', 'fit_decline', 'fit_decline_cum', 'forecast',
    'fit_ratio', 'ratio_forecast',
    'DeclineResult', 'ForecastResult', 'RatioResult',
]

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Union
from numpy.typing import ArrayLike

from pyrestoolbox.shared_fns import convert_to_numpy, process_output, ransac_linreg
from pyrestoolbox._accelerator import RUST_AVAILABLE as _RUST_AVAILABLE
if _RUST_AVAILABLE:
    from pyrestoolbox import _native as _rust


# ---------------------------------------------------------------------------
# Named constants — extracted from inline magic numbers per CLAUDE.md rules.
# ---------------------------------------------------------------------------

# Arps hyperbolic b-grid search (Arps 1945; boundary-dominated b in (0, 1))
# 0.05-0.95 avoids exponential (b=0) / harmonic (b=1) degeneracies.
_B_GRID_MIN = 0.05
_B_GRID_MAX = 0.96
_B_GRID_STEP = 0.01

# Modified hyperbolic (Robertson 1988, SPE-18731-MS) curve_fit settings.
# dterm is fitted as f * di with f in (0, 1) so the terminal decline can
# never reach the initial decline during optimisation.
_MH_FIT_B_MIN = 0.05
_MH_FIT_B_MAX = 2.0            # petbox-dca MH prior upper bound
_MH_FIT_F_MIN = 1e-3
_MH_FIT_F_MAX = 0.999
_MH_FIT_QI_FACTOR = 5.0        # qi bounds = first rate / and * factor
_MH_FIT_DI_SPAN = 1e3          # di bounds = log-slope guess / and * span
_MH_FIT_DI_GUESS_MIN = 1e-6    # Floor for slope-derived di initial guess
_MH_P0_B = 1.0                 # Initial b guess (petbox-dca deterministic fit)
_MH_P0_F = 0.1                 # Initial dterm/di fraction guess

# Two-segment hyperbolic (sharp-transition form of the transient hyperbolic,
# Fulford & Blasingame 2013, SPE-167242-MS) curve_fit settings. b2 is fitted
# as r * b1 (r in (0, 1]) so the second-segment b-factor can never exceed the
# first. qi/di/b1 bounds are shared with the MH fitter constants above.
_HYP2_P0_B1 = 2.0              # Initial b1 guess (petbox-dca THM: bi = 2.0)
_HYP2_P0_R = 0.5               # Initial b2/b1 ratio (THM p0: bf=1.0 / bi=2.0)
_HYP2_FIT_R_MIN = 0.01
_HYP2_P0_TELF_IDX = 3          # Initial telf guess = t[len(t) // idx] (THM p0)

# Duong (2011) cumulative integration grid
_DUONG_TRAP_LB = 0.001      # Lower trap bound to avoid t=0 singularity
_DUONG_GRID_MIN = 500       # Minimum trap grid points per integration
_DUONG_GRID_DENSITY = 10    # Additional points per unit t

# Duong curve_fit bounds and initial guess (Duong 2011, Eq. 5)
_DUONG_BOUNDS_LO = (0.0, 0.01, 1.001)
_DUONG_BOUNDS_HI_QI_FACTOR = 5.0   # Upper qi bound = first q * factor
_DUONG_BOUNDS_HI_A = 10.0
_DUONG_BOUNDS_HI_M = 3.0
_DUONG_P0_A = 1.0              # Initial guess a
_DUONG_P0_M = 1.2              # Initial guess m

# scipy.curve_fit iteration cap (covers Duong + logistic fits)
_CURVE_FIT_MAXFEV = 5000

# Numerical floors / near-zero guards
_HYPER_INNER_FLOOR = 1e-10     # Floor for (1 - di*Np*exp/qi)^(1/exp) argument
_ZERO_DIV_EPS = 1e-30          # General division-by-zero guard

# Logistic ratio curve_fit (unconventional GOR/WOR trending)
_LOGISTIC_P0_RMAX_INIT = 1.5   # Initial Rmax guess = max(ratio) * factor
_LOGISTIC_P0_TC = 0.01         # Initial carryover constant
_LOGISTIC_P0_ALPHA = 10.0      # Initial curvature
_LOGISTIC_BOUNDS_LO = (1e-8, 1e-3)  # Lower (tc, alpha) bounds
_LOGISTIC_BOUNDS_RMAX_FACTOR = 5.0  # Upper Rmax = max(ratio) * factor
_LOGISTIC_BOUNDS_HI = 1e6      # Upper tc / alpha bound


def _build_decline_result(method, q_obs, q_pred, **params):
    """DRY helper: compute R-squared + residuals, return a DeclineResult.

    Parameters
    ----------
    method : str
    q_obs, q_pred : np.ndarray  (observed and predicted rates)
    **params : scalar DeclineResult fields (qi, di, b, a, m, dterm)
    """
    residuals = q_obs - q_pred
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((q_obs - np.mean(q_obs)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return DeclineResult(
        method=method,
        r_squared=r2,
        residuals=residuals,
        **params,
    )


def _fmt_array(a, name):
    if a is None:
        return f"{name}=None"
    arr = np.asarray(a)
    return f"{name}=ndarray(shape={arr.shape}, dtype={arr.dtype})"


@dataclass
class DeclineResult:
    """Result from decline curve fitting.

    Attributes
    ----------
    method : str
        Decline model ('exponential', 'harmonic', 'hyperbolic', 'mh', 'hyp2',
        'duong').
    qi : float
        Initial rate.
    di : float
        Initial decline rate (1/time). Not used for Duong.
    b : float
        Arps b-factor (first-segment b1 for 'hyp2'). 0 for exponential,
        1 for harmonic. Not used for Duong.
    a : float
        Duong 'a' parameter. 0 for Arps models.
    m : float
        Duong 'm' parameter. 0 for Arps models.
    dterm : float
        Terminal decline rate (1/time) for 'mh' and 'hyp2'. 0 for other models.
    b2 : float
        Second-segment b-factor for 'hyp2'. 0 for other models.
    telf : float
        Transition time for 'hyp2'. 0 for other models.
    r_squared : float
        Coefficient of determination of the fit.
    residuals : np.ndarray
        Residual array (observed - predicted).
    uptime_mean : float, optional
        Mean uptime fraction from fit_decline_cum() when t_calendar provided.
    uptime_history : np.ndarray, optional
        Per-interval uptime fractions from fit_decline_cum() when t_calendar provided.
    """
    method: str
    qi: float
    di: float = 0.0
    b: float = 0.0
    a: float = 0.0
    m: float = 0.0
    dterm: float = 0.0
    b2: float = 0.0
    telf: float = 0.0
    r_squared: float = 0.0
    residuals: np.ndarray = field(default_factory=lambda: np.array([]))
    uptime_mean: Optional[float] = None
    uptime_history: Optional[np.ndarray] = None

    def __repr__(self):
        scalars = (
            f"method={self.method!r}, qi={self.qi}, di={self.di}, b={self.b}, "
            f"a={self.a}, m={self.m}, dterm={self.dterm}, "
            f"b2={self.b2}, telf={self.telf}, "
            f"r_squared={self.r_squared}, uptime_mean={self.uptime_mean}"
        )
        return (
            f"DeclineResult({scalars}, "
            f"{_fmt_array(self.residuals, 'residuals')}, "
            f"{_fmt_array(self.uptime_history, 'uptime_history')})"
        )


@dataclass
class ForecastResult:
    """Result from production forecast.

    Attributes
    ----------
    t : np.ndarray
        Time array.
    q : np.ndarray
        Rate array.
    Qcum : np.ndarray
        Cumulative production array.
    eur : float
        Estimated ultimate recovery (cumulative at end of forecast).
    secondary : dict, optional
        Secondary phase forecasts keyed by name. Each value is a dict with
        'ratio', 'rate', and 'cum' arrays.
    """
    t: np.ndarray
    q: np.ndarray
    Qcum: np.ndarray
    eur: float
    secondary: Optional[dict] = None

    def __repr__(self):
        sec = 'None' if self.secondary is None else f"{{{', '.join(self.secondary.keys())}}}"
        return (
            f"ForecastResult("
            f"{_fmt_array(self.t, 't')}, "
            f"{_fmt_array(self.q, 'q')}, "
            f"{_fmt_array(self.Qcum, 'Qcum')}, "
            f"eur={self.eur}, secondary={sec})"
        )


@dataclass
class RatioResult:
    """Result from ratio fitting.

    Attributes
    ----------
    method : str
        Ratio model ('linear', 'exponential', 'power', 'logistic').
    a : float
        Primary parameter (intercept / coefficient / Rmax for logistic).
    b : float
        Slope / exponent / growth rate.
    c : float
        Logistic offset parameter (only used for logistic model).
    domain : str
        'time' or 'cum' — tells forecast() which x-axis to use.
    r_squared : float
        Coefficient of determination of the fit.
    residuals : np.ndarray
        Residual array (observed - predicted).
    """
    method: str
    a: float
    b: float
    c: float = 0.0
    domain: str = 'cum'
    r_squared: float = 0.0
    residuals: np.ndarray = field(default_factory=lambda: np.array([]))

    def __repr__(self):
        return (
            f"RatioResult(method={self.method!r}, a={self.a}, b={self.b}, "
            f"c={self.c}, domain={self.domain!r}, r_squared={self.r_squared}, "
            f"{_fmt_array(self.residuals, 'residuals')})"
        )


def arps_rate(qi: float, di: float, b: float, t: ArrayLike) -> Union[float, np.ndarray]:
    """Arps decline rate.

    Parameters
    ----------
    qi : float or array
        Initial rate (volume/time).
    di : float
        Initial decline rate (1/time). Must be > 0.
    b : float
        Arps b-factor. b=0: exponential, 0<b<1: hyperbolic, b=1: harmonic,
        b>1: super-hyperbolic (transient flow, common in unconventionals).
    t : float or array
        Time.

    Returns
    -------
    float or np.ndarray
        Rate at time t.
    """
    t, is_list = convert_to_numpy(t)
    if qi <= 0:
        raise ValueError(f"qi must be positive, got {qi}")
    if di <= 0:
        raise ValueError(f"di must be positive, got {di}")
    if b < 0:
        raise ValueError(f"b must be non-negative, got {b}")

    if b == 0:
        q = qi * np.exp(-di * t)
    elif b == 1:
        q = qi / (1.0 + di * t)
    else:
        q = qi / (1.0 + b * di * t) ** (1.0 / b)

    return process_output(q, is_list)


def arps_cum(qi: float, di: float, b: float, t: ArrayLike) -> Union[float, np.ndarray]:
    """Arps cumulative production.

    Parameters
    ----------
    qi : float or array
        Initial rate (volume/time).
    di : float
        Initial decline rate (1/time). Must be > 0.
    b : float
        Arps b-factor. b=0: exponential, 0<b<1: hyperbolic, b=1: harmonic,
        b>1: super-hyperbolic (transient flow; cumulative is unbounded as
        t increases unless paired with a terminal decline - see mh_cum).
    t : float or array
        Time.

    Returns
    -------
    float or np.ndarray
        Cumulative production at time t.
    """
    t, is_list = convert_to_numpy(t)
    if qi <= 0:
        raise ValueError(f"qi must be positive, got {qi}")
    if di <= 0:
        raise ValueError(f"di must be positive, got {di}")
    if b < 0:
        raise ValueError(f"b must be non-negative, got {b}")

    if b == 0:
        Qcum = (qi / di) * (1.0 - np.exp(-di * t))
    elif b == 1:
        Qcum = (qi / di) * np.log(1.0 + di * t)
    else:
        Qcum = (qi / ((1.0 - b) * di)) * (1.0 - (1.0 + b * di * t) ** (-(1.0 - b) / b))

    return process_output(Qcum, is_list)


def duong_rate(qi: float, a: float, m: float, t: ArrayLike) -> Union[float, np.ndarray]:
    """Duong decline rate for unconventional reservoirs.

    q(t) = qi * t^(-m) * exp(a/(1-m) * (t^(1-m) - 1))

    Parameters
    ----------
    qi : float
        Rate at t=1 (volume/time).
    a : float
        Duong 'a' parameter. Must be > 0.
    m : float
        Duong 'm' parameter. Must be > 1.
    t : float or array
        Time (must be > 0).

    Returns
    -------
    float or np.ndarray
        Rate at time t.
    """
    t, is_list = convert_to_numpy(t)
    if a <= 0:
        raise ValueError(f"Duong parameter 'a' must be positive, got {a}")
    if m <= 1:
        raise ValueError(f"Duong parameter 'm' must be > 1, got {m}")
    if np.any(t <= 0):
        raise ValueError("Time must be positive for Duong model")

    q = qi * t ** (-m) * np.exp(a / (1.0 - m) * (t ** (1.0 - m) - 1.0))
    return process_output(q, is_list)


def duong_cum(qi: float, a: float, m: float, t: ArrayLike) -> Union[float, np.ndarray]:
    """Duong cumulative production via trapezoidal integration.

    Parameters
    ----------
    qi : float
        Rate at t=1 (volume/time).
    a : float
        Duong 'a' parameter.
    m : float
        Duong 'm' parameter.
    t : float or array
        Time (must be > 0).

    Returns
    -------
    float or np.ndarray
        Cumulative production at time t.
    """
    t, is_list = convert_to_numpy(t)
    if a <= 0:
        raise ValueError(f"Duong parameter 'a' must be positive, got {a}")
    if m <= 1:
        raise ValueError(f"Duong parameter 'm' must be > 1, got {m}")
    if np.any(t <= 0):
        raise ValueError("Time must be positive for Duong model")

    # Generate fine time grid for integration. Lower bound must be strictly
    # less than ti to give np.linspace an ascending range, otherwise trapezoid
    # integrates over a descending axis and returns negative cumulative.
    results = np.zeros_like(t, dtype=float)
    for i, ti in enumerate(t):
        lower = min(_DUONG_TRAP_LB, ti * _DUONG_TRAP_LB)
        t_fine = np.linspace(lower, ti, max(_DUONG_GRID_MIN, int(ti * _DUONG_GRID_DENSITY)))
        q_fine = qi * t_fine ** (-m) * np.exp(a / (1.0 - m) * (t_fine ** (1.0 - m) - 1.0))
        results[i] = np.trapezoid(q_fine, t_fine)

    return process_output(results, is_list)


def eur(qi: float, di: float, b: float, q_min: float) -> float:
    """Estimated ultimate recovery for Arps decline.

    Parameters
    ----------
    qi : float
        Initial rate.
    di : float
        Initial decline rate (1/time).
    b : float
        Arps b-factor. b > 1 is permitted (EUR remains finite for any
        q_min > 0), but see mh_eur for the two-segment form usually used
        with transient-flow b-factors.
    q_min : float
        Economic limit rate.

    Returns
    -------
    float
        EUR (cumulative production when rate reaches q_min).
    """
    if q_min >= qi:
        raise ValueError(f"q_min ({q_min}) must be less than qi ({qi})")
    if di <= 0:
        raise ValueError(f"di must be positive, got {di}")
    if b < 0:
        raise ValueError(f"b must be non-negative, got {b}")

    # Solve for time when q(t) = q_min
    if b == 0:
        t_end = -np.log(q_min / qi) / di
    elif b == 1:
        t_end = (qi / q_min - 1.0) / di
    else:
        t_end = ((qi / q_min) ** b - 1.0) / (b * di)

    return float(arps_cum(qi, di, b, t_end))


def _mh_switch(qi, di, b, dterm):
    """Hyperbolic-to-exponential switch point for the modified hyperbolic model.

    Returns (t_sw, q_sw, N_sw) where the nominal decline D(t) = di/(1+b*di*t)
    first equals dterm, or None when no terminal segment applies (dterm=0
    disabled, or b=0 where the nominal decline is constant at di and never
    declines to dterm).
    """
    if dterm < 0:
        raise ValueError(f"dterm must be non-negative, got {dterm}")
    if dterm == 0:
        return None
    if dterm >= di:
        raise ValueError(f"dterm ({dterm}) must be less than di ({di})")
    if b == 0:
        return None
    t_sw = (1.0 / dterm - 1.0 / di) / b
    q_sw = float(arps_rate(qi, di, b, t_sw))
    N_sw = float(arps_cum(qi, di, b, t_sw))
    return t_sw, q_sw, N_sw


def mh_rate(qi: float, di: float, t: ArrayLike, b: float = 2.0,
            dterm: float = 0.0) -> Union[float, np.ndarray]:
    """Modified hyperbolic (two-segment Arps) decline rate.

    Hyperbolic decline with b-factor b until the nominal decline
    D(t) = di / (1 + b*di*t) falls to dterm, exponential decline at dterm
    thereafter (Robertson 1988, SPE-18731-MS). Rate and nominal decline are
    both continuous at the switch. With dterm=0 (or b=0) the curve is a
    single Arps segment identical to arps_rate.

    Parameters
    ----------
    qi : float
        Initial rate (volume/time).
    di : float
        Initial decline rate (1/time). Must be > 0.
    t : float or array
        Time.
    b : float
        First-segment Arps b-factor (default 2.0, transient flow). Values
        above 1 are the usual application; any b >= 0 is accepted.
    dterm : float
        Terminal nominal decline rate (1/time). Must be < di. Default 0
        disables the terminal segment.

    Returns
    -------
    float or np.ndarray
        Rate at time t.
    """
    t, is_list = convert_to_numpy(t)
    sw = _mh_switch(qi, di, b, dterm)
    q = np.asarray(arps_rate(qi, di, b, t), dtype=float)
    if sw is not None:
        t_sw, q_sw, _ = sw
        late = t > t_sw
        q[late] = q_sw * np.exp(-dterm * (t[late] - t_sw))
    return process_output(q, is_list)


def mh_cum(qi: float, di: float, t: ArrayLike, b: float = 2.0,
           dterm: float = 0.0) -> Union[float, np.ndarray]:
    """Modified hyperbolic (two-segment Arps) cumulative production.

    Analytic piecewise integral of mh_rate: Arps cumulative to the switch
    time, exponential-segment cumulative beyond it.

    Parameters
    ----------
    qi : float
        Initial rate (volume/time).
    di : float
        Initial decline rate (1/time). Must be > 0.
    t : float or array
        Time.
    b : float
        First-segment Arps b-factor (default 2.0, transient flow).
    dterm : float
        Terminal nominal decline rate (1/time). Must be < di. Default 0
        disables the terminal segment.

    Returns
    -------
    float or np.ndarray
        Cumulative production at time t.
    """
    t, is_list = convert_to_numpy(t)
    sw = _mh_switch(qi, di, b, dterm)
    if sw is None:
        Qcum = np.asarray(arps_cum(qi, di, b, t), dtype=float)
    else:
        t_sw, q_sw, N_sw = sw
        Qcum = np.asarray(arps_cum(qi, di, b, np.minimum(t, t_sw)), dtype=float)
        late = t > t_sw
        Qcum[late] = N_sw + (q_sw / dterm) * (1.0 - np.exp(-dterm * (t[late] - t_sw)))
    return process_output(Qcum, is_list)


def mh_eur(qi: float, di: float, q_min: float, b: float = 2.0,
           dterm: float = 0.0) -> float:
    """Estimated ultimate recovery for modified hyperbolic decline.

    Parameters
    ----------
    qi : float
        Initial rate.
    di : float
        Initial decline rate (1/time).
    q_min : float
        Economic limit rate. Must be > 0 when b >= 1 and dterm = 0
        (the single-segment cumulative is unbounded).
    b : float
        First-segment Arps b-factor (default 2.0, transient flow).
    dterm : float
        Terminal nominal decline rate (1/time). Must be < di. Default 0
        disables the terminal segment.

    Returns
    -------
    float
        EUR (cumulative production when rate reaches q_min).
    """
    if q_min >= qi:
        raise ValueError(f"q_min ({q_min}) must be less than qi ({qi})")
    sw = _mh_switch(qi, di, b, dterm)
    if sw is None:
        if q_min <= 0 and b >= 1:
            raise ValueError(
                "EUR is unbounded for b >= 1 with no terminal decline; "
                "specify dterm or a positive q_min")
        return eur(qi, di, b, q_min)
    t_sw, q_sw, N_sw = sw
    if q_min >= q_sw:
        # Abandonment reached within the hyperbolic segment
        return eur(qi, di, b, q_min)
    # Exponential segment: integral from q_sw down to q_min is (q_sw - q_min)/dterm
    return N_sw + (q_sw - q_min) / dterm


def _hyp2_switch(qi, di, b1, b2, telf):
    """Segment-2 anchor for the two-segment hyperbolic model.

    Returns (q_sw, D_sw, N_sw) at t = telf: the first-segment rate, nominal
    decline and cumulative. Starting segment 2 from (q_sw, D_sw) makes rate
    and log-slope (nominal decline) both continuous at telf by construction.
    """
    if telf < 0:
        raise ValueError(f"telf must be non-negative, got {telf}")
    if b2 > b1:
        raise ValueError(f"b2 ({b2}) must not exceed b1 ({b1})")
    q_sw = float(arps_rate(qi, di, b1, telf))
    D_sw = di / (1.0 + b1 * di * telf)
    N_sw = float(arps_cum(qi, di, b1, telf))
    return q_sw, D_sw, N_sw


def hyp2_rate(qi: float, di: float, t: ArrayLike, telf: float, b1: float = 2.0,
              b2: float = 0.5, dterm: float = 0.0) -> Union[float, np.ndarray]:
    """Two-segment hyperbolic decline rate with transition at telf.

    Hyperbolic decline with b-factor b1 (default 2.0, transient linear flow)
    until the specified transition time telf, then a second hyperbolic segment
    with b-factor b2 whose parameters are set so that rate and nominal decline
    (log-slope) are identical to the first segment at telf. This is the
    sharp-transition form of the transient hyperbolic model (Fulford and
    Blasingame 2013, SPE-167242-MS). An optional terminal decline dterm
    converts the second segment to exponential once its nominal decline falls
    to dterm (as in mh_rate).

    Parameters
    ----------
    qi : float
        Initial rate (volume/time).
    di : float
        Initial nominal decline rate (1/time). Must be > 0.
    t : float or array
        Time.
    telf : float
        Transition time (end of linear flow). Same time units as t.
    b1 : float
        First-segment Arps b-factor (default 2.0, transient flow).
    b2 : float
        Second-segment Arps b-factor (default 0.5, boundary-dominated).
        Must not exceed b1.
    dterm : float
        Terminal nominal decline rate (1/time) applied to the second segment.
        Must be < the nominal decline at telf. Default 0 disables it.

    Returns
    -------
    float or np.ndarray
        Rate at time t.
    """
    t, is_list = convert_to_numpy(t)
    q_sw, D_sw, _ = _hyp2_switch(qi, di, b1, b2, telf)
    _mh_switch(q_sw, D_sw, b2, dterm)  # validate dterm against segment-2 decline
    q = np.asarray(arps_rate(qi, di, b1, t), dtype=float)
    late = t > telf
    if np.any(late):
        q[late] = mh_rate(q_sw, D_sw, t[late] - telf, b=b2, dterm=dterm)
    return process_output(q, is_list)


def hyp2_cum(qi: float, di: float, t: ArrayLike, telf: float, b1: float = 2.0,
             b2: float = 0.5, dterm: float = 0.0) -> Union[float, np.ndarray]:
    """Two-segment hyperbolic cumulative production.

    Analytic piecewise integral of hyp2_rate: Arps cumulative to telf, then
    the modified hyperbolic cumulative of the second segment re-anchored at
    (q_sw, D_sw).

    Parameters
    ----------
    Same as hyp2_rate.

    Returns
    -------
    float or np.ndarray
        Cumulative production at time t.
    """
    t, is_list = convert_to_numpy(t)
    q_sw, D_sw, N_sw = _hyp2_switch(qi, di, b1, b2, telf)
    _mh_switch(q_sw, D_sw, b2, dterm)  # validate dterm against segment-2 decline
    Qcum = np.asarray(arps_cum(qi, di, b1, np.minimum(t, telf)), dtype=float)
    late = t > telf
    if np.any(late):
        Qcum[late] = N_sw + mh_cum(q_sw, D_sw, t[late] - telf, b=b2, dterm=dterm)
    return process_output(Qcum, is_list)


def hyp2_eur(qi: float, di: float, q_min: float, telf: float, b1: float = 2.0,
             b2: float = 0.5, dterm: float = 0.0) -> float:
    """Estimated ultimate recovery for two-segment hyperbolic decline.

    Parameters
    ----------
    qi : float
        Initial rate.
    di : float
        Initial nominal decline rate (1/time).
    q_min : float
        Economic limit rate. Must be > 0 when b2 >= 1 and dterm = 0
        (the second-segment cumulative is unbounded).
    telf : float
        Transition time (end of linear flow).
    b1 : float
        First-segment Arps b-factor (default 2.0, transient flow).
    b2 : float
        Second-segment Arps b-factor (default 0.5). Must not exceed b1.
    dterm : float
        Terminal nominal decline rate (1/time) applied to the second segment.
        Default 0 disables it.

    Returns
    -------
    float
        EUR (cumulative production when rate reaches q_min).
    """
    q_sw, D_sw, N_sw = _hyp2_switch(qi, di, b1, b2, telf)
    if q_min >= q_sw:
        # Abandonment reached within the first segment
        return eur(qi, di, b1, q_min)
    return N_sw + mh_eur(q_sw, D_sw, q_min, b=b2, dterm=dterm)


def _d1_from_stage1_cum(qi, q2, b1, cum1):
    """First-segment nominal decline that produces cumulative cum1 while the
    rate falls from qi to q2 (inversion of the q-form Arps cumulative)."""
    if b1 == 0:
        return (qi - q2) / cum1
    if b1 == 1:
        return qi * np.log(qi / q2) / cum1
    return qi ** b1 * (qi ** (1.0 - b1) - q2 ** (1.0 - b1)) / (cum1 * (1.0 - b1))


def _telf_from_switch_rate(qi, q2, b1, d1):
    """Time at which the first segment's rate falls to q2."""
    if b1 == 0:
        return np.log(qi / q2) / d1
    return ((qi / q2) ** b1 - 1.0) / (b1 * d1)


def _solve_d1(f):
    """Root of f(d1) = EUR(d1) - target. EUR decreases monotonically in d1,
    so bracket by expanding upward from 1.0 then apply brentq."""
    from scipy.optimize import brentq
    lo, hi = 1e-12, 1.0
    while f(hi) > 0.0 and hi < 1e15:
        hi *= 10.0
    if f(hi) > 0.0:
        raise RuntimeError("hyp2_from_eur: could not bracket a solution for di")
    return brentq(f, lo, hi)


def hyp2_from_eur(eur_target: float, qi: float, q_min: float, b1: float = 2.0,
                  b2: float = 0.5, telf: Optional[float] = None,
                  d_sw: Optional[float] = None,
                  eur_frac: Optional[float] = None,
                  rate_frac: Optional[float] = None) -> 'DeclineResult':
    """Solve two-segment hyperbolic decline parameters for a target EUR.

    Iterates the first-segment nominal decline di so that the two-segment
    hyperbolic (see hyp2_rate: rate and log-slope continuous at the
    transition) delivers eur_target when the rate reaches q_min. The
    transition is specified by exactly one of:

    - telf: transition at a fixed time
    - d_sw: transition when the nominal decline falls to d_sw (1/time)
    - eur_frac: transition after this fraction of eur_target is produced
    - rate_frac: transition when the rate falls to rate_frac * qi

    Port of the two-stage type-curve generator logic (Burgoyne 2016-2018,
    'Simple DCA Generator-2Stage'), with bisection replaced by brentq.

    Parameters
    ----------
    eur_target : float
        Target EUR (volume) at abandonment rate q_min.
    qi : float
        Initial rate (volume/time).
    q_min : float
        Abandonment rate. Must satisfy 0 < q_min < qi.
    b1 : float
        First-segment Arps b-factor (default 2.0, transient flow). Must be > 0.
    b2 : float
        Second-segment Arps b-factor (default 0.5). Must not exceed b1.
    telf, d_sw, eur_frac, rate_frac : float, optional
        Transition specification - supply exactly one.

    Returns
    -------
    DeclineResult
        With method='hyp2', the solved di, b=b1, b2 and resolved telf; works
        directly with forecast(), hyp2_rate() and hyp2_eur(). For the d_sw
        specification, if the solved decline never reaches d_sw before
        abandonment the result is a single-segment fit (method='hyperbolic').
    """
    if eur_target <= 0:
        raise ValueError(f"eur_target must be positive, got {eur_target}")
    if qi <= 0:
        raise ValueError(f"qi must be positive, got {qi}")
    if not 0 < q_min < qi:
        raise ValueError(f"q_min must be in (0, qi), got {q_min}")
    if b1 <= 0:
        raise ValueError(f"b1 must be positive, got {b1}")
    if b2 < 0 or b2 > b1:
        raise ValueError(f"b2 must be in [0, b1], got {b2}")
    specs = [s is not None for s in (telf, d_sw, eur_frac, rate_frac)]
    if sum(specs) != 1:
        raise ValueError("Specify exactly one of telf, d_sw, eur_frac, rate_frac")

    if telf is not None:
        if telf <= 0:
            raise ValueError(f"telf must be positive, got {telf}")
        d1 = _solve_d1(lambda d: hyp2_eur(qi, d, q_min, telf, b1, b2) - eur_target)
        return DeclineResult(method='hyp2', qi=qi, di=d1, b=b1, b2=b2, telf=telf)

    if d_sw is not None:
        if d_sw <= 0:
            raise ValueError(f"d_sw must be positive, got {d_sw}")

        def eur_at(d1):
            if d1 <= d_sw:
                # Decline starts at or below d_sw: single segment throughout
                return eur(qi, d1, b1, q_min)
            t_sw = (d1 / d_sw - 1.0) / (b1 * d1)
            return hyp2_eur(qi, d1, q_min, t_sw, b1, b2)

        d1 = _solve_d1(lambda d: eur_at(d) - eur_target)
        if d1 <= d_sw:
            return DeclineResult(method='hyperbolic', qi=qi, di=d1, b=b1)
        return DeclineResult(method='hyp2', qi=qi, di=d1, b=b1, b2=b2,
                             telf=(d1 / d_sw - 1.0) / (b1 * d1))

    if rate_frac is not None:
        if not 0 < rate_frac < 1:
            raise ValueError(f"rate_frac must be in (0, 1), got {rate_frac}")
        q2 = qi * rate_frac
        if q2 <= q_min:
            raise ValueError(
                f"Switch rate qi*rate_frac ({q2}) must exceed q_min ({q_min})")
        # With the switch rate fixed, both segment cumulatives scale as 1/di,
        # so EUR(di) = EUR(di=1)/di and the solution is direct.
        d1 = hyp2_eur(qi, 1.0, q_min, _telf_from_switch_rate(qi, q2, b1, 1.0),
                      b1, b2) / eur_target
        return DeclineResult(method='hyp2', qi=qi, di=d1, b=b1, b2=b2,
                             telf=_telf_from_switch_rate(qi, q2, b1, d1))

    # eur_frac specification: first-stage cumulative is fixed; root-find the
    # switch rate q2, with di recovered from the stage-1 cumulative inversion.
    from scipy.optimize import brentq
    if not 0 < eur_frac < 1:
        raise ValueError(f"eur_frac must be in (0, 1), got {eur_frac}")
    cum1 = eur_target * eur_frac

    def eur_err(q2):
        d1 = _d1_from_stage1_cum(qi, q2, b1, cum1)
        d2 = d1 * (q2 / qi) ** b1  # nominal decline at the switch (slope match)
        return cum1 + eur(q2, d2, b2, q_min) - eur_target

    q2 = brentq(eur_err, q_min * (1.0 + 1e-9), qi * (1.0 - 1e-9))
    d1 = _d1_from_stage1_cum(qi, q2, b1, cum1)
    return DeclineResult(method='hyp2', qi=qi, di=d1, b=b1, b2=b2,
                         telf=_telf_from_switch_rate(qi, q2, b1, d1))


def _fit_exponential(t, q):
    """Fit exponential decline via log-linear RANSAC regression."""
    ln_q = np.log(q)
    slope, intercept, _ = ransac_linreg(t, ln_q)
    di = -slope
    qi = np.exp(intercept)
    if di <= 0:
        return None
    q_pred = qi * np.exp(-di * t)
    return _build_decline_result('exponential', q, q_pred, qi=qi, di=di, b=0.0)


def _fit_harmonic(t, q):
    """Fit harmonic decline via 1/q RANSAC linear regression."""
    inv_q = 1.0 / q
    slope, intercept, _ = ransac_linreg(t, inv_q)
    # 1/q = 1/qi + di/qi * t  =>  slope = di/qi, intercept = 1/qi
    if intercept <= 0:
        return None
    qi = 1.0 / intercept
    di = slope * qi
    if di <= 0:
        return None
    q_pred = qi / (1.0 + di * t)
    return _build_decline_result('harmonic', q, q_pred, qi=qi, di=di, b=1.0)


def _fit_hyperbolic(t, q):
    """Fit hyperbolic decline via linearized grid search over b + RANSAC.

    For a given b, q^(-b) = qi^(-b) + qi^(-b)*b*Di * t is linear in t.
    Grid search over b finds the best R-squared, recovering qi and di
    algebraically from intercept and slope.
    """
    if _RUST_AVAILABLE:
        try:
            qi, di, b, r2 = _rust.fit_hyperbolic_rust(t.tolist(), q.tolist())
            q_pred = qi / (1.0 + b * di * t) ** (1.0 / b)
            return _build_decline_result('hyperbolic', q, q_pred, qi=qi, di=di, b=b)
        except (ImportError, AttributeError):
            pass

    best_r2 = -np.inf
    best_result = None

    for b_trial in np.arange(_B_GRID_MIN, _B_GRID_MAX, _B_GRID_STEP):
        # Transform: Y = q^(-b) is linear in t
        Y = q ** (-b_trial)
        slope, intercept, _ = ransac_linreg(t, Y)

        # intercept = qi^(-b), slope = qi^(-b) * b * di
        if intercept <= 0:
            continue
        qi = intercept ** (-1.0 / b_trial)
        if qi <= 0:
            continue
        di = slope / (intercept * b_trial)
        if di <= 0:
            continue

        # Compute R-squared in original rate space
        q_pred = qi / (1.0 + b_trial * di * t) ** (1.0 / b_trial)
        result = _build_decline_result('hyperbolic', q, q_pred, qi=qi, di=di, b=b_trial)
        if result.r_squared > best_r2:
            best_r2 = result.r_squared
            best_result = result

    return best_result


def _fit_duong(t, q):
    """Fit Duong decline via scipy curve_fit."""
    from scipy.optimize import curve_fit

    def duong_func(t, qi, a, m):
        return qi * t ** (-m) * np.exp(a / (1.0 - m) * (t ** (1.0 - m) - 1.0))

    # Filter t > 0
    mask = t > 0
    if np.sum(mask) < 3:
        return None
    t_f, q_f = t[mask], q[mask]

    try:
        popt, _ = curve_fit(
            duong_func, t_f, q_f,
            p0=[q_f[0], _DUONG_P0_A, _DUONG_P0_M],
            bounds=(
                list(_DUONG_BOUNDS_LO),
                [q_f[0] * _DUONG_BOUNDS_HI_QI_FACTOR,
                 _DUONG_BOUNDS_HI_A, _DUONG_BOUNDS_HI_M],
            ),
            maxfev=_CURVE_FIT_MAXFEV,
        )
        qi, a, m = popt
        q_pred_valid = duong_func(t_f, qi, a, m)
        return _build_decline_result('duong', q_f, q_pred_valid, qi=qi, a=a, m=m)
    except (RuntimeError, ValueError):
        return None


def _fit_mh(t, q):
    """Fit modified hyperbolic decline via scipy curve_fit.

    dterm is reparameterised as f * di (0 < f < 1) so the terminal decline
    stays below the initial decline for every trial parameter set. The di
    initial guess comes from a robust log-rate slope, as in petbox-dca.
    """
    from scipy.optimize import curve_fit

    def mh_func(t, qi, di, b, f):
        return np.asarray(mh_rate(qi, di, t, b=b, dterm=f * di), dtype=float)

    ln_q = np.log(q)
    slope, _, _ = ransac_linreg(t, ln_q)
    d0 = max(-slope, _MH_FIT_DI_GUESS_MIN)
    q0 = q[0]

    try:
        popt, _ = curve_fit(
            mh_func, t, q,
            p0=[q0, d0, _MH_P0_B, _MH_P0_F],
            bounds=(
                [q0 / _MH_FIT_QI_FACTOR, d0 / _MH_FIT_DI_SPAN,
                 _MH_FIT_B_MIN, _MH_FIT_F_MIN],
                [q0 * _MH_FIT_QI_FACTOR, d0 * _MH_FIT_DI_SPAN,
                 _MH_FIT_B_MAX, _MH_FIT_F_MAX],
            ),
            maxfev=_CURVE_FIT_MAXFEV,
        )
        qi, di, b, f = popt
        q_pred = mh_func(t, *popt)
        return _build_decline_result('mh', q, q_pred, qi=qi, di=di, b=b,
                                     dterm=f * di)
    except (RuntimeError, ValueError):
        return None


def _fit_hyp2(t, q):
    """Fit two-segment hyperbolic decline via scipy curve_fit.

    b2 is reparameterised as r * b1 (0 < r <= 1) so the second-segment
    b-factor can never exceed the first. Initial guesses follow the
    petbox-dca THM fitter (b1 = 2, bf/bi = 0.5, telf at the first third of
    the data).
    """
    from scipy.optimize import curve_fit

    def hyp2_func(t, qi, di, b1, r, telf):
        return np.asarray(hyp2_rate(qi, di, t, telf, b1=b1, b2=r * b1),
                          dtype=float)

    ln_q = np.log(q)
    slope, _, _ = ransac_linreg(t, ln_q)
    d0 = max(-slope, _MH_FIT_DI_GUESS_MIN)
    q0 = q[0]
    telf0 = t[len(t) // _HYP2_P0_TELF_IDX]

    try:
        popt, _ = curve_fit(
            hyp2_func, t, q,
            p0=[q0, d0, _HYP2_P0_B1, _HYP2_P0_R, telf0],
            bounds=(
                [q0 / _MH_FIT_QI_FACTOR, d0 / _MH_FIT_DI_SPAN,
                 _MH_FIT_B_MIN, _HYP2_FIT_R_MIN, 0.0],
                [q0 * _MH_FIT_QI_FACTOR, d0 * _MH_FIT_DI_SPAN,
                 _MH_FIT_B_MAX, 1.0, t[-1]],
            ),
            maxfev=_CURVE_FIT_MAXFEV,
        )
        qi, di, b1, r, telf = popt
        q_pred = hyp2_func(t, *popt)
        return _build_decline_result('hyp2', q, q_pred, qi=qi, di=di, b=b1,
                                     b2=r * b1, telf=telf)
    except (RuntimeError, ValueError):
        return None


def _fit_exponential_cum(Np, q):
    """Fit exponential decline via q-vs-Np RANSAC linear regression.

    q = qi - di * Np  =>  slope = -di, intercept = qi
    """
    slope, intercept, _ = ransac_linreg(Np, q)
    di = -slope
    qi = intercept
    if di <= 0 or qi <= 0:
        return None
    q_pred = qi - di * Np
    return _build_decline_result('exponential', q, q_pred, qi=qi, di=di, b=0.0)


def _fit_harmonic_cum(Np, q):
    """Fit harmonic decline via log-linear q-vs-Np RANSAC regression.

    q = qi * exp(-di/qi * Np)  =>  ln(q) = ln(qi) - (di/qi) * Np
    """
    ln_q = np.log(q)
    slope, intercept, _ = ransac_linreg(Np, ln_q)
    # slope = -di/qi, intercept = ln(qi)
    qi = np.exp(intercept)
    di = -slope * qi
    if di <= 0 or qi <= 0:
        return None
    q_pred = qi * np.exp(-di / qi * Np)
    return _build_decline_result('harmonic', q, q_pred, qi=qi, di=di, b=1.0)


def _fit_hyperbolic_cum(Np, q):
    """Fit hyperbolic decline via linearized grid search over b + RANSAC.

    For a given b, Np = A + B * q^(1-b) where:
      A = qi / ((1-b)*di)
      B = -qi^b / ((1-b)*di)
    Recover qi and di algebraically from regression coefficients.
    """
    if _RUST_AVAILABLE:
        try:
            qi, di, b, r2 = _rust.fit_hyperbolic_cum_rust(Np.tolist(), q.tolist())
            exp = 1.0 - b
            inner = np.maximum(1.0 - exp * di * Np / qi, _HYPER_INNER_FLOOR)
            q_pred = qi * inner ** (1.0 / exp)
            return _build_decline_result('hyperbolic', q, q_pred, qi=qi, di=di, b=b)
        except (ImportError, AttributeError):
            pass

    best_r2 = -np.inf
    best_result = None

    for b_trial in np.arange(_B_GRID_MIN, _B_GRID_MAX, _B_GRID_STEP):
        exp = 1.0 - b_trial
        # Transform: Np = A + B * q^(1-b)
        X = q ** exp
        slope, intercept, _ = ransac_linreg(X, Np)

        # A = intercept = qi / ((1-b)*di)
        # B = slope = -qi^b / ((1-b)*di)
        # B/A = -qi^b / qi = -qi^(b-1) = -1/qi^(1-b)
        if abs(slope) < _ZERO_DIV_EPS or abs(intercept) < _ZERO_DIV_EPS:
            continue
        ratio = slope / intercept  # = -1/qi^(1-b)
        if ratio >= 0:
            continue  # Need ratio < 0
        qi_exp = -1.0 / ratio  # qi^(1-b)
        if qi_exp <= 0:
            continue
        qi = qi_exp ** (1.0 / exp)
        if qi <= 0:
            continue
        di = qi / (exp * intercept)
        if di <= 0:
            continue

        # Compute R-squared in original rate space
        inner = np.maximum(1.0 - exp * di * Np / qi, _HYPER_INNER_FLOOR)
        q_pred = qi * inner ** (1.0 / exp)
        result = _build_decline_result('hyperbolic', q, q_pred, qi=qi, di=di, b=b_trial)
        if result.r_squared > best_r2:
            best_r2 = result.r_squared
            best_result = result

    return best_result


def fit_decline_cum(Np: ArrayLike, q: ArrayLike, method: str = 'best',
                    t_calendar: Optional[ArrayLike] = None,
                    Np_start: Optional[float] = None,
                    Np_end: Optional[float] = None) -> 'DeclineResult':
    """Fit a decline model to rate-vs-cumulative data.

    Eliminates time from the Arps equations to fit q as a function of Np.
    The returned qi and di are identical to time-domain parameters, so the
    result works directly with arps_rate() and forecast().

    Parameters
    ----------
    Np : array-like
        Cumulative production array.
    q : array-like
        Rate array (must be > 0).
    method : str
        'exponential', 'harmonic', 'hyperbolic', or 'best' (default).
        Passing 'duong' raises ValueError (no analytical q-vs-Np form);
        'mh' is not supported in the cumulative domain.
    t_calendar : array-like, optional
        Calendar time stamps corresponding to each data point. When provided,
        uptime fractions are inferred by comparing calendar-average rates to
        fitted capacity rates.
    Np_start : float, optional
        Start of fitting window on cumulative axis (inclusive).
    Np_end : float, optional
        End of fitting window on cumulative axis (inclusive).

    Returns
    -------
    DeclineResult
        With uptime_mean and uptime_history populated when t_calendar is given.
    """
    Np = np.asarray(Np, dtype=float)
    q = np.asarray(q, dtype=float)

    if Np_start is not None or Np_end is not None:
        mask = np.ones(len(Np), dtype=bool)
        if Np_start is not None:
            mask &= (Np >= Np_start)
        if Np_end is not None:
            mask &= (Np <= Np_end)
        Np = Np[mask]
        q = q[mask]
        if t_calendar is not None:
            t_calendar = np.asarray(t_calendar, dtype=float)[mask]
        if len(Np) == 0:
            raise ValueError("No data points within the specified Np_start/Np_end window")
        Np = Np - Np[0]  # Shift so window starts at Np=0

    if len(Np) != len(q):
        raise ValueError(f"Np and q must have same length, got {len(Np)} and {len(q)}")
    if len(Np) < 3:
        raise ValueError("Need at least 3 data points for fitting")
    if np.any(q <= 0):
        raise ValueError("All rate values must be positive")

    if method == 'duong':
        raise ValueError("Duong model has no analytical rate-vs-cumulative form")

    fitters = {
        'exponential': _fit_exponential_cum,
        'harmonic': _fit_harmonic_cum,
        'hyperbolic': _fit_hyperbolic_cum,
    }

    if method == 'best':
        results = []
        for name, fitter in fitters.items():
            result = fitter(Np, q)
            if result is not None:
                results.append(result)
        if not results:
            raise RuntimeError("All cumulative decline fitting methods failed")
        best = max(results, key=lambda r: r.r_squared)
    else:
        if method not in fitters:
            raise ValueError(f"Unknown method '{method}'. Choose from: {list(fitters.keys())} or 'best'")
        best = fitters[method](Np, q)
        if best is None:
            raise RuntimeError(f"Cumulative decline fitting with method '{method}' failed")

    # Uptime inference
    if t_calendar is not None:
        t_calendar = np.asarray(t_calendar, dtype=float)
        if len(t_calendar) != len(q):
            raise ValueError(f"t_calendar must have same length as q, got {len(t_calendar)} and {len(q)}")
        n = len(q)
        uptime_arr = np.zeros(n - 1)
        for i in range(n - 1):
            dt_cal = t_calendar[i + 1] - t_calendar[i]
            if dt_cal <= 0:
                uptime_arr[i] = 1.0
                continue
            dNp = Np[i + 1] - Np[i]
            calendar_rate = dNp / dt_cal
            capacity_rate = (q[i] + q[i + 1]) / 2.0
            if capacity_rate <= 0:
                uptime_arr[i] = 1.0
            else:
                uptime_arr[i] = np.clip(calendar_rate / capacity_rate, 0.0, 1.0)
        best.uptime_mean = float(np.mean(uptime_arr))
        best.uptime_history = uptime_arr

    return best


def fit_decline(t: ArrayLike, q: ArrayLike, method: str = 'best',
                t_start: Optional[float] = None,
                t_end: Optional[float] = None) -> 'DeclineResult':
    """Fit a decline model to production data.

    Parameters
    ----------
    t : array-like
        Time array.
    q : array-like
        Rate array (must be > 0).
    method : str
        'exponential', 'harmonic', 'hyperbolic', 'duong', 'mh' (modified
        hyperbolic, hyperbolic-to-exponential), 'hyp2' (two-segment
        hyperbolic), or 'best' (default). 'best' tries exponential, harmonic,
        hyperbolic and Duong and returns the one with highest R-squared;
        'mh' and 'hyp2' are explicit opt-in only, since their extra
        parameters would dominate the R-squared comparison.
    t_start : float, optional
        Start of fitting window (inclusive). Data before t_start is excluded.
    t_end : float, optional
        End of fitting window (inclusive). Data after t_end is excluded.

    Returns
    -------
    DeclineResult
    """
    t = np.asarray(t, dtype=float)
    q = np.asarray(q, dtype=float)

    if t_start is not None or t_end is not None:
        mask = np.ones(len(t), dtype=bool)
        if t_start is not None:
            mask &= (t >= t_start)
        if t_end is not None:
            mask &= (t <= t_end)
        t = t[mask]
        q = q[mask]
        if len(t) == 0:
            raise ValueError("No data points within the specified t_start/t_end window")
        t = t - t[0]  # Shift so window starts at t=0

    if len(t) != len(q):
        raise ValueError(f"t and q must have same length, got {len(t)} and {len(q)}")
    if len(t) < 3:
        raise ValueError("Need at least 3 data points for fitting")
    if np.any(q <= 0):
        raise ValueError("All rate values must be positive")

    fitters = {
        'exponential': _fit_exponential,
        'harmonic': _fit_harmonic,
        'hyperbolic': _fit_hyperbolic,
        'duong': _fit_duong,
        'mh': _fit_mh,
        'hyp2': _fit_hyp2,
    }

    if method == 'best':
        results = []
        for name, fitter in fitters.items():
            if name in ('mh', 'hyp2'):
                continue  # explicit opt-in only; see docstring
            result = fitter(t, q)
            if result is not None:
                results.append(result)
        if not results:
            raise RuntimeError("All decline fitting methods failed")
        # The Duong fitter drops t <= 0 points, so its stored r_squared is
        # computed on a different subset to the Arps fitters. Compare all
        # methods on the common t > 0 subset so 'best' is apples-to-apples.
        # Fitted parameters and stored r_squared values are unchanged.
        if np.any(t <= 0):
            t_c = t[t > 0]
            q_c = q[t > 0]

            def _common_r2(r):
                if r.method == 'duong':
                    q_pred = np.asarray(duong_rate(r.qi, r.a, r.m, t_c), dtype=float)
                else:
                    q_pred = np.asarray(arps_rate(r.qi, r.di, r.b, t_c), dtype=float)
                ss_res = np.sum((q_c - q_pred) ** 2)
                ss_tot = np.sum((q_c - np.mean(q_c)) ** 2)
                return 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

            return max(results, key=_common_r2)
        return max(results, key=lambda r: r.r_squared)
    else:
        if method not in fitters:
            raise ValueError(f"Unknown method '{method}'. Choose from: {list(fitters.keys())} or 'best'")
        result = fitters[method](t, q)
        if result is None:
            raise RuntimeError(f"Decline fitting with method '{method}' failed")
        return result


def _fit_ratio_linear(x, ratio):
    """Fit R = a + b*x via RANSAC linear regression."""
    b, a, _ = ransac_linreg(x, ratio)
    r_pred = a + b * x
    ss_res = np.sum((ratio - r_pred) ** 2)
    ss_tot = np.sum((ratio - np.mean(ratio)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return RatioResult(
        method='linear', a=a, b=b,
        r_squared=r2, residuals=ratio - r_pred,
    )


def _fit_ratio_exponential(x, ratio):
    """Fit R = a * exp(b*x) via log-linear RANSAC regression."""
    if np.any(ratio <= 0):
        return None
    ln_r = np.log(ratio)
    b, ln_a, _ = ransac_linreg(x, ln_r)
    a = np.exp(ln_a)
    r_pred = a * np.exp(b * x)
    ss_res = np.sum((ratio - r_pred) ** 2)
    ss_tot = np.sum((ratio - np.mean(ratio)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return RatioResult(
        method='exponential', a=a, b=b,
        r_squared=r2, residuals=ratio - r_pred,
    )


def _fit_ratio_power(x, ratio):
    """Fit R = a * x^b via log-log RANSAC regression."""
    if np.any(x <= 0) or np.any(ratio <= 0):
        return None
    ln_x = np.log(x)
    ln_r = np.log(ratio)
    b, ln_a, _ = ransac_linreg(ln_x, ln_r)
    a = np.exp(ln_a)
    r_pred = a * x ** b
    ss_res = np.sum((ratio - r_pred) ** 2)
    ss_tot = np.sum((ratio - np.mean(ratio)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return RatioResult(
        method='power', a=a, b=b,
        r_squared=r2, residuals=ratio - r_pred,
    )


def _fit_ratio_logistic(x, ratio):
    """Fit R = Rmax / (1 + c * exp(-b*x)) via curve_fit."""
    from scipy.optimize import curve_fit

    def logistic_func(x, Rmax, b, c):
        return Rmax / (1.0 + c * np.exp(-b * x))

    try:
        Rmax_guess = np.max(ratio) * _LOGISTIC_P0_RMAX_INIT
        popt, _ = curve_fit(
            logistic_func, x, ratio,
            p0=[Rmax_guess, _LOGISTIC_P0_TC, _LOGISTIC_P0_ALPHA],
            bounds=(
                [0.0, _LOGISTIC_BOUNDS_LO[0], _LOGISTIC_BOUNDS_LO[1]],
                [Rmax_guess * _LOGISTIC_BOUNDS_RMAX_FACTOR,
                 _LOGISTIC_BOUNDS_HI, _LOGISTIC_BOUNDS_HI],
            ),
            maxfev=_CURVE_FIT_MAXFEV,
        )
        Rmax, b, c = popt
        r_pred = logistic_func(x, Rmax, b, c)
        ss_res = np.sum((ratio - r_pred) ** 2)
        ss_tot = np.sum((ratio - np.mean(ratio)) ** 2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        return RatioResult(
            method='logistic', a=Rmax, b=b, c=c,
            r_squared=r2, residuals=ratio - r_pred,
        )
    except (RuntimeError, ValueError):
        return None


def fit_ratio(x: ArrayLike, ratio: ArrayLike,
              method: str = 'best', domain: str = 'cum') -> 'RatioResult':
    """Fit a ratio model (e.g. GOR, WOR) to data.

    Parameters
    ----------
    x : array-like
        Independent variable (cumulative production or time).
    ratio : array-like
        Ratio values (e.g. GOR, WOR). Must be > 0 for exponential/power.
    method : str
        'linear', 'exponential', 'power', 'logistic', or 'best' (default).
    domain : str
        'cum' or 'time' — stored in result for use by forecast().

    Returns
    -------
    RatioResult
    """
    x = np.asarray(x, dtype=float)
    ratio = np.asarray(ratio, dtype=float)

    if len(x) != len(ratio):
        raise ValueError(f"x and ratio must have same length, got {len(x)} and {len(ratio)}")
    if len(x) < 3:
        raise ValueError("Need at least 3 data points for fitting")

    fitters = {
        'linear': _fit_ratio_linear,
        'exponential': _fit_ratio_exponential,
        'power': _fit_ratio_power,
        'logistic': _fit_ratio_logistic,
    }

    if method == 'best':
        results = []
        for name, fitter in fitters.items():
            result = fitter(x, ratio)
            if result is not None:
                results.append(result)
        if not results:
            raise RuntimeError("All ratio fitting methods failed")
        best = max(results, key=lambda r: r.r_squared)
        best.domain = domain
        return best
    else:
        if method not in fitters:
            raise ValueError(f"Unknown method '{method}'. Choose from: {list(fitters.keys())} or 'best'")
        result = fitters[method](x, ratio)
        if result is None:
            raise RuntimeError(f"Ratio fitting with method '{method}' failed")
        result.domain = domain
        return result


def ratio_forecast(result: 'RatioResult', x: ArrayLike) -> np.ndarray:
    """Evaluate a fitted ratio model at given x values.

    Parameters
    ----------
    result : RatioResult
        Fitted ratio result from fit_ratio().
    x : float or array-like
        Independent variable values to evaluate at.

    Returns
    -------
    float or np.ndarray
        Ratio values at x.
    """
    x, is_list = convert_to_numpy(x)

    if result.method == 'linear':
        r = result.a + result.b * x
    elif result.method == 'exponential':
        r = result.a * np.exp(result.b * x)
    elif result.method == 'power':
        r = result.a * x ** result.b
    elif result.method == 'logistic':
        r = result.a / (1.0 + result.c * np.exp(-result.b * x))
    else:
        raise ValueError(f"Unknown ratio method '{result.method}'")

    return process_output(r, is_list)


def forecast(result: 'DeclineResult', t_end: float, dt: float = 1.0,
             q_min: float = 0.0, uptime: float = 1.0,
             ratios: Optional[dict] = None) -> 'ForecastResult':
    """Generate a rate and cumulative forecast from a fitted decline model.

    Parameters
    ----------
    result : DeclineResult
        Fitted decline result from fit_decline() or fit_decline_cum().
    t_end : float
        End time for the forecast.
    dt : float
        Time step (default 1.0).
    q_min : float
        Economic limit rate (default 0, no cutoff).
    uptime : float
        Uptime fraction (0 to 1, default 1.0). Calendar-effective rate is
        capacity rate multiplied by uptime. Default 1.0 preserves existing
        behavior.
    ratios : dict, optional
        Maps names to RatioResult objects. For each entry, secondary phase
        rate and cumulative are computed. Domain='cum' evaluates ratio against
        cumulative production; domain='time' evaluates against time.

    Returns
    -------
    ForecastResult
    """
    if dt <= 0:
        raise ValueError(f"dt must be positive, got {dt}")
    if t_end <= 0:
        raise ValueError(f"t_end must be positive, got {t_end}")
    if not (0.0 < uptime <= 1.0):
        raise ValueError(f"uptime must be in (0, 1], got {uptime}")
    t = np.arange(dt, t_end + dt / 2, dt)

    if result.method == 'duong':
        q_capacity = np.asarray(duong_rate(result.qi, result.a, result.m, t), dtype=float)
    elif result.method == 'mh':
        q_capacity = np.asarray(mh_rate(result.qi, result.di, t, b=result.b,
                                        dterm=result.dterm), dtype=float)
    elif result.method == 'hyp2':
        q_capacity = np.asarray(hyp2_rate(result.qi, result.di, t, result.telf,
                                          b1=result.b, b2=result.b2,
                                          dterm=result.dterm), dtype=float)
    else:
        q_capacity = np.asarray(arps_rate(result.qi, result.di, result.b, t), dtype=float)

    # Apply uptime
    q = q_capacity * uptime

    # Apply economic limit
    if q_min > 0:
        mask = q >= q_min
        if not np.any(mask):
            t = t[:1]
            q = q[:1]
            q_capacity = q_capacity[:1]
        else:
            last_idx = np.where(mask)[0][-1] + 1
            t = t[:last_idx]
            q = q[:last_idx]
            q_capacity = q_capacity[:last_idx]

    # Analytic cumulative (right-rectangle cumsum biases Qcum/EUR low).
    # Uptime scales cumulative exactly as it scales rate.
    if result.method == 'duong':
        Qcum = np.asarray(duong_cum(result.qi, result.a, result.m, t), dtype=float) * uptime
    elif result.method == 'mh':
        Qcum = np.asarray(mh_cum(result.qi, result.di, t, b=result.b,
                                 dterm=result.dterm), dtype=float) * uptime
    elif result.method == 'hyp2':
        Qcum = np.asarray(hyp2_cum(result.qi, result.di, t, result.telf,
                                   b1=result.b, b2=result.b2,
                                   dterm=result.dterm), dtype=float) * uptime
    else:
        Qcum = np.asarray(arps_cum(result.qi, result.di, result.b, t), dtype=float) * uptime

    # Secondary phase ratios
    secondary = None
    if ratios is not None:
        secondary = {}
        for name, rr in ratios.items():
            if rr.domain == 'cum':
                x_eval = Qcum
            else:
                x_eval = t
            R = np.asarray(ratio_forecast(rr, x_eval), dtype=float)
            sec_rate = q * R
            sec_cum = np.cumsum(sec_rate * dt)
            secondary[name] = {'ratio': R, 'rate': sec_rate, 'cum': sec_cum}

    return ForecastResult(
        t=t, q=q, Qcum=Qcum,
        eur=float(Qcum[-1]) if len(Qcum) > 0 else 0.0,
        secondary=secondary,
    )
