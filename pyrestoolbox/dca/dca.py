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
    'arps_rate', 'arps_cum', 'duong_rate', 'duong_cum',
    'eur', 'fit_decline', 'fit_decline_cum', 'forecast',
    'fit_ratio', 'ratio_forecast',
    'DeclineResult', 'ForecastResult', 'RatioResult',
]

import numpy as np
from dataclasses import dataclass, field
from typing import Optional

from pyrestoolbox.shared_fns import convert_to_numpy, process_output, ransac_linreg
from pyrestoolbox._accelerator import RUST_AVAILABLE as _RUST_AVAILABLE
if _RUST_AVAILABLE:
    from pyrestoolbox import _native as _rust


@dataclass
class DeclineResult:
    """Result from decline curve fitting.

    Attributes
    ----------
    method : str
        Decline model ('exponential', 'harmonic', 'hyperbolic', 'duong').
    qi : float
        Initial rate.
    di : float
        Initial decline rate (1/time). Not used for Duong.
    b : float
        Arps b-factor. 0 for exponential, 1 for harmonic. Not used for Duong.
    a : float
        Duong 'a' parameter. 0 for Arps models.
    m : float
        Duong 'm' parameter. 0 for Arps models.
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
    r_squared: float = 0.0
    residuals: np.ndarray = field(default_factory=lambda: np.array([]))
    uptime_mean: Optional[float] = None
    uptime_history: Optional[np.ndarray] = None


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


def arps_rate(qi, di, b, t):
    """Arps decline rate.

    Parameters
    ----------
    qi : float or array
        Initial rate (volume/time).
    di : float
        Initial decline rate (1/time). Must be > 0.
    b : float
        Arps b-factor. b=0: exponential, 0<b<1: hyperbolic, b=1: harmonic.
    t : float or array
        Time.

    Returns
    -------
    float or np.ndarray
        Rate at time t.
    """
    t, is_list = convert_to_numpy(t)
    if di <= 0:
        raise ValueError(f"di must be positive, got {di}")
    if b < 0 or b > 1:
        raise ValueError(f"b must be between 0 and 1, got {b}")

    if b == 0:
        q = qi * np.exp(-di * t)
    elif b == 1:
        q = qi / (1.0 + di * t)
    else:
        q = qi / (1.0 + b * di * t) ** (1.0 / b)

    return process_output(q, is_list)


def arps_cum(qi, di, b, t):
    """Arps cumulative production.

    Parameters
    ----------
    qi : float or array
        Initial rate (volume/time).
    di : float
        Initial decline rate (1/time). Must be > 0.
    b : float
        Arps b-factor. b=0: exponential, 0<b<1: hyperbolic, b=1: harmonic.
    t : float or array
        Time.

    Returns
    -------
    float or np.ndarray
        Cumulative production at time t.
    """
    t, is_list = convert_to_numpy(t)
    if di <= 0:
        raise ValueError(f"di must be positive, got {di}")
    if b < 0 or b > 1:
        raise ValueError(f"b must be between 0 and 1, got {b}")

    if b == 0:
        Qcum = (qi / di) * (1.0 - np.exp(-di * t))
    elif b == 1:
        Qcum = (qi / di) * np.log(1.0 + di * t)
    else:
        Qcum = (qi / ((1.0 - b) * di)) * ((1.0 + b * di * t) ** (1.0 - 1.0 / b) - 1.0)
        # Equivalent to: qi^b / ((1-b)*di) * (qi^(1-b) - q(t)^(1-b))
        # Simplified: (qi / ((1-b)*di)) * (1 - (1+b*di*t)^((b-1)/b))
        # Using: (1+b*di*t)^(1-1/b) - 1 = (1+b*di*t)^((b-1)/b) - 1
        # Note sign: cumulative should be positive
        # Correct formula: Qcum = qi / ((1-b)*di) * (1 - (1+b*di*t)^(-(1-b)/b))
        # Let's recalculate properly
        Qcum = (qi / ((1.0 - b) * di)) * (1.0 - (1.0 + b * di * t) ** (-(1.0 - b) / b))

    return process_output(Qcum, is_list)


def duong_rate(qi, a, m, t):
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


def duong_cum(qi, a, m, t):
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

    # Generate fine time grid for integration
    results = np.zeros_like(t, dtype=float)
    for i, ti in enumerate(t):
        t_fine = np.linspace(0.001, ti, max(500, int(ti * 10)))
        q_fine = qi * t_fine ** (-m) * np.exp(a / (1.0 - m) * (t_fine ** (1.0 - m) - 1.0))
        results[i] = np.trapezoid(q_fine, t_fine)

    return process_output(results, is_list)


def eur(qi, di, b, q_min):
    """Estimated ultimate recovery for Arps decline.

    Parameters
    ----------
    qi : float
        Initial rate.
    di : float
        Initial decline rate (1/time).
    b : float
        Arps b-factor (0 to 1).
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
    if b < 0 or b > 1:
        raise ValueError(f"b must be between 0 and 1, got {b}")

    # Solve for time when q(t) = q_min
    if b == 0:
        t_end = -np.log(q_min / qi) / di
    elif b == 1:
        t_end = (qi / q_min - 1.0) / di
    else:
        t_end = ((qi / q_min) ** b - 1.0) / (b * di)

    return float(arps_cum(qi, di, b, t_end))


def _fit_exponential(t, q):
    """Fit exponential decline via log-linear RANSAC regression."""
    ln_q = np.log(q)
    slope, intercept, _ = ransac_linreg(t, ln_q)
    di = -slope
    qi = np.exp(intercept)
    if di <= 0:
        return None
    q_pred = qi * np.exp(-di * t)
    ss_res = np.sum((q - q_pred) ** 2)
    ss_tot = np.sum((q - np.mean(q)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return DeclineResult(
        method='exponential', qi=qi, di=di, b=0.0,
        r_squared=r2, residuals=q - q_pred,
    )


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
    ss_res = np.sum((q - q_pred) ** 2)
    ss_tot = np.sum((q - np.mean(q)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return DeclineResult(
        method='harmonic', qi=qi, di=di, b=1.0,
        r_squared=r2, residuals=q - q_pred,
    )


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
            return DeclineResult(
                method='hyperbolic', qi=qi, di=di, b=b,
                r_squared=r2, residuals=q - q_pred,
            )
        except Exception:
            pass

    best_r2 = -np.inf
    best_result = None

    for b_trial in np.arange(0.05, 0.96, 0.01):
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
        ss_res = np.sum((q - q_pred) ** 2)
        ss_tot = np.sum((q - np.mean(q)) ** 2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

        if r2 > best_r2:
            best_r2 = r2
            best_result = DeclineResult(
                method='hyperbolic', qi=qi, di=di, b=b_trial,
                r_squared=r2, residuals=q - q_pred,
            )

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
        popt, _ = curve_fit(duong_func, t_f, q_f, p0=[q_f[0], 1.0, 1.2],
                            bounds=([0, 0.01, 1.001], [q_f[0] * 5, 10.0, 3.0]),
                            maxfev=5000)
        qi, a, m = popt
        q_pred_full = np.full_like(q, np.nan, dtype=float)
        q_pred_full[mask] = duong_func(t_f, qi, a, m)
        q_pred_valid = q_pred_full[mask]
        ss_res = np.sum((q_f - q_pred_valid) ** 2)
        ss_tot = np.sum((q_f - np.mean(q_f)) ** 2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        return DeclineResult(
            method='duong', qi=qi, a=a, m=m,
            r_squared=r2, residuals=q_f - q_pred_valid,
        )
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
    ss_res = np.sum((q - q_pred) ** 2)
    ss_tot = np.sum((q - np.mean(q)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return DeclineResult(
        method='exponential', qi=qi, di=di, b=0.0,
        r_squared=r2, residuals=q - q_pred,
    )


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
    ss_res = np.sum((q - q_pred) ** 2)
    ss_tot = np.sum((q - np.mean(q)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return DeclineResult(
        method='harmonic', qi=qi, di=di, b=1.0,
        r_squared=r2, residuals=q - q_pred,
    )


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
            inner = np.maximum(1.0 - exp * di * Np / qi, 1e-10)
            q_pred = qi * inner ** (1.0 / exp)
            return DeclineResult(
                method='hyperbolic', qi=qi, di=di, b=b,
                r_squared=r2, residuals=q - q_pred,
            )
        except Exception:
            pass

    best_r2 = -np.inf
    best_result = None

    for b_trial in np.arange(0.05, 0.96, 0.01):
        exp = 1.0 - b_trial
        # Transform: Np = A + B * q^(1-b)
        X = q ** exp
        slope, intercept, _ = ransac_linreg(X, Np)

        # A = intercept = qi / ((1-b)*di)
        # B = slope = -qi^b / ((1-b)*di)
        # B/A = -qi^b / qi = -qi^(b-1) = -1/qi^(1-b)
        if abs(slope) < 1e-30 or abs(intercept) < 1e-30:
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
        inner = 1.0 - exp * di * Np / qi
        inner = np.maximum(inner, 1e-10)
        q_pred = qi * inner ** (1.0 / exp)
        ss_res = np.sum((q - q_pred) ** 2)
        ss_tot = np.sum((q - np.mean(q)) ** 2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

        if r2 > best_r2:
            best_r2 = r2
            best_result = DeclineResult(
                method='hyperbolic', qi=qi, di=di, b=b_trial,
                r_squared=r2, residuals=q - q_pred,
            )

    return best_result


def fit_decline_cum(Np, q, method='best', t_calendar=None, Np_start=None, Np_end=None):
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
        Passing 'duong' raises ValueError (no analytical q-vs-Np form).
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


def fit_decline(t, q, method='best', t_start=None, t_end=None):
    """Fit a decline model to production data.

    Parameters
    ----------
    t : array-like
        Time array.
    q : array-like
        Rate array (must be > 0).
    method : str
        'exponential', 'harmonic', 'hyperbolic', 'duong', or 'best' (default).
        'best' tries all four and returns the one with highest R-squared.
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
    }

    if method == 'best':
        results = []
        for name, fitter in fitters.items():
            result = fitter(t, q)
            if result is not None:
                results.append(result)
        if not results:
            raise RuntimeError("All decline fitting methods failed")
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
        Rmax_guess = np.max(ratio) * 1.5
        popt, _ = curve_fit(logistic_func, x, ratio,
                            p0=[Rmax_guess, 0.01, 10.0],
                            bounds=([0, 1e-8, 1e-3], [Rmax_guess * 5, 10.0, 1e6]),
                            maxfev=5000)
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


def fit_ratio(x, ratio, method='best', domain='cum'):
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


def ratio_forecast(result, x):
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


def forecast(result, t_end, dt=1.0, q_min=0.0, uptime=1.0, ratios=None):
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
    t = np.arange(dt, t_end + dt / 2, dt)

    if result.method == 'duong':
        q_capacity = np.array([float(duong_rate(result.qi, result.a, result.m, ti)) for ti in t])
    else:
        q_capacity = np.array([float(arps_rate(result.qi, result.di, result.b, ti)) for ti in t])

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

    Qcum = np.cumsum(q * dt)

    # Secondary phase ratios
    secondary = None
    if ratios is not None:
        secondary = {}
        for name, rr in ratios.items():
            if rr.domain == 'cum':
                x_eval = Qcum
            else:
                x_eval = t
            R = np.array([float(ratio_forecast(rr, xi)) for xi in x_eval])
            sec_rate = q * R
            sec_cum = np.cumsum(sec_rate * dt)
            secondary[name] = {'ratio': R, 'rate': sec_rate, 'cum': sec_cum}

    return ForecastResult(
        t=t, q=q, Qcum=Qcum,
        eur=float(Qcum[-1]) if len(Qcum) > 0 else 0.0,
        secondary=secondary,
    )
