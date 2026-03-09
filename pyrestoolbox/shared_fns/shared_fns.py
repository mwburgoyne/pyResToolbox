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

Shared utility functions used across pyResToolbox modules.

Functions
---------
bisect_solve        Generic bisection root-finder
convert_to_numpy    Normalize scalar/list/array inputs to numpy arrays
process_output      Convert results back to scalar or array based on input type
check_2_inputs      Validate matched float/list inputs
validate_pe_inputs  Centralized petroleum engineering input validation
halley_solve_cubic  Halley's method cubic equation solver
ransac_linreg       RANSAC-robust linear regression
"""

__all__ = [
    'bisect_solve', 'convert_to_numpy', 'process_output',
    'check_2_inputs', 'validate_pe_inputs', 'halley_solve_cubic',
    'ransac_linreg',
]

import numpy as np
from math import ceil
from typing import Union, List, Tuple

def bisect_solve(args, f, xmin, xmax, rtol):
    if xmin > xmax:
        xmin, xmax = xmax, xmin
    err_hi = f(args, xmax)
    err_lo = f(args, xmin)
    if err_hi * err_lo > 0:
        raise ValueError(f"bisect_solve: root is not bracketed between xmin={xmin} and xmax={xmax} (f(xmin)={err_lo}, f(xmax)={err_hi})")
    iternum = 0
    err_mid = 1e6
    while abs(err_mid) > rtol:
        mid_val = (xmax + xmin) / 2
        err_mid = f(args, mid_val)
        iternum += 1
        if iternum > 99:
            raise RuntimeError("bisect_solve: failed to converge after 99 iterations")
        if (err_hi * err_mid < 0):  # Solution point must be higher than current mid_val case
            xmin = mid_val
            err_lo = err_mid
            mid_val = (mid_val + xmax) / 2
        else:
            xmax = mid_val  # Otherwise must be lower than current mid_val case
            err_hi = err_mid
    return mid_val
    
def convert_to_numpy(input_data):
    # Convert input data to a numpy array ensuring it is always sizeable
    if isinstance(input_data, (np.ndarray, list, tuple)):
        return np.atleast_1d(input_data), True
    else:
        # Scalar input
        return np.atleast_1d(input_data), False

                
def process_output(input_data, is_list):
    # Check if input_data is a numpy array
    if isinstance(input_data, np.ndarray):
        if not is_list:
            # Return the single element if it's a single-element array
            return input_data.item()
        else:
            # Return the array itself if it's larger
            return input_data
    elif isinstance(input_data, list):
        if not is_list:
            # If it's a single-element list, return the element
            return input_data[0]
        else:
            # If it's a larger list, convert it to a numpy array and return
            return np.array(input_data)
    else:
        # If it's a single value (not in a list or array), return it directly
        return float(input_data)
        
def halley_solve_cubic(c2, c1, c0, flag=1, max_iter=50, tol=1e-12):
    """Solve Z^3 + c2*Z^2 + c1*Z + c0 = 0 via Halley's method.

    flag=1: return max real root (vapor)
    flag=-1: return min real root (liquid)
    flag=0: return all real roots as array
    Returns None if not converged (caller should use fallback).
    """
    # Find one root via Halley iteration from inflection point
    Z = -c2 / 3.0
    f = Z**3 + c2 * Z**2 + c1 * Z + c0
    if f < 0:
        Z = Z + 1.0  # Start on vapor side

    converged = False
    for _ in range(max_iter):
        f = Z**3 + c2 * Z**2 + c1 * Z + c0
        fp = 3.0 * Z**2 + 2.0 * c2 * Z + c1
        fpp = 6.0 * Z + 2.0 * c2
        if abs(fp) < 1e-30:
            break
        dZ = f / fp
        denom = fp - 0.5 * dZ * fpp
        if abs(denom) < 1e-30:
            break
        dZ = f / denom
        Z -= dZ
        if abs(dZ) < tol:
            converged = True
            break

    if not converged:
        return None

    # Check residual
    f = Z**3 + c2 * Z**2 + c1 * Z + c0
    if abs(f) > 1e-6:
        return None

    r1 = Z

    # Factor out converged root via synthetic division: Z^2 + e1*Z + e0
    e1 = c2 + r1
    e0 = c1 + r1 * e1

    roots = [r1]

    # Solve depressed quadratic
    disc = e1**2 - 4.0 * e0
    if disc >= 0:
        sqrt_disc = np.sqrt(disc)
        r2 = (-e1 + sqrt_disc) / 2.0
        r3 = (-e1 - sqrt_disc) / 2.0

        # Refine r2 with one Halley step (Michelsen pattern)
        for r in (r2, r3):
            f = r**3 + c2 * r**2 + c1 * r + c0
            fp = 3.0 * r**2 + 2.0 * c2 * r + c1
            fpp = 6.0 * r + 2.0 * c2
            if abs(fp) > 1e-30:
                dZ = f / fp
                denom = fp - 0.5 * dZ * fpp
                if abs(denom) > 1e-30:
                    r -= f / denom
            roots.append(r)

    roots = np.array(roots)

    if flag == 0:
        return roots
    elif flag == 1:
        return float(np.max(roots))
    elif flag == -1:
        return float(np.min(roots))
    return roots


def validate_pe_inputs(p=None, degf=None, sg=None, co2=None, h2s=None, n2=None, h2=None):
    """Validate common petroleum engineering inputs.
    Checks physical constraints: positive pressure, temperature above absolute zero,
    positive specific gravity, non-negative mole fractions summing to <= 1.0.
    Raises ValueError with descriptive message on failure.
    """
    if p is not None:
        p_arr = np.atleast_1d(p)
        if np.any(p_arr <= 0):
            raise ValueError(f"Pressure must be positive, got min value: {np.min(p_arr)}")
    if degf is not None:
        degf_arr = np.atleast_1d(degf)
        if np.any(degf_arr <= -459.67):
            raise ValueError(f"Temperature must be above absolute zero (-459.67 degF), got: {np.min(degf_arr)}")
    if sg is not None and sg <= 0:
        raise ValueError(f"Specific gravity must be positive, got: {sg}")
    fracs = {'co2': co2, 'h2s': h2s, 'n2': n2, 'h2': h2}
    frac_sum = 0
    for name, val in fracs.items():
        if val is not None:
            if val < 0:
                raise ValueError(f"Mole fraction {name} must be non-negative, got: {val}")
            frac_sum += val
    if frac_sum > 1.0:
        raise ValueError(f"Sum of non-hydrocarbon mole fractions ({frac_sum}) exceeds 1.0")


def check_2_inputs(x: Union[float, List[float]], y: Union[float, List[float]]) -> bool:
    """ Check inputs that need to be matched, either both float or both lists/arrays of same length"""
    # Check if both are scalars
    if isinstance(x, (int, float)) and isinstance(y, (int, float)):
        return True
    # Check if both are arrays or lists of same length
    if isinstance(x, (list, np.ndarray)) and isinstance(y, (list, np.ndarray)):
        return len(x) == len(y)
    return False


def ransac_linreg(x: np.ndarray, y: np.ndarray,
                  n_iter: int = 200, threshold_sigma: float = 3.0,
                  seed: int = 42,
                  through_origin: bool = False) -> Tuple[float, float, np.ndarray]:
    """RANSAC-robust linear regression.

    Fits y = slope * x + intercept (or y = slope * x if through_origin).
    Uses MAD-based inlier threshold for outlier robustness. With clean data
    the MAD threshold becomes very large, all points are inliers, and the
    result matches ordinary least squares exactly.

    Parameters
    ----------
    x : np.ndarray
        Independent variable.
    y : np.ndarray
        Dependent variable.
    n_iter : int
        Number of random sample iterations (default 200).
    threshold_sigma : float
        Inlier threshold in MAD-derived sigma units (default 3.0).
    seed : int
        Random seed for reproducibility (default 42).
    through_origin : bool
        If True, fit y = slope * x with zero intercept (default False).

    Returns
    -------
    slope : float
        Fitted slope.
    intercept : float
        Fitted intercept (0.0 if through_origin).
    inlier_mask : np.ndarray
        Boolean array identifying inlier points.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)

    if n < 2:
        raise ValueError("Need at least 2 data points for linear regression")

    min_samples = 1 if through_origin else 2
    min_inliers = max(min_samples, ceil(0.5 * n))

    def _fit_ols(x_s, y_s):
        """OLS fit on subset. Returns (slope, intercept)."""
        if through_origin:
            # slope = sum(x*y) / sum(x^2)
            denom = np.sum(x_s * x_s)
            if abs(denom) < 1e-30:
                return 0.0, 0.0
            return float(np.sum(x_s * y_s) / denom), 0.0
        else:
            n_s = len(x_s)
            sx = np.sum(x_s)
            sy = np.sum(y_s)
            sxx = np.sum(x_s * x_s)
            sxy = np.sum(x_s * y_s)
            denom = n_s * sxx - sx * sx
            if abs(denom) < 1e-30:
                return 0.0, float(np.mean(y_s))
            slope = (n_s * sxy - sx * sy) / denom
            intercept = (sy - slope * sx) / n_s
            return float(slope), float(intercept)

    # Compute MAD-based threshold from full-data OLS residuals
    slope_full, intercept_full = _fit_ols(x, y)
    residuals_full = y - (slope_full * x + intercept_full)
    mad = np.median(np.abs(residuals_full - np.median(residuals_full)))
    sigma_est = 1.4826 * mad  # MAD to std-dev for normal distribution
    if sigma_est < 1e-15:
        # Near-perfect fit (clean data) — all points are inliers
        return slope_full, intercept_full, np.ones(n, dtype=bool)

    threshold = threshold_sigma * sigma_est

    rng = np.random.RandomState(seed)
    best_n_inliers = 0
    best_inlier_mask = np.ones(n, dtype=bool)
    best_slope = slope_full
    best_intercept = intercept_full

    for _ in range(n_iter):
        idx = rng.choice(n, size=min_samples, replace=False)
        slope_t, intercept_t = _fit_ols(x[idx], y[idx])
        residuals_t = np.abs(y - (slope_t * x + intercept_t))
        inlier_mask = residuals_t < threshold
        n_inliers = int(np.sum(inlier_mask))

        if n_inliers > best_n_inliers and n_inliers >= min_inliers:
            # Refit on all inliers
            slope_r, intercept_r = _fit_ols(x[inlier_mask], y[inlier_mask])
            best_slope = slope_r
            best_intercept = intercept_r
            best_n_inliers = n_inliers
            best_inlier_mask = inlier_mask

    if best_n_inliers < min_inliers:
        # RANSAC couldn't find a good consensus — fall back to full OLS
        return slope_full, intercept_full, np.ones(n, dtype=bool)

    return best_slope, best_intercept, best_inlier_mask