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
"""

from enum import Enum
import numpy as np
import numpy.typing as npt
from typing import Union, List, Tuple

def bisect_solve(args, f, xmin, xmax, rtol):
    err_hi = f(args, xmax)
    err_lo = f(args, xmin)
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


def check_2_inputs(x: Union[float, List[float]], y: Union[float, List[float]]) -> bool:
    """ Check inputs that need to be matched, either both float or both lists/arrays of same length"""
    # Check if both are scalars
    if isinstance(x, (int, float)) and isinstance(y, (int, float)):
        return True
    # Check if both are arrays or lists of same length
    if isinstance(x, (list, np.ndarray)) and isinstance(y, (list, np.ndarray)):
        return len(x) == len(y)
    return False