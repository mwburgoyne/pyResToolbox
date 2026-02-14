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
        
def check_2_inputs(x: Union[float, List[float]], y: Union[float, List[float]]) -> bool:
    """ Check inputs that need to be matched, either both float or both lists/arrays of same length"""
    # Check if both are scalars
    if isinstance(x, (int, float)) and isinstance(y, (int, float)):
        return True
    # Check if both are arrays or lists of same length
    if isinstance(x, (list, np.ndarray)) and isinstance(y, (list, np.ndarray)):
        return len(x) == len(y)
    return False