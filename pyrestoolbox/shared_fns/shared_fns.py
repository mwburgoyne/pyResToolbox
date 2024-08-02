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
from typing import Union, List

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
            print("Could not solve via bisection")
            sys.exit()
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
    if isinstance(input_data, np.ndarray):
        # Input is already a numpy array, just return it
        return input_data
    else:
        # Convert list, tuple, scalar, or other types to numpy array
        # Ensuring even scalars become arrays with one element
        return np.atleast_1d(input_data)
        
def process_input(input_data):
    # Check if input_data is a numpy array
    if isinstance(input_data, np.ndarray):
        if input_data.size == 1:
            # Return the single element if it's a single-element array
            return input_data.item()
        else:
            # Return the array itself if it's larger
            return input_data
    elif isinstance(input_data, list):
        if len(input_data) == 1:
            # If it's a single-element list, return the element
            return input_data[0]
        else:
            # If it's a larger list, convert it to a numpy array and return
            return np.array(input_data)
    else:
        # If it's a single value (not in a list or array), return it directly
        return input_data
        
def check_2_inputs(x: Union[float, List[float]], y: Union[float, List[float]]) -> bool:
    """ Check inputs that need to be matched, either both float or both lists of same length"""
    # Check if both are floats
    if isinstance(x, float) and isinstance(y, float):
        return True
    # Check if both are lists of floats
    elif isinstance(x, list) and isinstance(y, list):
        if all(isinstance(item, float) for item in x) and all(isinstance(item, float) for item in y):
            # Check if lists have the same length
            return len(x) == len(y)
    # If neither condition is met, return False
    return False