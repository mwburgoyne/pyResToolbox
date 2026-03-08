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

Sensitivity analysis framework for parameter sweeps and tornado charts.

Functions
---------
sweep       Vary one parameter across a range, collecting results
tornado     Compute tornado-chart sensitivities for multiple parameters

Classes
-------
SweepResult     Result from sweep()
TornadoEntry    Single parameter sensitivity result
TornadoResult   Result from tornado()
"""

__all__ = [
    'sweep', 'tornado',
    'SweepResult', 'TornadoEntry', 'TornadoResult',
]

from dataclasses import dataclass, field
from typing import Callable, Dict, Any, List, Optional


@dataclass
class SweepResult:
    """Result from a parameter sweep.

    Attributes
    ----------
    param : str
        Name of the varied parameter.
    values : list
        Parameter values used.
    results : list
        Function results at each parameter value.
    """
    param: str
    values: list
    results: list


@dataclass
class TornadoEntry:
    """Single parameter sensitivity for tornado chart.

    Attributes
    ----------
    param : str
        Parameter name.
    low_value : float
        Low parameter value used.
    high_value : float
        High parameter value used.
    low_result : float
        Result at low parameter value.
    high_result : float
        Result at high parameter value.
    sensitivity : float
        Absolute relative sensitivity |high_result - low_result| / |base_result|.
    """
    param: str
    low_value: float
    high_value: float
    low_result: float
    high_result: float
    sensitivity: float


@dataclass
class TornadoResult:
    """Result from tornado sensitivity analysis.

    Attributes
    ----------
    base_result : float
        Result at base parameter values.
    entries : list of TornadoEntry
        Per-parameter sensitivities, sorted by decreasing sensitivity.
    """
    base_result: float
    entries: List[TornadoEntry] = field(default_factory=list)


def _extract_result(result, result_key):
    """Extract a scalar from a function result using result_key."""
    if result_key is None:
        return result
    if isinstance(result, dict):
        return result[result_key]
    return getattr(result, result_key)


def sweep(func: Callable, base_kwargs: Dict[str, Any],
          vary_param: str, vary_values: list,
          result_key: Optional[str] = None) -> SweepResult:
    """Vary one parameter across a range, collecting results.

    Parameters
    ----------
    func : callable
        Function to evaluate.
    base_kwargs : dict
        Base keyword arguments for the function.
    vary_param : str
        Name of the parameter to vary.
    vary_values : list
        Values to assign to vary_param.
    result_key : str, optional
        Key or attribute to extract from each result. If None, the raw
        result is stored.

    Returns
    -------
    SweepResult
    """
    if vary_param not in base_kwargs:
        raise ValueError(f"vary_param '{vary_param}' not found in base_kwargs")

    results = []
    for val in vary_values:
        kwargs = dict(base_kwargs)
        kwargs[vary_param] = val
        raw = func(**kwargs)
        results.append(_extract_result(raw, result_key))

    return SweepResult(param=vary_param, values=list(vary_values), results=results)


def tornado(func: Callable, base_kwargs: Dict[str, Any],
            ranges: Dict[str, tuple],
            result_key: Optional[str] = None) -> TornadoResult:
    """Compute tornado-chart sensitivities for multiple parameters.

    Parameters
    ----------
    func : callable
        Function to evaluate.
    base_kwargs : dict
        Base keyword arguments for the function.
    ranges : dict
        Mapping of parameter names to (low, high) tuples.
    result_key : str, optional
        Key or attribute to extract a scalar from each result.

    Returns
    -------
    TornadoResult
        Contains base_result and entries sorted by decreasing sensitivity.
    """
    base_raw = func(**base_kwargs)
    base_val = _extract_result(base_raw, result_key)
    base_result = float(base_val)

    entries = []
    for param, (lo, hi) in ranges.items():
        if param not in base_kwargs:
            raise ValueError(f"Parameter '{param}' not found in base_kwargs")

        kwargs_lo = dict(base_kwargs)
        kwargs_lo[param] = lo
        lo_result = float(_extract_result(func(**kwargs_lo), result_key))

        kwargs_hi = dict(base_kwargs)
        kwargs_hi[param] = hi
        hi_result = float(_extract_result(func(**kwargs_hi), result_key))

        if abs(base_result) > 0:
            sensitivity = abs(hi_result - lo_result) / abs(base_result)
        else:
            sensitivity = abs(hi_result - lo_result)

        entries.append(TornadoEntry(
            param=param,
            low_value=lo, high_value=hi,
            low_result=lo_result, high_result=hi_result,
            sensitivity=sensitivity,
        ))

    entries.sort(key=lambda e: e.sensitivity, reverse=True)
    return TornadoResult(base_result=base_result, entries=entries)
