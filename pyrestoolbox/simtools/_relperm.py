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

Relative permeability utilities: Corey/LET/Jerauld models, ECLIPSE table
generation (SWOF/SGOF/SGWFN) and model fitting.
"""

from typing import Union

import numpy as np
import numpy.typing as npt
import pandas as pd
from tabulate import tabulate

from pyrestoolbox.classes import kr_family, kr_table


def LET(s: np.ndarray, L: float, E: float, T: float) -> np.ndarray:
    """
    Returns LET Relative Permeability curve - Lomeland, F.; Ebeltoft, E.; Thomas, W.H. (2005).

    Input:
    s: Normalized saturation of the phase (np.array)
    L, E, T: Three correlation parameters 'L', 'E', and 'T'.

    Output:
    Relative permeability of the phase (np.array)
    """
    denom = s ** L + E * (1 - s) ** T
    denom = np.where(np.abs(denom) < 1e-30, 1e-30, denom)
    return np.clip(s ** L / denom, 0.0, 1.0)


def corey(s: np.ndarray, n: float) -> np.ndarray:
    """
    Returns Corey Relative Permeability curve.

    Input:
    s: Normalized saturation of the phase (np.array)
    n: Corey curve exponent.

    Output:
    Relative permeability of the phase (np.array)
    """
    return s ** n


def jerauld(s: np.ndarray, a: float, b: float) -> np.ndarray:
    """
    Returns Jerauld (Arco) Relative Permeability curve.

    Two-parameter model: kr = ((1+b) * S^a) / (1 + b * S^c)
    where c = a * (1 + 1/b).

    Input:
    s: Normalized saturation of the phase (np.array)
    a, b: Jerauld correlation parameters (both must be > 0).

    Output:
    Relative permeability of the phase (np.array)
    """
    b = max(b, 1e-10)
    c = a * (1 + 1 / b)
    denom = 1 + b * s ** c
    denom = np.where(np.abs(denom) < 1e-30, 1e-30, denom)
    return np.clip(((1 + b) * s ** a) / denom, 0.0, 1.0)


def is_let_physical(s: np.ndarray, L: float, E: float, T: float) -> bool:
    """
    Checks if a LET curve is physical (monotonically increasing with no
    spurious inflection-point reversals in concavity).

    Input:
    s: Normalized saturation array (should be sorted, e.g. np.linspace(0, 1, 100))
    L, E, T: LET correlation parameters.

    Output:
    True if the curve passes physicality checks, False otherwise.
    """
    s = np.sort(s)
    kr = LET(s, L, E, T)
    dkr = np.diff(kr)
    if np.any(dkr < -1e-10):
        return False
    ddkr = np.diff(dkr)
    for i in range(1, len(ddkr)):
        if np.sign(ddkr[i]) - np.sign(ddkr[i - 1]) > 0:
            return False
    return True


def _normalize_saturation(s_actual, s_min, s_max):
    """Normalize actual saturation to [0, 1] range."""
    span = s_max - s_min
    if span <= 0:
        return np.zeros_like(s_actual)
    return np.clip((s_actual - s_min) / span, 0, 1)


def _apply_kr_model(sn, krfamily, n_exp, L, E, T, a_jer, b_jer):
    """Apply the selected kr model to normalized saturation.

    sn: Normalized saturation array [0, 1]
    krfamily: kr_family enum instance
    n_exp: Corey exponent
    L, E, T: LET parameters
    a_jer, b_jer: Jerauld parameters
    Returns: kr array
    """
    fname = krfamily.name
    if fname == "COR":
        return corey(sn, n_exp)
    elif fname == "LET":
        return LET(sn, L, E, T)
    elif fname == "JER":
        return jerauld(sn, a_jer, b_jer)
    else:
        raise ValueError(f"Unknown kr family: {fname}")

def _build_kr_table(rows, krfamily, grid_lo, grid_hi, anchors, sat_col,
                    curves, export, keyword):
    """Shared builder for SWOF / SGOF / SGWFN tables.

    rows: Total number of table rows to return
    krfamily: kr_family enum instance
    grid_lo, grid_hi: Saturation range spanned by the regular grid
    anchors: Saturation values that must appear as rows (grid_lo must be one)
    sat_col: Name of the saturation column
    curves: List of (col_name, norm_lo, norm_hi, decreasing, krmax, params)
            where params = (n_exp, L, E, T, a_jer, b_jer)
    export: If True, write an ECLIPSE include file named <keyword>.INC
    keyword: ECLIPSE keyword (SWOF / SGOF / SGWFN)
    """
    unique_anchors = sorted(set(anchors))
    # The grid contributes ndiv points: grid_lo included, grid_hi excluded.
    # Every anchor other than grid_lo adds one extra row, so size the grid
    # to honour the requested total row count exactly.
    ndiv = max(rows - (len(unique_anchors) - 1), 1)
    sn = np.linspace(0.0, 1.0, ndiv, endpoint=False)
    s = sn * (grid_hi - grid_lo) + grid_lo
    s = np.array(sorted(set(list(s) + unique_anchors)))

    kr_df = pd.DataFrame()
    kr_df[sat_col] = s
    for col_name, norm_lo, norm_hi, decreasing, krmax, params in curves:
        n_exp, L, E, T, a_jer, b_jer = params
        s_norm = _normalize_saturation(s, norm_lo, norm_hi)
        if decreasing:
            s_norm = 1 - s_norm
        kr_df[col_name] = krmax * _apply_kr_model(
            s_norm, krfamily, n_exp, L, E, T, a_jer, b_jer)

    if export:
        df = kr_df.set_index(sat_col)
        headings = ["-- " + sat_col] + [c[0] for c in curves]
        fileout = keyword + "\n" + tabulate(df, headings) + "\n/"
        with open(keyword + ".INC", "w") as text_file:
            text_file.write(fileout)
    return kr_df


def rel_perm_table(
    rows: int,
    krtable: kr_table = kr_table.SWOF,
    krfamily: kr_family = kr_family.COR,
    kromax: float = 1,
    krgmax: float = 1,
    krwmax: float = 1,
    swc: float = 0,
    swcr: float = 0,
    sorg: float = 0,
    sorw: float = 0,
    sgcr: float = 0,
    no: float = 1,
    nw: float = 1,
    ng: float = 1,
    Lw: float = 1,
    Ew: float = 1,
    Tw: float = 1,
    Lo: float = 1,
    Eo: float = 1,
    To: float = 1,
    Lg: float = 1,
    Eg: float = 1,
    Tg: float = 1,
    aw: float = 1,
    bw: float = 1,
    ao: float = 1,
    bo: float = 1,
    ag: float = 1,
    bg: float = 1,
    export: bool = False,
) -> pd.DataFrame:
    """ Returns ECLIPSE styled relative permeability tables
        Users need only define parameters relevant to their table / family selection
        rows: Integer value specifying the number of table rows desired
        krtable: A string or kr_table Enum class that specifies one of three table type choices;
                   SWOF: Water / Oil table
                   SGOF: Gas / Oil table
                   SGWFN: Gas / Water table
        krfamily: A string or kr_family Enum class that specifies one of three curve function choices;
                   COR: Corey Curve function
                   LET: LET Relative permeability function
                   JER: Jerauld (Arco) two-parameter model
        kromax: Max Kr relative to oil. Default value = 1
        krgmax: Max Kr relative to gas. Default value = 1
        krwmax: Max Kr relative to water. Default value = 1
        swc: Minimum water saturation. Default value = 0
        swcr: Maximum water saturation for imobile water. Default value = 0
        sorg: Maximum oil saturation relative to gas for imobile oil. Default value = 0
        sorw: Maximum oil saturation relative to water for imobile oil. Default value = 0
        sgcr: Maximum gas saturation for immobile gas (honoured in both SGOF and SGWFN tables). Default value = 0
        no, nw, ng: Corey exponents to oil, water and gas respectively. Default values = 1
        Lw, Ew, Tw: LET exponents to water. Default values = 1
        Lo, Eo, To: LET exponents to oil. Default values = 1
        Lg, Eg, Tg: LET exponents to gas. Default values = 1
        aw, bw: Jerauld parameters for water. Default values = 1
        ao, bo: Jerauld parameters for oil. Default values = 1
        ag, bg: Jerauld parameters for gas. Default values = 1
        export: Boolean value that controls whether an include file with same name as krtable is created. Default: False
    """

    if isinstance(krtable, str):
        try:
            krtable = kr_table[krtable.upper()]
        except KeyError:
            raise ValueError(f"Incorrect table type specified: '{krtable}'")
    if isinstance(krfamily, str):
        try:
            krfamily = kr_family[krfamily.upper()]
        except KeyError:
            raise ValueError(f"Incorrect krfamily specified: '{krfamily}'")

    if rows < 2:
        raise ValueError(f"rows must be at least 2, got {rows}")

    # Consistency checks
    errors = []
    swcr = max(swc, swcr)
    if sorg + sgcr + swc >= 1:
        errors.append("sorg+sgcr+swc must be less than 1")
    if sorw + swcr >= 1:
        errors.append("sorw+swcr must be less than 1")
    if errors:
        raise ValueError("Saturation consistency check failure: " + "; ".join(errors))

    fname = krfamily.name

    # Validate that required parameters are available for the chosen family
    if fname == "COR":
        if krtable.name == "SWOF" and not (no > 0 and nw > 0):
            raise ValueError("Corey SWOF requires no > 0 and nw > 0")
        if krtable.name == "SGOF" and not (no > 0 and ng > 0):
            raise ValueError("Corey SGOF requires no > 0 and ng > 0")
        if krtable.name == "SGWFN" and not (ng > 0 and nw > 0):
            raise ValueError("Corey SGWFN requires ng > 0 and nw > 0")
    elif fname == "LET":
        if krtable.name == "SWOF" and not (Lw > 0 and Ew > 0 and Tw > 0 and Lo > 0 and Eo > 0 and To > 0):
            raise ValueError("LET SWOF requires Lw, Ew, Tw, Lo, Eo, To all > 0")
        if krtable.name == "SGOF" and not (Lg > 0 and Eg > 0 and Tg > 0 and Lo > 0 and Eo > 0 and To > 0):
            raise ValueError("LET SGOF requires Lg, Eg, Tg, Lo, Eo, To all > 0")
        if krtable.name == "SGWFN" and not (Lw > 0 and Ew > 0 and Tw > 0 and Lg > 0 and Eg > 0 and Tg > 0):
            raise ValueError("LET SGWFN requires Lw, Ew, Tw, Lg, Eg, Tg all > 0")
    elif fname == "JER":
        if krtable.name == "SWOF" and not (ao > 0 and bo > 0 and aw > 0 and bw > 0):
            raise ValueError("Jerauld SWOF requires ao, bo, aw, bw all > 0")
        if krtable.name == "SGOF" and not (ao > 0 and bo > 0 and ag > 0 and bg > 0):
            raise ValueError("Jerauld SGOF requires ao, bo, ag, bg all > 0")
        if krtable.name == "SGWFN" and not (ag > 0 and bg > 0 and aw > 0 and bw > 0):
            raise ValueError("Jerauld SGWFN requires ag, bg, aw, bw all > 0")

    water_params = (nw, Lw, Ew, Tw, aw, bw)
    oil_params = (no, Lo, Eo, To, ao, bo)
    gas_params = (ng, Lg, Eg, Tg, ag, bg)

    if krtable.name == "SWOF":
        # Grid spans [swcr, 1-sorw]; water mobile above swcr, oil immobile above 1-sorw
        return _build_kr_table(
            rows, krfamily, grid_lo=swcr, grid_hi=1 - sorw,
            anchors=[swc, swcr, 1 - sorw, 1], sat_col="Sw",
            curves=[
                ("Krwo", swcr, 1 - sorw, False, krwmax, water_params),
                ("Krow", swc, 1 - sorw, True, kromax, oil_params),
            ],
            export=export, keyword="SWOF")

    elif krtable.name == "SGOF":
        # Grid spans [0, 1-swc-sorg]; gas mobile above sgcr (anchor row at sgcr, krg=0)
        return _build_kr_table(
            rows, krfamily, grid_lo=0, grid_hi=1 - swc - sorg,
            anchors=[0, sgcr, 1 - swc - sorg], sat_col="Sg",
            curves=[
                ("Krgo", sgcr, 1 - swc - sorg, False, krgmax, gas_params),
                ("Krog", 0, 1 - swc - sorg, True, kromax, oil_params),
            ],
            export=export, keyword="SGOF")

    elif krtable.name == "SGWFN":
        # Grid spans [0, 1-swc]; gas mobile above sgcr (anchor row at sgcr, krg=0)
        return _build_kr_table(
            rows, krfamily, grid_lo=0, grid_hi=1 - swc,
            anchors=[0, sgcr, 1 - swc], sat_col="Sg",
            curves=[
                ("Krgw", sgcr, 1 - swc, False, krgmax, gas_params),
                ("Krwg", 0, 1 - swc, True, krwmax, water_params),
            ],
            export=export, keyword="SGWFN")

    raise ValueError("Check that you have specified table type as SWOF, SGOF or SGWFN")


def fit_rel_perm(
    sw: Union[npt.ArrayLike, list],
    kr: Union[npt.ArrayLike, list],
    krfamily: kr_family = kr_family.COR,
    krmax: float = 1.0,
    sw_min: float = 0.0,
    sw_max: float = 1.0,
) -> dict:
    """ Fits a relative permeability model (Corey, LET, or Jerauld) to measured kr data.

        Uses scipy.optimize.least_squares to find the best-fit parameters by minimizing
        the sum of squared residuals between the model and measured data.

        sw: Array of water saturation (or gas saturation) values corresponding to kr data.
        kr: Array of measured relative permeability values.
        krfamily: A string or kr_family Enum class specifying the model to fit;
                  COR: Corey (1 parameter: n)
                  LET: LET (3 parameters: L, E, T)
                  JER: Jerauld (2 parameters: a, b)
        krmax: Maximum kr endpoint value. The model is scaled by this value. Default = 1.0
        sw_min: Minimum saturation endpoint (critical saturation). Used to normalize sw to [0,1]. Default = 0.0
        sw_max: Maximum saturation endpoint. Used to normalize sw to [0,1]. Default = 1.0

        Returns: Dictionary with keys:
            'params': Dict of fitted parameter names and values (e.g. {'n': 2.5} for Corey)
            'krmax': The krmax used
            'sw_min': The sw_min used
            'sw_max': The sw_max used
            'residuals': Array of residuals at the solution
            'ssq': Sum of squared residuals
            'kr_pred': Array of predicted kr values at the input sw points
    """
    from scipy.optimize import least_squares

    sw = np.asarray(sw, dtype=float)
    kr = np.asarray(kr, dtype=float)

    if len(sw) != len(kr):
        raise ValueError(f"sw and kr must have same length (got {len(sw)} and {len(kr)})")
    if len(sw) < 2:
        raise ValueError("Need at least 2 data points for fitting")

    if isinstance(krfamily, str):
        try:
            krfamily = kr_family[krfamily.upper()]
        except KeyError:
            raise ValueError(f"Unknown kr family: '{krfamily}'")

    # Normalize saturation
    sn = _normalize_saturation(sw, sw_min, sw_max)

    # Scale measured kr by krmax
    kr_scaled = kr / max(krmax, 1e-30)

    fname = krfamily.name

    if fname == "COR":
        def residuals(p):
            n = p[0]
            return corey(sn, n) - kr_scaled
        x0 = [2.0]
        bounds = ([0.1], [20.0])
    elif fname == "LET":
        def residuals(p):
            L, E, T = p
            return LET(sn, L, E, T) - kr_scaled
        x0 = [2.0, 1.0, 2.0]
        bounds = ([0.01, 0.01, 0.01], [50.0, 50.0, 50.0])
    elif fname == "JER":
        def residuals(p):
            a, b = p
            return jerauld(sn, a, b) - kr_scaled
        x0 = [2.0, 1.0]
        bounds = ([0.1, 0.01], [20.0, 100.0])
    else:
        raise ValueError(f"Unknown kr family: {fname}")

    result = least_squares(residuals, x0, bounds=bounds)

    # Build prediction at input points
    kr_pred = krmax * (result.fun + kr_scaled)

    # Package results
    if fname == "COR":
        params = {"n": result.x[0]}
    elif fname == "LET":
        params = {"L": result.x[0], "E": result.x[1], "T": result.x[2]}
    elif fname == "JER":
        params = {"a": result.x[0], "b": result.x[1]}

    return {
        "params": params,
        "krmax": krmax,
        "sw_min": sw_min,
        "sw_max": sw_max,
        "residuals": result.fun,
        "ssq": float(np.sum(result.fun ** 2)),
        "kr_pred": kr_pred,
    }


def fit_rel_perm_best(
    sw: Union[npt.ArrayLike, list],
    kr: Union[npt.ArrayLike, list],
    krmax: float = 1.0,
    sw_min: float = 0.0,
    sw_max: float = 1.0,
) -> dict:
    """ Fits all three relative permeability models (Corey, LET, Jerauld) to measured data
        and returns the best fit (lowest sum of squared residuals).

        sw: Array of water saturation (or gas saturation) values corresponding to kr data.
        kr: Array of measured relative permeability values.
        krmax: Maximum kr endpoint value. Default = 1.0
        sw_min: Minimum saturation endpoint (critical saturation). Default = 0.0
        sw_max: Maximum saturation endpoint. Default = 1.0

        Returns: Dictionary with same keys as fit_rel_perm, plus:
            'family': Name of the best-fit model ('COR', 'LET', or 'JER')
            'all_results': Dict of all three fit results keyed by family name
    """
    all_results = {}
    for fam in [kr_family.COR, kr_family.LET, kr_family.JER]:
        all_results[fam.name] = fit_rel_perm(sw, kr, krfamily=fam,
                                              krmax=krmax, sw_min=sw_min, sw_max=sw_max)

    best_name = min(all_results, key=lambda k: all_results[k]["ssq"])
    best = all_results[best_name].copy()
    best["family"] = best_name
    best["all_results"] = all_results
    return best
