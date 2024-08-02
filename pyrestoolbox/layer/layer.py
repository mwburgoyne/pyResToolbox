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

import sys
from collections import Counter
import glob
from enum import Enum
import pkg_resources

import numpy as np
import numpy.typing as npt
import pandas as pd
from tabulate import tabulate

from pyrestoolbox.shared_fns import bisect_solve

def lorenz2b(lorenz: float, lrnz_method: str = "EXP") -> float:
    """ Returns B-factor that characterizes the Lorenz function
        Lorenz: Lorenz coefficient (0-1)
        lrnz_method: The method of calculation for the Lorenz coefficient
                Must be 'EXP' (Exponential) or 'LANG' (Langmuir).
                Defaults to EXP if undefined
                Background on Exponential formulation can be found in https://www.linkedin.com/pulse/loving-lorenz-new-life-old-parameter-mark-burgoyne/
                For Langmuir formulation; SumKh = Phih * VL / (Phih + PL)
                Lorenz = (VL - PL * VL * np.log(VL) + PL * VL * np.log(PL) - 0.5) * 2
                Where PL = 1 / B and VL = PL + 1
    """
    method = lrnz_method.upper()
    if method != "EXP" and method != "LANG":
        print('Method must be "LANG" or "EXP"')
        sys.exit()

    if lorenz < 0.000333:
        B = 2 / 1000
        if method == "LANG":
            B = 1 / 1000
        return B
    if lorenz > 0.997179125528914:
        B = 709
        if method == "LANG":
            B = 25000
        return B

    # Set bookends for B
    hi = 709
    if method == "LANG":
        hi = 25000
    lo = 0.000001
    args = (lorenz, method)

    def LorenzErr(args, B):
        lorenz, method = args
        B = max(B, 0.000001)
        if method == "EXP":
            B = min(B, 709)
            err = 2 * ((1 / (np.exp(B) - 1)) - (1 / B)) + 1 - lorenz
        else:
            B = min(B, 25000)
            PL = 1 / B
            VL = PL + 1
            err = (
                VL - PL * VL * np.log(VL) + PL * VL * np.log(PL) - 0.5
            ) * 2 - lorenz
        return err

    rtol = 0.0000001
    return bisect_solve(args, LorenzErr, lo, hi, rtol)


def lorenzfromb(B: float, lrnz_method: str = "EXP") -> float:
    """ Returns Lorenz coefficient that corresponds to a Beta value
        B: The B-Factor (positive float)
        lrnz_method: The method of calculation for the Lorenz coefficient
                Must be 'EXP' or 'LANG'.
                Defaults to Exponential if undefined
                Background on Exponential formulation can be found in https://www.linkedin.com/pulse/loving-lorenz-new-life-old-parameter-mark-burgoyne/
                For Langmuir formulation; SumKh = Phih * VL / (Phih + PL)
                Lorenz = (VL - PL * VL * np.log(VL) + PL * VL * np.log(PL) - 0.5) * 2
                Where PL = 1 / B and VL = PL + 1
    """
    method = lrnz_method.upper()
    B = max(B, 0.000001)
    if method == "LANG":
        B = min(B, 25000)
        PL = 1 / B
        VL = PL + 1
        L = (VL - PL * VL * np.log(VL) + PL * VL * np.log(PL) - 0.5) * 2
    else:
        B = min(B, 709)
        L = 2 * (1 / (np.exp(B) - 1) - (1 / B)) + 1
    return L


def lorenz_from_flow_fraction(
    kh_frac: float, phih_frac: float, lrnz_method: str = "EXP"
) -> float:
    """ Returns Lorenz coefficient consistent with observed flow fraction from a phi_h fraction
        kh_frac: (0 - 1). Fraction of total flow from best quality reservoir phi_h
        phih_frac: (0 - 1). phi_h fraction that delivers the observed kh_fraction of flow
        lrnz_method: The method of calculation for the Lorenz coefficient
                Must be 'EXP' or 'LANG'.
                Defaults to Exponential if undefined
                Background on Exponential formulation can be found in https://www.linkedin.com/pulse/loving-lorenz-new-life-old-parameter-mark-burgoyne/
                For Langmuir formulation; SumKh = Phih * VL / (Phih + PL)
                Lorenz = (VL - PL * VL * np.log(VL) + PL * VL * np.log(PL) - 0.5) * 2
                Where PL = 1 / B and VL = PL + 1
    """
    method = lrnz_method.upper()
    if kh_frac <= phih_frac:  #
        print("kh fraction should always be greater than phi_h fraction")
        return 0.001
    if kh_frac >= 1:
        print("kh Fraction must be less than 1")
        return 0.001

    # If Langmuir method, can explicitly calculate B
    if method == "LANG":
        x = phih_frac
        y = kh_frac
        B = (y - x) / (x * (1 - y))
        return lorenzfromb(B, method)

    # Set bookends and first guess of B
    hi = 709
    lo = 0.000001
    args = (kh_frac, phih_frac, method)

    def BErr(args, B):
        kh_frac, phih_frac, method = args
        method = method.upper()
        B = max(B, 0.000001)
        if method == "EXP":
            B = min(B, 709)
            err = (1 - np.exp(-B * phih_frac)) / (1 - np.exp(-B)) - kh_frac
        else:
            B = min(B, 25000)
            PL = 1 / B
            VL = PL + 1
            err = (VL * phih_frac) / (PL + phih_frac) - kh_frac
        return err

    rtol = 0.0000001
    B = bisect_solve(args, BErr, lo, hi, rtol)
    return lorenzfromb(B, method)


def lorenz_2_flow_frac(
    lorenz: float, phih_frac: float, lrnz_method: str = "EXP", B: float = -1
) -> float:
    """ Returns expected flow fraction from the best phi_h fraction, with a specified Lorenz coefficient

        lorenz: (0-1) Lorenz hetrogeneity factor
        phih_frac: (0 - 1). Best phi_h fraction
        lrnz_method: The method of calculation for the Lorenz coefficient
                Must be 'EXP' or 'LANG'.
                Defaults to Exponential if undefined
                Background on Exponential formulation can be found in https://www.linkedin.com/pulse/loving-lorenz-new-life-old-parameter-mark-burgoyne/
                For Langmuir formulation; SumKh = Phih * VL / (Phih + PL)
                Lorenz = (VL - PL * VL * np.log(VL) + PL * VL * np.log(PL) - 0.5) * 2
                Where PL = 1 / B and VL = PL + 1
        B: Factor that characterizes the Lorenz function for the given method. Will calculate if only lorenz variable defined
        lorenz: Lorenz coefficient (0-1). If B is provided, will ignore this parameter to be more efficient. If not, will calculate B from this parameter.
    """

    method = lrnz_method.upper()
    if B < 0 and lorenz < 0:
        print("Must define either B or lorenz parameters")
        sys.exit()

    if B < 0:  # Need to calculate B
        B = lorenz2b(lorenz=lorenz, lrnz_method=lrnz_method)

    B = max(B, 0.000001)
    if method == "EXP":
        B = min(B, 709)
        fraction = (1 - np.exp(-B * phih_frac)) / (1 - np.exp(-B))
    else:
        B = min(B, 25000)
        PL = 1 / B
        VL = PL + 1
        fraction = (VL * phih_frac) / (PL + phih_frac)
    return fraction


def lorenz_2_layers(
    lorenz: float,
    k_avg: float,
    nlayers: int = 1,
    shuffle: bool = False,
    lrnz_method: str = "EXP",
    B: float = -1,
    phi_h_fracs: list = [],
) -> np.ndarray:
    """ Returns np.array of permeability values honoring a specified average permeability (assuming equal thickness layers unless list of phi_h_fracs is provided), with degree of heterogeneity consistant with specified Lorenz coefficient and method

        If B is left default, then it will be calculated. If B is explictly specified > 0, then it will be used instead of the provided lorenz coefficient so as to eliminate repetitive solving for B.


        lorenz: Lorenz coefficient (0-1). If B is provided, will igonore this parameter to be more efficient. If not, will calculate B from this parameter.
        nlayers: The number of permeability layers desired (>1 needed unless a list of phi_h_fracs is supplied)
        kavg: The average permeability of all the layers (assuming equal thickness)
        shuffle: Boolean flag to determine whether to return the permeability array in decreasing order (False), or random order (True). Default False. Will be reset to False if user defined phi_h_fracs are supplied
        lrnz_method: The method of calculation for the Lorenz coefficient
                Must be 'EXP' or 'LANG'.
                Defaults to Exponential if undefined
                Background on Exponential formulation can be found in https://www.linkedin.com/pulse/loving-lorenz-new-life-old-parameter-mark-burgoyne/
                For Langmuir formulation; SumKh = Phih * VL / (Phih + PL)
                Lorenz = (VL - PL * VL * np.log(VL) + PL * VL * np.log(PL) - 0.5) * 2
                Where PL = 1 / B and VL = PL + 1
        B: Factor that characterizes the Lorenz function for the given method. Will calculate if only lorenz variable defined
        phi_h_fracs: Optional ability to specify a sorted list of phi_h fractions to get permeabilities for. If this list does not add to unity, then one additional layer permeability will be returned. The list needs to be in sorted order of best flow capacity to worst

    """
    if nlayers <= 1:
        if len(phi_h_fracs) < 2:
            return np.array([k_avg])

    method = lrnz_method.upper()

    if B < 0:  # Need to calculate B
        B = lorenz2b(lorenz=lorenz, lrnz_method=lrnz_method)

    B = max(B, 0.000001)
    if method == "EXP":
        B = min(B, 709)
    else:
        B = min(B, 25000)

    user_layers = False
    if len(phi_h_fracs) > 1:
        user_layers = True
        if sum(phi_h_fracs) > 1:
            phi_h_fracs = [x / sum(phi_h_fracs) for x in phi_h_fracs]
        if sum(phi_h_fracs) < 1:
            phi_h_fracs.append(1 - sum(phi_h_fracs))
        phih = (
            [0]
            + [sum(phi_h_fracs[: i + 1]) for i in range(len(phi_h_fracs) - 1)]
            + [1.0]
        )
        nlayers = len(phi_h_fracs)
    else:
        phih = np.arange(0, 1 + 1 / (nlayers), 1 / (nlayers))
        phi_h_fracs = np.array([1 / nlayers for i in range(len(phih) - 1)])
    sumkh = []

    for layer in phih:
        if method == "EXP":
            sumkh.append((1 - np.exp(-B * layer)) / (1 - np.exp(-B)))
        else:
            PL = 1 / B
            VL = PL + 1
            sumkh.append((VL * layer) / (PL + layer))

    kh = (
        np.array([sumkh[i] - sumkh[i - 1] for i in range(1, len(sumkh))])
        * k_avg
    )
    k = kh / np.array(phi_h_fracs)
    if shuffle:
        if not user_layers:
            np.random.shuffle(k)
    return k
