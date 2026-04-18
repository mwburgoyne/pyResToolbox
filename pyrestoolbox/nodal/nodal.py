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

Nodal analysis: VLP, IPR, and operating point calculations.

Classes
-------
WellSegment         Single wellbore segment with uniform geometry and deviation
Completion          Wellbore completion (legacy or multi-segment) for VLP calculations
Reservoir           Reservoir description for IPR calculations
NodalResult         Dict subclass returned by outflow_curve, ipr_curve, operating_point

Functions
---------
fbhp                Flowing bottom hole pressure via VLP correlation (HB, WG, GRAY, BB)
outflow_curve       VLP outflow curve (BHP vs rate)
ipr_curve           IPR inflow curve (rate vs Pwf)
operating_point     VLP/IPR intersection via bisection
"""

__all__ = [
    'WellSegment', 'Completion', 'Reservoir', 'NodalResult',
    'fbhp', 'fthp', 'outflow_curve', 'ipr_curve', 'operating_point',
]

import math
import warnings
from typing import Optional

import numpy as np

from pyrestoolbox.classes import vlp_method, class_dic
from pyrestoolbox.validate import validate_methods, validate_choice
from pyrestoolbox.shared_fns import bisect_solve, validate_pe_inputs
from pyrestoolbox.constants import (BAR_TO_PSI, PSI_TO_BAR, degc_to_degf, degf_to_degc,
                                    M_TO_FT, FT_TO_M, MM_TO_IN, IN_TO_MM,
                                    SM3_PER_SM3_TO_SCF_PER_STB, SCF_PER_STB_TO_SM3_PER_SM3,
                                    SM3_PER_SM3_TO_STB_PER_MSCF, STB_PER_MSCF_TO_SM3_PER_SM3,
                                    SM3_PER_SM3_TO_STB_PER_MMSCF, STB_PER_MMSCF_TO_SM3_PER_SM3,
                                    M3_TO_BBL, BBL_TO_M3, MSCF_TO_SM3, SM3_TO_MSCF,
                                    MMSCF_TO_SM3, SM3_TO_MMSCF,
                                    D_PER_SM3_TO_D_PER_MSCF, D_PER_MSCF_TO_D_PER_SM3,
                                    STB_TO_SM3, SM3_TO_STB)
import pyrestoolbox.gas as gas
import pyrestoolbox.oil as oil
from pyrestoolbox._accelerator import RUST_AVAILABLE as _RUST_AVAILABLE, rust_accelerated
if _RUST_AVAILABLE:
    from pyrestoolbox import _native as _rust


# Unit conversion constants used in VLP inner loops
_DYNCM_TO_LBFFT = 6.852e-5        # dyne/cm -> lbf/ft
_CP_TO_LBFTS = 6.7197e-4          # cP -> lb/(ft*s)
_LB_TO_KG = 0.453592              # lb -> kg
_LBFT3_TO_KGM3 = 16.01846         # lb/ft3 -> kg/m3
_R_GAS = 10.732                   # Gas constant, psia*ft3/(lb-mol*degR)


class NodalResult(dict):
    """Dict subclass returned by outflow_curve, ipr_curve, and operating_point.

    Supports both dict-style access (result['rate']) and attribute-style
    access (result.rate). Fully backward compatible with existing code
    that treats results as plain dicts.
    """
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(f"NodalResult has no attribute '{key}'")

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


# ============================================================================
#  Constants
# ============================================================================

_G_FT = 32.174          # ft/s^2
_GC = 32.174            # lbm*ft/(lbf*s^2)
_G_SI = 9.80665         # m/s^2
_P_ATM_PA = 101325.0    # Pa
_PSI_TO_PA = 6894.757
_FT_TO_M = 0.3048
_IN_TO_M = 0.0254
_LBFT3_TO_KGM3 = 16.01846
_CP_TO_PAS = 0.001
_DYNECM_TO_NM = 0.001   # dyne/cm -> N/m
_SCF_STB_TO_M3M3 = 0.178108
_LN10 = math.log(10.0)
_MW_AIR = 28.97


# ============================================================================
#  Standard Fluid & Unit Constants
# ============================================================================

_RHO_AIR_STC = 0.0765           # Air density at standard conditions (lb/scf)
_RHO_FW = 62.4                  # Freshwater density (lb/ft³)
_FT3_PER_BBL = 5.615            # ft³ per barrel
_SEC_PER_DAY = 86400.0          # seconds per day
_IN2_PER_FT2 = 144.0            # in² per ft²
_FW_GRAD = 0.433                # Freshwater pressure gradient (psi/ft)
_GAS_COL_K = 0.01875            # Gas column static pressure constant

# ============================================================================
#  IFT Correlation Constants
# ============================================================================

# Baker & Swerdloff (1956) dead oil IFT
_BS_TEMP_A = 1.11591            # Temperature coefficient (intercept)
_BS_TEMP_B = -0.00305           # Temperature coefficient (slope, per degC)
_BS_API_A = 38.085              # API coefficient (intercept)
_BS_API_B = -0.259              # API coefficient (slope, per API)

# Firoozabadi & Ramey (1988) gas-oil IFT
_FR_GOR_THRESH = 50.0           # GOR threshold (m³/m³) between low/high regimes
_FR_LOW_A = 0.02549             # Low-GOR regime coefficient
_FR_LOW_B = 1.0157              # Low-GOR regime exponent
_FR_HIGH_A = 32.0436            # High-GOR regime coefficient
_FR_HIGH_B = -1.1367            # High-GOR regime exponent

# Jennings & Newman (1971) gas-water IFT
_JN_74_A = 75.0                 # 74°F coefficient (intercept)
_JN_74_B = -1.108               # 74°F coefficient (pressure slope)
_JN_74_C = 0.349                # 74°F coefficient (pressure exponent)
_JN_280_A = 53.0                # 280°F coefficient (intercept)
_JN_280_B = -0.1048             # 280°F coefficient (pressure slope)
_JN_280_C = 0.637               # 280°F coefficient (pressure exponent)
_JN_T_LO = 74.0                 # Lower reference temperature (°F)
_JN_T_HI = 280.0                # Upper reference temperature (°F)

# Default IFT values (dyne/cm)
_IFT_OW_DEFAULT = 26.0          # Oil-water IFT default
_IFT_NOPHASE_DEFAULT = 20.0     # Default IFT when no phases flow

# ============================================================================
#  Serghides (1984) Friction Factor — Colebrook-White Constants
# ============================================================================

_CW_ROUGH = 3.7                 # Roughness denominator
_CW_RE = 12.0                   # Reynolds number term (laminar iteration)
_CW_TRANS = 2.51                # Transition term coefficient

# ============================================================================
#  HB (Hagedorn-Brown 1965) Holdup Constants
# ============================================================================

# Dimensionless group coefficients (also used in BB Payne correction)
_DG_VEL = 1.938                 # Velocity number coefficient
_DG_DIAM = 120.872              # Pipe diameter number coefficient
_DG_VISC = 0.15726              # Viscosity number coefficient

# CNL polynomial coefficients (5-term)
_CNL_C0 = -2.69851
_CNL_C1 = 0.1584095
_CNL_C2 = -0.5509976
_CNL_C3 = 0.5478492
_CNL_C4 = -0.1219458

# YL/NSI holdup polynomial coefficients (5-term)
_YL_C0 = -0.10306578
_YL_C1 = 0.617774
_YL_C2 = -0.632946
_YL_C3 = 0.29598
_YL_C4 = -0.0401

# S-correction polynomial coefficients (5-term)
_SC_C0 = 0.9116257
_SC_C1 = -4.821756
_SC_C2 = 1232.25
_SC_C3 = -22253.58
_SC_C4 = 116174.3

# S-correction thresholds
_SC_F2_CLAMP = 0.012            # f2 lower clamp
_SC_F2_FLOOR = 0.001            # f2 floor where S=1.0

# f1 exponents
_F1_VG_EXP = 0.575              # Gas velocity number exponent
_F1_P_ATM = 14.7                # Atmospheric pressure reference (psia)

# f2 exponents
_F2_VISC_EXP = 0.38             # Viscosity number exponent
_F2_DIAM_EXP = 2.14             # Diameter number exponent

# Orkiszewski (1967) bubble flow
_ORK_VS = 0.8                   # Bubble slip velocity (ft/s)
_ORK_LB_A = 1.071               # Bubble flow boundary coefficient
_ORK_LB_B = -0.2218             # Bubble flow boundary velocity coefficient
_ORK_LB_MIN = 0.13              # Minimum bubble flow boundary

# HB Reynolds/friction dimensional constants
_HB_RE_K = 96778.0              # Reynolds number dimensional constant
_HB_FRIC_K = 7.413e10           # Friction gradient dimensional constant

# ============================================================================
#  BB (Beggs & Brill 1973) Constants
# ============================================================================

# Flow pattern boundaries (L1-L4): coefficient, exponent
_BB_L1_A, _BB_L1_B = 316.0, 0.302
_BB_L2_A, _BB_L2_B = 0.0009252, -2.4684
_BB_L3_A, _BB_L3_B = 0.10, -1.4516
_BB_L4_A, _BB_L4_B = 0.5, -6.738

# Horizontal holdup (a, b, c per regime)
_BB_HL_SEG = (0.98, 0.4846, 0.0868)
_BB_HL_INT = (0.845, 0.5351, 0.0173)
_BB_HL_DIS = (1.065, 0.5824, 0.0609)

# Inclination correction (e, f, g, h per regime)
_BB_IC_SEG = (0.011, -3.7680, 3.5390, -1.6140)
_BB_IC_INT = (2.960, 0.3050, -0.4473, 0.0978)

# Payne et al. (1979) upward flow correction
_BB_PAYNE = 0.924

# Friction ratio S-factor polynomial
_BB_SF_C0 = -0.0523
_BB_SF_C1 = 3.182
_BB_SF_C2 = -0.8725
_BB_SF_C3 = 0.01853

# Friction S-factor clamp bounds
_BB_SF_LO = -5.0
_BB_SF_HI = 5.0

# ============================================================================
#  Gray (1978) Holdup & Roughness Constants
# ============================================================================

_GRAY_HL_A = 0.0814             # Holdup coefficient a
_GRAY_HL_B = 0.0554             # Holdup coefficient b (also used in f_e)
_GRAY_HL_C = 730.0              # Holdup coefficient c

_GRAY_ROUGH_K = 28.5            # Effective roughness coefficient

# ============================================================================
#  WG (Woldesemayat-Ghajar 2007) Drift-Flux Constants
# ============================================================================

_WG_DRIFT_K = 2.9               # Zuber-Findlay drift velocity coefficient
_WG_INCL_K = 1.22               # Inclination correction coefficient

# WG Chisholm C values (Lockhart-Martinelli)
_WG_CHIS_TT = 20.0              # Turbulent-turbulent
_WG_CHIS_LT = 12.0              # Laminar-turbulent
_WG_CHIS_TL = 10.0              # Turbulent-laminar
_WG_CHIS_LL = 5.0               # Laminar-laminar

# ============================================================================
#  Vogel / Darcy IPR Constants
# ============================================================================

_DARCY_K = 0.00708              # Darcy radial flow constant
_DIETZ_CORR = 0.75              # Dietz shape factor correction
_VOGEL_AOF_DENOM = 1.8          # Vogel AOF denominator
_VOGEL_LIN = 0.2                # Vogel linear coefficient
_VOGEL_QUAD = 0.8               # Vogel quadratic coefficient


def _log10(x):
    if isinstance(x, complex) or x <= 0:
        return -30.0
    return math.log(x) / _LN10


def _clamp(val, lo, hi):
    if isinstance(val, complex):
        val = abs(val)
    return max(lo, min(hi, val))


# ============================================================================
#  Data Classes: WellSegment, Completion and Reservoir
# ============================================================================

class WellSegment:
    """Single wellbore segment with uniform geometry and deviation.

        md: Measured depth of this segment (ft | m)
        id: Internal diameter (inches | mm)
        deviation: Deviation from vertical (degrees). 0=vertical, 90=horizontal. Defaults to 0
        roughness: Pipe roughness (inches | mm). Defaults to 0.0006 in (0.01524 mm)
        metric: If True, inputs in metric units (m, mm). Default False (ft, inches).
    """
    def __init__(self, md, id, deviation=0, roughness=None, metric=False):
        if roughness is None:
            roughness = 0.01524 if metric else 0.0006  # 0.0006 in = 0.01524 mm
        _md_in, _id_in, _rough_in = md, id, roughness
        _unit_len = 'm' if metric else 'ft'
        _unit_dia = 'mm' if metric else 'in'
        if metric:
            md = md * M_TO_FT
            id = id * MM_TO_IN
            roughness = roughness * MM_TO_IN
        if md <= 0:
            raise ValueError(f"Measured depth md must be positive, got {_md_in} {_unit_len}")
        if id <= 0:
            raise ValueError(f"Internal diameter id must be positive, got {_id_in} {_unit_dia}")
        if roughness < 0:
            raise ValueError(f"Roughness must be non-negative, got {_rough_in} {_unit_dia}")
        if not (0 <= deviation <= 90):
            raise ValueError(f"Deviation must be between 0 and 90 degrees, got {deviation}")
        self.md = md  # stored in ft
        self.id = id  # stored in inches
        self.deviation = deviation
        self.roughness = roughness  # stored in inches

    @property
    def tvd(self):
        """True vertical depth contribution of this segment (ft)."""
        return self.md * math.cos(math.radians(self.deviation))

    @property
    def theta(self):
        """Angle from horizontal (radians) for multiphase correlations."""
        return math.pi / 2.0 - math.radians(self.deviation)


class Completion:
    """ Wellbore completion description for VLP calculations.

        Can be constructed in two ways:

        Legacy mode (positional):
            tid: Tubing ID (inches | mm)
            length: Tubing length (ft | m) - from wellhead to tubing shoe
            tht: Tubing head (wellhead) temperature (degF | degC)
            bht: Bottom hole temperature (degF | degC)
            rough: Tubing roughness (inches | mm). Defaults to 0.0006 in (0.01524 mm)
            cid: Casing ID (inches | mm) below tubing shoe. Defaults to 0 (no casing section)
            crough: Casing roughness (inches | mm). Defaults to 0.0006 in
            mpd: Mid-perforation depth (ft | m). Defaults to length (no casing section)

        Segment mode (keyword):
            segments: List of WellSegment objects defining the wellbore
            tht: Tubing head (wellhead) temperature (degF | degC)
            bht: Bottom hole temperature (degF | degC)

        metric: If True, inputs in metric units (m, mm, degC). Default False.
                Note: WellSegments passed in segment mode should already be constructed
                with their own metric flag. The Completion metric flag only converts
                legacy mode dimensions and temperatures.
    """
    def __init__(self, tid=None, length=None, tht=None, bht=None, rough=None,
                 cid=0, crough=None, mpd=0, segments=None, metric=False):
        # Set roughness defaults appropriate for unit system
        if rough is None:
            rough = 0.01524 if metric else 0.0006  # 0.0006 in = 0.01524 mm
        if crough is None:
            crough = 0.01524 if metric else 0.0006

        self._metric = metric

        if segments is not None:
            # Segment mode — segments already store oilfield units internally
            if tht is None or bht is None:
                raise ValueError("tht and bht are required when using segments")
            if not segments:
                raise ValueError("segments list must not be empty")
            self._segments = list(segments)
            if metric:
                tht = degc_to_degf(tht)
                bht = degc_to_degf(bht)
            self.tht = tht  # stored in degF
            self.bht = bht  # stored in degF
            # Set legacy attributes from first segment for compatibility
            self.tid = self._segments[0].id
            self.rough = self._segments[0].roughness
            self.length = self._segments[0].md
            self.cid = 0
            self.crough = crough if not metric else crough * MM_TO_IN
            self.mpd = self.total_md
        elif tid is not None and length is not None:
            # Legacy mode
            if tht is None or bht is None:
                raise ValueError("tht and bht are required")
            if metric:
                tid = tid * MM_TO_IN
                length = length * M_TO_FT
                tht = degc_to_degf(tht)
                bht = degc_to_degf(bht)
                rough = rough * MM_TO_IN
                cid = cid * MM_TO_IN if cid > 0 else 0
                crough = crough * MM_TO_IN
                mpd = mpd * M_TO_FT if mpd > 0 else 0
            self.tid = tid  # stored in inches
            self.length = length  # stored in ft
            self.tht = tht  # stored in degF
            self.bht = bht  # stored in degF
            self.rough = rough  # stored in inches
            self.cid = cid  # stored in inches
            self.crough = crough  # stored in inches
            self.mpd = mpd if mpd > 0 else length  # stored in ft
            # Build segments from legacy parameters
            segs = [WellSegment(md=length, id=tid, deviation=0, roughness=rough)]
            if cid > 0 and self.mpd > length:
                segs.append(WellSegment(md=self.mpd - length, id=cid,
                                        deviation=0, roughness=crough))
            self._segments = segs
        else:
            raise ValueError("Must provide either (tid, length, tht, bht) or (segments, tht, bht)")

    @property
    def segments(self):
        """List of WellSegment objects defining the wellbore."""
        return self._segments

    @property
    def total_md(self):
        """Total measured depth across all segments (ft)."""
        return sum(s.md for s in self._segments)

    @property
    def total_tvd(self):
        """Total true vertical depth across all segments (ft)."""
        return sum(s.tvd for s in self._segments)

    @property
    def has_casing_section(self):
        return self.cid > 0 and self.mpd > self.length

    @property
    def casing_length(self):
        if self.has_casing_section:
            return self.mpd - self.length
        return 0.0

    @property
    def tubing_end_temperature(self):
        if self.has_casing_section:
            frac = self.length / self.mpd
            return self.tht + (self.bht - self.tht) * frac
        return self.bht

    def geometry_at_md(self, md):
        """Return wellbore geometry at a given measured depth.

        Args:
            md: Measured depth from surface (ft | m, matching Completion unit system)

        Returns:
            dict with keys: 'md', 'tvd', 'id', 'deviation', 'roughness'
            Values in oilfield units (ft, inches) or metric (m, mm) matching
            Completion construction.

        Raises:
            ValueError: if md < 0 or md > total_md
        """
        # Convert metric input to internal oilfield units
        md_ft = md * M_TO_FT if self._metric else md
        total = self.total_md
        if md_ft < 0 or md_ft > total + 1e-9:
            if self._metric:
                raise ValueError(
                    f"md must be between 0 and {total * FT_TO_M:.4f} m, got {md}"
                )
            raise ValueError(
                f"md must be between 0 and {total:.4f} ft, got {md}"
            )
        # Clamp to total_md for floating-point edge case
        md_ft = min(md_ft, total)

        md_traversed = 0.0
        tvd_traversed = 0.0
        for seg in self._segments:
            if md_traversed + seg.md >= md_ft - 1e-12:
                delta_md = md_ft - md_traversed
                tvd = tvd_traversed + delta_md * math.cos(math.radians(seg.deviation))
                if self._metric:
                    return {
                        'md': md_ft * FT_TO_M,
                        'tvd': tvd * FT_TO_M,
                        'id': seg.id * IN_TO_MM,
                        'deviation': seg.deviation,
                        'roughness': seg.roughness * IN_TO_MM,
                    }
                return {
                    'md': md_ft,
                    'tvd': tvd,
                    'id': seg.id,
                    'deviation': seg.deviation,
                    'roughness': seg.roughness,
                }
            md_traversed += seg.md
            tvd_traversed += seg.tvd
        # Should not reach here due to validation above
        raise RuntimeError("Failed to locate md within segments")  # pragma: no cover

    def profile(self):
        """Return wellbore profile at all segment boundaries.

        Returns:
            pandas DataFrame with columns: 'MD', 'TVD', 'Deviation', 'ID', 'Roughness'
            One row per node: surface (top of first segment), bottom of each
            segment, and top of each subsequent segment at crossover points.
            Units match Completion construction (oilfield or metric).
        """
        import pandas as pd

        mds, tvds, devs, ids, roughs = [], [], [], [], []
        md_cum = 0.0
        tvd_cum = 0.0

        for i, seg in enumerate(self._segments):
            if i == 0:
                # Top of first segment (surface)
                mds.append(md_cum)
                tvds.append(tvd_cum)
                devs.append(seg.deviation)
                ids.append(seg.id)
                roughs.append(seg.roughness)
            # Bottom of this segment
            md_cum += seg.md
            tvd_cum += seg.tvd
            mds.append(md_cum)
            tvds.append(tvd_cum)
            devs.append(seg.deviation)
            ids.append(seg.id)
            roughs.append(seg.roughness)
            # If there's a next segment, add its top at the same depth
            if i + 1 < len(self._segments):
                next_seg = self._segments[i + 1]
                mds.append(md_cum)
                tvds.append(tvd_cum)
                devs.append(next_seg.deviation)
                ids.append(next_seg.id)
                roughs.append(next_seg.roughness)

        if self._metric:
            mds = [m * FT_TO_M for m in mds]
            tvds = [t * FT_TO_M for t in tvds]
            ids = [d * IN_TO_MM for d in ids]
            roughs = [r * IN_TO_MM for r in roughs]

        return pd.DataFrame({
            'MD': mds,
            'TVD': tvds,
            'Deviation': devs,
            'ID': ids,
            'Roughness': roughs,
        })


class Reservoir:
    """ Reservoir description for IPR calculations.

        pr: Reservoir pressure (psia | barsa)
        degf: Reservoir temperature (degF | degC)
        k: Permeability (mD) — same in both unit systems
        h: Net pay thickness (ft | m)
        re: Drainage radius (ft | m)
        rw: Wellbore radius (ft | m)
        S: Skin factor. Defaults to 0
        D: Non-Darcy coefficient (day/Mscf | day/sm3 for gas, 0 for oil). Defaults to 0
        metric: If True, inputs in metric units (barsa, degC, m, day/sm3). Default False.
    """
    def __init__(self, pr, degf, k, h, re, rw, S=0, D=0, metric=False):
        # Preserve original user-facing values for error messages
        _pr_in, _degf_in, _h_in, _re_in, _rw_in = pr, degf, h, re, rw
        _unit_len = 'm' if metric else 'ft'
        if metric:
            pr = pr * BAR_TO_PSI
            degf = degc_to_degf(degf)
            h = h * M_TO_FT
            re = re * M_TO_FT
            rw = rw * M_TO_FT
            if D != 0:
                D = D * D_PER_SM3_TO_D_PER_MSCF
        validate_pe_inputs(p=pr, degf=degf)
        if k <= 0:
            raise ValueError(f"Permeability k must be positive, got {k}")
        if h <= 0:
            raise ValueError(f"Net pay thickness h must be positive, got {_h_in} {_unit_len}")
        if rw <= 0:
            raise ValueError(f"Wellbore radius rw must be positive, got {_rw_in} {_unit_len}")
        if re <= rw:
            raise ValueError(
                f"Drainage radius re ({_re_in} {_unit_len}) must be greater than "
                f"wellbore radius rw ({_rw_in} {_unit_len})"
            )
        self.pr = pr  # stored in psia
        self.degf = degf  # stored in degF
        self.k = k  # mD (same in both systems)
        self.h = h  # stored in ft
        self.re = re  # stored in ft
        self.rw = rw  # stored in ft
        self.S = S
        self.D = D  # stored in day/Mscf


# ============================================================================
#  Simplified PVT helpers (standalone, no GasPVT/OilPVT dependency)
#  Used inside VLP segment loops for performance
# ============================================================================

def _sutton_tc_pc(sg):
    """Sutton critical properties from gas SG."""
    tpc = 169.2 + 349.5 * sg - 74.0 * sg * sg
    ppc = 756.8 - 131.0 * sg - 3.6 * sg * sg
    return tpc, ppc


def _z_factor(sg, temp_f, press_psia, tc=None, pc=None):
    """Z-factor via Hall-Yarborough (1973)."""
    if press_psia < 1.0:
        return 1.0
    if tc is None or pc is None:
        tc, pc = _sutton_tc_pc(sg)
    tr = (temp_f + 459.67) / tc
    pr = press_psia / pc

    t_inv = 1.0 / tr
    a = -0.06125 * t_inv * math.exp(-1.2 * (1.0 - t_inv) ** 2)
    b = 14.76 * t_inv - 9.76 * t_inv ** 2 + 4.58 * t_inv ** 3
    c = 90.7 * t_inv - 242.2 * t_inv ** 2 + 42.4 * t_inv ** 3
    d = 2.18 + 2.82 * t_inv

    y = _clamp(0.0125 * pr * t_inv, 1e-10, 0.9)

    for _ in range(50):
        y2, y3, y4 = y * y, y * y * y, y ** 4
        one_m_y = 1.0 - y
        if abs(one_m_y) < 1e-15:
            break
        fy = ((y + y2 + y3 - y4) / (one_m_y ** 3) + a * pr - b * y2 + c * y ** d)
        dfy = ((1.0 + 4.0 * y + 4.0 * y2 - 4.0 * y3 + y4) / (one_m_y ** 4) -
               2.0 * b * y + c * d * y ** (d - 1.0))
        if abs(dfy) < 1e-30:
            break
        dy = fy / dfy
        y -= dy
        y = _clamp(y, 1e-10, 0.99)
        if abs(dy) < 1e-12:
            break

    if abs(y) < 1e-30:
        return 1.0
    z = -a * pr / y
    return _clamp(z, 0.1, 5.0)


def _gas_viscosity(sg, temp_f, press_psia, z=None, tc=None, pc=None):
    """Gas viscosity via Lee-Gonzalez-Eakin (1966) (cP)."""
    if z is None:
        z = _z_factor(sg, temp_f, press_psia, tc=tc, pc=pc)
    temp_r = temp_f + 459.67
    mw = _MW_AIR * sg
    rho_gcc = (_MW_AIR * sg * press_psia / (z * _R_GAS * temp_r)) / 62.428

    k_val = ((9.4 + 0.02 * mw) * temp_r ** 1.5 / (209.0 + 19.0 * mw + temp_r))
    x_val = 3.5 + 986.0 / temp_r + 0.01 * mw
    y_val = 2.4 - 0.2 * x_val
    exp_arg = x_val * rho_gcc ** y_val
    if exp_arg > 500:
        exp_arg = 500  # Cap to avoid overflow
    return max(1e-4 * k_val * math.exp(exp_arg), 1e-6)


def _gas_density(sg, temp_f, press_psia, z=None, tc=None, pc=None):
    """Gas density in lb/ft^3."""
    if z is None:
        z = _z_factor(sg, temp_f, press_psia, tc=tc, pc=pc)
    return _MW_AIR * sg * press_psia / (z * _R_GAS * (temp_f + 459.67))


def _water_viscosity(press_psia, temp_f, salinity=0.0):
    """Simplified water viscosity (cP)."""
    if temp_f < 32:
        temp_f = 32.0
    a_coeff = -3.79418 + 604.129 / (139.18 + temp_f)
    mu_w = 10.0 ** a_coeff
    if salinity > 0:
        mu_w *= (1.0 + 0.02 * salinity)
    mu_w *= (1.0 + 5e-4 * (press_psia - 14.7) / 1000.0)
    return max(mu_w, 0.01)


# ============================================================================
#  Oil PVT helpers (standalone for VLP segment loops)
# ============================================================================

def _standing_rs(gsg, press_psia, temp_f, api):
    """Standing (1947) solution GOR estimate (scf/STB)."""
    exponent = 0.00091 * temp_f - 0.0125 * api
    denom = 18.2 * 10.0 ** exponent
    if denom <= 0:
        return 0.0
    return gsg * max(press_psia / denom, 0.0) ** 1.2048


def _velarde_rs(sgsp, api, temp_f, pb, rsb, press_psia):
    """Velarde (1997) solution GOR vs pressure (scf/STB)."""
    if press_psia >= pb:
        return rsb
    a1 = (0.000000973 * sgsp ** 1.672608 * api ** 0.92987 *
          temp_f ** 0.247235 * max(pb - 14.7, 0.1) ** 1.056052)
    a2 = (0.022339 * sgsp ** (-1.004750) * api ** 0.337711 *
          temp_f ** 0.132795 * max(pb - 14.7, 0.1) ** 0.302065)
    a3 = (0.725167 * sgsp ** (-1.485480) * api ** (-0.164741) *
          temp_f ** (-0.09133) * max(pb - 14.7, 0.1) ** 0.047094)
    prr = max((press_psia - 14.7) / max(pb - 14.7, 0.1), 0.0)
    return (a1 * prr ** a2 + (1.0 - a1) * prr ** a3) * rsb


def _oil_density_mccain(rs, sgsp, sgsto, press_psia, temp_f):
    """Reservoir oil density via McCain iterative method (lb/ft^3)."""
    rhopo1 = 52.8 - 0.01 * rs
    for _ in range(100):
        rhoa2 = (-49.893 + 85.0149 * sgsp - 3.70373 * sgsp * rhopo1 +
                 0.0479818 * sgsp * rhopo1 ** 2 + 2.98914 * rhopo1 -
                 0.035688 * rhopo1 ** 2)
        denom = 73.71 + rs * sgsp / rhoa2 if abs(rhoa2) > 1e-30 else 73.71
        rhopo2 = (rs * sgsp + 4600.0 * sgsto) / denom
        if abs(rhopo1 - rhopo2) < 0.0003:
            break
        rhopo1 = rhopo2

    rhopo = rhopo2
    drhop = ((0.167 + 16.181 * 10.0 ** (-0.0425 * rhopo)) *
             (press_psia / 1000.0) -
             0.01 * (0.299 + 263.0 * 10.0 ** (-0.0603 * rhopo)) *
             (press_psia / 1000.0) ** 2)
    rhobs = rhopo + drhop
    dt = max(temp_f - 60.0, 0.001)
    drhot = ((0.00302 + 1.505 * rhobs ** (-0.951)) * dt ** 0.938 -
             (0.0216 - 0.0233 * 10.0 ** (-0.0161 * rhobs)) * dt ** 0.475)
    return rhobs - drhot


def _oil_viscosity_full(sgsp, api, temp_f, rsb, pb, press_psia, vis_frac=1.0, rsb_frac=1.0):
    """Oil viscosity (Beggs-Robinson + undersaturated correction) (cP)."""
    rs = _velarde_rs(sgsp, api, temp_f, pb, rsb / rsb_frac, press_psia) * rsb_frac
    c = 10.0 ** (3.0324 - 0.02023 * api) * temp_f ** (-1.163)
    mu_od = 10.0 ** c - 1.0
    a = 10.715 * (rs + 100.0) ** (-0.515)
    b = 5.44 * (rs + 150.0) ** (-0.338)
    mu_or = a * mu_od ** b

    if press_psia > pb:
        ab = 10.715 * (rsb + 100.0) ** (-0.515)
        bb = 5.44 * (rsb + 150.0) ** (-0.338)
        mu_orb = ab * mu_od ** bb
        log_mu = _log10(mu_orb)
        aa = (-1.0146 + 1.3322 * log_mu - 0.4876 * log_mu ** 2 -
              1.15036 * log_mu ** 3)
        mu_or = mu_orb + 0.0013449 * (press_psia - pb) * 10.0 ** aa

    return max(mu_or, 0.001) * vis_frac


# ============================================================================
#  IFT Correlations
# ============================================================================

def _dead_oil_ift(api, temp_f):
    temp_c = (temp_f - 32.0) / 1.8
    a = _BS_TEMP_A + _BS_TEMP_B * temp_c
    return max(a * (_BS_API_A + _BS_API_B * api), 1.0)


def _gas_oil_ift(api, temp_f, rs_scf_stb):
    sigma_od = _dead_oil_ift(api, temp_f)
    if rs_scf_stb <= 0:
        return sigma_od
    rs_m3m3 = rs_scf_stb * _SCF_STB_TO_M3M3
    if rs_m3m3 < _FR_GOR_THRESH:
        ratio = 1.0 / (1.0 + _FR_LOW_A * rs_m3m3 ** _FR_LOW_B)
    else:
        ratio = _FR_HIGH_A * rs_m3m3 ** _FR_HIGH_B
    return max(sigma_od * _clamp(ratio, 0.0, 1.0), 1.0)


def _gas_water_ift(press_psia, temp_f):
    p = max(press_psia, 14.7)
    sw74 = _JN_74_A + _JN_74_B * p ** _JN_74_C
    sw280 = _JN_280_A + _JN_280_B * p ** _JN_280_C
    if temp_f <= _JN_T_LO:
        return max(sw74, 1.0)
    elif temp_f >= _JN_T_HI:
        return max(sw280, 1.0)
    return max(sw74 + (sw280 - sw74) * (temp_f - _JN_T_LO) / (_JN_T_HI - _JN_T_LO), 1.0)


def _interfacial_tension(p_avg, temp_f, api, rs_scf_stb,
                         m_flow_o, m_flow_w, m_flow_g):
    sigma_og = _gas_oil_ift(api, temp_f, rs_scf_stb)
    sigma_wg = _gas_water_ift(p_avg, temp_f)
    sigma_ow = _IFT_OW_DEFAULT
    denom = m_flow_o * m_flow_w + m_flow_o * m_flow_g + m_flow_g * m_flow_w
    if denom <= 0:
        return _IFT_NOPHASE_DEFAULT
    return max(((m_flow_o * m_flow_w) * sigma_ow +
                (m_flow_o * m_flow_g) * sigma_og +
                (m_flow_g * m_flow_w) * sigma_wg) / denom, 1.0)


# ============================================================================
#  Serghides Friction Factor
# ============================================================================

def _serghides_fanning(re, eps_d):
    if re < 1:
        return 0.0
    if re <= 2100:
        return 16.0 / re
    a = -2.0 * _log10(eps_d / _CW_ROUGH + _CW_RE / re)
    b = -2.0 * _log10(eps_d / _CW_ROUGH + _CW_TRANS * a / re)
    c = -2.0 * _log10(eps_d / _CW_ROUGH + _CW_TRANS * b / re)
    diff = c - 2.0 * b + a
    if abs(diff) < 1e-30:
        return 0.005
    f_darcy = (a - (b - a) ** 2 / diff) ** (-2)
    return f_darcy / 4.0


# ============================================================================
#  Segment and Static Column Helpers
# ============================================================================

def _calc_segments(length, min_seg_ft=100.0):
    ndiv = 100
    seg_len = length / ndiv
    if seg_len < min_seg_ft:
        ndiv = int(length / min_seg_ft) + 1
        if ndiv < 1:
            ndiv = 1
    return ndiv


def _condensate_dropout(cgr, qg_mmscfd, p_avg, pr, osg, qw_bwpd, wsg):
    if pr > 14.7:
        cgr_local = cgr * max(0.0, pr - p_avg) / (pr - 14.7)
    else:
        cgr_local = cgr
    qo_local = cgr_local * qg_mmscfd
    ql_local = qo_local + qw_bwpd
    if ql_local > 0:
        lsg_local = (qo_local * osg + qw_bwpd * wsg) / ql_local
    else:
        lsg_local = osg  # No liquid; value won't be used meaningfully
    return cgr_local, qo_local, ql_local, lsg_local


def _condensate_vis(pr, cgr_local, gsg, api, temp_f, p_avg, oil_vis):
    if pr > 14.7 and cgr_local > 0.01:
        return _oil_viscosity_full(gsg, api, temp_f, rsb=0.0, pb=14.7,
                                   press_psia=p_avg)
    return oil_vis


def _static_gas_column_pressure(thp, length, tht, bht, gsg, theta=math.pi / 2.0):
    tc, pc = _sutton_tc_pc(gsg)
    n_seg = 50
    d_len = length / n_seg
    sin_theta = math.sin(theta)
    p = thp
    for i in range(n_seg):
        frac = (i + 0.5) / n_seg
        temp_local = tht + (bht - tht) * frac
        temp_r = temp_local + 459.67
        zee = _z_factor(gsg, temp_local, max(p, 14.7), tc=tc, pc=pc)
        rho_gas = _MW_AIR * gsg * p / (zee * _R_GAS * temp_r)
        p += rho_gas / _IN2_PER_FT2 * d_len * sin_theta
    return p


def _static_oil_column_pressure(thp, length, tht, bht, wc, wsg,
                                 api, sgsp, pb, rsb, rsb_scale=1.0,
                                 theta=math.pi / 2.0):
    sgsto = 141.5 / (api + 131.5)
    rsb_for_calc = rsb / rsb_scale
    sin_theta = math.sin(theta)
    n_seg = 50
    d_len = length / n_seg
    p = thp
    for i in range(n_seg):
        frac = (i + 0.5) / n_seg
        temp_local = tht + (bht - tht) * frac
        if p < pb:
            rs = _velarde_rs(sgsp, api, temp_local, pb, rsb_for_calc, p) * rsb_scale
            rho_oil = _oil_density_mccain(rs, sgsp, sgsto, p, temp_local)
        else:
            rho_oil = _oil_density_mccain(rsb, sgsp, sgsto, pb, temp_local)
        oil_sg_local = rho_oil / _RHO_FW
        mix_sg = (1.0 - wc) * oil_sg_local + wc * wsg
        p += _FW_GRAD * mix_sg * d_len * sin_theta
    return p


# ============================================================================
#  1. HAGEDORN-BROWN VLP
# ============================================================================

def _hb_gradient_gas(s):
    """Hagedorn-Brown gradient: polynomial holdup, Orkiszewski bubble flow, Serghides friction."""
    ul = s['v_sl']
    ugas = s['v_sg']
    ift = s['sigma']
    rho_l = s['rho_l']

    nvl = _DG_VEL * ul * (rho_l / ift) ** 0.25
    nvg = _DG_VEL * ugas * (rho_l / ift) ** 0.25
    nd = _DG_DIAM * s['diam_ft'] * (rho_l / ift) ** 0.5

    mul = s['mu_l'] if s['ql_loc'] > 0 else 1.0

    nl = _DG_VISC * mul * (1.0 / (rho_l * ift ** 3)) ** 0.25
    if nl <= 0:
        nl = 1e-8

    x1 = _log10(nl) + 3.0
    cnl = 10.0 ** (_CNL_C0 + _CNL_C1 * x1 + _CNL_C2 * x1 * x1 +
                   _CNL_C3 * x1 ** 3 + _CNL_C4 * x1 ** 4)

    nvg_safe = max(nvg, 1e-10)
    f1 = nvl * s['p_avg'] ** 0.1 * cnl / (nvg_safe ** _F1_VG_EXP * _F1_P_ATM ** 0.1 * nd)

    lf1 = _log10(f1) + 6.0
    ylonsi = (_YL_C0 + _YL_C1 * lf1 + _YL_C2 * lf1 ** 2 +
              _YL_C3 * lf1 ** 3 + _YL_C4 * lf1 ** 4)
    ylonsi = max(ylonsi, 0.0)

    f2 = nvg * nl ** _F2_VISC_EXP / nd ** _F2_DIAM_EXP
    idex = 1.0 if f2 >= _SC_F2_CLAMP else -1.0
    f2 = (1.0 - idex) / 2.0 * _SC_F2_CLAMP + (1.0 + idex) / 2.0 * f2

    si = (_SC_C0 + _SC_C1 * f2 + _SC_C2 * f2 * f2 +
          _SC_C3 * f2 ** 3 + _SC_C4 * f2 ** 4)
    if f2 <= _SC_F2_FLOOR:
        si = 1.0

    yl = _clamp(si * ylonsi, 0.0, 1.0)

    # Minimum holdup from mass fraction
    rho_g_sg = s['rho_g'] / 62.37
    mflow_lpd = s['mflow_l'] * _SEC_PER_DAY
    mflow_gpd = s['mflow_g'] * _SEC_PER_DAY
    mflow_total_pd = mflow_lpd + mflow_gpd
    mass_frac_liq = mflow_lpd / mflow_total_pd if mflow_total_pd > 0 else 0.0
    min_l = (mass_frac_liq * rho_g_sg / (rho_g_sg + s['lsg_loc'])
             if (rho_g_sg + s['lsg_loc']) > 0 else 0.0)
    if yl < min_l:
        yl = min_l

    # Orkiszewski bubble flow correction
    vm = ugas + ul
    vs = _ORK_VS
    lb = _ORK_LB_A + _ORK_LB_B * vm * vm / s['diam_ft']
    if lb < _ORK_LB_MIN:
        lb = _ORK_LB_MIN
    if vm > 0:
        b_ratio = ugas / vm
        if b_ratio < lb:
            disc = (1.0 + vm / vs) ** 2 - 4.0 * ugas / vs
            if disc >= 0:
                yl = 1.0 - 0.5 * (1.0 + vm / vs - math.sqrt(disc))
                yl = _clamp(yl, 0.0, 1.0)

    # Reynolds and friction (mass-flow basis)
    mflow_pd = mflow_total_pd
    nre = _HB_RE_K * (mflow_pd / _SEC_PER_DAY) / (s['diam_ft'] *
          mul ** yl * s['mu_g'] ** (1.0 - yl))
    if nre < 2100:
        nre = 2100.0

    eond = s['rough'] / s['tid']
    f = _serghides_fanning(nre, eond)

    rho_avg = yl * rho_l + (1.0 - yl) * s['rho_g']

    fric_term = f * mflow_pd ** 2 / _HB_FRIC_K / s['diam_ft'] ** 5 / rho_avg
    sign = -1.0 if s['injection'] else 1.0
    return (rho_avg * math.sin(s['theta']) + sign * fric_term) / _IN2_PER_FT2


@rust_accelerated('hb_fbhp_gas_rust')
def _hb_fbhp_gas(thp, api, gsg, tid, rough, length, tht, bht,
                  wsg, qg_mmscfd, cgr, qw_bwpd, oil_vis,
                  injection=False, pr=0.0, theta=math.pi / 2.0):
    return _segment_march_gas(thp, api, gsg, tid, rough, length, tht, bht,
                              wsg, qg_mmscfd, cgr, qw_bwpd, oil_vis,
                              injection, pr, theta, _hb_gradient_gas)


@rust_accelerated('hb_fbhp_oil_rust')
def _hb_fbhp_oil(thp, api, gsg, tid, rough, length, tht, bht,
                  wsg, qt_stbpd, gor, wc, pb, rsb, sgsp,
                  rsb_scale=1.0, injection=False, theta=math.pi / 2.0,
                  vis_frac=1.0, rsb_frac=1.0):
    return _segment_march_oil(thp, api, gsg, tid, rough, length, tht, bht,
                              wsg, qt_stbpd, gor, wc, pb, rsb, sgsp,
                              rsb_scale, injection, theta, _hb_gradient_gas,
                              vis_frac, rsb_frac)


# ============================================================================
#  2. WOLDESEMAYAT-GHAJAR VLP
# ============================================================================

def _wg_void_fraction(u_sg, u_sl, rho_g, rho_l, sigma, diam, theta, p_sys):
    if u_sg <= 0:
        return 0.0
    if u_sl <= 0:
        return 1.0
    if rho_l <= rho_g or sigma <= 0 or diam <= 0 or p_sys <= 0:
        return 0.5
    u_m = u_sg + u_sl
    lam = u_sg / u_m
    density_ratio = (rho_g / rho_l) ** 0.1
    co = lam * (1.0 + (u_sl / u_sg) ** density_ratio)
    buoy = (_G_SI * diam * sigma * (1.0 + math.cos(theta)) *
            (rho_l - rho_g) / (rho_l * rho_l))
    u_gm_base = _WG_DRIFT_K * max(buoy, 0.0) ** 0.25
    p_exp = _P_ATM_PA / p_sys
    incl_factor = (_WG_INCL_K + _WG_INCL_K * math.sin(theta)) ** p_exp
    u_gm = u_gm_base * incl_factor
    denom = co * u_m + u_gm
    if denom <= 0:
        return 0.0
    return _clamp(u_sg / denom, 0.0, 1.0)


def _wg_friction_gradient_lm(m_flow_g, m_flow_l, rho_g, rho_l,
                              mu_g, mu_l, diam, rough):
    if diam <= 0:
        return 0.0
    area_si = math.pi * diam * diam / 4.0
    if area_si <= 0:
        return 0.0
    m_total = m_flow_g + m_flow_l
    if m_total < 1e-30:
        return 0.0
    mass_flux = m_total / area_si
    x = m_flow_g / m_total
    e_over_d = rough / diam
    if x <= 0:
        re_l = mass_flux * diam / mu_l
        f_l = _serghides_fanning(re_l, e_over_d)
        return 2.0 * f_l * mass_flux ** 2 / (diam * rho_l)
    if x >= 1.0:
        re_g = mass_flux * diam / mu_g
        f_g = _serghides_fanning(re_g, e_over_d)
        return 2.0 * f_g * mass_flux ** 2 / (diam * rho_g)
    g_l = mass_flux * (1.0 - x)
    g_g = mass_flux * x
    re_l = g_l * diam / mu_l
    re_g = g_g * diam / mu_g
    f_l = _serghides_fanning(re_l, e_over_d)
    f_g = _serghides_fanning(re_g, e_over_d)
    dpdz_l = 2.0 * f_l * g_l ** 2 / (diam * rho_l)
    dpdz_g = 2.0 * f_g * g_g ** 2 / (diam * max(rho_g, 1e-30))
    x_param = math.sqrt(dpdz_l / max(dpdz_g, 1e-30))
    if x_param < 1e-30:
        return dpdz_g
    liq_turb = re_l >= 2100
    gas_turb = re_g >= 2100
    if liq_turb and gas_turb:
        chisholm_c = _WG_CHIS_TT
    elif not liq_turb and gas_turb:
        chisholm_c = _WG_CHIS_LT
    elif liq_turb and not gas_turb:
        chisholm_c = _WG_CHIS_TL
    else:
        chisholm_c = _WG_CHIS_LL
    phi2_l = 1.0 + chisholm_c / x_param + 1.0 / (x_param * x_param)
    return phi2_l * dpdz_l


def _wg_gradient_gas(s):
    """Woldesemayat-Ghajar gradient: drift-flux void fraction, Chisholm friction (SI)."""
    # Convert oilfield state to SI
    diam_m = s['tid'] * _IN_TO_M
    rough_m = s['rough'] * _IN_TO_M
    area_m2 = math.pi * diam_m ** 2 / 4.0
    rho_g_kgm3 = s['rho_g'] * _LBFT3_TO_KGM3
    rho_l_kgm3 = s['rho_l'] * _LBFT3_TO_KGM3
    mflow_g_kgs = s['mflow_g'] * _LB_TO_KG
    mflow_l_kgs = s['mflow_l'] * _LB_TO_KG
    u_sg = mflow_g_kgs / (rho_g_kgm3 * area_m2)
    u_sl = mflow_l_kgs / (rho_l_kgm3 * area_m2)
    sigma_nm = s['sigma'] * _DYNECM_TO_NM
    p_sys_pa = s['p_avg'] * _PSI_TO_PA

    alpha_g = _wg_void_fraction(u_sg, u_sl, rho_g_kgm3, rho_l_kgm3,
                                 sigma_nm, diam_m, s['theta'], p_sys_pa)
    rho_mix = alpha_g * rho_g_kgm3 + (1.0 - alpha_g) * rho_l_kgm3
    dpdz_hydro = rho_mix * _G_SI * math.sin(s['theta'])

    mu_g_pas = s['mu_g'] * _CP_TO_PAS
    mu_l_pas = s['mu_l'] * _CP_TO_PAS
    dpdz_fric = _wg_friction_gradient_lm(
        mflow_g_kgs, mflow_l_kgs, rho_g_kgm3, rho_l_kgm3,
        mu_g_pas, mu_l_pas, diam_m, rough_m)

    sign = -1.0 if s['injection'] else 1.0
    dpdz_pam = dpdz_hydro + sign * dpdz_fric
    return dpdz_pam / (_PSI_TO_PA / _FT_TO_M)


@rust_accelerated('wg_fbhp_gas_rust')
def _wg_fbhp_gas(thp, api, gsg, tid, rough, length, tht, bht,
                  wsg, qg_mmscfd, cgr, qw_bwpd, oil_vis,
                  injection=False, pr=0.0, theta=math.pi / 2.0):
    return _segment_march_gas(thp, api, gsg, tid, rough, length, tht, bht,
                              wsg, qg_mmscfd, cgr, qw_bwpd, oil_vis,
                              injection, pr, theta, _wg_gradient_gas)


@rust_accelerated('wg_fbhp_oil_rust')
def _wg_fbhp_oil(thp, api, gsg, tid, rough, length, tht, bht,
                  wsg, qt_stbpd, gor, wc, pb, rsb, sgsp,
                  rsb_scale=1.0, injection=False, theta=math.pi / 2.0,
                  vis_frac=1.0, rsb_frac=1.0):
    return _segment_march_oil(thp, api, gsg, tid, rough, length, tht, bht,
                              wsg, qt_stbpd, gor, wc, pb, rsb, sgsp,
                              rsb_scale, injection, theta, _wg_gradient_gas,
                              vis_frac, rsb_frac)


# ============================================================================
#  3. GRAY VLP
# ============================================================================

def _gray_liquid_holdup(v_m, rho_l, rho_g, sigma, diam, lambda_l, p_sys):
    if lambda_l <= 0:
        return 0.0
    if lambda_l >= 1.0:
        return 1.0
    if v_m < 1e-10:
        return lambda_l

    sigma_lbf_ft = sigma * _DYNCM_TO_LBFFT
    if sigma_lbf_ft <= 0:
        return lambda_l

    rho_ns = rho_l * lambda_l + rho_g * (1.0 - lambda_l)
    if rho_ns <= 0:
        return lambda_l

    rv = lambda_l
    drho = max(rho_l - rho_g, 0.1)
    n1 = v_m ** 2 * rho_ns / (_G_FT * diam * drho)
    n2 = _G_FT * diam ** 2 * drho / sigma_lbf_ft

    a = _GRAY_HL_A * (1.0 - _GRAY_HL_B * math.log(1.0 + _GRAY_HL_C * rv / (1.0 + rv)))
    b = _GRAY_HL_B

    if n1 <= 1e-15:
        return lambda_l

    argument = (a + b * n2) / n1
    if argument > 0:
        f_e = -_GRAY_HL_B * math.log(argument)
    else:
        f_e = 0.0

    hl = 1.0 - (1.0 - lambda_l) * math.exp(f_e)
    return _clamp(hl, lambda_l, 1.0)


def _gray_effective_roughness(rough_dry, sigma, rho_ns, v_m, lambda_l):
    sigma_lbf_ft = sigma * _DYNCM_TO_LBFFT
    if v_m < 1e-10 or rho_ns <= 0 or sigma_lbf_ft <= 0:
        return rough_dry
    ke = rho_ns * v_m * v_m
    r1 = _GRAY_ROUGH_K * sigma_lbf_ft / max(ke, 1e-10)
    r2 = ke * lambda_l / sigma_lbf_ft
    film = r1 * (1.0 - math.exp(-r2))
    return rough_dry + max(film, 0.0)


# ============================================================================
#  Shared Gas Segment March
# ============================================================================

def _segment_march_gas(thp, api, gsg, tid, rough, length, tht, bht,
                       wsg, qg_mmscfd, cgr, qw_bwpd, oil_vis,
                       injection, pr, theta, gradient_fn):
    """Shared segment march for all gas VLP methods.

    gradient_fn(s) -> dpdz (psi/ft)
        Called at each pressure iteration with a dict containing all
        computed PVT/flow properties for the current segment step.
    """
    tc, pc = _sutton_tc_pc(gsg)
    osg = 141.5 / (api + 131.5)
    total_mass = (_RHO_AIR_STC * gsg * qg_mmscfd * 1e6 +
                  osg * _RHO_FW * cgr * qg_mmscfd * _FT3_PER_BBL +
                  wsg * _RHO_FW * qw_bwpd * _FT3_PER_BBL)
    if total_mass < 1e-6 or qg_mmscfd < 0.001:
        if qg_mmscfd >= 0.001:
            warnings.warn(
                f"Gas rate {qg_mmscfd:.4f} MMscf/d with near-zero total mass; "
                "using static gas column pressure.",
                RuntimeWarning, stacklevel=3
            )
        return _static_gas_column_pressure(thp, length, tht, bht, gsg, theta=theta)

    diam_ft = tid / 12.0
    rough_ft = rough / 12.0
    area = math.pi * diam_ft ** 2 / 4.0

    ndiv = _calc_segments(length)
    seg_len = length / ndiv

    mflow_g = _RHO_AIR_STC * gsg * qg_mmscfd * 1e6 / _SEC_PER_DAY
    mflow_w = wsg * _RHO_FW * qw_bwpd * _FT3_PER_BBL / _SEC_PER_DAY

    p_psia = thp

    for i in range(1, ndiv + 1):
        frac = (i - 0.5) / ndiv
        temp_f = tht + (bht - tht) * frac
        temp_r = temp_f + 459.67
        p_est = p_psia

        for _it in range(2):
            p_avg = max((p_psia + p_est) / 2.0, 14.7)

            cgr_loc, qo_loc, ql_loc, lsg_loc = _condensate_dropout(
                cgr, qg_mmscfd, p_avg, pr, osg, qw_bwpd, wsg)

            mflow_o = osg * _RHO_FW * qo_loc * _FT3_PER_BBL / _SEC_PER_DAY
            mflow_l = mflow_o + mflow_w
            mflow_total = mflow_g + mflow_l

            oil_vis_loc = _condensate_vis(pr, cgr_loc, gsg, api,
                                          temp_f, p_avg, oil_vis)

            zee = _z_factor(gsg, temp_f, p_avg, tc=tc, pc=pc)
            mu_g = _gas_viscosity(gsg, temp_f, p_avg, z=zee, tc=tc, pc=pc)

            rho_g = _MW_AIR * gsg * p_avg / (zee * _R_GAS * temp_r)
            rho_l = lsg_loc * _RHO_FW

            v_sg = (mflow_g / max(rho_g, 1e-10)) / area
            v_sl = (mflow_l / max(rho_l, 1e-10)) / area
            v_m = v_sg + v_sl
            lambda_l = v_sl / v_m if v_m > 1e-10 else 0.0
            rho_ns = rho_l * lambda_l + rho_g * (1.0 - lambda_l)

            water_visc = _water_viscosity(p_avg, temp_f)
            mu_l = ((qo_loc * oil_vis_loc + qw_bwpd * water_visc) / ql_loc
                    if ql_loc > 0 else oil_vis_loc)

            rs_est = _standing_rs(gsg, p_avg, temp_f, api)
            sigma = _interfacial_tension(p_avg, temp_f, api, rs_est,
                                         mflow_o * _LB_TO_KG,
                                         mflow_w * _LB_TO_KG,
                                         mflow_g * _LB_TO_KG)

            s = {
                'p_avg': p_avg, 'temp_f': temp_f, 'temp_r': temp_r,
                'zee': zee, 'mu_g': mu_g, 'rho_g': rho_g, 'rho_l': rho_l,
                'rho_ns': rho_ns, 'v_sg': v_sg, 'v_sl': v_sl, 'v_m': v_m,
                'lambda_l': lambda_l, 'sigma': sigma, 'mu_l': mu_l,
                'diam_ft': diam_ft, 'rough_ft': rough_ft, 'area': area,
                'injection': injection, 'theta': theta,
                'mflow_g': mflow_g, 'mflow_l': mflow_l,
                'mflow_o': mflow_o, 'mflow_w': mflow_w,
                'mflow_total': mflow_total,
                'qo_loc': qo_loc, 'ql_loc': ql_loc, 'qw_bwpd': qw_bwpd,
                'oil_vis_loc': oil_vis_loc, 'lsg_loc': lsg_loc,
                'osg': osg, 'gsg': gsg, 'tid': tid, 'rough': rough,
                'rs_est': rs_est, 'tc': tc, 'pc': pc,
                'qg_mmscfd': qg_mmscfd, 'api': api,
            }

            dpdz = gradient_fn(s)
            p_est = p_psia + dpdz * seg_len

        p_psia = p_est

    return p_psia


# ============================================================================
#  Shared Oil Segment March
# ============================================================================

def _segment_march_oil(thp, api, gsg, tid, rough, length, tht, bht,
                       wsg, qt_stbpd, gor, wc, pb, rsb, sgsp,
                       rsb_scale, injection, theta, gradient_fn,
                       vis_frac=1.0, rsb_frac=1.0):
    """Shared segment march for all oil VLP methods.

    gradient_fn(s) -> dpdz (psi/ft)
        Called at each pressure iteration with a dict containing all
        computed PVT/flow properties for the current segment step.
    """
    tc, pc = _sutton_tc_pc(gsg)
    wc_adj = max(wc, 1e-9)
    if qt_stbpd < 1e-7:
        return _static_oil_column_pressure(
            thp, length, tht, bht, wc_adj, wsg, api, sgsp, pb, rsb, rsb_scale,
            theta=theta)

    qo = qt_stbpd * (1.0 - wc)
    qw = qt_stbpd * wc
    osg = 141.5 / (api + 131.5)
    rsb_for_calc = rsb / rsb_scale

    diam_ft = tid / 12.0
    rough_ft = rough / 12.0
    area = math.pi * diam_ft ** 2 / 4.0

    ndiv = _calc_segments(length)
    seg_len = length / ndiv

    p_psia = thp

    for i in range(1, ndiv + 1):
        frac = (i - 0.5) / ndiv
        temp_f = tht + (bht - tht) * frac
        temp_r = temp_f + 459.67
        p_est = p_psia

        for _it in range(2):
            p_avg = max((p_psia + p_est) / 2.0, 14.7)

            rs_local = _velarde_rs(sgsp, api, temp_f, pb, rsb_for_calc, p_avg) * rsb_scale
            free_gas = max(gor - rs_local, 0.0)
            qg_mmscfd = max(free_gas * qo / 1e6, 1e-9)

            oil_vis_seg = _oil_viscosity_full(sgsp, api, temp_f, rsb, pb, p_avg,
                                              vis_frac, rsb_frac)
            rho_oil = _oil_density_mccain(rs_local, sgsp, osg,
                                           min(p_avg, pb), temp_f)

            zee = _z_factor(gsg, temp_f, p_avg, tc=tc, pc=pc)
            mu_g = _gas_viscosity(gsg, temp_f, p_avg, z=zee, tc=tc, pc=pc)

            rho_g = _MW_AIR * gsg * p_avg / (zee * _R_GAS * temp_r)

            mflow_o = osg * _RHO_FW * qo * _FT3_PER_BBL / _SEC_PER_DAY
            mflow_w = wsg * _RHO_FW * qw * _FT3_PER_BBL / _SEC_PER_DAY
            mflow_g = _RHO_AIR_STC * gsg * qg_mmscfd * 1e6 / _SEC_PER_DAY
            mflow_l = mflow_o + mflow_w
            mflow_total = mflow_g + mflow_l

            if mflow_total < 1e-15:
                tavg_r = (tht + bht) / 2.0 + 459.67
                p_psia *= math.exp(_GAS_COL_K * gsg * seg_len / (zee * tavg_r))
                break

            ql = qo + qw
            rho_w = wsg * _RHO_FW
            rho_l = (qo * rho_oil + qw * rho_w) / ql if ql > 0 else rho_oil
            lsg = (qo * osg + qw * wsg) / ql if ql > 0 else osg

            v_sg = (mflow_g / max(rho_g, 1e-10)) / area
            v_sl = (mflow_l / max(rho_l, 1e-10)) / area
            v_m = v_sg + v_sl
            lambda_l = v_sl / v_m if v_m > 1e-10 else 0.0
            rho_ns = rho_l * lambda_l + rho_g * (1.0 - lambda_l)

            water_visc = _water_viscosity(p_avg, temp_f)
            mu_l = ((qo * oil_vis_seg + qw * water_visc) / ql
                    if ql > 0 else oil_vis_seg)

            sigma = _interfacial_tension(p_avg, temp_f, api, rs_local,
                                         mflow_o * _LB_TO_KG,
                                         mflow_w * _LB_TO_KG,
                                         mflow_g * _LB_TO_KG)

            s = {
                'p_avg': p_avg, 'temp_f': temp_f, 'temp_r': temp_r,
                'zee': zee, 'mu_g': mu_g, 'rho_g': rho_g, 'rho_l': rho_l,
                'rho_ns': rho_ns, 'v_sg': v_sg, 'v_sl': v_sl, 'v_m': v_m,
                'lambda_l': lambda_l, 'sigma': sigma, 'mu_l': mu_l,
                'diam_ft': diam_ft, 'rough_ft': rough_ft, 'area': area,
                'injection': injection, 'theta': theta,
                'mflow_g': mflow_g, 'mflow_l': mflow_l,
                'mflow_o': mflow_o, 'mflow_w': mflow_w,
                'mflow_total': mflow_total,
                'qo_loc': qo, 'ql_loc': ql, 'qw_bwpd': qw,
                'oil_vis_loc': oil_vis_seg, 'lsg_loc': lsg,
                'osg': osg, 'gsg': gsg, 'tid': tid, 'rough': rough,
                'rs_est': rs_local, 'tc': tc, 'pc': pc,
                'qg_mmscfd': qg_mmscfd, 'api': api,
            }

            dpdz = gradient_fn(s)
            p_est = p_psia + dpdz * seg_len

        p_psia = p_est

    return p_psia


# ============================================================================
#  Gradient Callbacks
# ============================================================================

def _gray_gradient_gas(s):
    """Gray method gradient: holdup via effective roughness, acceleration term."""
    hl = _gray_liquid_holdup(s['v_m'], s['rho_l'], s['rho_g'], s['sigma'],
                             s['diam_ft'], s['lambda_l'], s['p_avg'])
    rho_s = s['rho_l'] * hl + s['rho_g'] * (1.0 - hl)

    eps_eff = _gray_effective_roughness(s['rough_ft'], s['sigma'], s['rho_ns'],
                                        s['v_m'], s['lambda_l'])
    eps_d = eps_eff / s['diam_ft']

    mu_ns = s['mu_l'] * s['lambda_l'] + s['mu_g'] * (1.0 - s['lambda_l'])
    mu_ns_lbfts = mu_ns * _CP_TO_LBFTS
    n_re = (s['rho_ns'] * s['v_m'] * s['diam_ft'] / mu_ns_lbfts
            if mu_ns_lbfts > 0 else 0.0)

    f_moody = 4.0 * _serghides_fanning(n_re, eps_d)

    gm = s['mflow_total'] / s['area']
    dpdz_hydro = rho_s * math.sin(s['theta']) / _IN2_PER_FT2
    dpdz_fric = (f_moody * gm ** 2 / (2.0 * _GC * s['diam_ft'] * s['rho_ns'] * _IN2_PER_FT2)
                 if s['rho_ns'] > 0 else 0.0)
    ek = gm * s['v_sg'] / (_GC * s['p_avg'] * _IN2_PER_FT2)
    denom = max(1.0 - ek, 0.1)

    sign = -1.0 if s['injection'] else 1.0
    return (dpdz_hydro + sign * dpdz_fric) / denom


@rust_accelerated('gray_fbhp_gas_rust')
def _gray_fbhp_gas(thp, api, gsg, tid, rough, length, tht, bht,
                    wsg, qg_mmscfd, cgr, qw_bwpd, oil_vis,
                    injection=False, pr=0.0, theta=math.pi / 2.0):
    return _segment_march_gas(thp, api, gsg, tid, rough, length, tht, bht,
                              wsg, qg_mmscfd, cgr, qw_bwpd, oil_vis,
                              injection, pr, theta, _gray_gradient_gas)


@rust_accelerated('gray_fbhp_oil_rust')
def _gray_fbhp_oil(thp, api, gsg, tid, rough, length, tht, bht,
                    wsg, qt_stbpd, gor, wc, pb, rsb, sgsp,
                    rsb_scale=1.0, injection=False, theta=math.pi / 2.0,
                    vis_frac=1.0, rsb_frac=1.0):
    return _segment_march_oil(thp, api, gsg, tid, rough, length, tht, bht,
                              wsg, qt_stbpd, gor, wc, pb, rsb, sgsp,
                              rsb_scale, injection, theta, _gray_gradient_gas,
                              vis_frac, rsb_frac)


# ============================================================================
#  4. BEGGS & BRILL VLP
# ============================================================================

_BB_SEGREGATED = 0
_BB_INTERMITTENT = 1
_BB_DISTRIBUTED = 2
_BB_TRANSITION = 3


def _bb_flow_pattern(froude, lambda_l):
    if lambda_l <= 0 or lambda_l >= 1.0:
        return _BB_DISTRIBUTED, 0.0
    l1 = _BB_L1_A * lambda_l ** _BB_L1_B
    l2 = _BB_L2_A * lambda_l ** _BB_L2_B
    l3 = _BB_L3_A * lambda_l ** _BB_L3_B
    l4 = _BB_L4_A * lambda_l ** _BB_L4_B
    if lambda_l < 0.01 and froude < l1:
        return _BB_SEGREGATED, 0.0
    if lambda_l >= 0.01 and froude < l2:
        return _BB_SEGREGATED, 0.0
    if lambda_l >= 0.01 and l2 <= froude <= l3:
        a = (l3 - froude) / (l3 - l2) if l3 > l2 else 0.5
        return _BB_TRANSITION, _clamp(a, 0.0, 1.0)
    if ((lambda_l >= 0.01 and froude > l3 and froude < l1) or
            (lambda_l < 0.01 and froude >= l1)):
        return _BB_INTERMITTENT, 0.0
    if lambda_l < 0.4 and froude >= l1:
        return _BB_DISTRIBUTED, 0.0
    if lambda_l >= 0.4 and froude > l4:
        return _BB_DISTRIBUTED, 0.0
    return _BB_INTERMITTENT, 0.0


def _bb_horizontal_holdup(lambda_l, froude, pattern):
    if lambda_l <= 0:
        return 0.0
    if lambda_l >= 1.0:
        return 1.0
    if froude <= 0:
        return lambda_l
    if pattern == _BB_SEGREGATED:
        a, b, c = _BB_HL_SEG
        hl0 = a * lambda_l ** b / froude ** c
    elif pattern == _BB_INTERMITTENT:
        a, b, c = _BB_HL_INT
        hl0 = a * lambda_l ** b / froude ** c
    elif pattern == _BB_DISTRIBUTED:
        a, b, c = _BB_HL_DIS
        hl0 = a * lambda_l ** b / froude ** c
    else:
        hl0 = lambda_l
    return max(hl0, lambda_l)


def _bb_inclination_correction(hl0, lambda_l, n_lv, froude, pattern,
                                theta=math.pi / 2.0):
    if lambda_l <= 0 or lambda_l >= 1.0:
        return hl0
    if pattern == _BB_SEGREGATED:
        e_p, f_p, g_p, h_p = _BB_IC_SEG
    elif pattern == _BB_INTERMITTENT:
        e_p, f_p, g_p, h_p = _BB_IC_INT
    elif pattern == _BB_DISTRIBUTED:
        return hl0
    else:
        return hl0
    arg = (e_p * lambda_l ** f_p *
           max(n_lv, 1e-10) ** g_p *
           max(froude, 1e-10) ** h_p)
    c_corr = (1.0 - lambda_l) * math.log(arg) if arg > 0 else 0.0
    c_corr = max(c_corr, 0.0)
    sin18 = math.sin(1.8 * theta)
    psi = 1.0 + c_corr * (sin18 - 0.333 * sin18 ** 3)
    return _clamp(hl0 * psi, lambda_l, 1.0)


def _bb_two_phase_friction(f_ns, lambda_l, hl_theta):
    if hl_theta <= 0:
        return f_ns
    y = lambda_l / (hl_theta ** 2)
    if y <= 0 or y == 1.0:
        s = 0.0
    elif 1.0 < y < 1.2:
        s = math.log(2.2 * y - 1.2)
    else:
        ln_y = math.log(y)
        denom = (_BB_SF_C0 + _BB_SF_C1 * ln_y + _BB_SF_C2 * ln_y ** 2 +
                 _BB_SF_C3 * ln_y ** 4)
        s = ln_y / denom if abs(denom) >= 1e-6 else 0.0
    s = _clamp(s, _BB_SF_LO, _BB_SF_HI)
    return f_ns * math.exp(s)


def _bb_gradient_gas(s):
    """Beggs & Brill gradient: flow pattern map, Payne correction, two-phase friction."""
    froude = s['v_m'] ** 2 / (_G_FT * s['diam_ft'])
    pattern, trans_a = _bb_flow_pattern(froude, s['lambda_l'])

    if pattern == _BB_TRANSITION:
        hl0_seg = _bb_horizontal_holdup(s['lambda_l'], froude, _BB_SEGREGATED)
        hl0_int = _bb_horizontal_holdup(s['lambda_l'], froude, _BB_INTERMITTENT)
        hl0 = trans_a * hl0_seg + (1.0 - trans_a) * hl0_int
    else:
        hl0 = _bb_horizontal_holdup(s['lambda_l'], froude, pattern)

    if pattern in (_BB_SEGREGATED, _BB_INTERMITTENT, _BB_TRANSITION):
        hl0 *= _BB_PAYNE
        hl0 = max(hl0, s['lambda_l'])

    sigma_lbf_ft = s['sigma'] * _DYNCM_TO_LBFFT
    n_lv = (_DG_VEL * s['v_sl'] * (s['rho_l'] / sigma_lbf_ft) ** 0.25
            if sigma_lbf_ft > 0 else 0.0)

    if pattern == _BB_TRANSITION:
        hl_seg = _bb_inclination_correction(
            _bb_horizontal_holdup(s['lambda_l'], froude, _BB_SEGREGATED) * _BB_PAYNE,
            s['lambda_l'], n_lv, froude, _BB_SEGREGATED, theta=s['theta'])
        hl_int = _bb_inclination_correction(
            _bb_horizontal_holdup(s['lambda_l'], froude, _BB_INTERMITTENT) * _BB_PAYNE,
            s['lambda_l'], n_lv, froude, _BB_INTERMITTENT, theta=s['theta'])
        hl_theta = trans_a * hl_seg + (1.0 - trans_a) * hl_int
    else:
        hl_theta = _bb_inclination_correction(
            hl0, s['lambda_l'], n_lv, froude, pattern, theta=s['theta'])

    hl_theta = _clamp(hl_theta, s['lambda_l'], 1.0)
    rho_s = s['rho_l'] * hl_theta + s['rho_g'] * (1.0 - hl_theta)

    mu_ns = s['mu_l'] * s['lambda_l'] + s['mu_g'] * (1.0 - s['lambda_l'])
    mu_ns_lbfts = mu_ns * _CP_TO_LBFTS
    n_re = (s['rho_ns'] * s['v_m'] * s['diam_ft'] / mu_ns_lbfts
            if mu_ns_lbfts > 0 else 0.0)

    eps_d = s['rough_ft'] / s['diam_ft']
    f_ns = _serghides_fanning(n_re, eps_d)
    f_tp = _bb_two_phase_friction(f_ns, s['lambda_l'], hl_theta)

    dpdz_hydro = rho_s * math.sin(s['theta']) / _IN2_PER_FT2
    dpdz_fric = (4.0 * f_tp * s['rho_ns'] * s['v_m'] ** 2 /
                 (2.0 * _GC * s['diam_ft'] * _IN2_PER_FT2) if s['rho_ns'] > 0 else 0.0)
    ek = rho_s * s['v_m'] * s['v_sg'] / (_GC * s['p_avg'] * _IN2_PER_FT2)
    denom = max(1.0 - ek, 0.1)

    sign = -1.0 if s['injection'] else 1.0
    return (dpdz_hydro + sign * dpdz_fric) / denom


@rust_accelerated('bb_fbhp_gas_rust')
def _bb_core_gas(thp, api, gsg, tid, rough, length, tht, bht,
                 wsg, qg_mmscfd, cgr, qw_bwpd, oil_vis,
                 injection=False, pr=0.0, theta=math.pi / 2.0):
    """Beggs & Brill core for gas wells."""
    return _segment_march_gas(thp, api, gsg, tid, rough, length, tht, bht,
                              wsg, qg_mmscfd, cgr, qw_bwpd, oil_vis,
                              injection, pr, theta, _bb_gradient_gas)


@rust_accelerated('bb_fbhp_oil_rust')
def _bb_core_oil(thp, api, gsg, tid, rough, length, tht, bht,
                 wsg, qt_stbpd, gor, wc, pb, rsb, sgsp,
                 rsb_scale=1.0, injection=False, theta=math.pi / 2.0,
                 vis_frac=1.0, rsb_frac=1.0):
    """Beggs & Brill core for oil wells."""
    return _segment_march_oil(thp, api, gsg, tid, rough, length, tht, bht,
                              wsg, qt_stbpd, gor, wc, pb, rsb, sgsp,
                              rsb_scale, injection, theta, _bb_gradient_gas,
                              vis_frac, rsb_frac)


# ============================================================================
#  Method Dispatch Dictionaries
# ============================================================================

_GAS_METHOD_DIC = {
    "HB": _hb_fbhp_gas,
    "WG": _wg_fbhp_gas,
    "GRAY": _gray_fbhp_gas,
    "BB": _bb_core_gas,
}

_OIL_METHOD_DIC = {
    "HB": _hb_fbhp_oil,
    "WG": _wg_fbhp_oil,
    "GRAY": _gray_fbhp_oil,
    "BB": _bb_core_oil,
}


# ============================================================================
#  Public API: fbhp
# ============================================================================

def fbhp(thp: float, completion: 'Completion', vlpmethod: str = 'WG', well_type: str = 'gas',
         gas_pvt=None, oil_pvt=None,
         qg_mmscfd: float = 0, cgr: float = 0, qw_bwpd: float = 0,
         oil_vis: float = 1.0, api: float = 45, pr: float = 0,
         qt_stbpd: float = 0, gor: float = 0, wc: float = 0,
         wsg: float = 1.07, injection: bool = False,
         gsg: float = 0.65, pb: float = 0, rsb: float = 0, sgsp: float = 0.65,
         metric: bool = False, return_profile: bool = False):
    """ Returns flowing bottom hole pressure (psia | barsa) using specified VLP correlation.

        thp: Tubing head pressure (psia | barsa)
        completion: Completion object describing the wellbore
        vlpmethod: VLP method - 'HB' (Hagedorn-Brown), 'WG' (Woldesemayat-Ghajar), 'GRAY', or 'BB' (Beggs & Brill)
        well_type: 'gas' or 'oil'

        Gas well parameters:
            qg_mmscfd: Gas rate (MMscf/d | sm3/d)
            cgr: Condensate-gas ratio (STB/MMscf | sm3/sm3)
            qw_bwpd: Water rate (STB/d | sm3/d)
            oil_vis: Oil (condensate) viscosity (cP). Defaults to 1.0
            api: Condensate API gravity. Defaults to 45
            pr: Reservoir pressure (psia | barsa) - for condensate dropout. 0 disables

        Oil well parameters:
            qt_stbpd: Total liquid rate (STB/d | sm3/d)
            gor: Gas-oil ratio (scf/STB | sm3/sm3)
            wc: Water cut (fraction 0-1)

        Common parameters:
            wsg: Water specific gravity. Defaults to 1.07
            injection: True for injection wells. Defaults to False
            gsg: Gas specific gravity (relative to air). Defaults to 0.65
            pb: Bubble point pressure (psia | barsa). Required for oil wells
            rsb: Solution GOR at Pb (scf/STB | sm3/sm3). Required for oil wells
            sgsp: Separator gas specific gravity. Defaults to 0.65
            metric: If True, inputs/outputs in Eclipse METRIC units. Default False.

        gas_pvt: GasPVT object (unused by VLP methods directly, reserved for future use)
        oil_pvt: OilPVT object. If provided for oil wells, extracts api, sgsp, pb, rsb from it
        return_profile: If True, returns a NodalResult with per-segment-boundary arrays
            ``md`` (cumulative MD), ``tvd`` (cumulative TVD), ``p`` (pressure at each
            boundary, including THP at the top and BHP at the bottom). Defaults to False
            (scalar BHP return).
    """
    if metric:
        thp = thp * BAR_TO_PSI
        if pr > 0:
            pr = pr * BAR_TO_PSI
        if well_type == 'gas':
            qg_mmscfd = qg_mmscfd * SM3_TO_MMSCF  # sm3/d -> MMscf/d
            if cgr > 0:
                cgr = cgr * SM3_PER_SM3_TO_STB_PER_MMSCF  # sm3/sm3 -> STB/MMscf
            if qw_bwpd > 0:
                qw_bwpd = qw_bwpd * SM3_TO_STB  # sm3/d -> STB/d
        else:
            if qt_stbpd > 0:
                qt_stbpd = qt_stbpd * SM3_TO_STB  # sm3/d -> STB/d
            if gor > 0:
                gor = gor * SM3_PER_SM3_TO_SCF_PER_STB  # sm3/sm3 -> scf/STB
            if pb > 0:
                pb = pb * BAR_TO_PSI
            if rsb > 0:
                rsb = rsb * SM3_PER_SM3_TO_SCF_PER_STB

    validate_pe_inputs(p=thp)
    validate_choice(well_type, ('gas', 'oil'), 'well_type')
    vlpmethod = validate_methods(["vlpmethod"], [vlpmethod])

    # Extract oil PVT parameters if provided (already in oilfield units from OilPVT)
    vis_frac = 1.0
    rsb_frac = 1.0
    if oil_pvt is not None and well_type == 'oil':
        api = oil_pvt.api
        sgsp = oil_pvt.sg_sp
        pb = oil_pvt.pb
        rsb = oil_pvt.rsb
        vis_frac = oil_pvt.vis_frac
        rsb_frac = oil_pvt.rsb_frac
        if oil_pvt.sg_g > 0:
            gsg = oil_pvt.sg_g

    def _run_section(thp_in, tid, rough, length, tht_seg, bht_seg, theta):
        if well_type == 'gas':
            return _GAS_METHOD_DIC[vlpmethod.name](
                thp=thp_in, api=api, gsg=gsg, tid=tid, rough=rough,
                length=length, tht=tht_seg, bht=bht_seg, wsg=wsg,
                qg_mmscfd=qg_mmscfd, cgr=cgr, qw_bwpd=qw_bwpd,
                oil_vis=oil_vis, injection=injection, pr=pr, theta=theta)
        else:
            return _OIL_METHOD_DIC[vlpmethod.name](
                thp=thp_in, api=api, gsg=gsg, tid=tid, rough=rough,
                length=length, tht=tht_seg, bht=bht_seg, wsg=wsg,
                qt_stbpd=qt_stbpd, gor=gor, wc=wc,
                pb=pb, rsb=rsb, sgsp=sgsp,
                rsb_scale=rsb_frac, injection=injection, theta=theta,
                vis_frac=vis_frac, rsb_frac=rsb_frac)

    # Loop over wellbore segments
    # Temperature interpolated over TVD (not MD) since geothermal
    # gradient is a function of vertical depth
    total_tvd = completion.total_tvd
    tvd_traversed = 0.0
    md_traversed = 0.0
    p_current = thp
    md_profile = [0.0]
    tvd_profile = [0.0]
    p_profile = [thp]

    for seg in completion.segments:
        frac_start = tvd_traversed / total_tvd if total_tvd > 0 else 0.0
        frac_end = (tvd_traversed + seg.tvd) / total_tvd if total_tvd > 0 else 1.0
        tht_seg = completion.tht + (completion.bht - completion.tht) * frac_start
        bht_seg = completion.tht + (completion.bht - completion.tht) * frac_end

        p_current = _run_section(p_current, seg.id, seg.roughness, seg.md,
                                 tht_seg, bht_seg, seg.theta)
        tvd_traversed += seg.tvd
        md_traversed += seg.md
        md_profile.append(md_traversed)
        tvd_profile.append(tvd_traversed)
        p_profile.append(p_current)

    if not math.isfinite(p_current):
        raise RuntimeError(f"VLP calculation produced non-finite BHP: {p_current}")

    if return_profile:
        md_arr = np.asarray(md_profile)
        tvd_arr = np.asarray(tvd_profile)
        p_arr = np.asarray(p_profile)
        if metric:
            md_arr = md_arr * FT_TO_M
            tvd_arr = tvd_arr * FT_TO_M
            p_arr = p_arr * PSI_TO_BAR
        return NodalResult({
            'md': md_arr,
            'tvd': tvd_arr,
            'p': p_arr,
            'bhp': p_arr[-1],
        })

    if metric:
        return p_current * PSI_TO_BAR

    return p_current


# ============================================================================
#  Public API: fthp
# ============================================================================

def fthp(bhp: float, completion: 'Completion', vlpmethod: str = 'WG', well_type: str = 'gas',
         gas_pvt=None, oil_pvt=None,
         qg_mmscfd: float = 0, cgr: float = 0, qw_bwpd: float = 0,
         oil_vis: float = 1.0, api: float = 45, pr: float = 0,
         qt_stbpd: float = 0, gor: float = 0, wc: float = 0,
         wsg: float = 1.07, injection: bool = False,
         gsg: float = 0.65, pb: float = 0, rsb: float = 0, sgsp: float = 0.65,
         metric: bool = False, thp_min: float = 14.7, thp_max: float = 20000.0,
         tol: float = 1e-3) -> float:
    """ Solves for tubing head pressure (psia | barsa) that produces the specified BHP
        under the given VLP correlation. Inverse of fbhp.

        bhp: Target flowing bottom hole pressure (psia | barsa)
        completion, vlpmethod, well_type, gas/oil flow parameters: same as fbhp()
        thp_min: Lower bracket for THP search (psia). Default atmospheric
        thp_max: Upper bracket for THP search (psia). Default 20,000 psi
        tol: Absolute convergence tolerance on THP (psia). Default 1e-3
        metric: If True, bhp and return are in barsa. Default False.
    """
    if metric:
        bhp_oilfield = bhp * BAR_TO_PSI
        thp_min_oilfield = thp_min * BAR_TO_PSI if thp_min != 14.7 else thp_min
        thp_max_oilfield = thp_max * BAR_TO_PSI if thp_max != 20000.0 else thp_max
    else:
        bhp_oilfield = bhp
        thp_min_oilfield = thp_min
        thp_max_oilfield = thp_max

    def _err(_args, thp_trial):
        calculated_bhp = fbhp(thp=thp_trial, completion=completion, vlpmethod=vlpmethod,
                              well_type=well_type, gas_pvt=gas_pvt, oil_pvt=oil_pvt,
                              qg_mmscfd=qg_mmscfd, cgr=cgr, qw_bwpd=qw_bwpd,
                              oil_vis=oil_vis, api=api, pr=pr,
                              qt_stbpd=qt_stbpd, gor=gor, wc=wc,
                              wsg=wsg, injection=injection,
                              gsg=gsg, pb=pb, rsb=rsb, sgsp=sgsp,
                              metric=False)
        return calculated_bhp - bhp_oilfield

    try:
        thp_solved = bisect_solve(None, _err, thp_min_oilfield, thp_max_oilfield, tol)
    except (RuntimeError, ValueError) as exc:
        raise RuntimeError(
            f"fthp could not bracket THP in [{thp_min_oilfield}, {thp_max_oilfield}] psi "
            f"for target BHP {bhp_oilfield} psi: {exc}"
        )

    if metric:
        return thp_solved * PSI_TO_BAR
    return thp_solved


# ============================================================================
#  Public API: outflow_curve
# ============================================================================

def outflow_curve(thp: float, completion: 'Completion', vlpmethod: str = 'WG',
                  well_type: str = 'gas',
                  gas_pvt=None, oil_pvt=None,
                  rates=None, n_points: int = 20, max_rate: Optional[float] = None,
                  cgr: float = 0, qw_bwpd: float = 0, oil_vis: float = 1.0,
                  api: float = 45, pr: float = 0,
                  gor: float = 0, wc: float = 0,
                  wsg: float = 1.07, injection: bool = False,
                  gsg: float = 0.65, pb: float = 0, rsb: float = 0, sgsp: float = 0.65,
                  metric: bool = False, n_rates: Optional[int] = None) -> 'NodalResult':
    """ Returns VLP outflow curve as a dictionary.

        thp: Tubing head pressure (psia | barsa)
        completion: Completion object
        vlpmethod: VLP method string
        well_type: 'gas' or 'oil'
        rates: List of rates to evaluate (MMscf/d | sm3/d for gas, STB/d | sm3/d for oil). If None, auto-generated
        n_points: Number of rate points if rates is None. Default 20.
        n_rates: Deprecated alias for n_points (kept for backward compatibility; takes precedence if both given).
        max_rate: Maximum rate for auto-generation
        Other parameters: Same as fbhp()
        metric: If True, inputs/outputs in Eclipse METRIC units. Default False.

        Returns:
            dict with keys:
                'rate': list of flow rates (MMscf/d for gas, STB/d for oil; sm3/d if metric)
                'rates': alias for 'rate' (kept for backward compatibility)
                'bhp': list of flowing BHP values (psia; barsa if metric) at each rate
    """
    validate_choice(well_type, ('gas', 'oil'), 'well_type')
    # Convert metric inputs to oilfield at the boundary
    if metric:
        thp = thp * BAR_TO_PSI
        if pr > 0:
            pr = pr * BAR_TO_PSI
        if well_type == 'gas':
            if cgr > 0:
                cgr = cgr * SM3_PER_SM3_TO_STB_PER_MMSCF
            if qw_bwpd > 0:
                qw_bwpd = qw_bwpd * SM3_TO_STB
        else:
            if gor > 0:
                gor = gor * SM3_PER_SM3_TO_SCF_PER_STB
            if pb > 0:
                pb = pb * BAR_TO_PSI
            if rsb > 0:
                rsb = rsb * SM3_PER_SM3_TO_SCF_PER_STB

    if oil_pvt is not None and well_type == 'oil':
        api = oil_pvt.api
        sgsp = oil_pvt.sg_sp
        pb = oil_pvt.pb
        rsb = oil_pvt.rsb
        if oil_pvt.sg_g > 0:
            gsg = oil_pvt.sg_g

    if rates is None:
        if max_rate is None:
            if metric:
                # Default max rates in metric units
                max_rate = 50.0 * MMSCF_TO_SM3 if well_type == 'gas' else 10000.0 * STB_TO_SM3
            else:
                max_rate = 50.0 if well_type == 'gas' else 10000.0
        if metric:
            min_rate = 0.01 * MMSCF_TO_SM3 if well_type == 'gas' else 1.0 * STB_TO_SM3
        else:
            min_rate = 0.01 if well_type == 'gas' else 1.0
        n = n_rates if n_rates is not None else n_points
        rates = list(np.linspace(min_rate, max_rate, n))

    bhp_list = []
    for rate in rates:
        if well_type == 'gas':
            # Convert rate to MMscf/d for internal fbhp call
            rate_mmscfd = rate * SM3_TO_MMSCF if metric else rate
            bhp_val = fbhp(thp=thp, completion=completion, vlpmethod=vlpmethod,
                           well_type='gas', qg_mmscfd=rate_mmscfd, cgr=cgr,
                           qw_bwpd=qw_bwpd, oil_vis=oil_vis, api=api, pr=pr,
                           wsg=wsg, injection=injection, gsg=gsg)
        else:
            # Convert rate to STB/d for internal fbhp call
            rate_stbpd = rate * SM3_TO_STB if metric else rate
            bhp_val = fbhp(thp=thp, completion=completion, vlpmethod=vlpmethod,
                           well_type='oil', qt_stbpd=rate_stbpd, gor=gor, wc=wc,
                           wsg=wsg, injection=injection, gsg=gsg,
                           pb=pb, rsb=rsb, sgsp=sgsp, api=api, oil_pvt=oil_pvt)
        # bhp_val is in psia (fbhp called without metric=True)
        if metric:
            bhp_val = bhp_val * PSI_TO_BAR
        bhp_list.append(bhp_val)

    rates_out = list(rates)
    return NodalResult({'rate': rates_out, 'rates': rates_out, 'bhp': bhp_list})


# ============================================================================
#  Public API: ipr_curve
# ============================================================================

def ipr_curve(reservoir: 'Reservoir', well_type: str = 'gas',
              gas_pvt=None, oil_pvt=None,
              n_points: int = 20, min_pwf: Optional[float] = None,
              wc: float = 0, wsg: float = 1.07, bo: float = 1.2,
              uo: float = 1.0, gsg: float = 0.65,
              metric: bool = False) -> 'NodalResult':
    """ Returns IPR curve as dict {'pwf': [...], 'rate': [...]}.

        reservoir: Reservoir object (constructed with matching metric flag)
        well_type: 'gas', 'oil', or 'water'
        gas_pvt: GasPVT object (required for gas wells if not using defaults)
        oil_pvt: OilPVT object (optional for oil wells)
        n_points: Number of pressure points
        min_pwf: Minimum flowing BHP (psia | barsa). Defaults to 14.7 psia (1.01325 barsa)
        wc: Water cut (fraction 0-1). For oil wells
        wsg: Water specific gravity
        bo: Oil FVF (rb/stb | rm3/sm3). Used for oil/water wells if oil_pvt not provided
        uo: Oil viscosity (cP). Used for oil/water wells if oil_pvt not provided
        gsg: Gas specific gravity. Used if gas_pvt not provided
        metric: If True, inputs/outputs in Eclipse METRIC units. Default False.
    """
    validate_choice(well_type, ('gas', 'oil', 'water'), 'well_type')
    if min_pwf is None:
        min_pwf = 1.01325 if metric else 14.7
    if metric:
        min_pwf = min_pwf * BAR_TO_PSI

    # Reservoir stores oilfield units internally (converted in constructor)
    pr = reservoir.pr
    degf = reservoir.degf
    k = reservoir.k
    h = reservoir.h
    re = reservoir.re
    rw = reservoir.rw
    S = reservoir.S
    D = reservoir.D

    pwf_list = list(np.linspace(min_pwf, pr, n_points))
    rate_list = []

    if well_type == 'gas':
        sg = gas_pvt.sg if gas_pvt else gsg
        co2 = gas_pvt.co2 if gas_pvt else 0
        h2s = gas_pvt.h2s if gas_pvt else 0
        n2 = gas_pvt.n2 if gas_pvt else 0
        h2_frac = gas_pvt.h2 if gas_pvt else 0

        for pwf in pwf_list:
            qg = gas.gas_rate_radial(
                k=k, h=h, pr=pr, pwf=pwf, r_w=rw, r_ext=re,
                degf=degf, S=S, D=D, sg=sg,
                co2=co2, h2s=h2s, n2=n2, h2=h2_frac)
            rate_list.append(float(qg))  # Mscf/d in oilfield

    elif well_type == 'oil':
        # Darcy above Pb, Vogel below Pb
        if oil_pvt is not None:
            pb = oil_pvt.pb
            rsb = oil_pvt.rsb
            rs_pr = oil_pvt.rs(pr, degf)
            bo_pr = oil_pvt.bo(pr, degf, rs=rs_pr)
            uo_pr = oil_pvt.viscosity(pr, degf, rs=rs_pr)
        else:
            pb = 1e6  # No bubble point specified - all Darcy
            bo_pr = bo
            uo_pr = uo

        J = _DARCY_K * k * h / (uo_pr * bo_pr * (np.log(re / rw) + S - _DIETZ_CORR))

        for pwf in pwf_list:
            if oil_pvt is not None and pr > pb:
                # Undersaturated: Darcy from Pr to Pb, Vogel below Pb
                q_darcy = J * (pr - max(pwf, pb))
                if pwf < pb:
                    qsat_max = J * pb / _VOGEL_AOF_DENOM
                    q_vogel = qsat_max * (1 - _VOGEL_LIN * (pwf / pb) - _VOGEL_QUAD * (pwf / pb) ** 2)
                    rate_list.append(q_darcy + q_vogel)
                else:
                    rate_list.append(q_darcy)
            elif oil_pvt is not None and pr <= pb:
                # Saturated: Vogel
                qsat_max = J * pr / _VOGEL_AOF_DENOM
                q = qsat_max * (1 - _VOGEL_LIN * (pwf / pr) - _VOGEL_QUAD * (pwf / pr) ** 2)
                rate_list.append(q)
            else:
                # Simple Darcy (no Pb info)
                rate_list.append(J * (pr - pwf))

    elif well_type == 'water':
        # Water injectivity
        J = _DARCY_K * k * h / (uo * bo * (np.log(re / rw) + S - _DIETZ_CORR))
        for pwf in pwf_list:
            rate_list.append(J * (pr - pwf))

    # Convert outputs to metric if requested
    if metric:
        pwf_list = [p * PSI_TO_BAR for p in pwf_list]
        if well_type == 'gas':
            rate_list = [r * MSCF_TO_SM3 for r in rate_list]  # Mscf/d -> sm3/d
        elif well_type in ('oil', 'water'):
            rate_list = [r * STB_TO_SM3 for r in rate_list]  # STB/d -> sm3/d

    return NodalResult({'pwf': pwf_list, 'rate': rate_list})


# ============================================================================
#  Public API: operating_point
# ============================================================================

def operating_point(thp: float, completion: 'Completion', reservoir: 'Reservoir',
                    vlpmethod: str = 'WG', well_type: str = 'gas',
                    gas_pvt=None, oil_pvt=None,
                    cgr: float = 0, qw_bwpd: float = 0, oil_vis: float = 1.0,
                    api: float = 45,
                    gor: float = 0, wc: float = 0,
                    wsg: float = 1.07, injection: bool = False,
                    gsg: float = 0.65, pb: float = 0, rsb: float = 0, sgsp: float = 0.65,
                    bo: float = 1.2, uo: float = 1.0,
                    n_points: int = 25,
                    metric: bool = False) -> 'NodalResult':
    """ Returns operating point as dict with 'rate', 'bhp', 'vlp', 'ipr', 'converged'.

        Finds intersection of VLP outflow curve and IPR inflow curve via bisection.

        thp: Tubing head pressure (psia | barsa)
        completion: Completion object
        reservoir: Reservoir object (constructed with matching metric flag)
        vlpmethod: VLP method string
        well_type: 'gas' or 'oil'
        injection: True for injection wells — forwarded to the internal fbhp and
            outflow_curve calls. Defaults to False.
        Other parameters: Same as fbhp() and ipr_curve()
        metric: If True, inputs/outputs in Eclipse METRIC units. Default False.

        Returns:
            rate: Operating rate (MMscf/d | sm3/d for gas, STB/d | sm3/d for oil)
            bhp: Operating BHP (psia | barsa)
            vlp: VLP outflow curve dict
            ipr: IPR inflow curve dict
            converged: True if VLP/IPR intersection was found; False if bisection failed
                and result is a fallback (rate=0, bhp=pr).
    """
    validate_choice(well_type, ('gas', 'oil'), 'well_type')
    # Convert metric inputs to oilfield at the boundary
    if metric:
        thp = thp * BAR_TO_PSI
        if well_type == 'gas':
            if cgr > 0:
                cgr = cgr * SM3_PER_SM3_TO_STB_PER_MMSCF
            if qw_bwpd > 0:
                qw_bwpd = qw_bwpd * SM3_TO_STB
        else:
            if gor > 0:
                gor = gor * SM3_PER_SM3_TO_SCF_PER_STB
            if pb > 0:
                pb = pb * BAR_TO_PSI
            if rsb > 0:
                rsb = rsb * SM3_PER_SM3_TO_SCF_PER_STB

    # Reservoir stores oilfield units internally (converted in constructor)
    pr = reservoir.pr

    if oil_pvt is not None and well_type == 'oil':
        api = oil_pvt.api
        sgsp = oil_pvt.sg_sp
        pb = oil_pvt.pb
        rsb = oil_pvt.rsb
        if oil_pvt.sg_g > 0:
            gsg = oil_pvt.sg_g

    # Get IPR curve (in oilfield units — no metric flag)
    ipr = ipr_curve(reservoir=reservoir, well_type=well_type,
                    gas_pvt=gas_pvt, oil_pvt=oil_pvt,
                    n_points=n_points, wc=wc, wsg=wsg,
                    bo=bo, uo=uo, gsg=gsg)

    # Find AOF (maximum rate at minimum pwf)
    aof = max(ipr['rate'])
    if aof <= 0:
        result = {'rate': 0.0, 'bhp': pr, 'vlp': {'rate': [], 'rates': [], 'bhp': []},
                  'ipr': ipr, 'converged': False}
        if metric:
            result['bhp'] = pr * PSI_TO_BAR
            result['ipr'] = _convert_ipr_to_metric(ipr, well_type)
        return NodalResult(result)

    # IPR gas rates are in Mscf/d; VLP uses MMscf/d. Scale factor:
    gas_scale = 1000.0 if well_type == 'gas' else 1.0

    # Define error function: VLP BHP - IPR BHP at given rate
    # rate is in IPR units (Mscf/d for gas, STB/d for oil)
    def _err(args, rate):
        if rate <= 0:
            return -1.0  # At zero rate, VLP BHP < IPR BHP (Pr)

        # VLP BHP (convert gas rate to MMscf/d)
        if well_type == 'gas':
            vlp_bhp = fbhp(thp=thp, completion=completion, vlpmethod=vlpmethod,
                           well_type='gas', qg_mmscfd=rate / gas_scale, cgr=cgr,
                           qw_bwpd=qw_bwpd, oil_vis=oil_vis, api=api, pr=pr,
                           wsg=wsg, injection=injection, gsg=gsg)
        else:
            vlp_bhp = fbhp(thp=thp, completion=completion, vlpmethod=vlpmethod,
                           well_type='oil', qt_stbpd=rate, gor=gor, wc=wc,
                           wsg=wsg, injection=injection, gsg=gsg, pb=pb, rsb=rsb, sgsp=sgsp,
                           api=api, oil_pvt=oil_pvt)

        # IPR BHP: interpolate from IPR curve
        ipr_rates = ipr['rate']
        ipr_pwfs = ipr['pwf']
        # Rate decreases as pwf increases, so interpolate
        ipr_bhp = np.interp(rate, ipr_rates[::-1], ipr_pwfs[::-1])

        return vlp_bhp - ipr_bhp

    # Bisect to find operating rate (in IPR units).
    # VLP curves for oil wells can be non-monotonic at very low rates (spurious
    # near-shut-in high BHP from holdup correlations), so a naive bisect between
    # min_rate and AOF may miss the true crossing. Scan the rate range first and
    # pick the sign-change closest to AOF — that's the physical operating point.
    min_rate = 0.1 * gas_scale if well_type == 'gas' else 1.0
    max_rate_search = aof * 0.999

    converged = True
    scan_rates = list(np.linspace(min_rate, max_rate_search, max(n_points, 25)))
    scan_errs = [_err(None, r) for r in scan_rates]

    # Find all sign changes; prefer the one at the highest rate
    brackets = [
        (scan_rates[i], scan_rates[i + 1])
        for i in range(len(scan_rates) - 1)
        if scan_errs[i] * scan_errs[i + 1] < 0
    ]

    if brackets:
        lo, hi = brackets[-1]
        try:
            op_rate = bisect_solve(None, _err, lo, hi, 1e-4)
        except (RuntimeError, ValueError):
            op_rate = 0.0
            converged = False
    else:
        op_rate = 0.0
        converged = False

    # Calculate operating BHP
    if well_type == 'gas' and op_rate > 0:
        op_bhp = fbhp(thp=thp, completion=completion, vlpmethod=vlpmethod,
                      well_type='gas', qg_mmscfd=op_rate / gas_scale, cgr=cgr,
                      qw_bwpd=qw_bwpd, oil_vis=oil_vis, api=api, pr=pr,
                      wsg=wsg, injection=injection, gsg=gsg)
    elif well_type == 'oil' and op_rate > 0:
        op_bhp = fbhp(thp=thp, completion=completion, vlpmethod=vlpmethod,
                      well_type='oil', qt_stbpd=op_rate, gor=gor, wc=wc,
                      wsg=wsg, injection=injection, gsg=gsg, pb=pb, rsb=rsb, sgsp=sgsp,
                      api=api, oil_pvt=oil_pvt)
    else:
        op_bhp = pr

    # Generate VLP curve for output (in VLP units: MMscf/d for gas, STB/d for oil)
    vlp_max = max_rate_search / gas_scale
    vlp_min = 0.1 if well_type == 'gas' else 1.0
    vlp_rates = list(np.linspace(vlp_min, vlp_max, n_points))
    vlp = outflow_curve(thp=thp, completion=completion, vlpmethod=vlpmethod,
                        well_type=well_type, rates=vlp_rates,
                        cgr=cgr, qw_bwpd=qw_bwpd, oil_vis=oil_vis, api=api,
                        pr=pr, gor=gor, wc=wc, wsg=wsg, injection=injection, gsg=gsg,
                        pb=pb, rsb=rsb, sgsp=sgsp, oil_pvt=oil_pvt)

    # Convert operating rate to VLP units for return
    op_rate_out = op_rate / gas_scale if well_type == 'gas' else op_rate

    # Convert outputs to metric
    if metric:
        op_bhp = op_bhp * PSI_TO_BAR
        if well_type == 'gas':
            op_rate_out = op_rate_out * MMSCF_TO_SM3  # MMscf/d -> sm3/d
        else:
            op_rate_out = op_rate_out * STB_TO_SM3  # STB/d -> sm3/d
        vlp = _convert_vlp_to_metric(vlp, well_type)
        ipr = _convert_ipr_to_metric(ipr, well_type)

    return NodalResult({'rate': op_rate_out, 'bhp': op_bhp,
                        'vlp': vlp, 'ipr': ipr,
                        'converged': converged})


def _convert_vlp_to_metric(vlp, well_type):
    """Convert VLP outflow curve dict from oilfield to metric units."""
    if well_type == 'gas':
        rates = [r * MMSCF_TO_SM3 for r in vlp['rate']]  # MMscf/d -> sm3/d
    else:
        rates = [r * STB_TO_SM3 for r in vlp['rate']]  # STB/d -> sm3/d
    bhps = [b * PSI_TO_BAR for b in vlp['bhp']]  # psia -> barsa
    return {'rate': rates, 'rates': rates, 'bhp': bhps}


def _convert_ipr_to_metric(ipr, well_type):
    """Convert IPR curve dict from oilfield to metric units."""
    pwfs = [p * PSI_TO_BAR for p in ipr['pwf']]  # psia -> barsa
    if well_type == 'gas':
        rates = [r * MSCF_TO_SM3 for r in ipr['rate']]  # Mscf/d -> sm3/d
    else:
        rates = [r * STB_TO_SM3 for r in ipr['rate']]  # STB/d -> sm3/d
    return {'pwf': pwfs, 'rate': rates}
