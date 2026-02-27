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

import math
import numpy as np

from pyrestoolbox.classes import vlp_method, class_dic
from pyrestoolbox.validate import validate_methods
from pyrestoolbox.shared_fns import bisect_solve
import pyrestoolbox.gas as gas
import pyrestoolbox.oil as oil

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

        md: Measured depth of this segment (ft)
        id: Internal diameter (inches)
        deviation: Deviation from vertical (degrees). 0=vertical, 90=horizontal. Defaults to 0
        roughness: Pipe roughness (inches). Defaults to 0.0006
    """
    def __init__(self, md, id, deviation=0, roughness=0.0006):
        self.md = md
        self.id = id
        self.deviation = deviation
        self.roughness = roughness

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
            tid: Tubing ID (inches)
            length: Tubing length (ft) - from wellhead to tubing shoe
            tht: Tubing head (wellhead) temperature (degF)
            bht: Bottom hole temperature (degF)
            rough: Tubing roughness (inches). Defaults to 0.0006
            cid: Casing ID (inches) below tubing shoe. Defaults to 0 (no casing section)
            crough: Casing roughness (inches). Defaults to 0.0006
            mpd: Mid-perforation depth (ft). Defaults to length (no casing section)

        Segment mode (keyword):
            segments: List of WellSegment objects defining the wellbore
            tht: Tubing head (wellhead) temperature (degF)
            bht: Bottom hole temperature (degF)
    """
    def __init__(self, tid=None, length=None, tht=None, bht=None, rough=0.0006,
                 cid=0, crough=0.0006, mpd=0, segments=None):
        if segments is not None:
            # Segment mode
            if tht is None or bht is None:
                raise ValueError("tht and bht are required when using segments")
            if not segments:
                raise ValueError("segments list must not be empty")
            self._segments = list(segments)
            self.tht = tht
            self.bht = bht
            # Set legacy attributes from first segment for compatibility
            self.tid = self._segments[0].id
            self.rough = self._segments[0].roughness
            self.length = self._segments[0].md
            self.cid = 0
            self.crough = crough
            self.mpd = self.total_md
        elif tid is not None and length is not None:
            # Legacy mode
            if tht is None or bht is None:
                raise ValueError("tht and bht are required")
            self.tid = tid
            self.length = length
            self.tht = tht
            self.bht = bht
            self.rough = rough
            self.cid = cid
            self.crough = crough
            self.mpd = mpd if mpd > 0 else length
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


class Reservoir:
    """ Reservoir description for IPR calculations.

        pr: Reservoir pressure (psia)
        degf: Reservoir temperature (degF)
        k: Permeability (mD)
        h: Net pay thickness (ft)
        re: Drainage radius (ft)
        rw: Wellbore radius (ft)
        S: Skin factor. Defaults to 0
        D: Non-Darcy coefficient (day/mscf for gas, 0 for oil). Defaults to 0
    """
    def __init__(self, pr, degf, k, h, re, rw, S=0, D=0):
        self.pr = pr
        self.degf = degf
        self.k = k
        self.h = h
        self.re = re
        self.rw = rw
        self.S = S
        self.D = D


# ============================================================================
#  Simplified PVT helpers (standalone, no GasPVT/OilPVT dependency)
#  Used inside VLP segment loops for performance
# ============================================================================

def _z_factor(sg, temp_f, press_psia):
    """Z-factor via Hall-Yarborough (1973)."""
    if press_psia < 1.0:
        return 1.0
    tpc = 169.2 + 349.5 * sg - 74.0 * sg * sg
    ppc = 756.8 - 131.0 * sg - 3.6 * sg * sg
    tr = (temp_f + 459.67) / tpc
    pr = press_psia / ppc

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


def _gas_viscosity(sg, temp_f, press_psia):
    """Gas viscosity via Lee-Gonzalez-Eakin (1966) (cP)."""
    z = _z_factor(sg, temp_f, press_psia)
    temp_r = temp_f + 459.67
    mw = _MW_AIR * sg
    rho_gcc = (_MW_AIR * sg * press_psia / (z * 10.732 * temp_r)) / 62.428

    k_val = ((9.4 + 0.02 * mw) * temp_r ** 1.5 / (209.0 + 19.0 * mw + temp_r))
    x_val = 3.5 + 986.0 / temp_r + 0.01 * mw
    y_val = 2.4 - 0.2 * x_val
    exp_arg = x_val * rho_gcc ** y_val
    if exp_arg > 500:
        exp_arg = 500  # Cap to avoid overflow
    return max(1e-4 * k_val * math.exp(exp_arg), 1e-6)


def _gas_density(sg, temp_f, press_psia, z=None):
    """Gas density in lb/ft^3."""
    if z is None:
        z = _z_factor(sg, temp_f, press_psia)
    return _MW_AIR * sg * press_psia / (z * 10.732 * (temp_f + 459.67))


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


def _oil_viscosity_full(sgsp, api, temp_f, rsb, pb, press_psia):
    """Oil viscosity (Beggs-Robinson + undersaturated correction) (cP)."""
    rs = _velarde_rs(sgsp, api, temp_f, pb, rsb, press_psia)
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

    return max(mu_or, 0.001)


# ============================================================================
#  IFT Correlations
# ============================================================================

def _dead_oil_ift(api, temp_f):
    temp_c = (temp_f - 32.0) / 1.8
    a = 1.11591 - 0.00305 * temp_c
    return max(a * (38.085 - 0.259 * api), 1.0)


def _gas_oil_ift(api, temp_f, rs_scf_stb):
    sigma_od = _dead_oil_ift(api, temp_f)
    if rs_scf_stb <= 0:
        return sigma_od
    rs_m3m3 = rs_scf_stb * _SCF_STB_TO_M3M3
    if rs_m3m3 < 50.0:
        ratio = 1.0 / (1.0 + 0.02549 * rs_m3m3 ** 1.0157)
    else:
        ratio = 32.0436 * rs_m3m3 ** (-1.1367)
    return max(sigma_od * _clamp(ratio, 0.0, 1.0), 1.0)


def _gas_water_ift(press_psia, temp_f):
    p = max(press_psia, 14.7)
    sw74 = 75.0 - 1.108 * p ** 0.349
    sw280 = 53.0 - 0.1048 * p ** 0.637
    if temp_f <= 74.0:
        return max(sw74, 1.0)
    elif temp_f >= 280.0:
        return max(sw280, 1.0)
    return max(sw74 + (sw280 - sw74) * (temp_f - 74.0) / (280.0 - 74.0), 1.0)


def _interfacial_tension(p_avg, temp_f, api, rs_scf_stb,
                         m_flow_o, m_flow_w, m_flow_g):
    sigma_og = _gas_oil_ift(api, temp_f, rs_scf_stb)
    sigma_wg = _gas_water_ift(p_avg, temp_f)
    sigma_ow = 26.0
    denom = m_flow_o * m_flow_w + m_flow_o * m_flow_g + m_flow_g * m_flow_w
    if denom <= 0:
        return 20.0
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
    a = -2.0 * _log10(eps_d / 3.7 + 12.0 / re)
    b = -2.0 * _log10(eps_d / 3.7 + 2.51 * a / re)
    c = -2.0 * _log10(eps_d / 3.7 + 2.51 * b / re)
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
    qo_local = max(cgr_local * qg_mmscfd, 1e-7)
    ql_local = qo_local + qw_bwpd
    lsg_local = (qo_local * osg + qw_bwpd * wsg) / ql_local
    return cgr_local, qo_local, ql_local, lsg_local


def _condensate_vis(pr, cgr_local, gsg, api, temp_f, p_avg, oil_vis):
    if pr > 14.7 and cgr_local > 0.01:
        return _oil_viscosity_full(gsg, api, temp_f, rsb=0.0, pb=14.7,
                                   press_psia=p_avg)
    return oil_vis


def _static_gas_column_pressure(thp, length, tht, bht, gsg, theta=math.pi / 2.0):
    n_seg = 50
    d_len = length / n_seg
    sin_theta = math.sin(theta)
    p = thp
    for i in range(n_seg):
        frac = (i + 0.5) / n_seg
        temp_local = tht + (bht - tht) * frac
        temp_r = temp_local + 459.67
        zee = _z_factor(gsg, temp_local, max(p, 14.7))
        rho_gas = _MW_AIR * gsg * p / (zee * 10.732 * temp_r)
        p += rho_gas / 144.0 * d_len * sin_theta
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
        oil_sg_local = rho_oil / 62.4
        mix_sg = (1.0 - wc) * oil_sg_local + wc * wsg
        p += 0.433 * mix_sg * d_len * sin_theta
    return p


# ============================================================================
#  1. HAGEDORN-BROWN VLP
# ============================================================================

def _hb_fbhp_gas(thp, api, gsg, tid, rough, length, tht, bht,
                  wsg, qg_mmscfd, cgr, qw_bwpd, oil_vis,
                  injection=False, pr=0.0, theta=math.pi / 2.0):
    osg = 141.5 / (api + 131.5)
    total_mass = (0.0765 * gsg * qg_mmscfd * 1e6 +
                  osg * 62.4 * cgr * qg_mmscfd * 5.615 +
                  wsg * 62.4 * qw_bwpd * 5.615)
    if total_mass < 1e-6 or qg_mmscfd < 0.05:
        return _static_gas_column_pressure(thp, length, tht, bht, gsg, theta=theta)

    area = math.pi * tid * tid / 4.0 / 144.0
    ndiv = _calc_segments(length)
    seg_len = length / ndiv

    depth = [0.0] * (ndiv + 1)
    dz = [0.0] * (ndiv + 1)
    incr = [0.0] * (ndiv + 1)
    temp = [0.0] * (ndiv + 1)
    p = [0.0] * (ndiv + 1)

    temp[0] = tht
    midpoint = 0.0

    for i in range(ndiv + 1):
        depth[i] = i * seg_len
        if i == 0:
            incr[i] = seg_len
        else:
            incr[i] = depth[i] - depth[i - 1]
            dz[i] = (incr[i - 1] + incr[i]) / 2.0
        if i == ndiv:
            incr[i] = 0.5 * seg_len
            dz[i] = incr[i - 1] + incr[i]
        midpoint += dz[i]
        temp[i] = (bht - tht) * midpoint / length + tht + 459.67

    p[0] = thp

    for i in range(1, ndiv + 1):
        temp_f_i = temp[i] - 459.67
        p_est = p[i - 1]

        for _iter in range(2):
            p_avg = max((p[i - 1] + p_est) / 2.0, 14.7)

            cgr_loc, qo_loc, ql_loc, lsg_loc = _condensate_dropout(
                cgr, qg_mmscfd, p_avg, pr, osg, qw_bwpd, wsg)
            ul = ql_loc * 5.615 / 86400.0 / area

            oil_vis_loc = _condensate_vis(pr, cgr_loc, gsg, api,
                                          temp_f_i, p_avg, oil_vis)

            zee = _z_factor(gsg, temp_f_i, p_avg)

            mflow_o = osg * 62.4 * qo_loc * 5.615
            mflow_w = wsg * 62.4 * qw_bwpd * 5.615
            mflow_g = 0.0765 * gsg * qg_mmscfd * 1e6

            water_visc = _water_viscosity(p_avg, temp_f_i)

            rs_est = _standing_rs(gsg, p_avg, temp_f_i, api)
            ift = _interfacial_tension(p_avg, temp_f_i, api, rs_est,
                                       mflow_o, mflow_w, mflow_g)

            ugas = (qg_mmscfd * 1e6 /
                    (p_avg * 35.3741 / (zee * temp[i])) /
                    86400.0 / area)
            um = ugas + ul

            nvl = 1.938 * ul * (62.4 * lsg_loc / ift) ** 0.25
            nvg = 1.938 * ugas * (62.4 * lsg_loc / ift) ** 0.25
            nd = 120.872 * tid / 12.0 * (62.4 * lsg_loc / ift) ** 0.5

            mul = 1.0
            if ql_loc > 0:
                mul = (qo_loc * oil_vis_loc + qw_bwpd * water_visc) / (qo_loc + qw_bwpd)

            nl = 0.15726 * mul * (1.0 / (62.4 * lsg_loc * ift ** 3)) ** 0.25
            if nl <= 0:
                nl = 1e-8

            x1 = _log10(nl) + 3.0
            cnl = 10.0 ** (-2.69851 + 0.1584095 * x1 - 0.5509976 * x1 * x1 +
                           0.5478492 * x1 ** 3 - 0.1219458 * x1 ** 4)

            nvg_safe = max(nvg, 1e-10)
            f1 = nvl * p_avg ** 0.1 * cnl / (nvg_safe ** 0.575 * 14.7 ** 0.1 * nd)

            lf1 = _log10(f1) + 6.0
            ylonsi = (-0.10306578 + 0.617774 * lf1 - 0.632946 * lf1 ** 2 +
                      0.29598 * lf1 ** 3 - 0.0401 * lf1 ** 4)
            ylonsi = max(ylonsi, 0.0)

            f2 = nvg * nl ** 0.38 / nd ** 2.14
            idex = 1.0 if f2 >= 0.012 else -1.0
            f2 = (1.0 - idex) / 2.0 * 0.012 + (1.0 + idex) / 2.0 * f2

            si = (0.9116257 - 4.821756 * f2 + 1232.25 * f2 * f2 -
                  22253.58 * f2 ** 3 + 116174.3 * f2 ** 4)
            if f2 <= 0.001:
                si = 1.0

            yl = _clamp(si * ylonsi, 0.0, 1.0)

            rho = 29.0 * gsg * p_avg / (temp[i] * zee * 10.732 * 62.37)

            mflow = lsg_loc * 62.4 * ql_loc * 5.615 + 0.0765 * gsg * qg_mmscfd * 1e6
            mass_frac_liq = (lsg_loc * 62.4 * ql_loc * 5.615) / mflow if mflow > 0 else 0
            min_l = mass_frac_liq * rho * 62.37 / (rho * 62.37 + lsg_loc * 62.37)
            if yl < min_l:
                yl = min_l

            # Orkiszewski bubble flow correction
            vm = ugas + ul
            vs = 0.8
            d_ft = tid / 12.0
            lb = 1.071 - 0.2218 * vm * vm / d_ft
            if lb < 0.13:
                lb = 0.13
            if vm > 0:
                b_ratio = ugas / vm
                if b_ratio < lb:
                    disc = (1.0 + vm / vs) ** 2 - 4.0 * ugas / vs
                    if disc >= 0:
                        yl = 1.0 - 0.5 * (1.0 + vm / vs - math.sqrt(disc))
                        yl = _clamp(yl, 0.0, 1.0)

            mug = _gas_viscosity(gsg, temp_f_i, p_avg)
            nre = 96778.0 * (mflow / 86400.0) / ((tid / 12.0) *
                  mul ** yl * mug ** (1.0 - yl))
            if nre < 2100:
                nre = 2100.0

            eond = rough / tid
            f = _serghides_fanning(nre, eond)

            rho_g = _MW_AIR * gsg * p_avg / (zee * 10.73 * temp[i])
            rho_avg = yl * lsg_loc * 62.4 + (1.0 - yl) * rho_g

            fric_term = f * mflow ** 2 / 7.413e10 / (tid / 12.0) ** 5 / rho_avg
            dpdz = (rho_avg * math.sin(theta) + (-fric_term if injection else fric_term)) / 144.0

            p_est = p[i - 1] + dpdz * incr[i]

        p[i] = p_est

    return p[ndiv]


def _hb_fbhp_oil(thp, api, gsg, tid, rough, length, tht, bht,
                  wsg, qt_stbpd, gor, wc, pb, rsb, sgsp,
                  rsb_scale=1.0, injection=False, theta=math.pi / 2.0):
    wc_adj = max(wc, 1e-9)
    if qt_stbpd < 1e-7:
        return _static_oil_column_pressure(
            thp, length, tht, bht, wc_adj, wsg, api, sgsp, pb, rsb, rsb_scale,
            theta=theta)

    qo = qt_stbpd * (1.0 - wc)
    qw = qt_stbpd * wc
    osg = 141.5 / (api + 131.5)
    sgsto = osg
    rsb_for_calc = rsb / rsb_scale

    area = math.pi * tid * tid / 4.0 / 144.0
    ndiv = _calc_segments(length)
    seg_len_ft = length / ndiv

    depth = [0.0] * (ndiv + 1)
    dz = [0.0] * (ndiv + 1)
    incr = [0.0] * (ndiv + 1)
    temp = [0.0] * (ndiv + 1)
    p_arr = [0.0] * (ndiv + 1)

    temp[0] = tht
    midpoint = 0.0

    for i in range(ndiv + 1):
        depth[i] = i * seg_len_ft
        if i == 0:
            incr[i] = seg_len_ft
        else:
            incr[i] = depth[i] - depth[i - 1]
            dz[i] = (incr[i - 1] + incr[i]) / 2.0
        if i == ndiv:
            incr[i] = 0.5 * seg_len_ft
            dz[i] = incr[i - 1] + incr[i]
        midpoint += dz[i]
        temp[i] = (bht - tht) * midpoint / length + tht + 459.67

    p_arr[0] = thp

    for i in range(1, ndiv + 1):
        temp_f_i = temp[i] - 459.67
        p_est = p_arr[i - 1]

        for _iter in range(2):
            p_avg = max((p_arr[i - 1] + p_est) / 2.0, 14.7)

            rs_local = _velarde_rs(sgsp, api, temp_f_i, pb, rsb_for_calc, p_avg) * rsb_scale
            free_gas = max(gor - rs_local, 0.0)
            qg_mmscfd = max(free_gas * qo / 1e6, 1e-9)

            oil_vis_seg = _oil_viscosity_full(sgsp, api, temp_f_i, rsb, pb, p_avg)
            rho_oil_local = _oil_density_mccain(rs_local, sgsp, sgsto,
                                                 min(p_avg, pb), temp_f_i)

            zee = _z_factor(gsg, temp_f_i, p_avg)
            mug = _gas_viscosity(gsg, temp_f_i, p_avg)

            mflow_o = osg * 62.4 * qo * 5.615
            mflow_w = wsg * 62.4 * qw * 5.615
            mflow_g = 0.0765 * gsg * qg_mmscfd * 1e6
            mflow = mflow_o + mflow_w + mflow_g

            if mflow < 1e-6:
                tavg_r = (tht + bht) / 2.0 + 459.67
                p_arr[i] = p_arr[i - 1] * math.exp(0.01875 * gsg * incr[i] / (zee * tavg_r))
                break

            ql = qo + qw
            lsg = (qo * osg + qw * wsg) / ql
            ul = ql * 5.615 / 86400.0 / area

            water_visc = _water_viscosity(p_avg, temp_f_i)
            ift = _interfacial_tension(p_avg, temp_f_i, api, rs_local,
                                       mflow_o, mflow_w, mflow_g)

            ugas = (qg_mmscfd * 1e6 /
                    (p_avg * 35.3741 / (zee * temp[i])) / 86400.0 / area)

            nvl = 1.938 * ul * (62.4 * lsg / ift) ** 0.25
            nvg = 1.938 * ugas * (62.4 * lsg / ift) ** 0.25
            nd = 120.872 * tid / 12.0 * (62.4 * lsg / ift) ** 0.5

            mul = 1.0
            if ql > 0:
                mul = (qo * oil_vis_seg + qw * water_visc) / ql

            nl = 0.15726 * mul * (1.0 / (62.4 * lsg * ift ** 3)) ** 0.25
            if nl <= 0:
                nl = 1e-8

            x1 = _log10(nl) + 3.0
            cnl = 10.0 ** (-2.69851 + 0.1584095 * x1 - 0.5509976 * x1 ** 2 +
                           0.5478492 * x1 ** 3 - 0.1219458 * x1 ** 4)

            nvg_safe = max(nvg, 1e-10)
            f1 = nvl * p_avg ** 0.1 * cnl / (nvg_safe ** 0.575 * 14.7 ** 0.1 * nd)

            lf1 = _log10(f1) + 6.0
            ylonsi = (-0.10306578 + 0.617774 * lf1 - 0.632946 * lf1 ** 2 +
                      0.29598 * lf1 ** 3 - 0.0401 * lf1 ** 4)
            ylonsi = max(ylonsi, 0.0)

            f2 = nvg * nl ** 0.38 / nd ** 2.14
            idex = 1.0 if f2 >= 0.012 else -1.0
            f2 = (1.0 - idex) / 2.0 * 0.012 + (1.0 + idex) / 2.0 * f2

            si = (0.9116257 - 4.821756 * f2 + 1232.25 * f2 ** 2 -
                  22253.58 * f2 ** 3 + 116174.3 * f2 ** 4)
            if f2 <= 0.001:
                si = 1.0

            yl = _clamp(si * ylonsi, 0.0, 1.0)

            rho_g = _MW_AIR * gsg * p_avg / (zee * 10.73 * temp[i])

            mass_frac_liq = (lsg * 62.4 * ql * 5.615) / mflow if mflow > 0 else 0
            min_l = mass_frac_liq * rho_g * 62.37 / (rho_g * 62.37 + lsg * 62.37)
            if yl < min_l:
                yl = min_l

            # Orkiszewski bubble flow correction
            vm = ugas + ul
            vs = 0.8
            d_ft = tid / 12.0
            lb_val = 1.071 - 0.2218 * vm * vm / d_ft
            if lb_val < 0.13:
                lb_val = 0.13
            if vm > 0:
                b_ratio = ugas / vm
                if b_ratio < lb_val:
                    disc = (1.0 + vm / vs) ** 2 - 4.0 * ugas / vs
                    if disc >= 0:
                        yl = 1.0 - 0.5 * (1.0 + vm / vs - math.sqrt(disc))
                        yl = _clamp(yl, 0.0, 1.0)

            nre = 96778.0 * (mflow / 86400.0) / ((tid / 12.0) *
                  mul ** yl * mug ** (1.0 - yl))
            if nre < 2100:
                nre = 2100.0

            eond = rough / tid
            f_fan = _serghides_fanning(nre, eond)

            rho_avg = yl * rho_oil_local + (1.0 - yl) * rho_g

            fric_term = f_fan * mflow ** 2 / 7.413e10 / (tid / 12.0) ** 5 / rho_avg
            dpdz = (rho_avg * math.sin(theta) + (-fric_term if injection else fric_term)) / 144.0

            p_est = p_arr[i - 1] + dpdz * incr[i]

        if p_arr[i] == 0:
            p_arr[i] = p_est

    return p_arr[ndiv]


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
    u_gm_base = 2.9 * max(buoy, 0.0) ** 0.25
    p_exp = _P_ATM_PA / p_sys
    incl_factor = (1.22 + 1.22 * math.sin(theta)) ** p_exp
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
    liq_turb = re_l >= 2100
    gas_turb = re_g >= 2100
    if liq_turb and gas_turb:
        chisholm_c = 20.0
    elif not liq_turb and gas_turb:
        chisholm_c = 12.0
    elif liq_turb and not gas_turb:
        chisholm_c = 10.0
    else:
        chisholm_c = 5.0
    phi2_l = 1.0 + chisholm_c / x_param + 1.0 / (x_param * x_param)
    return phi2_l * dpdz_l


def _wg_fbhp_gas(thp, api, gsg, tid, rough, length, tht, bht,
                  wsg, qg_mmscfd, cgr, qw_bwpd, oil_vis,
                  injection=False, pr=0.0, theta=math.pi / 2.0):
    osg = 141.5 / (api + 131.5)
    total_mass = (0.0765 * gsg * qg_mmscfd * 1e6 +
                  osg * 62.4 * cgr * qg_mmscfd * 5.615 +
                  wsg * 62.4 * qw_bwpd * 5.615)
    if total_mass < 1e-6 or qg_mmscfd < 0.05:
        return _static_gas_column_pressure(thp, length, tht, bht, gsg, theta=theta)

    diam_m = tid * _IN_TO_M
    rough_m = rough * _IN_TO_M
    area_m2 = math.pi * diam_m ** 2 / 4.0

    m_flow_g_kgs = 0.0765 * gsg * qg_mmscfd * 1e6 * 0.453592 / 86400.0
    m_flow_w_kgs = wsg * 62.4 * qw_bwpd * 5.615 * 0.453592 / 86400.0

    ndiv = _calc_segments(length)
    seg_len_ft = length / ndiv

    p_psia = thp

    for i in range(1, ndiv + 1):
        frac = (i - 0.5) / ndiv
        temp_f_i = tht + (bht - tht) * frac
        p_est = p_psia

        for _it in range(2):
            p_avg = (p_psia + p_est) / 2.0
            if p_avg < 14.7:
                continue

            cgr_loc, qo_loc, ql_loc, lsg_loc = _condensate_dropout(
                cgr, qg_mmscfd, p_avg, pr, osg, qw_bwpd, wsg)

            m_flow_o_kgs = osg * 62.4 * qo_loc * 5.615 * 0.453592 / 86400.0
            m_flow_l_kgs = m_flow_o_kgs + m_flow_w_kgs
            rho_l_kgm3 = lsg_loc * 62.4 * _LBFT3_TO_KGM3
            u_sl = m_flow_l_kgs / (rho_l_kgm3 * area_m2)

            oil_vis_loc = _condensate_vis(pr, cgr_loc, gsg, api,
                                           temp_f_i, p_avg, oil_vis)

            zee = _z_factor(gsg, temp_f_i, p_avg)
            mu_g_cp = _gas_viscosity(gsg, temp_f_i, p_avg)

            temp_r = temp_f_i + 459.67
            rho_g_lbft3 = _MW_AIR * gsg * p_avg / (zee * 10.732 * temp_r)
            rho_g_kgm3 = rho_g_lbft3 * _LBFT3_TO_KGM3

            qg_actual_m3s = m_flow_g_kgs / rho_g_kgm3
            u_sg = qg_actual_m3s / area_m2

            water_visc = _water_viscosity(p_avg, temp_f_i)
            if ql_loc > 0:
                mu_l = (qo_loc * oil_vis_loc + qw_bwpd * water_visc) / ql_loc * _CP_TO_PAS
            else:
                mu_l = oil_vis_loc * _CP_TO_PAS

            rs_est = _standing_rs(gsg, p_avg, temp_f_i, api)
            ift = _interfacial_tension(p_avg, temp_f_i, api, rs_est,
                                       m_flow_o_kgs, m_flow_w_kgs, m_flow_g_kgs)
            sigma_nm = ift * _DYNECM_TO_NM
            p_sys_pa = p_avg * _PSI_TO_PA

            alpha_g = _wg_void_fraction(u_sg, u_sl, rho_g_kgm3, rho_l_kgm3,
                                         sigma_nm, diam_m, theta, p_sys_pa)
            liq_holdup = 1.0 - alpha_g

            rho_mix = alpha_g * rho_g_kgm3 + liq_holdup * rho_l_kgm3
            dpdz_hydro = rho_mix * _G_SI * math.sin(theta)

            mu_g_pas = mu_g_cp * _CP_TO_PAS
            dpdz_fric = _wg_friction_gradient_lm(
                m_flow_g_kgs, m_flow_l_kgs, rho_g_kgm3, rho_l_kgm3,
                mu_g_pas, mu_l, diam_m, rough_m)

            dpdz_total_pam = dpdz_hydro + (-dpdz_fric if injection else dpdz_fric)
            dpdz_psift = dpdz_total_pam / (_PSI_TO_PA / _FT_TO_M)

            p_est = p_psia + dpdz_psift * seg_len_ft

        p_psia = p_est

    return p_psia


def _wg_fbhp_oil(thp, api, gsg, tid, rough, length, tht, bht,
                  wsg, qt_stbpd, gor, wc, pb, rsb, sgsp,
                  rsb_scale=1.0, injection=False, theta=math.pi / 2.0):
    wc_adj = max(wc, 1e-9)
    if qt_stbpd < 1e-7:
        return _static_oil_column_pressure(
            thp, length, tht, bht, wc_adj, wsg, api, sgsp, pb, rsb, rsb_scale,
            theta=theta)

    qo = qt_stbpd * (1.0 - wc)
    qw = qt_stbpd * wc
    osg = 141.5 / (api + 131.5)
    rsb_for_calc = rsb / rsb_scale

    diam_m = tid * _IN_TO_M
    rough_m = rough * _IN_TO_M
    area_m2 = math.pi * diam_m ** 2 / 4.0

    ndiv = _calc_segments(length)
    seg_len_ft = length / ndiv

    p_psia = thp

    for i in range(1, ndiv + 1):
        frac = (i - 0.5) / ndiv
        temp_f_i = tht + (bht - tht) * frac
        p_est = p_psia

        for _it in range(2):
            p_avg = max((p_psia + p_est) / 2.0, 14.7)

            rs_local = _velarde_rs(sgsp, api, temp_f_i, pb, rsb_for_calc, p_avg) * rsb_scale
            free_gas = max(gor - rs_local, 0.0)
            qg_mmscfd = max(free_gas * qo / 1e6, 1e-9)

            oil_vis_seg = _oil_viscosity_full(sgsp, api, temp_f_i, rsb, pb, p_avg)
            rho_oil_lbft3 = _oil_density_mccain(rs_local, sgsp, osg,
                                                  min(p_avg, pb), temp_f_i)
            rho_oil_kgm3 = rho_oil_lbft3 * _LBFT3_TO_KGM3

            zee = _z_factor(gsg, temp_f_i, p_avg)
            mu_g_cp = _gas_viscosity(gsg, temp_f_i, p_avg)

            temp_r = temp_f_i + 459.67
            rho_g_lbft3 = _MW_AIR * gsg * p_avg / (zee * 10.732 * temp_r)
            rho_g_kgm3 = rho_g_lbft3 * _LBFT3_TO_KGM3

            m_flow_o_kgs = osg * 62.4 * qo * 5.615 * 0.453592 / 86400.0
            m_flow_w_kgs = wsg * 62.4 * qw * 5.615 * 0.453592 / 86400.0
            m_flow_g_kgs = 0.0765 * gsg * qg_mmscfd * 1e6 * 0.453592 / 86400.0
            m_flow_l_kgs = m_flow_o_kgs + m_flow_w_kgs
            m_flow_total = m_flow_g_kgs + m_flow_l_kgs

            if m_flow_total < 1e-15:
                tavg_r = (tht + bht) / 2.0 + 459.67
                p_psia *= math.exp(0.01875 * gsg * seg_len_ft / (zee * tavg_r))
                break

            ql = qo + qw
            rho_w_kgm3 = wsg * 62.4 * _LBFT3_TO_KGM3
            rho_l_kgm3 = ((qo * rho_oil_kgm3 + qw * rho_w_kgm3) / ql
                           if ql > 0 else rho_oil_kgm3)

            u_sg = m_flow_g_kgs / (rho_g_kgm3 * area_m2)
            u_sl = m_flow_l_kgs / (rho_l_kgm3 * area_m2)

            water_visc = _water_viscosity(p_avg, temp_f_i)
            mu_l_pas = ((qo * oil_vis_seg + qw * water_visc) / ql * _CP_TO_PAS
                        if ql > 0 else oil_vis_seg * _CP_TO_PAS)

            ift = _interfacial_tension(p_avg, temp_f_i, api, rs_local,
                                       m_flow_o_kgs, m_flow_w_kgs, m_flow_g_kgs)
            sigma_nm = ift * _DYNECM_TO_NM
            p_sys_pa = p_avg * _PSI_TO_PA

            alpha_g = _wg_void_fraction(u_sg, u_sl, rho_g_kgm3, rho_l_kgm3,
                                         sigma_nm, diam_m, theta, p_sys_pa)
            liq_holdup = 1.0 - alpha_g

            rho_mix = alpha_g * rho_g_kgm3 + liq_holdup * rho_l_kgm3
            dpdz_hydro = rho_mix * _G_SI * math.sin(theta)

            mu_g_pas = mu_g_cp * _CP_TO_PAS
            dpdz_fric = _wg_friction_gradient_lm(
                m_flow_g_kgs, m_flow_l_kgs, rho_g_kgm3, rho_l_kgm3,
                mu_g_pas, mu_l_pas, diam_m, rough_m)

            dpdz_total_pam = dpdz_hydro + (-dpdz_fric if injection else dpdz_fric)
            dpdz_psift = dpdz_total_pam / (_PSI_TO_PA / _FT_TO_M)

            p_est = p_psia + dpdz_psift * seg_len_ft

        p_psia = p_est

    return p_psia


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

    sigma_lbf_ft = sigma * 6.852e-5
    if sigma_lbf_ft <= 0:
        return lambda_l

    rho_ns = rho_l * lambda_l + rho_g * (1.0 - lambda_l)
    if rho_ns <= 0:
        return lambda_l

    rv = lambda_l
    drho = max(rho_l - rho_g, 0.1)
    n1 = v_m ** 2 * rho_ns / (_G_FT * diam * drho)
    n2 = _G_FT * diam ** 2 * drho / sigma_lbf_ft

    a = 0.0814 * (1.0 - 0.0554 * math.log(1.0 + 730.0 * rv / (1.0 + rv)))
    b = 0.0554

    if n1 <= 1e-15:
        return lambda_l

    argument = (a + b * n2) / n1
    if argument > 0:
        f_e = -0.0554 * math.log(argument)
    else:
        f_e = 0.0

    hl = 1.0 - (1.0 - lambda_l) * math.exp(f_e)
    return _clamp(hl, lambda_l, 1.0)


def _gray_effective_roughness(rough_dry, sigma, rho_ns, v_m, lambda_l):
    sigma_lbf_ft = sigma * 6.852e-5
    if v_m < 1e-10 or rho_ns <= 0 or sigma_lbf_ft <= 0:
        return rough_dry
    ke = rho_ns * v_m * v_m
    r1 = 28.5 * sigma_lbf_ft / max(ke, 1e-10)
    r2 = ke * lambda_l / sigma_lbf_ft
    film = r1 * (1.0 - math.exp(-r2))
    return rough_dry + max(film, 0.0)


def _gray_fbhp_gas(thp, api, gsg, tid, rough, length, tht, bht,
                    wsg, qg_mmscfd, cgr, qw_bwpd, oil_vis,
                    injection=False, pr=0.0, theta=math.pi / 2.0):
    osg = 141.5 / (api + 131.5)
    total_mass = (0.0765 * gsg * qg_mmscfd * 1e6 +
                  osg * 62.4 * cgr * qg_mmscfd * 5.615 +
                  wsg * 62.4 * qw_bwpd * 5.615)
    if total_mass < 1e-6 or qg_mmscfd < 0.05:
        return _static_gas_column_pressure(thp, length, tht, bht, gsg, theta=theta)

    diam_ft = tid / 12.0
    rough_ft = rough / 12.0
    area = math.pi * diam_ft ** 2 / 4.0

    ndiv = _calc_segments(length)
    seg_len = length / ndiv

    m_flow_g_lbs = 0.0765 * gsg * qg_mmscfd * 1e6 / 86400.0
    m_flow_w_lbs = wsg * 62.4 * qw_bwpd * 5.615 / 86400.0

    p_psia = thp

    for i in range(1, ndiv + 1):
        frac = (i - 0.5) / ndiv
        temp_f_i = tht + (bht - tht) * frac
        p_est = p_psia

        for _it in range(2):
            p_avg = (p_psia + p_est) / 2.0
            if p_avg < 14.7:
                continue

            cgr_loc, qo_loc, ql_loc, lsg_loc = _condensate_dropout(
                cgr, qg_mmscfd, p_avg, pr, osg, qw_bwpd, wsg)

            m_flow_o_lbs = osg * 62.4 * qo_loc * 5.615 / 86400.0
            m_flow_l_lbs = m_flow_o_lbs + m_flow_w_lbs
            m_flow_total_lbs = m_flow_g_lbs + m_flow_l_lbs

            oil_vis_loc = _condensate_vis(pr, cgr_loc, gsg, api,
                                           temp_f_i, p_avg, oil_vis)

            zee = _z_factor(gsg, temp_f_i, p_avg)
            mu_g_cp = _gas_viscosity(gsg, temp_f_i, p_avg)

            temp_r = temp_f_i + 459.67
            rho_g = _MW_AIR * gsg * p_avg / (zee * 10.732 * temp_r)
            rho_l = lsg_loc * 62.4

            v_sg = (m_flow_g_lbs / max(rho_g, 1e-10)) / area
            v_sl = (m_flow_l_lbs / max(rho_l, 1e-10)) / area
            v_m = v_sg + v_sl
            lambda_l = v_sl / v_m if v_m > 1e-10 else 0.0
            rho_ns = rho_l * lambda_l + rho_g * (1.0 - lambda_l)

            rs_est = _standing_rs(gsg, p_avg, temp_f_i, api)
            sigma = _interfacial_tension(p_avg, temp_f_i, api, rs_est,
                                         m_flow_o_lbs * 0.453592,
                                         m_flow_w_lbs * 0.453592,
                                         m_flow_g_lbs * 0.453592)

            hl = _gray_liquid_holdup(v_m, rho_l, rho_g, sigma, diam_ft,
                                      lambda_l, p_avg)
            rho_s = rho_l * hl + rho_g * (1.0 - hl)

            eps_eff = _gray_effective_roughness(rough_ft, sigma, rho_ns,
                                                v_m, lambda_l)
            eps_d = eps_eff / diam_ft

            water_visc = _water_viscosity(p_avg, temp_f_i)
            mu_l_cp = ((qo_loc * oil_vis_loc + qw_bwpd * water_visc) / ql_loc
                       if ql_loc > 0 else oil_vis_loc)
            mu_ns = mu_l_cp * lambda_l + mu_g_cp * (1.0 - lambda_l)
            mu_ns_lbfts = mu_ns * 6.7197e-4
            n_re = rho_ns * v_m * diam_ft / mu_ns_lbfts if mu_ns_lbfts > 0 else 0.0

            f_fanning = _serghides_fanning(n_re, eps_d)
            f_moody = 4.0 * f_fanning

            gm = m_flow_total_lbs / area

            dpdz_hydro = rho_s * math.sin(theta) / 144.0
            dpdz_fric = (f_moody * gm ** 2 / (2.0 * _GC * diam_ft * rho_ns * 144.0)
                         if rho_ns > 0 else 0.0)
            ek = gm * v_sg / (_GC * p_avg * 144.0)
            denom_val = max(1.0 - ek, 0.1)

            dpdz_total = (dpdz_hydro + (-dpdz_fric if injection else dpdz_fric)) / denom_val
            p_est = p_psia + dpdz_total * seg_len

        p_psia = p_est

    return p_psia


def _gray_fbhp_oil(thp, api, gsg, tid, rough, length, tht, bht,
                    wsg, qt_stbpd, gor, wc, pb, rsb, sgsp,
                    rsb_scale=1.0, injection=False, theta=math.pi / 2.0):
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
        temp_f_i = tht + (bht - tht) * frac
        p_est = p_psia

        for _it in range(2):
            p_avg = max((p_psia + p_est) / 2.0, 14.7)

            rs_local = _velarde_rs(sgsp, api, temp_f_i, pb, rsb_for_calc, p_avg) * rsb_scale
            free_gas = max(gor - rs_local, 0.0)
            qg_mmscfd = max(free_gas * qo / 1e6, 1e-9)

            oil_vis_seg = _oil_viscosity_full(sgsp, api, temp_f_i, rsb, pb, p_avg)
            rho_oil_local = _oil_density_mccain(rs_local, sgsp, osg,
                                                 min(p_avg, pb), temp_f_i)

            zee = _z_factor(gsg, temp_f_i, p_avg)
            mu_g_cp = _gas_viscosity(gsg, temp_f_i, p_avg)

            temp_r = temp_f_i + 459.67
            rho_g = _MW_AIR * gsg * p_avg / (zee * 10.732 * temp_r)

            m_flow_o_lbs = osg * 62.4 * qo * 5.615 / 86400.0
            m_flow_w_lbs = wsg * 62.4 * qw * 5.615 / 86400.0
            m_flow_g_lbs = 0.0765 * gsg * qg_mmscfd * 1e6 / 86400.0
            m_flow_l_lbs = m_flow_o_lbs + m_flow_w_lbs
            m_flow_total_lbs = m_flow_g_lbs + m_flow_l_lbs

            if m_flow_total_lbs < 1e-15:
                tavg_r = (tht + bht) / 2.0 + 459.67
                p_psia *= math.exp(0.01875 * gsg * seg_len / (zee * tavg_r))
                break

            ql = qo + qw
            rho_w = wsg * 62.4
            rho_l = (qo * rho_oil_local + qw * rho_w) / ql if ql > 0 else rho_oil_local

            v_sg = (m_flow_g_lbs / max(rho_g, 1e-10)) / area
            v_sl = (m_flow_l_lbs / max(rho_l, 1e-10)) / area
            v_m = v_sg + v_sl
            lambda_l = v_sl / v_m if v_m > 1e-10 else 0.0
            rho_ns = rho_l * lambda_l + rho_g * (1.0 - lambda_l)

            water_visc = _water_viscosity(p_avg, temp_f_i)
            mu_l_cp = ((qo * oil_vis_seg + qw * water_visc) / ql
                       if ql > 0 else oil_vis_seg)

            ift = _interfacial_tension(p_avg, temp_f_i, api, rs_local,
                                       m_flow_o_lbs * 0.453592,
                                       m_flow_w_lbs * 0.453592,
                                       m_flow_g_lbs * 0.453592)

            hl = _gray_liquid_holdup(v_m, rho_l, rho_g, ift, diam_ft,
                                      lambda_l, p_avg)
            rho_s = rho_l * hl + rho_g * (1.0 - hl)

            eps_eff = _gray_effective_roughness(rough_ft, ift, rho_ns,
                                                v_m, lambda_l)
            eps_d = eps_eff / diam_ft

            mu_ns = mu_l_cp * lambda_l + mu_g_cp * (1.0 - lambda_l)
            mu_ns_lbfts = mu_ns * 6.7197e-4
            n_re = rho_ns * v_m * diam_ft / mu_ns_lbfts if mu_ns_lbfts > 0 else 0.0

            f_moody = 4.0 * _serghides_fanning(n_re, eps_d)

            gm = m_flow_total_lbs / area
            dpdz_hydro = rho_s * math.sin(theta) / 144.0
            dpdz_fric = (f_moody * gm ** 2 / (2.0 * _GC * diam_ft * rho_ns * 144.0)
                         if rho_ns > 0 else 0.0)
            ek = gm * v_sg / (_GC * p_avg * 144.0)
            denom_val = max(1.0 - ek, 0.1)

            dpdz_total = (dpdz_hydro + (-dpdz_fric if injection else dpdz_fric)) / denom_val
            p_est = p_psia + dpdz_total * seg_len

        p_psia = p_est

    return p_psia


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
    l1 = 316.0 * lambda_l ** 0.302
    l2 = 0.0009252 * lambda_l ** (-2.4684)
    l3 = 0.10 * lambda_l ** (-1.4516)
    l4 = 0.5 * lambda_l ** (-6.738)
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
        hl0 = 0.98 * lambda_l ** 0.4846 / froude ** 0.0868
    elif pattern == _BB_INTERMITTENT:
        hl0 = 0.845 * lambda_l ** 0.5351 / froude ** 0.0173
    elif pattern == _BB_DISTRIBUTED:
        hl0 = 1.065 * lambda_l ** 0.5824 / froude ** 0.0609
    else:
        hl0 = lambda_l
    return max(hl0, lambda_l)


def _bb_inclination_correction(hl0, lambda_l, n_lv, froude, pattern,
                                theta=math.pi / 2.0):
    if lambda_l <= 0 or lambda_l >= 1.0:
        return hl0
    if pattern == _BB_SEGREGATED:
        e_p, f_p, g_p, h_p = 0.011, -3.7680, 3.5390, -1.6140
    elif pattern == _BB_INTERMITTENT:
        e_p, f_p, g_p, h_p = 2.960, 0.3050, -0.4473, 0.0978
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
        denom = (-0.0523 + 3.182 * ln_y - 0.8725 * ln_y ** 2 +
                 0.01853 * ln_y ** 4)
        s = ln_y / denom if abs(denom) >= 1e-6 else 0.0
    s = _clamp(s, -5.0, 5.0)
    return f_ns * math.exp(s)


def _bb_core_gas(thp, api, gsg, tid, rough, length, tht, bht,
                 wsg, qg_mmscfd, cgr, qw_bwpd, oil_vis,
                 injection, pr, theta=math.pi / 2.0):
    """Beggs & Brill core for gas wells."""
    osg = 141.5 / (api + 131.5)
    total_mass = (0.0765 * gsg * qg_mmscfd * 1e6 +
                  osg * 62.4 * cgr * qg_mmscfd * 5.615 +
                  wsg * 62.4 * qw_bwpd * 5.615)
    if total_mass < 1e-6 or qg_mmscfd < 0.05:
        return _static_gas_column_pressure(thp, length, tht, bht, gsg, theta=theta)

    diam_ft = tid / 12.0
    rough_ft = rough / 12.0
    area = math.pi * diam_ft ** 2 / 4.0
    eps_d = rough_ft / diam_ft

    ndiv = _calc_segments(length)
    seg_len = length / ndiv

    m_flow_g_lbs = 0.0765 * gsg * qg_mmscfd * 1e6 / 86400.0
    m_flow_w_lbs = wsg * 62.4 * qw_bwpd * 5.615 / 86400.0

    p_psia = thp

    for i in range(1, ndiv + 1):
        frac = (i - 0.5) / ndiv
        temp_f_i = tht + (bht - tht) * frac
        p_est = p_psia

        for _it in range(2):
            p_avg = (p_psia + p_est) / 2.0
            if p_avg < 14.7:
                continue

            cgr_loc, qo_loc, ql_loc, lsg_loc = _condensate_dropout(
                cgr, qg_mmscfd, p_avg, pr, osg, qw_bwpd, wsg)

            m_flow_o_lbs = osg * 62.4 * qo_loc * 5.615 / 86400.0
            m_flow_l_lbs = m_flow_o_lbs + m_flow_w_lbs

            oil_vis_loc = _condensate_vis(pr, cgr_loc, gsg, api,
                                           temp_f_i, p_avg, oil_vis)

            zee = _z_factor(gsg, temp_f_i, p_avg)
            mu_g_cp = _gas_viscosity(gsg, temp_f_i, p_avg)

            temp_r = temp_f_i + 459.67
            rho_g = _MW_AIR * gsg * p_avg / (zee * 10.732 * temp_r)
            rho_l = lsg_loc * 62.4

            v_sg = (m_flow_g_lbs / max(rho_g, 1e-10)) / area
            v_sl = (m_flow_l_lbs / max(rho_l, 1e-10)) / area
            v_m = v_sg + v_sl
            lambda_l = v_sl / v_m if v_m > 1e-10 else 0.0
            rho_ns = rho_l * lambda_l + rho_g * (1.0 - lambda_l)

            froude = v_m ** 2 / (_G_FT * diam_ft)
            pattern, trans_a = _bb_flow_pattern(froude, lambda_l)

            if pattern == _BB_TRANSITION:
                hl0_seg = _bb_horizontal_holdup(lambda_l, froude, _BB_SEGREGATED)
                hl0_int = _bb_horizontal_holdup(lambda_l, froude, _BB_INTERMITTENT)
                hl0 = trans_a * hl0_seg + (1.0 - trans_a) * hl0_int
            else:
                hl0 = _bb_horizontal_holdup(lambda_l, froude, pattern)

            if pattern in (_BB_SEGREGATED, _BB_INTERMITTENT, _BB_TRANSITION):
                hl0 *= 0.924
                hl0 = max(hl0, lambda_l)

            rs_est = _standing_rs(gsg, p_avg, temp_f_i, api)
            sigma = _interfacial_tension(p_avg, temp_f_i, api, rs_est,
                                         m_flow_o_lbs * 0.453592,
                                         m_flow_w_lbs * 0.453592,
                                         m_flow_g_lbs * 0.453592)
            sigma_lbf_ft = sigma * 6.852e-5
            n_lv = (1.938 * v_sl * (rho_l / sigma_lbf_ft) ** 0.25
                    if sigma_lbf_ft > 0 else 0.0)

            if pattern == _BB_TRANSITION:
                hl_seg = _bb_inclination_correction(
                    _bb_horizontal_holdup(lambda_l, froude, _BB_SEGREGATED) * 0.924,
                    lambda_l, n_lv, froude, _BB_SEGREGATED, theta=theta)
                hl_int = _bb_inclination_correction(
                    _bb_horizontal_holdup(lambda_l, froude, _BB_INTERMITTENT) * 0.924,
                    lambda_l, n_lv, froude, _BB_INTERMITTENT, theta=theta)
                hl_theta = trans_a * hl_seg + (1.0 - trans_a) * hl_int
            else:
                hl_theta = _bb_inclination_correction(
                    hl0, lambda_l, n_lv, froude, pattern, theta=theta)

            hl_theta = _clamp(hl_theta, lambda_l, 1.0)
            rho_s = rho_l * hl_theta + rho_g * (1.0 - hl_theta)

            water_visc = _water_viscosity(p_avg, temp_f_i)
            mu_l_cp = ((qo_loc * oil_vis_loc + qw_bwpd * water_visc) / ql_loc
                       if ql_loc > 0 else oil_vis_loc)
            mu_ns = mu_l_cp * lambda_l + mu_g_cp * (1.0 - lambda_l)
            mu_ns_lbfts = mu_ns * 6.7197e-4
            n_re = rho_ns * v_m * diam_ft / mu_ns_lbfts if mu_ns_lbfts > 0 else 0.0

            f_ns = _serghides_fanning(n_re, eps_d)
            f_tp = _bb_two_phase_friction(f_ns, lambda_l, hl_theta)

            dpdz_hydro = rho_s * math.sin(theta) / 144.0
            dpdz_fric = (4.0 * f_tp * rho_ns * v_m ** 2 /
                         (2.0 * _GC * diam_ft * 144.0) if rho_ns > 0 else 0.0)
            ek = rho_s * v_m * v_sg / (_GC * p_avg * 144.0)
            denom_val = max(1.0 - ek, 0.1)

            dpdz_total = (dpdz_hydro + (-dpdz_fric if injection else dpdz_fric)) / denom_val
            p_est = p_psia + dpdz_total * seg_len

        p_psia = p_est

    return p_psia


def _bb_core_oil(thp, api, gsg, tid, rough, length, tht, bht,
                 wsg, qt_stbpd, gor, wc, pb, rsb, sgsp,
                 rsb_scale, injection, theta=math.pi / 2.0):
    """Beggs & Brill core for oil wells."""
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
    eps_d = rough_ft / diam_ft

    ndiv = _calc_segments(length)
    seg_len = length / ndiv

    p_psia = thp

    for i in range(1, ndiv + 1):
        frac = (i - 0.5) / ndiv
        temp_f_i = tht + (bht - tht) * frac
        p_est = p_psia

        for _it in range(2):
            p_avg = max((p_psia + p_est) / 2.0, 14.7)

            rs_local = _velarde_rs(sgsp, api, temp_f_i, pb, rsb_for_calc, p_avg) * rsb_scale
            free_gas = max(gor - rs_local, 0.0)
            qg_mmscfd = max(free_gas * qo / 1e6, 1e-9)

            oil_vis_seg = _oil_viscosity_full(sgsp, api, temp_f_i, rsb, pb, p_avg)
            rho_oil_local = _oil_density_mccain(rs_local, sgsp, osg,
                                                 min(p_avg, pb), temp_f_i)

            zee = _z_factor(gsg, temp_f_i, p_avg)
            mu_g_cp = _gas_viscosity(gsg, temp_f_i, p_avg)

            temp_r = temp_f_i + 459.67
            rho_g = _MW_AIR * gsg * p_avg / (zee * 10.732 * temp_r)

            m_flow_o_lbs = osg * 62.4 * qo * 5.615 / 86400.0
            m_flow_w_lbs = wsg * 62.4 * qw * 5.615 / 86400.0
            m_flow_g_lbs = 0.0765 * gsg * qg_mmscfd * 1e6 / 86400.0
            m_flow_l_lbs = m_flow_o_lbs + m_flow_w_lbs
            m_flow_total_lbs = m_flow_g_lbs + m_flow_l_lbs

            if m_flow_total_lbs < 1e-15:
                tavg_r = (tht + bht) / 2.0 + 459.67
                p_psia *= math.exp(0.01875 * gsg * seg_len / (zee * tavg_r))
                break

            ql = qo + qw
            rho_w = wsg * 62.4
            rho_l = (qo * rho_oil_local + qw * rho_w) / ql if ql > 0 else rho_oil_local

            v_sg = (m_flow_g_lbs / max(rho_g, 1e-10)) / area
            v_sl = (m_flow_l_lbs / max(rho_l, 1e-10)) / area
            v_m = v_sg + v_sl
            lambda_l = v_sl / v_m if v_m > 1e-10 else 0.0
            rho_ns = rho_l * lambda_l + rho_g * (1.0 - lambda_l)

            froude = v_m ** 2 / (_G_FT * diam_ft)
            pattern, trans_a = _bb_flow_pattern(froude, lambda_l)

            if pattern == _BB_TRANSITION:
                hl0_seg = _bb_horizontal_holdup(lambda_l, froude, _BB_SEGREGATED)
                hl0_int = _bb_horizontal_holdup(lambda_l, froude, _BB_INTERMITTENT)
                hl0 = trans_a * hl0_seg + (1.0 - trans_a) * hl0_int
            else:
                hl0 = _bb_horizontal_holdup(lambda_l, froude, pattern)

            if pattern in (_BB_SEGREGATED, _BB_INTERMITTENT, _BB_TRANSITION):
                hl0 *= 0.924
                hl0 = max(hl0, lambda_l)

            ift = _interfacial_tension(p_avg, temp_f_i, api, rs_local,
                                       m_flow_o_lbs * 0.453592,
                                       m_flow_w_lbs * 0.453592,
                                       m_flow_g_lbs * 0.453592)
            sigma_lbf_ft = ift * 6.852e-5
            n_lv = (1.938 * v_sl * (rho_l / sigma_lbf_ft) ** 0.25
                    if sigma_lbf_ft > 0 else 0.0)

            if pattern == _BB_TRANSITION:
                hl_seg = _bb_inclination_correction(
                    _bb_horizontal_holdup(lambda_l, froude, _BB_SEGREGATED) * 0.924,
                    lambda_l, n_lv, froude, _BB_SEGREGATED, theta=theta)
                hl_int = _bb_inclination_correction(
                    _bb_horizontal_holdup(lambda_l, froude, _BB_INTERMITTENT) * 0.924,
                    lambda_l, n_lv, froude, _BB_INTERMITTENT, theta=theta)
                hl_theta = trans_a * hl_seg + (1.0 - trans_a) * hl_int
            else:
                hl_theta = _bb_inclination_correction(
                    hl0, lambda_l, n_lv, froude, pattern, theta=theta)

            hl_theta = _clamp(hl_theta, lambda_l, 1.0)
            rho_s = rho_l * hl_theta + rho_g * (1.0 - hl_theta)

            water_visc = _water_viscosity(p_avg, temp_f_i)
            mu_l_cp = ((qo * oil_vis_seg + qw * water_visc) / ql
                       if ql > 0 else oil_vis_seg)
            mu_ns = mu_l_cp * lambda_l + mu_g_cp * (1.0 - lambda_l)
            mu_ns_lbfts = mu_ns * 6.7197e-4
            n_re = rho_ns * v_m * diam_ft / mu_ns_lbfts if mu_ns_lbfts > 0 else 0.0

            f_ns = _serghides_fanning(n_re, eps_d)
            f_tp = _bb_two_phase_friction(f_ns, lambda_l, hl_theta)

            dpdz_hydro = rho_s * math.sin(theta) / 144.0
            dpdz_fric = (4.0 * f_tp * rho_ns * v_m ** 2 /
                         (2.0 * _GC * diam_ft * 144.0) if rho_ns > 0 else 0.0)
            ek = rho_s * v_m * v_sg / (_GC * p_avg * 144.0)
            denom_val = max(1.0 - ek, 0.1)

            dpdz_total = (dpdz_hydro + (-dpdz_fric if injection else dpdz_fric)) / denom_val
            p_est = p_psia + dpdz_total * seg_len

        p_psia = p_est

    return p_psia


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

def fbhp(thp, completion, vlpmethod='HB', well_type='gas',
         gas_pvt=None, oil_pvt=None,
         qg_mmscfd=0, cgr=0, qw_bwpd=0, oil_vis=1.0, api=45, pr=0,
         qt_stbpd=0, gor=0, wc=0,
         wsg=1.07, injection=False,
         gsg=0.65, pb=0, rsb=0, sgsp=0.65):
    """ Returns flowing bottom hole pressure (psia) using specified VLP correlation.

        thp: Tubing head pressure (psia)
        completion: Completion object describing the wellbore
        vlpmethod: VLP method - 'HB' (Hagedorn-Brown), 'WG' (Woldesemayat-Ghajar), 'GRAY', or 'BB' (Beggs & Brill)
        well_type: 'gas' or 'oil'

        Gas well parameters:
            qg_mmscfd: Gas rate (MMscf/d)
            cgr: Condensate-gas ratio (STB/MMscf)
            qw_bwpd: Water rate (STB/d)
            oil_vis: Oil (condensate) viscosity (cP). Defaults to 1.0
            api: Condensate API gravity. Defaults to 45
            pr: Reservoir pressure (psia) - for condensate dropout. 0 disables

        Oil well parameters:
            qt_stbpd: Total liquid rate (STB/d)
            gor: Gas-oil ratio (scf/STB)
            wc: Water cut (fraction 0-1)

        Common parameters:
            wsg: Water specific gravity. Defaults to 1.07
            injection: True for injection wells. Defaults to False
            gsg: Gas specific gravity (relative to air). Defaults to 0.65
            pb: Bubble point pressure (psia). Required for oil wells
            rsb: Solution GOR at Pb (scf/STB). Required for oil wells
            sgsp: Separator gas specific gravity. Defaults to 0.65

        gas_pvt: GasPVT object (unused by VLP methods directly, reserved for future use)
        oil_pvt: OilPVT object. If provided for oil wells, extracts api, sgsp, pb, rsb from it
    """
    vlpmethod = validate_methods(["vlpmethod"], [vlpmethod])

    # Extract oil PVT parameters if provided
    if oil_pvt is not None and well_type == 'oil':
        api = oil_pvt.api
        sgsp = oil_pvt.sg_sp
        pb = oil_pvt.pb
        rsb = oil_pvt.rsb
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
                rsb_scale=1.0, injection=injection, theta=theta)

    # Loop over wellbore segments
    total_md = completion.total_md
    md_traversed = 0.0
    p_current = thp

    for seg in completion.segments:
        frac_start = md_traversed / total_md if total_md > 0 else 0.0
        frac_end = (md_traversed + seg.md) / total_md if total_md > 0 else 1.0
        tht_seg = completion.tht + (completion.bht - completion.tht) * frac_start
        bht_seg = completion.tht + (completion.bht - completion.tht) * frac_end

        p_current = _run_section(p_current, seg.id, seg.roughness, seg.md,
                                 tht_seg, bht_seg, seg.theta)
        md_traversed += seg.md

    return p_current


# ============================================================================
#  Public API: outflow_curve
# ============================================================================

def outflow_curve(thp, completion, vlpmethod='HB', well_type='gas',
                  gas_pvt=None, oil_pvt=None,
                  rates=None, n_rates=20, max_rate=None,
                  cgr=0, qw_bwpd=0, oil_vis=1.0, api=45, pr=0,
                  gor=0, wc=0,
                  wsg=1.07, injection=False,
                  gsg=0.65, pb=0, rsb=0, sgsp=0.65):
    """ Returns VLP outflow curve as dict {'rates': [...], 'bhp': [...]}.

        thp: Tubing head pressure (psia)
        completion: Completion object
        vlpmethod: VLP method string
        well_type: 'gas' or 'oil'
        rates: List of rates to evaluate (MMscf/d for gas, STB/d for oil). If None, auto-generated
        n_rates: Number of rate points if rates is None
        max_rate: Maximum rate for auto-generation. If None, defaults to 50 MMscf/d (gas) or 10000 STB/d (oil)
        Other parameters: Same as fbhp()
    """
    if oil_pvt is not None and well_type == 'oil':
        api = oil_pvt.api
        sgsp = oil_pvt.sg_sp
        pb = oil_pvt.pb
        rsb = oil_pvt.rsb
        if oil_pvt.sg_g > 0:
            gsg = oil_pvt.sg_g

    if rates is None:
        if max_rate is None:
            max_rate = 50.0 if well_type == 'gas' else 10000.0
        rates = list(np.linspace(0.01 if well_type == 'gas' else 1.0,
                                 max_rate, n_rates))

    bhp_list = []
    for rate in rates:
        if well_type == 'gas':
            bhp_val = fbhp(thp=thp, completion=completion, vlpmethod=vlpmethod,
                           well_type='gas', qg_mmscfd=rate, cgr=cgr,
                           qw_bwpd=qw_bwpd, oil_vis=oil_vis, api=api, pr=pr,
                           wsg=wsg, injection=injection, gsg=gsg)
        else:
            bhp_val = fbhp(thp=thp, completion=completion, vlpmethod=vlpmethod,
                           well_type='oil', qt_stbpd=rate, gor=gor, wc=wc,
                           wsg=wsg, injection=injection, gsg=gsg,
                           pb=pb, rsb=rsb, sgsp=sgsp, api=api, oil_pvt=oil_pvt)
        bhp_list.append(bhp_val)

    return {'rates': list(rates), 'bhp': bhp_list}


# ============================================================================
#  Public API: ipr_curve
# ============================================================================

def ipr_curve(reservoir, well_type='gas', gas_pvt=None, oil_pvt=None,
              n_points=20, min_pwf=14.7,
              wc=0, wsg=1.07, bo=1.2, uo=1.0, gsg=0.65):
    """ Returns IPR curve as dict {'pwf': [...], 'rate': [...]}.

        reservoir: Reservoir object
        well_type: 'gas', 'oil', or 'water'
        gas_pvt: GasPVT object (required for gas wells if not using defaults)
        oil_pvt: OilPVT object (optional for oil wells)
        n_points: Number of pressure points
        min_pwf: Minimum flowing BHP (psia). Defaults to 14.7
        wc: Water cut (fraction 0-1). For oil wells
        wsg: Water specific gravity
        bo: Oil FVF (rb/stb). Used for oil/water wells if oil_pvt not provided
        uo: Oil viscosity (cP). Used for oil/water wells if oil_pvt not provided
        gsg: Gas specific gravity. Used if gas_pvt not provided
    """
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
            rate_list.append(float(qg))

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

        J = 0.00708 * k * h / (uo_pr * bo_pr * (np.log(re / rw) + S - 0.75))

        for pwf in pwf_list:
            if oil_pvt is not None and pr > pb:
                # Undersaturated: Darcy from Pr to Pb, Vogel below Pb
                q_darcy = J * (pr - max(pwf, pb))
                if pwf < pb:
                    qsat_max = J * pb / 1.8
                    q_vogel = qsat_max * (1 - 0.2 * (pwf / pb) - 0.8 * (pwf / pb) ** 2)
                    rate_list.append(q_darcy + q_vogel)
                else:
                    rate_list.append(q_darcy)
            elif oil_pvt is not None and pr <= pb:
                # Saturated: Vogel
                qsat_max = J * pr / 1.8
                q = qsat_max * (1 - 0.2 * (pwf / pr) - 0.8 * (pwf / pr) ** 2)
                rate_list.append(q)
            else:
                # Simple Darcy (no Pb info)
                rate_list.append(J * (pr - pwf))

    elif well_type == 'water':
        # Water injectivity
        J = 0.00708 * k * h / (uo * bo * (np.log(re / rw) + S - 0.75))
        for pwf in pwf_list:
            rate_list.append(J * (pr - pwf))

    return {'pwf': pwf_list, 'rate': rate_list}


# ============================================================================
#  Public API: operating_point
# ============================================================================

def operating_point(thp, completion, reservoir,
                    vlpmethod='HB', well_type='gas',
                    gas_pvt=None, oil_pvt=None,
                    cgr=0, qw_bwpd=0, oil_vis=1.0, api=45,
                    gor=0, wc=0,
                    wsg=1.07, gsg=0.65, pb=0, rsb=0, sgsp=0.65,
                    bo=1.2, uo=1.0,
                    n_points=25):
    """ Returns operating point as dict with 'rate', 'bhp', 'vlp', 'ipr'.

        Finds intersection of VLP outflow curve and IPR inflow curve via bisection.

        thp: Tubing head pressure (psia)
        completion: Completion object
        reservoir: Reservoir object
        vlpmethod: VLP method string
        well_type: 'gas' or 'oil'
        Other parameters: Same as fbhp() and ipr_curve()
    """
    pr = reservoir.pr

    if oil_pvt is not None and well_type == 'oil':
        api = oil_pvt.api
        sgsp = oil_pvt.sg_sp
        pb = oil_pvt.pb
        rsb = oil_pvt.rsb
        if oil_pvt.sg_g > 0:
            gsg = oil_pvt.sg_g

    # Get IPR curve
    ipr = ipr_curve(reservoir=reservoir, well_type=well_type,
                    gas_pvt=gas_pvt, oil_pvt=oil_pvt,
                    n_points=n_points, wc=wc, wsg=wsg,
                    bo=bo, uo=uo, gsg=gsg)

    # Find AOF (maximum rate at minimum pwf)
    aof = max(ipr['rate'])
    if aof <= 0:
        return {'rate': 0.0, 'bhp': pr, 'vlp': {'rates': [], 'bhp': []},
                'ipr': ipr}

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
                           wsg=wsg, gsg=gsg)
        else:
            vlp_bhp = fbhp(thp=thp, completion=completion, vlpmethod=vlpmethod,
                           well_type='oil', qt_stbpd=rate, gor=gor, wc=wc,
                           wsg=wsg, gsg=gsg, pb=pb, rsb=rsb, sgsp=sgsp,
                           api=api, oil_pvt=oil_pvt)

        # IPR BHP: interpolate from IPR curve
        ipr_rates = ipr['rate']
        ipr_pwfs = ipr['pwf']
        # Rate decreases as pwf increases, so interpolate
        ipr_bhp = np.interp(rate, ipr_rates[::-1], ipr_pwfs[::-1])

        return vlp_bhp - ipr_bhp

    # Bisect to find operating rate (in IPR units)
    min_rate = 0.01 * gas_scale if well_type == 'gas' else 1.0
    max_rate_search = aof * 0.999

    try:
        op_rate = bisect_solve(None, _err, min_rate, max_rate_search, 1e-4)
    except (RuntimeError, ValueError):
        op_rate = 0.0

    # Calculate operating BHP
    if well_type == 'gas' and op_rate > 0:
        op_bhp = fbhp(thp=thp, completion=completion, vlpmethod=vlpmethod,
                      well_type='gas', qg_mmscfd=op_rate / gas_scale, cgr=cgr,
                      qw_bwpd=qw_bwpd, oil_vis=oil_vis, api=api, pr=pr,
                      wsg=wsg, gsg=gsg)
    elif well_type == 'oil' and op_rate > 0:
        op_bhp = fbhp(thp=thp, completion=completion, vlpmethod=vlpmethod,
                      well_type='oil', qt_stbpd=op_rate, gor=gor, wc=wc,
                      wsg=wsg, gsg=gsg, pb=pb, rsb=rsb, sgsp=sgsp,
                      api=api, oil_pvt=oil_pvt)
    else:
        op_bhp = pr

    # Generate VLP curve for output (in VLP units: MMscf/d for gas, STB/d for oil)
    vlp_max = max_rate_search / gas_scale
    vlp_min = 0.01 if well_type == 'gas' else 1.0
    vlp_rates = list(np.linspace(vlp_min, vlp_max, n_points))
    vlp = outflow_curve(thp=thp, completion=completion, vlpmethod=vlpmethod,
                        well_type=well_type, rates=vlp_rates,
                        cgr=cgr, qw_bwpd=qw_bwpd, oil_vis=oil_vis, api=api,
                        pr=pr, gor=gor, wc=wc, wsg=wsg, gsg=gsg,
                        pb=pb, rsb=rsb, sgsp=sgsp, oil_pvt=oil_pvt)

    # Convert operating rate to VLP units for return
    op_rate_out = op_rate / gas_scale if well_type == 'gas' else op_rate

    return {'rate': op_rate_out, 'bhp': op_bhp, 'vlp': vlp, 'ipr': ipr}
