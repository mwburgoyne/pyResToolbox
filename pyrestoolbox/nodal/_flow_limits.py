"""
Flow limits and wellhead chokes: Turner/Coleman critical gas rate for liquid
loading, API RP 14E erosional velocity, and the Gilbert-family critical-flow
choke correlations.
"""

import math
import warnings
from dataclasses import dataclass

import numpy as np

import pyrestoolbox.gas as gas
from pyrestoolbox.constants import (R, psc, tscr, degF2R, MW_AIR, WDEN, CUFTperBBL,
                                    BAR_TO_PSI, PSI_TO_BAR, degc_to_degf,
                                    MM_TO_IN, FT_TO_M, LBCUFT_TO_KGM3,
                                    KGM3_TO_LBCUFT, MSCF_TO_SM3, SM3_TO_MSCF,
                                    STB_TO_SM3, SM3_TO_STB,
                                    SM3_PER_SM3_TO_SCF_PER_STB)
from pyrestoolbox.shared_fns import convert_to_numpy, process_output, validate_pe_inputs
from pyrestoolbox.validate import validate_choice

__all__ = ['critical_gas_rate', 'erosional_velocity', 'choke_flow', 'choke_flow_sachdeva']

# Critical (droplet-removal) velocity, v = K [sigma (rho_l - rho_g)]^0.25 / rho_g^0.5,
# v in ft/s, sigma in dyne/cm, densities in lbm/ft3.
# Coleman, Clay, McCurdy & Norris (1991), JPT 43(3) 329, SPE-20280-PA: K = 1.593 (no uplift).
# Turner, Hubbard & Dukler (1969), JPT 21(11) 1475, SPE-2198-PA: the same droplet model
# with the 20% uplift Turner applied to fit his wellhead data, K = 1.2 x 1.593 = 1.912.
_COLEMAN_K = 1.593
_TURNER_UPLIFT = 1.2
# Liquid properties Turner et al. (1969) used for each liquid phase
_TURNER_WATER = (67.0, 60.0)       # (density lbm/ft3, interfacial tension dyne/cm)
_TURNER_CONDENSATE = (45.0, 20.0)

# API RP 14E (5th ed., 1991), Section 2.5: erosional velocity Ve = C / sqrt(rho_m) (ft/s,
# rho_m lbm/ft3); mixture density (Eq. 2.14) and minimum flow area (Eq. 2.15).
_RP14E_LIQ = 12409.0               # rho_m numerator, liquid term
_RP14E_GAS = 2.7                   # rho_m numerator, gas term
_RP14E_DEN = 198.7                 # rho_m denominator, liquid term
_RP14E_AREA_LIQ = 9.35             # min area, in2 per 1000 bbl/d, liquid term
_RP14E_AREA_GAS = 21.25            # min area, gas-term divisor

# Gilbert-family critical-flow chokes, p_up = A q GLR^B / S^C with p_up in psig,
# q gross liquid STB/d, GLR scf/STB, S bean size in 64ths of an inch. Coefficients as
# tabulated with GLR in scf/STB (Guo, Lyons & Ghalambor, Petroleum Production
# Engineering, 2007, Table 5.1); Gilbert's 10.00 is his 435 with GLR in Mscf/STB.
# Gilbert (1954), API Drilling and Production Practice 126; Ros (1960), Shell report;
# Baxendell (1958), API Drilling and Production Practice 145; Achong (1961), Shell report.
_CHOKE_COEFFS = {
    'GILBERT': (10.00, 0.546, 1.89),
    'ROS': (17.40, 0.500, 2.00),
    'BAXENDELL': (9.56, 0.546, 1.93),
    'ACHONG': (3.82, 0.650, 1.88),
}
# Gilbert (1954): the correlation holds while upstream pressure is at least 1.7 times
# downstream (critical flow), i.e. p_down / p_up <= 1 / 1.7.
_CHOKE_CRITICAL_RATIO = 1.0 / 1.7

# Sachdeva, Schmidt, Brill & Blais (1986), Two-Phase Flow Through Chokes, SPE-15657-MS.
# Discharge coefficient: 0.75 for a choke in a housing with an elbow upstream (usual in
# the field), 0.85 for undisturbed upstream flow.
_SACHDEVA_CD = 0.75
# Liquid specific heat, Btu/(lbm degF), the value of the paper's Fig. 5 base case
_SACHDEVA_CL = 0.8
# Ideal-gas heat-capacity ratio Cp/Cv of natural gas (the model's adiabatic exponent)
_SACHDEVA_K = 1.3
_GC = 32.174                         # lbm ft / (lbf s^2)
_R_BTU = 1.98588                     # Btu / (lbmol degR)
_KJ_PER_KG_K_PER_BTU_LBM_F = 4.1868

_SEC_PER_DAY = 86400.0
_SQIN_PER_SQFT = 144.0


def _gas_density_z(p, degf, sg, co2, h2s, n2, h2, zmethod, cmethod):
    """Z and in-situ gas density (lbm/ft3) at field-unit p (psia), degf."""
    z = np.asarray(gas.gas_z(p=p, sg=sg, degf=degf, zmethod=zmethod, cmethod=cmethod,
                             co2=co2, h2s=h2s, n2=n2, h2=h2), dtype=float)
    rho = p * sg * MW_AIR / (z * R * (degf + degF2R))
    return z, rho


def _std_rate_mscfd(v, area_sqft, p, z, degf):
    """Standard gas rate (Mscf/d) for in-situ velocity v (ft/s) through area_sqft."""
    return v * area_sqft * _SEC_PER_DAY * (p / psc) * (tscr / (z * (degf + degF2R))) / 1000.0


def critical_gas_rate(p, degf, sg, tid, method='TURNER', liquid='water',
                      rho_liq=None, ift=None, co2=0, h2s=0, n2=0, h2=0,
                      zmethod='DAK', cmethod='PMC', metric=False):
    """Minimum gas rate to lift liquid droplets (liquid-loading onset).

    Droplet model v_c = K [ift (rho_liq - rho_g)]^0.25 / rho_g^0.5, with K = 1.912
    for Turner et al. (1969) (their 20% uplift on the derived 1.593) and K = 1.593 for
    Coleman et al. (1991), who found the uplift unnecessary at the low wellhead
    pressures (below about 500 psia) of their data. Evaluate at the wellhead and at
    the bottom of the tubing; the larger rate controls.

    p: Flowing pressure (psia | barsa). Scalar or array
    degf: Flowing temperature (deg F | deg C)
    sg: Gas SG relative to air
    tid: Tubing inside diameter (in | mm)
    method: 'TURNER' (default) or 'COLEMAN'
    liquid: 'water' (default; 67 lbm/ft3, 60 dyne/cm) or 'condensate' (45 lbm/ft3,
        20 dyne/cm), the values Turner et al. used. Sets defaults for rho_liq and ift
    rho_liq: Liquid density, overrides the preset (lbm/ft3 | kg/m3)
    ift: Gas-liquid interfacial tension, overrides the preset (dyne/cm | mN/m)
    co2, h2s, n2, h2: Mole fractions in the gas
    zmethod, cmethod: Z-factor and critical-property methods, as for gas.gas_z
    metric: If True, inputs and outputs in Eclipse METRIC units

    Returns dict with 'v_crit' (ft/s | m/s), 'q_crit' (Mscf/d | sm3/d), 'z' and
    'rho_gas' (lbm/ft3 | kg/m3), each scalar or array following p.
    """
    validate_choice(method, ('TURNER', 'COLEMAN'), 'method')
    validate_choice(liquid, ('water', 'condensate'), 'liquid')
    p, is_list = convert_to_numpy(p)
    if metric:
        p = p * BAR_TO_PSI
        degf = degc_to_degf(degf)
        tid = tid * MM_TO_IN
        if rho_liq is not None:
            rho_liq = rho_liq * KGM3_TO_LBCUFT
    validate_pe_inputs(p=p, degf=degf, sg=sg, co2=co2, h2s=h2s, n2=n2, h2=h2)
    if tid <= 0:
        raise ValueError(f"tid must be positive, got {tid}")

    rho_default, ift_default = _TURNER_WATER if liquid == 'water' else _TURNER_CONDENSATE
    rho_liq = rho_default if rho_liq is None else rho_liq
    ift = ift_default if ift is None else ift
    if ift <= 0:
        raise ValueError(f"ift must be positive, got {ift}")

    z, rho_g = _gas_density_z(p, degf, sg, co2, h2s, n2, h2, zmethod, cmethod)
    if np.any(rho_liq <= rho_g):
        raise ValueError("Liquid density must exceed gas density; no droplet "
                         "model applies at these conditions")
    k = _COLEMAN_K * (_TURNER_UPLIFT if method == 'TURNER' else 1.0)
    v = k * (ift * (rho_liq - rho_g)) ** 0.25 / np.sqrt(rho_g)
    area = math.pi * tid ** 2 / 4.0 / _SQIN_PER_SQFT
    q = _std_rate_mscfd(v, area, p, z, degf)

    if metric:
        v, q, rho_g = v * FT_TO_M, q * MSCF_TO_SM3, rho_g * LBCUFT_TO_KGM3
    return {'v_crit': process_output(v, is_list), 'q_crit': process_output(q, is_list),
            'z': process_output(z, is_list), 'rho_gas': process_output(rho_g, is_list)}


def _erosional_inputs(p, degf, glr, tid, metric):
    """Convert erosional_velocity inputs to field units and check glr and tid."""
    if metric:
        p, degf = p * BAR_TO_PSI, degc_to_degf(degf)
        tid = None if tid is None else tid * MM_TO_IN
        glr = None if glr is None else glr * SM3_PER_SM3_TO_SCF_PER_STB
    if glr is not None and glr < 0:
        raise ValueError(f"glr must be non-negative, got {glr}")
    if tid is not None and tid <= 0:
        raise ValueError(f"tid must be positive, got {tid}")
    return p, degf, glr, tid


def erosional_velocity(p, degf, sg, glr=None, sg_liq=1.0, c=100.0, tid=None,
                       co2=0, h2s=0, n2=0, h2=0, zmethod='DAK', cmethod='PMC',
                       metric=False):
    """API RP 14E erosional velocity, Ve = C / sqrt(rho_m).

    rho_m is the RP 14E no-slip mixture density (Eq. 2.14) when glr is given, or the
    gas density for dry gas (glr=None). C is the empirical constant in field units,
    (lbm/ft s^2)^0.5, whichever unit system is used: RP 14E gives 100 for continuous
    service and 125 for intermittent; 150-200 is used for solids-free, non-corrosive
    service. The constant is a screening limit, not a measured erosion rate.

    p: Flowing pressure (psia | barsa). Scalar or array
    degf: Flowing temperature (deg F | deg C)
    sg: Gas SG relative to air
    glr: Gas-liquid ratio at standard conditions (scf/STB | sm3/sm3); None for dry gas
    sg_liq: Liquid SG relative to water (default 1.0)
    c: RP 14E C-factor (default 100)
    tid: Optional pipe inside diameter (in | mm). When given, also returns the
        maximum rate at Ve: liquid rate from RP 14E Eq. 2.15 when glr is given,
        gas rate otherwise
    co2, h2s, n2, h2: Mole fractions in the gas
    zmethod, cmethod: Z-factor and critical-property methods, as for gas.gas_z
    metric: If True, inputs and outputs in Eclipse METRIC units

    Returns dict with 've' (ft/s | m/s), 'rho_mix' (lbm/ft3 | kg/m3) and 'z', plus
    'q_liq_max' (STB/d | sm3/d) or 'q_gas_max' (Mscf/d | sm3/d) when tid is given.
    """
    p, is_list = convert_to_numpy(p)
    p, degf, glr, tid = _erosional_inputs(p, degf, glr, tid, metric)
    validate_pe_inputs(p=p, degf=degf, sg=sg, co2=co2, h2s=h2s, n2=n2, h2=h2)
    if c <= 0 or sg_liq <= 0:
        raise ValueError(f"c and sg_liq must be positive, got {c} and {sg_liq}")

    z, rho_g = _gas_density_z(p, degf, sg, co2, h2s, n2, h2, zmethod, cmethod)
    t_r = degf + degF2R
    if glr is None:
        rho_m = rho_g
    else:
        rho_m = ((_RP14E_LIQ * sg_liq * p + _RP14E_GAS * glr * sg * p)
                 / (_RP14E_DEN * p + glr * t_r * z))
    ve = c / np.sqrt(rho_m)

    out = {}
    if tid is not None:
        area_sqin = math.pi * tid ** 2 / 4.0
        if glr is None:
            q = _std_rate_mscfd(ve, area_sqin / _SQIN_PER_SQFT, p, z, degf)
            out['q_gas_max'] = q * MSCF_TO_SM3 if metric else q
        else:
            area_per_kbpd = (_RP14E_AREA_LIQ + z * glr * t_r / (_RP14E_AREA_GAS * p)) / ve
            q = 1000.0 * area_sqin / area_per_kbpd
            out['q_liq_max'] = q * STB_TO_SM3 if metric else q
    if metric:
        ve, rho_m = ve * FT_TO_M, rho_m * LBCUFT_TO_KGM3
    out = {'ve': ve, 'rho_mix': rho_m, 'z': z, **out}
    return {k: process_output(v, is_list) for k, v in out.items()}


def _choke_inputs(glr, p_up, q_liq, size, p_down, metric):
    """Check the choke inputs and convert them to field units (psia, STB/d, scf/STB)."""
    if sum(x is None for x in (p_up, q_liq, size)) != 1:
        raise ValueError("Give exactly two of p_up, q_liq and size")
    if metric:
        glr = glr * SM3_PER_SM3_TO_SCF_PER_STB
        p_up = None if p_up is None else p_up * BAR_TO_PSI
        p_down = None if p_down is None else p_down * BAR_TO_PSI
        q_liq = None if q_liq is None else q_liq * SM3_TO_STB
    if glr <= 0:
        raise ValueError(f"glr must be positive, got {glr}")
    for name, val in (('q_liq', q_liq), ('size', size)):
        if val is not None and val <= 0:
            raise ValueError(f"{name} must be positive, got {val}")
    if p_up is not None and p_up <= psc:
        raise ValueError(f"p_up must exceed atmospheric ({psc} psia), got {p_up} psia")
    if p_down is not None:
        validate_pe_inputs(p=p_down)
    return glr, p_up, q_liq, size, p_down


def choke_flow(glr, p_up=None, q_liq=None, size=None, method='GILBERT',
               p_down=None, metric=False):
    """Gilbert-family critical-flow choke: p_up = A q GLR^B / S^C.

    Give exactly two of p_up, q_liq and size; the third is solved. The correlations
    are empirical fits to critical (choked) flow of oil wells and do not see the
    downstream pressure. When p_down is given and p_down / p_up exceeds 1/1.7
    (Gilbert's criterion, absolute pressures), the flow is subcritical and the
    result is not valid; a warning is issued and 'critical' is False.

    glr: Gross gas-liquid ratio (scf/STB | sm3/sm3)
    p_up: Upstream (wellhead) pressure (psia | barsa)
    q_liq: Gross liquid rate (STB/d | sm3/d)
    size: Bean size in 64ths of an inch (both unit systems)
    method: 'GILBERT' (default), 'ROS', 'BAXENDELL' or 'ACHONG'
    p_down: Optional downstream pressure for the critical-flow check (psia | barsa)
    metric: If True, inputs and outputs in Eclipse METRIC units

    Returns dict with 'p_up' (psia | barsa), 'q_liq' (STB/d | sm3/d), 'size' (64ths)
    and 'critical' (True, False, or None when p_down is not given).
    """
    method = validate_choice(str(method).upper(), tuple(_CHOKE_COEFFS), 'method')
    glr, p_up, q_liq, size, p_down = _choke_inputs(glr, p_up, q_liq, size, p_down, metric)

    a, b, c = _CHOKE_COEFFS[method]
    if p_up is None:
        p_up = a * q_liq * glr ** b / size ** c + psc
    elif q_liq is None:
        q_liq = (p_up - psc) * size ** c / (a * glr ** b)
    else:
        size = (a * q_liq * glr ** b / (p_up - psc)) ** (1.0 / c)

    critical = None
    if p_down is not None:
        critical = p_down / p_up <= _CHOKE_CRITICAL_RATIO
        if not critical:
            warnings.warn(f"p_down/p_up = {p_down / p_up:.3f} exceeds the critical-flow "
                          f"limit {_CHOKE_CRITICAL_RATIO:.3f}; the {method} correlation "
                          f"assumes critical flow and is not valid here")
    if metric:
        p_up, q_liq = p_up * PSI_TO_BAR, q_liq * STB_TO_SM3
    return {'p_up': p_up, 'q_liq': q_liq, 'size': size, 'critical': critical}


@dataclass
class _ChokeFluid:
    """Upstream two-phase state for the Sachdeva model, field units."""
    x: float          # free-gas mass quality
    vg1: float        # upstream gas specific volume, ft3/lbm
    vl: float         # liquid specific volume, ft3/lbm
    n: float          # polytropic exponent, Eq. (A.18)
    k: float
    mass_per_stb: float   # lbm of liquid plus free gas per STB of liquid


def _sachdeva_fluid(p1, degf, sg, glr, sg_liq, k, cl, zkw):
    """Quality, specific volumes and polytropic exponent at upstream conditions."""
    mw = sg * MW_AIR
    z = float(gas.gas_z(p=p1, sg=sg, degf=degf, **zkw))
    m_liq = CUFTperBBL * sg_liq * WDEN
    m_gas = glr * psc * mw / (R * tscr)
    x = m_gas / (m_gas + m_liq)
    cv = _R_BTU / (mw * (k - 1.0))
    cp = k * cv
    n = 1.0 + x * (cp - cv) / (x * cv + (1.0 - x) * cl)
    return _ChokeFluid(x=x, vg1=z * R * (degf + degF2R) / (p1 * mw),
                       vl=1.0 / (sg_liq * WDEN), n=n, k=k, mass_per_stb=m_gas + m_liq)


def _sachdeva_yc(fl):
    """Critical pressure ratio from Eq. (2)/(A.25), solved for y (V_G2 depends on y)."""
    from scipy.optimize import brentq
    k, x, vl, vg1, n = fl.k, fl.x, fl.vl, fl.vg1, fl.n
    kk = k / (k - 1.0)

    def resid(y):
        a = (1.0 - x) * vl / (x * vg1 * y ** (-1.0 / k))
        num = kk + (1.0 - x) * vl * (1.0 - y) / (x * vg1)
        return y - (num / (kk + n / 2.0 + n * a + n / 2.0 * a * a)) ** kk
    return brentq(resid, 1e-9, 1.0 - 1e-12, xtol=1e-14)


def _sachdeva_flux(y, p1, fl, cd):
    """Mass flux at the throat, lbm/(ft2 s), Eq. (4)/(A.27) with Eqs. (6)-(7)."""
    k, x = fl.k, fl.x
    vg2 = fl.vg1 * y ** (-1.0 / k)
    rho_m2 = 1.0 / (x * vg2 + (1.0 - x) * fl.vl)
    bracket = ((1.0 - x) * (1.0 - y) * fl.vl
               + x * k / (k - 1.0) * (fl.vg1 - y * vg2))
    return cd * math.sqrt(2.0 * _GC * _SQIN_PER_SQFT * p1 * rho_m2 ** 2 * bracket)


def _sachdeva_rate(p1, p2, size, fluid_at, cd):
    """Liquid rate (STB/d) through a size/64 in choke, critical ratio and regime."""
    fl = fluid_at(p1)
    yc = _sachdeva_yc(fl)
    y = max(p2 / p1, yc)
    area = math.pi / 4.0 * (size / 64.0 / 12.0) ** 2
    q = _sachdeva_flux(y, p1, fl, cd) * area / fl.mass_per_stb * _SEC_PER_DAY
    return q, yc, fl


def _sachdeva_inputs(p_up, p_down, q_liq, glr, degf, cl, metric):
    """Field-unit conversion and checks for choke_flow_sachdeva."""
    if metric:
        p_up = None if p_up is None else p_up * BAR_TO_PSI
        p_down = p_down * BAR_TO_PSI
        q_liq = None if q_liq is None else q_liq * SM3_TO_STB
        glr = glr * SM3_PER_SM3_TO_SCF_PER_STB
        degf = degc_to_degf(degf)
        cl = None if cl is None else cl / _KJ_PER_KG_K_PER_BTU_LBM_F
    cl = _SACHDEVA_CL if cl is None else cl
    if glr <= 0:
        raise ValueError(f"glr must be positive (the model needs free gas), got {glr}")
    if q_liq is not None and q_liq <= 0:
        raise ValueError(f"q_liq must be positive, got {q_liq}")
    if cl <= 0:
        raise ValueError(f"cl must be positive, got {cl}")
    validate_pe_inputs(p=p_down, degf=degf)
    if p_up is not None:
        validate_pe_inputs(p=p_up)
        if p_up <= p_down:
            raise ValueError(f"p_up ({p_up:.6g} psia) must exceed p_down ({p_down:.6g} psia)")
    return p_up, p_down, q_liq, glr, degf, cl


def choke_flow_sachdeva(p_down, degf, glr, p_up=None, q_liq=None, size=None, sg=0.65,
                        sg_liq=1.0, cd=_SACHDEVA_CD, k=_SACHDEVA_K, cl=None,
                        co2=0, h2s=0, n2=0, h2=0, zmethod='DAK', cmethod='PMC',
                        metric=False):
    """Sachdeva et al. (1986) mechanistic two-phase choke model, critical and subcritical.

    Give exactly two of p_up, q_liq and size; the third is solved. The critical
    pressure ratio y_c comes from Eq. (2) (iterated, as V_G2 depends on y); when
    p_down / p_up <= y_c the flow is critical and the throat ratio is y_c, otherwise
    it is subcritical at the actual ratio. The mass flux is Eq. (4) with the
    no-slip throat density of Eq. (7) and the polytropic exponent of Eq. (A.18).
    Free gas only: dissolved gas and liquid shrinkage are ignored, as in the paper,
    so glr is the free-gas ratio at upstream conditions and q_liq the liquid rate at
    them. The gas expands adiabatically with an ideal-gas exponent k (not the
    real-gas Cp/Cv); Cv per lbm is R / (MW (k - 1)).

    p_down: Downstream pressure (psia | barsa)
    degf: Upstream temperature (deg F | deg C)
    glr: Free gas-liquid ratio (scf/STB | sm3/sm3), > 0
    p_up: Upstream pressure (psia | barsa)
    q_liq: Liquid rate (STB/d | sm3/d)
    size: Choke diameter in 64ths of an inch (both unit systems)
    sg: Gas SG relative to air
    sg_liq: Liquid SG relative to water (default 1.0)
    cd: Discharge coefficient: 0.75 (default) for a choke in a housing with an
        upstream elbow, 0.85 for undisturbed upstream flow (Sachdeva et al.)
    k: Gas heat-capacity ratio Cp/Cv (default 1.3)
    cl: Liquid specific heat (Btu/(lbm degF) | kJ/(kg K)); default 0.8 Btu/(lbm degF),
        the paper's base case. Water is about 1.0, crude oil about 0.5
    co2, h2s, n2, h2: Mole fractions in the gas
    zmethod, cmethod: Z-factor and critical-property methods, as for gas.gas_z
    metric: If True, inputs and outputs in Eclipse METRIC units

    Returns dict with 'p_up' (psia | barsa), 'q_liq' (STB/d | sm3/d), 'q_gas'
    (Mscf/d | sm3/d), 'size' (64ths), 'critical' (bool), 'y_crit' (critical
    downstream/upstream ratio), 'y' (actual ratio) and 'x' (free-gas mass quality).
    """
    if sum(v is None for v in (p_up, q_liq, size)) != 1:
        raise ValueError("Give exactly two of p_up, q_liq and size")
    if size is not None and size <= 0:
        raise ValueError(f"size must be positive, got {size}")
    if not 0 < cd <= 1.0 or k <= 1.0 or sg_liq <= 0:
        raise ValueError(f"Need 0 < cd <= 1, k > 1 and sg_liq > 0, got {cd}, {k}, {sg_liq}")
    p_up, p_down, q_liq, glr, degf, cl = _sachdeva_inputs(
        p_up, p_down, q_liq, glr, degf, cl, metric)
    validate_pe_inputs(sg=sg, co2=co2, h2s=h2s, n2=n2, h2=h2)
    zkw = dict(zmethod=zmethod, cmethod=cmethod, co2=co2, h2s=h2s, n2=n2, h2=h2)

    def fluid_at(p1):
        return _sachdeva_fluid(p1, degf, sg, glr, sg_liq, k, cl, zkw)

    if q_liq is None:
        q_liq, yc, fl = _sachdeva_rate(p_up, p_down, size, fluid_at, cd)
    elif size is None:
        q_one, yc, fl = _sachdeva_rate(p_up, p_down, 1.0, fluid_at, cd)
        size = math.sqrt(q_liq / q_one)          # rate scales with throat area
    else:
        # Bracket probes reach pressures far from the answer; their Z-method range
        # warnings are noise. The final evaluation below warns if the answer is out.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            p_up = _sachdeva_p_up(p_down, q_liq, size, fluid_at, cd)
        q_liq, yc, fl = _sachdeva_rate(p_up, p_down, size, fluid_at, cd)

    out = {'p_up': p_up, 'q_liq': q_liq, 'q_gas': q_liq * glr / 1000.0, 'size': size,
           'critical': p_down / p_up <= yc, 'y_crit': yc, 'y': p_down / p_up, 'x': fl.x}
    if metric:
        out['p_up'] *= PSI_TO_BAR
        out['q_liq'] *= STB_TO_SM3
        out['q_gas'] *= MSCF_TO_SM3
    return out


def _sachdeva_p_up(p_down, q_liq, size, fluid_at, cd):
    """Upstream pressure delivering q_liq; the rate rises with p_up at fixed p_down."""
    from scipy.optimize import brentq

    def resid(p1):
        return _sachdeva_rate(p1, p_down, size, fluid_at, cd)[0] - q_liq
    lo, hi = p_down * (1.0 + 1e-9), 2.0 * p_down
    while resid(hi) < 0:
        if hi > 1e5:
            raise ValueError(f"No upstream pressure below 100,000 psia delivers {q_liq:.6g} STB/d")
        hi *= 2.0
    return brentq(resid, lo, hi, xtol=1e-10, rtol=1e-12)
