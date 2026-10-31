"""
Steady-state, isothermal gas flow in a single pipeline segment: the General Flow
equation with Colebrook-White or AGA transmission factors, and the Weymouth,
Panhandle A, Panhandle B and IGT equations, each with the elevation correction.

Equations and constants from Menon, E.S. (2005), Gas Pipeline Hydraulics, CRC Press
(doi 10.1201/9781420038224), Chapter 2, in USCS form: Q in scf/d, D in inches, L in
miles, pressures psia, temperatures deg R, viscosity lb/ft-s. Equation numbers below
are those of Menon's Chapter 2 as reproduced in his course notes P06-001.
"""

import math
import warnings
from dataclasses import dataclass
from typing import Callable, Optional

import pyrestoolbox.gas as gas
from pyrestoolbox.constants import (psc, tscr, degF2R, BAR_TO_PSI, PSI_TO_BAR,
                                    degc_to_degf, MM_TO_IN, IN_TO_MM, FT_TO_M, M_TO_FT,
                                    MSCF_TO_SM3, SM3_TO_MSCF)
from pyrestoolbox.shared_fns import validate_pe_inputs
from pyrestoolbox.validate import validate_choice

__all__ = ['gas_pipeline']

# General Flow equation with transmission factor F = 2 / sqrt(f_darcy), Eq. (2.15):
# Q = 38.77 F (Tb/Pb) [(P1^2 - e^s P2^2) / (G Tf Le Z)]^0.5 D^2.5
_GENERAL_K = 38.77
# Empirical equations, Q = E K (Tb/Pb)^a [(P1^2 - e^s P2^2) / (G^g Tf Le Z^z mu^u)]^m D^n
# as (K, a, g, z, u, m, n). Weymouth Eq. (2.30); Panhandle A Eq. (2.34); Panhandle B
# Eq. (2.36); IGT Eq. (2.38), which carries viscosity (lb/ft-s) and no Z.
_EMPIRICAL = {
    'WEYMOUTH': (433.5, 1.0, 1.0, 1.0, 0.0, 0.5, 2.667),
    'PANHANDLE_A': (435.87, 1.0788, 0.8539, 1.0, 0.0, 0.5394, 2.6182),
    'PANHANDLE_B': (737.0, 1.02, 0.961, 1.0, 0.0, 0.51, 2.53),
    'IGT': (136.9, 1.0, 0.8, 0.0, 0.2, 0.555, 2.667),
}
_METHODS = ('COLEBROOK', 'AGA') + tuple(_EMPIRICAL)

# Reynolds number Re = 0.0004778 (Pb/Tb) G Q / (mu D), Eq. (2.18)
_RE_K = 0.0004778
# Laminar below Re = 2000, f = 64 / Re, Eq. (2.20)
_LAMINAR_RE = 2000.0
_LAMINAR_K = 64.0
# Colebrook-White, 1/sqrt(f) = -2 log10(e / (3.7 D) + 2.51 / (Re sqrt(f))), Eq. (2.21)
_CW_ROUGH = 3.7
_CW_SMOOTH = 2.51
# AGA NB-13: rough-pipe F = 4 log10(3.7 D / e), Eq. (2.22); partially turbulent
# F = 4 Df log10(Re / (1.4125 Ft)), Eq. (2.23), with the Von Karman smooth-pipe factor
# Ft = 4 log10(Re / Ft) - 0.6, Eq. (2.24)
_AGA_PT = 1.4125
_AGA_SMOOTH_OFFSET = 0.6
# Elevation adjustment s = 0.0375 G (H2 - H1) / (Tf Z), H in ft, Eq. (2.12)
_ELEV_K = 0.0375

# Bare steel absolute roughness, 700 micro-inches (Menon, after Eq. 2.21)
_BARE_STEEL_ROUGHNESS = 0.0007
_FT_PER_MILE = 5280.0
# cP to lb/(ft s): 1e-3 Pa s x 0.3048 m/ft / 0.45359237 kg/lb
_CP_TO_LB_FT_S = 1.0e-3 * FT_TO_M / 0.45359237
_ITER_TOL = 1.0e-12
_MAX_ITER = 100


@dataclass
class _Pipe:
    """Field-unit pipeline description passed through the internal solvers."""
    d: float                # inside diameter, in
    l_mi: float             # length, miles
    tf: float               # flowing temperature, deg R
    sg: float
    method: str
    e: float = 1.0          # pipeline efficiency
    rough: float = 0.0007   # absolute roughness, in
    df: float = 0.96        # AGA drag factor
    dz: float = 0.0         # outlet minus inlet elevation, ft
    pb: float = psc
    tb: float = tscr
    f_fixed: Optional[float] = None   # General Flow equation with a given Darcy f


def _pavg(p1, p2):
    """Average pipe pressure, Eq. (2.25)."""
    return 2.0 / 3.0 * (p1 + p2 - p1 * p2 / (p1 + p2))


def _elevation(pipe, z):
    """Elevation parameter s and equivalent length Le (miles), Eqs. (2.11)-(2.12)."""
    s = _ELEV_K * pipe.sg * pipe.dz / (pipe.tf * z)
    le = pipe.l_mi if abs(s) < 1e-12 else pipe.l_mi * math.expm1(s) / s
    return s, le


def _colebrook_f(re, rel_rough):
    """Darcy friction factor from Colebrook-White by fixed-point iteration."""
    inv = 7.0
    for _ in range(_MAX_ITER):
        new = -2.0 * math.log10(rel_rough / _CW_ROUGH + _CW_SMOOTH * inv / re)
        if abs(new - inv) < _ITER_TOL * new:
            break
        inv = new
    return 1.0 / new ** 2


def _von_karman_ft(re):
    """Smooth-pipe transmission factor Ft = 4 log10(Re / Ft) - 0.6, Eq. (2.24)."""
    ft = 20.0
    for _ in range(_MAX_ITER):
        new = 4.0 * math.log10(re / ft) - _AGA_SMOOTH_OFFSET
        if abs(new - ft) < _ITER_TOL * new:
            break
        ft = new
    return new


def _aga_transmission(re, pipe):
    """AGA transmission factor: the smaller of the rough-pipe and partially
    turbulent values, Eqs. (2.22)-(2.24)."""
    f_rough = 4.0 * math.log10(_CW_ROUGH * pipe.d / pipe.rough)
    f_partial = 4.0 * pipe.df * math.log10(re / (_AGA_PT * _von_karman_ft(re)))
    return min(f_rough, f_partial)


def _transmission(re, pipe):
    """Transmission factor F = 2 / sqrt(f) for the General Flow equation."""
    if pipe.f_fixed is not None:
        return 2.0 / math.sqrt(pipe.f_fixed)
    if re <= _LAMINAR_RE:
        return 2.0 / math.sqrt(_LAMINAR_K / re)
    if pipe.method == 'AGA':
        return _aga_transmission(re, pipe)
    return 2.0 / math.sqrt(_colebrook_f(re, pipe.rough / pipe.d))


def _reynolds(q, pipe, mu):
    return _RE_K * (pipe.pb / pipe.tb) * pipe.sg * abs(q) / (mu * pipe.d)


def _general_q(f_trans, drive, pipe, le, z):
    """General Flow equation (2.15) rate (scf/d) for drive = P1^2 - e^s P2^2."""
    root = math.copysign(math.sqrt(abs(drive) / (pipe.sg * pipe.tf * le * z)), drive)
    return pipe.e * _GENERAL_K * f_trans * (pipe.tb / pipe.pb) * root * pipe.d ** 2.5


def _empirical_q(drive, pipe, le, z, mu):
    k, a, g, ze, u, m, n = _EMPIRICAL[pipe.method]
    denom = pipe.sg ** g * pipe.tf * le * z ** ze * mu ** u
    return (pipe.e * k * (pipe.tb / pipe.pb) ** a
            * math.copysign((abs(drive) / denom) ** m, drive) * pipe.d ** n)


def _rate(p1, p2, pipe, props):
    """Standard rate (scf/d, negative for reverse flow) and the flow state, for
    upstream p1 and downstream p2 (psia). props(pavg) returns (Z, mu lb/ft-s)."""
    pavg = _pavg(p1, p2)
    z, mu = props(pavg)
    s, le = _elevation(pipe, z)
    drive = p1 ** 2 - math.exp(s) * p2 ** 2
    if drive == 0.0:
        return 0.0, dict(pavg=pavg, z=z, mu=mu, s=s, re=0.0, f_trans=math.nan)
    if pipe.method in _EMPIRICAL and pipe.f_fixed is None:
        q = _empirical_q(drive, pipe, le, z, mu)
    else:
        f_trans = _transmission(1.0e7, pipe)
        for _ in range(_MAX_ITER):
            q = _general_q(f_trans, drive, pipe, le, z)
            new = _transmission(_reynolds(q, pipe, mu), pipe)
            if abs(new - f_trans) < _ITER_TOL * new:
                break
            f_trans = new
    re = _reynolds(q, pipe, mu)
    # Transmission factor the rate implies in the General Flow equation (for the
    # empirical equations, their equivalent F)
    f_trans = q / _general_q(1.0, drive, pipe, le, z)
    return q, dict(pavg=pavg, z=z, mu=mu, s=s, re=re, f_trans=f_trans)


def _bracket_root(fn, lo, hi, grow_hi, limit):
    """Brent root of fn on [lo, hi], expanding hi by grow_hi until the sign changes."""
    from scipy.optimize import brentq
    f_lo, f_hi = fn(lo), fn(hi)
    while f_lo * f_hi > 0:
        if hi >= limit:
            return None
        hi = min(hi * grow_hi, limit)
        f_hi = fn(hi)
    return brentq(fn, lo, hi, xtol=1e-10, rtol=1e-12)


_P_FLOOR = 1.0   # psia, lowest pressure the solvers try
_P_CEIL = 1.0e5  # psia


def _solve_p2(p1, q, pipe, props):
    """Downstream pressure delivering q (scf/d) from upstream p1."""
    p2 = _bracket_root(lambda p: _rate(p1, p, pipe, props)[0] - q,
                       _P_FLOOR, p1, 2.0, _P_CEIL)
    if p2 is None or _rate(p1, _P_FLOOR, pipe, props)[0] < q:
        raise ValueError(f"Rate {q / 1000:.6g} Mscf/d exceeds the line's capacity from "
                         f"p1 = {p1:.6g} psia")
    return p2


def _solve_p1(p2, q, pipe, props):
    """Upstream pressure needed to deliver q (scf/d) at downstream p2."""
    def resid(p):
        return _rate(p, p2, pipe, props)[0] - q
    lo = p2
    while resid(lo) > 0 and lo > _P_FLOOR:   # downhill line flowing at p1 = p2
        lo = max(lo / 2.0, _P_FLOOR)
    p1 = _bracket_root(resid, lo, max(2.0 * p2, lo * 2.0), 2.0, _P_CEIL)
    if p1 is None:
        raise ValueError(f"No upstream pressure below {_P_CEIL:.0f} psia delivers "
                         f"{q / 1000:.6g} Mscf/d")
    return p1


def _to_field(metric, p1, p2, qg, degf, d, length, roughness, dz):
    if not metric:
        return p1, p2, qg, degf, d, length, roughness, dz
    conv = (lambda v, k: None if v is None else v * k)
    return (conv(p1, BAR_TO_PSI), conv(p2, BAR_TO_PSI), conv(qg, SM3_TO_MSCF),
            degc_to_degf(degf), d * MM_TO_IN, length * M_TO_FT, roughness * MM_TO_IN,
            dz * M_TO_FT)


def _check_inputs(p1, p2, qg, d, length, efficiency, roughness, drag_factor):
    if sum(v is None for v in (p1, p2, qg)) != 1:
        raise ValueError("Give exactly two of p1, p2 and qg; the missing one is solved for")
    if d <= 0 or length <= 0:
        raise ValueError(f"d and length must be positive, got {d} and {length}")
    if not 0 < efficiency <= 1.0:
        raise ValueError(f"efficiency must be in (0, 1], got {efficiency}")
    if roughness <= 0 or not 0 < drag_factor <= 1.0:
        raise ValueError(f"roughness must be positive and drag_factor in (0, 1], got "
                         f"{roughness} and {drag_factor}")
    if qg is not None and qg <= 0:
        raise ValueError(f"qg must be positive, got {qg}")


def _gas_props(degf, sg, co2, h2s, n2, h2, zmethod, cmethod) -> Callable:
    """props(pavg) -> (Z, viscosity lb/ft-s) at the flowing temperature."""
    kw = dict(sg=sg, degf=degf, zmethod=zmethod, cmethod=cmethod,
              co2=co2, h2s=h2s, n2=n2, h2=h2)

    def props(pavg):
        z = float(gas.gas_z(p=pavg, **kw))
        mu = float(gas.gas_ug(p=pavg, **kw)) * _CP_TO_LB_FT_S
        return z, mu
    return props


def gas_pipeline(d, length, degf, sg=0.65, p1=None, p2=None, qg=None,
                 method='COLEBROOK', efficiency=1.0, roughness=None,
                 drag_factor=0.96, dz=0.0, co2=0, h2s=0, n2=0, h2=0,
                 zmethod='DAK', cmethod='PMC', metric=False):
    """Steady isothermal gas flow in one pipeline segment; solves for whichever one
    of p1, p2 and qg is not given.

    Equations as in Menon (2005), Gas Pipeline Hydraulics, Chapter 2. Z (and
    viscosity, where used) are evaluated at the flowing temperature and the average
    pressure 2/3 (p1 + p2 - p1 p2 / (p1 + p2)). Standard conditions are the
    library's 14.696 psia and 60 deg F, not Menon's 14.7 psia. The segment is
    treated as one length at one temperature: for long lines with large pressure
    drop, split the line and chain the calls (outlet of one is the inlet of the next).

    d: Pipe inside diameter (in | mm)
    length: Segment length (ft | m)
    degf: Average flowing gas temperature (deg F | deg C)
    sg: Gas SG relative to air
    p1, p2: Upstream and downstream pressure (psia | barsa)
    qg: Gas rate (Mscf/d | sm3/d)
    method: 'COLEBROOK' (default; General Flow equation, Colebrook-White friction,
        64/Re below Re 2000), 'AGA' (General Flow equation, AGA NB-13 fully and
        partially turbulent transmission factors, the smaller governing),
        'WEYMOUTH', 'PANHANDLE_A', 'PANHANDLE_B' or 'IGT'
    efficiency: Pipeline efficiency E, multiplies the rate (default 1.0; Panhandle
        and Weymouth are usually run at 0.88-0.95)
    roughness: Absolute pipe roughness (in | mm), COLEBROOK and AGA. Default
        0.0007 in | 0.01778 mm (700 micro-inches, bare steel)
    drag_factor: AGA drag factor Df (default 0.96, average bend index, bare steel)
    dz: Outlet elevation minus inlet elevation (ft | m). Default 0
    co2, h2s, n2, h2: Mole fractions in the gas
    zmethod, cmethod: Z-factor and critical-property methods, as for gas.gas_z
    metric: If True, inputs and outputs in Eclipse METRIC units

    Returns dict with 'p1', 'p2' (psia | barsa), 'qg' (Mscf/d | sm3/d), 'pavg'
    (psia | barsa), 'z' at pavg, 'reynolds', 'transmission_factor' (F the rate
    implies in the General Flow equation; for the empirical equations, their
    equivalent F), 'friction_factor' (Darcy, 4 / F^2) and 's' (elevation parameter).
    """
    validate_choice(method, _METHODS, 'method')
    if roughness is None:
        roughness = _BARE_STEEL_ROUGHNESS * (IN_TO_MM if metric else 1.0)
    _check_inputs(p1, p2, qg, d, length, efficiency, roughness, drag_factor)
    p1, p2, qg, degf, d, length, roughness, dz = _to_field(
        metric, p1, p2, qg, degf, d, length, roughness, dz)
    for p in (p1, p2):
        if p is not None:
            validate_pe_inputs(p=p)
    validate_pe_inputs(degf=degf, sg=sg, co2=co2, h2s=h2s, n2=n2, h2=h2)

    pipe = _Pipe(d=d, l_mi=length / _FT_PER_MILE, tf=degf + degF2R, sg=sg,
                 method=method, e=efficiency, rough=roughness, df=drag_factor, dz=dz)
    props = _gas_props(degf, sg, co2, h2s, n2, h2, zmethod, cmethod)
    if qg is None:
        q, state = _rate(p1, p2, pipe, props)
        if q <= 0:
            raise ValueError(f"No forward flow: p1 = {p1:.6g} psia does not exceed "
                             f"p2 = {p2:.6g} psia plus the elevation head")
        qg = q / 1000.0
    else:
        # Bracket probes reach pressures far from the answer; their Z-method range
        # warnings are noise. The final state below warns if the answer itself is out.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            if p2 is None:
                p2 = _solve_p2(p1, qg * 1000.0, pipe, props)
            else:
                p1 = _solve_p1(p2, qg * 1000.0, pipe, props)
        state = _rate(p1, p2, pipe, props)[1]

    out = {'p1': p1, 'p2': p2, 'qg': qg, 'pavg': state['pavg'], 'z': state['z'],
           'reynolds': state['re'], 'transmission_factor': state['f_trans'],
           'friction_factor': 4.0 / state['f_trans'] ** 2, 's': state['s']}
    if metric:
        for k in ('p1', 'p2', 'pavg'):
            out[k] *= PSI_TO_BAR
        out['qg'] *= MSCF_TO_SM3
    return out
