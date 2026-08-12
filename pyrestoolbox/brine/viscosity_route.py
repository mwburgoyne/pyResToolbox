"""Dispatcher for the BRINE VISCOSITY leg.

    mu_brine(T, P, salts) = mu_water(T, P) * salt ratio * [pressure factor]

then the per-gas correction in `garcia_mixing` multiplies that. This module
covers the gas-free brine only.

DEFAULT since 2026-07-26: **`'jones_dole'` with the IAPWS-2008 water leg and the
Kestin pressure factor.** Mark's call, taken on capability and reservoir-
temperature accuracy after the salt term was measured against Kestin.

    'jones_dole'  Appelo & Parkhurst's ion-additive modified Jones-Dole model
                  as implemented in PHREEQC. **MULTI-SALT** over 15 ions.
                  Reproduces PHREEQC to 0.0015% mean over 145 cases.
    'mao_duan'    Mao-Duan (2009), NaCl only. The pre-2026-07-26 default,
                  retained for reproducibility. It RAISES on a non-NaCl brine
                  rather than silently pretending the salt is NaCl.

**SALINITY WITHOUT A SPECIES IS NaCl.** `brine_viscosity(T, P, m=3.0)` and
`salts={'NaCl': 3.0}` are the same call, matching the density side and the
`S`-as-weight-fraction-NaCl convention used throughout this project.

WHY JONES-DOLE, WITH THE NUMBERS. Against Kestin's measured NaCl viscosity
(`Papers/50`) the two salt terms TIE - 0.415% mean for Jones-Dole, 0.443% for
Mao-Duan, both inside Kestin's own +/-0.5%. Jones-Dole wins at reservoir
temperature (0.408% against 0.754% over 125-150 degC) and loses in the cold
concentrated corner. What decided it is capability, measured on Kestin's KCl
companion paper (`Papers/51`):

| KCl, 25-150 degC, 0.1-35 MPa, 0.5-5 molal | mean | max |
|---|---|---|
| ion-additive Jones-Dole (K+ from pitzer.dat) | **0.77%** | 4.6% |
| a NaCl-only leg, forced to treat KCl as NaCl | **15.9%** | 51.2% |

THE PRESSURE FACTOR. Neither salt model carries any pressure dependence, and the
measured one is not small: Kestin's NaCl ratio moves +3.2% from 0.1 to 35 MPa at
20 degC / 4 molal, the dependence's sign reversing above ~91 degC (molality-
independent; KCl crosses at ~68 degC). The shipped ratio is therefore
multiplied by Kestin's factor normalised to 0.1 MPa, evaluated at the brine's
IONIC STRENGTH as an NaCl-equivalent molality. That is a first-order,
electrolyte-generic correction rather than an NaCl-specific one: NaCl and KCl
agree on its sign everywhere and on its size to about 0.6 pp (at 25 degC / 4
molal, +2.74% for NaCl against +2.09% for KCl). It improves both -

| vs the measured ratio | mean, as-is | mean, grafted | max, as-is | max, grafted |
|---|---|---|---|---|
| NaCl (Papers/50) | 0.365% | **0.299%** | 2.57% | **0.78%** |
| KCl (Papers/51) | 0.973% | **0.767%** | 6.57% | **4.62%** |

and it is CLAMPED, not extrapolated, outside Kestin's 20-150 degC / 0.1-35 MPa /
6 molal calibration box.

WHAT IS NOT VALIDATED: ion MIXING. Every measured brine scored here is a single
salt, because no measured mixed-brine viscosity was found. Mixed-brine viscosity
rests on reproducing PHREEQC, which is verification of the implementation, not
validation of the model.
"""

DEFAULT_ROUTE = 'jones_dole'
ROUTES = ('jones_dole', 'mao_duan')

# Kestin's calibration box for the pressure factor. Outside it the factor is
# held at the boundary rather than extrapolated.
_P_FACTOR_T_C = (20.0, 150.0)
_P_FACTOR_P_MPA = (0.1, 35.0)
_P_FACTOR_M_MAX = 6.0


def _composition(composition=None, salts=None, m=None, S=None):
    """Resolve any of the four salinity inputs into {ion: molality}.

    A bare molality or weight fraction means NaCl - the same default the
    density side uses.
    """
    given = [x for x in (composition, salts, m, S) if x is not None]
    if len(given) != 1:
        raise ValueError(
            "supply exactly one of composition={ion: molality}, "
            "salts={salt: molality}, m=<NaCl molality> or S=<NaCl weight "
            "fraction>")

    if composition is not None:
        return dict(composition)
    if salts is not None:
        from .appelo_volumes import ions_from_salts
        return ions_from_salts(salts)
    if S is not None:
        from .appelo_volumes import molality_from_salinity
        m = molality_from_salinity(S)
    if m <= 0.0:
        return {}
    return {'Na+': float(m), 'Cl-': float(m)}


def pressure_factor(T_K, P_MPa, composition):
    """
    Kestin's measured pressure dependence of the salt ratio, normalised to
    0.1 MPa, evaluated at the brine's ionic strength as NaCl-equivalent
    molality. Returns 1.0 for pure water.
    """
    from .appelo_volumes import ionic_strength
    from . import kestin_nacl_viscosity as kestin

    I = ionic_strength(composition)
    if I <= 0.0:
        return 1.0

    t = min(max(T_K - 273.15, _P_FACTOR_T_C[0]), _P_FACTOR_T_C[1])
    p = min(max(P_MPa, _P_FACTOR_P_MPA[0]), _P_FACTOR_P_MPA[1])
    m = min(I, _P_FACTOR_M_MAX)
    return kestin.salt_ratio(t, p, m) / kestin.salt_ratio(t, 0.1, m)


def salt_ratio(T_K, P_MPa, composition=None, salts=None, m=None, S=None,
               route=None, pressure_term=True):
    """mu(brine)/mu(water) at the same T and P."""
    comp = _composition(composition, salts, m, S)
    route = route or DEFAULT_ROUTE
    if route not in ROUTES:
        raise ValueError(f"route must be one of {ROUTES}, got {route!r}")
    if not comp:
        return 1.0

    if route == 'jones_dole':
        from .jones_dole_viscosity import viscosity_ratio
        ratio = viscosity_ratio(T_K, P_MPa, comp)
    else:
        non_nacl = {k: v for k, v in comp.items()
                    if k not in ('Na+', 'Cl-') and v > 0}
        if non_nacl:
            raise ValueError(
                f"route='mao_duan' is NaCl only and cannot represent "
                f"{sorted(non_nacl)}. Treating them as NaCl costs 15.9% mean "
                f"on KCl (Papers/51); use the default 'jones_dole' route.")
        from ._lib_maoduan_viscosity import maoduan_ratio
        ratio = maoduan_ratio(T_K, comp.get('Na+', 0.0))

    if pressure_term:
        ratio *= pressure_factor(T_K, P_MPa, comp)
    return ratio


def mu_water(T_K, P_MPa, route=None):
    """Pure-water viscosity in mPa s (== cP).

    IAPWS-2008 for the default route - the international standard, and the
    correlation PHREEQC itself uses, which the Jones-Dole electrostatic term
    scales with. Mao-Duan's own water fit for `route='mao_duan'`.
    """
    route = route or DEFAULT_ROUTE
    if route == 'jones_dole':
        from .jones_dole_viscosity import mu_water as _mu
        return _mu(T_K, P_MPa)
    from ._lib_maoduan_viscosity import mu_water_maoduan
    return mu_water_maoduan(T_K, P_MPa)


def brine_viscosity(T_K, P_MPa, composition=None, salts=None, m=None, S=None,
                    route=None, pressure_term=True):
    """
    Gas-free brine viscosity in mPa s (== cP).

    >>> brine_viscosity(373.15, 30.0, m=3.0)                    # NaCl assumed
    >>> brine_viscosity(373.15, 30.0, S=0.15)                   # weight fraction
    >>> brine_viscosity(373.15, 30.0, salts={'NaCl': 2.0, 'CaCl2': 1.0})
    """
    route = route or DEFAULT_ROUTE
    return mu_water(T_K, P_MPa, route) * salt_ratio(
        T_K, P_MPa, composition, salts, m, S, route, pressure_term)
