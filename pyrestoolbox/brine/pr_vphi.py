"""V2inf of dissolved gases from the Soreide-Whitson modified PR EOS + VSHIFT.

An alternative to `plyasunov_model.V2_inf`, built 2026-07-25. It is NOT a
replacement: Plyasunov remains the reference and covers 8 gases to 573 K. This
route exists because it is far more parsimonious and because it repairs the one
V_phi defect the project has pinned (the H2S temperature trend).

THE EXACT RELATION. For any pressure-explicit EOS the partial molar volume is

    Vbar_2 = -(dP/dn2)_{T,V,n1} / (dP/dV)_{T,n}

and at x2 -> 0 that is V2inf. No flash and no root selection for the solute: the
volume root is the water-rich LIQUID root, and the gas enters only through b2
and the cross term a12 = sqrt(a1 a2)(1 - kij). The solute never has a root of
its own, so there is nothing to "force" into a liquid state.

WHAT TO SAY WHEN CHALLENGED ON PR'S POOR WATER PROPERTIES. Two separate points,
and the first is a scope point rather than a numerical one:

  1. PR's water volumetrics never reach the delivered answer. This module
     supplies V_phi ONLY. The solvent density rho1 in the Garcia mixing rule
     comes from Spivey (`brine_properties.rho_brine`), never from this EOS. So
     PR's water error is confined to an intermediate quantity that is then
     calibrated against measured V_phi.

  2. A constant shift absorbs an OFFSET but cannot fix a SHAPE, so the question
     that matters is whether V_phi's T and P dependence survives. The answer is
     empirical: one constant s per gas reproduces Hnedkovsky densimetry to 0.5%
     mean absolute error for H2S, 0.9% for CO2 and 1.4% for CH4 across
     298-473 K. Lead with that measurement, not with a story about errors
     cancelling.

Do NOT justify the temperature cap by citing PR's 15-17% water density error -
that error is a FLAT systematic factor (density ratio PR/IF97 is 0.847 at 298 K
and 0.820 at 473 K) and is precisely what the shift absorbs. The cap exists for
two other reasons, given at `T_MAX`.

VOLUME SHIFT. Peneloux translation v = v_EOS - sum(n_i c_i) gives, exactly,

    V2inf(shifted) = V2inf(EOS) - c2,      c2 = s2 * b2

so a constant s2 is a pure offset with no temperature or pressure shape. The
shift is stored dimensionless (s = c/b, the standard VSHIFT convention) and
b2 is the gas's PR co-volume.

ONE PARAMETER PER GAS, NOT TWO. A held-out test (`fit_pr_vshift.py`) showed
s(T) = s0 + s1/T overfits: on the 21 Murphy & Gaines H2S points it scored 1.3%
mean absolute error against 0.6% for a constant s. Constant s is used.

CALIBRATION BASIS: DIRECT VOLUMETRIC MEASUREMENT ONLY. See `fit_pr_vshift.py`
for the set and the exclusions. Indirect determinations (solubility-versus-
pressure, electrochemical) and V_phi inverted from brine density are excluded
from the fit, the latter because they scatter five to seven times wider.

PROVENANCE. Developed in ~/projects/brine_props (`code/pr_vphi_model.py`,
calibrated by `code/fit_pr_vshift.py`) and ported here 2026-07-25. The S&W alpha,
the aqueous kij correlations and the critical properties are taken from this
package's own `_lib_vle_engine`, not re-copied.
"""
import numpy as np

from . import _lib_vle_engine as _SW
from ..plyasunov.iapws_if97 import p_sat_if97

R = 8.31446261815324          # cm3.MPa/(mol.K), so volumes come out in cm3/mol

# PR constants (Peng-Robinson 1976)
_OMEGA_A = 0.45723553
_OMEGA_B = 0.07779607

# Validity envelope. Two reasons for the 473.15 K ceiling, neither of them the
# 15-17% water density error (that is a flat offset and the shift absorbs it):
#   (a) the calibration data stops there. Hnedkovsky's top point is 473.15 K, so
#       anything above is unvalidated extrapolation of a fitted shift.
#   (b) V2inf diverges at the SOLVENT critical point, and a cubic EOS gets the
#       near-critical behaviour wrong in a way no constant offset can repair.
#       Measured on this implementation, the PR/IF97 compressibility ratio for
#       water is 0.27 at 298 K, 0.84 at 473 K, crosses 1.0 near 520 K and reaches
#       1.32 by 613 K - i.e. a growing SHAPE error, not an offset. It shows up
#       downstream as H2S drifting 9% by 523 K and 36% by 573 K.
T_MIN, T_MAX = 273.15, 623.15
T_CALIBRATED_MAX = 473.15   # top of the calibration data (Hnedkovsky)
T_VOUCHED_MAX = 450.0       # ~350 degF; the only one of the three that is an
                            # accuracy claim. T_MAX is where the arithmetic
                            # stops (IF97 Region 1); T_CALIBRATED_MAX is where
                            # the fitted shift stops being fitted.
P_MAX = 100.0

SUPPORTED = ('CH4', 'CO2', 'H2S', 'N2', 'H2', 'C2H6', 'C3H8', 'NC4H10')

# This library spells butane 'NC4H10'; the S&W component tables spell it
# 'nC4H10'. That mismatch is why butane looked absent from the component
# set until 2026-07-25 - it was there all along under the other spelling.
_SW_NAME = {'NC4H10': 'nC4H10'}

# WATER ALPHA INSENSITIVITY, measured 2026-07-25 and worth knowing before
# porting. The S&W framework uses `alpha_water_soreide` for the sw_original and
# dropin frameworks but Mathias-Copeman (`alpha_water_mc3`) otherwise, and the
# ResToolbox3 Dart path uses MC3. Swapping the water alpha moves V2inf_raw by
# only 0.00 to 0.05 cm3/mol over 298-473 K, i.e. delta-s of 0.0002 to 0.0018,
# under 2% of the shift values and negligible against the 0.5-1.4% fit error.
# So ONE VSHIFT set is valid for both water alphas and the ports do not need
# separate calibrations. (The alpha enters only a1, and its effect largely
# cancels in the ratio of derivatives.)

# Dimensionless volume shifts s = c2/b2, fitted to direct volumetric
# measurements only by `fit_pr_vshift.py`. Regenerate with that script if any
# calibration point changes; the values below are pinned in validation.py.
VSHIFT = {
    'CH4':  -0.109632,
    'CO2':  -0.037913,
    'H2S':  -0.078975,
    'N2':   -0.176510,
    'H2':   -0.177625,
    'C2H6': -0.073142,
    # C3H8 added 2026-07-25. NOT fitted to a densimetric data set, because none
    # exists for propane in water. Set from the only two direct 298 K
    # determinations available, Moore (1982) 70.7 and Zhou & Battino (2001)
    # 75.0 cm3/mol, which imply s = -0.074763 and -0.151164; the adopted value
    # is their mean, giving V_phi = 72.85 at 298 K. Those two disagree by 6.1%,
    # which is the honest uncertainty on this gas and is far wider than for any
    # other. It replaces a fallback that sat BELOW both of them (66.99).
    'C3H8': -0.112963,
    # NC4H10 added 2026-07-25 from Moore (1982) Table I, 76.6 +/- 0.1 cm3/mol
    # (2 runs; his stated overall imprecision is +/-1.5). Two oddities that
    # should not be smoothed over: the shift is POSITIVE, alone among the
    # eight, and it breaks Moore's own homologous series (his CH2 increments
    # run +18.4, +17.8, then +5.9). Adopted because it is the only direct
    # measurement of this quantity, but it is a single lab and two runs, and
    # C4 sits outside the five-gas validated scope.
    'NC4H10': +0.110924,
}

# Gases whose shift rests on 298 K data alone, so their temperature behaviour is
# the EOS's unaided prediction and is NOT calibrated. Reported, not hidden.
UNCALIBRATED_IN_T = ('N2', 'H2', 'C2H6', 'C3H8', 'NC4H10')


def _ab(species, T, m_nacl=0.0):
    """PR a(T) and b for one species, cm6.MPa/mol2 and cm3/mol."""
    c = _SW.COMPONENTS[_SW_NAME.get(species, species)]
    Pc = c.Pc / 1e6
    Tr = T / c.Tc
    if species == 'H2O':
        alpha = _SW.alpha_water_soreide(Tr, m_nacl)
    else:
        alpha = _SW.alpha_standard_pr(Tr, c.omega)
    return (_OMEGA_A * R ** 2 * c.Tc ** 2 / Pc * alpha,
            _OMEGA_B * R * c.Tc / Pc)


def b_covolume(gas):
    """PR co-volume of the gas, cm3/mol - the normaliser for the VSHIFT."""
    return _ab(gas, 300.0)[1]


def _P_mix(T, V, n1, n2, gas, kij, m_nacl):
    """PR pressure of the binary water + gas mixture, MPa. V in cm3, n in mol."""
    a1, b1 = _ab('H2O', T, m_nacl)
    a2, b2 = _ab(gas, T)
    n = n1 + n2
    x1, x2 = n1 / n, n2 / n
    a12 = np.sqrt(a1 * a2) * (1.0 - kij)
    am = x1 * x1 * a1 + 2.0 * x1 * x2 * a12 + x2 * x2 * a2
    bm = x1 * b1 + x2 * b2
    return (n * R * T / (V - n * bm)
            - n * n * am / (V * V + 2.0 * n * bm * V - n * n * bm * bm))


def v_water_liquid(T, P, m_nacl=0.0):
    """Smallest real PR root for water/brine, cm3/mol.

    Above water's critical temperature the cubic has one real root; below it the
    smallest root is the liquid. Taking the minimum root is continuous across
    both cases, which is what the derivative below needs.
    """
    a, b = _ab('H2O', T, m_nacl)
    A = a * P / (R * T) ** 2
    B = b * P / (R * T)
    roots = np.roots([1.0, -(1.0 - B), A - 3.0 * B ** 2 - 2.0 * B,
                      -(A * B - B ** 2 - B ** 3)])
    real = roots[np.abs(roots.imag) < 1e-9].real
    real = real[real > B]
    if real.size == 0:
        raise ValueError(f'No admissible PR liquid root for water at {T} K, {P} MPa')
    return float(real.min()) * R * T / P


def _check(gas, T, P):
    if gas not in SUPPORTED:
        raise ValueError(
            f'{gas!r} has no fitted VSHIFT. Supported: {SUPPORTED}. Use '
            f'plyasunov_model.V2_inf for the other gases.')
    if not (T_MIN <= T <= T_MAX):
        raise ValueError(
            f'T = {T} K is outside the calibrated range {T_MIN}-{T_MAX} K. '
            f'Above {T_MAX} K the fitted shift is unvalidated (the calibration '
            f'data stops there) and the near-critical shape error grows; use '
            f'plyasunov_model.V2_inf instead.')
    if not (0.0 < P <= P_MAX):
        raise ValueError(f'P = {P} MPa is outside 0 to {P_MAX} MPa.')
    psat = float(p_sat_if97(T))
    if P < psat:
        raise ValueError(
            f'P = {P} MPa is below the water saturation pressure {psat:.4f} MPa '
            f'at {T} K, so the solvent is vapour and V2inf is undefined.')


def V2_inf_raw(gas, T, P, m_nacl=0.0):
    """Unshifted V2inf from the exact EOS relation, cm3/mol. No range guard."""
    # The two derivatives of V2 = -(dP/dn2)/(dP/dV) are evaluated ANALYTICALLY
    # at n1 = 1, n2 = 0, V = v_w (closed forms: manuscript Appendix B).
    # Switched from finite differences 2026-07-31; the difference evaluation
    # agrees to 4e-5 cm3/mol, its truncation error.
    v = v_water_liquid(T, P, m_nacl)
    kij = _SW.get_kij_aq(_SW_NAME.get(gas, gas), T, m_nacl)
    a1, b1 = _ab('H2O', T, m_nacl)
    a2, b2 = _ab(gas, T)
    a12 = np.sqrt(a1 * a2) * (1.0 - kij)
    D = v * v + 2.0 * b1 * v - b1 * b1
    dPdn2 = (R * T / (v - b1) + R * T * b2 / (v - b1) ** 2
             - (2.0 * a12 * D - 2.0 * a1 * b2 * (v - b1)) / D ** 2)
    dPdV = -R * T / (v - b1) ** 2 + 2.0 * a1 * (v + b1) / D ** 2
    return -dPdn2 / dPdV


def V2_inf(gas, T, P, m_nacl=0.0, s=None):
    """Apparent molar volume of the dissolved gas at infinite dilution, cm3/mol.

    Args:
        gas: one of SUPPORTED
        T: K, within T_MIN..T_MAX and above the water saturation pressure
        P: MPa, up to P_MAX
        m_nacl: NaCl molality. Enters the S&W water alpha and kij only; there is
            NO salinity term on the volume shift. The salt effect on V_phi
            is applied one level up (vphi_route, as a relative fraction since
            2026-07-30), never here, so there is no double counting. Leave at
            0 for the calibrated behaviour.
        s: override the dimensionless volume shift (for refitting).

    Returns:
        V2inf in cm3/mol.
    """
    _check(gas, T, P)
    shift = VSHIFT[gas] if s is None else s
    return V2_inf_raw(gas, T, P, m_nacl) - shift * b_covolume(gas)


# `plyasunov_model` exposes V_phi as an alias for V2_inf; mirror that so the two
# modules are drop-in comparable in scripts.
V_phi = V2_inf


def _demo():
    print('S&W PR + VSHIFT, V2inf in cm3/mol (Plyasunov in brackets)')
    from ..plyasunov import V2_inf as V_ply
    print(f"\n  {'gas':<6}{'s':>10}{'b':>8}  " +
          ''.join(f'{t:>18.0f} K' for t in (298, 348, 398, 448)))
    for g in SUPPORTED:
        cells = []
        for T in (298.15, 348.15, 398.15, 448.15):
            cells.append(f'{V2_inf(g, T, 20.0):7.2f} '
                         f'[{float(V_ply(g, T, 20.0)):6.2f}]')
        flag = ' *' if g in UNCALIBRATED_IN_T else ''
        print(f'  {g:<6}{VSHIFT[g]:10.4f}{b_covolume(g):8.2f}  '
              + ''.join(f'{c:>20}' for c in cells) + flag)
    print('\n  * temperature behaviour is the EOS unaided: these gases have no '
          'T-resolved\n    calibration data, only a 298 K cluster.')


__all__ = ['V2_inf', 'V_phi', 'V2_inf_raw', 'b_covolume', 'v_water_liquid',
           'VSHIFT', 'SUPPORTED', 'UNCALIBRATED_IN_T', 'T_MIN', 'T_MAX', 'P_MAX']
