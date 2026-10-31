"""
Gas well deliverability tests: Rawlins-Schellhardt back-pressure (C, n) and
Houpeurt laminar-inertial-turbulent (LIT, a and b) fits to multi-rate test data,
in pressure-squared or real-gas pseudopressure, with absolute open flow (AOF).

Rawlins & Schellhardt (1935), Back-Pressure Data on Natural-Gas Wells and Their
Application to Production Practices, USBM Monograph 7: q = C (pr^2 - pwf^2)^n.
Houpeurt (1959), Revue de l'IFP 14: dpsi = a q + b q^2.
Cullender (1955), Trans. AIME 204, 137: isochronal testing. Katz et al. (1959),
Handbook of Natural Gas Engineering: modified isochronal testing.
Test procedures and analysis as in Lee & Wattenbarger, Gas Reservoir Engineering
(SPE Textbook Vol. 5, 1996), Chapter 7.
"""

import math
import warnings
from dataclasses import dataclass
from typing import Optional

import numpy as np

from pyrestoolbox.constants import psc, BAR_TO_PSI, PSI_TO_BAR, MSCF_TO_SM3, SM3_TO_MSCF, degc_to_degf
from pyrestoolbox.shared_fns import validate_pe_inputs
from pyrestoolbox.validate import validate_choice

# Back-pressure exponent bounds: n = 1 is laminar (Darcy) flow, n = 0.5 fully
# turbulent; a fitted n outside them indicates non-stabilised or bad test data.
_N_MIN, _N_MAX = 0.5, 1.0


@dataclass
class DeliverabilityResult:
    """Fitted gas deliverability equation.

    method: 'CN' (q = C dp^n) or 'LIT' (dp = a q + b q^2)
    pseudo: True when dp is the pseudopressure difference, False for p^2
    c, n: back-pressure coefficient and exponent (CN only, else None)
    a, b: laminar and inertial-turbulent coefficients (LIT only, else None)
    pr: reservoir pressure the AOF refers to (psia | barsa)
    aof: absolute open flow, rate at pwf = 14.696 psia (Mscf/d | sm3/d)
    r_squared: fit of the points that set the slope (log-log for CN, dp/q vs q for LIT)
    test_type: 'flow_after_flow', 'isochronal' or 'modified_isochronal'
    metric: unit system of every dimensional field. dp is in psia^2 or psia^2/cP
        (barsa^2 or barsa^2/cP) and q in Mscf/d (sm3/d)
    """
    method: str
    pseudo: bool
    c: Optional[float]
    n: Optional[float]
    a: Optional[float]
    b: Optional[float]
    pr: float
    aof: float
    r_squared: float
    test_type: str
    metric: bool

    def rate(self, pwf, pr=None):
        """Gas rate (Mscf/d | sm3/d) from the fitted equation at pwf (psia | barsa).

        pr defaults to the reservoir pressure of the fit; pass a lower pr to move
        the curve with depletion (the coefficients are held fixed).
        """
        pr = self.pr if pr is None else pr
        scale = BAR_TO_PSI if self.metric else 1.0
        p_hi, p_lo = pr * scale, np.atleast_1d(np.asarray(pwf, dtype=float)) * scale
        if np.any(p_lo > p_hi) or np.any(p_lo <= 0):
            raise ValueError("pwf must be positive and not above pr")
        d = np.array([_dp(p_hi, p, self.pseudo, self._dp_args) for p in p_lo])
        if self.metric:
            d = d * PSI_TO_BAR ** 2
        if self.method == 'CN':
            q = self.c * d ** self.n
        else:
            q = np.array([_lit_rate(self.a, self.b, x) for x in d])
        return float(q[0]) if np.ndim(pwf) == 0 else q


def _r_squared(y, y_fit):
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 1.0 - np.sum((y - y_fit) ** 2) / ss_tot if ss_tot > 0 else 1.0


def _lit_rate(a, b, dp):
    """Positive root of b q^2 + a q - dp = 0."""
    if b == 0:
        return dp / a
    return (-a + math.sqrt(a * a + 4.0 * b * dp)) / (2.0 * b)


def gas_deliverability(pwf, qg, pr, degf=None, sg=0.65, method='LIT', pseudo=True,
                       pwf_stab=None, qg_stab=None, pr_stab=None,
                       co2=0, h2s=0, n2=0, h2=0, zmethod='DAK', cmethod='PMC',
                       metric=False):
    """Fit a gas well deliverability equation to multi-rate test data.

    Flow-after-flow test: give pwf and qg for each stabilised flow period and a scalar
    pr. Isochronal or modified isochronal test: pwf and qg are the equal-duration
    (transient) flow periods, which set the slope (n, or b); the extended stabilised
    flow period (pwf_stab, qg_stab) then sets the deliverability level (C, or a).
    For a modified isochronal test, pass pr as an array of the shut-in pressure
    before each flow period, and pr_stab as the average reservoir pressure for
    the stabilised point and the AOF.

    The pseudopressure form (pseudo=True, default) is valid at any pressure; the
    pressure-squared form assumes constant mu z, reasonable below about 2000 psia.

    pwf: Flowing bottomhole pressures (psia | barsa), one per flow period
    qg: Gas rates (Mscf/d | sm3/d), one per flow period
    pr: Reservoir (static) pressure (psia | barsa). Scalar, or an array of the
        preceding shut-in pressure per flow period (modified isochronal)
    degf: Reservoir temperature (deg F | deg C). Required when pseudo=True
    sg: Gas SG relative to air (default 0.65)
    method: 'LIT' (default, Houpeurt dp = a q + b q^2) or 'CN' (Rawlins-Schellhardt
        q = C dp^n)
    pseudo: True (default) for pseudopressure differences, False for p^2
    pwf_stab, qg_stab: Extended stabilised flow point for isochronal tests
        (psia | barsa, Mscf/d | sm3/d)
    pr_stab: Reservoir pressure for the stabilised point and the AOF (psia | barsa).
        Defaults to pr when pr is a scalar; required when pr is an array
    co2, h2s, n2, h2: Mole fractions in the gas
    zmethod, cmethod: Z-factor and critical-property methods, as for gas_dmp
    metric: If True, inputs and outputs in Eclipse METRIC units

    Returns DeliverabilityResult. Its rate(pwf) method evaluates the fitted equation.
    """
    method = validate_choice(str(method).upper(), ('LIT', 'CN'), 'method')
    pwf = np.atleast_1d(np.asarray(pwf, dtype=float))
    qg = np.atleast_1d(np.asarray(qg, dtype=float))
    pr_arr = np.asarray(pr, dtype=float)
    test_type = _test_type(pr_arr.ndim, pwf_stab, qg_stab, pr_stab)
    if pr_arr.ndim == 0:
        pr_stab = float(pr) if pr_stab is None else pr_stab
        pr_arr = np.full_like(pwf, float(pr_arr))
    if not (len(pwf) == len(qg) == len(pr_arr)):
        raise ValueError("pwf, qg and per-period pr must have the same length")
    if len(pwf) < 2:
        raise ValueError("Need at least two flow periods to fit a deliverability slope")
    if pseudo and degf is None:
        raise ValueError("degf is required for the pseudopressure form")
    stab = test_type != 'flow_after_flow'
    pwf, qg, pr_arr, pr_stab, pwf_stab, qg_stab, degf = _field_points(
        pwf, qg, pr_arr, pr_stab, pwf_stab, qg_stab, degf, metric)
    validate_pe_inputs(p=np.concatenate([pwf, pr_arr, [pr_stab]]), sg=sg,
                       co2=co2, h2s=h2s, n2=n2, h2=h2)
    if degf is not None:
        validate_pe_inputs(degf=degf)

    dp_args = (degf, sg, co2, h2s, n2, h2, zmethod, cmethod)
    d = np.array([_dp(hi, lo, pseudo, dp_args) for hi, lo in zip(pr_arr, pwf)])
    d_aof = _dp(pr_stab, psc, pseudo, dp_args)
    stab_pt = (_dp(pr_stab, pwf_stab, pseudo, dp_args), qg_stab) if stab else None

    if method == 'CN':
        c, n, r2 = _fit_cn(d, qg, stab_pt)
        a = b = None
        aof = c * d_aof ** n
    else:
        a, b, r2 = _fit_lit(d, qg, stab_pt)
        c = n = None
        aof = _lit_rate(a, b, d_aof)

    if metric:
        c, a, b = _coeffs_to_metric(c, n, a, b)
        pr_stab, aof = pr_stab * PSI_TO_BAR, aof * MSCF_TO_SM3
    result = DeliverabilityResult(method=method, pseudo=pseudo, c=_opt_float(c), n=_opt_float(n),
                                  a=_opt_float(a), b=_opt_float(b), pr=float(pr_stab),
                                  aof=float(aof), r_squared=float(r2), test_type=test_type,
                                  metric=metric)
    # Fluid and methods for rate(); a plain attribute, so not serialised as a field
    result._dp_args = dp_args
    return result


def _field_points(pwf, qg, pr_arr, pr_stab, pwf_stab, qg_stab, degf, metric):
    """Convert the test points to field units and check rates and drawdowns."""
    stab = pwf_stab is not None
    if metric:
        pwf, pr_arr, qg = pwf * BAR_TO_PSI, pr_arr * BAR_TO_PSI, qg * SM3_TO_MSCF
        pr_stab = pr_stab * BAR_TO_PSI
        if stab:
            pwf_stab, qg_stab = pwf_stab * BAR_TO_PSI, qg_stab * SM3_TO_MSCF
        degf = None if degf is None else degc_to_degf(degf)
    q_all, pwf_all, pref_all = qg, pwf, pr_arr
    if stab:
        q_all, pwf_all = np.append(qg, qg_stab), np.append(pwf, pwf_stab)
        pref_all = np.append(pr_arr, pr_stab)
    if np.any(q_all <= 0):
        raise ValueError("Gas rates must be positive")
    if np.any(pwf_all >= pref_all):
        raise ValueError("Each pwf must be below its reservoir (shut-in) pressure")
    return pwf, qg, pr_arr, pr_stab, pwf_stab, qg_stab, degf


def _test_type(pr_ndim, pwf_stab, qg_stab, pr_stab):
    """Classify the test from which inputs were given, and check they are complete."""
    if (pwf_stab is None) != (qg_stab is None):
        raise ValueError("Give both pwf_stab and qg_stab, or neither")
    stab = pwf_stab is not None
    if pr_ndim == 0:
        return 'isochronal' if stab else 'flow_after_flow'
    if pr_stab is None:
        raise ValueError("pr_stab is required when pr is given per flow period")
    if not stab:
        raise ValueError("A modified isochronal test needs the stabilised point "
                         "(pwf_stab, qg_stab)")
    return 'modified_isochronal'


def _dp(p_hi, p_lo, pseudo, dp_args):
    """Drawdown term: pseudopressure difference (psia^2/cP) or p^2 difference (psia^2)."""
    if not pseudo:
        return p_hi ** 2 - p_lo ** 2
    from pyrestoolbox.gas.gas import gas_dmp
    degf, sg, co2, h2s, n2, h2, zmethod, cmethod = dp_args
    return gas_dmp(p_lo, p_hi, degf, sg, zmethod=zmethod, cmethod=cmethod,
                   co2=co2, h2s=h2s, n2=n2, h2=h2)


def _fit_cn(d, qg, stab_pt):
    """Back-pressure fit: n from log q vs log dp; C from the stabilised point if given."""
    n, log_c = np.polyfit(np.log(d), np.log(qg), 1)
    r2 = _r_squared(np.log(qg), log_c + n * np.log(d))
    c = math.exp(log_c) if stab_pt is None else stab_pt[1] / stab_pt[0] ** n
    if not _N_MIN <= n <= _N_MAX:
        warnings.warn(f"Back-pressure exponent n = {n:.3f} is outside [0.5, 1]; "
                      f"the test data are likely not stabilised or inconsistent")
    return c, n, r2


def _fit_lit(d, qg, stab_pt):
    """LIT fit: b and a from dp/q vs q; a from the stabilised point if given."""
    b, a = np.polyfit(qg, d / qg, 1)
    r2 = _r_squared(d / qg, a + b * qg)
    if stab_pt is not None:
        a = (stab_pt[0] - b * stab_pt[1] ** 2) / stab_pt[1]
    if a <= 0:
        raise ValueError(f"Fitted laminar coefficient a = {a:.4g} is not positive; "
                         f"the test data cannot be represented by the LIT equation")
    if b < 0:
        warnings.warn(f"Fitted turbulence coefficient b = {b:.4g} is negative; "
                      f"the test data are likely not stabilised or inconsistent")
    return a, b, r2


def _coeffs_to_metric(c, n, a, b):
    """Field-unit coefficients to metric: dp in barsa^2 (or barsa^2/cP), q in sm3/d."""
    dp_conv = PSI_TO_BAR ** 2
    if c is not None:
        c = c * MSCF_TO_SM3 / dp_conv ** n
    if a is not None:
        a, b = a * dp_conv / MSCF_TO_SM3, b * dp_conv / MSCF_TO_SM3 ** 2
    return c, a, b


def _opt_float(x):
    return None if x is None else float(x)
