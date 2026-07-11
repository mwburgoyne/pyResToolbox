#!/usr/bin/env python3
"""Tests for the DCA (Decline Curve Analysis) module."""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
from pyrestoolbox import dca


# ======================== arps_rate tests ========================

def test_arps_rate_exponential():
    """Exponential decline: q = qi * exp(-di*t)"""
    q = dca.arps_rate(1000, 0.1, 0, 10)
    assert abs(q - 1000 * np.exp(-1.0)) < 0.01


def test_arps_rate_harmonic():
    """Harmonic decline: q = qi / (1 + di*t)"""
    q = dca.arps_rate(1000, 0.1, 1.0, 10)
    assert abs(q - 1000 / 2.0) < 0.01


def test_arps_rate_hyperbolic():
    """Hyperbolic decline: q = qi / (1 + b*di*t)^(1/b)"""
    q = dca.arps_rate(1000, 0.1, 0.5, 10)
    expected = 1000 / (1.0 + 0.5 * 0.1 * 10) ** (1.0 / 0.5)
    assert abs(q - expected) < 0.01


def test_arps_rate_vectorized():
    """Test with array input."""
    t = [0, 5, 10]
    q = dca.arps_rate(1000, 0.1, 0, t)
    assert isinstance(q, np.ndarray)
    assert len(q) == 3
    assert abs(q[0] - 1000) < 0.01


def test_arps_rate_bad_di():
    try:
        dca.arps_rate(1000, -0.1, 0, 10)
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_arps_rate_bad_b():
    try:
        dca.arps_rate(1000, 0.1, -0.5, 10)
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_arps_rate_super_hyperbolic():
    """b > 1 (transient flow) is permitted: q = qi / (1 + b*di*t)^(1/b)"""
    q = dca.arps_rate(1000, 0.1, 2.0, 10)
    expected = 1000 / (1.0 + 2.0 * 0.1 * 10) ** 0.5
    assert abs(q - expected) < 0.01


# ======================== arps_cum tests ========================

def test_arps_cum_exponential():
    """Exponential cumulative: Qcum = qi/di * (1 - exp(-di*t))"""
    Qcum = dca.arps_cum(1000, 0.1, 0, 10)
    expected = (1000 / 0.1) * (1.0 - np.exp(-1.0))
    assert abs(Qcum - expected) < 0.1


def test_arps_cum_harmonic():
    """Harmonic cumulative: Qcum = qi/di * ln(1 + di*t)"""
    Qcum = dca.arps_cum(1000, 0.1, 1.0, 10)
    expected = (1000 / 0.1) * np.log(2.0)
    assert abs(Qcum - expected) < 0.1


def test_arps_cum_hyperbolic():
    """Hyperbolic cumulative should be positive and between exp and harmonic."""
    Qcum_exp = dca.arps_cum(1000, 0.1, 0, 10)
    Qcum_hyp = dca.arps_cum(1000, 0.1, 0.5, 10)
    Qcum_har = dca.arps_cum(1000, 0.1, 1.0, 10)
    assert Qcum_hyp > 0
    assert Qcum_exp < Qcum_hyp < Qcum_har


def test_arps_cum_at_zero():
    """Cumulative at t=0 should be 0."""
    assert abs(dca.arps_cum(1000, 0.1, 0, 0)) < 1e-10
    assert abs(dca.arps_cum(1000, 0.1, 0.5, 0)) < 1e-10
    assert abs(dca.arps_cum(1000, 0.1, 1.0, 0)) < 1e-10


# ======================== duong tests ========================

def test_duong_rate_basic():
    q = dca.duong_rate(1000, 1.5, 1.2, 1.0)
    # At t=1: q = qi * 1^(-m) * exp(a/(1-m)*(1-1)) = qi
    assert abs(q - 1000) < 0.01


def test_duong_rate_declining():
    """Rate should decrease with time for suitable parameters."""
    q1 = dca.duong_rate(1000, 1.5, 1.2, 10.0)
    q50 = dca.duong_rate(1000, 1.5, 1.2, 50.0)
    assert q50 < q1


def test_duong_rate_bad_a():
    try:
        dca.duong_rate(1000, -1.0, 1.2, 1.0)
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_duong_rate_bad_m():
    try:
        dca.duong_rate(1000, 1.5, 0.5, 1.0)
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_duong_cum_positive():
    Qcum = dca.duong_cum(1000, 1.5, 1.2, 10.0)
    assert Qcum > 0


# ======================== eur tests ========================

def test_eur_exponential():
    """EUR with exponential decline."""
    e = dca.eur(1000, 0.1, 0, 10)
    # Total recovery when rate drops to 10
    assert e > 0
    # Should be less than qi/di (infinite time recovery)
    assert e < 1000 / 0.1


def test_eur_harmonic():
    e = dca.eur(1000, 0.1, 1.0, 10)
    assert e > 0


def test_eur_bad_qmin():
    try:
        dca.eur(1000, 0.1, 0, 2000)  # q_min > qi
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_eur_super_hyperbolic_finite():
    """EUR with b > 1 is finite for q_min > 0."""
    e = dca.eur(1000, 0.005, 2.0, 10)
    t_end = ((1000 / 10) ** 2.0 - 1.0) / (2.0 * 0.005)
    assert abs(e - dca.arps_cum(1000, 0.005, 2.0, t_end)) < 1e-6


# ======================== modified hyperbolic tests ========================

# Shared MH reference case: qi=1000, di=0.005/day nominal, b=2, dterm=0.0005/day
# Switch when D(t) = dterm: t_sw = (1/dterm - 1/di)/b = 900 days
_MH_CASE = dict(qi=1000.0, di=0.005, b=2.0, dterm=0.0005)
_MH_T_SW = 900.0


def test_mh_rate_matches_arps_before_switch():
    c = _MH_CASE
    t = 400.0
    q = dca.mh_rate(c['qi'], c['di'], t, b=c['b'], dterm=c['dterm'])
    assert abs(q - dca.arps_rate(c['qi'], c['di'], c['b'], t)) < 1e-9


def test_mh_rate_exponential_after_switch():
    c = _MH_CASE
    q_sw = dca.arps_rate(c['qi'], c['di'], c['b'], _MH_T_SW)
    t = 2000.0
    q = dca.mh_rate(c['qi'], c['di'], t, b=c['b'], dterm=c['dterm'])
    expected = q_sw * np.exp(-c['dterm'] * (t - _MH_T_SW))
    assert abs(q - expected) < 1e-9


def test_mh_rate_continuous_at_switch():
    c = _MH_CASE
    q_lo = dca.mh_rate(c['qi'], c['di'], _MH_T_SW - 1e-6, b=c['b'], dterm=c['dterm'])
    q_hi = dca.mh_rate(c['qi'], c['di'], _MH_T_SW + 1e-6, b=c['b'], dterm=c['dterm'])
    assert abs(q_lo - q_hi) < 1e-5


def test_mh_rate_default_b_is_two():
    c = _MH_CASE
    q = dca.mh_rate(c['qi'], c['di'], 400.0, dterm=c['dterm'])
    assert abs(q - dca.mh_rate(c['qi'], c['di'], 400.0, b=2.0, dterm=c['dterm'])) < 1e-12


def test_mh_rate_no_dterm_is_single_segment():
    c = _MH_CASE
    q = dca.mh_rate(c['qi'], c['di'], 5000.0, b=c['b'])
    assert abs(q - dca.arps_rate(c['qi'], c['di'], c['b'], 5000.0)) < 1e-9


def test_mh_rate_vectorized():
    c = _MH_CASE
    q = dca.mh_rate(c['qi'], c['di'], [100.0, _MH_T_SW, 5000.0], b=c['b'], dterm=c['dterm'])
    assert isinstance(q, np.ndarray)
    assert len(q) == 3
    assert np.all(np.diff(q) < 0)


def test_mh_cum_matches_numeric_integral():
    c = _MH_CASE
    t_end = 3.0 * _MH_T_SW
    t = np.linspace(1e-9, t_end, 200001)
    q = dca.mh_rate(c['qi'], c['di'], t, b=c['b'], dterm=c['dterm'])
    numeric = np.trapezoid(q, t)
    analytic = dca.mh_cum(c['qi'], c['di'], t_end, b=c['b'], dterm=c['dterm'])
    assert abs(numeric / analytic - 1.0) < 1e-6


def test_mh_cum_continuous_at_switch():
    c = _MH_CASE
    n_lo = dca.mh_cum(c['qi'], c['di'], _MH_T_SW - 1e-6, b=c['b'], dterm=c['dterm'])
    n_hi = dca.mh_cum(c['qi'], c['di'], _MH_T_SW + 1e-6, b=c['b'], dterm=c['dterm'])
    assert abs(n_lo - n_hi) < 1e-2


def test_mh_eur_exponential_segment():
    """q_min below the switch rate: EUR = N_sw + (q_sw - q_min)/dterm."""
    c = _MH_CASE
    q_sw = dca.arps_rate(c['qi'], c['di'], c['b'], _MH_T_SW)
    n_sw = dca.arps_cum(c['qi'], c['di'], c['b'], _MH_T_SW)
    e = dca.mh_eur(c['qi'], c['di'], 5.0, b=c['b'], dterm=c['dterm'])
    assert abs(e - (n_sw + (q_sw - 5.0) / c['dterm'])) < 1e-6


def test_mh_eur_hyperbolic_segment():
    """q_min above the switch rate: identical to single-segment eur()."""
    c = _MH_CASE
    e = dca.mh_eur(c['qi'], c['di'], 600.0, b=c['b'], dterm=c['dterm'])
    assert abs(e - dca.eur(c['qi'], c['di'], c['b'], 600.0)) < 1e-9


def test_mh_bad_dterm():
    try:
        dca.mh_rate(1000, 0.005, 100.0, b=2.0, dterm=0.006)  # dterm >= di
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_fit_mh_recovers_parameters():
    """Noiseless MH data: fit recovers qi, di, b, dterm."""
    t = np.arange(15, 2900, 30.0)
    q = dca.mh_rate(1200, 0.008, t, b=1.6, dterm=0.0006)
    result = dca.fit_decline(t, q, method='mh')
    assert result.method == 'mh'
    assert abs(result.qi - 1200) < 1.0
    assert abs(result.di - 0.008) < 1e-4
    assert abs(result.b - 1.6) < 0.01
    assert abs(result.dterm - 0.0006) < 1e-5
    assert result.r_squared > 0.9999


def test_fit_best_excludes_mh():
    """'best' never returns the mh model (explicit opt-in only)."""
    t = np.arange(15, 2900, 30.0)
    q = dca.mh_rate(1200, 0.008, t, b=1.6, dterm=0.0006)
    result = dca.fit_decline(t, q, method='best')
    assert result.method != 'mh'


def test_forecast_mh():
    """forecast() dispatches on method='mh' with analytic cumulative."""
    c = _MH_CASE
    res = dca.DeclineResult(method='mh', qi=c['qi'], di=c['di'], b=c['b'], dterm=c['dterm'])
    fc = dca.forecast(res, t_end=3000, dt=10)
    assert abs(fc.q[-1] - dca.mh_rate(c['qi'], c['di'], fc.t[-1], b=c['b'], dterm=c['dterm'])) < 1e-9
    assert abs(fc.eur - dca.mh_cum(c['qi'], c['di'], fc.t[-1], b=c['b'], dterm=c['dterm'])) < 1e-6
    assert np.all(np.diff(fc.Qcum) > 0)


# ======================== two-segment hyperbolic (hyp2) tests ========================

# Reference case from 'Simple DCA Generator-2Stage v3.xlsm' (Time Specified
# sheet): EUR 1.5 BCF, qi 1500 mscf/d, qa 100 mscf/d, t2 24 months, b1=2,
# b2=0.5. Workbook times are months * 30.4 days/month; declines nominal /day.
_H2 = dict(qi=1500.0, di=0.005201588333851846, b1=2.0, b2=0.5, telf=24 * 30.4)


def test_hyp2_rate_matches_arps_before_switch():
    c = _H2
    q = dca.hyp2_rate(c['qi'], c['di'], 300.0, c['telf'], b1=c['b1'], b2=c['b2'])
    assert abs(q - dca.arps_rate(c['qi'], c['di'], c['b1'], 300.0)) < 1e-9


def test_hyp2_rate_continuous_at_switch():
    c = _H2
    q_lo = dca.hyp2_rate(c['qi'], c['di'], c['telf'] - 1e-6, c['telf'], b1=c['b1'], b2=c['b2'])
    q_hi = dca.hyp2_rate(c['qi'], c['di'], c['telf'] + 1e-6, c['telf'], b1=c['b1'], b2=c['b2'])
    assert abs(q_lo - q_hi) < 1e-5


def test_hyp2_slope_continuous_at_switch():
    """Nominal decline (log-slope) is identical either side of telf."""
    c = _H2
    h, eps = 1e-4, 1e-6

    def logq(t):
        return np.log(dca.hyp2_rate(c['qi'], c['di'], t, c['telf'], b1=c['b1'], b2=c['b2']))

    D_lo = -(logq(c['telf'] - eps) - logq(c['telf'] - eps - h)) / h
    D_hi = -(logq(c['telf'] + eps + h) - logq(c['telf'] + eps)) / h
    assert abs(D_lo - D_hi) / D_lo < 1e-4


def test_hyp2_second_segment_form():
    """Beyond telf the curve is hyperbolic b2 anchored at (q_sw, D_sw)."""
    c = _H2
    q_sw = dca.arps_rate(c['qi'], c['di'], c['b1'], c['telf'])
    D_sw = c['di'] / (1.0 + c['b1'] * c['di'] * c['telf'])
    t = 2500.0
    q = dca.hyp2_rate(c['qi'], c['di'], t, c['telf'], b1=c['b1'], b2=c['b2'])
    expected = dca.arps_rate(q_sw, D_sw, c['b2'], t - c['telf'])
    assert abs(q - expected) < 1e-9
    # Workbook q2i cross-check (Time Specified sheet D9)
    assert abs(q_sw - 511.78869771853385) / 511.78869771853385 < 1e-9


def test_hyp2_cum_matches_numeric_integral():
    c = _H2
    t_end = 3.0 * c['telf']
    t = np.linspace(1e-9, t_end, 200001)
    q = dca.hyp2_rate(c['qi'], c['di'], t, c['telf'], b1=c['b1'], b2=c['b2'])
    numeric = np.trapezoid(q, t)
    analytic = dca.hyp2_cum(c['qi'], c['di'], t_end, c['telf'], b1=c['b1'], b2=c['b2'])
    assert abs(numeric / analytic - 1.0) < 1e-6


def test_hyp2_eur_both_segments():
    c = _H2
    # Abandonment in first segment: identical to single-segment eur
    e1 = dca.hyp2_eur(c['qi'], c['di'], 800.0, c['telf'], b1=c['b1'], b2=c['b2'])
    assert abs(e1 - dca.eur(c['qi'], c['di'], c['b1'], 800.0)) < 1e-9
    # Abandonment in second segment: workbook target EUR (1.5 BCF = 1.5e6 mscf)
    e2 = dca.hyp2_eur(c['qi'], c['di'], 100.0, c['telf'], b1=c['b1'], b2=c['b2'])
    assert abs(e2 - 1.5e6) / 1.5e6 < 1e-5


def test_hyp2_dterm_tail():
    """Optional terminal decline applies to the second segment (mh form)."""
    c = _H2
    dterm = 1e-4
    D_sw = c['di'] / (1.0 + c['b1'] * c['di'] * c['telf'])
    q_sw = dca.arps_rate(c['qi'], c['di'], c['b1'], c['telf'])
    t = 8000.0
    q = dca.hyp2_rate(c['qi'], c['di'], t, c['telf'], b1=c['b1'], b2=c['b2'], dterm=dterm)
    expected = dca.mh_rate(q_sw, D_sw, t - c['telf'], b=c['b2'], dterm=dterm)
    assert abs(q - expected) < 1e-9


def test_hyp2_b2_exceeds_b1_raises():
    try:
        dca.hyp2_rate(1000, 0.005, 100.0, 500.0, b1=0.5, b2=2.0)
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_hyp2_from_eur_telf_mode():
    """Workbook 'Time Specified' sheet regression."""
    r = dca.hyp2_from_eur(1.5e6, 1500.0, 100.0, b1=2.0, b2=0.5, telf=24 * 30.4)
    assert r.method == 'hyp2'
    assert abs(r.di - 0.005201588333851846) / 0.005201588333851846 < 1e-4
    assert abs(dca.hyp2_eur(r.qi, r.di, 100.0, r.telf, r.b, r.b2) - 1.5e6) / 1.5e6 < 1e-9


def test_hyp2_from_eur_d_sw_mode():
    """Workbook 'Switch Decline Specified' sheet. The VBA freezes t2 from its
    initial guess (delivering 8.811 BCF vs the 8.8 target when evaluated
    self-consistently); the brentq port refreshes t2, so di differs from the
    workbook cache by ~0.24% and the EUR round-trip is exact."""
    r = dca.hyp2_from_eur(8.8e6, 10130.0, 100.0, b1=2.0, b2=0.5, d_sw=0.1 / 365.0)
    assert r.method == 'hyp2'
    assert abs(r.di - 0.025095477342955762) / 0.025095477342955762 < 5e-3
    # Switch is exactly at D = d_sw and the EUR round-trip is exact
    assert abs(r.di / (1.0 + r.b * r.di * r.telf) - 0.1 / 365.0) / (0.1 / 365.0) < 1e-9
    assert abs(dca.hyp2_eur(r.qi, r.di, 100.0, r.telf, r.b, r.b2) - 8.8e6) / 8.8e6 < 1e-9


def test_hyp2_from_eur_eur_frac_mode():
    """Workbook 'EUR Fraction Specified' sheet regression."""
    r = dca.hyp2_from_eur(8.8e6, 10130.0, 100.0, b1=2.0, b2=0.5, eur_frac=0.3)
    assert abs(r.di - 0.009009606094717059) / 0.009009606094717059 < 1e-4
    assert abs(r.telf / 30.4 - 18.63724102605943) / 18.63724102605943 < 1e-4
    assert abs(dca.hyp2_eur(r.qi, r.di, 100.0, r.telf, r.b, r.b2) - 8.8e6) / 8.8e6 < 1e-6


def test_hyp2_from_eur_rate_frac_mode():
    """Workbook 'Rate Fraction Specified' sheet regression (direct 1/di scaling)."""
    r = dca.hyp2_from_eur(8.8e6, 10130.0, 100.0, b1=2.0, b2=0.5, rate_frac=0.25)
    assert abs(r.di - 0.010832538161681804) / 0.010832538161681804 < 1e-9
    assert abs(dca.arps_rate(r.qi, r.di, r.b, r.telf) - 2532.5) / 2532.5 < 1e-9


def test_hyp2_from_eur_spec_validation():
    for kwargs in [dict(), dict(telf=100.0, d_sw=0.001), dict(rate_frac=0.001)]:
        try:
            dca.hyp2_from_eur(1e6, 1000.0, 100.0, **kwargs)
            assert False, f"Should raise ValueError for {kwargs}"
        except ValueError:
            pass


def test_fit_hyp2_recovers_parameters():
    """Noiseless hyp2 data: fit recovers qi, di, b1, b2, telf."""
    t = np.arange(15, 2900, 30.0)
    q = dca.hyp2_rate(1500.0, 0.0052, t, 730.0, b1=1.8, b2=0.5)
    result = dca.fit_decline(t, q, method='hyp2')
    assert result.method == 'hyp2'
    assert abs(result.qi - 1500.0) / 1500.0 < 0.01
    assert abs(result.di - 0.0052) / 0.0052 < 0.01
    assert abs(result.b - 1.8) / 1.8 < 0.02
    assert abs(result.b2 - 0.5) / 0.5 < 0.02
    assert abs(result.telf - 730.0) / 730.0 < 0.02
    assert result.r_squared > 0.9999


def test_fit_best_excludes_hyp2():
    t = np.arange(15, 2900, 30.0)
    q = dca.hyp2_rate(1500.0, 0.0052, t, 730.0, b1=1.8, b2=0.5)
    result = dca.fit_decline(t, q, method='best')
    assert result.method not in ('hyp2', 'mh')


def test_forecast_hyp2():
    """forecast() dispatches on method='hyp2' with analytic cumulative."""
    c = _H2
    res = dca.DeclineResult(method='hyp2', qi=c['qi'], di=c['di'], b=c['b1'],
                            b2=c['b2'], telf=c['telf'])
    fc = dca.forecast(res, t_end=4000, dt=10)
    assert abs(fc.q[-1] - dca.hyp2_rate(c['qi'], c['di'], fc.t[-1], c['telf'],
                                        b1=c['b1'], b2=c['b2'])) < 1e-9
    assert abs(fc.eur - dca.hyp2_cum(c['qi'], c['di'], fc.t[-1], c['telf'],
                                     b1=c['b1'], b2=c['b2'])) < 1e-6
    assert np.all(np.diff(fc.Qcum) > 0)


# ======================== fit_decline tests ========================

def test_fit_exponential():
    """Fit exponential decline and recover parameters."""
    t = np.arange(1, 51, dtype=float)
    q_true = 1000 * np.exp(-0.05 * t)
    result = dca.fit_decline(t, q_true, method='exponential')
    assert result.method == 'exponential'
    assert abs(result.qi - 1000) < 50
    assert abs(result.di - 0.05) < 0.01
    assert result.r_squared > 0.99


def test_fit_harmonic():
    """Fit harmonic decline."""
    t = np.arange(1, 51, dtype=float)
    q_true = 1000 / (1 + 0.05 * t)
    result = dca.fit_decline(t, q_true, method='harmonic')
    assert result.method == 'harmonic'
    assert abs(result.qi - 1000) / 1000 < 0.05
    assert result.r_squared > 0.99


def test_fit_hyperbolic():
    """Fit hyperbolic decline."""
    t = np.arange(1, 51, dtype=float)
    q_true = 1000 / (1 + 0.5 * 0.1 * t) ** (1.0 / 0.5)
    result = dca.fit_decline(t, q_true, method='hyperbolic')
    assert result.method == 'hyperbolic'
    assert result.r_squared > 0.95


def test_fit_best():
    """Best fit should pick the right model for exponential data."""
    t = np.arange(1, 51, dtype=float)
    q_true = 1000 * np.exp(-0.05 * t)
    result = dca.fit_decline(t, q_true, method='best')
    assert result.r_squared > 0.95


def test_fit_too_few_points():
    try:
        dca.fit_decline([1, 2], [100, 90])
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_fit_bad_method():
    try:
        dca.fit_decline([1, 2, 3], [100, 90, 80], method='invalid')
        assert False, "Should raise ValueError"
    except ValueError:
        pass


# ======================== forecast tests ========================

def test_forecast_basic():
    result = dca.DeclineResult(method='exponential', qi=1000, di=0.05, b=0)
    fc = dca.forecast(result, t_end=100, dt=1.0)
    assert isinstance(fc, dca.ForecastResult)
    assert len(fc.t) == 100
    assert len(fc.q) == 100
    assert len(fc.Qcum) == 100
    assert fc.eur > 0
    # Rate should be declining
    assert fc.q[-1] < fc.q[0]


def test_forecast_with_qlimit():
    result = dca.DeclineResult(method='exponential', qi=1000, di=0.1, b=0)
    fc = dca.forecast(result, t_end=200, dt=1.0, q_min=100)
    # Should truncate before t_end
    assert len(fc.t) < 200
    assert fc.q[-1] >= 100


def test_forecast_cumulative_increasing():
    result = dca.DeclineResult(method='exponential', qi=1000, di=0.05, b=0)
    fc = dca.forecast(result, t_end=50, dt=1.0)
    # Cumulative should be monotonically increasing
    assert np.all(np.diff(fc.Qcum) > 0)


# ======================== fit_decline_cum tests ========================

def test_fit_cum_exponential_roundtrip():
    """Generate exponential q(t) -> compute Np -> fit_decline_cum -> recover qi, di."""
    qi, di = 1000.0, 0.05
    t = np.arange(1, 101, dtype=float)
    q = qi * np.exp(-di * t)
    Np = np.array([float(dca.arps_cum(qi, di, 0, ti)) for ti in t])
    result = dca.fit_decline_cum(Np, q, method='exponential')
    assert result.method == 'exponential'
    assert abs(result.qi - qi) / qi < 0.02
    assert abs(result.di - di) / di < 0.02
    assert result.r_squared > 0.99


def test_fit_cum_harmonic_roundtrip():
    """Generate harmonic q(t) -> compute Np -> fit_decline_cum -> recover qi, di."""
    qi, di = 800.0, 0.08
    t = np.arange(1, 101, dtype=float)
    q = qi / (1.0 + di * t)
    Np = np.array([float(dca.arps_cum(qi, di, 1.0, ti)) for ti in t])
    result = dca.fit_decline_cum(Np, q, method='harmonic')
    assert result.method == 'harmonic'
    assert abs(result.qi - qi) / qi < 0.05
    assert abs(result.di - di) / di < 0.05
    assert result.r_squared > 0.99


def test_fit_cum_hyperbolic_roundtrip():
    """Generate hyperbolic q(t) -> compute Np -> fit_decline_cum -> recover qi, di."""
    qi, di, b = 1000.0, 0.1, 0.5
    t = np.arange(1, 101, dtype=float)
    q = np.array([float(dca.arps_rate(qi, di, b, ti)) for ti in t])
    Np = np.array([float(dca.arps_cum(qi, di, b, ti)) for ti in t])
    result = dca.fit_decline_cum(Np, q, method='hyperbolic')
    assert result.method == 'hyperbolic'
    assert result.r_squared > 0.95


def test_fit_cum_best_selection():
    """Best method picks the right model for exponential data."""
    qi, di = 1000.0, 0.05
    t = np.arange(1, 101, dtype=float)
    q = qi * np.exp(-di * t)
    Np = np.array([float(dca.arps_cum(qi, di, 0, ti)) for ti in t])
    result = dca.fit_decline_cum(Np, q, method='best')
    assert result.r_squared > 0.95


def test_fit_cum_duong_rejected():
    """Duong should raise ValueError for cumulative fitting."""
    try:
        dca.fit_decline_cum([1, 2, 3], [100, 90, 80], method='duong')
        assert False, "Should raise ValueError"
    except ValueError as e:
        assert 'Duong' in str(e)


def test_fit_cum_bad_inputs():
    """Too few points and mismatched lengths raise ValueError."""
    try:
        dca.fit_decline_cum([1, 2], [100, 90])
        assert False, "Should raise ValueError"
    except ValueError:
        pass
    try:
        dca.fit_decline_cum([1, 2, 3], [100, 90])
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_fit_cum_uptime_inference():
    """Generate data with known uptime gaps, verify uptime_mean ≈ 0.8."""
    qi, di = 1000.0, 0.05
    uptime_true = 0.8
    # Generate capacity rates at integer times
    t = np.arange(1, 51, dtype=float)
    q_capacity = qi * np.exp(-di * t)
    # Actual production in each interval: q_capacity * uptime * dt
    # Cumulative = sum of actual production
    dt = 1.0
    q_actual = q_capacity * uptime_true
    Np = np.cumsum(q_actual * dt)
    # Calendar time = t (each step is 1 unit of calendar time)
    t_calendar = t
    result = dca.fit_decline_cum(Np, q_capacity, t_calendar=t_calendar)
    assert result.uptime_mean is not None
    assert abs(result.uptime_mean - uptime_true) < 0.05
    assert result.uptime_history is not None
    assert len(result.uptime_history) == len(t) - 1


def test_fit_cum_to_forecast_roundtrip():
    """fit_decline_cum result works with forecast() for cumulative consistency."""
    qi, di = 1000.0, 0.05
    t = np.arange(1, 101, dtype=float)
    q = qi * np.exp(-di * t)
    Np = np.array([float(dca.arps_cum(qi, di, 0, ti)) for ti in t])
    result = dca.fit_decline_cum(Np, q, method='exponential')
    fc = dca.forecast(result, t_end=100, dt=1.0)
    assert fc.eur > 0
    assert len(fc.t) == 100
    # Cumulative should be monotonically increasing
    assert np.all(np.diff(fc.Qcum) > 0)


# ======================== fit_ratio tests ========================

def test_fit_ratio_linear():
    """Fit linear ratio R = a + b*x."""
    x = np.arange(1, 51, dtype=float)
    ratio = 0.5 + 0.02 * x
    result = dca.fit_ratio(x, ratio, method='linear')
    assert result.method == 'linear'
    assert abs(result.a - 0.5) < 0.01
    assert abs(result.b - 0.02) < 0.001
    assert result.r_squared > 0.99


def test_fit_ratio_exponential():
    """Fit exponential ratio R = a * exp(b*x)."""
    x = np.arange(1, 51, dtype=float)
    ratio = 0.5 * np.exp(0.02 * x)
    result = dca.fit_ratio(x, ratio, method='exponential')
    assert result.method == 'exponential'
    assert abs(result.a - 0.5) / 0.5 < 0.02
    assert abs(result.b - 0.02) / 0.02 < 0.02
    assert result.r_squared > 0.99


def test_fit_ratio_power():
    """Fit power ratio R = a * x^b."""
    x = np.arange(1, 51, dtype=float)
    ratio = 2.0 * x ** 0.5
    result = dca.fit_ratio(x, ratio, method='power')
    assert result.method == 'power'
    assert abs(result.a - 2.0) / 2.0 < 0.02
    assert abs(result.b - 0.5) / 0.5 < 0.02
    assert result.r_squared > 0.99


def test_fit_ratio_logistic():
    """Fit logistic ratio R = Rmax / (1 + c*exp(-b*x))."""
    x = np.arange(1, 101, dtype=float)
    Rmax, b_true, c_true = 10.0, 0.1, 50.0
    ratio = Rmax / (1.0 + c_true * np.exp(-b_true * x))
    result = dca.fit_ratio(x, ratio, method='logistic')
    assert result.method == 'logistic'
    assert abs(result.a - Rmax) / Rmax < 0.1
    assert result.r_squared > 0.95


def test_fit_ratio_best():
    """Best method picks the right model."""
    x = np.arange(1, 51, dtype=float)
    ratio = 0.5 + 0.02 * x
    result = dca.fit_ratio(x, ratio, method='best')
    assert result.r_squared > 0.95


def test_ratio_forecast_evaluation():
    """ratio_forecast evaluates a fitted model correctly."""
    rr = dca.RatioResult(method='linear', a=1.0, b=0.5, domain='cum')
    r = dca.ratio_forecast(rr, 10.0)
    assert abs(r - 6.0) < 0.01
    # Array input
    r_arr = dca.ratio_forecast(rr, [5.0, 10.0])
    assert abs(r_arr[0] - 3.5) < 0.01
    assert abs(r_arr[1] - 6.0) < 0.01


def test_fit_ratio_bad_inputs():
    """Too few points and mismatched lengths raise ValueError."""
    try:
        dca.fit_ratio([1, 2], [0.5, 0.6])
        assert False, "Should raise ValueError"
    except ValueError:
        pass
    try:
        dca.fit_ratio([1, 2, 3], [0.5, 0.6])
        assert False, "Should raise ValueError"
    except ValueError:
        pass


# ======================== extended forecast tests ========================

def test_forecast_uptime_reduces_eur():
    """uptime=0.8 should give ~80% of full-uptime EUR."""
    result = dca.DeclineResult(method='exponential', qi=1000, di=0.05, b=0)
    fc_full = dca.forecast(result, t_end=100, dt=1.0, uptime=1.0)
    fc_80 = dca.forecast(result, t_end=100, dt=1.0, uptime=0.8)
    assert abs(fc_80.eur / fc_full.eur - 0.8) < 0.01


def test_forecast_ratio_cum_domain():
    """Ratio with cum domain evaluated against Qcum."""
    result = dca.DeclineResult(method='exponential', qi=1000, di=0.05, b=0)
    rr = dca.RatioResult(method='linear', a=0.5, b=0.001, domain='cum')
    fc = dca.forecast(result, t_end=50, dt=1.0, ratios={'GOR': rr})
    assert fc.secondary is not None
    assert 'GOR' in fc.secondary
    assert len(fc.secondary['GOR']['ratio']) == len(fc.t)
    assert len(fc.secondary['GOR']['rate']) == len(fc.t)
    assert len(fc.secondary['GOR']['cum']) == len(fc.t)
    # Ratio should increase with cumulative
    assert fc.secondary['GOR']['ratio'][-1] > fc.secondary['GOR']['ratio'][0]


def test_forecast_ratio_time_domain():
    """Ratio with time domain evaluated against t."""
    result = dca.DeclineResult(method='exponential', qi=1000, di=0.05, b=0)
    rr = dca.RatioResult(method='linear', a=0.5, b=0.01, domain='time')
    fc = dca.forecast(result, t_end=50, dt=1.0, ratios={'WOR': rr})
    assert fc.secondary is not None
    assert 'WOR' in fc.secondary
    # WOR ratio at t=1 should be ~0.51, at t=50 should be ~1.0
    assert abs(fc.secondary['WOR']['ratio'][0] - 0.51) < 0.01
    assert abs(fc.secondary['WOR']['ratio'][-1] - 1.0) < 0.01


def test_forecast_multiple_ratios():
    """Multiple ratios (GOR + WOR) simultaneously."""
    result = dca.DeclineResult(method='exponential', qi=1000, di=0.05, b=0)
    gor = dca.RatioResult(method='linear', a=0.5, b=0.001, domain='cum')
    wor = dca.RatioResult(method='linear', a=0.1, b=0.01, domain='time')
    fc = dca.forecast(result, t_end=50, dt=1.0, ratios={'GOR': gor, 'WOR': wor})
    assert fc.secondary is not None
    assert 'GOR' in fc.secondary
    assert 'WOR' in fc.secondary


def test_forecast_backward_compat():
    """No new params produces identical behavior to original."""
    result = dca.DeclineResult(method='exponential', qi=1000, di=0.05, b=0)
    fc = dca.forecast(result, t_end=50, dt=1.0)
    assert fc.secondary is None
    assert len(fc.t) == 50
    assert fc.eur > 0
    # Rate should be declining
    assert fc.q[-1] < fc.q[0]


# ======================== fit_decline windowing tests ========================

def test_fit_decline_window_start_only():
    """Window t_start=20 on exponential data. qi ≈ rate at t=20, di ≈ 0.05."""
    t = np.arange(1, 101, dtype=float)
    q = 1000 * np.exp(-0.05 * t)
    result = dca.fit_decline(t, q, method='exponential', t_start=20)
    expected_qi = 1000 * np.exp(-0.05 * 20)  # 367.879...
    assert abs(result.qi - expected_qi) / expected_qi < 0.01
    assert abs(result.di - 0.05) / 0.05 < 0.01
    assert result.r_squared > 0.99


def test_fit_decline_window_both():
    """Window t_start=20, t_end=60. qi ≈ rate at t=20, R² > 0.99."""
    t = np.arange(1, 101, dtype=float)
    q = 1000 * np.exp(-0.05 * t)
    result = dca.fit_decline(t, q, method='exponential', t_start=20, t_end=60)
    expected_qi = 1000 * np.exp(-0.05 * 20)
    assert abs(result.qi - expected_qi) / expected_qi < 0.01
    assert abs(result.di - 0.05) / 0.05 < 0.01
    assert result.r_squared > 0.99


def test_fit_decline_window_end_only():
    """Window t_end=50 only. Shift by t[0]. R² > 0.99."""
    t = np.arange(1, 101, dtype=float)
    q = 1000 * np.exp(-0.05 * t)
    result = dca.fit_decline(t, q, method='exponential', t_end=50)
    # t[0]=1, so qi ≈ rate at t=1
    expected_qi = 1000 * np.exp(-0.05 * 1)
    assert abs(result.qi - expected_qi) / expected_qi < 0.01
    assert result.r_squared > 0.99


def test_fit_decline_window_too_few():
    """Window leaves < 3 points. ValueError."""
    t = np.arange(1, 101, dtype=float)
    q = 1000 * np.exp(-0.05 * t)
    try:
        dca.fit_decline(t, q, method='exponential', t_start=99, t_end=100)
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_fit_decline_window_empty():
    """Window outside data range. ValueError."""
    t = np.arange(1, 101, dtype=float)
    q = 1000 * np.exp(-0.05 * t)
    try:
        dca.fit_decline(t, q, method='exponential', t_start=200)
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_fit_decline_window_backward_compat():
    """No t_start/t_end. Same result as before."""
    t = np.arange(1, 51, dtype=float)
    q = 1000 * np.exp(-0.05 * t)
    result = dca.fit_decline(t, q, method='exponential')
    assert result.method == 'exponential'
    assert abs(result.qi - 1000.0000000000007) / 1000.0 < 1e-4
    assert abs(result.di - 0.05000000000000006) / 0.05 < 1e-4
    assert abs(result.r_squared - 1.0) < 1e-6


def test_fit_decline_window_best():
    """Windowed with method='best'. Valid result."""
    t = np.arange(1, 101, dtype=float)
    q = 1000 * np.exp(-0.05 * t)
    result = dca.fit_decline(t, q, method='best', t_start=20, t_end=60)
    assert result.r_squared > 0.99
    assert result.qi > 0


def test_fit_decline_cum_window_basic():
    """Cumulative window: Np_start at midpoint. R² > 0.95."""
    qi, di = 1000.0, 0.05
    t = np.arange(1, 101, dtype=float)
    q = qi * np.exp(-di * t)
    Np = np.array([float(dca.arps_cum(qi, di, 0, ti)) for ti in t])
    Np_mid = Np[49]
    result = dca.fit_decline_cum(Np, q, method='exponential', Np_start=Np_mid)
    assert abs(result.di - 0.05) / 0.05 < 0.05
    assert result.r_squared > 0.95


def test_fit_decline_cum_window_both():
    """Cumulative window: Np_start and Np_end. R² > 0.95."""
    qi, di = 1000.0, 0.05
    t = np.arange(1, 101, dtype=float)
    q = qi * np.exp(-di * t)
    Np = np.array([float(dca.arps_cum(qi, di, 0, ti)) for ti in t])
    result = dca.fit_decline_cum(Np, q, method='exponential',
                                  Np_start=Np[20], Np_end=Np[70])
    assert abs(result.di - 0.05) / 0.05 < 0.05
    assert result.r_squared > 0.95


def test_fit_decline_cum_window_with_tcalendar():
    """Windowed fit_decline_cum with t_calendar. uptime_history length correct."""
    qi, di = 1000.0, 0.05
    t = np.arange(1, 101, dtype=float)
    q = qi * np.exp(-di * t)
    Np = np.array([float(dca.arps_cum(qi, di, 0, ti)) for ti in t])
    t_calendar = t.copy()
    result = dca.fit_decline_cum(Np, q, method='exponential',
                                  Np_start=Np[20], Np_end=Np[70],
                                  t_calendar=t_calendar)
    assert result.uptime_history is not None
    # Window from index 20 to 70 = 51 points, uptime_history = n-1 = 50
    assert len(result.uptime_history) == 50


# ======================== RANSAC outlier robustness tests ========================

def test_fit_exponential_outliers():
    """Exponential fit with outliers should still recover parameters.
    RANSAC recovers correct qi/di even though R-squared (computed against
    all points including outliers) may be low."""
    t = np.arange(1, 101, dtype=float)
    q = 1000 * np.exp(-0.05 * t)
    # Add 5 large outliers
    rng = np.random.RandomState(123)
    outlier_idx = rng.choice(len(t), 5, replace=False)
    q_noisy = q.copy()
    q_noisy[outlier_idx] = q[outlier_idx] * rng.uniform(3, 10, 5)
    result = dca.fit_decline(t, q_noisy, method='exponential')
    assert abs(result.qi - 1000) / 1000 < 0.1
    assert abs(result.di - 0.05) / 0.05 < 0.1


def test_fit_harmonic_outliers():
    """Harmonic fit with outliers should still recover parameters."""
    t = np.arange(1, 101, dtype=float)
    q = 1000 / (1 + 0.05 * t)
    q_noisy = q.copy()
    q_noisy[10] = q[10] * 5
    q_noisy[50] = q[50] * 8
    result = dca.fit_decline(t, q_noisy, method='harmonic')
    assert abs(result.qi - 1000) / 1000 < 0.15


def test_fit_hyperbolic_outliers():
    """Hyperbolic fit with outliers should still recover reasonable b.
    Rate outliers in the transformed q^(-b) space are compressed, so
    parameter recovery is less precise than exponential/harmonic cases."""
    t = np.arange(1, 101, dtype=float)
    b_true = 0.5
    q = 1000 / (1 + b_true * 0.1 * t) ** (1.0 / b_true)
    q_noisy = q.copy()
    q_noisy[5] = q[5] * 6
    q_noisy[30] = q[30] * 4
    q_noisy[70] = q[70] * 5
    result = dca.fit_decline(t, q_noisy, method='hyperbolic')
    assert 0.2 < result.b < 0.9
    assert result.qi > 0
    assert result.di > 0


def test_fit_ratio_linear_outliers():
    """Linear ratio fit with outliers should recover slope and intercept."""
    x = np.arange(1, 51, dtype=float)
    ratio = 0.5 + 0.02 * x
    ratio_noisy = ratio.copy()
    ratio_noisy[5] = 20.0  # big outlier
    ratio_noisy[25] = -5.0  # big outlier
    result = dca.fit_ratio(x, ratio_noisy, method='linear')
    assert abs(result.a - 0.5) < 0.2
    assert abs(result.b - 0.02) < 0.005


def test_fit_cum_exponential_outliers():
    """Cumulative exponential fit with outliers should recover di."""
    qi, di = 1000.0, 0.05
    t = np.arange(1, 101, dtype=float)
    q = qi * np.exp(-di * t)
    Np = np.array([float(dca.arps_cum(qi, di, 0, ti)) for ti in t])
    q_noisy = q.copy()
    q_noisy[20] = q[20] * 5
    q_noisy[60] = q[60] * 4
    result = dca.fit_decline_cum(Np, q_noisy, method='exponential')
    assert abs(result.di - 0.05) / 0.05 < 0.15


def test_fit_hyperbolic_cum_outliers():
    """Cumulative hyperbolic fit with outliers should recover reasonable b."""
    qi, di, b = 1000.0, 0.1, 0.5
    t = np.arange(1, 101, dtype=float)
    q = np.array([float(dca.arps_rate(qi, di, b, ti)) for ti in t])
    Np = np.array([float(dca.arps_cum(qi, di, b, ti)) for ti in t])
    q_noisy = q.copy()
    q_noisy[15] = q[15] * 4
    q_noisy[45] = q[45] * 3
    result = dca.fit_decline_cum(Np, q_noisy, method='hyperbolic')
    assert 0.1 < result.b < 0.9


# ======================== forecast analytic cumulative regression ========================

def test_forecast_eur_matches_analytic_eur():
    """Regression: forecast now uses analytic cumulative (arps_cum/duong_cum)
    instead of right-rectangle cumsum, which biased EUR low by ~half a
    timestep of production. Truncated forecast EUR must match eur() to
    within grid resolution."""
    result = dca.DeclineResult(method='exponential', qi=1000, di=0.05, b=0)
    fc = dca.forecast(result, t_end=500, dt=1.0, q_min=50)
    analytic = dca.eur(qi=1000, di=0.05, b=0, q_min=50)
    assert abs(fc.eur - analytic) / analytic < 0.01

    # Untruncated case: Qcum at t equals arps_cum exactly
    fc2 = dca.forecast(result, t_end=50, dt=1.0)
    expected = float(dca.arps_cum(1000, 0.05, 0, 50))
    assert abs(fc2.eur - expected) / expected < 1e-12


def test_fit_best_common_subset_with_t0():
    """Regression: with a t=0 point present, method='best' compares all
    models on the common t > 0 subset (the Duong fitter drops t <= 0
    internally, so its stored r_squared covers a different subset)."""
    t = np.arange(0, 51, dtype=float)
    q = 1000 * np.exp(-0.05 * t)
    result = dca.fit_decline(t, q, method='best')
    assert result.method == 'exponential'
    assert result.r_squared > 0.999
