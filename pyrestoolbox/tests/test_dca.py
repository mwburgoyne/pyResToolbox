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
        dca.arps_rate(1000, 0.1, 1.5, 10)
        assert False, "Should raise ValueError"
    except ValueError:
        pass


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
