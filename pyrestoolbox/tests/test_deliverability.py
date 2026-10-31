"""Tests for gas.gas_deliverability (Rawlins-Schellhardt C,n and Houpeurt LIT fits)."""

import warnings

import numpy as np
import pytest

from pyrestoolbox import gas
from pyrestoolbox.constants import PSI_TO_BAR, MSCF_TO_SM3, degf_to_degc

# Ahmed, Reservoir Engineering Handbook, 4th ed. (2010), Example 8-2: three-point
# conventional (flow-after-flow) test, pr = 1952 psia
_PWF = [1700.0, 1500.0, 1300.0]
_QG = [2624.6, 4154.7, 5425.1]
_PR = 1952.0


def test_back_pressure_matches_ahmed_example_8_2():
    r = gas.gas_deliverability(_PWF, _QG, _PR, pseudo=False, method='CN')
    # Book (graphical): n = 0.87, C = 0.0169 Mscf/d/psi^2n, AOF 8980 Mscf/d at pwf = 0
    assert r.n == pytest.approx(0.87, rel=5e-3)
    assert r.c == pytest.approx(0.0169, rel=1e-2)
    assert r.aof == pytest.approx(8980, rel=1e-2)
    # Book IPR points (Mscf/d) at pwf 1800, 1600, 1000, 500 psia
    assert np.allclose(r.rate([1800, 1600, 1000, 500]), [1720, 3406, 6891, 8465], rtol=1e-2)


def test_lit_pressure_squared_matches_ahmed_example_8_2():
    r = gas.gas_deliverability(_PWF, _QG, _PR, pseudo=False, method='LIT')
    # Book reads a = 318 and b = 0.01333 off a hand-drawn line; least squares on the
    # same three points gives a = 313.7, b = 0.01441. AOF 8763 Mscf/d, IPR at 1000 psia 6862.
    assert r.a == pytest.approx(318, rel=2e-2)
    assert r.b == pytest.approx(0.01333, rel=0.1)
    assert r.aof == pytest.approx(8763, rel=1.5e-2)
    assert r.rate(1000) == pytest.approx(6862, rel=1e-2)
    assert r.test_type == 'flow_after_flow'


def test_exact_recovery_flow_after_flow_pseudopressure():
    a, b, pr, degf, sg = 2.0e4, 1.5, 3500.0, 200.0, 0.7
    pwf = np.array([3200.0, 2900.0, 2500.0, 2000.0])
    d = np.array([gas.gas_dmp(p, pr, degf, sg) for p in pwf])
    q = (-a + np.sqrt(a * a + 4 * b * d)) / (2 * b)
    r = gas.gas_deliverability(pwf, q, pr, degf=degf, sg=sg)
    assert (r.a, r.b) == (pytest.approx(a, rel=1e-9), pytest.approx(b, rel=1e-7))
    assert r.r_squared == pytest.approx(1.0, abs=1e-12)
    assert np.allclose(r.rate(pwf), q, rtol=1e-9)


def test_isochronal_slope_from_transient_level_from_stabilised():
    # Transient points share n (and b) with the stabilised curve but sit at a higher rate
    pr, n, c_t, c_s = 2500.0, 0.8, 0.05, 0.02
    pwf = np.array([2300.0, 2100.0, 1900.0])
    q = c_t * (pr ** 2 - pwf ** 2) ** n
    q_s = c_s * (pr ** 2 - 1800.0 ** 2) ** n
    r = gas.gas_deliverability(pwf, q, pr, pseudo=False, method='CN',
                               pwf_stab=1800.0, qg_stab=q_s)
    assert r.test_type == 'isochronal'
    assert (r.n, r.c) == (pytest.approx(n, rel=1e-10), pytest.approx(c_s, rel=1e-9))
    assert r.aof == pytest.approx(c_s * (pr ** 2 - 14.696 ** 2) ** n, rel=1e-9)


def test_modified_isochronal_uses_per_period_shut_in_pressure():
    a_t, a_s, b = 400.0, 900.0, 0.02
    pws = np.array([2500.0, 2480.0, 2455.0])
    pwf = np.array([2350.0, 2200.0, 2000.0])
    d = pws ** 2 - pwf ** 2
    q = (-a_t + np.sqrt(a_t ** 2 + 4 * b * d)) / (2 * b)
    pr_stab, pwf_s = 2500.0, 1900.0
    d_s = pr_stab ** 2 - pwf_s ** 2
    q_s = (-a_s + np.sqrt(a_s ** 2 + 4 * b * d_s)) / (2 * b)
    r = gas.gas_deliverability(pwf, q, pws, pseudo=False, pwf_stab=pwf_s, qg_stab=q_s,
                               pr_stab=pr_stab)
    assert r.test_type == 'modified_isochronal'
    assert (r.a, r.b) == (pytest.approx(a_s, rel=1e-9), pytest.approx(b, rel=1e-7))


def test_metric_matches_field():
    f = gas.gas_deliverability(_PWF, _QG, _PR, degf=180, sg=0.7)
    m = gas.gas_deliverability(np.array(_PWF) * PSI_TO_BAR, np.array(_QG) * MSCF_TO_SM3,
                               _PR * PSI_TO_BAR, degf=degf_to_degc(180), sg=0.7, metric=True)
    assert m.aof == pytest.approx(f.aof * MSCF_TO_SM3, rel=1e-9)
    assert m.rate(1000 * PSI_TO_BAR) == pytest.approx(f.rate(1000) * MSCF_TO_SM3, rel=1e-9)
    fc = gas.gas_deliverability(_PWF, _QG, _PR, pseudo=False, method='CN')
    mc = gas.gas_deliverability(np.array(_PWF) * PSI_TO_BAR, np.array(_QG) * MSCF_TO_SM3,
                                _PR * PSI_TO_BAR, pseudo=False, method='CN', metric=True)
    assert mc.n == pytest.approx(fc.n, rel=1e-12)
    assert mc.rate(500 * PSI_TO_BAR) == pytest.approx(fc.rate(500) * MSCF_TO_SM3, rel=1e-9)


def test_deliverability_rejects_and_warns():
    with pytest.raises(ValueError, match="below its reservoir"):
        gas.gas_deliverability([2000, 1500], [1000, 2000], 1900, pseudo=False)
    with pytest.raises(ValueError, match="degf is required"):
        gas.gas_deliverability(_PWF, _QG, _PR)
    with pytest.raises(ValueError, match="pr_stab"):
        gas.gas_deliverability(_PWF, _QG, [1952, 1950, 1948], pseudo=False)
    with pytest.raises(ValueError, match="both"):
        gas.gas_deliverability(_PWF, _QG, _PR, pseudo=False, pwf_stab=1000)
    with pytest.warns(UserWarning, match="outside"):
        # Rate rising faster than drawdown: n > 1
        gas.gas_deliverability([1900, 1800], [100, 1000], 1952, pseudo=False, method='CN')
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        gas.gas_deliverability(_PWF, _QG, _PR, pseudo=False, method='CN')
