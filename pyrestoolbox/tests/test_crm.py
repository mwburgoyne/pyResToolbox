"""Capacitance-resistance models (Sayarpour et al. 2009, JPSE 69, 227). The paper
does not tabulate its rate histories, so the fits are pinned by exact parameter
recovery from synthetic data and by the paper's closed-form solutions."""

import numpy as np
import pytest

from pyrestoolbox import crm
from pyrestoolbox.crm.crm import _simulate

T = np.arange(60, dtype=float) * 30.0


def _field(seed=0, ni=3, n_p=4, bhp=False):
    rng = np.random.default_rng(seed)
    inj = rng.uniform(500, 1500, (len(T), ni))
    f = np.array([[0.4, 0.3, 0.2, 0.1], [0.1, 0.5, 0.2, 0.2], [0.25, 0.25, 0.25, 0.2]])[:ni, :n_p]
    tau = np.array([40.0, 90.0, 200.0, 15.0])[:n_p]
    q0 = np.array([800.0, 600.0, 400.0, 900.0])[:n_p]
    J = np.array([2.0, 5.0, 1.0, 3.0])[:n_p] if bhp else None
    p = 1000 + np.cumsum(rng.normal(0, 20, (len(T), n_p)), axis=0) if bhp else None
    prod = _simulate(np.diff(T), inj, None if p is None else np.diff(p, axis=0), f, tau, q0, J)
    return inj, prod, p, f, tau, q0, J


def test_recursion_equals_paper_eq5_closed_form():
    # Eq. (5): q(tn) = q0 e^(-(tn-t0)/tau)
    #   + sum_k e^(-(tn-tk)/tau) (1 - e^(-dt_k/tau)) [sum_i f_i I_i(k) - J tau dp(k)/dt_k]
    inj, prod, p, f, tau, q0, J = _field(bhp=True)
    j, n = 1, len(T) - 1
    dt = np.diff(T)
    terms = [np.exp(-(T[n] - T[k]) / tau[j]) * (1 - np.exp(-dt[k - 1] / tau[j]))
             * (inj[k] @ f[:, j] - J[j] * tau[j] * (p[k, j] - p[k - 1, j]) / dt[k - 1])
             for k in range(1, n + 1)]
    closed = q0[j] * np.exp(-(T[n] - T[0]) / tau[j]) + sum(terms)
    assert prod[n, j] == pytest.approx(closed, rel=1e-12)


def test_crmt_matches_eq8():
    # Eq. (8): q(tk) = q(tk-1) e^(-dt/tau) + I(k) (1 - e^(-dt/tau)), field-lumped
    inj, _, _, _, _, _, _ = _field()
    q = [700.0]
    for k in range(1, len(T)):
        a = np.exp(-30.0 / 120.0)
        q.append(q[-1] * a + inj[k].sum() * (1 - a))
    r = crm.crm_fit(T, inj, np.array(q), model='CRMT')
    assert r.tau[0] == pytest.approx(120.0, rel=1e-6)
    assert r.q0[0] == pytest.approx(700.0, rel=1e-6)
    assert r.f.tolist() == [[1.0]] and r.J is None


def test_crmp_exact_recovery_with_bhp():
    inj, prod, p, f, tau, q0, J = _field(bhp=True)
    r = crm.crm_fit(T, inj, prod, bhp=p)
    np.testing.assert_allclose(r.f, f, atol=1e-6)
    np.testing.assert_allclose(r.tau, tau, rtol=1e-6)
    np.testing.assert_allclose(r.J, J, rtol=1e-6)
    np.testing.assert_allclose(r.q0, q0, rtol=1e-5)
    assert r.r_squared.min() > 1 - 1e-10
    np.testing.assert_allclose(r.pore_volume(ct=1e-5), tau * J / 1e-5, rtol=1e-6)


def test_connectivity_constraint_binds():
    # Injector 0 sends all its water to the producers; with noise the unconstrained
    # fit would exceed one, Eq. (20) caps it
    rng = np.random.default_rng(5)
    inj, _, _, f, tau, q0, _ = _field()
    f = f.copy()
    f[0] = f[0] / f[0].sum()
    prod = _simulate(np.diff(T), inj, None, f, tau, q0, None)
    prod = prod * (1 + 0.03 * rng.standard_normal(prod.shape))
    r = crm.crm_fit(T, inj, prod)
    assert r.f.sum(axis=1).max() <= 1 + 1e-9
    assert r.f.min() >= 0
    np.testing.assert_allclose(r.f, f, atol=0.08)


def test_predict_extends_history():
    inj, prod, _, f, tau, q0, _ = _field()
    r = crm.crm_fit(T[:40], inj[:40], prod[:40])
    np.testing.assert_allclose(r.predict(T, inj), prod, rtol=1e-6)
    with pytest.raises(ValueError, match="bhp"):
        r.predict(T, inj, bhp=np.ones_like(prod))
    with pytest.raises(ValueError, match="pore volume|Pore volume"):
        r.pore_volume(1e-5)


def test_oil_cut_fit_recovers_power_law():
    cwi = np.geomspace(1e5, 1e7, 20)
    fo = 1 / (1 + 3e-4 * cwi ** 0.7)
    r = crm.crm_oil_cut_fit(cwi, fo)
    assert r['alpha'] == pytest.approx(3e-4, rel=1e-9)
    assert r['beta'] == pytest.approx(0.7, rel=1e-9)
    assert r['n_points'] == 20


def test_input_errors():
    inj, prod, _, _, _, _, _ = _field()
    with pytest.raises(ValueError, match="rows"):
        crm.crm_fit(T, inj[:-1], prod)
    with pytest.raises(ValueError, match="increasing"):
        crm.crm_fit(T[::-1], inj, prod)
    with pytest.raises(ValueError, match="model"):
        crm.crm_fit(T, inj, prod, model='CRMIP')
