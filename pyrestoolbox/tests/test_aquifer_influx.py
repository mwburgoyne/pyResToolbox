"""Tests for matbal.aquifer_influx (Carter-Tracy), aquifer_constants and the gas_matbal p/z curvature diagnostic."""

import numpy as np
import pytest

from pyrestoolbox import matbal
from pyrestoolbox.simtools import influence_tables

# Ahmed, Reservoir Engineering Handbook, 4th ed. (2010), Example 10-9 monthly table:
# infinite-acting aquifer, B = 20.4 bbl/psi, tD and total pressure drop per month.
_TD = [0, 30.0892, 60.1784, 90.2676, 120.357, 150.446, 180.535, 210.624, 240.713, 270.802,
       300.891, 330.981, 361.070, 391.159, 421.248, 451.337, 481.426, 511.516, 541.071,
       571.130, 601.190, 631.249, 661.309, 691.369, 721.428]
_DP = [0, 1.06, 2.31, 3.81, 5.56, 7.55, 9.79, 12.27, 15.00, 17.98, 21.20, 24.67, 28.38,
       32.34, 36.55, 41.00, 45.70, 50.64, 55.74, 61.16, 66.84, 72.75, 78.92, 85.32, 91.98]
_WE_BOOK = {6: 7057.3, 12: 31591.5, 18: 81578.8, 24: 165119.0}


def test_carter_tracy_matches_ahmed_example_10_9():
    r = matbal.aquifer_influx(_TD, 2500 - np.array(_DP), B=20.4, td_per_t=1.0)
    # The book rounds pD to two decimals (2.11 at tD 30 against 2.124); the gap falls
    # from 1.9% at month 1 to 0.2% at month 24
    for month, we in _WE_BOOK.items():
        assert r['We'][month] == pytest.approx(we, rel=8e-3)
    assert r['We'][0] == 0.0


def test_pd_and_derivative_match_influence_tables():
    r = matbal.aquifer_influx([0, 1, 10, 100], [3000, 2990, 2950, 2900], B=10, td_per_t=1.0, red=20)
    td, pds = influence_tables([20], min_td=1, max_td=100, n_incr=2, M=7)
    assert np.allclose(r['pD'][1:], pds[0], rtol=1e-5)
    # Chatas (1953) infinite-acting pD at tD = 100 is 2.723 (reD = 20 has felt the
    # boundary by then: 2.742)
    inf = matbal.aquifer_influx([0, 100], [3000, 2990], B=10, td_per_t=1.0)
    assert inf['pD'][1] == pytest.approx(2.723, rel=2e-4)
    # Derivative against a central difference of the exact table
    _, p_lo_hi = influence_tables([20], min_td=9.9, max_td=10.1, n_incr=2, M=7)
    fd = (p_lo_hi[0][2] - p_lo_hi[0][0]) / 0.2
    assert r['dpD'][2] == pytest.approx(fd, rel=1e-4)


def test_influx_is_linear_in_B_and_unit_homogeneous():
    t = np.arange(0, 721, 30.0)
    p = 2500 - 0.12 * t
    a = matbal.aquifer_influx(t, p, B=10.0, td_per_t=0.5)
    b = matbal.aquifer_influx(t, p, B=20.0, td_per_t=0.5)
    assert np.allclose(b['We'], 2 * a['We'], rtol=1e-12)
    # Same history in bar with B per bar gives the same influx volume
    from pyrestoolbox.constants import PSI_TO_BAR
    m = matbal.aquifer_influx(t, p * PSI_TO_BAR, B=10.0 / PSI_TO_BAR, td_per_t=0.5)
    assert np.allclose(m['We'], a['We'], rtol=1e-12)


def test_finite_aquifer_gives_less_influx_than_infinite():
    t = np.arange(0, 3001, 30.0)
    p = 3000 - 0.1 * t
    inf = matbal.aquifer_influx(t, p, B=15, td_per_t=1.0)
    with pytest.warns(UserWarning, match="pseudo-steady"):
        fin = matbal.aquifer_influx(t, p, B=15, td_per_t=1.0, red=10)
    assert fin['We'][-1] < 0.5 * inf['We'][-1]


def test_aquifer_influx_rejects_bad_inputs():
    with pytest.raises(ValueError, match="increasing"):
        matbal.aquifer_influx([0, 2, 1], [3000, 2990, 2980], 10, 1)
    with pytest.raises(ValueError, match="greater than 1"):
        matbal.aquifer_influx([0, 1, 2], [3000, 2990, 2980], 10, 1, red=1.0)
    with pytest.raises(ValueError, match="positive"):
        matbal.aquifer_influx([0, 1, 2], [3000, 2990, 2980], 0, 1)


def test_pz_curvature_sign_and_value():
    gp = np.linspace(0, 40, 8)
    table = {'p': [100, 5000], 'Z': [1.0, 1.0]}  # Z = 1, so p/z = p
    straight = matbal.gas_matbal(4000 - 50 * gp, gp, 200, pvt_table=table)
    assert abs(straight.pz_curvature) < 1e-10
    c = 0.5
    up = matbal.gas_matbal(4000 - 50 * gp + c * gp ** 2, gp, 200, pvt_table=table)
    fall = abs(-50 * 40 + c * 40 ** 2)
    assert up.pz_curvature == pytest.approx(c * 40 ** 2 / 4 / fall, rel=1e-9)
    down = matbal.gas_matbal(4000 - 50 * gp - c * gp ** 2, gp, 200, pvt_table=table)
    assert down.pz_curvature < 0
    assert np.isnan(matbal.gas_matbal([4000, 3500, 3000], [0, 10, 20], 200,
                                      pvt_table=table).pz_curvature)


def test_aquifer_constants_textbook_field_factors():
    # Field forms (Ahmed, Dake): B = 1.119 phi ct h ro^2 (theta/360) rb/psi,
    # tD per day = 0.0063283 k / (phi mu ct ro^2); 1.119 and 0.0063283 are rounded.
    c = matbal.aquifer_constants(phi=0.25, ct=7e-6, h=100, ro=3000, k=50, mu=0.5)
    assert c['B'] == pytest.approx(1.119 * 0.25 * 7e-6 * 100 * 3000 ** 2, rel=2e-4)
    assert c['td_per_t'] == pytest.approx(0.0063283 * 50 / (0.25 * 0.5 * 7e-6 * 3000 ** 2), rel=1e-5)
    half = matbal.aquifer_constants(phi=0.25, ct=7e-6, h=100, ro=3000, k=50, mu=0.5, theta=180)
    assert half['B'] == pytest.approx(c['B'] / 2, rel=1e-12)
    assert half['td_per_t'] == c['td_per_t']


def test_aquifer_constants_metric_matches_field():
    from pyrestoolbox.constants import FT_TO_M, PSI_TO_BAR, BBL_TO_M3
    f = matbal.aquifer_constants(phi=0.25, ct=7e-6, h=100, ro=3000, k=50, mu=0.5)
    m = matbal.aquifer_constants(phi=0.25, ct=7e-6 / PSI_TO_BAR, h=100 * FT_TO_M,
                                 ro=3000 * FT_TO_M, k=50, mu=0.5, metric=True)
    assert m['td_per_t'] == pytest.approx(f['td_per_t'], rel=1e-9)
    # metric (bar, m) Darcy constant is the ECLIPSE METRIC 0.0085270 per day
    assert m['td_per_t'] == pytest.approx(0.0085270 * 50 / (0.25 * 0.5 * (7e-6 / PSI_TO_BAR)
                                                           * (3000 * FT_TO_M) ** 2), rel=1e-5)
    t = [0, 182.5, 365, 547.5, 730]
    p = np.array([2500, 2490, 2472, 2444, 2408])
    we_f = matbal.aquifer_influx(t, p, **f)['We']
    we_m = matbal.aquifer_influx(t, p * PSI_TO_BAR, **m)['We']
    # rm3 vs rb differ only by the library's 5.61458 ft3/bbl rounding (~6e-7)
    np.testing.assert_allclose(we_m / BBL_TO_M3, we_f, rtol=2e-6)


def test_aquifer_constants_rejects_bad_inputs():
    good = dict(phi=0.25, ct=7e-6, h=100, ro=3000, k=50, mu=0.5)
    for key, bad in (('phi', 0), ('phi', 1.5), ('ct', -1e-6), ('h', 0), ('ro', float('nan')),
                     ('k', 0), ('mu', -0.5)):
        with pytest.raises(ValueError):
            matbal.aquifer_constants(**{**good, key: bad})
    for theta in (0, 400):
        with pytest.raises(ValueError):
            matbal.aquifer_constants(**good, theta=theta)
