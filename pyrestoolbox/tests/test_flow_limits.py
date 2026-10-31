"""Tests for nodal critical gas rate, RP 14E erosional velocity and Gilbert-family chokes."""

import math
import warnings

import numpy as np
import pytest

from pyrestoolbox import nodal
from pyrestoolbox.constants import (BAR_TO_PSI, PSI_TO_BAR, MSCF_TO_SM3, FT_TO_M,
                                    STB_TO_SM3, SCF_PER_STB_TO_SM3_PER_SM3, degf_to_degc)


def _droplet_k(r, rho_l, ift):
    """Recover K in v = K [ift (rho_l - rho_g)]^0.25 / rho_g^0.5 from a result."""
    rg = r['rho_gas']
    return r['v_crit'] * math.sqrt(rg) / (ift * (rho_l - rg)) ** 0.25


@pytest.mark.parametrize("liquid,rho_l,ift,simplified", [
    # Coleman et al. (1991) simplified forms: v = 4.434 (67 - rho_g)^0.25 / rho_g^0.5
    # for water and 3.369 (45 - rho_g)^0.25 / rho_g^0.5 for condensate
    ('water', 67.0, 60.0, 4.434),
    ('condensate', 45.0, 20.0, 3.369),
])
def test_coleman_matches_published_simplified_constants(liquid, rho_l, ift, simplified):
    r = nodal.critical_gas_rate(500, 120, 0.6, 2.441, method='COLEMAN', liquid=liquid)
    k = _droplet_k(r, rho_l, ift)
    assert k * ift ** 0.25 == pytest.approx(simplified, rel=2e-4)


def test_turner_is_coleman_with_20_percent_uplift():
    t = nodal.critical_gas_rate([200, 1500], 150, 0.7, 1.995)
    c = nodal.critical_gas_rate([200, 1500], 150, 0.7, 1.995, method='COLEMAN')
    assert np.allclose(t['v_crit'] / c['v_crit'], 1.2, rtol=1e-12)
    assert _droplet_k({'v_crit': t['v_crit'][0], 'rho_gas': t['rho_gas'][0]},
                      67.0, 60.0) == pytest.approx(1.912, rel=1e-3)


def test_critical_rate_matches_hand_conversion():
    # q (MMscf/d) = 3.06 p v A / (T z), A in ft2: the field-unit form in Lea et al.,
    # Gas Well Deliquification (2008). 3.06 is rounded from 3.0552.
    p, degf, tid = 500.0, 120.0, 2.441
    r = nodal.critical_gas_rate(p, degf, 0.6, tid)
    area = math.pi * tid ** 2 / 4 / 144
    q_hand = 3.06 * p * r['v_crit'] * area / ((degf + 459.67) * r['z']) * 1000
    assert r['q_crit'] == pytest.approx(q_hand, rel=2e-3)
    assert r['q_crit'] == pytest.approx(1123.2543, rel=1e-6)


def test_critical_rate_overrides_and_metric():
    r = nodal.critical_gas_rate(800, 180, 0.65, 2.992, rho_liq=50.0, ift=30.0)
    assert _droplet_k(r, 50.0, 30.0) == pytest.approx(1.912, rel=1e-3)
    m = nodal.critical_gas_rate(800 * PSI_TO_BAR, degf_to_degc(180), 0.65, 2.992 * 25.4,
                                rho_liq=50.0 * 16.01846337, ift=30.0, metric=True)
    assert m['v_crit'] == pytest.approx(r['v_crit'] * FT_TO_M, rel=1e-9)
    assert m['q_crit'] == pytest.approx(r['q_crit'] * MSCF_TO_SM3, rel=1e-9)


def test_critical_rate_rejects_bad_inputs():
    with pytest.raises(ValueError):
        nodal.critical_gas_rate(500, 120, 0.6, 2.441, method='LEA')
    with pytest.raises(ValueError):
        nodal.critical_gas_rate(500, 120, 0.6, 0.0)
    with pytest.raises(ValueError, match="exceed gas density"):
        nodal.critical_gas_rate(500, 120, 0.6, 2.441, rho_liq=1.0)


def test_rp14e_mixture_density_limits():
    # No gas: 12409 S / 198.7 = 62.45 S lbm/ft3
    r = nodal.erosional_velocity(1000, 150, 0.7, glr=0.0, sg_liq=0.85)
    assert r['rho_mix'] == pytest.approx(12409 * 0.85 / 198.7, rel=1e-12)
    assert r['ve'] == pytest.approx(100 / math.sqrt(r['rho_mix']), rel=1e-12)
    # All gas: the RP 14E gas term 2.7 Sg P / (T Z) is the real-gas density
    # (2.7 rounds 28.97 / 10.7316 = 2.6995)
    g = nodal.erosional_velocity(2000, 120, 0.65)
    w = nodal.erosional_velocity(2000, 120, 0.65, glr=1e9)
    assert w['rho_mix'] == pytest.approx(g['rho_mix'], rel=3e-4)


def test_rp14e_max_rates():
    p, degf, glr, tid = 500.0, 120.0, 1000.0, 2.441
    r = nodal.erosional_velocity(p, degf, 0.65, glr=glr, tid=tid, c=125)
    # RP 14E Eq. 2.15: A (in2 per 1000 bbl/d) = (9.35 + Z R T / (21.25 P)) / Ve
    a_req = (9.35 + r['z'] * glr * (degf + 459.67) / (21.25 * p)) / r['ve']
    assert r['q_liq_max'] == pytest.approx(1000 * math.pi * tid ** 2 / 4 / a_req, rel=1e-12)
    g = nodal.erosional_velocity(2000, 120, 0.65, tid=tid)
    c = nodal.critical_gas_rate(2000, 120, 0.65, tid)
    # Same area and standard-rate conversion: rate scales with velocity
    assert g['q_gas_max'] / c['q_crit'] == pytest.approx(g['ve'] / c['v_crit'], rel=1e-12)


def test_rp14e_metric():
    f = nodal.erosional_velocity(1500, 140, 0.7, glr=800, tid=2.441)
    m = nodal.erosional_velocity(1500 * PSI_TO_BAR, degf_to_degc(140), 0.7,
                                 glr=800 * SCF_PER_STB_TO_SM3_PER_SM3, tid=2.441 * 25.4,
                                 metric=True)
    assert m['ve'] == pytest.approx(f['ve'] * FT_TO_M, rel=1e-9)
    assert m['q_liq_max'] == pytest.approx(f['q_liq_max'] * STB_TO_SM3, rel=1e-9)


def test_gilbert_hand_calculation():
    # p_up (psig) = 10.00 q GLR^0.546 / S^1.89
    r = nodal.choke_flow(1000, q_liq=1000, size=32)
    assert r['p_up'] == pytest.approx(10 * 1000 * 1000 ** 0.546 / 32 ** 1.89 + 14.696,
                                      rel=1e-12)
    # 10.00 is Gilbert's 435 with GLR in Mscf/STB
    assert 435 / 1000 ** 0.546 == pytest.approx(10.0, rel=5e-3)


@pytest.mark.parametrize("method", ['GILBERT', 'ROS', 'BAXENDELL', 'ACHONG'])
def test_choke_round_trip(method):
    p = nodal.choke_flow(800, q_liq=1500, size=24, method=method)['p_up']
    assert nodal.choke_flow(800, p_up=p, size=24, method=method)['q_liq'] == pytest.approx(1500, rel=1e-12)
    assert nodal.choke_flow(800, p_up=p, q_liq=1500, method=method)['size'] == pytest.approx(24, rel=1e-12)


def test_choke_critical_check_and_metric():
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        assert nodal.choke_flow(1000, p_up=1000, size=32, p_down=500)['critical'] is True
    with pytest.warns(UserWarning, match="critical-flow"):
        assert nodal.choke_flow(1000, p_up=1000, size=32, p_down=700)['critical'] is False
    f = nodal.choke_flow(1000, q_liq=1000, size=32)
    m = nodal.choke_flow(1000 * SCF_PER_STB_TO_SM3_PER_SM3, q_liq=1000 * STB_TO_SM3, size=32,
                         metric=True)
    assert m['p_up'] * BAR_TO_PSI == pytest.approx(f['p_up'], rel=1e-9)
    with pytest.raises(ValueError, match="exactly two"):
        nodal.choke_flow(1000, p_up=500)
    with pytest.raises(ValueError, match="atmospheric"):
        nodal.choke_flow(1000, p_up=10, size=32)


# --- Sachdeva et al. (1986) mechanistic choke ------------------------------------

def _fig5_fluid(x, p1, degf, sg_liq):
    from pyrestoolbox.nodal._flow_limits import _ChokeFluid
    from pyrestoolbox.constants import R, WDEN
    k, cp, cl = 1.4, 0.24, 0.8                  # Fig. 5 base case (air, kerosene-like)
    cv = cp / k
    n = 1 + x * (cp - cv) / (x * cv + (1 - x) * cl)
    return _ChokeFluid(x=x, vg1=R * (degf + 459.67) / (p1 * 28.97), vl=1 / (sg_liq * WDEN),
                       n=n, k=k, mass_per_stb=1.0)


@pytest.mark.parametrize("curve, p1, degf, sg_liq, x, y_read", [
    # Points read off Sachdeva et al. Fig. 5 (legend column headers are swapped in the
    # paper: the text makes B, C the higher P1 and D the higher T1)
    ('A', 80, 100, 0.9, 0.00029, 0.18), ('A', 80, 100, 0.9, 0.002, 0.39),
    ('A', 80, 100, 0.9, 0.019, 0.583), ('A', 80, 100, 0.9, 0.096, 0.612),
    ('A', 80, 100, 0.9, 0.58, 0.595), ('B', 400, 100, 0.9, 0.01, 0.40),
    ('B', 400, 100, 0.9, 0.041, 0.54), ('C', 800, 100, 0.9, 0.0015, 0.148),
    ('C', 800, 100, 0.9, 0.01, 0.30), ('D', 80, 200, 0.9, 0.00029, 0.21),
    ('E', 80, 100, 0.8, 0.00029, 0.17)])
def test_sachdeva_critical_ratio_matches_fig5(curve, p1, degf, sg_liq, x, y_read):
    from pyrestoolbox.nodal._flow_limits import _sachdeva_yc
    assert _sachdeva_yc(_fig5_fluid(x, p1, degf, sg_liq)) == pytest.approx(y_read, abs=0.025)


def test_sachdeva_pure_gas_limit_is_eq3():
    from pyrestoolbox.nodal._flow_limits import _sachdeva_yc
    k = 1.4
    assert _sachdeva_yc(_fig5_fluid(1 - 1e-9, 80, 100, 0.9)) == pytest.approx(
        (2 / (k + 1)) ** (k / (k - 1)), abs=1e-6)


def test_sachdeva_round_trips_regimes_and_scaling():
    kw = dict(p_down=100, degf=120, glr=800, sg=0.7, sg_liq=0.85)
    r = nodal.choke_flow_sachdeva(p_up=400, size=32, **kw)
    assert r['critical'] and r['y'] < r['y_crit']
    assert nodal.choke_flow_sachdeva(q_liq=r['q_liq'], size=32, **kw)['p_up'] == pytest.approx(400, rel=1e-9)
    assert nodal.choke_flow_sachdeva(q_liq=r['q_liq'], p_up=400, **kw)['size'] == pytest.approx(32, rel=1e-9)
    # Rate scales with throat area; in critical flow it does not see p_down
    assert nodal.choke_flow_sachdeva(p_up=400, size=16, **kw)['q_liq'] == pytest.approx(r['q_liq'] / 4)
    kw2 = dict(kw, p_down=150)
    assert nodal.choke_flow_sachdeva(p_up=400, size=32, **kw2)['q_liq'] == pytest.approx(r['q_liq'])
    # Subcritical: lower rate, continuous at the boundary
    sub = nodal.choke_flow_sachdeva(p_up=400, size=32, **dict(kw, p_down=350))
    assert not sub['critical'] and sub['q_liq'] < r['q_liq']
    edge = nodal.choke_flow_sachdeva(p_up=400, size=32, **dict(kw, p_down=400 * r['y_crit'] * (1 + 1e-9)))
    assert not edge['critical'] and edge['q_liq'] == pytest.approx(r['q_liq'], rel=1e-6)
    assert r['q_gas'] == pytest.approx(r['q_liq'] * 0.8)


def test_sachdeva_metric_matches_field():
    from pyrestoolbox.constants import PSI_TO_BAR, SCF_PER_STB_TO_SM3_PER_SM3, STB_TO_SM3
    f = nodal.choke_flow_sachdeva(p_down=100, degf=120, glr=800, p_up=400, size=32)
    m = nodal.choke_flow_sachdeva(p_down=100 * PSI_TO_BAR, degf=(120 - 32) / 1.8,
                                  glr=800 * SCF_PER_STB_TO_SM3_PER_SM3, p_up=400 * PSI_TO_BAR,
                                  size=32, cl=0.8 * 4.1868, metric=True)
    assert m['q_liq'] == pytest.approx(f['q_liq'] * STB_TO_SM3, rel=1e-9)
    assert m['y_crit'] == pytest.approx(f['y_crit'], rel=1e-9)


def test_sachdeva_errors():
    with pytest.raises(ValueError, match="exactly two"):
        nodal.choke_flow_sachdeva(p_down=100, degf=120, glr=800, p_up=400)
    with pytest.raises(ValueError, match="free gas"):
        nodal.choke_flow_sachdeva(p_down=100, degf=120, glr=0, p_up=400, size=32)
    with pytest.raises(ValueError, match="exceed p_down"):
        nodal.choke_flow_sachdeva(p_down=500, degf=120, glr=800, p_up=400, size=32)
