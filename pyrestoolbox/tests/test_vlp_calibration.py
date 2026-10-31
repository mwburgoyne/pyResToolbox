"""Tests for nodal.vlp_calibrate: multiplier recovery from synthetic well tests,
identifiability warnings, leave-one-out, PI scale and metric units."""

import warnings

import numpy as np
import pytest

from pyrestoolbox import nodal

_C = nodal.Completion(tid=2.441, length=2000, tht=100, bht=130)
_COMMON = dict(gor=600, pb=2800, rsb=600, api=33, wc=0.3)


def _tests(rates, gm=1.0, fm=1.0, thp=200, method='BB'):
    return [dict(thp=thp, qt_stbpd=q,
                 bhp=nodal.fbhp(thp=thp, completion=_C, vlpmethod=method, well_type='oil',
                                qt_stbpd=q, gravity_mult=gm, friction_mult=fm, **_COMMON))
            for q in rates]


def _calibrate(tests, **kw):
    kw = dict(dict(vlpmethod='BB', well_type='oil', common=_COMMON), **kw)
    return nodal.vlp_calibrate(tests, _C, **kw)


def test_recovers_friction_multiplier_with_loo():
    r = _calibrate(_tests([300, 900, 1800, 3000], fm=1.25), fit='friction')
    assert abs(r['friction_mult'] - 1.25) < 1e-3 and r['gravity_mult'] == 1.0
    assert r['rms'] < 0.05 and r['loo_rms'] < 0.05
    assert r['identifiable'] and r['at_bound'] == [] and len(r['std_err']) == 1


def test_recovers_gravity_multiplier():
    r = _calibrate(_tests([300, 1500, 3000], gm=0.92), fit='gravity', loo=False)
    assert abs(r['gravity_mult'] - 0.92) < 1e-3 and r['friction_mult'] == 1.0


def test_recovers_both_multipliers():
    r = _calibrate(_tests([300, 1200, 3000], gm=0.95, fm=1.2), fit='both', loo=False)
    assert abs(r['gravity_mult'] - 0.95) < 1e-3 and abs(r['friction_mult'] - 1.2) < 1e-3
    assert r['identifiable'] and abs(r['correlation']) < 0.95


def test_equal_rates_cannot_separate_gravity_and_friction():
    """Tests at one rate fix one combination of the multipliers, not both."""
    with pytest.warns(UserWarning, match="correlated"):
        r = _calibrate(_tests([1000, 1000, 1000], gm=0.95, fm=1.2), fit='both', loo=False)
    assert not r['identifiable']


def test_truth_outside_bounds_is_flagged():
    with pytest.warns(UserWarning, match="ended on its bound"):
        r = _calibrate(_tests([300, 1500, 3000], fm=2.5), fit='friction', loo=False)
    assert r['at_bound'] == ['friction_mult'] and r['friction_mult'] == pytest.approx(2.0)


def test_pi_scale_recovered():
    res = nodal.Reservoir(pr=3200, degf=210, k=20, h=40, re=1000, rw=0.35)
    ipr = dict(wc=0.3, bo=1.3, uo=0.8)
    tests = []
    for t in _tests([500, 1000, 1500], fm=1.1):
        q = nodal.ipr_curve(res, well_type='oil', n_points=2, min_pwf=t['bhp'], **ipr)['rate'][0]
        tests.append(dict(t, qt_stbpd=0.7 * q))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r = _calibrate(tests, reservoir=res, ipr_kwargs=ipr, loo=False)
    assert r['pi_scale'] == pytest.approx(0.7, rel=1e-9)
    assert np.allclose(r['ipr_rate_residual'], 0, atol=1e-6)


def test_metric_gas():
    c = nodal.Completion(tid=62.0, length=1000, tht=40, bht=60, metric=True)
    tests = [dict(thp=30, qg_mscfd=q, qw_bwpd=3,
                  bhp=nodal.fbhp(thp=30, completion=c, vlpmethod='HB', well_type='gas',
                                 qg_mscfd=q, qw_bwpd=3, friction_mult=1.3, metric=True))
             for q in (40000, 90000, 160000)]
    r = nodal.vlp_calibrate(tests, c, vlpmethod='HB', well_type='gas', metric=True)
    assert abs(r['friction_mult'] - 1.3) < 1e-3 and r['rms'] < 1e-3


def test_input_validation():
    with pytest.raises(ValueError, match="bounds"):
        _calibrate(_tests([500]), bounds=(2.0, 0.5))
    with pytest.raises(ValueError, match="'bhp'"):
        _calibrate([dict(thp=200, qt_stbpd=500)])
    with pytest.raises(ValueError):
        _calibrate(_tests([500]), fit='pi')
