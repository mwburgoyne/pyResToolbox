"""Gas pipeline equations, pinned to Menon (2005), Gas Pipeline Hydraulics, Chapter 2
worked examples (internal functions with Menon's fixed Z and 14.7 psia / 520 R base)
and to round-trip and unit consistency through the public API."""

import math

import pytest

from pyrestoolbox import pipeline as pl
from pyrestoolbox.constants import IN_TO_MM, FT_TO_M, PSI_TO_BAR, BAR_TO_PSI, MSCF_TO_SM3
from pyrestoolbox.pipeline.pipeline import (_Pipe, _rate, _solve_p1, _solve_p2, _reynolds,
                                            _colebrook_f, _von_karman_ft, _aga_transmission)


def _const(z, mu=8.0e-6):
    return lambda pavg: (z, mu)


def test_menon_example_7_general_flow():
    # Book prints 380.99 MMscf/d; its arithmetic used Tb = 540 R (the "60+460 = 540"
    # slip in its variable list). With Tb = 520 R the same equation gives 366.88.
    pipe = _Pipe(d=19.25, l_mi=10, tf=530, sg=0.6, method='COLEBROOK', pb=14.7, tb=520,
                 f_fixed=0.02)
    assert _rate(1014.7, 814.7, pipe, _const(0.85))[0] / 1e6 == pytest.approx(366.88, abs=0.01)
    pipe.tb = 540
    assert _rate(1014.7, 814.7, pipe, _const(0.85))[0] / 1e6 == pytest.approx(380.99, abs=0.01)


def test_menon_example_8_reynolds_and_colebrook():
    pipe = _Pipe(d=23.0, l_mi=1, tf=530, sg=0.65, method='COLEBROOK', pb=14.73, tb=520)
    re = _reynolds(250e6, pipe, 8.0e-6)
    assert re == pytest.approx(11953115, rel=1e-6)
    f = _colebrook_f(re, 700e-6 / 23.0)
    assert round(f, 4) == 0.0101
    assert round(2 / math.sqrt(f), 2) == 19.95


def test_menon_example_9_aga():
    pipe = _Pipe(d=23.0, l_mi=1, tf=530, sg=0.59, method='AGA', pb=14.7, tb=520,
                 rough=0.00075, df=0.96)
    re = _reynolds(250e6, pipe, 8.0e-6)
    assert re == pytest.approx(10827653, rel=1e-6)
    ft = _von_karman_ft(re)
    assert round(ft, 2) == 22.16
    assert round(4 * 0.96 * math.log10(re / (1.4125 * ft)), 2) == 21.27   # partially turbulent
    assert round(_aga_transmission(re, pipe), 2) == 20.22                  # rough pipe governs


@pytest.mark.parametrize("method, p2_book", [('PANHANDLE_A', 970.81), ('PANHANDLE_B', 971.81),
                                             ('WEYMOUTH', 946.24)])
def test_menon_example_11_empirical(method, p2_book):
    pipe = _Pipe(d=17.5, l_mi=20, tf=530, sg=0.6, method=method, e=0.95, pb=14.7, tb=520)
    assert _solve_p2(1014.7, 150e6, pipe, _const(0.85)) == pytest.approx(p2_book, abs=0.01)


def test_menon_example_12_upstream_pressure():
    pipe = _Pipe(d=15.5, l_mi=100, tf=530, sg=0.6, method='COLEBROOK', pb=14.7, tb=520,
                 f_fixed=4.0 / 20.0 ** 2)
    assert _solve_p1(1014.7, 100e6, pipe, _const(0.844)) == pytest.approx(1195.14, abs=0.01)


@pytest.mark.parametrize("method", ['COLEBROOK', 'AGA', 'WEYMOUTH', 'PANHANDLE_A',
                                    'PANHANDLE_B', 'IGT'])
def test_round_trip_all_methods(method):
    kw = dict(d=17.5, length=105600, degf=70, sg=0.6, method=method, dz=300)
    r = pl.gas_pipeline(p1=1014.7, qg=150000, **kw)
    assert pl.gas_pipeline(p1=1014.7, p2=r['p2'], **kw)['qg'] == pytest.approx(150000, rel=1e-8)
    assert pl.gas_pipeline(p2=r['p2'], qg=150000, **kw)['p1'] == pytest.approx(1014.7, rel=1e-8)
    assert r['friction_factor'] == pytest.approx(4 / r['transmission_factor'] ** 2)
    assert r['reynolds'] > 1e6


def test_elevation_and_flat_limit():
    kw = dict(d=17.5, length=105600, degf=70, sg=0.6, p1=1014.7, qg=150000)
    flat = pl.gas_pipeline(**kw)['p2']
    assert pl.gas_pipeline(dz=1e-9, **kw)['p2'] == pytest.approx(flat, rel=1e-9)
    up, down = pl.gas_pipeline(dz=1000, **kw), pl.gas_pipeline(dz=-1000, **kw)
    assert up['p2'] < flat < down['p2']
    assert up['s'] > 0 > down['s']
    # Downhill line: p1 solve must find the root even though p1 = p2 already flows
    assert pl.gas_pipeline(p2=down['p2'], qg=150000, d=17.5, length=105600, degf=70, sg=0.6,
                           dz=-1000)['p1'] == pytest.approx(1014.7, rel=1e-8)


def test_metric_matches_field():
    field = pl.gas_pipeline(d=17.5, length=105600, degf=70, sg=0.6, p1=1014.7, qg=150000,
                            dz=200)
    m = pl.gas_pipeline(d=17.5 * IN_TO_MM, length=105600 * FT_TO_M, degf=(70 - 32) / 1.8,
                        sg=0.6, p1=1014.7 * PSI_TO_BAR, qg=150000 * MSCF_TO_SM3,
                        dz=200 * FT_TO_M, metric=True)
    assert m['p2'] * BAR_TO_PSI == pytest.approx(field['p2'], rel=1e-10)
    assert m['reynolds'] == pytest.approx(field['reynolds'], rel=1e-10)


def test_laminar_branch():
    r = pl.gas_pipeline(d=2, length=1000, degf=70, sg=0.6, p1=100, qg=0.5)
    assert r['reynolds'] < 2000
    assert r['friction_factor'] == pytest.approx(64 / r['reynolds'], rel=1e-8)


def test_errors():
    kw = dict(d=2, length=1e5, degf=70, sg=0.6)
    with pytest.raises(ValueError, match="exactly two"):
        pl.gas_pipeline(p1=100, **kw)
    with pytest.raises(ValueError, match="capacity"):
        pl.gas_pipeline(p1=100, qg=1e5, **kw)
    with pytest.raises(ValueError, match="No forward flow"):
        pl.gas_pipeline(p1=100, p2=100, dz=10, **kw)
    with pytest.raises(ValueError, match="method"):
        pl.gas_pipeline(p1=100, p2=90, method='FRITZSCHE', **kw)
