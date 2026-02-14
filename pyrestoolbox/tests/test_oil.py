#!/usr/bin/env python3
"""
Validation tests for oil module.
Run with: PYTHONPATH=/home/mark/projects python3 -m pytest tests/ -v
Or standalone: PYTHONPATH=/home/mark/projects python3 tests/test_oil.py
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import pyrestoolbox.oil as oil

RTOL = 0.02

# =============================================================================
# Physical reasonableness tests
# =============================================================================

def test_oil_sg_api_roundtrip():
    """API ↔ SG conversion roundtrip"""
    for api in [20, 30, 40, 50]:
        sg = oil.oil_sg(api)
        api_back = oil.oil_api(sg)
        assert abs(api_back - api) < 1e-10, f"Roundtrip failed: API {api} → SG {sg} → API {api_back}"

def test_oil_sg_ranges():
    """SG should be in physical range for typical oils"""
    assert 0.8 < oil.oil_sg(30) < 0.9
    assert 0.7 < oil.oil_sg(50) < 0.8

def test_oil_pbub_standing():
    """Standing bubble point pressure"""
    pb = oil.oil_pbub(api=35, degf=200, rsb=500, sg_g=0.75, pbmethod='STAN')
    assert isinstance(pb, float)
    assert 1000 < pb < 6000, f"Pb={pb} psia outside expected range for these inputs"

def test_oil_pbub_valko_mccain():
    """Valko-McCain bubble point pressure"""
    pb = oil.oil_pbub(api=35, degf=200, rsb=500, sg_sp=0.75, pbmethod='VALMC')
    assert isinstance(pb, float)
    assert 1000 < pb < 6000, f"Pb={pb} psia outside expected range"

def test_oil_pbub_velarde():
    """Velarde bubble point pressure"""
    pb = oil.oil_pbub(api=35, degf=200, rsb=500, sg_sp=0.75, pbmethod='VELAR')
    assert isinstance(pb, float)
    assert 1000 < pb < 6000, f"Pb={pb} psia outside expected range"

def test_oil_pbub_methods_agree():
    """Different Pb methods should be within ~20% for typical inputs"""
    kwargs = dict(api=35, degf=200, rsb=500, sg_sp=0.75, sg_g=0.75)
    pb_stan = oil.oil_pbub(**kwargs, pbmethod='STAN')
    pb_valmc = oil.oil_pbub(**kwargs, pbmethod='VALMC')
    pb_velar = oil.oil_pbub(**kwargs, pbmethod='VELAR')
    pbs = [pb_stan, pb_valmc, pb_velar]
    mean_pb = np.mean(pbs)
    for pb in pbs:
        assert abs(pb - mean_pb) / mean_pb < 0.25, f"Pb methods diverge too much: {pbs}"

def test_oil_pbub_increases_with_rsb():
    """Pb should increase with increasing Rs"""
    pb1 = oil.oil_pbub(api=35, degf=200, rsb=300, sg_sp=0.75, pbmethod='VALMC')
    pb2 = oil.oil_pbub(api=35, degf=200, rsb=600, sg_sp=0.75, pbmethod='VALMC')
    assert pb2 > pb1, f"Pb should increase with Rs: Pb(300)={pb1}, Pb(600)={pb2}"

def test_oil_rs_bub_velarde():
    """Rsb from Velarde at bubble point"""
    rsb = oil.oil_rs_bub(api=35, degf=200, pb=3000, sg_sp=0.75, rsmethod='VELAR')
    assert isinstance(rsb, float)
    assert 100 < rsb < 2000, f"Rsb={rsb} outside expected range"

def test_oil_rs_bub_roundtrip():
    """Pb → Rsb → Pb roundtrip using Velarde"""
    pb_input = 3000
    rsb = oil.oil_rs_bub(api=35, degf=200, pb=pb_input, sg_sp=0.75, rsmethod='VELAR')
    pb_back = oil.oil_pbub(api=35, degf=200, rsb=rsb, sg_sp=0.75, pbmethod='VELAR')
    assert abs(pb_back - pb_input) / pb_input < 0.05, \
        f"Roundtrip: Pb={pb_input} → Rs={rsb} → Pb={pb_back}"

def test_oil_rs_below_pb():
    """Rs should decrease below Pb"""
    rsb = 500
    pb = 3000
    rs_at_2000 = oil.oil_rs(api=35, degf=200, sg_sp=0.75, p=2000, pb=pb, rsb=rsb)
    rs_at_1000 = oil.oil_rs(api=35, degf=200, sg_sp=0.75, p=1000, pb=pb, rsb=rsb)
    assert rs_at_2000 < rsb, f"Rs at 2000 should be less than Rsb"
    assert rs_at_1000 < rs_at_2000, f"Rs should decrease with pressure below Pb"

def test_oil_rs_above_pb():
    """Rs should equal Rsb above Pb"""
    rsb = 500
    pb = 3000
    rs_above = oil.oil_rs(api=35, degf=200, sg_sp=0.75, p=4000, pb=pb, rsb=rsb)
    assert abs(rs_above - rsb) < 1e-6, f"Rs above Pb should equal Rsb: {rs_above} vs {rsb}"

def test_oil_bo_at_pb():
    """Bo at bubble point should be > 1.0"""
    sg_o = oil.oil_sg(35)
    bo = oil.oil_bo(p=3000, pb=3000, degf=200, rs=500, rsb=500, sg_o=sg_o, sg_g=0.75)
    assert bo > 1.0, f"Bo at Pb should be > 1.0, got {bo}"
    assert bo < 3.0, f"Bo at Pb = {bo} seems too high"

def test_oil_bo_above_pb_decreases():
    """Bo above Pb should decrease with increasing pressure"""
    sg_o = oil.oil_sg(35)
    bo_at_pb = oil.oil_bo(p=3000, pb=3000, degf=200, rs=500, rsb=500, sg_o=sg_o, sg_g=0.75)
    bo_above = oil.oil_bo(p=5000, pb=3000, degf=200, rs=500, rsb=500, sg_o=sg_o, sg_g=0.75)
    assert bo_above < bo_at_pb, f"Bo should decrease above Pb: Bo(Pb)={bo_at_pb}, Bo(5000)={bo_above}"

def test_oil_deno():
    """Oil density at reservoir conditions"""
    rho = oil.oil_deno(p=3000, degf=200, rs=500, rsb=500, sg_sp=0.75, api=35)
    assert isinstance(rho, float)
    assert 30 < rho < 60, f"Oil density = {rho} lb/cuft, outside expected range"

def test_oil_deno_above_pb():
    """Oil density above Pb should increase with pressure"""
    rho_pb = oil.oil_deno(p=3000, degf=200, rs=500, rsb=500, sg_sp=0.75, pb=3000, api=35)
    rho_above = oil.oil_deno(p=5000, degf=200, rs=500, rsb=500, sg_sp=0.75, pb=3000, api=35)
    assert rho_above > rho_pb, f"Density should increase above Pb: {rho_pb} vs {rho_above}"

def test_oil_viso():
    """Oil viscosity at typical conditions"""
    uo = oil.oil_viso(p=3000, api=35, degf=200, pb=3000, rs=500)
    assert isinstance(uo, float)
    assert 0.1 < uo < 20, f"Oil viscosity = {uo} cP, outside expected range"

def test_oil_viso_below_pb():
    """Oil viscosity should increase below Pb (less gas in solution)"""
    uo_pb = oil.oil_viso(p=3000, api=35, degf=200, pb=3000, rs=500)
    uo_low = oil.oil_viso(p=1000, api=35, degf=200, pb=3000, rs=200)
    assert uo_low > uo_pb, f"Viscosity should increase below Pb: {uo_pb} vs {uo_low}"

def test_oil_viso_above_pb():
    """Oil viscosity should increase above Pb"""
    uo_pb = oil.oil_viso(p=3000, api=35, degf=200, pb=3000, rs=500)
    uo_above = oil.oil_viso(p=5000, api=35, degf=200, pb=3000, rs=500)
    assert uo_above > uo_pb, f"Viscosity should increase above Pb: {uo_pb} vs {uo_above}"

def test_oil_co():
    """Oil compressibility should be positive"""
    co = oil.oil_co(p=3000, api=35, degf=200, sg_sp=0.75, pb=3000, rsb=500)
    assert isinstance(co, float)
    assert co > 0, f"Oil compressibility should be positive, got {co}"
    assert co < 0.001, f"Oil compressibility = {co} seems too high"

def test_oil_rate_radial():
    """Radial oil rate"""
    q = oil.oil_rate_radial(k=100, h=50, pr=3000, pwf=2000, r_w=0.35, r_ext=1000, uo=1.0, bo=1.3)
    assert isinstance(q, (float, np.floating, np.ndarray))
    q_val = float(q) if isinstance(q, np.ndarray) else q
    assert q_val > 0, "Flow rate should be positive"

def test_oil_rate_linear():
    """Linear oil rate"""
    q = oil.oil_rate_linear(k=100, pr=3000, pwf=2000, area=10000, length=5000, uo=1.0, bo=1.3)
    assert isinstance(q, (float, np.floating, np.ndarray))
    q_val = float(q) if isinstance(q, np.ndarray) else q
    assert q_val > 0

def test_oil_twu_props():
    """Twu critical property correlations"""
    sg, tb, tc, pc, vc = oil.oil_twu_props(mw=200, ja=0.2)
    assert 0.7 < sg < 1.0, f"SG = {sg}"
    assert 500 < tb < 1500, f"Tb = {tb} R"
    assert 800 < tc < 2000, f"Tc = {tc} R"
    assert 100 < pc < 500, f"Pc = {pc} psia"
    assert vc > 0, f"Vc = {vc}"

def test_make_bot_og():
    """Black oil table generation"""
    result = oil.make_bot_og(
        pi=4000, api=35, degf=200, sg_g=0.75,
        pmax=5000, pb=3000, rsb=500, nrows=10
    )
    assert 'bot' in result
    assert 'deno' in result
    assert 'deng' in result
    df = result['bot']
    assert len(df) >= 10
    assert 'Pressure (psia)' in df.columns
    assert 'Rs (mscf/stb)' in df.columns
    assert 'Bo (rb/stb)' in df.columns

def test_sg_evolved_gas():
    """Evolved gas SG"""
    sg = oil.sg_evolved_gas(p=2000, degf=200, rsb=500, api=35, sg_sp=0.75)
    assert isinstance(sg, float)
    assert sg >= 0.75, "Evolved gas SG should be >= separator gas SG"

def test_oil_rs_st():
    """Stock tank GOR"""
    rst = oil.oil_rs_st(psp=100, degf_sp=80, api=35)
    assert isinstance(rst, float)
    assert rst >= 0, "Stock tank GOR should be non-negative"

# =============================================================================
# Regression baselines
# =============================================================================

_BASELINES = None

def capture_oil_baselines():
    baselines = {}
    baselines['pb_valmc'] = oil.oil_pbub(api=35, degf=200, rsb=500, sg_sp=0.75, pbmethod='VALMC')
    baselines['pb_velar'] = oil.oil_pbub(api=35, degf=200, rsb=500, sg_sp=0.75, pbmethod='VELAR')
    baselines['rsb_velar'] = oil.oil_rs_bub(api=35, degf=200, pb=3000, sg_sp=0.75, rsmethod='VELAR')
    baselines['bo_3000'] = oil.oil_bo(p=3000, pb=3000, degf=200, rs=500, rsb=500, sg_o=oil.oil_sg(35), sg_g=0.75)
    baselines['deno_3000'] = oil.oil_deno(p=3000, degf=200, rs=500, rsb=500, sg_sp=0.75, api=35)
    baselines['uo_3000'] = oil.oil_viso(p=3000, api=35, degf=200, pb=3000, rs=500)
    return baselines

def get_baselines():
    global _BASELINES
    if _BASELINES is None:
        _BASELINES = capture_oil_baselines()
    return _BASELINES

def test_regression_pb_valmc():
    b = get_baselines()
    pb = oil.oil_pbub(api=35, degf=200, rsb=500, sg_sp=0.75, pbmethod='VALMC')
    assert abs(pb - b['pb_valmc']) / b['pb_valmc'] < 1e-6

def test_regression_bo():
    b = get_baselines()
    bo = oil.oil_bo(p=3000, pb=3000, degf=200, rs=500, rsb=500, sg_o=oil.oil_sg(35), sg_g=0.75)
    assert abs(bo - b['bo_3000']) / b['bo_3000'] < 1e-6

def test_regression_deno():
    b = get_baselines()
    rho = oil.oil_deno(p=3000, degf=200, rs=500, rsb=500, sg_sp=0.75, api=35)
    assert abs(rho - b['deno_3000']) / b['deno_3000'] < 1e-6


if __name__ == '__main__':
    print("=" * 70)
    print("OIL MODULE VALIDATION TESTS")
    print("=" * 70)

    tests = [v for k, v in globals().items() if k.startswith('test_')]
    passed = 0
    failed = 0
    errors = []

    for test in tests:
        try:
            test()
            passed += 1
            print(f"  PASS: {test.__name__}")
        except Exception as e:
            failed += 1
            errors.append((test.__name__, str(e)))
            print(f"  FAIL: {test.__name__}: {e}")

    print(f"\n{'=' * 70}")
    print(f"Results: {passed} passed, {failed} failed out of {passed + failed}")

    if _BASELINES:
        print(f"\nCaptured Baselines:")
        for k, v in _BASELINES.items():
            print(f"  {k}: {v}")

    print("=" * 70)
    sys.exit(1 if failed > 0 else 0)
