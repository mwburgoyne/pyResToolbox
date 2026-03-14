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

# Hard-coded regression baselines - frozen reference values
# If any of these change, it means a code modification has altered computational results.
# Do NOT update these values without explicit review and approval.
# =============================================================================

# Frozen baselines captured 2026-02-27
_FROZEN_BASELINES = {
    'pb_valmc': 2332.9599554806437,
    'pb_velar': 2477.58071878937,
    'bo_3000': 1.286959379509055,
    'deno_3000': 45.143156554713826,
    'uo_3000': 0.5657565566131671,
    'co_4500': 1.6230543665347515e-05,
    'rs_bub_velar': 1872.666133282599,
}

def test_regression_pb_valmc():
    pb = oil.oil_pbub(api=35, degf=200, rsb=500, sg_sp=0.75, pbmethod='VALMC')
    expected = _FROZEN_BASELINES['pb_valmc']
    assert abs(pb - expected) / expected < 1e-6, f"Pb VALMC changed: {pb} vs frozen {expected}"

def test_regression_bo():
    bo = oil.oil_bo(p=3000, pb=3000, degf=200, rs=500, rsb=500, sg_o=oil.oil_sg(35), sg_g=0.75)
    expected = _FROZEN_BASELINES['bo_3000']
    assert abs(bo - expected) / expected < 1e-6, f"Bo changed: {bo} vs frozen {expected}"

def test_regression_deno():
    rho = oil.oil_deno(p=3000, degf=200, rs=500, rsb=500, sg_sp=0.75, api=35)
    expected = _FROZEN_BASELINES['deno_3000']
    assert abs(rho - expected) / expected < 1e-6, f"Oil density changed: {rho} vs frozen {expected}"

def test_regression_viso():
    """Oil viscosity must match frozen baseline"""
    uo = oil.oil_viso(p=3000, api=35, degf=200, pb=3000, rs=500)
    expected = _FROZEN_BASELINES['uo_3000']
    assert abs(uo - expected) / expected < 1e-6, f"Oil viscosity changed: {uo} vs frozen {expected}"

def test_regression_co():
    """Oil undersaturated compressibility must match frozen baseline"""
    co = oil.oil_co(p=4500, api=47, degf=180, sg_sp=0.72, rsb=2750)
    expected = _FROZEN_BASELINES['co_4500']
    assert abs(co - expected) / expected < 1e-4, f"Oil Co changed: {co} vs frozen {expected}"

def test_regression_pb_velar():
    """Velarde Pb must match frozen baseline"""
    pb = oil.oil_pbub(api=35, degf=200, rsb=500, sg_sp=0.75, pbmethod='VELAR')
    expected = _FROZEN_BASELINES['pb_velar']
    assert abs(pb - expected) / expected < 1e-6, f"Pb VELAR changed: {pb} vs frozen {expected}"

# =============================================================================
# oil_harmonize and vis_frac tests
# =============================================================================

def test_oil_harmonize_returns_4_tuple():
    """oil_harmonize returns (pb, rsb, rsb_frac, vis_frac)"""
    result = oil.oil_harmonize(pb=3500, degf=175, api=38, sg_g=0.68)
    assert len(result) == 4, f"Expected 4-tuple, got {len(result)}"
    pb, rsb, rsb_frac, vis_frac = result
    assert pb == 3500
    assert rsb > 0
    assert rsb_frac == 1.0
    assert vis_frac == 1.0

def test_oil_harmonize_vis_frac_no_target():
    """vis_frac == 1.0 when uo_target is not specified"""
    _pb, _rsb, _frac, vis_frac = oil.oil_harmonize(
        pb=3500, rsb=1200, degf=175, api=38, sg_sp=0.68, sg_g=0.68
    )
    assert vis_frac == 1.0

def test_oil_harmonize_vis_frac_with_target():
    """vis_frac should scale viscosity to match target"""
    pb, rsb, rsb_frac, vis_frac = oil.oil_harmonize(
        pb=3000, rsb=500, degf=200, api=35, sg_sp=0.75, sg_g=0.75,
        uo_target=1.0, p_uo=3000
    )
    # Verify: oil_viso at p_uo * vis_frac should equal uo_target
    rs_at_p = oil.oil_rs(api=35, degf=200, sg_sp=0.75, p=3000, pb=pb,
                         rsb=rsb / rsb_frac) * rsb_frac
    uo_corr = oil.oil_viso(p=3000, api=35, degf=200, pb=pb, rs=rs_at_p)
    assert abs(uo_corr * vis_frac - 1.0) < 1e-10, \
        f"Scaled viscosity should equal target 1.0, got {uo_corr * vis_frac}"

def test_oil_harmonize_pb_rsb_compat():
    """Deprecated oil_harmonize_pb_rsb returns 3-tuple"""
    result = oil.oil_harmonize_pb_rsb(pb=3500, degf=175, api=38, sg_g=0.68)
    assert len(result) == 3, f"Expected 3-tuple, got {len(result)}"
    pb, rsb, rsb_frac = result
    assert pb == 3500
    assert rsb > 0
    assert rsb_frac == 1.0

def test_oilpvt_vis_frac():
    """OilPVT(vis_frac=2.0).viscosity() should be 2x unscaled"""
    opvt_base = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    opvt_scaled = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500, vis_frac=2.0)
    uo_base = opvt_base.viscosity(2000, 180)
    uo_scaled = opvt_scaled.viscosity(2000, 180)
    assert abs(uo_scaled - 2.0 * uo_base) < 1e-10, \
        f"vis_frac=2.0 should double viscosity: {uo_scaled} vs {2.0 * uo_base}"

def test_oilpvt_vis_frac_default():
    """Default vis_frac=1.0 leaves results unchanged"""
    opvt = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    assert opvt.vis_frac == 1.0
    uo = opvt.viscosity(2000, 180)
    uo_direct = oil.oil_viso(p=2000, api=35, degf=180, pb=2500,
                             rs=oil.oil_rs(api=35, degf=180, sg_sp=0.65, p=2000,
                                           pb=2500, rsb=500))
    assert abs(uo - uo_direct) < 1e-10

def test_oilpvt_vis_frac_validation():
    """OilPVT should reject vis_frac <= 0"""
    try:
        oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500, vis_frac=0)
        assert False, "Should have raised ValueError"
    except ValueError:
        pass
    try:
        oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500, vis_frac=-1.0)
        assert False, "Should have raised ValueError"
    except ValueError:
        pass

# =============================================================================
# OilPVT rsb_frac and from_harmonize tests
# =============================================================================

def test_oilpvt_rsb_frac_default():
    """Default rsb_frac=1.0 leaves rs unchanged"""
    opvt = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    assert opvt.rsb_frac == 1.0
    rs = opvt.rs(2000, 180)
    rs_direct = oil.oil_rs(api=35, degf=180, sg_sp=0.65, p=2000, pb=2500, rsb=500)
    assert abs(rs - rs_direct) < 1e-10

def test_oilpvt_rsb_frac_stored():
    """rsb_frac is stored and accessible on OilPVT"""
    opvt = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500, rsb_frac=0.85)
    assert opvt.rsb_frac == 0.85
    # Rs should still be computable
    rs = opvt.rs(2000, 180)
    assert rs > 0

def test_oilpvt_rsb_frac_at_pb():
    """With rsb_frac, Rs at Pb should still equal rsb (the physical value)"""
    rsb = 500
    rsb_frac = 0.85
    opvt = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=rsb, rsb_frac=rsb_frac)
    rs_at_pb = opvt.rs(2500, 180)
    # Rs at Pb: oil_rs returns rsb/rsb_frac (which is rsb at the unscaled pb),
    # then multiplied by rsb_frac, giving back rsb
    assert abs(rs_at_pb - rsb) < 1.0, \
        f"Rs at Pb should be ~rsb={rsb}, got {rs_at_pb}"

def test_oilpvt_rsb_frac_validation():
    """OilPVT should reject rsb_frac <= 0"""
    try:
        oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500, rsb_frac=0)
        assert False, "Should have raised ValueError"
    except ValueError:
        pass

def test_oilpvt_from_harmonize_pb_only():
    """from_harmonize with pb only should compute rsb"""
    opvt = oil.OilPVT.from_harmonize(degf=200, api=35, sg_g=0.75, pb=3000)
    assert opvt.pb == 3000
    assert opvt.rsb > 0
    assert opvt.rsb_frac == 1.0
    assert opvt.vis_frac == 1.0
    # Rs at Pb should equal rsb
    rs_at_pb = opvt.rs(3000, 200)
    assert abs(rs_at_pb - opvt.rsb) < 1.0

def test_oilpvt_from_harmonize_both():
    """from_harmonize with both pb and rsb should set rsb_frac"""
    opvt = oil.OilPVT.from_harmonize(degf=200, api=35, sg_sp=0.75, sg_g=0.75,
                                      pb=3000, rsb=500)
    assert opvt.pb == 3000
    assert opvt.rsb == 500
    # rsb_frac may or may not be 1.0 depending on correlation match

def test_oilpvt_from_harmonize_with_viscosity():
    """from_harmonize with uo_target should set vis_frac"""
    opvt = oil.OilPVT.from_harmonize(degf=200, api=35, sg_sp=0.75, sg_g=0.75,
                                      pb=3000, rsb=500, uo_target=1.0, p_uo=3000)
    assert opvt.vis_frac != 1.0
    # Viscosity at p_uo should be close to target
    uo = opvt.viscosity(3000, 200)
    assert abs(uo - 1.0) < 0.01, f"Viscosity should be ~1.0, got {uo}"


# =============================================================================
# OilPVT auto-harmonization tests
# =============================================================================

def test_oilpvt_auto_harmonize_pb_only():
    """OilPVT with degf and no rsb should calculate rsb from correlation"""
    opvt = oil.OilPVT(api=35, sg_sp=0.75, pb=3000, degf=200)
    assert opvt.pb == 3000
    assert opvt.rsb > 0
    assert opvt.rsb_frac == 1.0
    assert opvt.vis_frac == 1.0
    # Rs at Pb should equal rsb
    rs_at_pb = opvt.rs(3000, 200)
    assert abs(rs_at_pb - opvt.rsb) < 1.0

def test_oilpvt_auto_harmonize_both():
    """OilPVT with degf and both pb and rsb should compute rsb_frac"""
    opvt = oil.OilPVT(api=35, sg_sp=0.75, pb=3000, rsb=500, degf=200)
    assert opvt.pb == 3000
    assert opvt.rsb == 500

def test_oilpvt_auto_harmonize_with_vis():
    """OilPVT with degf, uo_target, and p_uo should compute vis_frac"""
    opvt = oil.OilPVT(api=35, sg_sp=0.75, pb=3000, rsb=500, degf=200,
                      uo_target=1.0, p_uo=3000)
    assert opvt.vis_frac != 1.0
    uo = opvt.viscosity(3000, 200)
    assert abs(uo - 1.0) < 0.01, f"Viscosity should be ~1.0, got {uo}"

def test_oilpvt_no_rsb_no_degf_raises():
    """OilPVT with rsb=0 and degf=0 should raise ValueError"""
    try:
        oil.OilPVT(api=35, sg_sp=0.75, pb=3000)
        assert False, "Should have raised ValueError"
    except ValueError:
        pass

def test_oilpvt_legacy_with_rsb():
    """OilPVT with rsb>0 and no degf should work as before"""
    opvt = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    assert opvt.pb == 2500
    assert opvt.rsb == 500
    rs = opvt.rs(2000, 180)
    assert rs > 0

def test_oilpvt_auto_harmonize_matches_from_harmonize():
    """OilPVT(degf=...) should produce same result as from_harmonize"""
    opvt1 = oil.OilPVT(api=35, sg_sp=0.75, pb=3000, degf=200, sg_g=0.75)
    opvt2 = oil.OilPVT.from_harmonize(degf=200, api=35, sg_sp=0.75, sg_g=0.75, pb=3000)
    assert abs(opvt1.rsb - opvt2.rsb) < 1e-6
    assert abs(opvt1.rsb_frac - opvt2.rsb_frac) < 1e-6
    assert abs(opvt1.vis_frac - opvt2.vis_frac) < 1e-6

# =============================================================================
# oil_rate_radial/linear with oil_pvt tests
# =============================================================================

def test_oil_rate_radial_with_pvt():
    """oil_rate_radial with oil_pvt should produce valid results"""
    opvt = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    q = oil.oil_rate_radial(k=100, h=50, pr=3000, pwf=2000, r_w=0.35, r_ext=1000,
                            oil_pvt=opvt, degf=180)
    q_val = float(q) if isinstance(q, np.ndarray) else q
    assert q_val > 0, f"Rate should be positive, got {q_val}"

def test_oil_rate_radial_pvt_no_degf_raises():
    """oil_rate_radial with oil_pvt but no degf should raise"""
    opvt = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    try:
        oil.oil_rate_radial(k=100, h=50, pr=3000, pwf=2000, r_w=0.35, r_ext=1000,
                            oil_pvt=opvt)
        assert False, "Should have raised ValueError"
    except ValueError:
        pass

def test_oil_rate_linear_with_pvt():
    """oil_rate_linear with oil_pvt should produce valid results"""
    opvt = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    q = oil.oil_rate_linear(k=100, pr=3000, pwf=2000, area=10000, length=5000,
                            oil_pvt=opvt, degf=180)
    q_val = float(q) if isinstance(q, np.ndarray) else q
    assert q_val > 0, f"Rate should be positive, got {q_val}"


# =============================================================================
# oil_co co_sat tests
# =============================================================================

def test_oil_co_sat_false_returns_float():
    """Default co_sat=False returns a scalar float, not a list"""
    co = oil.oil_co(p=3000, api=35, degf=200, sg_sp=0.75, pb=3500, rsb=500)
    assert isinstance(co, (float, np.floating)), f"Expected float, got {type(co)}"

def test_oil_co_sat_true_returns_list():
    """co_sat=True returns a 2-element list"""
    result = oil.oil_co(p=3000, api=35, degf=200, sg_sp=0.75, pb=3500, rsb=500, co_sat=True)
    assert isinstance(result, list), f"Expected list, got {type(result)}"
    assert len(result) == 2, f"Expected 2 elements, got {len(result)}"

def test_oil_co_sat_below_pb():
    """Below Pb, co_sat > co_usat (gas evolution adds compressibility)"""
    result = oil.oil_co(p=2000, api=47, degf=180, sg_sp=0.72, rsb=2750, pb=4945, co_sat=True)
    co_usat, co_sat_val = result
    assert co_usat > 0, f"co_usat should be positive, got {co_usat}"
    assert co_sat_val > 0, f"co_sat should be positive, got {co_sat_val}"
    assert co_sat_val > co_usat, f"Below Pb, co_sat ({co_sat_val}) should exceed co_usat ({co_usat})"

def test_oil_co_sat_above_pb():
    """Above Pb, co_sat == co_usat (no gas evolution)"""
    result = oil.oil_co(p=6000, api=47, degf=180, sg_sp=0.72, rsb=2750, pb=4945, co_sat=True)
    co_usat, co_sat_val = result
    assert abs(co_usat - co_sat_val) < 1e-15, f"Above Pb co_usat ({co_usat}) should equal co_sat ({co_sat_val})"

def test_oil_co_sat_false_unchanged():
    """co_sat=False must return the same value as before (backward compat with frozen baseline)"""
    co = oil.oil_co(p=4500, api=47, degf=180, sg_sp=0.72, rsb=2750)
    expected = _FROZEN_BASELINES['co_4500']
    assert abs(co - expected) / expected < 1e-4, f"co_sat=False changed: {co} vs frozen {expected}"

def test_oil_co_sat_usat_matches_default():
    """co_usat from co_sat=True matches the scalar from co_sat=False"""
    co_scalar = oil.oil_co(p=2000, api=47, degf=180, sg_sp=0.72, rsb=2750, pb=4945)
    co_list = oil.oil_co(p=2000, api=47, degf=180, sg_sp=0.72, rsb=2750, pb=4945, co_sat=True)
    assert abs(co_scalar - co_list[0]) < 1e-15, f"co_usat mismatch: {co_scalar} vs {co_list[0]}"


# =============================================================================
# oil_bt tests
# =============================================================================

def test_oil_bt_above_pb():
    """Above Pb, Bt should equal Bo (no free gas)"""
    bt = oil.oil_bt(p=6000, api=47, degf=180, sg_sp=0.72, rsb=2750, pb=4945)
    # Compute Bo for comparison
    sg_g, sg_sp2 = oil.check_sgs(sg_g=0, sg_sp=0.72)
    rs = oil.oil_rs(api=47, degf=180, sg_sp=0.72, p=6000, pb=4945, rsb=2750)
    bo = oil.oil_bo(p=6000, pb=4945, degf=180, rs=rs, rsb=2750, sg_sp=sg_sp2, sg_g=sg_g, sg_o=oil.oil_sg(47))
    assert abs(bt - bo) < 1e-10, f"Above Pb, Bt ({bt}) should equal Bo ({bo})"

def test_oil_bt_below_pb():
    """Below Pb, Bt > Bo (includes free gas volume)"""
    bt = oil.oil_bt(p=2000, api=47, degf=180, sg_sp=0.72, rsb=2750, pb=4945)
    sg_g, sg_sp2 = oil.check_sgs(sg_g=0, sg_sp=0.72)
    rs = oil.oil_rs(api=47, degf=180, sg_sp=0.72, p=2000, pb=4945, rsb=2750)
    bo = oil.oil_bo(p=2000, pb=4945, degf=180, rs=rs, rsb=2750, sg_sp=sg_sp2, sg_g=sg_g, sg_o=oil.oil_sg(47))
    assert bt > bo, f"Below Pb, Bt ({bt}) should exceed Bo ({bo})"

def test_oil_bt_positive():
    """Bt should always be positive"""
    bt = oil.oil_bt(p=1000, api=35, degf=200, sg_sp=0.75, rsb=500, pb=3000)
    assert bt > 0, f"Bt should be positive, got {bt}"

def test_oil_bt_increases_as_p_drops():
    """Bt should increase as pressure drops below Pb (more gas evolves)"""
    bt_3000 = oil.oil_bt(p=3000, api=47, degf=180, sg_sp=0.72, rsb=2750, pb=4945)
    bt_2000 = oil.oil_bt(p=2000, api=47, degf=180, sg_sp=0.72, rsb=2750, pb=4945)
    bt_1000 = oil.oil_bt(p=1000, api=47, degf=180, sg_sp=0.72, rsb=2750, pb=4945)
    assert bt_2000 > bt_3000, f"Bt at 2000 ({bt_2000}) should exceed Bt at 3000 ({bt_3000})"
    assert bt_1000 > bt_2000, f"Bt at 1000 ({bt_1000}) should exceed Bt at 2000 ({bt_2000})"


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

    print("=" * 70)
    sys.exit(1 if failed > 0 else 0)
