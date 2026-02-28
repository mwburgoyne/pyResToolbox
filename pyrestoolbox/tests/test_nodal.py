#!/usr/bin/env python3
"""
Tests for pyRestoolbox nodal module (VLP/IPR/Nodal Analysis).
"""

import sys
import os
import math

# Ensure project parent is on path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(project_root))

from pyrestoolbox.gas import GasPVT, gas_z, gas_ug, gas_den, gas_bg
from pyrestoolbox.oil import OilPVT, oil_rs, oil_bo, oil_deno, oil_viso
from pyrestoolbox.nodal import (
    WellSegment, Completion, Reservoir, fbhp, outflow_curve, ipr_curve, operating_point
)


# ============================================================================
#  GasPVT Tests
# ============================================================================

def test_gas_pvt_z_matches_direct():
    """GasPVT.z() must match gas_z() for same inputs."""
    gpvt = GasPVT(sg=0.75, co2=0.01, h2s=0.02)
    z_obj = gpvt.z(2000, 200)
    z_direct = gas_z(2000, 0.75, 200, co2=0.01, h2s=0.02)
    assert abs(z_obj - z_direct) < 1e-10, f"Z mismatch: {z_obj} vs {z_direct}"


def test_gas_pvt_viscosity_matches_direct():
    """GasPVT.viscosity() must match gas_ug() for same inputs."""
    gpvt = GasPVT(sg=0.65)
    vis_obj = gpvt.viscosity(3000, 250)
    vis_direct = gas_ug(3000, 0.65, 250)
    assert abs(vis_obj - vis_direct) < 1e-10, f"Viscosity mismatch: {vis_obj} vs {vis_direct}"


def test_gas_pvt_density_matches_direct():
    """GasPVT.density() must match gas_den() for same inputs."""
    gpvt = GasPVT(sg=0.7)
    den_obj = gpvt.density(1500, 180)
    den_direct = gas_den(1500, 0.7, 180)
    assert abs(den_obj - den_direct) < 1e-10, f"Density mismatch: {den_obj} vs {den_direct}"


def test_gas_pvt_bg_matches_direct():
    """GasPVT.bg() must match gas_bg() for same inputs."""
    gpvt = GasPVT(sg=0.65)
    bg_obj = gpvt.bg(2000, 200)
    bg_direct = gas_bg(2000, 0.65, 200)
    assert abs(bg_obj - bg_direct) < 1e-10, f"Bg mismatch: {bg_obj} vs {bg_direct}"


def test_gas_pvt_h2_auto_selects_bns():
    """GasPVT with h2>0 should auto-select BNS method."""
    gpvt = GasPVT(sg=0.5, h2=0.1)
    assert gpvt.zmethod.name in ('BNS', 'BUR'), f"Expected BNS, got {gpvt.zmethod.name}"
    assert gpvt.cmethod.name in ('BNS', 'BUR'), f"Expected BNS, got {gpvt.cmethod.name}"


def test_gas_pvt_precomputed_tc_pc():
    """GasPVT should pre-compute tc and pc in constructor."""
    gpvt = GasPVT(sg=0.75)
    assert gpvt.tc > 0, "tc should be positive"
    assert gpvt.pc > 0, "pc should be positive"


# ============================================================================
#  OilPVT Tests
# ============================================================================

def test_oil_pvt_rs_matches_direct():
    """OilPVT.rs() must match oil_rs() for same inputs."""
    opvt = OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    rs_obj = opvt.rs(2000, 180)
    rs_direct = oil_rs(api=35, degf=180, sg_sp=0.65, p=2000, pb=2500, rsb=500)
    assert abs(rs_obj - rs_direct) < 1e-6, f"Rs mismatch: {rs_obj} vs {rs_direct}"


def test_oil_pvt_bo_matches_direct():
    """OilPVT.bo() must match oil_bo() for same inputs."""
    opvt = OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    p, degf = 2000, 180
    rs = opvt.rs(p, degf)
    bo_obj = opvt.bo(p, degf, rs=rs)
    bo_direct = oil_bo(p=p, pb=2500, degf=degf, rs=rs, rsb=500,
                       sg_o=opvt.sg_o, sg_sp=0.65)
    assert abs(bo_obj - bo_direct) < 1e-6, f"Bo mismatch: {bo_obj} vs {bo_direct}"


def test_oil_pvt_viscosity_matches_direct():
    """OilPVT.viscosity() must match oil_viso() for same inputs."""
    opvt = OilPVT(api=30, sg_sp=0.7, pb=3000, rsb=600)
    p, degf = 2500, 200
    rs = opvt.rs(p, degf)
    vis_obj = opvt.viscosity(p, degf, rs=rs)
    vis_direct = oil_viso(p=p, api=30, degf=degf, pb=3000, rs=rs)
    assert abs(vis_obj - vis_direct) < 1e-8, f"Viscosity mismatch: {vis_obj} vs {vis_direct}"


def test_oil_pvt_above_pb():
    """OilPVT.rs() above Pb should return rsb."""
    opvt = OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    rs = opvt.rs(3000, 200)
    assert abs(rs - 500) < 1e-6, f"Rs above Pb should equal rsb: {rs}"


def test_oil_pvt_sg_o():
    """OilPVT should correctly compute sg_o from API."""
    opvt = OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    expected = 141.5 / (35 + 131.5)
    assert abs(opvt.sg_o - expected) < 1e-10, f"sg_o mismatch: {opvt.sg_o} vs {expected}"


# ============================================================================
#  Completion Tests
# ============================================================================

def test_completion_basic():
    """Basic Completion properties."""
    c = Completion(tid=2.441, length=10000, tht=100, bht=200)
    assert c.tid == 2.441
    assert c.length == 10000
    assert c.mpd == 10000  # Default mpd = length
    assert not c.has_casing_section
    assert c.casing_length == 0


def test_completion_with_casing():
    """Completion with casing section below tubing."""
    c = Completion(tid=2.441, length=8000, tht=100, bht=250,
                   cid=4.892, crough=0.001, mpd=10000)
    assert c.has_casing_section
    assert c.casing_length == 2000
    expected_temp = 100 + (250 - 100) * 8000 / 10000
    assert abs(c.tubing_end_temperature - expected_temp) < 0.01


def test_completion_no_casing_if_mpd_equals_length():
    """No casing section when mpd == length even if cid > 0."""
    c = Completion(tid=2.441, length=10000, tht=100, bht=200,
                   cid=4.892, mpd=10000)
    assert not c.has_casing_section


# ============================================================================
#  Reservoir Tests
# ============================================================================

def test_reservoir_basic():
    """Basic Reservoir properties."""
    r = Reservoir(pr=3000, degf=200, k=10, h=50, re=1500, rw=0.35)
    assert r.pr == 3000
    assert r.S == 0
    assert r.D == 0


# ============================================================================
#  VLP Gas Well Tests
# ============================================================================

_GAS_COMP = Completion(tid=2.441, length=10000, tht=100, bht=200)
_GAS_PARAMS = dict(
    thp=500, completion=_GAS_COMP, well_type='gas',
    qg_mmscfd=5.0, cgr=10, qw_bwpd=10, oil_vis=1.0, api=45, gsg=0.65
)


def test_vlp_gas_hb_bhp_gt_thp():
    """HB gas: BHP should exceed THP for production."""
    bhp = fbhp(vlpmethod='HB', **_GAS_PARAMS)
    assert bhp > 500, f"BHP {bhp} should exceed THP 500"


def test_vlp_gas_wg_bhp_gt_thp():
    """WG gas: BHP should exceed THP for production."""
    bhp = fbhp(vlpmethod='WG', **_GAS_PARAMS)
    assert bhp > 500, f"BHP {bhp} should exceed THP 500"


def test_vlp_gas_gray_bhp_gt_thp():
    """Gray gas: BHP should exceed THP for production."""
    bhp = fbhp(vlpmethod='GRAY', **_GAS_PARAMS)
    assert bhp > 500, f"BHP {bhp} should exceed THP 500"


def test_vlp_gas_bb_bhp_gt_thp():
    """BB gas: BHP should exceed THP for production."""
    bhp = fbhp(vlpmethod='BB', **_GAS_PARAMS)
    assert bhp > 500, f"BHP {bhp} should exceed THP 500"


def test_vlp_gas_bhp_increases_with_rate():
    """Gas BHP should increase with increasing rate (friction effect)."""
    bhp_low = fbhp(vlpmethod='HB', **{**_GAS_PARAMS, 'qg_mmscfd': 2.0})
    bhp_high = fbhp(vlpmethod='HB', **{**_GAS_PARAMS, 'qg_mmscfd': 10.0})
    assert bhp_high > bhp_low, f"BHP should increase with rate: {bhp_low} vs {bhp_high}"


def test_vlp_gas_zero_rate_static_column():
    """Zero gas rate should return static column pressure."""
    c = Completion(tid=2.441, length=10000, tht=100, bht=200)
    bhp = fbhp(vlpmethod='HB', thp=500, completion=c, well_type='gas',
               qg_mmscfd=0.0, cgr=0, qw_bwpd=0, oil_vis=1.0, api=45, gsg=0.65)
    # Static column should give BHP > THP
    assert bhp > 500, f"Static BHP {bhp} should exceed THP 500"
    # Gas column ~ 0.01-0.03 psi/ft * 10000 ft + 500 ~ 600-800 psia
    assert bhp < 1000, f"Static gas BHP {bhp} too high"


def test_vlp_gas_methods_within_range():
    """All 4 gas VLP methods should produce BHP within reasonable range of each other."""
    bhps = []
    for method in ['HB', 'WG', 'GRAY', 'BB']:
        bhps.append(fbhp(vlpmethod=method, **_GAS_PARAMS))
    # All should be within 500-3000 for this case
    for bhp in bhps:
        assert 500 < bhp < 3000, f"Gas BHP {bhp} outside reasonable range"


# ============================================================================
#  VLP Oil Well Tests
# ============================================================================

_OIL_COMP = Completion(tid=2.441, length=8000, tht=100, bht=180)
_OIL_PARAMS = dict(
    thp=200, completion=_OIL_COMP, well_type='oil',
    qt_stbpd=1000, gor=800, wc=0.3, pb=2500, rsb=500,
    sgsp=0.65, gsg=0.65, api=35
)


def test_vlp_oil_hb_bhp_gt_thp():
    """HB oil: BHP should exceed THP for production."""
    bhp = fbhp(vlpmethod='HB', **_OIL_PARAMS)
    assert bhp > 200, f"BHP {bhp} should exceed THP 200"


def test_vlp_oil_wg_bhp_gt_thp():
    """WG oil: BHP should exceed THP."""
    bhp = fbhp(vlpmethod='WG', **_OIL_PARAMS)
    assert bhp > 200, f"BHP {bhp} should exceed THP 200"


def test_vlp_oil_gray_bhp_gt_thp():
    """Gray oil: BHP should exceed THP."""
    bhp = fbhp(vlpmethod='GRAY', **_OIL_PARAMS)
    assert bhp > 200, f"BHP {bhp} should exceed THP 200"


def test_vlp_oil_bb_bhp_gt_thp():
    """BB oil: BHP should exceed THP."""
    bhp = fbhp(vlpmethod='BB', **_OIL_PARAMS)
    assert bhp > 200, f"BHP {bhp} should exceed THP 200"


def test_vlp_oil_bhp_increases_with_rate():
    """Oil BHP should increase with increasing rate."""
    bhp_low = fbhp(vlpmethod='HB', **{**_OIL_PARAMS, 'qt_stbpd': 500})
    bhp_high = fbhp(vlpmethod='HB', **{**_OIL_PARAMS, 'qt_stbpd': 3000})
    assert bhp_high > bhp_low, f"BHP should increase with rate: {bhp_low} vs {bhp_high}"


def test_vlp_oil_zero_rate_static_column():
    """Zero oil rate should return static column pressure."""
    bhp = fbhp(vlpmethod='HB', **{**_OIL_PARAMS, 'qt_stbpd': 0.0})
    assert bhp > 200, f"Static BHP {bhp} should exceed THP 200"


def test_vlp_oil_methods_within_range():
    """All 4 oil VLP methods should produce BHP within reasonable range."""
    bhps = []
    for method in ['HB', 'WG', 'GRAY', 'BB']:
        bhps.append(fbhp(vlpmethod=method, **_OIL_PARAMS))
    for bhp in bhps:
        assert 200 < bhp < 5000, f"Oil BHP {bhp} outside reasonable range"


# ============================================================================
#  Two-Section (Casing) Tests
# ============================================================================

def test_vlp_casing_section_increases_bhp():
    """Adding a casing section below tubing should increase BHP."""
    comp_tubing_only = Completion(tid=2.441, length=8000, tht=100, bht=200)
    comp_with_casing = Completion(tid=2.441, length=8000, tht=100, bht=250,
                                  cid=4.892, mpd=10000)
    bhp_tubing = fbhp(thp=500, completion=comp_tubing_only, vlpmethod='HB',
                      well_type='gas', qg_mmscfd=5.0, gsg=0.65, cgr=10,
                      qw_bwpd=10, api=45, oil_vis=1.0)
    bhp_casing = fbhp(thp=500, completion=comp_with_casing, vlpmethod='HB',
                      well_type='gas', qg_mmscfd=5.0, gsg=0.65, cgr=10,
                      qw_bwpd=10, api=45, oil_vis=1.0)
    assert bhp_casing > bhp_tubing, f"Casing BHP {bhp_casing} should exceed tubing-only {bhp_tubing}"


# ============================================================================
#  Outflow Curve Tests
# ============================================================================

def test_outflow_curve_gas_returns_dict():
    """outflow_curve should return dict with rates and bhp."""
    c = Completion(tid=2.441, length=10000, tht=100, bht=200)
    result = outflow_curve(thp=500, completion=c, vlpmethod='HB',
                           well_type='gas', n_rates=5, max_rate=20.0,
                           gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    assert 'rates' in result
    assert 'bhp' in result
    assert len(result['rates']) == 5
    assert len(result['bhp']) == 5


def test_outflow_curve_gas_bhp_increasing():
    """Gas outflow curve BHP should generally increase with rate at moderate rates."""
    c = Completion(tid=2.441, length=10000, tht=100, bht=200)
    rates = [2.0, 5.0, 10.0, 15.0, 20.0]
    result = outflow_curve(thp=500, completion=c, vlpmethod='HB',
                           well_type='gas', rates=rates,
                           gsg=0.65, cgr=0, qw_bwpd=0, api=45, oil_vis=1.0)
    # BHP should increase monotonically at moderate-to-high rates
    for i in range(1, len(result['bhp'])):
        assert result['bhp'][i] > result['bhp'][i-1], \
            f"BHP should increase: rate {rates[i]}: {result['bhp'][i]} <= {result['bhp'][i-1]}"


def test_outflow_curve_oil_returns_dict():
    """outflow_curve for oil should return correct structure."""
    c = Completion(tid=2.441, length=8000, tht=100, bht=180)
    result = outflow_curve(thp=200, completion=c, vlpmethod='HB',
                           well_type='oil', n_rates=5, max_rate=5000.0,
                           gsg=0.65, gor=800, wc=0.3, pb=2500, rsb=500,
                           sgsp=0.65, api=35)
    assert len(result['rates']) == 5
    assert len(result['bhp']) == 5


# ============================================================================
#  IPR Curve Tests
# ============================================================================

def test_ipr_gas_rate_decreases_with_pwf():
    """Gas IPR: rate should decrease as pwf increases toward Pr."""
    r = Reservoir(pr=3000, degf=200, k=10, h=50, re=1500, rw=0.35)
    gpvt = GasPVT(sg=0.65)
    result = ipr_curve(reservoir=r, well_type='gas', gas_pvt=gpvt, n_points=10)
    assert 'pwf' in result
    assert 'rate' in result
    # At Pr, rate should be ~0
    assert abs(result['rate'][-1]) < 1.0, f"Rate at Pr should be ~0, got {result['rate'][-1]}"
    # Rate should be highest at lowest Pwf
    assert result['rate'][0] > result['rate'][-1], "Rate should be highest at min Pwf"


def test_ipr_oil_darcy():
    """Oil IPR without Pb info should give linear Darcy IPR."""
    r = Reservoir(pr=3000, degf=200, k=50, h=30, re=1000, rw=0.35)
    result = ipr_curve(reservoir=r, well_type='oil', n_points=5, bo=1.2, uo=2.0)
    # Linear: rate should be proportional to (Pr - Pwf)
    assert result['rate'][-1] < 1.0, f"Rate at Pr should be ~0"
    assert result['rate'][0] > 0, "Rate at min Pwf should be positive"


def test_ipr_oil_with_vogel():
    """Oil IPR with OilPVT and Pb should use Vogel below Pb."""
    r = Reservoir(pr=2000, degf=180, k=50, h=30, re=1000, rw=0.35)
    opvt = OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    result = ipr_curve(reservoir=r, well_type='oil', oil_pvt=opvt, n_points=10)
    # Since Pr < Pb, Vogel applies - rate should be positive at min Pwf
    assert result['rate'][0] > 0, "Rate at min Pwf should be positive"


def test_ipr_gas_with_default_sg():
    """Gas IPR should work without explicit GasPVT."""
    r = Reservoir(pr=3000, degf=200, k=10, h=50, re=1500, rw=0.35)
    result = ipr_curve(reservoir=r, well_type='gas', n_points=5, gsg=0.65)
    assert len(result['rate']) == 5


# ============================================================================
#  Operating Point Tests
# ============================================================================

def test_operating_point_gas():
    """Gas operating point should return valid intersection."""
    c = Completion(tid=2.441, length=10000, tht=100, bht=200)
    r = Reservoir(pr=3000, degf=200, k=10, h=50, re=1500, rw=0.35)
    gpvt = GasPVT(sg=0.65)
    result = operating_point(thp=500, completion=c, reservoir=r,
                             vlpmethod='HB', well_type='gas',
                             gas_pvt=gpvt, gsg=0.65, cgr=0,
                             qw_bwpd=0, api=45, oil_vis=1.0)
    assert 'rate' in result
    assert 'bhp' in result
    assert result['rate'] > 0, f"Operating rate should be positive: {result['rate']}"
    assert result['bhp'] > 500, f"Operating BHP should exceed THP: {result['bhp']}"
    assert result['bhp'] < 3000, f"Operating BHP should be less than Pr: {result['bhp']}"


def test_operating_point_oil():
    """Oil operating point should return valid intersection."""
    c = Completion(tid=2.441, length=8000, tht=100, bht=180)
    r = Reservoir(pr=3000, degf=180, k=50, h=30, re=1000, rw=0.35)
    opvt = OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    result = operating_point(thp=200, completion=c, reservoir=r,
                             vlpmethod='HB', well_type='oil',
                             oil_pvt=opvt, gor=800, wc=0.3, gsg=0.65)
    assert result['rate'] > 0, f"Operating rate should be positive: {result['rate']}"
    assert result['bhp'] > 200, f"Operating BHP should exceed THP: {result['bhp']}"


def test_operating_point_returns_curves():
    """Operating point should include VLP and IPR curves."""
    c = Completion(tid=2.441, length=10000, tht=100, bht=200)
    r = Reservoir(pr=3000, degf=200, k=10, h=50, re=1500, rw=0.35)
    result = operating_point(thp=500, completion=c, reservoir=r,
                             vlpmethod='HB', well_type='gas',
                             gsg=0.65, cgr=0, qw_bwpd=0, api=45, oil_vis=1.0)
    assert 'vlp' in result
    assert 'ipr' in result
    assert len(result['vlp']['rates']) > 0
    assert len(result['ipr']['pwf']) > 0


# ============================================================================
#  Injection Test
# ============================================================================

def test_vlp_gas_injection():
    """Injection should give lower BHP than production at same rate."""
    c = Completion(tid=2.441, length=5000, tht=100, bht=150)
    bhp_prod = fbhp(thp=1000, completion=c, vlpmethod='HB', well_type='gas',
                    qg_mmscfd=2.0, gsg=0.65, cgr=0, qw_bwpd=0,
                    api=45, oil_vis=1.0, injection=False)
    bhp_inj = fbhp(thp=1000, completion=c, vlpmethod='HB', well_type='gas',
                   qg_mmscfd=2.0, gsg=0.65, cgr=0, qw_bwpd=0,
                   api=45, oil_vis=1.0, injection=True)
    # Injection friction opposes gravity, so BHP should be lower
    assert bhp_inj < bhp_prod, f"Injection BHP {bhp_inj} should be less than production {bhp_prod}"


# ============================================================================
#  OilPVT integration with fbhp
# ============================================================================

def test_fbhp_oil_with_oil_pvt_object():
    """fbhp() should accept OilPVT and extract parameters."""
    c = Completion(tid=2.441, length=8000, tht=100, bht=180)
    opvt = OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    bhp = fbhp(thp=200, completion=c, vlpmethod='HB', well_type='oil',
               qt_stbpd=1000, gor=800, wc=0.3, oil_pvt=opvt, gsg=0.65)
    assert bhp > 200, f"BHP {bhp} should exceed THP 200"


# ============================================================================
#  Bad input Tests
# ============================================================================

def test_vlp_method_validation():
    """Invalid VLP method should raise ValueError."""
    c = Completion(tid=2.441, length=10000, tht=100, bht=200)
    try:
        fbhp(thp=500, completion=c, vlpmethod='INVALID', well_type='gas',
             qg_mmscfd=5.0, gsg=0.65)
        assert False, "Should have raised ValueError"
    except ValueError:
        pass


# ============================================================================
#  WellSegment and Multi-Segment Tests
# ============================================================================

def test_wellsegment_vertical():
    """Vertical WellSegment has TVD == MD and theta == pi/2."""
    seg = WellSegment(md=10000, id=2.441, deviation=0)
    assert abs(seg.tvd - 10000.0) < 0.01
    assert abs(seg.theta - math.pi / 2.0) < 1e-10


def test_wellsegment_horizontal():
    """Horizontal WellSegment has TVD ~0 and theta ~0."""
    seg = WellSegment(md=5000, id=2.441, deviation=90)
    assert abs(seg.tvd) < 0.01
    assert abs(seg.theta) < 1e-10


def test_wellsegment_deviated():
    """45-degree deviation gives TVD = MD * cos(45)."""
    seg = WellSegment(md=10000, id=2.441, deviation=45)
    expected_tvd = 10000 * math.cos(math.radians(45))
    assert abs(seg.tvd - expected_tvd) < 0.01
    expected_theta = math.pi / 4.0
    assert abs(seg.theta - expected_theta) < 1e-10


def test_completion_segments_mode():
    """Completion constructed with segments list."""
    segs = [WellSegment(md=5000, id=2.441), WellSegment(md=3000, id=4.0)]
    c = Completion(segments=segs, tht=100, bht=250)
    assert len(c.segments) == 2
    assert abs(c.total_md - 8000.0) < 0.01
    assert abs(c.total_tvd - 8000.0) < 0.01
    assert c.tht == 100
    assert c.bht == 250


def test_completion_legacy_builds_segments():
    """Legacy Completion constructor builds internal segments list."""
    c = Completion(tid=2.441, length=10000, tht=100, bht=200)
    assert len(c.segments) == 1
    assert abs(c.segments[0].md - 10000) < 0.01
    assert abs(c.segments[0].id - 2.441) < 0.001


def test_completion_legacy_casing_builds_two_segments():
    """Legacy Completion with casing builds two segments."""
    c = Completion(tid=2.441, length=8000, tht=100, bht=250,
                   cid=6.0, crough=0.001, mpd=10000)
    assert len(c.segments) == 2
    assert abs(c.segments[0].md - 8000) < 0.01
    assert abs(c.segments[1].md - 2000) < 0.01
    assert abs(c.segments[1].id - 6.0) < 0.01


def test_single_segment_matches_legacy_gas():
    """Single vertical WellSegment gives same BHP as legacy Completion for all 4 methods."""
    legacy = Completion(tid=2.441, length=10000, tht=100, bht=200)
    seg_comp = Completion(segments=[WellSegment(md=10000, id=2.441)], tht=100, bht=200)
    for method in ['HB', 'WG', 'GRAY', 'BB']:
        bhp_legacy = fbhp(thp=500, completion=legacy, vlpmethod=method,
                          well_type='gas', qg_mmscfd=5, cgr=10, qw_bwpd=10,
                          gsg=0.65, wsg=1.07)
        bhp_seg = fbhp(thp=500, completion=seg_comp, vlpmethod=method,
                       well_type='gas', qg_mmscfd=5, cgr=10, qw_bwpd=10,
                       gsg=0.65, wsg=1.07)
        assert abs(bhp_legacy - bhp_seg) < 0.01, \
            f"Method {method}: legacy={bhp_legacy:.4f} vs seg={bhp_seg:.4f}"


def test_single_segment_matches_legacy_oil():
    """Single vertical WellSegment gives same BHP as legacy Completion for oil wells."""
    legacy = Completion(tid=2.441, length=10000, tht=120, bht=220)
    seg_comp = Completion(segments=[WellSegment(md=10000, id=2.441)], tht=120, bht=220)
    for method in ['HB', 'WG', 'GRAY', 'BB']:
        bhp_legacy = fbhp(thp=200, completion=legacy, vlpmethod=method,
                          well_type='oil', qt_stbpd=1000, gor=500, wc=0.3,
                          gsg=0.65, wsg=1.07, api=35, pb=3000, rsb=500, sgsp=0.65)
        bhp_seg = fbhp(thp=200, completion=seg_comp, vlpmethod=method,
                       well_type='oil', qt_stbpd=1000, gor=500, wc=0.3,
                       gsg=0.65, wsg=1.07, api=35, pb=3000, rsb=500, sgsp=0.65)
        assert abs(bhp_legacy - bhp_seg) < 0.01, \
            f"Method {method}: legacy={bhp_legacy:.4f} vs seg={bhp_seg:.4f}"


def test_deviated_less_than_vertical_gas():
    """45-degree deviated well gives lower BHP than vertical (less hydrostatic)."""
    vert = Completion(segments=[WellSegment(md=10000, id=2.441, deviation=0)],
                      tht=100, bht=200)
    dev = Completion(segments=[WellSegment(md=10000, id=2.441, deviation=45)],
                     tht=100, bht=200)
    for method in ['HB', 'GRAY', 'BB']:
        bhp_vert = fbhp(thp=500, completion=vert, vlpmethod=method,
                         well_type='gas', qg_mmscfd=5, cgr=10, qw_bwpd=10,
                         gsg=0.65, wsg=1.07)
        bhp_dev = fbhp(thp=500, completion=dev, vlpmethod=method,
                        well_type='gas', qg_mmscfd=5, cgr=10, qw_bwpd=10,
                        gsg=0.65, wsg=1.07)
        assert bhp_dev < bhp_vert, \
            f"Method {method}: deviated BHP {bhp_dev:.2f} should be < vertical {bhp_vert:.2f}"


def test_horizontal_near_thp_gas():
    """Horizontal well BHP should be close to THP (minimal hydrostatic, mostly friction)."""
    horiz = Completion(segments=[WellSegment(md=5000, id=2.441, deviation=90)],
                       tht=150, bht=150)
    vert = Completion(segments=[WellSegment(md=5000, id=2.441, deviation=0)],
                      tht=150, bht=150)
    bhp_horiz = fbhp(thp=500, completion=horiz, vlpmethod='HB',
                     well_type='gas', qg_mmscfd=5, cgr=10, qw_bwpd=10,
                     gsg=0.65, wsg=1.07)
    bhp_vert = fbhp(thp=500, completion=vert, vlpmethod='HB',
                    well_type='gas', qg_mmscfd=5, cgr=10, qw_bwpd=10,
                    gsg=0.65, wsg=1.07)
    # Horizontal should have significantly lower BHP than vertical (no hydrostatic head)
    assert bhp_horiz < bhp_vert, \
        f"Horizontal BHP {bhp_horiz:.2f} should be < vertical BHP {bhp_vert:.2f}"
    assert bhp_horiz > 490, f"Horizontal BHP {bhp_horiz:.2f} should not be below THP"


def test_multi_segment_gas_all_methods():
    """Multi-segment completion runs without error for all methods."""
    segs = [
        WellSegment(md=3000, id=2.441, deviation=0),    # vertical
        WellSegment(md=4000, id=2.441, deviation=30),   # build section
        WellSegment(md=3000, id=2.441, deviation=60),   # near-horizontal
    ]
    comp = Completion(segments=segs, tht=100, bht=250)
    for method in ['HB', 'WG', 'GRAY', 'BB']:
        bhp = fbhp(thp=500, completion=comp, vlpmethod=method,
                   well_type='gas', qg_mmscfd=5, cgr=10, qw_bwpd=10,
                   gsg=0.65, wsg=1.07)
        assert bhp > 500, f"Method {method}: BHP {bhp:.2f} should exceed THP"


def test_multi_segment_oil_all_methods():
    """Multi-segment completion runs without error for oil wells."""
    segs = [
        WellSegment(md=5000, id=2.441, deviation=0),
        WellSegment(md=3000, id=4.0, deviation=0, roughness=0.001),
    ]
    comp = Completion(segments=segs, tht=120, bht=220)
    for method in ['HB', 'WG', 'GRAY', 'BB']:
        bhp = fbhp(thp=200, completion=comp, vlpmethod=method,
                   well_type='oil', qt_stbpd=1000, gor=500, wc=0.3,
                   gsg=0.65, wsg=1.07, api=35, pb=3000, rsb=500, sgsp=0.65)
        assert bhp > 200, f"Method {method}: BHP {bhp:.2f} should exceed THP"


# ============================================================================
#  geometry_at_md() Tests
# ============================================================================

def test_geometry_at_md_single_segment():
    """Single vertical segment: geometry at surface, midpoint, and bottom."""
    c = Completion(tid=2.441, length=10000, tht=100, bht=200)
    # Surface
    g0 = c.geometry_at_md(0)
    assert abs(g0['md']) < 1e-10
    assert abs(g0['tvd']) < 1e-10
    assert abs(g0['id'] - 2.441) < 1e-6
    assert g0['deviation'] == 0
    # Midpoint
    g5 = c.geometry_at_md(5000)
    assert abs(g5['md'] - 5000) < 1e-6
    assert abs(g5['tvd'] - 5000) < 1e-6  # vertical: tvd == md
    assert abs(g5['id'] - 2.441) < 1e-6
    # Bottom
    g10 = c.geometry_at_md(10000)
    assert abs(g10['md'] - 10000) < 1e-6
    assert abs(g10['tvd'] - 10000) < 1e-6


def test_geometry_at_md_multi_segment():
    """3-segment deviated well: geometry within each segment and at crossovers."""
    segs = [
        WellSegment(md=3000, id=2.441, deviation=0),    # vertical
        WellSegment(md=4000, id=2.441, deviation=45),   # build section
        WellSegment(md=3000, id=4.0, deviation=60),     # near-horizontal
    ]
    c = Completion(segments=segs, tht=100, bht=250)
    # Within first segment (md=1500)
    g1 = c.geometry_at_md(1500)
    assert abs(g1['tvd'] - 1500) < 0.1  # vertical
    assert abs(g1['id'] - 2.441) < 1e-6
    assert g1['deviation'] == 0
    # At crossover point (md=3000 — end of first, start of second)
    g2 = c.geometry_at_md(3000)
    assert abs(g2['tvd'] - 3000) < 0.1
    assert abs(g2['id'] - 2.441) < 1e-6
    assert g2['deviation'] == 0  # still in first segment
    # Within second segment (md=5000 — 2000 ft into 45-degree segment)
    g3 = c.geometry_at_md(5000)
    expected_tvd = 3000 + 2000 * math.cos(math.radians(45))
    assert abs(g3['tvd'] - expected_tvd) < 0.1
    assert g3['deviation'] == 45
    # Within third segment (md=8000 — 1000 ft into 60-degree segment)
    g4 = c.geometry_at_md(8000)
    tvd_seg1 = 3000
    tvd_seg2 = 4000 * math.cos(math.radians(45))
    tvd_seg3_partial = 1000 * math.cos(math.radians(60))
    expected_tvd = tvd_seg1 + tvd_seg2 + tvd_seg3_partial
    assert abs(g4['tvd'] - expected_tvd) < 0.1
    assert abs(g4['id'] - 4.0) < 1e-6
    assert g4['deviation'] == 60
    # At total_md (md=10000)
    g5 = c.geometry_at_md(10000)
    total_tvd = tvd_seg1 + tvd_seg2 + 3000 * math.cos(math.radians(60))
    assert abs(g5['tvd'] - total_tvd) < 0.1


def test_geometry_at_md_out_of_range():
    """md < 0 and md > total_md should raise ValueError."""
    c = Completion(tid=2.441, length=10000, tht=100, bht=200)
    try:
        c.geometry_at_md(-1)
        assert False, "Should have raised ValueError for md < 0"
    except ValueError:
        pass
    try:
        c.geometry_at_md(10001)
        assert False, "Should have raised ValueError for md > total_md"
    except ValueError:
        pass


def test_geometry_at_md_metric():
    """Metric Completion: input and output should be in metric units."""
    segs = [WellSegment(md=3000, id=62, deviation=0, metric=True),
            WellSegment(md=2000, id=100, deviation=45, metric=True)]
    c = Completion(segments=segs, tht=50, bht=90, metric=True)
    # Query at 1000 m (within first segment)
    g = c.geometry_at_md(1000)
    assert abs(g['md'] - 1000) < 0.1
    assert abs(g['tvd'] - 1000) < 0.1  # vertical
    assert abs(g['id'] - 62) < 0.1  # mm
    assert g['deviation'] == 0
    # Query at 4000 m (within second segment, 1000 m into it)
    g2 = c.geometry_at_md(4000)
    expected_tvd = 3000 + 1000 * math.cos(math.radians(45))
    assert abs(g2['tvd'] - expected_tvd) < 0.5
    assert abs(g2['id'] - 100) < 0.1
    assert g2['deviation'] == 45


# ============================================================================
#  profile() Tests
# ============================================================================

def test_profile_single_segment():
    """Single segment profile has 2 rows: top and bottom."""
    c = Completion(tid=2.441, length=10000, tht=100, bht=200)
    df = c.profile()
    assert len(df) == 2
    assert abs(df['MD'].iloc[0]) < 1e-10
    assert abs(df['TVD'].iloc[0]) < 1e-10
    assert abs(df['MD'].iloc[1] - 10000) < 1e-6
    assert abs(df['TVD'].iloc[1] - 10000) < 1e-6
    assert abs(df['ID'].iloc[0] - 2.441) < 1e-6
    assert abs(df['ID'].iloc[1] - 2.441) < 1e-6
    assert 'Roughness' in df.columns
    assert abs(df['Roughness'].iloc[0] - 0.0006) < 1e-8
    assert abs(df['Roughness'].iloc[1] - 0.0006) < 1e-8


def test_profile_multi_segment():
    """3-segment well: crossover rows have same MD/TVD but different ID/deviation."""
    segs = [
        WellSegment(md=3000, id=2.441, deviation=0),
        WellSegment(md=4000, id=2.441, deviation=45),
        WellSegment(md=3000, id=4.0, deviation=60),
    ]
    c = Completion(segments=segs, tht=100, bht=250)
    df = c.profile()
    # 3 segments: top + (bottom + crossover_top)*2 + bottom = 1 + 2 + 2 + 1 = 6 rows
    # Actually: top(0), bottom(0)/top(1) = 2 rows, bottom(1)/top(2) = 2 rows, bottom(2) = 1 row
    # Total = 1 + 2 + 2 + 1 = 6? Let's count:
    #   i=0: top of seg0, bottom of seg0, top of seg1 = 3 rows
    #   i=1: bottom of seg1, top of seg2 = 2 rows
    #   i=2: bottom of seg2 = 1 row
    # Total = 3 + 2 + 1 = 6
    assert len(df) == 6, f"Expected 6 rows, got {len(df)}"
    # Check crossover at md=3000: rows 1 and 2 should have same MD/TVD but different deviation
    assert abs(df['MD'].iloc[1] - df['MD'].iloc[2]) < 1e-6
    assert abs(df['TVD'].iloc[1] - df['TVD'].iloc[2]) < 1e-6
    assert df['Deviation'].iloc[1] == 0   # end of segment 0
    assert df['Deviation'].iloc[2] == 45  # start of segment 1
    # Check crossover at md=7000: rows 3 and 4
    assert abs(df['MD'].iloc[3] - df['MD'].iloc[4]) < 1e-6
    assert df['Deviation'].iloc[3] == 45  # end of segment 1
    assert df['Deviation'].iloc[4] == 60  # start of segment 2
    assert abs(df['ID'].iloc[4] - 4.0) < 1e-6


def test_profile_legacy_with_casing():
    """Legacy completion with casing section shows tubing->casing transition."""
    c = Completion(tid=2.441, length=8000, tht=100, bht=250,
                   cid=6.0, crough=0.001, mpd=10000)
    df = c.profile()
    # 2 segments -> top, bottom/crossover_top, bottom = 4 rows
    assert len(df) == 4, f"Expected 4 rows, got {len(df)}"
    # Crossover at 8000 ft
    assert abs(df['MD'].iloc[1] - 8000) < 1e-6
    assert abs(df['ID'].iloc[1] - 2.441) < 1e-6  # tubing
    assert abs(df['MD'].iloc[2] - 8000) < 1e-6
    assert abs(df['ID'].iloc[2] - 6.0) < 1e-6    # casing
    # Bottom at 10000 ft
    assert abs(df['MD'].iloc[3] - 10000) < 1e-6
    assert abs(df['ID'].iloc[3] - 6.0) < 1e-6
    # Roughness: tubing=0.0006, casing=0.001
    assert abs(df['Roughness'].iloc[0] - 0.0006) < 1e-8
    assert abs(df['Roughness'].iloc[2] - 0.001) < 1e-8


def test_profile_metric():
    """Metric Completion profile should have values in metric units."""
    segs = [WellSegment(md=3000, id=62, deviation=0, metric=True),
            WellSegment(md=2000, id=100, deviation=45, metric=True)]
    c = Completion(segments=segs, tht=50, bht=90, metric=True)
    df = c.profile()
    # Check first row is surface in metric
    assert abs(df['MD'].iloc[0]) < 1e-6
    assert abs(df['TVD'].iloc[0]) < 1e-6
    assert abs(df['ID'].iloc[0] - 62) < 0.1  # mm
    # Check last row total MD ~ 5000 m
    assert abs(df['MD'].iloc[-1] - 5000) < 0.5
    # Check ID of second segment ~ 100 mm
    assert abs(df['ID'].iloc[-1] - 100) < 0.1
    # Roughness should be in mm (default 0.01524 mm)
    assert abs(df['Roughness'].iloc[0] - 0.01524) < 0.001


if __name__ == '__main__':
    import traceback
    tests = [(k, v) for k, v in list(globals().items()) if k.startswith('test_') and callable(v)]
    passed, failed = 0, 0
    for name, func in sorted(tests):
        try:
            func()
            passed += 1
            print(f"  PASS: {name}")
        except Exception as e:
            failed += 1
            print(f"  FAIL: {name}")
            print(f"        {e}")
            traceback.print_exc()
    print(f"\n{passed} passed, {failed} failed out of {passed + failed}")
