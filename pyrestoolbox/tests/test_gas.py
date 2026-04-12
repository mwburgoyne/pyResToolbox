#!/usr/bin/env python3
"""
Validation tests for gas module.
Run with: PYTHONPATH=/home/mark/projects python3 -m pytest tests/ -v
Or standalone: PYTHONPATH=/home/mark/projects python3 tests/test_gas.py
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import pyrestoolbox.gas as gas

RTOL = 0.02  # 2% relative tolerance for correlation-based results
ATOL_Z = 0.005  # Absolute tolerance for Z-factor (near 1.0)

# =============================================================================
# Baseline values captured from current implementation
# These serve as regression anchors - if values change, tests alert us
# =============================================================================

def test_gas_z_dak_single():
    """DAK Z-factor for typical gas at moderate conditions"""
    z = gas.gas_z(p=2000, sg=0.75, degf=200, zmethod='DAK', cmethod='PMC')
    assert isinstance(z, float), f"Expected float, got {type(z)}"
    assert 0.5 < z < 1.2, f"Z={z} outside physical bounds"
    # Capture baseline
    baseline = gas.gas_z(p=2000, sg=0.75, degf=200, zmethod='DAK', cmethod='PMC')
    assert abs(z - baseline) < 1e-10, "DAK Z-factor not reproducible"

def test_gas_z_dak_array():
    """DAK Z-factor for array of pressures"""
    pressures = np.array([500, 1000, 2000, 4000, 8000])
    z = gas.gas_z(p=pressures, sg=0.75, degf=200, zmethod='DAK', cmethod='PMC')
    assert isinstance(z, np.ndarray), f"Expected ndarray, got {type(z)}"
    assert len(z) == len(pressures)
    # Z should decrease then increase with pressure (typical behavior)
    assert all(0.3 < zi < 1.5 for zi in z), f"Z values outside physical bounds: {z}"

def test_gas_z_hy():
    """Hall-Yarborough Z-factor"""
    z = gas.gas_z(p=2000, sg=0.75, degf=200, zmethod='HY', cmethod='PMC')
    z_dak = gas.gas_z(p=2000, sg=0.75, degf=200, zmethod='DAK', cmethod='PMC')
    assert isinstance(z, float)
    assert 0.5 < z < 1.2
    # HY and DAK should agree within ~3% for typical conditions
    assert abs(z - z_dak) / z_dak < 0.03, f"HY={z} vs DAK={z_dak}, diff > 3%"

def test_gas_z_wyw():
    """Wang-Ye-Wu Z-factor"""
    z = gas.gas_z(p=2000, sg=0.75, degf=200, zmethod='WYW', cmethod='PMC')
    z_dak = gas.gas_z(p=2000, sg=0.75, degf=200, zmethod='DAK', cmethod='PMC')
    assert isinstance(z, float)
    assert 0.5 < z < 1.2
    assert abs(z - z_dak) / z_dak < 0.05, f"WYW={z} vs DAK={z_dak}, diff > 5%"

def test_gas_z_bns():
    """BNS (Peng-Robinson EOS) Z-factor"""
    z = gas.gas_z(p=2000, sg=0.75, degf=200, zmethod='BNS', cmethod='BNS')
    assert isinstance(z, float)
    assert 0.5 < z < 1.2, f"BNS Z={z} outside physical bounds"

def test_gas_z_with_h2():
    """Z-factor with hydrogen forces BNS method"""
    z = gas.gas_z(p=2000, sg=0.75, degf=200, h2=0.1)
    assert isinstance(z, float)
    assert 0.5 < z < 1.5

def test_gas_z_with_impurities():
    """Z-factor with CO2, H2S, N2 impurities"""
    z = gas.gas_z(p=2000, sg=0.8, degf=200, co2=0.05, h2s=0.03, n2=0.02)
    assert isinstance(z, float)
    assert 0.3 < z < 1.5

def test_gas_z_low_pressure():
    """Z at very low pressure should approach 1.0"""
    z = gas.gas_z(p=14.7, sg=0.75, degf=200)
    assert abs(z - 1.0) < 0.05, f"Z at 14.7 psia = {z}, should be ~1.0"

def test_gas_tc_pc_pmc():
    """PMC critical property correlation"""
    tc, pc = gas.gas_tc_pc(sg=0.75, cmethod='PMC')
    assert 300 < tc < 500, f"Tc={tc} outside expected range"
    assert 500 < pc < 800, f"Pc={pc} outside expected range"

def test_gas_tc_pc_sut():
    """Sutton critical property correlation"""
    tc, pc = gas.gas_tc_pc(sg=0.75, cmethod='SUT')
    assert 300 < tc < 500, f"Tc={tc} outside expected range"
    assert 500 < pc < 800, f"Pc={pc} outside expected range"

def test_gas_tc_pc_user_override():
    """User-specified critical properties should be returned as-is"""
    tc, pc = gas.gas_tc_pc(sg=0.75, tc=400, pc=650)
    assert tc == 400
    assert pc == 650

def test_gas_viscosity():
    """Gas viscosity at typical conditions"""
    ug = gas.gas_ug(p=2000, sg=0.75, degf=200)
    assert isinstance(ug, float)
    assert 0.005 < ug < 0.1, f"Gas viscosity = {ug} cP, outside expected range"

def test_gas_viscosity_array():
    """Gas viscosity for array of pressures"""
    pressures = np.array([500, 1000, 2000, 4000])
    ug = gas.gas_ug(p=pressures, sg=0.75, degf=200)
    assert isinstance(ug, np.ndarray)
    assert len(ug) == len(pressures)
    # Viscosity should increase with pressure
    for i in range(len(ug) - 1):
        assert ug[i + 1] > ug[i] * 0.9, "Viscosity should generally increase with pressure"

def test_gas_viscosity_bns():
    """LBC viscosity model via BNS method"""
    ug = gas.gas_ug(p=2000, sg=0.75, degf=200, zmethod='BNS', cmethod='BNS')
    assert isinstance(ug, float)
    assert 0.005 < ug < 0.1, f"BNS gas viscosity = {ug} cP, outside expected range"

def test_gas_compressibility():
    """Gas compressibility at typical conditions"""
    cg = gas.gas_cg(p=2000, sg=0.75, degf=200)
    assert isinstance(cg, float)
    assert 0 < cg < 0.01, f"Cg = {cg}, outside expected range"

def test_gas_bg():
    """Gas FVF"""
    bg = gas.gas_bg(p=2000, sg=0.75, degf=200)
    assert isinstance(bg, float)
    assert 0 < bg < 0.1, f"Bg = {bg} rcf/scf, outside expected range"

def test_gas_bg_decreases_with_pressure():
    """Bg should decrease with increasing pressure"""
    pressures = np.array([500, 1000, 2000, 4000])
    bg = gas.gas_bg(p=pressures, sg=0.75, degf=200)
    for i in range(len(bg) - 1):
        assert bg[i + 1] < bg[i], f"Bg not decreasing: {bg}"

def test_gas_density():
    """Gas density"""
    rho = gas.gas_den(p=2000, sg=0.75, degf=200)
    assert isinstance(rho, float)
    assert 1 < rho < 30, f"Gas density = {rho} lb/cuft, outside expected range"

def test_gas_density_increases_with_pressure():
    """Gas density should increase with pressure"""
    pressures = np.array([500, 1000, 2000, 4000])
    rho = gas.gas_den(p=pressures, sg=0.75, degf=200)
    for i in range(len(rho) - 1):
        assert rho[i + 1] > rho[i], f"Density not increasing: {rho}"

def test_gas_ponz2p():
    """P/Z to P conversion roundtrip"""
    p_original = 3000
    z = gas.gas_z(p=p_original, sg=0.75, degf=200)
    ponz = p_original / z
    p_recovered = gas.gas_ponz2p(poverz=ponz, sg=0.75, degf=200)
    assert abs(p_recovered - p_original) / p_original < 0.001, \
        f"P/Z roundtrip: original={p_original}, recovered={p_recovered}"

def test_gas_grad2sg():
    """Gradient to SG inversion"""
    sg_input = 0.75
    p = 3000.0
    degf = 200.0
    rho = gas.gas_den(p=p, sg=sg_input, degf=degf)
    grad = float(rho) / 144  # psi/ft
    sg_recovered = gas.gas_grad2sg(grad=grad, p=p, degf=degf)
    assert abs(sg_recovered - sg_input) / sg_input < 0.01, \
        f"Grad2SG: input={sg_input}, recovered={sg_recovered}"

def test_gas_dmp():
    """Delta pseudo pressure integration"""
    dmp = gas.gas_dmp(p1=1000, p2=3000, degf=200, sg=0.75)
    assert dmp > 0, "Pseudopressure integral should be positive when p2 > p1"
    assert isinstance(dmp, float)

def test_gas_dmp_zero():
    """Delta pseudo pressure should be zero when p1 = p2"""
    dmp = gas.gas_dmp(p1=2000, p2=2000, degf=200, sg=0.75)
    assert dmp == 0

def test_gas_rate_radial():
    """Radial gas rate calculation"""
    qg = gas.gas_rate_radial(k=10, h=50, pr=3000, pwf=1500, r_w=0.35, r_ext=1000, degf=200, sg=0.75)
    assert isinstance(qg, (float, np.floating)), f"Expected float, got {type(qg)}"
    assert qg > 0, "Flow rate should be positive when Pr > Pwf"

def test_gas_rate_radial_reverse():
    """Radial gas rate with reversed pressures should give negative rate"""
    qg = gas.gas_rate_radial(k=10, h=50, pr=1500, pwf=3000, r_w=0.35, r_ext=1000, degf=200, sg=0.75)
    assert qg < 0, "Flow rate should be negative when Pr < Pwf"

def test_gas_rate_linear():
    """Linear gas rate calculation"""
    qg = gas.gas_rate_linear(k=10, pr=3000, pwf=1500, area=10000, length=5000, degf=200, sg=0.75)
    assert isinstance(qg, (float, np.floating))
    assert qg > 0

def test_gas_water_content():
    """Saturated water content in gas"""
    wc = gas.gas_water_content(p=2000, degf=200)
    assert isinstance(wc, float)
    assert wc > 0, "Water content should be positive"

def test_gas_fws_sg():
    """Full wellstream SG from gas condensate"""
    sg = gas.gas_fws_sg(sg_g=0.75, cgr=50, api_st=50)
    assert isinstance(sg, float)
    assert sg > 0.75, "FWS SG should be higher than separator gas SG"

# =============================================================================
# Hard-coded regression baselines - frozen reference values
# If any of these change, it means a code modification has altered computational results.
# Do NOT update these values without explicit review and approval.
# =============================================================================

# Frozen baselines captured 2026-02-27
_FROZEN_BASELINES = {
    'z_dak_2000_075_200': 0.8682651934424434,
    'z_hy_2000_075_200': 0.8677386022105378,
    'z_wyw_2000_075_200': 0.8570731189820472,
    'z_bns_2000_075_200': 0.8457226741039529,
    'ug_dak_2000_075_200': 0.017374276181391906,
    'bg_dak_2000_075_200': 0.008098799120908043,
    'cg_dak_2000_075_200': 0.0005200298285282161,
    'den_dak_2000_075_200': 7.069636415667095,
    'tc_pmc_075': (385.66677319794064, 653.255423805908),
    'tc_sut_075': (389.7, 656.525),
    'dmp_1000_2000_hy_sut': 213690308.9907268,
    'dmp_0_4000_dak': 919822133.7011306,
    'qg_radial': 2078.9101970773477,
}

def test_regression_z_dak():
    z = gas.gas_z(p=2000, sg=0.75, degf=200, zmethod='DAK', cmethod='PMC')
    expected = _FROZEN_BASELINES['z_dak_2000_075_200']
    assert abs(z - expected) / expected < 1e-8, f"DAK Z-factor changed: {z} vs frozen {expected}"

def test_regression_z_hy():
    z = gas.gas_z(p=2000, sg=0.75, degf=200, zmethod='HY', cmethod='PMC')
    expected = _FROZEN_BASELINES['z_hy_2000_075_200']
    assert abs(z - expected) / expected < 1e-8, f"HY Z-factor changed: {z} vs frozen {expected}"

def test_regression_z_wyw():
    z = gas.gas_z(p=2000, sg=0.75, degf=200, zmethod='WYW', cmethod='PMC')
    expected = _FROZEN_BASELINES['z_wyw_2000_075_200']
    assert abs(z - expected) / expected < 1e-8, f"WYW Z-factor changed: {z} vs frozen {expected}"

def test_regression_ug():
    ug = gas.gas_ug(p=2000, sg=0.75, degf=200)
    expected = _FROZEN_BASELINES['ug_dak_2000_075_200']
    assert abs(ug - expected) / expected < 1e-6, f"Gas viscosity changed: {ug} vs frozen {expected}"

def test_regression_bg():
    bg = gas.gas_bg(p=2000, sg=0.75, degf=200)
    expected = _FROZEN_BASELINES['bg_dak_2000_075_200']
    assert abs(bg - expected) / expected < 1e-6, f"Gas Bg changed: {bg} vs frozen {expected}"

def test_regression_z_bns():
    z = gas.gas_z(p=2000, sg=0.75, degf=200, zmethod='BNS', cmethod='BNS')
    expected = _FROZEN_BASELINES['z_bns_2000_075_200']
    assert abs(z - expected) / expected < 1e-8, f"BNS Z-factor changed: {z} vs frozen {expected}"

def test_regression_dmp():
    """Delta pseudopressure must match frozen baseline"""
    dmp = gas.gas_dmp(p1=1000, p2=2000, degf=185, sg=0.78, zmethod='HY', cmethod='SUT', n2=0.05, co2=0.1, h2s=0.02)
    expected = _FROZEN_BASELINES['dmp_1000_2000_hy_sut']
    assert abs(dmp - expected) / expected < 1e-6, f"gas_dmp changed: {dmp} vs frozen {expected}"

def test_regression_dmp_wide_range():
    """Delta pseudopressure 0-4000 psia must match frozen baseline"""
    dmp = gas.gas_dmp(p1=0, p2=4000, sg=0.75, degf=200)
    expected = _FROZEN_BASELINES['dmp_0_4000_dak']
    assert abs(dmp - expected) / expected < 1e-6, f"gas_dmp wide range changed: {dmp} vs frozen {expected}"

def test_regression_gas_rate():
    """Gas rate radial must match frozen baseline"""
    qg = gas.gas_rate_radial(k=5, h=50, pr=2000, pwf=750, r_w=0.3, r_ext=1500, degf=180, sg=0.75, D=0.01, S=5)
    expected = _FROZEN_BASELINES['qg_radial']
    assert abs(qg - expected) / expected < 1e-6, f"gas_rate_radial changed: {qg} vs frozen {expected}"

def test_regression_tc_pc():
    """Critical properties must match frozen baselines"""
    tc, pc = gas.gas_tc_pc(sg=0.75, cmethod='PMC')
    tc_exp, pc_exp = _FROZEN_BASELINES['tc_pmc_075']
    assert abs(tc - tc_exp) / tc_exp < 1e-8, f"PMC Tc changed: {tc} vs frozen {tc_exp}"
    assert abs(pc - pc_exp) / pc_exp < 1e-8, f"PMC Pc changed: {pc} vs frozen {pc_exp}"

def test_regression_cg():
    """Gas compressibility must match frozen baseline"""
    cg = gas.gas_cg(p=2000, sg=0.75, degf=200)
    expected = _FROZEN_BASELINES['cg_dak_2000_075_200']
    assert abs(cg - expected) / expected < 1e-6, f"gas_cg changed: {cg} vs frozen {expected}"

def test_regression_den():
    """Gas density must match frozen baseline"""
    den = gas.gas_den(p=2000, sg=0.75, degf=200)
    expected = _FROZEN_BASELINES['den_dak_2000_075_200']
    assert abs(den - expected) / expected < 1e-6, f"gas_den changed: {den} vs frozen {expected}"

# =============================================================================
# List input tests (convert_to_numpy fix)
# =============================================================================

def test_gas_z_list_input():
    """gas_z should accept Python lists, not just numpy arrays"""
    z_list = gas.gas_z(p=[500, 1000, 2000], sg=0.75, degf=200)
    z_array = gas.gas_z(p=np.array([500, 1000, 2000]), sg=0.75, degf=200)
    assert isinstance(z_list, np.ndarray), f"Expected ndarray from list input, got {type(z_list)}"
    np.testing.assert_allclose(z_list, z_array, rtol=1e-10)

def test_gas_bg_list_input():
    """gas_bg should accept Python lists"""
    bg = gas.gas_bg(p=[1000, 2000, 3000], sg=0.75, degf=200)
    assert isinstance(bg, np.ndarray)
    assert len(bg) == 3

# =============================================================================
# BNS implementation validation (against reference bns.py from SPE-229932-MS)
# Reference values computed from github.com/mwburgoyne/5_Component_PengRobinson_Z-Factor
# =============================================================================

def test_bns_z_pure_gas():
    """BNS Z-factor for pure hydrocarbon gas should match reference implementation"""
    z = gas.gas_z(p=2000, sg=0.75, degf=200, zmethod='BNS', cmethod='BNS')
    assert abs(z - 0.8457226741) < 1e-6, f"BNS Z={z}, expected 0.8457 (reference)"

def test_bns_z_with_co2():
    """BNS Z-factor with CO2 impurity should match reference"""
    z = gas.gas_z(p=2000, sg=0.8, degf=200, zmethod='BNS', cmethod='BNS', co2=0.05)
    assert abs(z - 0.8338808668) < 1e-6, f"BNS Z={z}, expected 0.8339 (reference)"

def test_bns_z_with_mixed_impurities():
    """BNS Z-factor with CO2+H2S+N2 should match reference"""
    z = gas.gas_z(p=2000, sg=0.8, degf=200, zmethod='BNS', cmethod='BNS', co2=0.05, h2s=0.03, n2=0.02)
    assert abs(z - 0.8409252043) < 1e-6, f"BNS Z={z}, expected 0.8409 (reference)"

def test_bns_z_with_h2():
    """BNS Z-factor with hydrogen should match reference"""
    z = gas.gas_z(p=2000, sg=0.75, degf=200, zmethod='BNS', cmethod='BNS', h2=0.1)
    assert abs(z - 0.8582347790) < 1e-6, f"BNS Z={z}, expected 0.8582 (reference)"

def test_bns_z_high_pressure():
    """BNS Z-factor at high pressure should match reference"""
    z = gas.gas_z(p=8000, sg=0.75, degf=200, zmethod='BNS', cmethod='BNS')
    assert abs(z - 1.2547712914) < 1e-6, f"BNS Z={z}, expected 1.2548 (reference)"

def test_bns_z_array_pressures():
    """BNS Z-factor should work correctly with array of pressures"""
    pressures = np.array([500, 2000, 8000])
    z = gas.gas_z(p=pressures, sg=0.75, degf=200, zmethod='BNS', cmethod='BNS')
    assert isinstance(z, np.ndarray)
    assert len(z) == 3
    # Verify each value matches scalar call
    for i, p in enumerate(pressures):
        z_scalar = gas.gas_z(p=p, sg=0.75, degf=200, zmethod='BNS', cmethod='BNS')
        assert abs(z[i] - z_scalar) < 1e-10, f"Array vs scalar mismatch at p={p}"

def test_bns_viscosity_pure_gas():
    """BNS LBC viscosity for pure gas should match reference"""
    ug = gas.gas_ug(p=2000, sg=0.75, degf=200, zmethod='BNS', cmethod='BNS')
    assert abs(ug - 0.0172603874) < 1e-6, f"BNS ug={ug}, expected 0.01726 (reference)"

def test_bns_viscosity_with_h2():
    """BNS LBC viscosity with hydrogen should match reference"""
    ug = gas.gas_ug(p=2000, sg=0.75, degf=200, zmethod='BNS', cmethod='BNS', h2=0.1)
    assert abs(ug - 0.0168458688) < 1e-6, f"BNS ug={ug}, expected 0.01685 (reference)"

def test_bns_viscosity_array_pressures():
    """BNS viscosity should return per-pressure results, not use wrong Z"""
    pressures = np.array([500, 2000, 8000])
    ug = gas.gas_ug(p=pressures, sg=0.75, degf=200, zmethod='BNS', cmethod='BNS')
    assert isinstance(ug, np.ndarray)
    assert len(ug) == 3
    # Verify each matches scalar call (catches bug where full Z array was passed to LBC)
    for i, p in enumerate(pressures):
        ug_scalar = gas.gas_ug(p=p, sg=0.75, degf=200, zmethod='BNS', cmethod='BNS')
        assert abs(ug[i] - ug_scalar) < 1e-10, f"Array vs scalar viscosity mismatch at p={p}"

def test_bns_z_subcritical_fugacity_selection():
    """BNS Z-factor must use fugacity-based root selection in the 3-root regime.
    Pure H2S at 80 degF is sub-critical (Tc_H2S = 672 R = 212 degF).
    The cubic has 3 real roots between ~330-605 psia; the thermodynamically
    stable root (lowest fugacity) switches from vapor to liquid at ~340 psia.
    Without fugacity selection the solver returns the metastable vapor root,
    producing a discontinuous Z-vs-P curve."""
    # Reference values from bns.py (SPE-229932-MS reference implementation)
    # Pressures chosen to span: vapor-only, fugacity crossover, liquid-only
    ref = {
        100:  0.9496046861917848,   # Vapor (single root)
        300:  0.83365711413685,     # Vapor (3 roots, vapor has min fugacity)
        350:  0.042605680488353545, # Liquid (3 roots, liquid has min fugacity)
        400:  0.04862660238021889,  # Liquid (3 roots, liquid has min fugacity)
        600:  0.07256104178906098,  # Liquid (3 roots merging)
        1000: 0.11977671955327113,  # Liquid (single root)
    }
    pressures = np.array(sorted(ref.keys()))
    z_arr = gas.gas_z(p=pressures, sg=0.8, degf=80, h2s=1.0, zmethod='BNS', cmethod='BNS')
    for i, p in enumerate(pressures):
        assert abs(z_arr[i] - ref[p]) < 1e-8, (
            f"BNS Z at p={p} psia: got {z_arr[i]:.10f}, expected {ref[p]:.10f}"
        )
    # Verify monotonic in the liquid region (no discontinuous jump)
    z_liquid = z_arr[pressures >= 350]
    assert np.all(np.diff(z_liquid) > 0), "Z should increase monotonically in liquid region"

def test_gas_water_content_salinity():
    """Water content should decrease with salinity"""
    wc_fresh = gas.gas_water_content(p=2000, degf=200, salinity=0)
    wc_saline = gas.gas_water_content(p=2000, degf=200, salinity=10)
    assert wc_saline < wc_fresh, f"Saline ({wc_saline}) should be less than fresh ({wc_fresh})"
    # Salinity correction factor at 10 wt%: 1 - 0.00492*10 - 0.00017672*100 = 0.93308
    expected_ratio = 1 - 0.00492 * 10 - 0.00017672 * 100
    actual_ratio = wc_saline / wc_fresh
    assert abs(actual_ratio - expected_ratio) < 0.001, f"Salinity correction ratio {actual_ratio} != {expected_ratio}"

# =============================================================================
# gas_rate_radial/linear with gas_pvt tests
# =============================================================================

def test_gas_rate_radial_with_pvt():
    """gas_rate_radial with gas_pvt should match manual params"""
    gpvt = gas.GasPVT(sg=0.75, co2=0.05)
    # Using gas_pvt
    q_pvt = gas.gas_rate_radial(k=10, h=50, pr=3000, pwf=1500, r_w=0.35, r_ext=1000,
                                degf=200, gas_pvt=gpvt)
    # Using manual params
    q_manual = gas.gas_rate_radial(k=10, h=50, pr=3000, pwf=1500, r_w=0.35, r_ext=1000,
                                   degf=200, sg=0.75, co2=0.05)
    assert abs(float(q_pvt) - float(q_manual)) / float(q_manual) < 1e-10, \
        f"gas_pvt result {q_pvt} != manual {q_manual}"

def test_gas_rate_linear_with_pvt():
    """gas_rate_linear with gas_pvt should match manual params"""
    gpvt = gas.GasPVT(sg=0.8)
    q_pvt = gas.gas_rate_linear(k=0.1, area=50, length=200, pr=2000, pwf=250,
                                degf=180, gas_pvt=gpvt)
    q_manual = gas.gas_rate_linear(k=0.1, area=50, length=200, pr=2000, pwf=250,
                                   degf=180, sg=0.8)
    assert abs(float(q_pvt) - float(q_manual)) / float(q_manual) < 1e-10, \
        f"gas_pvt result {q_pvt} != manual {q_manual}"


# =============================================================================
# Gas hydrate prediction tests
# =============================================================================

# Frozen regression baselines captured 2026-03-05
_HYDRATE_BASELINES = {
    'motiee_hft_1000_065': 96.20360674625366,     # Motiee HFT at 1000 psia, sg=0.65
    'motiee_hfp_60_065': 149.60397973632814,       # Motiee HFP at 60 degF, sg=0.65
    'towler_hft_1000_065': 62.918902535978695,      # Towler HFT at 1000 psia, sg=0.65
    'meoh_depression_25wt': 17.958375,              # Østergaard MEOH 25wt% depression (degF)
    'meoh_inhibited_hft': 91.02848313045367,        # Inhibited HFT at 2000 psia, sg=0.7, MEOH 25wt%
}

def test_hydrate_frozen_baselines():
    """Frozen regression baselines for hydrate calculations"""
    r1 = gas.gas_hydrate(p=1000, degf=60, sg=0.65, hydmethod='MOTIEE')
    assert abs(r1.hft - _HYDRATE_BASELINES['motiee_hft_1000_065']) < 1e-10, \
        f"Motiee HFT: {r1.hft} != {_HYDRATE_BASELINES['motiee_hft_1000_065']}"
    assert abs(r1.hfp - _HYDRATE_BASELINES['motiee_hfp_60_065']) < 1e-6, \
        f"Motiee HFP: {r1.hfp} != {_HYDRATE_BASELINES['motiee_hfp_60_065']}"

    r2 = gas.gas_hydrate(p=1000, degf=60, sg=0.65)  # default is now TOWLER
    assert abs(r2.hft - _HYDRATE_BASELINES['towler_hft_1000_065']) < 1e-10, \
        f"Towler HFT: {r2.hft} != {_HYDRATE_BASELINES['towler_hft_1000_065']}"

    r3 = gas.gas_hydrate(p=2000, degf=80, sg=0.7, hydmethod='MOTIEE', inhibitor_type='MEOH', inhibitor_wt_pct=25)
    assert abs(r3.inhibitor_depression - _HYDRATE_BASELINES['meoh_depression_25wt']) < 1e-10, \
        f"MEOH depression: {r3.inhibitor_depression} != {_HYDRATE_BASELINES['meoh_depression_25wt']}"
    assert abs(r3.inhibited_hft - _HYDRATE_BASELINES['meoh_inhibited_hft']) < 1e-10, \
        f"MEOH inhibited HFT: {r3.inhibited_hft} != {_HYDRATE_BASELINES['meoh_inhibited_hft']}"

def test_hydrate_hft_increases_with_pressure():
    """HFT should increase monotonically with pressure"""
    pressures = [200, 500, 1000, 2000, 4000, 8000]
    hfts = [gas.gas_hydrate(p=p, degf=60, sg=0.7).hft for p in pressures]
    for i in range(1, len(hfts)):
        assert hfts[i] > hfts[i-1], \
            f"HFT not monotonic: {hfts[i]:.2f} <= {hfts[i-1]:.2f} at P={pressures[i]} vs {pressures[i-1]}"

def test_hydrate_hfp_round_trip():
    """HFT at P -> HFP at that HFT should recover P within 1 psia"""
    for p_test in [300, 500, 1000, 2000, 5000]:
        r1 = gas.gas_hydrate(p=p_test, degf=60, sg=0.7)
        r2 = gas.gas_hydrate(p=p_test, degf=r1.hft, sg=0.7)
        assert abs(r2.hfp - p_test) < 1.0, \
            f"HFP round-trip: P={p_test}, HFT={r1.hft:.2f}, recovered HFP={r2.hfp:.2f}"

def test_hydrate_subcooling_and_window():
    """Subcooling and in_hydrate_window should be consistent"""
    # Operating below HFT -> in window
    r_cold = gas.gas_hydrate(p=2000, degf=60, sg=0.7)
    assert r_cold.in_hydrate_window, "Should be in hydrate window at 60F, 2000 psia"
    assert r_cold.subcooling > 0, "Subcooling should be positive in hydrate window"

    # Operating above HFT -> not in window
    r_hot = gas.gas_hydrate(p=200, degf=150, sg=0.6)
    assert not r_hot.in_hydrate_window, "Should NOT be in hydrate window at 150F, 200 psia"
    assert r_hot.subcooling < 0, "Subcooling should be negative outside hydrate window"

def test_hydrate_inhibitor_depression():
    """Inhibitor should reduce HFT; depression should increase with concentration"""
    r_base = gas.gas_hydrate(p=1000, degf=60, sg=0.65)
    r_inh = gas.gas_hydrate(p=1000, degf=60, sg=0.65, inhibitor_type='MEG', inhibitor_wt_pct=30)
    assert r_inh.inhibited_hft < r_base.hft, "Inhibited HFT should be lower"
    assert r_inh.inhibitor_depression > 0, "Depression should be positive"

    # Higher concentration -> more depression
    r_low = gas.gas_hydrate(p=1000, degf=60, sg=0.65, inhibitor_type='MEG', inhibitor_wt_pct=10)
    r_high = gas.gas_hydrate(p=1000, degf=60, sg=0.65, inhibitor_type='MEG', inhibitor_wt_pct=40)
    assert r_high.inhibitor_depression > r_low.inhibitor_depression, \
        "Higher concentration should give more depression"

def test_hydrate_required_concentration_round_trip():
    """Required concentration should round-trip through depression calculation"""
    r = gas.gas_hydrate(p=1000, degf=60, sg=0.65, inhibitor_type='MEG', inhibitor_wt_pct=20)
    # Compute what 20wt% gives as depression, then check required_concentration for that depression
    from pyrestoolbox.gas.gas import _ostergaard_depression, _required_concentration
    from pyrestoolbox.classes import inhibitor as inh_enum
    dep_c = _ostergaard_depression(20, inh_enum.MEG)
    conc = _required_concentration(dep_c, inh_enum.MEG)
    assert abs(conc - 20.0) < 1e-4, f"Round-trip concentration: {conc} != 20.0"

def test_hydrate_no_inhibitor():
    """Without inhibitor, inhibited_hft should be NaN and depression should be 0"""
    r = gas.gas_hydrate(p=1000, degf=60, sg=0.65)
    assert np.isnan(r.inhibited_hft), "inhibited_hft should be NaN without inhibitor"
    assert r.inhibitor_depression == 0, "depression should be 0 without inhibitor"
    assert r.required_inhibitor_wt_pct == 0, "required wt% should be 0 without inhibitor"

def test_hydrate_metric():
    """Metric results should match oilfield results after conversion"""
    r_f = gas.gas_hydrate(p=1000, degf=60, sg=0.65)
    r_m = gas.gas_hydrate(p=1000 * 0.0689475729, degf=(60 - 32) * 5 / 9, sg=0.65, metric=True)
    # HFT: convert field degF to degC
    hft_degc = (r_f.hft - 32) * 5 / 9
    assert abs(r_m.hft - hft_degc) < 0.01, f"Metric HFT {r_m.hft:.4f} != {hft_degc:.4f}"
    # HFP: convert field psia to barsa
    hfp_bar = r_f.hfp * 0.0689475729
    assert abs(r_m.hfp - hfp_bar) < 0.01, f"Metric HFP {r_m.hfp:.4f} != {hfp_bar:.4f}"

def test_hydrate_metric_with_inhibitor():
    """Metric inhibitor results should be consistent with oilfield"""
    r_f = gas.gas_hydrate(p=2000, degf=80, sg=0.7, inhibitor_type='MEG', inhibitor_wt_pct=20)
    r_m = gas.gas_hydrate(p=2000 * 0.0689475729, degf=(80 - 32) * 5 / 9, sg=0.7,
                          inhibitor_type='MEG', inhibitor_wt_pct=20, metric=True)
    # Depression: field degF delta -> degC delta
    dep_degc = r_f.inhibitor_depression * 5 / 9
    assert abs(r_m.inhibitor_depression - dep_degc) < 0.01, \
        f"Metric depression {r_m.inhibitor_depression:.4f} != {dep_degc:.4f}"
    # Required wt% should be identical (unit-independent)
    assert abs(r_m.required_inhibitor_wt_pct - r_f.required_inhibitor_wt_pct) < 0.01, \
        f"Metric required wt% {r_m.required_inhibitor_wt_pct:.2f} != {r_f.required_inhibitor_wt_pct:.2f}"

def test_hydrate_towler_vs_motiee():
    """Both methods should give physically reasonable HFT; may differ numerically"""
    r_m = gas.gas_hydrate(p=2000, degf=80, sg=0.7, hydmethod='MOTIEE')
    r_t = gas.gas_hydrate(p=2000, degf=80, sg=0.7, hydmethod='TOWLER')
    # Both should be reasonable (between 30 and 200 degF at 2000 psia)
    assert 30 < r_m.hft < 200, f"Motiee HFT out of range: {r_m.hft}"
    assert 30 < r_t.hft < 200, f"Towler HFT out of range: {r_t.hft}"

def test_hydrate_all_inhibitors():
    """All five inhibitor types should work and give positive depression"""
    for inh_name in ['MEOH', 'MEG', 'DEG', 'TEG', 'ETOH']:
        r = gas.gas_hydrate(p=1000, degf=60, sg=0.65, inhibitor_type=inh_name, inhibitor_wt_pct=15)
        assert r.inhibitor_depression > 0, f"{inh_name}: depression should be positive"
        assert r.inhibited_hft < r.hft, f"{inh_name}: inhibited HFT should be < uninhibited"

def test_hydrate_bad_inputs():
    """Invalid inputs should raise ValueError"""
    bad_cases = [
        dict(p=-100, degf=60, sg=0.65),          # negative pressure
        dict(p=1000, degf=60, sg=-0.5),           # negative SG
        dict(p=1000, degf=60, sg=0.65, inhibitor_wt_pct=-5),   # negative wt%
        dict(p=1000, degf=60, sg=0.65, inhibitor_wt_pct=100),  # wt% = 100
        dict(p=1000, degf=60, sg=0.65, hydmethod='INVALID'),   # bad method
    ]
    for kwargs in bad_cases:
        try:
            gas.gas_hydrate(**kwargs)
            raise AssertionError(f"Expected ValueError for {kwargs}")
        except ValueError:
            pass

def test_hydrate_zero_inhibitor_wt_pct():
    """Inhibitor type with 0 wt% should give no depression"""
    r = gas.gas_hydrate(p=1000, degf=60, sg=0.65, inhibitor_type='MEG', inhibitor_wt_pct=0)
    assert r.inhibitor_depression == 0, "Depression should be 0 at 0 wt%"
    assert r.inhibited_hft == r.hft, "Inhibited HFT should equal uninhibited at 0 wt%"


# =============================================================================
# Water balance, capping, injection rate, reservoir P/T tests
# =============================================================================

# Frozen regression baselines for water balance and injection (captured 2026-03-06)
_HYDRATE_NEW_BASELINES = {
    'danesh_wc_op_1000_60_065': 0.051307948560405346,    # Danesh vaporized water at operating
    'sw_wc_op_2000_80_070_co2': 0.07027556282207788,     # SoreideWhitson vaporized at operating
    'wc_res_3000_200': 0.821279572319976,                 # Vaporized at reservoir P=3000,T=200
    'condensed_3000_200_to_1000_60': 0.7699716237595706,  # Condensed between res→op
    'meg_mass_rate_res_to_op': 576.8222533381149,         # MEG injection lb/MMscf
    'meg_vol_rate_res_to_op': 62.26896381181326,          # MEG injection gal/MMscf
    'meoh_mass_rate_res_3000_200': 91.91442238644271,     # MEOH capped injection lb/MMscf
}

def test_hydrate_water_balance_no_reservoir():
    """Without reservoir P,T: vaporized_res = vaporized_op, condensed = 0"""
    r = gas.gas_hydrate(p=1000, degf=60, sg=0.65)
    wc_direct = gas.gas_water_content(p=1000, degf=60, salinity=0)
    assert abs(r.water_vaporized_op - wc_direct) < 1e-10, \
        f"vaporized_op {r.water_vaporized_op} != gas_water_content {wc_direct}"
    assert r.water_vaporized_res == r.water_vaporized_op, \
        "Without reservoir P,T, vaporized_res should equal vaporized_op"
    assert r.water_condensed == 0.0, "No condensation without reservoir P,T"
    assert r.free_water == 0.0, "No free water without additional_water"
    assert r.total_liquid_water == 0.0, "No liquid water without reservoir or free water"

def test_hydrate_water_balance_with_reservoir():
    """With reservoir P,T: condensed = vaporized_res - vaporized_op"""
    r = gas.gas_hydrate(p=1000, degf=60, sg=0.65, p_res=3000, degf_res=200)
    assert r.water_vaporized_res > r.water_vaporized_op, \
        "Reservoir (hotter) should have more vaporized water"
    expected_condensed = r.water_vaporized_res - r.water_vaporized_op
    assert abs(r.water_condensed - expected_condensed) < 1e-15, \
        f"condensed {r.water_condensed} != res-op {expected_condensed}"
    assert r.total_liquid_water == r.water_condensed + r.free_water, \
        "total_liquid = condensed + free_water"

def test_hydrate_water_balance_with_free_water():
    """Free water (additional_water) adds to total liquid, not to vaporized"""
    r = gas.gas_hydrate(p=1000, degf=60, sg=0.65, p_res=3000, degf_res=200,
                         additional_water=0.5)
    assert r.free_water == 0.5, f"free_water should echo additional_water input"
    assert abs(r.total_liquid_water - (r.water_condensed + 0.5)) < 1e-15, \
        "total_liquid = condensed + free"

def test_hydrate_water_content_soreide_whitson():
    """Water content uses SoreideWhitson when composition provided"""
    r = gas.gas_hydrate(p=2000, degf=80, sg=0.7, co2=0.1)
    wc_danesh = gas.gas_water_content(p=2000, degf=80, salinity=0)
    assert r.water_vaporized_op != wc_danesh, \
        "With CO2 composition, should use SoreideWhitson, not Danesh"
    assert r.water_vaporized_op > 0, "Vaporized water should be positive"

def test_hydrate_water_balance_frozen_baselines():
    """Frozen baselines for water balance fields"""
    r1 = gas.gas_hydrate(p=1000, degf=60, sg=0.65)
    assert abs(r1.water_vaporized_op - _HYDRATE_NEW_BASELINES['danesh_wc_op_1000_60_065']) < 1e-10

    r2 = gas.gas_hydrate(p=2000, degf=80, sg=0.7, co2=0.1)
    assert abs(r2.water_vaporized_op - _HYDRATE_NEW_BASELINES['sw_wc_op_2000_80_070_co2']) < 1e-10

    r3 = gas.gas_hydrate(p=1000, degf=60, sg=0.65, p_res=3000, degf_res=200)
    assert abs(r3.water_vaporized_res - _HYDRATE_NEW_BASELINES['wc_res_3000_200']) < 1e-10
    assert abs(r3.water_condensed - _HYDRATE_NEW_BASELINES['condensed_3000_200_to_1000_60']) < 1e-10

def test_hydrate_meoh_capping():
    """MEOH should cap at 25wt% for high-subcooling scenario"""
    r = gas.gas_hydrate(p=2000, degf=80, sg=0.7, hydmethod='MOTIEE', inhibitor_type='MEOH', inhibitor_wt_pct=25)
    assert r.required_inhibitor_wt_pct == 25.0, \
        f"MEOH required should be capped at 25.0, got {r.required_inhibitor_wt_pct}"
    assert r.max_inhibitor_wt_pct == 25.0, f"Max should be 25.0, got {r.max_inhibitor_wt_pct}"
    assert r.inhibitor_underdosed is True, "Should be underdosed"

def test_hydrate_meg_no_capping():
    """MEG should NOT cap at moderate conditions"""
    # At 1000 psia, 90F with MOTIEE — small subcooling, MEG required should be < 70%
    r = gas.gas_hydrate(p=1000, degf=90, sg=0.65, hydmethod='MOTIEE', inhibitor_type='MEG')
    assert r.required_inhibitor_wt_pct < r.max_inhibitor_wt_pct, \
        f"MEG required {r.required_inhibitor_wt_pct} should be < max {r.max_inhibitor_wt_pct}"
    assert r.inhibitor_underdosed is False, "Should NOT be underdosed"

def test_hydrate_all_max_wt_pct():
    """All 5 inhibitors report correct max_inhibitor_wt_pct"""
    expected = {'MEOH': 25.0, 'MEG': 70.0, 'DEG': 70.0, 'TEG': 50.0, 'ETOH': 30.0}
    for inh_name, exp_max in expected.items():
        r = gas.gas_hydrate(p=1000, degf=60, sg=0.65, inhibitor_type=inh_name)
        assert r.max_inhibitor_wt_pct == exp_max, \
            f"{inh_name}: max_wt_pct {r.max_inhibitor_wt_pct} != {exp_max}"

def test_hydrate_injection_rate_algebra():
    """Injection rate matches hand calculation based on total_liquid_water"""
    # Need reservoir conditions for non-zero condensation and injection rate
    r = gas.gas_hydrate(p=2000, degf=60, sg=0.7, inhibitor_type='MEOH',
                         p_res=3000, degf_res=200)
    # Hand calculation: inhibitor treats total liquid water
    total_liquid = float(r.total_liquid_water)
    liquid_mass_lb = total_liquid * 350.2  # lb/MMscf
    w_frac = r.required_inhibitor_wt_pct / 100.0
    expected_mass = liquid_mass_lb * w_frac / (1.0 - w_frac)
    expected_vol = expected_mass / (0.791 * 8.34540445)
    assert abs(r.inhibitor_mass_rate - expected_mass) < 1e-8, \
        f"Mass rate {r.inhibitor_mass_rate} != hand calc {expected_mass}"
    assert abs(r.inhibitor_vol_rate - expected_vol) < 1e-8, \
        f"Vol rate {r.inhibitor_vol_rate} != hand calc {expected_vol}"

def test_hydrate_injection_rate_frozen_baselines():
    """Frozen baselines for injection rates"""
    r = gas.gas_hydrate(p=1000, degf=60, sg=0.65, hydmethod='MOTIEE', inhibitor_type='MEG',
                         p_res=3000, degf_res=200)
    assert abs(r.inhibitor_mass_rate - _HYDRATE_NEW_BASELINES['meg_mass_rate_res_to_op']) < 1e-8, \
        f"Mass rate: {r.inhibitor_mass_rate} != {_HYDRATE_NEW_BASELINES['meg_mass_rate_res_to_op']}"
    assert abs(r.inhibitor_vol_rate - _HYDRATE_NEW_BASELINES['meg_vol_rate_res_to_op']) < 1e-8, \
        f"Vol rate: {r.inhibitor_vol_rate} != {_HYDRATE_NEW_BASELINES['meg_vol_rate_res_to_op']}"

    # MEOH capped with reservoir
    r2 = gas.gas_hydrate(p=2000, degf=60, sg=0.7, hydmethod='MOTIEE', inhibitor_type='MEOH',
                          p_res=3000, degf_res=200)
    assert abs(r2.inhibitor_mass_rate - _HYDRATE_NEW_BASELINES['meoh_mass_rate_res_3000_200']) < 1e-8

def test_hydrate_injection_rate_increases_with_additional_water():
    """Injection rate should increase with additional_water (free water)"""
    r0 = gas.gas_hydrate(p=2000, degf=80, sg=0.7, hydmethod='MOTIEE', inhibitor_type='MEOH',
                          p_res=3000, degf_res=200)
    r1 = gas.gas_hydrate(p=2000, degf=80, sg=0.7, hydmethod='MOTIEE', inhibitor_type='MEOH',
                          p_res=3000, degf_res=200, additional_water=1.0)
    assert r1.inhibitor_mass_rate > r0.inhibitor_mass_rate, \
        "Mass rate should increase with additional_water"
    assert r1.free_water == 1.0, "free_water should reflect additional_water input"

def test_hydrate_injection_rate_zero_outside_window():
    """Injection rate should be 0 when outside hydrate window"""
    r = gas.gas_hydrate(p=200, degf=150, sg=0.6, inhibitor_type='MEG',
                         inhibitor_wt_pct=20, p_res=1000, degf_res=200)
    assert not r.in_hydrate_window, "Should be outside hydrate window"
    assert r.inhibitor_mass_rate == 0.0, f"Mass rate should be 0, got {r.inhibitor_mass_rate}"
    assert r.inhibitor_vol_rate == 0.0, f"Vol rate should be 0, got {r.inhibitor_vol_rate}"

def test_hydrate_injection_rate_zero_no_inhibitor():
    """Injection rate should be 0 when no inhibitor specified"""
    r = gas.gas_hydrate(p=2000, degf=80, sg=0.7, p_res=3000, degf_res=200)
    assert r.inhibitor_mass_rate == 0.0, "Mass rate should be 0 without inhibitor"
    assert r.inhibitor_vol_rate == 0.0, "Vol rate should be 0 without inhibitor"

def test_hydrate_injection_rate_zero_no_liquid():
    """Injection rate should be 0 when no liquid water (no reservoir, no free water)"""
    r = gas.gas_hydrate(p=2000, degf=80, sg=0.7, inhibitor_type='MEOH')
    assert r.total_liquid_water == 0.0, "No liquid without reservoir or free water"
    assert r.inhibitor_mass_rate == 0.0, "No injection rate without liquid water"

def test_hydrate_injection_rate_free_water_only():
    """Injection rate should work with free water even without reservoir conditions"""
    r = gas.gas_hydrate(p=1000, degf=60, sg=0.65, inhibitor_type='MEG',
                         additional_water=0.5)
    assert r.water_condensed == 0.0, "No condensation without reservoir P,T"
    assert r.free_water == 0.5
    assert r.total_liquid_water == 0.5
    assert r.inhibitor_mass_rate > 0, "Should have injection rate from free water"

def test_hydrate_reservoir_pt_does_not_affect_hft():
    """Reservoir P,T affects water balance but not hydrate assessment"""
    r_default = gas.gas_hydrate(p=1000, degf=60, sg=0.65)
    r_res = gas.gas_hydrate(p=1000, degf=60, sg=0.65, p_res=3000, degf_res=200)
    assert abs(r_res.hft - r_default.hft) < 1e-10, \
        "HFT should not change with reservoir P,T"
    assert abs(r_res.hfp - r_default.hfp) < 1e-10, \
        "HFP should not change with reservoir P,T"

def test_hydrate_metric_water_balance():
    """Metric water balance fields should be consistent with oilfield conversion"""
    from pyrestoolbox.constants import STB_PER_MMSCF_TO_SM3_PER_SM3
    r_f = gas.gas_hydrate(p=1000, degf=60, sg=0.65, p_res=3000, degf_res=200)
    r_m = gas.gas_hydrate(p=1000 * 0.0689475729, degf=(60 - 32) * 5 / 9, sg=0.65,
                           p_res=3000 * 0.0689475729, degf_res=(200 - 32) * 5 / 9, metric=True)
    for field in ['water_vaporized_res', 'water_vaporized_op', 'water_condensed', 'total_liquid_water']:
        f_val = getattr(r_f, field)
        m_val = getattr(r_m, field)
        expected = f_val * STB_PER_MMSCF_TO_SM3_PER_SM3
        assert abs(m_val - expected) / expected < 1e-4, \
            f"Metric {field}: {m_val} != converted {expected}"

def test_hydrate_metric_injection_rate():
    """Metric injection rates should be consistent with oilfield conversion"""
    from pyrestoolbox.constants import LB_PER_MMSCF_TO_KG_PER_SM3, GAL_PER_MMSCF_TO_L_PER_SM3
    r_f = gas.gas_hydrate(p=1000, degf=60, sg=0.65, inhibitor_type='MEG',
                           p_res=3000, degf_res=200)
    r_m = gas.gas_hydrate(p=1000 * 0.0689475729, degf=(60 - 32) * 5 / 9, sg=0.65,
                           inhibitor_type='MEG', p_res=3000 * 0.0689475729,
                           degf_res=(200 - 32) * 5 / 9, metric=True)
    mass_kg = r_f.inhibitor_mass_rate * LB_PER_MMSCF_TO_KG_PER_SM3
    vol_l = r_f.inhibitor_vol_rate * GAL_PER_MMSCF_TO_L_PER_SM3
    assert abs(r_m.inhibitor_mass_rate - mass_kg) / mass_kg < 1e-4, \
        f"Metric mass rate {r_m.inhibitor_mass_rate} != converted {mass_kg}"
    assert abs(r_m.inhibitor_vol_rate - vol_l) / vol_l < 1e-4, \
        f"Metric vol rate {r_m.inhibitor_vol_rate} != converted {vol_l}"

def test_hydrate_composition_validation():
    """Composition inputs should be validated"""
    # Negative mole fraction
    for param in ['co2', 'h2s', 'n2', 'h2']:
        try:
            gas.gas_hydrate(p=1000, degf=60, sg=0.65, **{param: -0.1})
            raise AssertionError(f"Expected ValueError for negative {param}")
        except ValueError:
            pass
    # Sum > 1
    try:
        gas.gas_hydrate(p=1000, degf=60, sg=0.65, co2=0.5, h2s=0.6)
        raise AssertionError("Expected ValueError for sum > 1")
    except ValueError:
        pass
    # Negative additional_water
    try:
        gas.gas_hydrate(p=1000, degf=60, sg=0.65, additional_water=-1)
        raise AssertionError("Expected ValueError for negative additional_water")
    except ValueError:
        pass

def test_hydrate_new_fields_no_inhibitor():
    """New fields should have sensible defaults when no inhibitor"""
    r = gas.gas_hydrate(p=1000, degf=60, sg=0.65)
    assert r.max_inhibitor_wt_pct == 0.0, "max should be 0 without inhibitor"
    assert r.inhibitor_underdosed is False, "underdosed should be False without inhibitor"
    assert r.water_vaporized_op > 0, "vaporized water at operating should be positive"
    assert r.inhibitor_mass_rate == 0.0, "mass rate should be 0 without inhibitor"
    assert r.inhibitor_vol_rate == 0.0, "vol rate should be 0 without inhibitor"


if __name__ == '__main__':
    print("=" * 70)
    print("GAS MODULE VALIDATION TESTS")
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
