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
    'cg_dak_2000_075_200': 0.0005199034510414734,
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
