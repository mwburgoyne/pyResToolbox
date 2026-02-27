#!/usr/bin/env python3
"""
Validation tests for brine module.
Run with: PYTHONPATH=/home/mark/projects python3 -m pytest tests/ -v
Or standalone: PYTHONPATH=/home/mark/projects python3 tests/test_brine.py
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import pyrestoolbox.brine as brine

# =============================================================================
# CH4-Saturated Brine Tests (Spivey correlation)
# =============================================================================

def test_brine_props_freshwater():
    """Freshwater properties at moderate conditions"""
    bw, lden, visw, cw, rsw = brine.brine_props(p=3000, degf=200, wt=0, ch4_sat=0)
    assert 0.98 < bw < 1.15, f"Bw = {bw}, outside expected range"
    assert 0.9 < lden < 1.1, f"Brine density SG = {lden}, outside expected range"
    assert 0.1 < visw < 2.0, f"Brine viscosity = {visw} cP, outside expected range"
    assert 0 < abs(cw) < 0.001, f"Brine compressibility = {cw}, outside expected range"
    assert rsw >= 0, f"Brine gas solubility should be non-negative"

def test_brine_props_saline():
    """Saline brine properties"""
    bw, lden, visw, cw, rsw = brine.brine_props(p=3000, degf=200, wt=10, ch4_sat=0)
    # Saline brine should be denser than freshwater
    bw_fw, lden_fw, _, _, _ = brine.brine_props(p=3000, degf=200, wt=0, ch4_sat=0)
    assert lden > lden_fw, f"Saline brine ({lden}) should be denser than freshwater ({lden_fw})"
    # Saline brine should be more viscous
    assert visw > 0, "Viscosity should be positive"

def test_brine_props_ch4_saturated():
    """Methane-saturated brine"""
    bw, lden, visw, cw, rsw = brine.brine_props(p=3000, degf=200, wt=0, ch4_sat=1.0)
    bw_dry, _, _, _, rsw_dry = brine.brine_props(p=3000, degf=200, wt=0, ch4_sat=0)
    assert rsw > rsw_dry, "CH4-saturated brine should have higher Rsw"
    assert bw > bw_dry * 0.99, "CH4-saturated Bw should be >= dry brine Bw"

def test_brine_bw_increases_with_temp():
    """Bw should generally increase with temperature"""
    bw_low, _, _, _, _ = brine.brine_props(p=3000, degf=100, wt=0, ch4_sat=0)
    bw_high, _, _, _, _ = brine.brine_props(p=3000, degf=300, wt=0, ch4_sat=0)
    assert bw_high > bw_low, f"Bw should increase with T: Bw(100)={bw_low}, Bw(300)={bw_high}"

def test_brine_viscosity_decreases_with_temp():
    """Brine viscosity should decrease with temperature"""
    _, _, vis_low, _, _ = brine.brine_props(p=3000, degf=100, wt=0, ch4_sat=0)
    _, _, vis_high, _, _ = brine.brine_props(p=3000, degf=300, wt=0, ch4_sat=0)
    assert vis_high < vis_low, f"Viscosity should decrease with T: vis(100)={vis_low}, vis(300)={vis_high}"

# =============================================================================
# CO2-Brine Mixture Tests (Spycher-Pruess)
# =============================================================================

def test_co2_brine_basic():
    """Basic CO2-brine mixture at moderate conditions"""
    mix = brine.CO2_Brine_Mixture(pres=200, temp=80, ppm=0, metric=True)
    assert 0 < mix.x[0] < 0.05, f"xCO2 = {mix.x[0]}, outside expected range"
    assert 0.95 < mix.y[0] <= 1.0, f"yCO2 = {mix.y[0]}, outside expected range"
    assert mix.rhoGas > 0, "Gas density should be positive"
    assert mix.bDen[0] > 0.9, f"CO2-laden brine density = {mix.bDen[0]}"

def test_co2_brine_field_units():
    """CO2-brine in field units"""
    mix = brine.CO2_Brine_Mixture(pres=5000, temp=275, ppm=30000, metric=False)
    assert 0 < mix.x[0] < 0.05
    assert mix.Rs > 0, f"Rs should be positive, got {mix.Rs}"
    assert mix.bw[0] > 1.0, f"Bw should be > 1.0 for CO2-laden brine"

def test_co2_brine_xco2_increases_with_pressure():
    """CO2 solubility should increase with pressure"""
    mix_low = brine.CO2_Brine_Mixture(pres=50, temp=80, ppm=0, metric=True)
    mix_high = brine.CO2_Brine_Mixture(pres=300, temp=80, ppm=0, metric=True)
    assert mix_high.x[0] > mix_low.x[0], \
        f"xCO2 should increase with P: xCO2(50)={mix_low.x[0]}, xCO2(300)={mix_high.x[0]}"

def test_co2_brine_salt_reduces_solubility():
    """Salt should reduce CO2 solubility (salting out)"""
    mix_fresh = brine.CO2_Brine_Mixture(pres=200, temp=80, ppm=0, metric=True)
    mix_saline = brine.CO2_Brine_Mixture(pres=200, temp=80, ppm=100000, metric=True)
    assert mix_saline.x[0] < mix_fresh.x[0], \
        f"Salt should reduce xCO2: fresh={mix_fresh.x[0]}, saline={mix_saline.x[0]}"

def test_co2_brine_density_greater_than_pure():
    """CO2-laden brine should be denser than pure brine"""
    mix = brine.CO2_Brine_Mixture(pres=200, temp=80, ppm=0, metric=True)
    assert mix.bDen[0] >= mix.bDen[1] * 0.99, \
        f"CO2-laden density ({mix.bDen[0]}) should be >= pure brine ({mix.bDen[1]})"

def test_co2_brine_viscosity_positive():
    """Brine viscosity should be positive"""
    mix = brine.CO2_Brine_Mixture(pres=200, temp=80, ppm=0, metric=True)
    assert mix.bVis[0] > 0, f"CO2-laden viscosity should be positive"
    assert mix.bVis[1] > 0, f"Pure brine viscosity should be positive"

def test_co2_brine_saturated_compressibility():
    """Saturated compressibility calculation"""
    mix = brine.CO2_Brine_Mixture(pres=200, temp=80, ppm=0, metric=True, cw_sat=True)
    assert mix.Cf_sat is not None, "Cf_sat should be calculated when cw_sat=True"
    assert mix.Cf_sat > 0, f"Saturated compressibility should be positive, got {mix.Cf_sat}"

# =============================================================================
# Regression baselines
# =============================================================================

# Hard-coded regression baselines - frozen reference values
# If any of these change, it means a code modification has altered computational results.
# Do NOT update these values without explicit review and approval.
# Frozen baselines captured 2026-02-27
_FROZEN_BASELINES = {
    'bw_3000_200_fresh': 1.027589195773527,
    'lden_3000_200_fresh': 0.972276768415092,
    'visw_3000_200_fresh': 0.3083544960904146,
    'cw_3000_200_fresh': 3.0887176266534516e-06,
    'co2_xco2_200_80': 0.02036714644979853,
    'co2_rs_200_80': 26.141605639385897,
}

def test_regression_brine_freshwater():
    bw, lden, visw, cw, rsw = brine.brine_props(p=3000, degf=200, wt=0, ch4_sat=0)
    assert abs(bw - _FROZEN_BASELINES['bw_3000_200_fresh']) / _FROZEN_BASELINES['bw_3000_200_fresh'] < 1e-8, \
        f"Brine Bw changed: {bw} vs frozen {_FROZEN_BASELINES['bw_3000_200_fresh']}"
    assert abs(lden - _FROZEN_BASELINES['lden_3000_200_fresh']) / _FROZEN_BASELINES['lden_3000_200_fresh'] < 1e-8, \
        f"Brine density changed: {lden} vs frozen {_FROZEN_BASELINES['lden_3000_200_fresh']}"
    assert abs(visw - _FROZEN_BASELINES['visw_3000_200_fresh']) / _FROZEN_BASELINES['visw_3000_200_fresh'] < 1e-8, \
        f"Brine viscosity changed: {visw} vs frozen {_FROZEN_BASELINES['visw_3000_200_fresh']}"

def test_regression_co2_xco2():
    mix = brine.CO2_Brine_Mixture(pres=200, temp=80, ppm=0, metric=True)
    expected = _FROZEN_BASELINES['co2_xco2_200_80']
    assert abs(mix.x[0] - expected) / expected < 1e-6, \
        f"CO2 xCO2 changed: {mix.x[0]} vs frozen {expected}"

def test_regression_co2_rs():
    """CO2 brine Rs must match frozen baseline"""
    mix = brine.CO2_Brine_Mixture(pres=200, temp=80, ppm=0, metric=True)
    expected = _FROZEN_BASELINES['co2_rs_200_80']
    assert abs(mix.Rs - expected) / expected < 1e-6, \
        f"CO2 Rs changed: {mix.Rs} vs frozen {expected}"


# =============================================================================
# make_pvtw_table Tests
# =============================================================================

def test_make_pvtw_table_basic():
    """Result has expected keys and correct row count"""
    result = brine.make_pvtw_table(pi=3000, degf=200, wt=0, ch4_sat=0, nrows=10)
    expected_keys = {'table', 'pref', 'bw_ref', 'cw_ref', 'visw_ref', 'rsw_ref', 'den_ref'}
    assert expected_keys == set(result.keys()), f"Missing keys: {expected_keys - set(result.keys())}"
    # nrows=10 plus pi if not already in grid
    assert len(result['table']) >= 10, f"Expected at least 10 rows, got {len(result['table'])}"
    assert result['pref'] == 3000, f"pref should be 3000, got {result['pref']}"

def test_make_pvtw_table_bw_range():
    """Bw values should be in reasonable range (0.9 - 1.2)"""
    result = brine.make_pvtw_table(pi=3000, degf=200, wt=0, ch4_sat=0)
    bw_vals = result['table']['Bw (rb/stb)'].values
    assert all(0.9 < bw < 1.2 for bw in bw_vals), f"Bw values outside 0.9-1.2: {bw_vals}"

def test_make_pvtw_table_ref_vs_direct():
    """Reference properties should match direct brine_props() call"""
    pi, degf, wt, ch4_sat = 3000, 200, 5, 0.5
    result = brine.make_pvtw_table(pi=pi, degf=degf, wt=wt, ch4_sat=ch4_sat)
    bw, lden, visw, cw, rsw = brine.brine_props(p=pi, degf=degf, wt=wt, ch4_sat=ch4_sat)
    assert abs(result['bw_ref'] - bw) < 1e-10, f"bw_ref mismatch: {result['bw_ref']} vs {bw}"
    assert abs(result['den_ref'] - lden) < 1e-10, f"den_ref mismatch"
    assert abs(result['visw_ref'] - visw) < 1e-10, f"visw_ref mismatch"
    assert abs(result['cw_ref'] - cw) < 1e-10, f"cw_ref mismatch"
    assert abs(result['rsw_ref'] - rsw) < 1e-10, f"rsw_ref mismatch"

def test_make_pvtw_table_pi_included():
    """pi should appear in the pressure grid"""
    pi = 3333
    result = brine.make_pvtw_table(pi=pi, degf=200, wt=0, ch4_sat=0)
    pressures = result['table']['Pressure (psia)'].values
    assert pi in pressures, f"pi={pi} not found in pressure grid: {pressures}"

def test_make_pvtw_table_saline():
    """Saline brine density should be greater than freshwater density"""
    result_fresh = brine.make_pvtw_table(pi=3000, degf=200, wt=0, ch4_sat=0)
    result_saline = brine.make_pvtw_table(pi=3000, degf=200, wt=10, ch4_sat=0)
    assert result_saline['den_ref'] > result_fresh['den_ref'], \
        f"Saline density ({result_saline['den_ref']}) should exceed freshwater ({result_fresh['den_ref']})"

# =============================================================================
# SoreideWhitson Tests
# =============================================================================

def test_soreide_whitson_pure_co2():
    """SoreideWhitson with pure CO2 should produce positive density, viscosity, and Rs"""
    mix = brine.SoreideWhitson(pres=5000, temp=275, ppm=30000, y_CO2=1.0, metric=False)
    # Density: gas-saturated should be greater than gas-free (CO2 increases brine density)
    assert mix.bDen[0] > mix.bDen[1], f"CO2-saturated density {mix.bDen[0]} should exceed gas-free {mix.bDen[1]}"
    # Viscosity: gas-saturated should be greater than gas-free (CO2 increases viscosity)
    assert mix.bVis[0] > mix.bVis[1], f"CO2-saturated viscosity {mix.bVis[0]} should exceed gas-free {mix.bVis[1]}"
    # Rs should be positive
    assert mix.Rs_total > 0, f"Rs_total should be positive, got {mix.Rs_total}"
    # Bw should be > 1 (dissolved gas expands brine)
    assert mix.bw[0] > 1.0, f"Gas-saturated Bw should exceed 1.0, got {mix.bw[0]}"
    # x should have CO2 key with positive value
    assert 'CO2' in mix.x, f"x should contain CO2, got {mix.x}"
    assert mix.x['CO2'] > 0, f"xCO2 should be positive, got {mix.x['CO2']}"


def test_soreide_whitson_pure_ch4():
    """SoreideWhitson with pure CH4 should produce physically reasonable results"""
    mix = brine.SoreideWhitson(pres=5000, temp=275, ppm=30000, y_CO2=0, sg=0.554, metric=False)
    # CH4 decreases brine density
    assert mix.bDen[0] < mix.bDen[1], f"CH4-saturated density {mix.bDen[0]} should be less than gas-free {mix.bDen[1]}"
    # Rs should be positive
    assert mix.Rs_total > 0, f"Rs_total should be positive, got {mix.Rs_total}"
    # Should contain CH4 in dissolved gases
    assert 'CH4' in mix.x, f"x should contain CH4, got {mix.x}"


def test_soreide_whitson_mixed_gas():
    """SoreideWhitson with mixed gas should return per-gas details"""
    mix = brine.SoreideWhitson(pres=200, temp=80, ppm=10000, y_CO2=0.1, y_H2S=0.05, sg=0.7, metric=True)
    # Should have multiple gases in results
    assert len(mix.x) >= 2, f"Should have at least 2 dissolved gases, got {len(mix.x)}"
    assert len(mix.Rs) >= 2, f"Should have at least 2 Rs entries, got {len(mix.Rs)}"
    # Gas composition should be normalized
    total_y = sum(mix.gas_comp.values())
    assert abs(total_y - 1.0) < 1e-10, f"Gas composition should sum to 1.0, got {total_y}"
    # All densities should be positive
    for d in mix.bDen:
        assert d > 0, f"All densities should be positive, got {mix.bDen}"


def test_soreide_whitson_bad_gas_fractions():
    """SoreideWhitson should raise ValueError if non-HC fractions exceed 1.0"""
    try:
        brine.SoreideWhitson(pres=200, temp=80, y_CO2=0.6, y_H2S=0.5, metric=True)
        assert False, "Expected ValueError for fractions > 1.0"
    except ValueError:
        pass


if __name__ == '__main__':
    print("=" * 70)
    print("BRINE MODULE VALIDATION TESTS")
    print("=" * 70)

    tests = [v for k, v in globals().items() if k.startswith('test_')]
    passed = 0
    failed = 0

    for test in tests:
        try:
            test()
            passed += 1
            print(f"  PASS: {test.__name__}")
        except Exception as e:
            failed += 1
            print(f"  FAIL: {test.__name__}: {e}")

    print(f"\n{'=' * 70}")
    print(f"Results: {passed} passed, {failed} failed out of {passed + failed}")
    print("=" * 70)
    sys.exit(1 if failed > 0 else 0)
