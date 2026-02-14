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

_BASELINES = None

def capture_brine_baselines():
    b = {}
    b['brine_3000_200_0'] = brine.brine_props(p=3000, degf=200, wt=0, ch4_sat=0)
    b['brine_3000_200_10'] = brine.brine_props(p=3000, degf=200, wt=10, ch4_sat=0)
    mix = brine.CO2_Brine_Mixture(pres=200, temp=80, ppm=0, metric=True)
    b['co2_xco2_200_80'] = mix.x[0]
    b['co2_rs_200_80'] = mix.Rs
    return b

def get_baselines():
    global _BASELINES
    if _BASELINES is None:
        _BASELINES = capture_brine_baselines()
    return _BASELINES

def test_regression_brine_freshwater():
    b = get_baselines()
    bw, lden, visw, cw, rsw = brine.brine_props(p=3000, degf=200, wt=0, ch4_sat=0)
    bw_b, lden_b, visw_b, cw_b, rsw_b = b['brine_3000_200_0']
    assert abs(bw - bw_b) < 1e-10
    assert abs(lden - lden_b) < 1e-10
    assert abs(visw - visw_b) < 1e-10

def test_regression_co2_xco2():
    b = get_baselines()
    mix = brine.CO2_Brine_Mixture(pres=200, temp=80, ppm=0, metric=True)
    assert abs(mix.x[0] - b['co2_xco2_200_80']) < 1e-8


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
