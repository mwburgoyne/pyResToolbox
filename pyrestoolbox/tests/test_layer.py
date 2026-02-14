#!/usr/bin/env python3
"""
Validation tests for layer module.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import pyrestoolbox.layer as layer

def test_lorenz2b_zero():
    """Near-zero Lorenz should give small B"""
    B = layer.lorenz2b(0.001)
    assert B > 0
    assert B < 1

def test_lorenz2b_high():
    """High Lorenz should give large B"""
    B = layer.lorenz2b(0.9)
    assert B > 5

def test_lorenz_roundtrip_exp():
    """Lorenz → B → Lorenz roundtrip (exponential)"""
    for L_input in [0.1, 0.3, 0.5, 0.7, 0.9]:
        B = layer.lorenz2b(L_input, 'EXP')
        L_back = layer.lorenzfromb(B, 'EXP')
        assert abs(L_back - L_input) < 0.001, \
            f"EXP roundtrip: L={L_input} → B={B} → L={L_back}"

def test_lorenz_roundtrip_lang():
    """Lorenz → B → Lorenz roundtrip (Langmuir)"""
    for L_input in [0.1, 0.3, 0.5, 0.7, 0.9]:
        B = layer.lorenz2b(L_input, 'LANG')
        L_back = layer.lorenzfromb(B, 'LANG')
        assert abs(L_back - L_input) < 0.001, \
            f"LANG roundtrip: L={L_input} → B={B} → L={L_back}"

def test_lorenz_2_layers_avg():
    """Generated layers should have correct average permeability"""
    k_avg = 100
    for lorenz in [0.2, 0.5, 0.8]:
        k = layer.lorenz_2_layers(lorenz=lorenz, k_avg=k_avg, nlayers=20)
        calculated_avg = np.mean(k)
        assert abs(calculated_avg - k_avg) / k_avg < 0.01, \
            f"Avg perm mismatch: L={lorenz}, avg={calculated_avg}, expected={k_avg}"

def test_lorenz_2_layers_sorted():
    """Layers should be in decreasing order by default"""
    k = layer.lorenz_2_layers(lorenz=0.5, k_avg=100, nlayers=10)
    for i in range(len(k) - 1):
        assert k[i] >= k[i + 1], "Layers should be sorted in decreasing order"

def test_lorenz_2_layers_single():
    """Single layer should return k_avg"""
    k = layer.lorenz_2_layers(lorenz=0.5, k_avg=100, nlayers=1)
    assert len(k) == 1
    assert k[0] == 100

def test_lorenz_2_layers_heterogeneity():
    """Higher Lorenz should give more heterogeneous distribution"""
    k_low = layer.lorenz_2_layers(lorenz=0.2, k_avg=100, nlayers=20)
    k_high = layer.lorenz_2_layers(lorenz=0.8, k_avg=100, nlayers=20)
    cv_low = np.std(k_low) / np.mean(k_low)
    cv_high = np.std(k_high) / np.mean(k_high)
    assert cv_high > cv_low, f"Higher Lorenz should give more variation: CV(0.2)={cv_low}, CV(0.8)={cv_high}"

def test_lorenz_from_flow_fraction():
    """Flow fraction to Lorenz"""
    L = layer.lorenz_from_flow_fraction(kh_frac=0.8, phih_frac=0.5)
    assert 0 < L < 1, f"Lorenz = {L}, should be between 0 and 1"

def test_lorenz_2_flow_frac_homogeneous():
    """Zero Lorenz should give kh_frac = phih_frac"""
    frac = layer.lorenz_2_flow_frac(lorenz=0.001, phih_frac=0.5)
    assert abs(frac - 0.5) < 0.05, f"Homogeneous: flow_frac = {frac}, expected ~0.5"

def test_lorenz_2_flow_frac_heterogeneous():
    """High Lorenz should give kh_frac >> phih_frac"""
    frac = layer.lorenz_2_flow_frac(lorenz=0.8, phih_frac=0.3)
    assert frac > 0.3, f"Heterogeneous: flow_frac = {frac}, should be > phih_frac=0.3"

def test_lorenz_flow_frac_roundtrip():
    """Flow fraction roundtrip"""
    L_input = 0.5
    phih = 0.4
    kh_frac = layer.lorenz_2_flow_frac(lorenz=L_input, phih_frac=phih)
    L_back = layer.lorenz_from_flow_fraction(kh_frac=kh_frac, phih_frac=phih)
    assert abs(L_back - L_input) < 0.01, \
        f"Flow frac roundtrip: L={L_input} → kh={kh_frac} → L={L_back}"


if __name__ == '__main__':
    print("=" * 70)
    print("LAYER MODULE VALIDATION TESTS")
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
