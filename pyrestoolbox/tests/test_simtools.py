#!/usr/bin/env python3
"""
Validation tests for simtools module (rel-perm and flash only;
IX extraction and influence tables require external files/long runtimes).
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import pyrestoolbox.simtools as simtools

# =============================================================================
# Relative Permeability Tables
# =============================================================================

def test_swof_corey():
    """SWOF table with Corey curves"""
    df = simtools.rel_perm_table(rows=20, krtable='SWOF', krfamily='COR',
                                  swc=0.2, swcr=0.2, sorw=0.2, no=2, nw=3)
    assert df is not None, "SWOF table should not be None"
    assert 'Sw' in df.columns
    assert 'Krwo' in df.columns
    assert 'Krow' in df.columns
    # Endpoint checks
    sw = df['Sw'].values
    assert sw[0] <= 0.2 + 0.001, f"Min Sw should be ~swc: {sw[0]}"
    assert sw[-1] >= 0.8 - 0.001, f"Max Sw should be ~1-sorw: {sw[-1]}"
    # Krw at Swc should be 0
    assert df['Krwo'].iloc[0] < 0.001
    # Kro at 1-Sorw should be 0
    assert df['Krow'].iloc[-1] < 0.001

def test_swof_let():
    """SWOF table with LET curves"""
    df = simtools.rel_perm_table(rows=20, krtable='SWOF', krfamily='LET',
                                  swc=0.2, swcr=0.2, sorw=0.2,
                                  Lw=2, Ew=1, Tw=2, Lo=2, Eo=1, To=2)
    assert df is not None
    assert len(df) >= 10

def test_sgof_corey():
    """SGOF table with Corey curves"""
    df = simtools.rel_perm_table(rows=20, krtable='SGOF', krfamily='COR',
                                  swc=0.2, sorg=0.1, no=2, ng=2)
    assert df is not None
    assert 'Sg' in df.columns
    assert 'Krgo' in df.columns
    assert 'Krog' in df.columns
    # Sg should start at 0
    assert df['Sg'].iloc[0] < 0.001
    # Krog at max Sg should be ~0
    assert df['Krog'].iloc[-1] < 0.01

def test_sgwfn_corey():
    """SGWFN table with Corey curves"""
    df = simtools.rel_perm_table(rows=20, krtable='SGWFN', krfamily='COR',
                                  swc=0.2, sgcr=0.05, nw=2, ng=2)
    assert df is not None
    assert 'Sg' in df.columns
    assert 'Krgw' in df.columns
    assert 'Krwg' in df.columns

def test_kr_values_bounded():
    """All Kr values should be between 0 and 1"""
    df = simtools.rel_perm_table(rows=20, krtable='SWOF', krfamily='COR',
                                  swc=0.15, swcr=0.15, sorw=0.25,
                                  no=2, nw=3, kromax=0.9, krwmax=0.5)
    assert all(df['Krwo'] >= -0.001), "Krw should be non-negative"
    assert all(df['Krwo'] <= 1.001), "Krw should be <= 1"
    assert all(df['Krow'] >= -0.001), "Kro should be non-negative"
    assert all(df['Krow'] <= 1.001), "Kro should be <= 1"

def test_kr_monotonic():
    """Krw should be monotonically increasing, Kro decreasing with Sw"""
    df = simtools.rel_perm_table(rows=30, krtable='SWOF', krfamily='COR',
                                  swc=0.2, swcr=0.2, sorw=0.2, no=2, nw=3)
    krw = df['Krwo'].values
    kro = df['Krow'].values
    for i in range(len(krw) - 1):
        assert krw[i + 1] >= krw[i] - 1e-10, f"Krw not monotonic at index {i}"
        assert kro[i + 1] <= kro[i] + 1e-10, f"Kro not monotonic at index {i}"

# =============================================================================
# LET and Corey functions directly
# =============================================================================

def test_let_endpoints():
    """LET function at s=0 should be 0, at s=1 should be 1"""
    s = np.array([0.0, 1.0])
    kr = simtools.LET(s, L=2, E=1, T=2)
    assert abs(kr[0]) < 1e-10, f"LET(0) = {kr[0]}, expected 0"
    assert abs(kr[1] - 1.0) < 1e-10, f"LET(1) = {kr[1]}, expected 1"

def test_corey_endpoints():
    """Corey function at s=0 should be 0, at s=1 should be 1"""
    s = np.array([0.0, 1.0])
    kr = simtools.corey(s, n=2)
    assert abs(kr[0]) < 1e-10
    assert abs(kr[1] - 1.0) < 1e-10

def test_corey_linear():
    """Corey with n=1 should give linear relationship"""
    s = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    kr = simtools.corey(s, n=1)
    np.testing.assert_allclose(kr, s, atol=1e-10)

# =============================================================================
# Rachford-Rice Solver
# =============================================================================

def test_rr_solver_basic():
    """Basic flash calculation"""
    zi = np.array([0.5, 0.3, 0.2])
    ki = np.array([5.0, 1.5, 0.1])
    N_it, yi, xi, V, L = simtools.rr_solver(zi, ki)
    assert 0 <= V <= 1, f"V = {V}, should be between 0 and 1"
    assert abs(V + L - 1.0) < 1e-10, f"V + L = {V + L}, should be 1.0"
    assert all(yi >= 0) and all(yi <= 1), f"yi out of range: {yi}"
    assert all(xi >= 0) and all(xi <= 1), f"xi out of range: {xi}"
    # Material balance check
    for i in range(len(zi)):
        mb = V * yi[i] + L * xi[i]
        assert abs(mb - zi[i]) < 0.01, f"Material balance error for component {i}: {mb} vs {zi[i]}"

def test_rr_solver_ki_check():
    """K-value consistency: yi/xi should approximate Ki"""
    zi = np.array([0.4, 0.3, 0.2, 0.1])
    ki = np.array([10.0, 2.0, 0.5, 0.01])
    N_it, yi, xi, V, L = simtools.rr_solver(zi, ki)
    for i in range(len(zi)):
        if xi[i] > 1e-10:
            ki_calc = yi[i] / xi[i]
            assert abs(ki_calc - ki[i]) / ki[i] < 0.05, \
                f"K-value mismatch for comp {i}: calc={ki_calc}, expected={ki[i]}"

def test_rr_solver_convergence():
    """Solver should converge in reasonable iterations"""
    zi = np.array([0.5, 0.3, 0.2])
    ki = np.array([5.0, 1.5, 0.1])
    N_it, _, _, _, _ = simtools.rr_solver(zi, ki)
    assert N_it < 50, f"Solver took {N_it} iterations, expected < 50"

def test_rr_solver_near_dew():
    """Near dew point conditions (small V)"""
    zi = np.array([0.1, 0.3, 0.6])
    ki = np.array([3.0, 1.1, 0.8])
    N_it, yi, xi, V, L = simtools.rr_solver(zi, ki)
    assert 0 <= V <= 1
    assert abs(V + L - 1.0) < 1e-10


if __name__ == '__main__':
    print("=" * 70)
    print("SIMTOOLS MODULE VALIDATION TESTS")
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
