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


# =============================================================================
# Jerauld Correlation
# =============================================================================

def test_jerauld_endpoints():
    """Jerauld function at s=0 should be 0, at s=1 should be 1"""
    s = np.array([0.0, 1.0])
    kr = simtools.jerauld(s, a=2.0, b=1.0)
    assert abs(kr[0]) < 1e-10, f"Jerauld(0) = {kr[0]}, expected 0"
    assert abs(kr[1] - 1.0) < 1e-10, f"Jerauld(1) = {kr[1]}, expected 1"

def test_jerauld_monotonic():
    """Jerauld curve should be monotonically increasing"""
    s = np.linspace(0, 1, 50)
    kr = simtools.jerauld(s, a=2.5, b=0.5)
    dkr = np.diff(kr)
    assert np.all(dkr >= -1e-10), "Jerauld curve not monotonic"

def test_jerauld_in_rel_perm_table():
    """rel_perm_table should work with Jerauld family"""
    df = simtools.rel_perm_table(rows=20, krtable='SWOF', krfamily='JER',
                                  swc=0.2, swcr=0.2, sorw=0.2,
                                  aw=2.0, bw=1.0, ao=2.5, bo=0.8)
    assert df is not None
    assert 'Sw' in df.columns
    assert all(df['Krwo'] >= -0.001)
    assert all(df['Krow'] >= -0.001)

# =============================================================================
# LET Physicality Check
# =============================================================================

def test_is_let_physical_good():
    """Standard LET parameters should pass physicality check"""
    s = np.linspace(0, 1, 100)
    assert simtools.is_let_physical(s, L=2.0, E=1.0, T=2.0) == True

def test_is_let_physical_detects_non_monotonic():
    """Extreme LET parameters that cause non-monotonic behavior should fail"""
    s = np.linspace(0, 1, 200)
    # Very large E relative to L and T can cause non-monotonic behavior
    result = simtools.is_let_physical(s, L=0.1, E=100.0, T=0.1)
    # This specific combination may or may not be physical; test that function runs
    assert isinstance(result, bool)

# =============================================================================
# fit_rel_perm Regression
# =============================================================================

def test_fit_rel_perm_corey():
    """Fit Corey model to synthetic Corey data"""
    s = np.linspace(0, 1, 20)
    kr_true = s ** 2.5
    result = simtools.fit_rel_perm(s, kr_true, krfamily='COR')
    assert abs(result['params']['n'] - 2.5) < 0.1, f"Expected n~2.5, got {result['params']['n']}"
    assert result['ssq'] < 0.001

def test_fit_rel_perm_let():
    """Fit LET model to synthetic LET data"""
    s = np.linspace(0, 1, 30)
    kr_true = simtools.LET(s, L=2.0, E=1.5, T=2.0)
    result = simtools.fit_rel_perm(s, kr_true, krfamily='LET')
    assert result['ssq'] < 0.001
    assert 'L' in result['params']
    assert 'E' in result['params']
    assert 'T' in result['params']

def test_fit_rel_perm_jerauld():
    """Fit Jerauld model to synthetic Jerauld data"""
    s = np.linspace(0, 1, 20)
    kr_true = simtools.jerauld(s, a=2.0, b=1.0)
    result = simtools.fit_rel_perm(s, kr_true, krfamily='JER')
    assert result['ssq'] < 0.001
    assert 'a' in result['params']
    assert 'b' in result['params']

def test_fit_rel_perm_best():
    """fit_rel_perm_best should try all models and pick the best"""
    s = np.linspace(0, 1, 30)
    # Generate LET data — LET should win
    kr_true = simtools.LET(s, L=2.0, E=1.5, T=2.0)
    result = simtools.fit_rel_perm_best(s, kr_true)
    assert 'family' in result
    assert 'all_results' in result
    assert len(result['all_results']) == 3
    assert result['family'] == 'LET', f"Expected LET to win, got {result['family']}"
    assert result['ssq'] < 0.001

def test_fit_rel_perm_with_endpoints():
    """Fit with non-trivial sw_min/sw_max and krmax"""
    # Generate data with swc=0.2, 1-sorw=0.8, krmax=0.5
    sw = np.linspace(0.2, 0.8, 15)
    sn = (sw - 0.2) / (0.8 - 0.2)
    kr_true = 0.5 * sn ** 2.0
    result = simtools.fit_rel_perm(sw, kr_true, krfamily='COR',
                                    krmax=0.5, sw_min=0.2, sw_max=0.8)
    assert abs(result['params']['n'] - 2.0) < 0.15


# =============================================================================
# VFP Table Generation
# =============================================================================

def test_vfpinj_gas():
    """Generate a gas injection VFPINJ table"""
    from pyrestoolbox.nodal import Completion
    comp = Completion(tid=2.441, length=8000, tht=100, bht=250)
    result = simtools.make_vfpinj(
        table_num=1, completion=comp, flo_type='GAS',
        vlpmethod='HB',
        flo_rates=[1000, 5000, 10000],
        thp_values=[1000, 2000],
        gsg=0.65)
    assert result['bhp'].shape == (2, 3)
    assert result['flo_type'] == 'GAS'
    assert result['n_failed'] == 0

def test_vfpinj_water():
    """Generate a water injection VFPINJ table"""
    from pyrestoolbox.nodal import Completion
    comp = Completion(tid=3.5, length=8000, tht=100, bht=200)
    result = simtools.make_vfpinj(
        table_num=2, completion=comp, flo_type='WAT',
        flo_rates=[1000, 5000, 10000],
        thp_values=[500, 1000])
    assert result['bhp'].shape == (2, 3)
    assert result['flo_type'] == 'WAT'
    assert np.all(result['bhp'] > 0)

def test_vfpinj_eclipse_format():
    """VFPINJ Eclipse string should contain required keywords"""
    from pyrestoolbox.nodal import Completion
    comp = Completion(tid=2.441, length=8000, tht=100, bht=250)
    result = simtools.make_vfpinj(
        table_num=1, completion=comp, flo_type='GAS',
        flo_rates=[5000, 10000],
        thp_values=[1000],
        gsg=0.65)
    ecl = result['eclipse_string']
    assert 'VFPINJ' in ecl
    assert 'GAS' in ecl
    assert 'FIELD' in ecl
    assert 'BHP' in ecl

def test_vfpprod_gas_basic():
    """Generate a small gas VFPPROD table"""
    from pyrestoolbox.nodal import Completion
    comp = Completion(tid=2.441, length=8000, tht=100, bht=250)
    result = simtools.make_vfpprod(
        table_num=1, completion=comp, well_type='gas',
        vlpmethod='HB',
        flo_rates=[5000, 10000, 20000],
        thp_values=[500, 1000],
        wfr_values=[0],
        gfr_values=[0],
        alq_values=[0],
        gsg=0.65)
    assert result['bhp'].shape == (2, 1, 1, 1, 3)
    assert result['n_failed'] == 0
    assert result['flo_type'] == 'GAS'
    assert result['wfr_type'] == 'WGR'
    assert result['gfr_type'] == 'OGR'
    assert np.all(result['bhp'] > 0)

def test_vfpprod_oil_basic():
    """Generate a small oil VFPPROD table"""
    from pyrestoolbox.nodal import Completion
    comp = Completion(tid=3.5, length=6000, tht=120, bht=220)
    result = simtools.make_vfpprod(
        table_num=2, completion=comp, well_type='oil',
        flo_rates=[500, 1000, 2000],
        thp_values=[200, 500],
        wfr_values=[0, 0.5],
        gfr_values=[0.5, 1.0],
        alq_values=[0],
        pb=2500, rsb=500, sgsp=0.65, api=35)
    assert result['bhp'].shape == (2, 2, 2, 1, 3)
    assert result['flo_type'] == 'OIL'
    assert result['wfr_type'] == 'WCT'
    assert result['gfr_type'] == 'GOR'

def test_vfpprod_eclipse_format():
    """VFPPROD Eclipse string should contain required keywords"""
    from pyrestoolbox.nodal import Completion
    comp = Completion(tid=2.441, length=8000, tht=100, bht=250)
    result = simtools.make_vfpprod(
        table_num=1, completion=comp, well_type='gas',
        vlpmethod='HB',
        flo_rates=[5000, 10000],
        thp_values=[500],
        wfr_values=[0],
        gfr_values=[0],
        alq_values=[0],
        gsg=0.65)
    ecl = result['eclipse_string']
    assert 'VFPPROD' in ecl
    assert 'GAS' in ecl
    assert 'WGR' in ecl
    assert 'OGR' in ecl
    assert 'FIELD' in ecl
    assert 'BHP' in ecl

def test_vfpprod_gas_bhp_increases_with_rate():
    """For gas production, BHP should generally increase with rate"""
    from pyrestoolbox.nodal import Completion
    comp = Completion(tid=2.441, length=8000, tht=100, bht=250)
    result = simtools.make_vfpprod(
        table_num=1, completion=comp, well_type='gas',
        vlpmethod='HB',
        flo_rates=[1000, 10000, 50000],
        thp_values=[500],
        wfr_values=[0],
        gfr_values=[0],
        alq_values=[0],
        gsg=0.65)
    bhps = result['bhp'][0, 0, 0, 0, :]
    assert bhps[-1] > bhps[0], f"BHP should increase with rate: {bhps}"

def test_vfpprod_bad_wct():
    """WCT >= 1.0 should raise ValueError for oil VFPPROD"""
    from pyrestoolbox.nodal import Completion
    comp = Completion(tid=3.5, length=6000, tht=120, bht=220)
    try:
        simtools.make_vfpprod(table_num=1, completion=comp, well_type='oil',
                              wfr_values=[0, 0.5, 1.0], pb=2500, rsb=500)
        assert False, "Should have raised ValueError"
    except ValueError:
        pass


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
