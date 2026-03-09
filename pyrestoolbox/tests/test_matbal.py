#!/usr/bin/env python3
"""Tests for the material balance module."""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
from pyrestoolbox import matbal


# ======================== gas_matbal tests ========================

def test_gas_matbal_basic():
    """Basic P/Z material balance with synthetic data."""
    # Create synthetic data: OGIP = 100 Bscf, pi = 5000 psia, T = 200 degF
    # P/Z linear decline: p/z = pi/zi * (1 - Gp/OGIP)
    from pyrestoolbox import gas
    pi = 5000
    degf = 200
    sg = 0.65
    zi = gas.gas_z(pi, degf, sg)
    ogip_true = 100.0  # Bscf

    # Generate data points along the P/Z line
    Gp = np.array([0, 10, 20, 30, 40, 50])
    pz_true = (pi / zi) * (1.0 - Gp / ogip_true)
    # Recover pressures: p = pz * z, but z depends on p, so iterate
    # Instead, use gas_matbal which handles this internally
    # For test, use p values that approximately match the P/Z line
    p = np.array([5000, 4500, 4000, 3500, 3000, 2500])

    result = matbal.gas_matbal(p, Gp, degf, sg=sg)
    assert isinstance(result, matbal.GasMatbalResult)
    assert result.ogip > 0
    assert result.r_squared > 0.9
    assert len(result.pz) == 6
    assert result.p_initial == 5000


def test_gas_matbal_perfect_linear():
    """Perfect linear P/Z data should give very high R-squared."""
    # Construct pressures where P/Z is exactly linear
    # Use a simple case: constant Z approximation
    from pyrestoolbox import gas
    degf = 200
    sg = 0.65
    # Use high-pressure range where Z varies less
    p = np.array([6000, 5500, 5000, 4500, 4000])
    Gp = np.array([0, 5, 10, 15, 20])

    result = matbal.gas_matbal(p, Gp, degf, sg=sg)
    assert result.ogip > 0
    assert result.r_squared > 0.95


def test_gas_matbal_z_initial():
    """Z-factor at initial pressure should be physically reasonable."""
    result = matbal.gas_matbal(
        [5000, 4000, 3000], [0, 10, 20], degf=200, sg=0.65
    )
    assert 0.5 < result.z_initial < 1.5


def test_gas_matbal_slope_negative():
    """Slope of P/Z vs Gp should be negative (declining P/Z)."""
    result = matbal.gas_matbal(
        [5000, 4000, 3000], [0, 10, 20], degf=200, sg=0.65
    )
    assert result.slope < 0


def test_gas_matbal_too_few_points():
    try:
        matbal.gas_matbal([5000], [0], degf=200)
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_gas_matbal_length_mismatch():
    try:
        matbal.gas_matbal([5000, 4000], [0], degf=200)
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_gas_matbal_bns_method():
    """Gas matbal with BNS method should work."""
    result = matbal.gas_matbal(
        [5000, 4000, 3000], [0, 10, 20], degf=200, sg=0.65,
        zmethod='BNS', cmethod='BNS'
    )
    assert result.ogip > 0


def test_gas_matbal_with_impurities():
    """Gas matbal with CO2 and N2."""
    result = matbal.gas_matbal(
        [5000, 4000, 3000], [0, 10, 20], degf=200, sg=0.70,
        co2=0.05, n2=0.03
    )
    assert result.ogip > 0


# ======================== gas_matbal Cole plot & aquifer tests ========================

def test_gas_matbal_cole_volumetric():
    """For volumetric reservoir, F/Et should be approximately constant (= OGIP)."""
    p = [3000, 2700, 2400, 2100, 1800]
    Gp = [0, 5, 12, 22, 35]
    result = matbal.gas_matbal(p, Gp, degf=200, sg=0.65)
    # F/Et values at indices 1: should be roughly constant
    cole = result.cole_F_over_Et[1:]
    assert all(np.isfinite(cole))
    # Relative spread should be small for volumetric data
    cole_mean = np.mean(cole)
    assert cole_mean > 0
    # Each value should be within 30% of the mean (P/Z gives non-perfect OGIP)
    for v in cole:
        assert abs(v - cole_mean) / cole_mean < 0.30


def test_gas_matbal_bg_populated():
    """Bg array should be correct length, all positive, increasing as p decreases."""
    result = matbal.gas_matbal(
        [5000, 4000, 3000, 2000], [0, 10, 20, 30], degf=200, sg=0.65
    )
    assert len(result.bg) == 4
    assert all(result.bg > 0)
    # Bg increases as pressure decreases
    for i in range(1, len(result.bg)):
        assert result.bg[i] > result.bg[i-1]


def test_gas_matbal_Et_initial_zero():
    """Et at initial pressure should be zero (Bg - Bgi = 0)."""
    result = matbal.gas_matbal(
        [5000, 4000, 3000], [0, 10, 20], degf=200, sg=0.65
    )
    assert result.Et[0] == 0.0


def test_gas_matbal_cole_initial_nan():
    """Cole F/Et at index 0 should be NaN."""
    result = matbal.gas_matbal(
        [5000, 4000, 3000], [0, 10, 20], degf=200, sg=0.65
    )
    assert np.isnan(result.cole_F_over_Et[0])


def test_gas_matbal_method_default_pz():
    """Default method should be 'pz' when We not provided."""
    result = matbal.gas_matbal(
        [5000, 4000, 3000], [0, 10, 20], degf=200, sg=0.65
    )
    assert result.method == 'pz'


def test_gas_matbal_with_wp():
    """F should include water term when Wp provided, making F > Gp*Bg."""
    p = [3000, 2700, 2400, 2100]
    Gp = [0, 5e6, 12e6, 22e6]  # scf
    Wp = [0, 100, 300, 600]     # STB
    result = matbal.gas_matbal(p, Gp, degf=200, sg=0.65, Wp=Wp, Bw=1.02)
    result_no_wp = matbal.gas_matbal(p, Gp, degf=200, sg=0.65)
    # F with Wp should be larger at non-zero Wp points
    for i in range(1, len(p)):
        assert result.F[i] > result_no_wp.F[i]


def test_gas_matbal_wp_length_mismatch():
    """Wp with wrong length should raise ValueError."""
    try:
        matbal.gas_matbal(
            [5000, 4000, 3000], [0, 10, 20], degf=200, Wp=[0, 100]
        )
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_gas_matbal_havlena_odeh():
    """With We provided, should use Havlena-Odeh and produce reasonable OGIP."""
    p = [3000, 2700, 2400, 2100, 1800]
    Gp = [0, 5e9, 12e9, 22e9, 35e9]  # scf
    # Synthetic We: small aquifer influx (rcf)
    We = [0, 1e6, 3e6, 6e6, 10e6]
    result = matbal.gas_matbal(p, Gp, degf=200, sg=0.65, We=We)
    assert result.method == 'havlena_odeh'
    assert result.ogip > 0


def test_gas_matbal_we_zero_matches_pz():
    """We all zeros should give HO ogip close to P/Z ogip."""
    p = [3000, 2700, 2400, 2100, 1800]
    Gp = [0, 5, 12, 22, 35]
    result_pz = matbal.gas_matbal(p, Gp, degf=200, sg=0.65)
    We_zero = [0, 0, 0, 0, 0]
    result_ho = matbal.gas_matbal(p, Gp, degf=200, sg=0.65, We=We_zero)
    assert result_ho.method == 'havlena_odeh'
    # Both OGIPs should be reasonably close
    assert abs(result_ho.ogip - result_pz.ogip) / abs(result_pz.ogip) < 0.15


def test_gas_matbal_we_length_mismatch():
    """We with wrong length should raise ValueError."""
    try:
        matbal.gas_matbal(
            [5000, 4000, 3000], [0, 10, 20], degf=200, We=[0, 1e6]
        )
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_gas_matbal_ho_pz_still_computed():
    """P/Z regression should still be populated when We is given."""
    p = [3000, 2700, 2400, 2100]
    Gp = [0, 5, 12, 22]
    We = [0, 1e5, 3e5, 6e5]
    result = matbal.gas_matbal(p, Gp, degf=200, sg=0.65, We=We)
    assert result.slope < 0
    assert result.intercept > 0
    assert result.r_squared > 0


def test_gas_matbal_backward_compat():
    """Existing call pattern should produce identical ogip/slope/intercept."""
    p = [3000, 2700, 2400, 2100, 1800]
    Gp = [0, 5, 12, 22, 35]
    result = matbal.gas_matbal(p, Gp, degf=200, sg=0.65)
    # These should match the documented values from matbal.rst
    assert abs(result.ogip - 87.602774253829) / 87.602774253829 < 1e-6
    assert abs(result.z_initial - 0.9163208839373836) / 0.9163208839373836 < 1e-6
    assert abs(result.r_squared - 0.9734794008096929) / 0.9734794008096929 < 1e-6


# ======================== oil_matbal tests ========================

def test_oil_matbal_basic():
    """Basic Havlena-Odeh oil material balance."""
    result = matbal.oil_matbal(
        p=[4000, 3500, 3000, 2500, 2000],
        Np=[0, 500000, 1200000, 2100000, 3200000],
        degf=200,
        api=35,
        sg_sp=0.65,
        pb=3000,
        rsb=500,
    )
    assert isinstance(result, matbal.OilMatbalResult)
    assert result.ooip > 0
    assert len(result.F) == 5
    assert len(result.Eo) == 5
    assert 'DDI' in result.drive_indices
    assert 'SDI' in result.drive_indices
    assert 'CDI' in result.drive_indices


def test_oil_matbal_above_pb():
    """All pressures above Pb: undersaturated depletion."""
    result = matbal.oil_matbal(
        p=[5000, 4800, 4600, 4400],
        Np=[0, 100000, 200000, 300000],
        degf=200,
        api=35,
        sg_sp=0.65,
        pb=3000,
        rsb=500,
        cf=3e-6,
        sw_i=0.2,
        cw=3e-6,
    )
    assert result.ooip > 0


def test_oil_matbal_with_gas_cap():
    """Oil matbal with gas cap (m > 0)."""
    result = matbal.oil_matbal(
        p=[4000, 3500, 3000, 2500],
        Np=[0, 300000, 700000, 1200000],
        degf=200,
        api=35,
        sg_sp=0.65,
        pb=3500,
        rsb=600,
        m=0.5,
    )
    assert result.ooip > 0


def test_oil_matbal_pvt_output():
    """Check that PVT outputs are populated."""
    result = matbal.oil_matbal(
        p=[4000, 3500, 3000],
        Np=[0, 500000, 1200000],
        degf=200,
        api=35,
        sg_sp=0.65,
        pb=3500,
        rsb=500,
    )
    assert 'Rs' in result.pvt
    assert 'Bo' in result.pvt
    assert 'Bg' in result.pvt
    assert len(result.pvt['Rs']) == 3
    # Bo should be > 1
    assert all(bo > 1.0 for bo in result.pvt['Bo'])


def test_oil_matbal_rsb_from_pb():
    """When rsb=0, it should be calculated from pb."""
    result = matbal.oil_matbal(
        p=[4000, 3500, 3000],
        Np=[0, 500000, 1200000],
        degf=200,
        api=35,
        sg_sp=0.65,
        pb=3000,
    )
    assert result.ooip > 0


def test_oil_matbal_pb_from_rsb():
    """When pb=0, it should be calculated from rsb."""
    result = matbal.oil_matbal(
        p=[4000, 3500, 3000],
        Np=[0, 500000, 1200000],
        degf=200,
        api=35,
        sg_sp=0.65,
        rsb=500,
    )
    assert result.ooip > 0


def test_oil_matbal_too_few_points():
    try:
        matbal.oil_matbal(
            p=[4000], Np=[0], degf=200, api=35, sg_sp=0.65, pb=3000
        )
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_oil_matbal_no_pb_rsb():
    """Should raise error when neither pb nor rsb specified."""
    try:
        matbal.oil_matbal(
            p=[4000, 3500], Np=[0, 500000], degf=200, api=35, sg_sp=0.65
        )
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_oil_matbal_drive_indices_sum():
    """Drive indices should approximately sum to 1 at each valid point."""
    result = matbal.oil_matbal(
        p=[4000, 3500, 3000, 2500],
        Np=[0, 500000, 1200000, 2100000],
        degf=200,
        api=35,
        sg_sp=0.65,
        pb=3500,
        rsb=500,
        m=0.3,
        cf=3e-6,
        sw_i=0.2,
        cw=3e-6,
    )
    ddi = result.drive_indices['DDI']
    sdi = result.drive_indices['SDI']
    cdi = result.drive_indices['CDI']
    # For non-zero points (skip initial), indices should sum to ~1
    for i in range(1, len(ddi)):
        total = ddi[i] + sdi[i] + cdi[i]
        if total > 0:
            assert abs(total - 1.0) < 0.05, f"Drive indices sum to {total} at step {i}"


# ======================== oil_matbal regression tests ========================

def test_oil_matbal_regress_m():
    """Regress m (gas cap ratio). regressed['m'] > 0, ooip > 0."""
    result = matbal.oil_matbal(
        p=[4000, 3500, 3000, 2500],
        Np=[0, 300000, 700000, 1200000],
        degf=200, api=35, sg_sp=0.65,
        pb=3500, rsb=600,
        m=0,
        regress={'m': (0, 2)}
    )
    assert result.regressed is not None
    assert 'm' in result.regressed
    assert result.regressed['m'] >= 0
    assert result.ooip > 0


def test_oil_matbal_regress_cf():
    """Regress cf. regressed['cf'] > 0."""
    result = matbal.oil_matbal(
        p=[4000, 3500, 3000, 2500],
        Np=[0, 500000, 1200000, 2100000],
        degf=200, api=35, sg_sp=0.65,
        pb=3500, rsb=500,
        cf=0, sw_i=0.2, cw=3e-6,
        regress={'cf': (1e-7, 50e-6)}
    )
    assert result.regressed is not None
    assert 'cf' in result.regressed
    assert result.regressed['cf'] >= 1e-7


def test_oil_matbal_regress_multiple():
    """Regress m and cf simultaneously."""
    result = matbal.oil_matbal(
        p=[4000, 3500, 3000, 2500],
        Np=[0, 300000, 700000, 1200000],
        degf=200, api=35, sg_sp=0.65,
        pb=3500, rsb=600,
        m=0, cf=0, sw_i=0.2, cw=3e-6,
        regress={'m': (0, 2), 'cf': (1e-7, 50e-6)}
    )
    assert result.regressed is not None
    assert 'm' in result.regressed
    assert 'cf' in result.regressed


def test_oil_matbal_regress_none_default():
    """No regress. regressed is None. Frozen baseline matches."""
    result = matbal.oil_matbal(
        p=[4000, 3500, 3000, 2500],
        Np=[0, 1e6, 3e6, 6e6],
        degf=220, api=35, sg_sp=0.75,
        pb=3500, rsb=500, cf=3e-6, sw_i=0.2, cw=3e-6
    )
    assert result.regressed is None
    assert abs(result.ooip - 82793519.84914012) / 82793519.84914012 < 1e-6


def test_oil_matbal_regress_bad_param():
    """regress={'invalid': (0, 1)}. ValueError."""
    try:
        matbal.oil_matbal(
            p=[4000, 3500, 3000],
            Np=[0, 500000, 1200000],
            degf=200, api=35, sg_sp=0.65,
            pb=3500, rsb=500,
            regress={'invalid': (0, 1)}
        )
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_oil_matbal_regress_bad_bounds():
    """regress={'m': (2, 0)}. ValueError (lower >= upper)."""
    try:
        matbal.oil_matbal(
            p=[4000, 3500, 3000],
            Np=[0, 500000, 1200000],
            degf=200, api=35, sg_sp=0.65,
            pb=3500, rsb=500,
            regress={'m': (2, 0)}
        )
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_oil_matbal_regress_bad_format():
    """regress={'m': 0.5}. ValueError (not a tuple)."""
    try:
        matbal.oil_matbal(
            p=[4000, 3500, 3000],
            Np=[0, 500000, 1200000],
            degf=200, api=35, sg_sp=0.65,
            pb=3500, rsb=500,
            regress={'m': 0.5}
        )
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_oil_matbal_regress_reduces_cv():
    """CV with regressed params < CV without."""
    # Without regression
    r_noreg = matbal.oil_matbal(
        p=[4000, 3500, 3000, 2500],
        Np=[0, 500000, 1200000, 2100000],
        degf=200, api=35, sg_sp=0.65,
        pb=3500, rsb=500,
        m=0, cf=0, sw_i=0.2, cw=3e-6,
    )
    denom_noreg = r_noreg.Eo + 0 * r_noreg.Eg + 1.0 * r_noreg.Efw
    valid_noreg = np.abs(denom_noreg) > 1e-30
    N_noreg = r_noreg.F[valid_noreg] / denom_noreg[valid_noreg]
    cv_noreg = np.std(N_noreg) / abs(np.mean(N_noreg))

    # With regression
    r_reg = matbal.oil_matbal(
        p=[4000, 3500, 3000, 2500],
        Np=[0, 500000, 1200000, 2100000],
        degf=200, api=35, sg_sp=0.65,
        pb=3500, rsb=500,
        m=0, cf=0, sw_i=0.2, cw=3e-6,
        regress={'m': (0, 2), 'cf': (1e-7, 50e-6)}
    )
    m_reg = r_reg.regressed['m']
    denom_reg = r_reg.Eo + m_reg * r_reg.Eg + (1.0 + m_reg) * r_reg.Efw
    valid_reg = np.abs(denom_reg) > 1e-30
    N_reg = r_reg.F[valid_reg] / denom_reg[valid_reg]
    cv_reg = np.std(N_reg) / abs(np.mean(N_reg))

    assert cv_reg < cv_noreg


# ======================== oil_matbal pvt_table tests ========================

def test_oil_matbal_pvt_table_roundtrip():
    """Build PVT table from correlations, compare ooip."""
    from pyrestoolbox import oil as oil_mod, gas as gas_mod
    from pyrestoolbox.constants import CUFTperBBL as CUFT

    p_survey = [4000, 3500, 3000, 2500]
    degf = 220
    api = 35
    sg_sp = 0.75
    pb = 3500
    rsb = 500
    sg_o = 141.5 / (api + 131.5)

    # Include exact survey pressures for perfect match
    p_table = np.sort(np.unique(np.concatenate([
        np.linspace(2000, 5000, 200),
        np.array(p_survey, dtype=float)
    ])))
    Rs_t = [oil_mod.oil_rs(api, degf, sg_sp, pi, pb=pb, rsb=rsb) for pi in p_table]
    Bo_t = [oil_mod.oil_bo(pi, pb, degf, rs, rsb, sg_o, sg_sp=sg_sp) for pi, rs in zip(p_table, Rs_t)]
    Bg_t = [gas_mod.gas_bg(pi, sg_sp, degf) / CUFT for pi in p_table]

    r_corr = matbal.oil_matbal(
        p=p_survey, Np=[0, 1e6, 3e6, 6e6], degf=degf,
        api=api, sg_sp=sg_sp, pb=pb, rsb=rsb, cf=3e-6, sw_i=0.2, cw=3e-6
    )
    r_table = matbal.oil_matbal(
        p=p_survey, Np=[0, 1e6, 3e6, 6e6], degf=degf,
        pb=pb, rsb=rsb, cf=3e-6, sw_i=0.2, cw=3e-6,
        pvt_table={'p': list(p_table), 'Rs': Rs_t, 'Bo': Bo_t, 'Bg': Bg_t}
    )
    assert abs(r_table.ooip - r_corr.ooip) / r_corr.ooip < 0.02


def test_oil_matbal_pvt_table_no_api():
    """pvt_table + api=0, sg_sp=0 should work."""
    from pyrestoolbox import oil as oil_mod, gas as gas_mod
    from pyrestoolbox.constants import CUFTperBBL as CUFT

    api, sg_sp, pb, rsb, degf = 35, 0.75, 3500, 500, 220
    sg_o = 141.5 / (api + 131.5)
    p_table = np.linspace(2000, 5000, 200)
    Rs_t = [oil_mod.oil_rs(api, degf, sg_sp, pi, pb=pb, rsb=rsb) for pi in p_table]
    Bo_t = [oil_mod.oil_bo(pi, pb, degf, rs, rsb, sg_o, sg_sp=sg_sp) for pi, rs in zip(p_table, Rs_t)]
    Bg_t = [gas_mod.gas_bg(pi, sg_sp, degf) / CUFT for pi in p_table]

    result = matbal.oil_matbal(
        p=[4000, 3500, 3000, 2500], Np=[0, 1e6, 3e6, 6e6], degf=degf,
        pb=pb, rsb=rsb, cf=3e-6, sw_i=0.2, cw=3e-6,
        pvt_table={'p': list(p_table), 'Rs': Rs_t, 'Bo': Bo_t, 'Bg': Bg_t}
    )
    assert result.ooip > 0


def test_oil_matbal_pvt_table_no_pvt_no_api():
    """No pvt_table and api=0 should raise ValueError."""
    try:
        matbal.oil_matbal(
            p=[4000, 3500, 3000], Np=[0, 500000, 1200000], degf=200,
            pb=3500, rsb=500
        )
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_oil_matbal_pvt_table_missing_key():
    """pvt_table missing 'Bg' should raise ValueError."""
    try:
        matbal.oil_matbal(
            p=[4000, 3500, 3000], Np=[0, 500000, 1200000], degf=200,
            pb=3500, rsb=500,
            pvt_table={'p': [2000, 3000, 4000], 'Rs': [300, 400, 500], 'Bo': [1.2, 1.3, 1.4]}
        )
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_oil_matbal_pvt_table_infer_pb():
    """No pb/rsb specified — inferred from table."""
    from pyrestoolbox import oil as oil_mod, gas as gas_mod
    from pyrestoolbox.constants import CUFTperBBL as CUFT

    api, sg_sp, pb, rsb, degf = 35, 0.75, 3500, 500, 220
    sg_o = 141.5 / (api + 131.5)
    p_table = np.sort(np.unique(np.concatenate([
        np.linspace(2000, 5000, 200),
        np.array([4000, 3500, 3000, 2500], dtype=float)
    ])))
    Rs_t = [oil_mod.oil_rs(api, degf, sg_sp, pi, pb=pb, rsb=rsb) for pi in p_table]
    Bo_t = [oil_mod.oil_bo(pi, pb, degf, rs, rsb, sg_o, sg_sp=sg_sp) for pi, rs in zip(p_table, Rs_t)]
    Bg_t = [gas_mod.gas_bg(pi, sg_sp, degf) / CUFT for pi in p_table]

    result = matbal.oil_matbal(
        p=[4000, 3500, 3000, 2500], Np=[0, 1e6, 3e6, 6e6], degf=degf,
        cf=3e-6, sw_i=0.2, cw=3e-6,
        pvt_table={'p': list(p_table), 'Rs': Rs_t, 'Bo': Bo_t, 'Bg': Bg_t}
    )
    assert result.ooip > 0


def test_oil_matbal_pvt_table_pressure_range():
    """Survey pressure outside table range should raise ValueError."""
    try:
        matbal.oil_matbal(
            p=[4000, 3500, 3000, 1000], Np=[0, 1e6, 3e6, 6e6], degf=220,
            pb=3500, rsb=500,
            pvt_table={'p': [2000, 3000, 4000], 'Rs': [300, 400, 500],
                       'Bo': [1.2, 1.3, 1.4], 'Bg': [0.001, 0.0008, 0.0006]}
        )
        assert False, "Should raise ValueError"
    except ValueError:
        pass


# ======================== gas_matbal pvt_table tests ========================

def test_gas_matbal_pvt_z():
    """Z table should give ogip matching correlation within 1%."""
    from pyrestoolbox import gas as gas_mod

    p_survey = [3000, 2700, 2400, 2100, 1800]
    degf = 200
    sg = 0.65

    p_table = np.linspace(1500, 3500, 200)
    Z_table = [gas_mod.gas_z(pi, sg, degf) for pi in p_table]

    r_corr = matbal.gas_matbal(p=p_survey, Gp=[0, 5, 12, 22, 35], degf=degf, sg=sg)
    r_z = matbal.gas_matbal(p=p_survey, Gp=[0, 5, 12, 22, 35], degf=degf, sg=sg,
                             pvt_table={'p': list(p_table), 'Z': Z_table})
    assert abs(r_z.ogip - r_corr.ogip) / abs(r_corr.ogip) < 0.01


def test_gas_matbal_pvt_bg():
    """Bg table should give ogip matching correlation."""
    from pyrestoolbox import gas as gas_mod

    p_survey = [3000, 2700, 2400, 2100, 1800]
    degf = 200
    sg = 0.65

    p_table = np.linspace(1500, 3500, 200)
    Bg_table = [gas_mod.gas_bg(pi, sg, degf) for pi in p_table]

    r_corr = matbal.gas_matbal(p=p_survey, Gp=[0, 5, 12, 22, 35], degf=degf, sg=sg)
    r_bg = matbal.gas_matbal(p=p_survey, Gp=[0, 5, 12, 22, 35], degf=degf, sg=sg,
                              pvt_table={'p': list(p_table), 'Bg': Bg_table})
    assert abs(r_bg.ogip - r_corr.ogip) / abs(r_corr.ogip) < 0.01


def test_gas_matbal_pvt_both_z_bg():
    """Both Z and Bg keys should raise ValueError."""
    try:
        matbal.gas_matbal(
            p=[3000, 2700, 2400], Gp=[0, 5, 12], degf=200,
            pvt_table={'p': [2000, 2500, 3000], 'Z': [0.9, 0.85, 0.8], 'Bg': [0.01, 0.008, 0.006]}
        )
        assert False, "Should raise ValueError"
    except ValueError:
        pass


def test_gas_matbal_pvt_backward_compat():
    """No pvt_table should give frozen baseline result."""
    result = matbal.gas_matbal(
        p=[3000, 2700, 2400, 2100, 1800],
        Gp=[0, 5, 12, 22, 35],
        degf=200, sg=0.65
    )
    assert abs(result.ogip - 87.602774253829) / 87.602774253829 < 1e-6
