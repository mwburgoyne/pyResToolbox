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
