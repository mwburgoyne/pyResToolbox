#!/usr/bin/env python3
"""
Tests for Rust acceleration layer.

Verifies:
1. Accelerator module loads and reports status correctly
2. Rust vs Python numerical equivalence for all accelerated functions
3. Graceful fallback behavior

Run with: PYTHONPATH=/home/mark/projects python3 -m pytest tests/test_rust_acceleration.py -v
"""

import os
import sys
import contextlib
import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from pyrestoolbox._accelerator import get_status, RUST_AVAILABLE

rust_required = pytest.mark.skipif(
    not RUST_AVAILABLE,
    reason="Rust extension not available"
)

RTOL_TIGHT = 1e-6
RTOL_MEDIUM = 1e-4
RTOL_LOOSE = 1e-3


@contextlib.contextmanager
def force_python():
    """Context manager that disables Rust in all modules, then restores."""
    import pyrestoolbox._accelerator as acc
    import pyrestoolbox.gas.gas as gas_mod
    import pyrestoolbox.oil._density as oil_mod
    import pyrestoolbox.brine.brine as brine_mod
    import pyrestoolbox.nodal.nodal as nodal_mod
    import pyrestoolbox.dca.dca as dca_mod
    import pyrestoolbox.matbal.matbal as matbal_mod

    modules = [
        (acc, 'RUST_AVAILABLE'),
        (gas_mod, 'RUST_AVAILABLE'),
        (oil_mod, '_RUST_AVAILABLE'),
        (brine_mod, '_RUST_AVAILABLE'),
        (nodal_mod, '_RUST_AVAILABLE'),
        (dca_mod, '_RUST_AVAILABLE'),
        (matbal_mod, '_RUST_AVAILABLE'),
    ]

    try:
        import pyrestoolbox.brine._lib_vle_engine as vle_mod
        modules.append((vle_mod, '_RUST_AVAILABLE'))
    except ImportError:
        pass

    saved = [(mod, attr, getattr(mod, attr)) for mod, attr in modules]

    for mod, attr, _ in saved:
        setattr(mod, attr, False)
    try:
        yield
    finally:
        for mod, attr, orig_val in saved:
            setattr(mod, attr, orig_val)


# =============================================================================
# Accelerator module tests
# =============================================================================

class TestAcceleratorStatus:
    def test_get_status_returns_dict(self):
        status = get_status()
        assert isinstance(status, dict)
        for key in ('rust_available', 'failure_reason', 'forced_python'):
            assert key in status

    def test_status_rust_available_is_bool(self):
        assert isinstance(get_status()['rust_available'], bool)

    def test_forced_python_reflects_env(self):
        env_val = os.environ.get("PYRESTOOLBOX_NO_RUST", "").strip()
        expected = env_val in ("1", "true", "yes")
        assert get_status()['forced_python'] == expected


# =============================================================================
# Gas module: Z-factor equivalence
# =============================================================================

@rust_required
class TestGasZFactorEquivalence:

    def _compare_z(self, **kwargs):
        import pyrestoolbox.gas as gas
        z_rust = gas.gas_z(**kwargs)
        with force_python():
            z_python = gas.gas_z(**kwargs)
        np.testing.assert_allclose(
            np.atleast_1d(z_rust), np.atleast_1d(z_python),
            rtol=RTOL_MEDIUM, atol=1e-6,
            err_msg=f"Z-factor mismatch for {kwargs}"
        )

    def test_dak_single(self):
        self._compare_z(p=2000, sg=0.75, degf=200, zmethod='DAK', cmethod='SUT')

    def test_dak_array(self):
        self._compare_z(p=[1000, 2000, 4000, 8000], sg=0.75, degf=200,
                         zmethod='DAK', cmethod='SUT')

    def test_hy_single(self):
        self._compare_z(p=3000, sg=0.65, degf=180, zmethod='HY', cmethod='SUT')

    def test_bns_single(self):
        self._compare_z(p=2000, sg=0.75, degf=200, zmethod='BNS')

    def test_bns_with_inerts(self):
        self._compare_z(p=3000, sg=0.75, degf=200, zmethod='BNS',
                         co2=0.05, h2s=0.02, n2=0.03)

    def test_bns_array(self):
        self._compare_z(p=[500, 1000, 3000, 6000, 10000], sg=0.7, degf=250,
                         zmethod='BNS')


# =============================================================================
# Gas module: Viscosity equivalence
# =============================================================================

@rust_required
class TestGasViscosityEquivalence:

    def _compare_ug(self, **kwargs):
        import pyrestoolbox.gas as gas
        ug_rust = gas.gas_ug(**kwargs)
        with force_python():
            ug_python = gas.gas_ug(**kwargs)
        np.testing.assert_allclose(
            np.atleast_1d(ug_rust), np.atleast_1d(ug_python),
            rtol=RTOL_MEDIUM, atol=1e-6,
            err_msg=f"Viscosity mismatch for {kwargs}"
        )

    def test_lge_viscosity(self):
        self._compare_ug(p=2000, sg=0.75, degf=200)

    def test_lbc_viscosity_bns(self):
        self._compare_ug(p=2000, sg=0.75, degf=200, zmethod='BNS')

    def test_viscosity_array(self):
        self._compare_ug(p=[1000, 3000, 5000], sg=0.7, degf=250)


# =============================================================================
# Gas module: Pseudopressure equivalence
# =============================================================================

@rust_required
class TestPseudopressureEquivalence:

    def _compare_dmp(self, **kwargs):
        import pyrestoolbox.gas as gas
        dmp_rust = gas.gas_dmp(**kwargs)
        with force_python():
            dmp_python = gas.gas_dmp(**kwargs)
        np.testing.assert_allclose(
            np.atleast_1d(dmp_rust), np.atleast_1d(dmp_python),
            rtol=RTOL_LOOSE, atol=10.0,
            err_msg=f"Pseudopressure mismatch for {kwargs}"
        )

    def test_dmp_single(self):
        self._compare_dmp(p1=14.7, p2=2000, sg=0.75, degf=200)

    def test_dmp_high_pressure(self):
        self._compare_dmp(p1=14.7, p2=5000, sg=0.7, degf=250)


# =============================================================================
# Oil module: Density equivalence
# =============================================================================

@rust_required
class TestOilEquivalence:

    def test_oil_density_below_pb(self):
        import pyrestoolbox.oil as oil
        kwargs = dict(p=2000, degf=200, rs=500, rsb=600, api=35,
                      sg_sp=0.8, sg_g=0.75, pb=2500, denomethod='SWMH')
        den_rust = oil.oil_deno(**kwargs)
        with force_python():
            den_python = oil.oil_deno(**kwargs)
        np.testing.assert_allclose(den_rust, den_python, rtol=RTOL_MEDIUM,
                                   err_msg="Oil density mismatch below Pb")

    def test_oil_density_above_pb(self):
        import pyrestoolbox.oil as oil
        kwargs = dict(p=4000, degf=200, rs=500, rsb=500, api=35,
                      sg_sp=0.8, sg_g=0.75, pb=2500, denomethod='SWMH')
        den_rust = oil.oil_deno(**kwargs)
        with force_python():
            den_python = oil.oil_deno(**kwargs)
        np.testing.assert_allclose(den_rust, den_python, rtol=RTOL_MEDIUM,
                                   err_msg="Oil density mismatch above Pb")


# =============================================================================
# Nodal module: VLP segment loop equivalence
# =============================================================================

@rust_required
class TestNodalVLPEquivalence:

    def _gas_completion(self):
        from pyrestoolbox.nodal import nodal
        from pyrestoolbox.gas.gas import GasPVT
        comp = nodal.Completion(tid=2.441, length=10000, tht=100, bht=200)
        pvt = GasPVT(sg=0.75, co2=0, h2s=0, n2=0)
        return pvt, comp

    def _oil_completion(self):
        from pyrestoolbox.nodal import nodal
        from pyrestoolbox.oil import OilPVT
        comp = nodal.Completion(tid=2.441, length=10000, tht=100, bht=200)
        pvt = OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
        return pvt, comp

    def _compare_fbhp_gas(self, vlpmethod, rtol=RTOL_LOOSE):
        from pyrestoolbox.nodal import nodal
        pvt, comp = self._gas_completion()
        kwargs = dict(thp=500, gas_pvt=pvt, completion=comp,
                      vlpmethod=vlpmethod, qg_mmscfd=5)
        result_rust = nodal.fbhp(**kwargs)
        with force_python():
            result_python = nodal.fbhp(**kwargs)
        np.testing.assert_allclose(
            result_rust, result_python,
            rtol=rtol, atol=1.0,
            err_msg=f"VLP {vlpmethod} gas mismatch"
        )

    def _compare_fbhp_oil(self, vlpmethod, rtol=RTOL_LOOSE):
        from pyrestoolbox.nodal import nodal
        pvt, comp = self._oil_completion()
        kwargs = dict(thp=200, oil_pvt=pvt, completion=comp,
                      vlpmethod=vlpmethod, qt_stbpd=1000, gor=500, wc=0.1)
        result_rust = nodal.fbhp(**kwargs)
        with force_python():
            result_python = nodal.fbhp(**kwargs)
        np.testing.assert_allclose(
            result_rust, result_python,
            rtol=rtol, atol=1.0,
            err_msg=f"VLP {vlpmethod} oil mismatch"
        )

    def test_hb_gas(self):
        self._compare_fbhp_gas('HB')

    def test_hb_oil(self):
        self._compare_fbhp_oil('HB')

    def test_wg_gas(self):
        self._compare_fbhp_gas('WG')

    def test_wg_oil(self):
        self._compare_fbhp_oil('WG')

    def test_gray_gas(self):
        self._compare_fbhp_gas('GRAY')

    def test_gray_oil(self):
        self._compare_fbhp_oil('GRAY')

    def test_bb_gas(self):
        self._compare_fbhp_gas('BB')

    def test_bb_oil(self):
        self._compare_fbhp_oil('BB')


# =============================================================================
# DCA module: Hyperbolic fitting equivalence
# =============================================================================

@rust_required
class TestDCAEquivalence:

    def test_fit_decline(self):
        import pyrestoolbox.dca as dca
        t = np.arange(1, 37, dtype=float)
        q = 1000 * (1 + 0.8 * 0.05 * t) ** (-1 / 0.8) + np.random.RandomState(42).normal(0, 10, len(t))
        result_rust = dca.fit_decline(t, q, method='hyperbolic')
        with force_python():
            result_python = dca.fit_decline(t, q, method='hyperbolic')
        np.testing.assert_allclose(result_rust.qi, result_python.qi, rtol=RTOL_LOOSE,
                                   err_msg="DCA qi mismatch")
        np.testing.assert_allclose(result_rust.di, result_python.di, rtol=RTOL_LOOSE,
                                   err_msg="DCA di mismatch")
        np.testing.assert_allclose(result_rust.b, result_python.b, atol=0.05,
                                   err_msg="DCA b mismatch")


# =============================================================================
# Brine module: CO2 solubility equivalence
# =============================================================================

@rust_required
class TestBrineCO2Equivalence:

    def test_co2_solubility(self):
        from pyrestoolbox.brine import brine

        # Rust-accelerated construction
        mix_rust = brine.CO2_Brine_Mixture(pres=5000, temp=275, ppm=30000, metric=False)

        with force_python():
            mix_python = brine.CO2_Brine_Mixture(pres=5000, temp=275, ppm=30000, metric=False)

        # Compare xCO2 and Rs
        np.testing.assert_allclose(mix_rust.x[0], mix_python.x[0], rtol=RTOL_TIGHT,
                                   err_msg="CO2 xCO2 mismatch")
        np.testing.assert_allclose(mix_rust.Rs, mix_python.Rs, rtol=RTOL_TIGHT,
                                   err_msg="CO2 Rs mismatch")


# =============================================================================
# Material balance: Objective function equivalence
# =============================================================================

@rust_required
class TestMatbalEquivalence:

    def test_oil_matbal(self):
        from pyrestoolbox.matbal import matbal
        p = [5000, 4500, 4000, 3500, 3000]
        Np = [0, 50000, 120000, 210000, 320000]
        kwargs = dict(p=p, Np=Np, api=35, sg_sp=0.8, sg_g=0.75,
                      degf=200, pb=3500, sw_i=0.25, cf=5e-6, cw=3e-6)
        result_rust = matbal.oil_matbal(**kwargs)
        with force_python():
            result_python = matbal.oil_matbal(**kwargs)
        np.testing.assert_allclose(result_rust.ooip, result_python.ooip, rtol=RTOL_MEDIUM,
                                   err_msg="Matbal OOIP mismatch")


# =============================================================================

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
