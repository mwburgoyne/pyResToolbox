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
    """Context manager that disables Rust in all modules, then restores.

    Iterates the central registry in _accelerator so new module-level
    Rust-flag snapshots are picked up automatically.
    """
    import importlib
    from pyrestoolbox._accelerator import RUST_FLAG_REGISTRY

    modules = []
    for module_path, attr in RUST_FLAG_REGISTRY:
        try:
            mod = importlib.import_module(module_path)
        except ImportError:
            continue
        modules.append((mod, attr))

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
# rust_accelerated decorator behaviour (pure Python, no Rust needed)
# =============================================================================

class _CountingFakeNative:
    """Fake _native module whose attribute misses are counted."""

    def __init__(self):
        object.__setattr__(self, 'lookup_count', 0)

    def __getattr__(self, name):
        self.lookup_count += 1  # lookup_count is in __dict__, no recursion
        raise AttributeError(name)


@contextlib.contextmanager
def _fake_rust(fake_module):
    """Temporarily install a fake Rust module into the accelerator."""
    import pyrestoolbox._accelerator as acc
    saved_avail, saved_mod = acc.RUST_AVAILABLE, acc._rust_module
    acc.RUST_AVAILABLE = True
    acc._rust_module = fake_module
    try:
        yield acc
    finally:
        acc.RUST_AVAILABLE = saved_avail
        acc._rust_module = saved_mod


class TestRustAcceleratedDecorator:

    def test_missing_rust_fn_falls_back_and_caches(self):
        """Nonexistent Rust fn: Python fallback works, lookup happens once."""
        fake = _CountingFakeNative()
        with _fake_rust(fake) as acc:
            @acc.rust_accelerated('does_not_exist_rust')
            def doubler(x):
                return x * 2

            assert doubler(3) == 6
            assert doubler(4) == 8
            assert fake.lookup_count == 1, "missing fn lookup not cached"

    def test_exception_inside_rust_fn_propagates(self):
        """Errors raised during Rust execution must not be masked by fallback."""
        class FakeNative:
            @staticmethod
            def boom_rust(x):
                raise RuntimeError("rust blew up")

        with _fake_rust(FakeNative()) as acc:
            @acc.rust_accelerated('boom_rust')
            def quiet(x):
                return x

            with pytest.raises(RuntimeError, match="rust blew up"):
                quiet(1)

    def test_attributeerror_inside_rust_fn_propagates(self):
        """Regression: AttributeError raised inside the Rust call used to be
        swallowed and silently rerouted to the Python fallback."""
        class FakeNative:
            @staticmethod
            def attr_err_rust(x):
                raise AttributeError("inner attribute failure")

        with _fake_rust(FakeNative()) as acc:
            @acc.rust_accelerated('attr_err_rust')
            def quiet(x):
                return x

            with pytest.raises(AttributeError, match="inner attribute failure"):
                quiet(1)

    def test_existing_rust_fn_is_used(self):
        class FakeNative:
            @staticmethod
            def triple_rust(x):
                return x * 3

        with _fake_rust(FakeNative()) as acc:
            @acc.rust_accelerated('triple_rust')
            def triple(x):
                return -1  # Python fallback should not run

            assert triple(2) == 6


# =============================================================================
# Rust-flag registry completeness
# =============================================================================

class TestRustFlagRegistry:

    def test_registry_covers_all_module_snapshots(self):
        """Every module-level RUST_AVAILABLE/_RUST_AVAILABLE snapshot in the
        package source must appear in RUST_FLAG_REGISTRY."""
        import re
        from pathlib import Path
        import pyrestoolbox
        from pyrestoolbox._accelerator import RUST_FLAG_REGISTRY

        registered = set(RUST_FLAG_REGISTRY)
        pkg_root = Path(pyrestoolbox.__file__).parent

        found = set()
        import_re = re.compile(
            r'^from pyrestoolbox\._accelerator import (.+)$', re.M)
        define_re = re.compile(r'^(_?RUST_AVAILABLE)\s*[:=]', re.M)

        for py in pkg_root.rglob('*.py'):
            rel = py.relative_to(pkg_root)
            if rel.parts[0] == 'tests':
                continue
            module_path = 'pyrestoolbox.' + '.'.join(rel.with_suffix('').parts)
            text = py.read_text(encoding='utf-8', errors='replace')
            for match in import_re.finditer(text):
                for item in match.group(1).split(','):
                    item = item.strip()
                    if not item.startswith('RUST_AVAILABLE'):
                        continue
                    parts = item.split(' as ')
                    bound = parts[1].strip() if len(parts) == 2 else parts[0].strip()
                    found.add((module_path, bound))
            for match in define_re.finditer(text):
                found.add((module_path, match.group(1)))

        missing = found - registered
        assert not missing, (
            f"Rust-flag snapshots missing from RUST_FLAG_REGISTRY: {sorted(missing)}"
        )

    def test_registry_entries_resolve(self):
        """Each registry entry must import and expose the named attribute."""
        import importlib
        from pyrestoolbox._accelerator import RUST_FLAG_REGISTRY
        for module_path, attr in RUST_FLAG_REGISTRY:
            mod = importlib.import_module(module_path)
            assert isinstance(getattr(mod, attr), bool), (
                f"{module_path}.{attr} is not a bool flag"
            )


# =============================================================================
# Version metadata
# =============================================================================

class TestVersionMetadata:

    def test_package_version_resolves(self):
        import re
        import pyrestoolbox
        v = pyrestoolbox.__version__
        assert isinstance(v, str)
        assert re.match(r'\d+\.\d+', v), f"unexpected version string: {v!r}"

    def test_unknown_attribute_raises(self):
        import pyrestoolbox
        with pytest.raises(AttributeError):
            pyrestoolbox.no_such_attribute_xyz

    @rust_required
    def test_rust_version_reported(self):
        status = get_status()
        assert 'rust_version' in status


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
        self._compare_z(p=2000, sg=0.75, degf=200, zmethod='BNS', cmethod='BNS')

    def test_bns_with_inerts(self):
        self._compare_z(p=3000, sg=0.75, degf=200, zmethod='BNS', cmethod='BNS',
                         co2=0.05, h2s=0.02, n2=0.03)

    def test_bns_array(self):
        self._compare_z(p=[500, 1000, 3000, 6000, 10000], sg=0.7, degf=250,
                         zmethod='BNS', cmethod='BNS')

    def test_bns_user_tc_pc_hc_only(self):
        # WS-8.5: Rust BNS path honors user HC Tc/Pc override (not mixture)
        self._compare_z(p=[1000, 3000, 6000], sg=0.75, degf=200,
                         zmethod='BNS', cmethod='BNS',
                         co2=0.1, h2s=0.05, tc=343.0, pc=667.8)

    def test_dak_user_tc_pc(self):
        # WS-8.5: Rust DAK+SUT batch path honors user mixture Tc/Pc
        self._compare_z(p=[1000, 3000, 6000], sg=0.75, degf=200,
                         zmethod='DAK', cmethod='SUT', tc=380.0, pc=670.0)


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
        self._compare_ug(p=2000, sg=0.75, degf=200, zmethod='BNS', cmethod='BNS')

    def test_lbc_viscosity_bns_user_tc_pc(self):
        # WS-8.5: Rust LBC path honors user HC Tc/Pc override
        self._compare_ug(p=[1000, 3000, 6000], sg=0.75, degf=200,
                          zmethod='BNS', cmethod='BNS',
                          co2=0.1, h2s=0.05, tc=343.0, pc=667.8)

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

    def _compare_ponz2p(self, **kwargs):
        import pyrestoolbox.gas as gas
        p_rust = gas.gas_ponz2p(**kwargs)
        with force_python():
            p_python = gas.gas_ponz2p(**kwargs)
        np.testing.assert_allclose(
            np.atleast_1d(p_rust), np.atleast_1d(p_python),
            rtol=RTOL_MEDIUM, atol=1.0,
            err_msg=f"gas_ponz2p mismatch for {kwargs}"
        )

    def test_ponz2p_single(self):
        self._compare_ponz2p(poverz=3000, sg=0.75, degf=200)

    def test_ponz2p_array(self):
        self._compare_ponz2p(poverz=[1000, 3000, 5000], sg=0.7, degf=180)

    def test_ponz2p_with_inerts(self):
        self._compare_ponz2p(poverz=3500, sg=0.8, degf=220, co2=0.1, n2=0.05)


@rust_required
class TestInfluenceTableEquivalence:
    """Van Everdingen-Hurst influence tables: Rust ≈ Python."""

    def test_influence_tables(self):
        import pyrestoolbox.simtools as simtools
        # Small grid keeps the (slow) Python ilt path fast.
        kwargs = dict(ReDs=[2.0, 5.0], min_td=0.1, max_td=50, n_incr=5, M=7)
        tD_r, pDs_r = simtools.influence_tables(**kwargs)
        with force_python():
            tD_p, pDs_p = simtools.influence_tables(**kwargs)
        np.testing.assert_allclose(np.asarray(tD_r), np.asarray(tD_p),
                                   rtol=1e-10, err_msg="influence_tables tD drift")
        np.testing.assert_allclose(np.asarray(pDs_r), np.asarray(pDs_p),
                                   rtol=RTOL_MEDIUM, atol=1e-4,
                                   err_msg="influence_tables pD drift Rust vs Python")


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

    def _compare_fbhp_gas(self, vlpmethod, rtol=RTOL_MEDIUM):
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

    def _compare_fbhp_oil(self, vlpmethod, rtol=RTOL_MEDIUM):
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
# Parity harness: broader nodal grid (VLP × THP × rate × geometry)
# =============================================================================
# Prophylactic — catches future one-sided edits between Python and Rust VLP
# implementations. Tolerance is tight; any drift above 1e-4 relative is a
# genuine divergence, not numerical noise. Grid spans representative gas/oil
# operating conditions including both vertical and deviated wellbores.


def _nodal_parity_check_gas(vlpmethod, thp, qg, completion, rtol):
    """Run both Rust and Python fbhp paths and assert agreement."""
    from pyrestoolbox.nodal import nodal
    kwargs = dict(
        thp=thp, completion=completion, vlpmethod=vlpmethod,
        well_type='gas', qg_mmscfd=qg, gsg=0.75, cgr=0, qw_bwpd=0,
        api=45, oil_vis=1.0,
    )
    result_rust = nodal.fbhp(**kwargs)
    with force_python():
        result_python = nodal.fbhp(**kwargs)
    np.testing.assert_allclose(
        result_rust, result_python,
        rtol=rtol, atol=1.0,
        err_msg=f"Nodal {vlpmethod} gas parity drift at thp={thp}, qg={qg}",
    )


def _nodal_parity_check_oil(vlpmethod, thp, qt, completion, rtol):
    from pyrestoolbox.nodal import nodal
    from pyrestoolbox.oil import OilPVT
    opvt = OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    kwargs = dict(
        thp=thp, completion=completion, vlpmethod=vlpmethod,
        well_type='oil', qt_stbpd=qt, oil_pvt=opvt,
        gor=500, wc=0.2, gsg=0.65,
    )
    result_rust = nodal.fbhp(**kwargs)
    with force_python():
        result_python = nodal.fbhp(**kwargs)
    np.testing.assert_allclose(
        result_rust, result_python,
        rtol=rtol, atol=1.0,
        err_msg=f"Nodal {vlpmethod} oil parity drift at thp={thp}, qt={qt}",
    )


def _build_vertical_completion():
    from pyrestoolbox.nodal import nodal
    return nodal.Completion(tid=2.441, length=8000, tht=100, bht=200)


def _build_deviated_completion():
    from pyrestoolbox.nodal import nodal
    segs = [
        nodal.WellSegment(md=4000, id=2.441, deviation=0),
        nodal.WellSegment(md=4000, id=2.441, deviation=45),
    ]
    return nodal.Completion(segments=segs, tht=100, bht=200)


@rust_required
@pytest.mark.parametrize('vlpmethod', ['HB', 'WG', 'GRAY', 'BB'])
@pytest.mark.parametrize('thp,qg', [
    (300, 1.0), (300, 5.0),
    (800, 1.0), (800, 5.0), (800, 15.0),
    (1500, 10.0),
])
def test_nodal_gas_parity_vertical(vlpmethod, thp, qg):
    """Gas-well VLP: Rust ≈ Python across methods × (THP, rate) grid, vertical."""
    comp = _build_vertical_completion()
    _nodal_parity_check_gas(vlpmethod, thp, qg, comp, rtol=RTOL_MEDIUM)


@rust_required
@pytest.mark.parametrize('vlpmethod', ['WG', 'BB'])  # HB and GRAY are vertical-only
@pytest.mark.parametrize('thp,qg', [(500, 3.0), (1200, 8.0)])
def test_nodal_gas_parity_deviated(vlpmethod, thp, qg):
    """Gas-well VLP: Rust ≈ Python on deviated wells (WG/BB only — HB/GRAY are vertical)."""
    comp = _build_deviated_completion()
    _nodal_parity_check_gas(vlpmethod, thp, qg, comp, rtol=RTOL_MEDIUM)


@rust_required
@pytest.mark.parametrize('vlpmethod', ['HB', 'WG', 'GRAY', 'BB'])
@pytest.mark.parametrize('thp,qt', [
    (100, 500), (100, 2000),
    (300, 500), (300, 2000), (300, 5000),
])
def test_nodal_oil_parity_vertical(vlpmethod, thp, qt):
    """Oil-well VLP: Rust ≈ Python across methods × (THP, rate) grid, vertical."""
    comp = _build_vertical_completion()
    _nodal_parity_check_oil(vlpmethod, thp, qt, comp, rtol=RTOL_MEDIUM)


@rust_required
@pytest.mark.parametrize('vlpmethod', ['WG', 'BB'])
@pytest.mark.parametrize('thp,qt', [(200, 1000), (500, 3000)])
def test_nodal_oil_parity_deviated(vlpmethod, thp, qt):
    """Oil-well VLP: Rust ≈ Python on deviated wells (WG/BB only)."""
    comp = _build_deviated_completion()
    _nodal_parity_check_oil(vlpmethod, thp, qt, comp, rtol=RTOL_MEDIUM)


# =============================================================================
# Parity harness: broader brine SoreideWhitson grid
# =============================================================================


def _sw_parity_check(pres, temp, ppm, y_CO2, y_H2S, y_N2, sg, rtol):
    """Compare SoreideWhitson Python and Rust paths on a given (P,T,comp) point."""
    from pyrestoolbox.brine import brine
    kwargs = dict(pres=pres, temp=temp, ppm=ppm, y_CO2=y_CO2, y_H2S=y_H2S,
                  y_N2=y_N2, sg=sg, metric=False)
    mix_rust = brine.SoreideWhitson(**kwargs)
    with force_python():
        mix_python = brine.SoreideWhitson(**kwargs)
    # Compare dissolved-gas mole fractions (primary flash output)
    for comp in mix_rust.x:
        np.testing.assert_allclose(
            mix_rust.x[comp], mix_python.x[comp], rtol=rtol, atol=1e-10,
            err_msg=f"SW x[{comp}] drift at P={pres}, T={temp}, ppm={ppm}",
        )
    # Compare water content and total Rs
    np.testing.assert_allclose(
        mix_rust.y_H2O, mix_python.y_H2O, rtol=rtol, atol=1e-10,
        err_msg=f"SW y_H2O drift at P={pres}, T={temp}",
    )
    np.testing.assert_allclose(
        mix_rust.Rs_total, mix_python.Rs_total, rtol=rtol, atol=1e-6,
        err_msg=f"SW Rs_total drift at P={pres}, T={temp}",
    )


@rust_required
@pytest.mark.parametrize('pres,temp', [
    (1000, 150),
    (3000, 200),
    (5000, 275),
    (8000, 300),
])
@pytest.mark.parametrize('ppm', [0, 30000, 100000])
def test_sw_parity_pure_co2(pres, temp, ppm):
    """Pure-CO2 brine flash: Rust ≈ Python across P/T/salinity grid."""
    _sw_parity_check(pres, temp, ppm, y_CO2=1.0, y_H2S=0, y_N2=0, sg=0.65,
                     rtol=RTOL_MEDIUM)


@rust_required
@pytest.mark.parametrize('pres,temp', [(2500, 200), (5000, 250)])
def test_sw_parity_sour_gas(pres, temp):
    """CO2+H2S+N2 mix in brine: Rust ≈ Python across the representative sour-gas case."""
    _sw_parity_check(pres, temp, ppm=50000, y_CO2=0.15, y_H2S=0.05, y_N2=0.03,
                     sg=0.75, rtol=RTOL_MEDIUM)


@rust_required
@pytest.mark.parametrize('pres,temp', [(1500, 180), (4000, 250)])
def test_sw_parity_natural_gas_only(pres, temp):
    """Pure natural gas in brine (HC-only): Rust ≈ Python."""
    _sw_parity_check(pres, temp, ppm=20000, y_CO2=0, y_H2S=0, y_N2=0,
                     sg=0.65, rtol=RTOL_MEDIUM)


@rust_required
@pytest.mark.parametrize('framework', ['dropin', 'sw_original'])
def test_sw_non_proposed_framework_not_downgraded(framework):
    """Regression: the Rust flash only implements the 'proposed' framework
    (MC-3 water alpha + proposed kij_AQ). 'dropin'/'sw_original' must take the
    Python path so they are not silently computed as 'proposed'.

    With Rust available, a non-proposed framework must (a) match its own
    forced-Python result and (b) differ from the 'proposed' result on a
    saline point where the framework choice matters.
    """
    from pyrestoolbox.brine import brine
    kwargs = dict(pres=3000, temp=200, ppm=50000, y_CO2=1.0, y_H2S=0, y_N2=0,
                  sg=0.65, metric=False)
    mix_rust = brine.SoreideWhitson(framework=framework, **kwargs)
    with force_python():
        mix_python = brine.SoreideWhitson(framework=framework, **kwargs)
    # (a) Rust-available path must not downgrade — equals forced-Python path.
    np.testing.assert_allclose(
        mix_rust.Rs_total, mix_python.Rs_total, rtol=RTOL_MEDIUM, atol=1e-6,
        err_msg=f"framework={framework} silently downgraded on Rust path",
    )
    # (b) The framework choice must actually change the answer vs 'proposed'.
    mix_proposed = brine.SoreideWhitson(framework='proposed', **kwargs)
    assert abs(mix_rust.Rs_total - mix_proposed.Rs_total) > 1e-4, (
        f"framework={framework} produced the 'proposed' result — gate ineffective"
    )


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
