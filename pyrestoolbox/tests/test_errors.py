"""
Error-case tests using pytest.raises.

Covers key validation paths across all modules.
Ensures bad inputs raise clear exceptions rather than silently producing garbage.
"""

import pytest
import numpy as np


# ── shared_fns / validate ──────────────────────────────────────────

class TestValidation:
    def test_validate_bad_method_string(self):
        from pyrestoolbox.validate import validate_methods
        with pytest.raises(ValueError, match="incorrect method"):
            validate_methods(["zmethod"], ["NOSUCH"])

    def test_validate_pe_negative_pressure(self):
        from pyrestoolbox.gas import gas
        with pytest.raises(ValueError, match="[Pp]ressure"):
            gas.gas_z(p=-100, sg=0.7, degf=200)

    def test_validate_pe_negative_sg(self):
        from pyrestoolbox.gas import gas
        with pytest.raises(ValueError, match="[Ss]pecific gravity"):
            gas.gas_z(p=2000, sg=-0.1, degf=200)

    def test_validate_pe_below_absolute_zero(self):
        from pyrestoolbox.gas import gas
        with pytest.raises(ValueError, match="[Tt]emperature"):
            gas.gas_z(p=2000, sg=0.7, degf=-500)

    def test_validate_pe_mole_fraction_negative(self):
        from pyrestoolbox.gas import gas
        with pytest.raises(ValueError, match="[Mm]ole fraction"):
            gas.gas_z(p=2000, sg=0.7, degf=200, co2=-0.1)

    def test_validate_pe_mole_fractions_exceed_one(self):
        from pyrestoolbox.gas import gas
        with pytest.raises(ValueError, match="[Ff]raction"):
            gas.gas_z(p=2000, sg=0.7, degf=200, co2=0.6, h2s=0.5)


# ── gas module ─────────────────────────────────────────────────────

class TestGasErrors:
    def test_sut_pure_inerts_raises(self):
        """SUT critical property method requires HC fraction > 0."""
        from pyrestoolbox.gas import gas
        with pytest.raises(ValueError, match="hydrocarbon"):
            gas.gas_tc_pc(sg=0.7, co2=0.5, h2s=0.3, n2=0.2, cmethod='SUT')

    def test_gas_rate_radial_bad_rw(self):
        from pyrestoolbox.gas import gas
        with pytest.raises(ValueError, match="r_w"):
            gas.gas_rate_radial(
                k=10, h=50, pr=4000, pwf=3000, r_w=0, r_ext=1000,
                degf=200, sg=0.7
            )

    def test_gas_rate_radial_re_less_than_rw(self):
        from pyrestoolbox.gas import gas
        with pytest.raises(ValueError, match="r_ext"):
            gas.gas_rate_radial(
                k=10, h=50, pr=4000, pwf=3000, r_w=0.5, r_ext=0.3,
                degf=200, sg=0.7
            )

    def test_gas_rate_linear_bad_length(self):
        from pyrestoolbox.gas import gas
        with pytest.raises(ValueError, match="[Ll]ength"):
            gas.gas_rate_linear(
                k=10, pr=4000, pwf=3000, area=10000, length=0,
                degf=200, sg=0.7
            )

    def test_gas_grad2sg_negative_gradient(self):
        from pyrestoolbox.gas import gas
        with pytest.raises(ValueError, match="[Gg]radient"):
            gas.gas_grad2sg(grad=-0.01, degf=200, p=3000)

    def test_hydrate_negative_pressure(self):
        from pyrestoolbox.gas import gas
        with pytest.raises(ValueError, match="[Pp]ressure"):
            gas.gas_hydrate(p=-100, degf=200, sg=0.7)

    def test_hydrate_bad_inhibitor_wt(self):
        from pyrestoolbox.gas import gas
        from pyrestoolbox.classes import inhibitor
        with pytest.raises(ValueError, match="[Ii]nhibitor"):
            gas.gas_hydrate(p=1000, degf=200, sg=0.7,
                            inhibitor_type=inhibitor.MEOH, inhibitor_wt_pct=110)


# ── oil module ─────────────────────────────────────────────────────

class TestOilErrors:
    def test_oil_rate_radial_missing_pvt(self):
        from pyrestoolbox.oil import oil
        with pytest.raises(ValueError):
            oil.oil_rate_radial(
                k=50, h=30, pr=4000, pwf=2000,
                r_w=0.5, r_ext=1000, degf=180
            )

    def test_oil_rate_radial_bad_rw(self):
        from pyrestoolbox.oil import oil
        with pytest.raises(ValueError, match="r_w"):
            oil.oil_rate_radial(
                k=50, h=30, pr=4000, pwf=2000,
                r_w=-1, r_ext=1000, uo=1.0, bo=1.2
            )

    def test_oil_harmonize_no_sg(self):
        """oil_harmonize requires at least one of sg_g or sg_sp."""
        from pyrestoolbox.oil import oil
        with pytest.raises(ValueError, match="sg_g.*sg_sp|sg_sp.*sg_g"):
            oil.oil_harmonize(api=35, degf=180, pb=2500, sg_g=0, sg_sp=0)

    def test_oil_rate_linear_bad_length(self):
        from pyrestoolbox.oil import oil
        with pytest.raises(ValueError, match="[Ll]ength"):
            oil.oil_rate_linear(
                k=50, pr=4000, pwf=2000,
                length=0, area=10000, uo=1.0, bo=1.2
            )


# ── brine module ───────────────────────────────────────────────────

class TestBrineErrors:
    def test_brine_props_negative_pressure(self):
        from pyrestoolbox.brine import brine
        with pytest.raises(ValueError, match="[Pp]ressure"):
            brine.brine_props(degf=200, p=-100)

    def test_brine_props_salt_out_of_range(self):
        from pyrestoolbox.brine import brine
        with pytest.raises(ValueError, match="[Ss]alt"):
            brine.brine_props(degf=200, p=3000, wt=100)

    def test_sw_negative_ppm(self):
        from pyrestoolbox.brine import brine
        with pytest.raises(ValueError, match="ppm"):
            brine.SoreideWhitson(pres=3000, temp=200, ppm=-1, metric=False)

    def test_sw_ppm_too_high(self):
        from pyrestoolbox.brine import brine
        with pytest.raises(ValueError, match="ppm"):
            brine.SoreideWhitson(pres=3000, temp=200, ppm=1_500_000, metric=False)

    def test_co2_brine_negative_pressure(self):
        from pyrestoolbox.brine import brine
        with pytest.raises(ValueError, match="[Pp]ressure"):
            brine.CO2_Brine_Mixture(pres=-100, temp=200, ppm=10000, metric=False)

    def test_co2_brine_gas_fractions_exceed_one(self):
        from pyrestoolbox.brine import brine
        with pytest.raises(ValueError, match="[Ff]raction"):
            brine.SoreideWhitson(pres=3000, temp=200, ppm=10000,
                                 y_CO2=0.6, y_H2S=0.5, metric=False)

    def test_sw_negative_pressure(self):
        from pyrestoolbox.brine import brine
        with pytest.raises(ValueError, match="[Pp]ressure"):
            brine.SoreideWhitson(pres=-100, temp=200, ppm=10000, metric=False)


# ── dca module ─────────────────────────────────────────────────────

class TestDCAErrors:
    def test_arps_rate_negative_qi(self):
        from pyrestoolbox.dca import dca
        with pytest.raises(ValueError, match="qi"):
            dca.arps_rate(qi=-100, di=0.1, b=0.5, t=np.array([1, 2, 3]))

    def test_arps_rate_negative_di(self):
        from pyrestoolbox.dca import dca
        with pytest.raises(ValueError, match="di"):
            dca.arps_rate(qi=1000, di=-0.1, b=0.5, t=np.array([1, 2, 3]))

    def test_arps_rate_b_out_of_range(self):
        from pyrestoolbox.dca import dca
        with pytest.raises(ValueError, match="b"):
            dca.arps_rate(qi=1000, di=0.1, b=1.5, t=np.array([1, 2, 3]))

    def test_arps_cum_negative_qi(self):
        from pyrestoolbox.dca import dca
        with pytest.raises(ValueError, match="qi"):
            dca.arps_cum(qi=-100, di=0.1, b=0.5, t=np.array([1, 2, 3]))

    def test_duong_rate_bad_a(self):
        from pyrestoolbox.dca import dca
        with pytest.raises(ValueError, match="'a'"):
            dca.duong_rate(qi=1000, a=-1, m=1.5, t=np.array([1, 2, 3]))

    def test_duong_rate_m_less_than_one(self):
        from pyrestoolbox.dca import dca
        with pytest.raises(ValueError, match="'m'"):
            dca.duong_rate(qi=1000, a=1.0, m=0.5, t=np.array([1, 2, 3]))

    def test_eur_qmin_exceeds_qi(self):
        from pyrestoolbox.dca import dca
        with pytest.raises(ValueError, match="q_min"):
            dca.eur(qi=100, di=0.1, b=0.5, q_min=200)

    def test_fit_decline_too_few_points(self):
        from pyrestoolbox.dca import dca
        with pytest.raises(ValueError, match="3 data points"):
            dca.fit_decline(t=np.array([1, 2]), q=np.array([100, 90]))

    def test_fit_decline_negative_rates(self):
        from pyrestoolbox.dca import dca
        with pytest.raises(ValueError, match="positive"):
            dca.fit_decline(
                t=np.array([1, 2, 3, 4]),
                q=np.array([100, 90, -10, 70])
            )

    def test_fit_decline_mismatched_lengths(self):
        from pyrestoolbox.dca import dca
        with pytest.raises(ValueError, match="same length"):
            dca.fit_decline(t=np.array([1, 2, 3]), q=np.array([100, 90]))


# ── matbal module ──────────────────────────────────────────────────

class TestMatbalErrors:
    def test_gas_matbal_negative_pressure(self):
        from pyrestoolbox.matbal import matbal
        with pytest.raises(ValueError, match="[Pp]ressure"):
            matbal.gas_matbal(
                p=np.array([-100, 2000]),
                Gp=np.array([0, 100]),
                degf=200,
                pvt_table={'p': [1000, 3000], 'Z': [0.9, 0.8]}
            )

    def test_gas_matbal_mismatched_arrays(self):
        from pyrestoolbox.matbal import matbal
        with pytest.raises(ValueError, match="same length"):
            matbal.gas_matbal(
                p=np.array([3000, 2500, 2000]),
                Gp=np.array([0, 50]),
                degf=200,
                pvt_table={'p': [1000, 3000], 'Z': [0.9, 0.8]}
            )

    def test_gas_matbal_too_few_points(self):
        from pyrestoolbox.matbal import matbal
        with pytest.raises(ValueError, match="2"):
            matbal.gas_matbal(
                p=np.array([3000]),
                Gp=np.array([0]),
                degf=200,
                pvt_table={'p': [1000, 3000], 'Z': [0.9, 0.8]}
            )

    def test_gas_matbal_both_z_and_bg(self):
        from pyrestoolbox.matbal import matbal
        with pytest.raises(ValueError, match="not both"):
            matbal.gas_matbal(
                p=np.array([3000, 2500]),
                Gp=np.array([0, 100]),
                degf=200,
                pvt_table={'p': [1000, 3000], 'Z': [0.9, 0.8], 'Bg': [0.005, 0.006]}
            )

    def test_gas_matbal_missing_p_in_pvt(self):
        from pyrestoolbox.matbal import matbal
        with pytest.raises(ValueError, match="'p'"):
            matbal.gas_matbal(
                p=np.array([3000, 2500]),
                Gp=np.array([0, 100]),
                degf=200,
                pvt_table={'Z': [0.9, 0.8]}
            )

    def test_oil_matbal_bad_regress_key(self):
        from pyrestoolbox.matbal import matbal
        with pytest.raises(ValueError, match="regress"):
            matbal.oil_matbal(
                p=np.array([4000, 3500, 3000]),
                Np=np.array([0, 100, 200]),
                api=35, sg_sp=0.8, degf=200, pb=3500,
                regress={'nosuch': (0, 1)}
            )


# ── nodal module ───────────────────────────────────────────────────

class TestNodalErrors:
    def test_wellsegment_negative_md(self):
        from pyrestoolbox.nodal.nodal import WellSegment
        with pytest.raises(ValueError, match="[Mm]easured depth|md"):
            WellSegment(md=-1, id=2.441)

    def test_wellsegment_negative_id(self):
        from pyrestoolbox.nodal.nodal import WellSegment
        with pytest.raises(ValueError, match="[Ii]nternal diameter|id"):
            WellSegment(md=8000, id=-1)

    def test_wellsegment_deviation_out_of_range(self):
        from pyrestoolbox.nodal.nodal import WellSegment
        with pytest.raises(ValueError, match="[Dd]eviation"):
            WellSegment(md=8000, id=2.441, deviation=95)

    def test_completion_missing_params(self):
        from pyrestoolbox.nodal.nodal import Completion
        with pytest.raises(ValueError, match="[Mm]ust provide"):
            Completion()

    def test_completion_empty_segments(self):
        from pyrestoolbox.nodal.nodal import Completion
        with pytest.raises(ValueError, match="empty"):
            Completion(segments=[], tht=100, bht=250)

    def test_reservoir_bad_pr(self):
        from pyrestoolbox.nodal.nodal import Reservoir
        with pytest.raises(ValueError, match="[Pp]ressure|pr"):
            Reservoir(pr=-100, degf=200, k=10, h=50, rw=0.328, re=1000)

    def test_reservoir_bad_k(self):
        from pyrestoolbox.nodal.nodal import Reservoir
        with pytest.raises(ValueError, match="[Pp]ermeability|k"):
            Reservoir(pr=4000, degf=200, k=-1, h=50, rw=0.328, re=1000)


# ── simtools module ────────────────────────────────────────────────

class TestSimtoolsErrors:
    def test_relperm_bad_table_type(self):
        from pyrestoolbox.simtools import simtools
        with pytest.raises(ValueError, match="[Tt]able|SWOF|SGOF"):
            simtools.rel_perm_table(
                rows=20, krtable='BADTYPE',
                swc=0.2, swcr=0.2, sorw=0.2, sorg=0.2,
                sgcr=0.05, no=2.0, nw=2.0, ng=2.0
            )

    def test_relperm_too_few_rows(self):
        from pyrestoolbox.simtools import simtools
        with pytest.raises(ValueError, match="rows"):
            simtools.rel_perm_table(
                rows=1, krtable='SWOF',
                swc=0.2, swcr=0.2, sorw=0.2, sorg=0.2,
                sgcr=0.05, no=2.0, nw=2.0, ng=2.0
            )

    def test_rr_solver_empty_arrays(self):
        from pyrestoolbox.simtools import simtools
        with pytest.raises(ValueError, match="empty"):
            simtools.rr_solver(zi=np.array([]), ki=np.array([]))

    def test_rr_solver_mismatched(self):
        from pyrestoolbox.simtools import simtools
        with pytest.raises(ValueError, match="same length"):
            simtools.rr_solver(zi=np.array([0.5, 0.5]), ki=np.array([2.0]))


# ── layer module ───────────────────────────────────────────────────

class TestLayerErrors:
    def test_lorenz_out_of_range(self):
        from pyrestoolbox.layer import layer
        with pytest.raises(ValueError, match="[Ll]orenz"):
            layer.lorenz2b(lorenz=1.5)

    def test_bad_method(self):
        from pyrestoolbox.layer import layer
        with pytest.raises(ValueError, match="[Mm]ethod"):
            layer.lorenz2b(lorenz=0.5, lrnz_method='BAD')


# ── library module ─────────────────────────────────────────────────

class TestLibraryErrors:
    def test_bad_component(self):
        from pyrestoolbox.library import library
        with pytest.raises(KeyError):
            library.prop('NOSUCH_COMPONENT', 'SRK', 'Tc')

    def test_bad_model(self):
        from pyrestoolbox.library import library
        with pytest.raises(ValueError, match="[Mm]odel"):
            library.prop('C1', 'NOSUCH_MODEL', 'Tc')


# ── correlation validity range warnings ────────────────────────────

class TestRangeWarnings:
    """Verify non-blocking warnings fire for out-of-range inputs."""

    def test_dak_tr_out_of_range(self):
        from pyrestoolbox.gas import gas
        with pytest.warns(UserWarning, match="DAK.*Tr.*outside"):
            gas.gas_z(p=500, sg=0.7, degf=800)  # Very high T -> Tr > 3

    def test_dak_ppr_out_of_range(self):
        from pyrestoolbox.gas import gas
        with pytest.warns(UserWarning, match="DAK.*Ppr.*outside"):
            gas.gas_z(p=50000, sg=0.7, degf=200)  # Very high P -> Ppr > 30

    def test_hy_tr_out_of_range(self):
        from pyrestoolbox.gas import gas
        with pytest.warns(UserWarning, match="HY.*Tr.*outside"):
            gas.gas_z(p=2000, sg=0.7, degf=800, zmethod='HY')

    def test_standing_pb_api_out_of_range(self):
        from pyrestoolbox.oil import oil
        with pytest.warns(UserWarning, match="Standing.*API.*outside"):
            oil.oil_pbub(api=5, degf=180, rsb=500, sg_g=0.8, pbmethod='STAN')

    def test_beggs_robinson_temp_out_of_range(self):
        from pyrestoolbox.oil import oil
        with pytest.warns(UserWarning, match="Beggs-Robinson.*outside"):
            oil.oil_viso(p=2000, api=35, degf=400, pb=2500, rs=500)

    def test_no_warning_in_range(self):
        """Normal-range inputs should not produce warnings."""
        import warnings as w
        from pyrestoolbox.gas import gas
        with w.catch_warnings():
            w.simplefilter("error")
            gas.gas_z(p=2000, sg=0.7, degf=200)
