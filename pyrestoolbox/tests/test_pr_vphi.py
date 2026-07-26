"""Tests for the S&W PR + VSHIFT route to dissolved-gas V2inf.

Ported from ~/projects/brine_props 2026-07-25 (the project was named
garcia_extension at the time). The reference values below are the brine_props
implementation, which is itself pinned against direct volumetric measurements
in that project's validation.py.
"""
import pytest

from pyrestoolbox.brine import pr_vphi
from pyrestoolbox.plyasunov import V2_inf as V_plyasunov


def test_vshift_values_pinned():
    """The shifts are load-bearing: refitting them requires regenerating the
    brine_props calibration (code/fit_pr_vshift.py) and updating both."""
    expected = {'CH4': -0.109632, 'CO2': -0.037913, 'H2S': -0.078975,
                'N2': -0.155288, 'H2': -0.177625, 'C2H6': -0.073142,
                # C3H8 is not from that calibration: it is the mean of two
                # direct 298 K determinations (see test below).
                'C3H8': -0.112963,
                # NC4H10 from Moore (1982) 76.6; the only positive shift.
                'NC4H10': +0.110924}
    assert pr_vphi.VSHIFT == pytest.approx(expected, abs=1e-9)


def test_peneloux_identity_is_exact():
    """V2inf = V2inf_raw - s*b, exactly. Guards the volume-shift algebra."""
    for gas in ('CO2', 'H2S', 'CH4'):
        raw = pr_vphi.V2_inf_raw(gas, 350.0, 20.0)
        shifted = pr_vphi.V2_inf(gas, 350.0, 20.0)
        implied = (raw - shifted) / pr_vphi.b_covolume(gas)
        assert implied == pytest.approx(pr_vphi.VSHIFT[gas], abs=1e-10)


@pytest.mark.parametrize('gas,T,P,expected', [
    ('CH4', 350.0, 20.0, 39.4072),
    ('CO2', 350.0, 20.0, 36.1979),
    ('H2S', 350.0, 20.0, 37.1617),
    ('N2', 350.0, 20.0, 37.1810),
    ('H2', 350.0, 20.0, 27.3687),
    ('C2H6', 350.0, 20.0, 55.5740),
])
def test_reference_values(gas, T, P, expected):
    assert pr_vphi.V2_inf(gas, T, P) == pytest.approx(expected, abs=1e-3)


@pytest.mark.parametrize('gas,T,P,measured', [
    # Hnedkovsky 1996 vibrating-tube densimetry, the calibration basis
    ('CH4', 373.15, 28.0, 40.7),
    ('CO2', 373.15, 20.0, 37.8),
    ('H2S', 373.15, 20.0, 39.0),
    ('CH4', 473.15, 28.0, 54.4),
    ('CO2', 473.15, 20.0, 50.0),
    ('H2S', 473.15, 20.0, 49.2),
])
def test_agrees_with_hnedkovsky_densimetry(gas, T, P, measured):
    v = pr_vphi.V2_inf(gas, T, P)
    assert abs(100 * (v / measured - 1)) < 3.0


def test_recovers_h2s_temperature_trend():
    """The point of carrying this route. Murphy & Gaines 1974 measured H2S V_phi
    RISING with temperature; the Plyasunov correlation is flat there."""
    lo, hi = 296.5, 310.0
    trend_pr = pr_vphi.V2_inf('H2S', hi, 0.5) - pr_vphi.V2_inf('H2S', lo, 0.5)
    trend_ply = (float(V_plyasunov('H2S', hi, 0.5))
                 - float(V_plyasunov('H2S', lo, 0.5)))
    assert trend_pr > 0.2
    assert abs(trend_ply) < 0.1


@pytest.mark.parametrize('T,P', [
    # 500 K is NO LONGER refused (the 473 K handover was removed 2026-07-25);
    # T_MAX is now the IF97 Region 1 ceiling.
    (624.0, 20.0),    # above T_MAX
    (250.0, 20.0),    # below T_MIN
    (350.0, 150.0),   # above P_MAX
    (450.0, 0.2),     # below the water saturation pressure
])
def test_range_guards_raise(T, P):
    with pytest.raises(ValueError):
        pr_vphi.V2_inf('CO2', T, P)


def test_unsupported_gas_raises():
    # C3H8 and NC4H10 both gained shifts 2026-07-25 and are supported now.
    with pytest.raises(ValueError):
        pr_vphi.V2_inf('C5H12', 350.0, 20.0)


def test_water_alpha_choice_is_immaterial():
    """The S&W framework uses alpha_water_soreide for sw_original/dropin and
    Mathias-Copeman otherwise. The shift is transferable: swapping the alpha
    moves V2inf by well under the 0.5-1.4% fit error, which is why one VSHIFT
    set serves every framework and the Dart port too."""
    from pyrestoolbox.brine import _lib_vle_engine as E
    for gas in ('CO2', 'CH4', 'H2S'):
        for T in (298.15, 373.15, 473.15):
            v = pr_vphi.V2_inf(gas, T, 20.0)
            # perturb: MC3 alpha instead of the S&W one
            orig = E.alpha_water_soreide
            try:
                E.alpha_water_soreide = lambda Tr, m=0.0: E.alpha_water_mc3(Tr)
                v_mc3 = pr_vphi.V2_inf(gas, T, 20.0)
            finally:
                E.alpha_water_soreide = orig
            assert abs(v_mc3 - v) < 0.10


def test_agrees_with_plyasunov_within_stated_spread():
    """Two independent routes to the same quantity should not disagree wildly.
    They differ most for the gases whose shift rests on 298 K data alone."""
    for gas in ('CO2', 'CH4', 'H2S'):
        for T in (300.0, 373.15, 450.0):
            a = pr_vphi.V2_inf(gas, T, 20.0)
            b = float(V_plyasunov(gas, T, 20.0))
            assert abs(100 * (a / b - 1)) < 8.0


# ---------------------------------------------------------------------------
# V_phi route dispatch (added 2026-07-25, when the delivered density default
# moved from Plyasunov to the S&W modified-PR route)
# ---------------------------------------------------------------------------
def test_vphi_route_dispatch_inside_box():
    """Inside the PR calibration box, 'auto' must equal 'pr'."""
    from pyrestoolbox.brine import vphi_route as vr

    for gas in ("CH4", "CO2", "H2S", "N2", "H2", "C2H6"):
        assert vr.route_used(gas, 373.15, 30.0) == "pr"
        assert vr.V_phi(gas, 373.15, 30.0) == vr.V_phi(gas, 373.15, 30.0, "pr")


def test_vphi_route_falls_back_only_where_it_must():
    """Only unsupported gases and undefined states fall back.

    The 473 K handover was REMOVED 2026-07-25: PR now runs to the IF97 Region 1
    ceiling, so there is no route change anywhere in the usable range.
    """
    from pyrestoolbox.brine import vphi_route as vr

    # NC4H10 alone has no fitted shift and is absent from the component set.
    assert vr.route_used("C5H12", 298.15, 20.0) == "plyasunov"
    # C3H8 gained a shift 2026-07-25 and must no longer fall back.
    assert vr.route_used("C3H8", 298.15, 20.0) == "pr"
    assert vr.route_used("CO2", 624.0, 30.0) == "plyasunov"    # past IF97 R1
    assert vr.route_used("CO2", 260.0, 20.0) == "plyasunov"    # below 273.15 K
    # No handover anywhere in between.
    for T in (300.0, 400.0, 450.0, 473.0, 523.15, 600.0, 620.0):
        assert vr.route_used("CO2", T, 30.0) == "pr", f"route changed at {T} K"


def test_vphi_routes_agree_inside_the_vouched_envelope():
    """Inside the range accuracy is claimed for, the route choice must not matter.

    There is no longer a handover to justify, but this pins the claim made in
    the deliverables: within 300-450 K the two routes agree on delivered brine
    density to 0.12%, so which one is used is not what limits accuracy there.
    """
    from pyrestoolbox.brine import vphi_route as vr
    from pyrestoolbox.plyasunov import gas_mw

    MW_WATER, S, P, T = 18.015268, 0.05, 20.0, 450.0
    rho1 = 928.0                       # kg/m3, gas-free brine at these conditions

    def rho(gas, v, x2):
        W = (1.0 - x2) * MW_WATER / (1.0 - S)
        return (W + x2 * gas_mw(gas)) / (W / (rho1 / 1000.0) + x2 * v) * 1000.0

    worst = 0.0
    for gas in ("CH4", "CO2", "H2S", "N2", "H2", "C2H6"):
        x2 = 0.02 if gas in ("CO2", "H2S") else 0.008
        a = rho(gas, vr.V_phi(gas, T, P, "plyasunov"), x2)
        b = rho(gas, vr.V_phi(gas, T, P, "pr"), x2)
        worst = max(worst, abs(100.0 * (b / a - 1.0)))
    assert worst < 0.20, f"route disagreement grew to {worst:.3f}%"


def test_sw_vphi_route_selectable_and_rejects_junk():
    """SoreideWhitson must expose the route and reject an unknown one."""
    import pytest

    from pyrestoolbox import brine

    kw = dict(pres=345.0, temp=100.0, ppm=50000, y_CO2=1.0, metric=True)
    auto = brine.SoreideWhitson(**kw)
    pr = brine.SoreideWhitson(**kw, vphi_route="pr")
    ply = brine.SoreideWhitson(**kw, vphi_route="plyasunov")

    assert auto.vphi_route == "auto"
    assert auto.bDen[0] == pr.bDen[0]              # inside the box, auto == pr
    assert abs(auto.bDen[0] - ply.bDen[0]) > 1e-9  # and really differs from old
    with pytest.raises(ValueError):
        brine.SoreideWhitson(**kw, vphi_route="nonsense")


def test_c3h8_shift_sits_between_its_two_anchors():
    """C3H8 is the weakest entry and must not drift silently.

    Its shift is not fitted to a data set: no densimetry for propane in water
    exists. It is the mean of the two direct 298 K determinations available,
    Moore (1982) 70.7 and Zhou & Battino (2001) 75.0 cm3/mol, which disagree by
    6.1%. That spread is the honest uncertainty on this gas.
    """
    from pyrestoolbox.brine import pr_vphi, vphi_route as vr

    v = vr.V_phi("C3H8", 298.15, 20.0, "pr")
    assert abs(v - 72.85) < 0.02, f"C3H8 V_phi moved to {v:.4f}"
    assert 70.7 < v < 75.0, "must sit inside the bracket its anchors define"
    # The superseded fallback sat below both anchors; that is the bar cleared.
    assert vr.V_phi("C3H8", 298.15, 20.0, "plyasunov") < 70.7
    # No temperature-resolved data at all, so it must stay declared.
    assert "C3H8" in pr_vphi.UNCALIBRATED_IN_T


def test_nc4h10_reproduces_moores_measurement():
    """Butane is the only gas needing a POSITIVE shift, and its evidence is thin.

    Moore (1982) Table I gives 76.6 +/- 0.1 cm3/mol from two runs, against his
    own stated overall imprecision of +/-1.5. It also breaks his own homologous
    series: his CH2 increments run +18.4, +17.8, then +5.9. It is adopted
    because it is the only direct measurement, not because it is well founded.
    """
    from pyrestoolbox.brine import pr_vphi, vphi_route as vr

    assert abs(vr.V_phi("NC4H10", 298.15, 20.0, "pr") - 76.6) < 0.02
    assert pr_vphi.VSHIFT["NC4H10"] > 0
    assert all(v < 0 for g, v in pr_vphi.VSHIFT.items() if g != "NC4H10")
    assert "NC4H10" in pr_vphi.UNCALIBRATED_IN_T
