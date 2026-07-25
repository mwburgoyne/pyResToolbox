"""Which V_phi route supplies the dissolved-gas molar volume.

DEFAULT CHANGED 2026-07-25: the delivered gas-saturated brine density uses the
Soreide-Whitson modified-PR route with one dimensionless volume shift per gas,
falling back to the Plyasunov A12inf model where PR is not calibrated. Plyasunov
was the default only because it was implemented first.

WHY PR IS PREFERRED. Against direct volumetric measurement PR wins on five of
six gases (in-sample MAE: CH4 1.4 vs 1.5%, H2S 0.5 vs 1.4%, N2 3.3 vs 4.9%,
H2 5.1 vs 6.2%, C2H6 3.0 vs 3.4%; it loses CO2, 0.9 vs 0.7%), and it is the only
route that reproduces the H2S temperature trend - Murphy and Gaines measured
+0.69 cm3/mol between their cold and warm groups, PR gives +0.53, Plyasunov is
flat at +0.02. It also costs one number per gas rather than 35 at eight
significant digits plus IAPWS-IF97.

NO TEMPERATURE HANDOVER (2026-07-25). An earlier version fell back to Plyasunov
above 473 K, the top of the calibration data. It no longer does: PR runs to the
IAPWS-IF97 Region 1 ceiling of 623.15 K. A mid-range handover bought nothing and
put a discontinuity in the delivered default.

That is not an accuracy claim. Accuracy is vouched for to about 300-350 degF
(420-450 K); above that neither route is calibrated and neither is claimed. The
practical envelope therefore sits inside the calibration data, and within it the
two routes agree on brine density to 0.12% at worst.

WHERE PR STILL CANNOT BE USED, and what happens instead:
  * C3H8 and NC4H10 have no fitted volume shift, so they always use Plyasunov.
  * Outside 273.15-623.15 K, above 100 MPa, or below the water saturation
    pressure, where V2inf is undefined because the solvent is vapour.

CAVEAT. N2, H2 and C2H6 have shifts fitted to a 298 K cluster only, so for those
three the temperature dependence of the default route is the EOS unaided rather
than calibrated. Plyasunov stays available via route='plyasunov'.
"""

from pyrestoolbox.plyasunov import V_phi as _plyasunov_V_phi

from . import pr_vphi as _pr

ROUTES = ("auto", "pr", "plyasunov")
DEFAULT_ROUTE = "auto"


def pr_applicable(gas: str, T: float, P: float) -> bool:
    """True if the PR + VSHIFT route is calibrated for this gas and state."""
    try:
        _pr._check(gas, T, P)
        return True
    except (ValueError, KeyError):
        return False


def route_used(gas: str, T: float, P: float, route: str = DEFAULT_ROUTE) -> str:
    """Which route V_phi would actually use. For reporting, never for logic."""
    if route in ("pr", "plyasunov"):
        return route
    return "pr" if pr_applicable(gas, T, P) else "plyasunov"


def V_phi(gas: str, T: float, P: float, route: str = DEFAULT_ROUTE,
          m_nacl: float = 0.0) -> float:
    """Dissolved-gas apparent molar volume at infinite dilution, cm3/mol.

    Args:
        gas: gas name ('CO2', 'CH4', 'H2S', 'N2', 'H2', 'C2H6', 'C3H8', 'NC4H10')
        T: temperature, K
        P: pressure, MPa
        route: 'auto' uses PR where calibrated and Plyasunov elsewhere; 'pr'
            forces PR and raises outside its box; 'plyasunov' forces A12inf.
        m_nacl: NaCl molality, passed to the PR route only, where it enters the
            S&W water alpha and kij. There is no salinity term on the volume
            shift itself; leave at 0 for the calibrated behaviour.
    """
    if route not in ROUTES:
        raise ValueError(f"route must be one of {ROUTES}, got {route!r}")
    if route == "plyasunov":
        return float(_plyasunov_V_phi(gas, T, P))
    if route == "pr":
        return float(_pr.V2_inf(gas, T, P, m_nacl))
    return (float(_pr.V2_inf(gas, T, P, m_nacl)) if pr_applicable(gas, T, P)
            else float(_plyasunov_V_phi(gas, T, P)))
