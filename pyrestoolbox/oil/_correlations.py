"""Core oil PVT correlations: Pb, Rs, Bo, viscosity."""

import warnings
import numpy as np

from pyrestoolbox.constants import (
    psc, BAR_TO_PSI, PSI_TO_BAR, degc_to_degf,
    SM3_PER_SM3_TO_SCF_PER_STB, SCF_PER_STB_TO_SM3_PER_SM3,
)
from pyrestoolbox.classes import pb_method, rs_method, bo_method, deno_method
from pyrestoolbox.validate import validate_methods
from pyrestoolbox.shared_fns import validate_pe_inputs

from ._constants import (
    _STAN_T_COEFF, _STAN_API_COEFF, _STAN_DENOM, _STAN_SG_EXP, _STAN_OFFSET,
    _VEL_X_C0, _VEL_X_T_EXP, _VEL_X_API_C, _VEL_X_API_EXP,
    _VEL_PB_MULT, _VEL_RS_EXP, _VEL_SG_EXP, _VEL_PB_OFFSET, _VEL_PB_POW,
    _VEL_RS_A, _VEL_RS_B, _VEL_RS_C,
    _VALMC_C, _VALMC_LNPB_A0, _VALMC_LNPB_A1, _VALMC_LNPB_A2,
    _BR_Z0, _BR_Z_API, _BR_T_EXP, _BR_A_MULT, _BR_A_EXP, _BR_A_OFFSET,
    _BR_B_MULT, _BR_B_EXP, _BR_B_OFFSET,
    _PF_POLY, _PF_P_COEFF,
    _BO_STAN_A, _BO_STAN_B, _BO_STAN_C, _BO_STAN_EXP,
    _BO_MC_WDEN, _BO_MC_RS_COEFF,
)
from ._utils import check_sgs, get_real_part, oil_sg
from ._density import oil_deno


def oil_pbub(
    api: float,
    degf: float,
    rsb: float,
    sg_g: float = 0,
    sg_sp: float = 0,
    pbmethod: pb_method = pb_method.VALMC,
    metric: bool = False,
) -> float:
    """ Returns bubble point pressure (psia | barsa) calculated with different correlations

        api: Stock tank oil density (deg API)
        degf: Reservoir Temperature (deg F | deg C)
        rsb: Oil solution gas volume at Pbub (scf/stb | sm3/sm3)
        sg_sp: Separator Gas specific Gravity (relative to air) <-- Required for Valko McCain & Velarde
        sg_g: Weighted average specific gravity of surface gas (relative to air). <-- Required for Standing
        pbmethod: A string or pb_method Enum class that specifies one of following calculation choices;
                   STAN: Standing Correlation (1947)
                   VALMC: Valko-McCain Correlation (2003) - https://www.sciencedirect.com/science/article/abs/pii/S0920410502003194
                   VELAR: Velarde, Blasingame & McCain (1997)
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, sm3/sm3). Defaults to False (FIELD)
    """
    if metric:
        degf = degc_to_degf(degf)
        rsb = rsb * SM3_PER_SM3_TO_SCF_PER_STB

    validate_pe_inputs(degf=degf)

    sg_g, sg_sp = check_sgs(sg_g=sg_g, sg_sp=sg_sp)

    pbmethod = validate_methods(["pbmethod"], [pbmethod])

    if pbmethod.name == "STAN":
        if rsb * api * sg_g * degf == 0:
            raise ValueError(
                "Need valid values for rs, api, sg_g and degf for Standing Pb calculation"
            )
    else:
        if rsb * api * sg_sp * degf == 0:
            raise ValueError(
                "Need valid values for rsb, api, sg_sp and degf for Velarde or Valko McCain Pb calculation"
            )

    # Correlation validity range warnings (non-blocking)
    if pbmethod.name == "STAN":
        if api < 16.5 or api > 63.8:
            warnings.warn(f"Standing Pb: API={api:.1f} outside calibration range [16.5, 63.8]", stacklevel=2)
        if degf < 100 or degf > 258:
            warnings.warn(f"Standing Pb: T={degf:.1f}\u00b0F outside calibration range [100, 258]", stacklevel=2)
    elif pbmethod.name == "VALMC":
        if api < 6 or api > 56.8:
            warnings.warn(f"Valko-McCain Pb: API={api:.1f} outside calibration range [6, 56.8]", stacklevel=2)
        if degf < 78 or degf > 330:
            warnings.warn(f"Valko-McCain Pb: T={degf:.1f}\u00b0F outside calibration range [78, 330]", stacklevel=2)

    def pbub_standing(
        api, degf, sg_g, rsb, sg_sp
    ) -> float:  # 1.63 in 'Oil & Gas Properties & Correlations' - http://dx.doi.org/10.1016/B978-0-12-803437-8.00001-4
        a = _STAN_T_COEFF * degf - _STAN_API_COEFF * api
        return (
            _STAN_DENOM * ((rsb / sg_g) ** _STAN_SG_EXP * 10 ** a - _STAN_OFFSET) + psc
        )  # Adding 14.7 as I suspect this is in psig

    def pbub_valko_mccain(api, degf, sg_g, rsb, sg_sp) -> float:
        extrap = False
        rsb2 = rsb
        if rsb <= 1: # Extrapolate to 14.696 psia at zero GOR
            extrap = True
            rsb = 1
        var = [np.log(rsb), api, sg_sp, degf]
        Zn = [sum([_VALMC_C[i][n] * var[n] ** i for i in range(4)]) for n in range(4)] # Eq 2-1
        Z = sum(Zn)
        lnpb = _VALMC_LNPB_A0 + _VALMC_LNPB_A1 * Z + _VALMC_LNPB_A2 * Z ** 2
        pb = np.exp(lnpb)

        if extrap:
            slope = (pb - 14.696)/(1 - 0)
            intercept = pb - 1*slope
            pb = slope * rsb2 + intercept

        return pb

    def pbub_velarde(api, degf, sg_g, rsb, sg_sp) -> float:
        x = _VEL_X_C0 * degf ** _VEL_X_T_EXP - _VEL_X_API_C * api ** _VEL_X_API_EXP

        rsb_lim = (_VEL_PB_OFFSET / (sg_sp ** _VEL_SG_EXP * 10 ** x))**(1/_VEL_RS_EXP) # If Rsb < than this value, then the term inside the pbp brackets of pbp goes negative and causes imaginary numbers when raised to a power
        if rsb < rsb_lim+1e-6:
            rsb = rsb_lim+1e-6

        pbp = (
            _VEL_PB_MULT
            * (rsb ** _VEL_RS_EXP * sg_sp ** _VEL_SG_EXP * 10 ** x - _VEL_PB_OFFSET)
            ** _VEL_PB_POW # psig
        )
        return pbp+psc # psia

    fn_dic = {
        "STAN": pbub_standing,
        "VALMC": pbub_valko_mccain,
        "VELAR": pbub_velarde,
    }

    result = fn_dic[pbmethod.name](
        api=api, degf=degf, sg_g=sg_g, rsb=rsb, sg_sp=sg_sp
    )
    if metric:
        return result * PSI_TO_BAR  # psia -> barsa
    return result

def oil_rs_bub(
    api: float,
    degf: float,
    pb: float,
    sg_g: float = 0,
    sg_sp: float = 0,
    rsmethod: rs_method = rs_method.VELAR,
    metric: bool = False,
) -> float:
    """ Returns Solution GOR (scf/stb | sm3/sm3) at bubble point pressure.
        Uses the inverse of the Bubble point pressure correlations, with the same method families
        Note: At low pressures, the VALMC method will fail (generally when Rsb < 10 scf/stb).
              The VALMC method will revert to the STAN method in these cases

        api: Stock tank oil density (deg API)
        degf: Reservoir Temperature (deg F | deg C)
        pb: Bubble point Pressure (psia | barsa)
        sg_sp: Separator Gas specific Gravity (relative to air) <-- Required for Valko McCain & Velarde
        sg_g: Weighted average specific gravity of surface gas (relative to air). <-- Required for Standing
        rsmethod: A string or pb_method Enum class that specifies one of following calculation choices. Note that VASBG is not available as it requires Pb as an input;
                   VELAR: Velarde, Blasingame & McCain (1997) - Default
                   STAN: Standing Correlation (1947), using form from https://www.sciencedirect.com/science/article/pii/B9780128034378000014
                   VALMC: Valko-McCain Correlation (2003) - https://www.sciencedirect.com/science/article/abs/pii/S0920410502003194
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, sm3/sm3). Defaults to False (FIELD)
    """
    if metric:
        degf = degc_to_degf(degf)
        pb = pb * BAR_TO_PSI

    validate_pe_inputs(degf=degf, p=pb)

    sg_g, sg_sp = check_sgs(sg_g=sg_g, sg_sp=sg_sp)

    pbmethod = 'VALMC' # Pbub calculations only needed with 'VALMC' rsmethod
    pbmethod, rsmethod = validate_methods(
        ["pbmethod", "rsmethod"], [pbmethod, rsmethod]
    )


    def rsbub_standing(api, degf, pb, sg_g, sg_sp) -> float:
        #print('Standing')
        a = _STAN_T_COEFF * degf - _STAN_API_COEFF * api  # Eq 1.64
        return sg_g * (((pb - psc) / _STAN_DENOM + _STAN_OFFSET) / 10 ** a) ** (
            1 / _STAN_SG_EXP
        )  # Eq 1.72 - Subtracting 14.7 as suspect this pressure in psig

    def rsbub_valko_mccain(api, degf, pb, sg_g, sg_sp) -> float:
        #print('Valko McCain')
        # Solve via iteration. First guess using Velarde Rsb, then simple Newton Iterations
        old_rsb = rsbub_velarde(api, degf, pb, sg_g, sg_sp)
        standing = False

        old_pbcalc = oil_pbub(degf=degf, api=api, sg_sp=sg_sp, rsb=old_rsb, pbmethod=pbmethod)
        old_err = old_pbcalc - pb
        new_rsb = old_rsb * pb / old_pbcalc
        i = 0
        new_err = 1000
        while abs(new_err) > 1e-5:
            i += 1
            new_pbcalc = oil_pbub(degf=degf,api=api,sg_sp=sg_sp,rsb=new_rsb,pbmethod=pbmethod)
            new_err = new_pbcalc - pb

            error_slope = (new_rsb - old_rsb) / (new_err - old_err)
            intcpt = new_rsb - error_slope * new_err

            old_err, old_rsb = new_err, new_rsb

            new_rsb = intcpt

            if (i > 100):
                import warnings
                warnings.warn("oil_rs_bub: Valko-McCain did not converge after 100 iterations", RuntimeWarning)
                return new_rsb
        return new_rsb

    def rsbub_velarde(api, degf, pb, sg_g, sg_sp) -> float:
        x = _VEL_X_C0 * degf ** _VEL_X_T_EXP - (_VEL_X_API_C * api ** _VEL_X_API_EXP) # Eq 14

        # Note that the Velarde approach predicts increasing Pbub with rs falling below 1 scf/stb, so below this we will linearly extrapolate to zero instead
        p_1scfstb = _VEL_PB_MULT*(sg_sp**_VEL_SG_EXP * 10**x - _VEL_PB_OFFSET)**_VEL_PB_POW # Eq 13
        psig = pb - psc

        if psig >= p_1scfstb:
            rsb = (-10**(-x)* sg_sp**(-_VEL_SG_EXP) *(-(psig/_VEL_PB_MULT)**(1/_VEL_PB_POW) - _VEL_PB_OFFSET))**(1/_VEL_RS_EXP) # Rearranged Eq 13
        else:
            slope = (1-0)/(p_1scfstb - 0)
            intercept = 1 - slope * p_1scfstb
            rsb = slope * psig + intercept
        return rsb

    fn_dic = {
        "STAN": rsbub_standing,
        "VALMC": rsbub_valko_mccain,
        "VELAR": rsbub_velarde,
    }

    rsbub = fn_dic[rsmethod.name](
        api=api, degf=degf, pb=pb, sg_g=sg_g, sg_sp=sg_sp
    )
    if np.isnan(rsbub):
        import warnings
        warnings.warn("oil_rs_bub: correlation returned NaN, returning 0", RuntimeWarning)
        return 0
    if metric:
        return rsbub * SCF_PER_STB_TO_SM3_PER_SM3  # scf/stb -> sm3/sm3
    return rsbub

def oil_rs(
    api: float,
    degf: float,
    sg_sp: float,
    p: float,
    pb: float = 0,
    rsb: float = 0,
    rsmethod: rs_method = rs_method.VELAR,
    pbmethod: pb_method = pb_method.VALMC,
    metric: bool = False,
) -> float:
    """ Returns solution gas oil ratio (scf/stb | sm3/sm3) calculated from different correlations. Either pb, rsb or both need to be specified.
        If one is missing, the other will be calculated from correlation.

        api: Stock tank oil density (deg API)
        degf: Reservoir Temperature (deg F | deg C)
        sg_sp: SG of separator gas
        p: Pressure of oil (psia | barsa)
        pb: Original bubble point pressure of oil (psia | barsa)
        rsb: Oil solution gas volume at bubblepoint pressure (scf/stb | sm3/sm3) <-- Required for Velarde, Blasingame & McCain
        rsmethod: A string or pb_method Enum class that specifies one of following calculation choices;
                   VELAR: Velarde, Blasingame & McCain (1997) - Default
                   STAN: Standing Correlation (1947), using form from https://www.sciencedirect.com/science/article/pii/B9780128034378000014
                   VALMC: Valko-McCain Correlation (2003)
        pbmethod: A string or pb_method Enum class that specifies one of following calculation choices;
                   STAN: Standing Correlation (1947)
                   VALMC: Valko-McCain Correlation (2003) - https://www.sciencedirect.com/science/article/abs/pii/S0920410502003194
                   VELAR: Velarde, Blasingame & McCain (1997) - Default
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, sm3/sm3). Defaults to False (FIELD)
    """
    if metric:
        degf = degc_to_degf(degf)
        p = p * BAR_TO_PSI
        if pb > 0:
            pb = pb * BAR_TO_PSI
        if rsb > 0:
            rsb = rsb * SM3_PER_SM3_TO_SCF_PER_STB

    validate_pe_inputs(p=p, degf=degf)
    if pb > 0:
        validate_pe_inputs(p=pb)

    pbmethod, rsmethod = validate_methods(
        ["pbmethod", "rsmethod"], [pbmethod, rsmethod]
    )

    sg_g, sg_sp = check_sgs(sg_g=0, sg_sp=sg_sp)

    if pb <= 0:  # Calculate Pb
        pb = oil_pbub(
            api=api, degf=degf, rsb=rsb, sg_sp=sg_sp, pbmethod=pbmethod
        )
    if rsb <= 0:  # Calculate rsb
        rsb = oil_rs_bub(
            api=api,
            degf=degf,
            pb=pb,
            sg_sp=sg_sp,
            rsmethod=rsmethod,
        )
        rsb = get_real_part(rsb)

    #print(rsb)

    if p >= pb:
        if metric:
            return rsb * SCF_PER_STB_TO_SM3_PER_SM3
        return rsb

    def Rs_velarde(
        api, degf, sg_g, sg_sp, p, pb, rsb
    ):  # Velarde, Blasingame & McCain (1997)
        # Equations 3.8a - 3.8f
        # Estimates Rs of depleting oil from separator oil observations
        pb = max(psc, pb)
        p = max(psc, p)

        if sg_sp * api * rsb == 0:
            raise ValueError(
                "Missing one of the required inputs: sg_sp, api, rsb, for the Velarde, Blasingame & McCain Rs calculation"
            )
        # Degenerate: pb at (or below) atmospheric means no dissolved gas can evolve.
        # Skip the correlation; the 0/0 in pr would otherwise return NaN silently.
        if pb - psc <= 0:
            return 0.0
        xs = [_VEL_RS_A, _VEL_RS_B, _VEL_RS_C]
        a = [
            x[0]
            * sg_sp ** x[1]
            * api ** x[2]
            * degf ** x[3]
            * (pb - psc) ** x[4]
            for x in xs
        ]
        pr = (p - psc) / (pb - psc)
        rsr = a[0] * pr ** a[1] + (1 - a[0]) * pr ** a[2] # Eq 3 from Velarde & Blasingame
        rs = rsb * rsr
        return get_real_part(rs)

    def rs_standing(api, degf, sg_g, sg_sp, p, pb, rsb):
        a = _STAN_T_COEFF * degf - _STAN_API_COEFF * api  # Eq 1.64
        return sg_g * (((p - psc) / _STAN_DENOM + _STAN_OFFSET) / 10 ** a) ** (
            1 / _STAN_SG_EXP
        )  # Eq 1.72 - Subtracting 14.7 as suspect this pressure in psig

    def rs_valko_mccain(api, degf, sg_g, sg_sp, p, pb, rsb):
        rsb_valko = oil_rs_bub(api, degf, pb, sg_g, sg_sp, rsmethod = 'VALMC') # Rsb from Valko-McCain approach
        rs_scaler = rsb / rsb_valko # Scalar to adjust calculated Valko McCain Rs back to be in line with the supplied Rsb
        return rs_scaler * oil_rs_bub(api, degf, p, sg_g, sg_sp, rsmethod = 'VALMC')


    fn_dic = {
        "VELAR": Rs_velarde,
        "STAN": rs_standing,
        "VALMC": rs_valko_mccain,
    }

    result = fn_dic[rsmethod.name](
        api=api, degf=degf, sg_g=sg_g, sg_sp=sg_sp, p=p, pb=pb, rsb=rsb
    )
    if metric:
        return result * SCF_PER_STB_TO_SM3_PER_SM3  # scf/stb -> sm3/sm3
    return result


def oil_bo(
    p: float,
    pb: float,
    degf: float,
    rs: float,
    rsb: float,
    sg_o: float,
    sg_g: float = 0,
    sg_sp: float = 0,
    bomethod: bo_method = bo_method.MCAIN,
    denomethod: deno_method = deno_method.SWMH,
    metric: bool = False,
) -> float:
    """ Returns oil formation volume factor (rb/stb | rm3/sm3) calculated with different correlations

        p: Pressure (psia | barsa)
        pb: Bubble point pressure (psia | barsa). Defaults to 1E6, and not used for densities below Pb. A valid value is required for density calculations above Pb
        degf: Reservoir Temperature (deg F | deg C)
        rs: Oil solution gas volume (scf/stb | sm3/sm3)
        rsb: Oil solution gas volume (scf/stb | sm3/sm3) at Pb
        sg_g: Weighted average specific gravity of surface gas (relative to air). If not known, it will be estimated from sg_sp
        sg_sp: Separator gas specific gravity (relative to air).
        sg_o: Stock tank oil specific gravity (SG relative to water).
        bomethod: A string or deno_method Enum class that specifies one of following calculation choices;
                   STAN: Standing Correlation
                   MCAIN: McCain approach, calculating from densities
        denomethod: A string or deno_method Enum class that specifies one of following calculation choices;
                   SWMH: Standing, Witte, McCain-Hill (1995) - Default
        metric: If True, input/output in Eclipse METRIC units (barsa, degC, sm3/sm3). Defaults to False (FIELD)
    """
    if metric:
        p = p * BAR_TO_PSI
        pb = pb * BAR_TO_PSI
        degf = degc_to_degf(degf)
        rs = rs * SM3_PER_SM3_TO_SCF_PER_STB
        rsb = rsb * SM3_PER_SM3_TO_SCF_PER_STB

    validate_pe_inputs(p=p, degf=degf)

    sg_g, sg_sp = check_sgs(sg_g=sg_g, sg_sp=sg_sp)
    denomethod, bomethod = validate_methods(
        ["denomethod", "bomethod"], [denomethod, bomethod]
    )

    def Bo_standing(p, pb, degf, rs, rsb, sg_sp, sg_g, sg_o):
        Bob = (
            _BO_STAN_A
            + _BO_STAN_B * (rs * (sg_g / sg_o) ** 0.5 + _BO_STAN_C * degf) ** _BO_STAN_EXP
        )
        return Bob

    def Bo_mccain(p, pb, degf, rs, rsb, sg_sp, sg_g, sg_o):
        rhor = oil_deno(
            p=p,
            degf=degf,
            rs=rs,
            rsb=rsb,
            sg_g=sg_g,
            sg_sp=sg_sp,
            pb=pb,
            sg_o=sg_o,
            denomethod=denomethod,
        )
        return (sg_o * _BO_MC_WDEN + _BO_MC_RS_COEFF * rs * sg_g) / rhor  # Eq 3.21

    fn_dic = {"STAN": Bo_standing, "MCAIN": Bo_mccain}

    return fn_dic[bomethod.name](p, pb, degf, rs, rsb, sg_sp, sg_g, sg_o)


def oil_viso(p: float, api: float, degf: float, pb: float, rs: float, metric: bool = False) -> float:
    """ Returns Oil Viscosity (cP) with Beggs-Robinson (1975) correlation at saturated pressures
        and Petrosky-Farshad (1995) at undersaturated pressures

        p: Pressure (psia | barsa)
        api: Stock tank oil density (deg API)
        degf: Reservoir Temperature (deg F | deg C)
        pb: Bubble point Pressure (psia | barsa)
        rs: Solution GOR (scf/stb | sm3/sm3)
        metric: If True, input in Eclipse METRIC units (barsa, degC, sm3/sm3). Defaults to False (FIELD)
    """
    if metric:
        p = p * BAR_TO_PSI
        degf = degc_to_degf(degf)
        pb = pb * BAR_TO_PSI
        rs = rs * SM3_PER_SM3_TO_SCF_PER_STB

    validate_pe_inputs(p=p, degf=degf)

    # Correlation validity range warnings (non-blocking)
    if api < 16 or api > 58:
        warnings.warn(f"Beggs-Robinson viscosity: API={api:.1f} outside calibration range [16, 58]", stacklevel=2)
    if degf < 70 or degf > 295:
        warnings.warn(f"Beggs-Robinson viscosity: T={degf:.1f}\u00b0F outside calibration range [70, 295]", stacklevel=2)

    def uo_br(p, api, degf, pb, rs):
        Z = _BR_Z0 + _BR_Z_API * api
        y = 10 ** Z
        X = y * degf ** _BR_T_EXP
        A = _BR_A_MULT * (rs + _BR_A_OFFSET) ** _BR_A_EXP
        B = _BR_B_MULT * (rs + _BR_B_OFFSET) ** _BR_B_EXP

        uod = max(10 ** X - 1, 1e-6)  # Guard against negative/zero dead-oil viscosity
        uor = A * uod ** B  # Eq 3.23c
        return max(uor, 0.01)  # Floor at 0.01 cP

    def uo_pf(p, api, degf, pb, rs):
        uob = uo_br(pb, api, degf, pb, rs)
        if uob <= 0:
            uob = 0.01  # Safety floor
        loguob = np.log(uob)
        A = (
            _PF_POLY[0]
            + _PF_POLY[1] * loguob
            + _PF_POLY[2] * loguob ** 2
            + _PF_POLY[3] * loguob ** 3
        )  # Eq 3.24b
        uor = uob + _PF_P_COEFF * (p - pb) * 10 ** A  # Eq 3.24a
        return max(uor, 0.01)  # Floor at 0.01 cP

    if p <= pb:
        return uo_br(p, api, degf, pb, rs)
    else:
        return uo_pf(p, api, degf, pb, rs)
