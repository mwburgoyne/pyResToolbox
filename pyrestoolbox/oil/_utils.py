"""Utility functions for the oil module."""

import numpy as np
from typing import Tuple

from ._constants import _API_NUMER, _API_DENOM
from pyrestoolbox.constants import PSI_TO_BAR


def get_real_part(value):
    if isinstance(value, complex):
        return value.real
    else:
        return value

def oil_sg(api_value: float) -> float:
    """ Returns oil specific gravity given API value of oil
        api_value: API value
    """
    return _API_NUMER / (api_value + _API_DENOM)

def oil_api(sg_value: float) -> float:
    """ Returns oil API given specific gravity value of oil
        sg_value: Specific gravity (relative to water)
    """
    return _API_NUMER / sg_value - _API_DENOM

def check_sgs(
    sg_g: float,
    sg_sp: float,
    rst: float = 5,
    rsp: float = 1000,
    sg_st: float = 1.15,
) -> Tuple:
    """ Function used to impute sg_g or sg_sp when one or the other is zero
        sg_g: The weighted average surface-gas specific gravity (sep gas + gas evolved from liquid after separation)
        sg_sp: Separator specific gas gravity

        Optional
            rst: Post separation GOR (scf/stb)
            rsp: Separator GOR (scf/stb)
            sg_st: Gas sg evolved from post separartor liquid (rel air)

    """


    if sg_g > 0 and sg_sp > 0:
        return sg_g, sg_sp
    if sg_g <= 0 and sg_sp > 0:  # Estimate sg_g from sg_sp
        sg_g = (sg_sp * rsp + sg_st * rst) / (rsp + rst)
    if sg_g > 0 and sg_sp <= 0:  # Estimate sg_sp from sg_g
        sg_sp = (sg_g * (rsp + rst) - (sg_st * rst)) / rsp
    if sg_g <= 0 and sg_sp <= 0:
        raise ValueError("At least one of sg_g or sg_sp must be positive")
    if sg_g < sg_sp:
        sg_sp = sg_g
    return (sg_g, sg_sp)

def oil_ja_sg(mw: float, ja: float) -> float:
    """ Returns liquid hydrocarbon specific gravity using Jacoby Aromaticity Factor relationship
        mw: Molecular weight of the liquid (g/gmole or lb/lb-mol)
        Ja: Varies between 0 (Paraffins) - 1 (Aromatic)
    """
    ja = min(1, ja)
    ja = max(0, ja)
    return 0.8468 - 15.8 / mw + ja * (0.2456 - 1.77 / mw)


def oil_twu_props(
    mw: float, ja: float = 0, sg: float = 0, damp: float = 1, metric: bool = False
) -> Tuple:
    """ Returns Tuple of sg, tb (degR | K), tc (DegR | K), pc (psia | barsa), Vc (ft3/lb-mol) using method from Twu (1984) correlations for petroleum liquids
        Modified with damping factor proposed by A. Zick between 0 (paraffin) and 1 (original Twu)
        Returns sg, tb (R | K), tc (R | K), pc (psia | barsa), vc (ft3/lbmol)

        mw: Molecular weight of the liquid hydrocarbon (g/g.mol / lb/lb.mol)
        ja: jacoby Aromaticity Factor relationship. Varies between 0 (Paraffins) - 1 (Aromatic). Defaults to zero if undefined
        sg: Specific gravity of the liquid (fraction relative to water density). Will use jacoby method to estimate sg from mw if undefined.
        damp: damping factor proposed by A. Zick between 0 (paraffin) and 1 (original Twu). Defaults to 1
        metric: If True, output Tc/Tb in Kelvin and Pc in barsa. Defaults to False (FIELD: Rankine, psia)
        Unless otherwise mentioned, all Twu equation references are from Whitson Monograph
    """
    if sg == 0:
        sg = oil_ja_sg(mw, ja)  # Use jacoby relationship to estimate sg if not specified
        #print('sg', sg)

    # Estimate boiling point
    # Return boiling point (deg R) and Paraffin properties
    def Twu_tb(mw, sg, damp=1):
        Mp_guess = mw  # Guess for paraffinic mw
        tb, tcp, pcp, vcp, sgp = paraffin_props(Mp_guess)
        d_err = mw - M(tb, sgp, sg, Mp_guess, damp)
        n_iter = 0
        while abs(d_err / mw) > 0.0001:
            n_iter += 1
            Mp_guess += d_err
            tb, tcp, pcp, vcp, sgp = paraffin_props(Mp_guess)
            d_err = mw - M(tb, sgp, sg, Mp_guess, damp)
            if n_iter > 100:
                raise RuntimeError(f"Twu algorithm did not converge for mw={mw}, ja={ja}, sg={sg}, damp={damp}")
        return tb, Mp_guess, tcp, pcp, vcp, sgp

    # Return mw from modified Eq 5.78 to take into account damping
    def M(tb, sgp, sg, Mp, damp):
        absx = abs(0.012342 - 0.328086 / tb ** 0.5)  # Just above Eq 5.78
        dsgM = (
            np.exp(5 * (sgp - sg)) - 1
        )  # Modified Eq 5.78 to take into account damping
        fm = dsgM * (
            absx + (-0.0175691 + 0.193168 / tb ** 0.5) * dsgM
        )  # Just above Eq 5.78
        M = np.exp(
            np.log(Mp) * (1 + 8 * damp * fm / (1 - 2 * fm) ** 2)
        )  # Modified Eq 5.78 to take into account damping
        return M

    def Twu_tc(tb, sgp, sg):
        tcp = (
            tb
            * (
                0.533272
                + 0.000191017 * tb
                + 0.0000000779681 * tb ** 2
                - 2.84376e-11 * tb ** 3
                + 95.9468 / (0.01 * tb) ** 13
            )
            ** -1
        )  # Eq 5.67
        dsgT = np.exp(5 * (sgp - sg)) - 1  # Eq 5.75
        ft = dsgT * (
            (-0.362456 / tb ** 0.5)
            + (0.0398285 - (0.948125 / tb ** 0.5)) * dsgT
        )  # Eq 5.75
        tc = tcp * ((1 + 2 * ft) / (1 - 2 * ft)) ** 2  # Eq 5.75
        return tc

    def Twu_vc(tb, tcp, sg, sgp):
        alpha = 1 - tb / tcp  # Eq 5.72
        vcp = (
            1
            - (
                0.419869
                - 0.505839 * alpha
                - 1.56436 * alpha ** 3
                - 9481.7 * alpha ** 14
            )
        ) ** -8  # Eq 5.69
        dsgV = np.exp(4 * (sgp ** 2 - sg ** 2)) - 1  # Eq 5.76
        f_v = dsgV * (
            (0.46659 / tb ** 0.5) + (-0.182421 + (3.01721 / tb ** 0.5)) * dsgV
        )  # Eq 5.76
        vc = vcp * ((1 + 2 * f_v) / (1 - 2 * f_v)) ** 2  # Eq 5.76
        return vc

    def Twu_pc(tb, sgp, sg, pcp, tc, tcp, vc, vcp):
        dsgp = np.exp(0.5 * (sgp - sg)) - 1  # Eq 5.77
        fp = dsgp * (
            (2.53262 - 46.1955 / tb ** 0.5 - 0.00127885 * tb)
            + (-11.4277 + 252.14 / tb ** 0.5 + 0.00230533 * tb) * dsgp
        )  # Eq 5.77
        pc = (
            pcp * (tc / tcp) * (vcp / vc) * ((1 + 2 * fp) / (1 - 2 * fp)) ** 2
        )  # Eq 5.77
        return pc

    def paraffin_props(Mp):
        theta = np.log(Mp)  # Eq 5.73
        tb = (
            np.exp(
                5.71419
                + 2.71579 * theta
                - 0.28659 * theta ** 2
                - 39.8544 / theta
                - 0.122488 / theta ** 2
            )
            - 24.7522 * theta
            + 35.3155 * theta ** 2
        )  # Eq 5.71
        tcp = (
            tb
            * (
                0.533272
                + 0.000191017 * tb
                + 0.0000000779681 * tb ** 2
                - 2.84376e-11 * tb ** 3
                + 95.9468 / (0.01 * tb) ** 13
            )
            ** -1
        )  # Eq. 5.67
        alpha = 1 - tb / tcp  # Eq 5.72
        pcp = (
            3.83354
            + 1.19629 * alpha ** 0.5
            + 34.8888 * alpha
            + 36.1952 * alpha ** 2
            + 104.193 * alpha ** 4
        ) ** 2  # Eq 5.68
        vcp = (
            1
            - (
                0.419869
                - 0.505839 * alpha
                - 1.56436 * alpha ** 3
                - 9481.7 * alpha ** 14
            )
        ) ** -8  # Eq 5.69
        sgp = (
            0.843593
            - 0.128624 * alpha
            - 3.36159 * alpha ** 3
            - 13749.5 * alpha ** 12
        )  # Eq 5.70
        return tb, tcp, pcp, vcp, sgp

    tb, Mp, tcp, pcp, vcp, sgp = Twu_tb(mw, sg, damp)
    #print(tb, Mp, tcp, pcp, vcp, sgp)
    tc = Twu_tc(tb, sgp, sg)

    vc = Twu_vc(tb, tcp, sg, sgp)
    pc = Twu_pc(tb, sgp, sg, pcp, tc, tcp, vc, vcp)
    #print(tc, pc, vc)
    if metric:
        tb = tb / 1.8  # Rankine -> Kelvin
        tc = tc / 1.8  # Rankine -> Kelvin
        pc = pc * PSI_TO_BAR  # psia -> barsa
    return (sg, tb, tc, pc, vc)
