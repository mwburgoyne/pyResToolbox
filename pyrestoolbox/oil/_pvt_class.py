"""OilPVT convenience class."""

import numpy as np

from pyrestoolbox.constants import (
    BAR_TO_PSI, degc_to_degf,
    SM3_PER_SM3_TO_SCF_PER_STB, SCF_PER_STB_TO_SM3_PER_SM3,
    LBCUFT_TO_KGM3,
)
from pyrestoolbox.validate import validate_methods

from ._utils import check_sgs, oil_sg
from ._correlations import oil_rs, oil_bo, oil_deno, oil_viso
from ._harmonize import oil_harmonize


class OilPVT:
    """ Oil PVT wrapper that stores oil characterization and method choices.
        Exposes rs(), bo(), density(), viscosity() methods delegating to oil_rs, oil_bo, oil_deno, oil_viso.

        api: Stock tank oil density (deg API)
        sg_sp: Separator gas specific gravity (relative to air)
        pb: Bubble point pressure (psia | barsa)
        rsb: Solution GOR at Pb (scf/stb | sm3/sm3). If 0 and degf > 0, calculated from pb via oil_harmonize()
        sg_g: Weighted average specific gravity of surface gas (relative to air). Estimated from sg_sp if not provided
        degf: Reservoir temperature (deg F | deg C). If > 0, triggers auto-harmonization via oil_harmonize()
        uo_target: Target oil viscosity (cP) at p_uo. Used with degf for auto-harmonization. Default 0 (no tuning)
        p_uo: Pressure at which uo_target was measured (psia | barsa). Required if uo_target > 0
        vis_frac: Viscosity scaling factor. All viscosity outputs are multiplied by this value. Default 1.0
        rsb_frac: Rs scaling factor from oil_harmonize(). Rs = oil_rs(rsb=rsb/rsb_frac) * rsb_frac. Default 1.0
        rsmethod: Method for Rs calculation. Defaults to 'VELAR'
        pbmethod: Method for Pb calculation. Defaults to 'VALMC'
        bomethod: Method for Bo calculation. Defaults to 'MCAIN'
        metric: If True, constructor inputs (pb, rsb) and method inputs/outputs use Eclipse METRIC units. Defaults to False (FIELD)
    """
    def __init__(self, api: float, sg_sp: float, pb: float, rsb: float = 0,
                 sg_g: float = 0, degf: float = 0,
                 uo_target: float = 0, p_uo: float = 0,
                 vis_frac: float = 1.0, rsb_frac: float = 1.0,
                 rsmethod: str = 'VELAR', pbmethod: str = 'VALMC',
                 bomethod: str = 'MCAIN',
                 metric: bool = False):
        self.api = api
        self.sg_sp = sg_sp
        self.metric = metric
        # Auto-harmonize when degf is provided
        if degf > 0:
            pb_h, rsb_h, rsb_frac_h, vis_frac_h = oil_harmonize(
                pb=pb, rsb=rsb, degf=degf, api=api, sg_sp=sg_sp, sg_g=sg_g,
                uo_target=uo_target, p_uo=p_uo,
                rsmethod=rsmethod, pbmethod=pbmethod, metric=metric)
            pb, rsb = pb_h, rsb_h
            rsb_frac, vis_frac = rsb_frac_h, vis_frac_h
        elif rsb <= 0:
            raise ValueError("Either rsb or degf must be specified")
        if vis_frac <= 0:
            raise ValueError("vis_frac must be positive")
        if rsb_frac <= 0:
            raise ValueError("rsb_frac must be positive")
        self.vis_frac = vis_frac
        self.rsb_frac = rsb_frac
        if metric:
            self.pb = pb * BAR_TO_PSI
            self.rsb = rsb * SM3_PER_SM3_TO_SCF_PER_STB
        else:
            self.pb = pb
            self.rsb = rsb
        self.sg_o = oil_sg(api)
        self.sg_g, self.sg_sp = check_sgs(sg_g=sg_g, sg_sp=sg_sp)
        self.rsmethod = validate_methods(["rsmethod"], [rsmethod])
        self.pbmethod = validate_methods(["pbmethod"], [pbmethod])
        self.bomethod = validate_methods(["bomethod"], [bomethod])

    @classmethod
    def from_harmonize(cls, degf: float, api: float, sg_sp: float = 0,
                       sg_g: float = 0, pb: float = 0, rsb: float = 0,
                       uo_target: float = 0, p_uo: float = 0,
                       rsmethod: str = 'VELAR', pbmethod: str = 'VELAR',
                       bomethod: str = 'MCAIN',
                       metric: bool = False) -> 'OilPVT':
        """Deprecated: use OilPVT(degf=...) directly.

        Convenience constructor that calls oil_harmonize() internally.
        """
        return cls(api=api, sg_sp=sg_sp, pb=pb, rsb=rsb, sg_g=sg_g, degf=degf,
                   uo_target=uo_target, p_uo=p_uo,
                   rsmethod=rsmethod, pbmethod=pbmethod, bomethod=bomethod,
                   metric=metric)

    def _rs_field(self, p_field, degf_field):
        """Internal: compute Rs in oilfield units applying rsb_frac scaling."""
        raw = oil_rs(api=self.api, degf=degf_field, sg_sp=self.sg_sp, p=p_field,
                     pb=self.pb, rsb=self.rsb / self.rsb_frac,
                     rsmethod=self.rsmethod, pbmethod=self.pbmethod)
        return raw * self.rsb_frac

    @staticmethod
    def _is_array(x):
        return isinstance(x, (list, tuple, np.ndarray))

    def rs(self, p, degf: float):
        """ Returns solution GOR (scf/stb | sm3/sm3) at pressure p (psia | barsa) and temperature degf (deg F | deg C).
            p can be a scalar, list, or array. degf must be scalar. """
        if self._is_array(p):
            return np.array([self.rs(pi, degf) for pi in p])
        if self.metric:
            p = p * BAR_TO_PSI
            degf = degc_to_degf(degf)
        result = self._rs_field(p, degf)
        if self.metric:
            return result * SCF_PER_STB_TO_SM3_PER_SM3
        return result

    def bo(self, p, degf: float, rs=None):
        """ Returns oil FVF (rb/stb | rm3/sm3) at pressure p (psia | barsa) and temperature degf (deg F | deg C).
            p can be a scalar, list, or array. degf must be scalar. """
        if self._is_array(p):
            rs_arr = rs if rs is None else list(rs)
            return np.array([self.bo(pi, degf, rs=None if rs is None else rs_arr[i])
                             for i, pi in enumerate(p)])
        if self.metric:
            p_field = p * BAR_TO_PSI
            degf_field = degc_to_degf(degf)
        else:
            p_field = p
            degf_field = degf
        if rs is None:
            rs_field = self._rs_field(p_field, degf_field)
        else:
            rs_field = rs * SM3_PER_SM3_TO_SCF_PER_STB if self.metric else rs
        return oil_bo(p=p_field, pb=self.pb, degf=degf_field, rs=rs_field, rsb=self.rsb,
                      sg_o=self.sg_o, sg_g=self.sg_g, sg_sp=self.sg_sp,
                      bomethod=self.bomethod)

    def density(self, p, degf: float, rs=None):
        """ Returns live oil density (lb/cuft | kg/m3) at pressure p (psia | barsa) and temperature degf (deg F | deg C).
            p can be a scalar, list, or array. degf must be scalar. """
        if self._is_array(p):
            rs_arr = rs if rs is None else list(rs)
            return np.array([self.density(pi, degf, rs=None if rs is None else rs_arr[i])
                             for i, pi in enumerate(p)])
        if self.metric:
            p_field = p * BAR_TO_PSI
            degf_field = degc_to_degf(degf)
        else:
            p_field = p
            degf_field = degf
        if rs is None:
            rs_field = self._rs_field(p_field, degf_field)
        else:
            rs_field = rs * SM3_PER_SM3_TO_SCF_PER_STB if self.metric else rs
        result = oil_deno(p=p_field, degf=degf_field, rs=rs_field, rsb=self.rsb,
                        sg_g=self.sg_g, sg_sp=self.sg_sp, pb=self.pb,
                        sg_o=self.sg_o)
        if self.metric:
            return result * LBCUFT_TO_KGM3
        return result

    def viscosity(self, p, degf: float, rs=None):
        """ Returns oil viscosity (cP) at pressure p (psia | barsa) and temperature degf (deg F | deg C).
            p can be a scalar, list, or array. degf must be scalar. """
        if self._is_array(p):
            rs_arr = rs if rs is None else list(rs)
            return np.array([self.viscosity(pi, degf, rs=None if rs is None else rs_arr[i])
                             for i, pi in enumerate(p)])
        if self.metric:
            p_field = p * BAR_TO_PSI
            degf_field = degc_to_degf(degf)
        else:
            p_field = p
            degf_field = degf
        if rs is None:
            rs_field = self._rs_field(p_field, degf_field)
        else:
            rs_field = rs * SM3_PER_SM3_TO_SCF_PER_STB if self.metric else rs
        return oil_viso(p=p_field, api=self.api, degf=degf_field, pb=self.pb, rs=rs_field) * self.vis_frac
