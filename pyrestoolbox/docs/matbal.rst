===================================
Material Balance
===================================

Material balance functions for estimating original hydrocarbons in place from pressure-production history. Gas reservoirs use P/Z linear regression with optional Havlena-Odeh aquifer influx correction. Oil reservoirs use the Havlena-Odeh method with drive index decomposition and optional parameter regression.


pyrestoolbox.matbal.gas_matbal
======================

.. code-block:: python

    gas_matbal(p, Gp, degf, sg=0.65, co2=0, h2s=0, n2=0, h2=0,
               Wp=None, Bw=1.0, We=None,
               zmethod='DAK', cmethod='PMC', metric=False,
               pvt_table=None) -> GasMatbalResult

P/Z gas material balance for OGIP estimation. Performs linear regression of P/Z vs cumulative gas production to determine original gas in place (OGIP = -intercept/slope). Optionally computes Cole plot diagnostics (F/Et vs Gp) and Havlena-Odeh regression when cumulative water influx (We) is provided.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - array-like
     - Reservoir pressures at each survey (psia | barsa). First value is initial pressure
   * - Gp
     - array-like
     - Cumulative gas production at each pressure survey. Same length as p. Units are user-defined (e.g. Bscf, MMscf) — OGIP will be in the same units. When Wp or We are provided, Gp should be in scf (or sm3 if metric) for dimensional consistency with Bg
   * - degf
     - float
     - Reservoir temperature (deg F | deg C)
   * - sg
     - float
     - Gas specific gravity (default 0.65)
   * - co2
     - float
     - CO2 mole fraction (default 0)
   * - h2s
     - float
     - H2S mole fraction (default 0)
   * - n2
     - float
     - N2 mole fraction (default 0)
   * - h2
     - float
     - H2 mole fraction (default 0)
   * - Wp
     - array-like, optional
     - Cumulative water production (STB | sm3). Same length as p
   * - Bw
     - float
     - Water FVF (rb/stb | rm3/sm3, default 1.0). Used when Wp provided
   * - We
     - array-like, optional
     - Cumulative water influx (rcf | rm3). Same length as p. When provided, Havlena-Odeh regression is used for OGIP
   * - zmethod
     - str
     - Z-factor method (default 'DAK')
   * - cmethod
     - str
     - Critical property method (default 'PMC')
   * - metric
     - bool
     - If True, p in barsa and degf in deg C (default False)
   * - pvt_table
     - dict, optional
     - Tabulated gas PVT. Either ``{'p': [...], 'Z': [...]}`` or ``{'p': [...], 'Bg': [...]}``. Providing both 'Z' and 'Bg' raises ValueError

.. list-table:: Returns (GasMatbalResult)
   :widths: 10 15 40
   :header-rows: 1

   * - Attribute
     - Type
     - Description
   * - ogip
     - float
     - Original gas in place (same volume units as Gp input)
   * - pz
     - np.ndarray
     - P/Z values at each pressure step
   * - gp
     - np.ndarray
     - Cumulative gas production at each step
   * - slope
     - float
     - Slope of the P/Z vs Gp regression line
   * - intercept
     - float
     - Intercept of the P/Z vs Gp regression line
   * - r_squared
     - float
     - Coefficient of determination of the linear fit
   * - p_initial
     - float
     - Initial reservoir pressure
   * - z_initial
     - float
     - Z-factor at initial pressure
   * - bg
     - np.ndarray
     - Gas FVF (rcf/scf | rm3/sm3) at each pressure step
   * - F
     - np.ndarray
     - Underground withdrawal at each step
   * - Et
     - np.ndarray
     - Gas expansion term (Bg - Bgi) at each step
   * - cole_F_over_Et
     - np.ndarray
     - Cole plot diagnostic: F/Et at each step (NaN at index 0)
   * - method
     - str
     - OGIP determination method: 'pz' or 'havlena_odeh'

P/Z Example:

.. code-block:: python

    >>> from pyrestoolbox import matbal
    >>> r = matbal.gas_matbal(
    ...     p=[3000, 2700, 2400, 2100, 1800],
    ...     Gp=[0, 5, 12, 22, 35],
    ...     degf=200, sg=0.65
    ... )
    >>> r.ogip
    87.602774253829
    >>> r.z_initial
    0.9163208839373836
    >>> r.r_squared
    0.9734794008096929

Cole Plot Example (volumetric diagnostic):

.. code-block:: python

    >>> from pyrestoolbox import matbal
    >>> r = matbal.gas_matbal(
    ...     p=[3000, 2700, 2400, 2100, 1800],
    ...     Gp=[0, 5, 12, 22, 35],
    ...     degf=200, sg=0.65
    ... )
    >>> r.method
    'pz'
    >>> r.cole_F_over_Et[1]
    53.34214259832056
    >>> r.cole_F_over_Et[2]
    62.64134359331348

Havlena-Odeh Example (with aquifer influx):

.. code-block:: python

    >>> from pyrestoolbox import matbal
    >>> r = matbal.gas_matbal(
    ...     p=[3000, 2700, 2400, 2100, 1800],
    ...     Gp=[0, 5e9, 12e9, 22e9, 35e9],
    ...     degf=200, sg=0.65,
    ...     We=[0, 5e6, 15e6, 35e6, 60e6]
    ... )
    >>> r.method
    'havlena_odeh'
    >>> r.ogip
    67036445117.070206

Tabulated PVT Example (Z-factor table):

.. code-block:: python

    >>> from pyrestoolbox import matbal, gas
    >>> p_table = list(range(1500, 3501, 50))
    >>> Z_table = [gas.gas_z(pi, 0.65, 200) for pi in p_table]
    >>> r = matbal.gas_matbal(
    ...     p=[3000, 2700, 2400, 2100, 1800],
    ...     Gp=[0, 5, 12, 22, 35],
    ...     degf=200, sg=0.65,
    ...     pvt_table={'p': p_table, 'Z': Z_table}
    ... )
    >>> r.ogip
    87.602774253829
    >>> r.z_initial
    0.9163208839373836


pyrestoolbox.matbal.oil_matbal
======================

.. code-block:: python

    oil_matbal(p, Np, degf, api=0, sg_sp=0, sg_g=0, pb=0, rsb=0,
               Rp=None, Wp=None, Wi=None, Gi=None,
               Bw=1.0, m=0, cf=0, sw_i=0, cw=0,
               rsmethod='VELAR', bomethod='MCAIN',
               zmethod='DAK', cmethod='PMC', metric=False,
               pvt_table=None, regress=None) -> OilMatbalResult

Havlena-Odeh oil material balance for OOIP estimation. Computes underground withdrawal (F), oil expansion (Eo), gas cap expansion (Eg), and formation/water compressibility (Efw) terms at each pressure step. OOIP is estimated as the mean of F/(Eo + m*Eg + (1+m)*Efw) across valid steps. Drive indices (DDI, SDI, CDI) are computed at each step.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - array-like
     - Reservoir pressures at each survey (psia | barsa). First value is initial pressure
   * - Np
     - array-like
     - Cumulative oil production (STB | sm3) at each step. First entry is typically 0
   * - degf
     - float
     - Reservoir temperature (deg F | deg C)
   * - api
     - float
     - Stock tank oil API gravity. Required when pvt_table is not provided
   * - sg_sp
     - float
     - Separator gas specific gravity. Required when pvt_table is not provided
   * - sg_g
     - float
     - Weighted average surface gas specific gravity (default 0 = use sg_sp)
   * - pb
     - float
     - Bubble point pressure (psia | barsa). If 0, calculated from rsb
   * - rsb
     - float
     - Solution GOR at bubble point (scf/stb | sm3/sm3). If 0, calculated from pb
   * - Rp
     - array-like, optional
     - Cumulative producing GOR (scf/stb | sm3/sm3) at each step. Default = Rs
   * - Wp
     - array-like, optional
     - Cumulative water production (STB | sm3). Default all zeros
   * - Wi
     - array-like, optional
     - Cumulative water injection (STB | sm3). Default all zeros
   * - Gi
     - array-like, optional
     - Cumulative gas injection (scf | sm3). Default all zeros
   * - Bw
     - float
     - Water FVF (rb/stb | rm3/sm3, default 1.0)
   * - m
     - float
     - Gas cap size ratio (default 0)
   * - cf
     - float
     - Formation compressibility (1/psi | 1/bar, default 0)
   * - sw_i
     - float
     - Initial water saturation (fraction, default 0)
   * - cw
     - float
     - Water compressibility (1/psi | 1/bar, default 0). See compressibility guidance below
   * - rsmethod
     - str
     - Solution GOR method (default 'VELAR')
   * - bomethod
     - str
     - Oil FVF method (default 'MCAIN')
   * - zmethod
     - str
     - Z-factor method for Bg (default 'DAK')
   * - cmethod
     - str
     - Critical property method for Bg (default 'PMC')
   * - metric
     - bool
     - If True, inputs/outputs in Eclipse METRIC units (default False)
   * - pvt_table
     - dict, optional
     - Tabulated oil PVT: ``{'p': [...], 'Rs': [...], 'Bo': [...], 'Bg': [...]}``. Units follow the metric flag. When provided, api and sg_sp are not required
   * - regress
     - dict, optional
     - Parameters to regress with bounds: ``{'m': (0, 2), 'cf': (1e-6, 10e-6)}``. Allowed keys: 'm', 'cf', 'cw', 'sw_i'. Minimizes coefficient of variation of OOIP estimates

.. list-table:: Returns (OilMatbalResult)
   :widths: 10 15 40
   :header-rows: 1

   * - Attribute
     - Type
     - Description
   * - ooip
     - float
     - Original oil in place (STB | sm3)
   * - F
     - np.ndarray
     - Underground withdrawal at each step
   * - Eo
     - np.ndarray
     - Oil expansion term at each step
   * - Eg
     - np.ndarray
     - Gas cap expansion term at each step
   * - Efw
     - np.ndarray
     - Formation and water compressibility term at each step
   * - drive_indices
     - dict
     - Drive index fractions at each step: 'DDI' (depletion), 'SDI' (segregation/gas cap), 'CDI' (compaction/water). Each maps to np.ndarray
   * - p
     - np.ndarray
     - Pressure array used
   * - pvt
     - dict
     - PVT properties at each step: 'Rs', 'Bo', 'Bg' arrays
   * - regressed
     - dict or None
     - Optimized parameter values when ``regress`` is used, else None

Examples:

.. code-block:: python

    >>> r = matbal.oil_matbal(
    ...     p=[4000, 3500, 3000, 2500],
    ...     Np=[0, 1e6, 3e6, 6e6],
    ...     degf=220, api=35, sg_sp=0.75,
    ...     pb=3500, rsb=500, cf=3e-6, sw_i=0.2, cw=3e-6
    ... )
    >>> r.ooip
    82793519.84914012
    >>> r.drive_indices['DDI'][1]
    0.7108509458427899

Regression Example (optimizing m and cf):

.. code-block:: python

    >>> r = matbal.oil_matbal(
    ...     p=[4000, 3500, 3000, 2500],
    ...     Np=[0, 500000, 1200000, 2100000],
    ...     degf=200, api=35, sg_sp=0.65,
    ...     pb=3500, rsb=500,
    ...     m=0, cf=0, sw_i=0.2, cw=3e-6,
    ...     regress={'m': (0, 2), 'cf': (1e-7, 50e-6)}
    ... )
    >>> r.regressed['m']
    2.0
    >>> r.regressed['cf']
    5e-05
    >>> r.ooip
    1468306.5302524716

Tabulated PVT Example:

.. code-block:: python

    >>> from pyrestoolbox import matbal, oil, gas
    >>> from pyrestoolbox.constants import CUFTperBBL
    >>> api, sg_sp, pb, rsb, degf = 35, 0.75, 3500, 500, 220
    >>> sg_o = 141.5 / (api + 131.5)
    >>> p_table = list(range(2000, 5001, 100))
    >>> Rs_t = [oil.oil_rs(api, degf, sg_sp, pi, pb=pb, rsb=rsb) for pi in p_table]
    >>> Bo_t = [oil.oil_bo(pi, pb, degf, rs, rsb, sg_o, sg_sp=sg_sp) for pi, rs in zip(p_table, Rs_t)]
    >>> Bg_t = [gas.gas_bg(pi, sg_sp, degf) / CUFTperBBL for pi in p_table]
    >>> r = matbal.oil_matbal(
    ...     p=[4000, 3500, 3000, 2500], Np=[0, 1e6, 3e6, 6e6], degf=degf,
    ...     pb=pb, rsb=rsb, cf=3e-6, sw_i=0.2, cw=3e-6,
    ...     pvt_table={'p': p_table, 'Rs': Rs_t, 'Bo': Bo_t, 'Bg': Bg_t}
    ... )
    >>> r.ooip
    82793519.84914012


.. note:: **Compressibility Guidance**

   The ``cw`` parameter controls how water compressibility is treated in the formation/water expansion term (Efw). Two approaches are common:

   - **Undersaturated (liquid-phase) compressibility**: Use ``cw_usat`` from ``brine.brine_props()``, which is the isothermal compressibility of brine at constant dissolved gas content. Appropriate when the material balance tracks free gas volumes separately, or when the water does not evolve dissolved gas over the depletion range.

   - **Saturated (pseudo) compressibility**: Use ``cw_sat`` from ``brine.brine_props()``, which includes the volume of differentially evolved gas. Appropriate when no separate water-gas tracking is performed. Dissolved gas in water can contribute compressibility as large as or larger than the liquid-phase value alone (Ramey, 1964).

   To compute ``cw`` from brine properties:

   .. code-block:: python

       >>> from pyrestoolbox import brine
       >>> props = brine.brine_props(p=4000, degf=220, wt=3)
       >>> cw_usat, cw_sat = props['cw']  # [undersaturated, saturated]

   Similarly, ``oil.oil_co(co_sat=True)`` returns ``[co_usat, co_sat]`` — the undersaturated value is the liquid-phase compressibility, while the saturated value includes gas evolution effects (Perrine's definition).


Class Objects
======================

.. list-table::
   :widths: 15 40
   :header-rows: 1

   * - Class
     - Description
   * - GasMatbalResult
     - Result from ``gas_matbal()``. See `Returns table above <#pyrestoolbox-matbal-gas-matbal>`__ for full attribute details
   * - OilMatbalResult
     - Result from ``oil_matbal()``. See `Returns table above <#pyrestoolbox-matbal-oil-matbal>`__ for full attribute details
