================
Material Balance
================

Material balance functions for estimating original hydrocarbons in place from pressure-production history. Gas reservoirs use P/Z linear regression with optional Havlena-Odeh aquifer influx correction. Oil reservoirs use the Havlena-Odeh method with drive index decomposition and optional parameter regression.


pyrestoolbox.matbal.gas_matbal
==============================

.. code-block:: python

    gas_matbal(p, Gp, degf, sg=0.65, co2=0, h2s=0, n2=0, h2=0,
               Wp=None, Bw=1.0, We=None,
               zmethod='DAK', cmethod='PMC', metric=False,
               pvt_table=None, he=0) -> GasMatbalResult

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
     - Cumulative gas production at each pressure survey. Same length as p. Units are user-defined (e.g. Bscf, MMscf) - OGIP will be in the same units. When Wp or We are provided, Gp should be in scf (or sm3 if metric) for dimensional consistency with Bg
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
   * - he
     - float
     - He mole fraction (default 0). Keyword, last in the signature (new in 3.8.3); like h2 it forces BNS
   * - Wp
     - array-like, optional
     - Cumulative water production (STB | sm3). Same length as p
   * - Bw
     - float
     - Water FVF (rb/stb | rm3/sm3, default 1.0). Used when Wp provided
   * - We
     - array-like, optional
     - Cumulative water influx (rcf | rm3). Same length as p. When provided, Havlena-Odeh regression is used for OGIP. Field-unit influx from ``aquifer_influx`` is in rb: multiply by 5.61458 to get rcf
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
   * - pz_curvature
     - float
     - Signed bend of P/Z vs Gp: midpoint offset of a quadratic fit from its chord, as a fraction of the P/Z decline. Positive (concave up) is consistent with pressure support (aquifer influx, a second compartment); negative with abnormal pressure (compaction) or a growing drainage area. Compare it with the scatter before trusting a straight-line OGIP. NaN with fewer than four points

P/Z Example:

.. code-block:: python

    >>> from pyrestoolbox import matbal
    >>> r = matbal.gas_matbal(
    ...     p=[3000, 2700, 2400, 2100, 1800],
    ...     Gp=[0, 5, 12, 22, 35],
    ...     degf=200, sg=0.65
    ... )
    >>> r.ogip
    87.60264634235111
    >>> r.z_initial
    0.9163202021564056
    >>> r.r_squared
    0.9734793606102145
    >>> round(r.pz_curvature, 4)  # concave up: the straight line is a poor model here
    0.1391

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
    53.342120020574896
    >>> r.cole_F_over_Et[2]
    62.641251784427226

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
    67036355572.60765

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
    87.60264634235111
    >>> r.z_initial
    0.9163202021564056


pyrestoolbox.matbal.oil_matbal
==============================

.. code-block:: python

    oil_matbal(p, Np, degf, api=0, sg_sp=0, sg_g=0, pb=0, rsb=0,
               Rp=None, Wp=None, Wi=None, Gi=None, We=None,
               Bw=1.0, m=0, cf=0, sw_i=0, cw=0,
               rsmethod='VELAR', bomethod='MCAIN',
               zmethod='DAK', cmethod='PMC', metric=False,
               pvt_table=None, regress=None) -> OilMatbalResult

Havlena-Odeh oil material balance for OOIP estimation. Computes underground withdrawal (F), oil expansion (Eo), gas cap expansion (Eg), and formation/water compressibility (Efw) terms at each pressure step. When cumulative aquifer influx ``We`` is supplied it is subtracted from the withdrawal first (``F - We = N·[Eo + m*Eg + (1+m)*Efw]``); OOIP is the least-squares slope of ``F - We`` against ``Et = Eo + m*Eg + (1+m)*Efw`` through the origin (Havlena and Odeh 1963), so late points with large Et carry most weight; earlier releases used the mean of the per-point ratios. Drive indices (DDI, SDI, CDI, WDI) use that fitted OOIP at each step.

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
   * - We
     - array-like, optional
     - Cumulative aquifer water influx in reservoir volume (rb | rm3). Subtracted from withdrawal before estimating OOIP and surfaced as the 'WDI' water-drive index. No metric conversion (already a reservoir volume). Default all zeros
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
     - Parameters to regress with bounds: ``{'m': (0, 2), 'cf': (1e-6, 10e-6)}``. Allowed keys: 'm', 'cf', 'cw', 'sw_i'. Least squares on F - We - N Et, with N the slope for each trial set. Raises a ``UserWarning`` when a parameter ends on a bound or is correlated with N above 0.99 (see ``regress_report``)

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
     - Drive index fractions of the withdrawal F at each step, from the fitted OOIP: 'DDI' (depletion), 'SDI' (segregation/gas cap), 'CDI' (compaction/expansion), 'WDI' (water/aquifer influx). Each maps to np.ndarray. Their sum is 1 only where the material balance closes, so a sum away from 1 flags a mismatch (Pletcher 2002). 'WDI' is all-zero when no ``We`` is supplied
   * - p
     - np.ndarray
     - Pressure array used
   * - pvt
     - dict
     - PVT properties at each step: 'Rs', 'Bo', 'Bg' arrays
   * - regressed
     - dict or None
     - Regressed parameter values when ``regress`` is used, else None
   * - ooip_points
     - np.ndarray
     - Per-point OOIP (F - We)/Et, NaN where Et is 0. Against F it is the Campbell plot: flat supports the model, rising points to support the model omits (an aquifer)
   * - regress_report
     - dict or None
     - With ``regress``: 'at_bound' (parameters that ended on a bound), 'max_correlation' (largest absolute correlation between N and the regressed parameters, from a finite-difference Jacobian), 'identifiable' (False when either check fires)

Examples:

.. code-block:: python

    >>> r = matbal.oil_matbal(
    ...     p=[4000, 3500, 3000, 2500],
    ...     Np=[0, 1e6, 3e6, 6e6],
    ...     degf=220, api=35, sg_sp=0.75,
    ...     pb=3500, rsb=500, cf=3e-6, sw_i=0.2, cw=3e-6
    ... )
    >>> r.ooip
    57244701.589534335
    >>> r.drive_indices['DDI'][1]
    0.3149045515580967

These inputs do not satisfy a single OOIP, and the diagnostics show it: the per-point estimates fall (a Campbell plot that declines) and the drive indices sum well below 1 early on:

.. code-block:: python

    >>> [round(x / 1e6, 1) for x in r.ooip_points[1:]]
    [129.2, 63.6, 55.6]
    >>> [round(float(x), 3) for x in sum(r.drive_indices.values())[1:]]
    [0.443, 0.9, 1.03]

Water-drive Example (aquifer influx ``We``):

.. code-block:: python

    >>> r = matbal.oil_matbal(
    ...     p=[4000, 3500, 3000, 2500],
    ...     Np=[0, 1e6, 3e6, 6e6],
    ...     degf=220, api=35, sg_sp=0.75,
    ...     pb=3500, rsb=500, cf=3e-6, sw_i=0.2, cw=3e-6,
    ...     We=[0, 3e5, 1.2e6, 3.0e6]
    ... )
    >>> r.ooip
    35077977.022075504
    >>> r.drive_indices['WDI'][-1]
    0.4055063920332815

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
    UserWarning: regressed m, cf ended on its bound: the data do not constrain it inside the bounds
    UserWarning: OOIP and the regressed parameters are correlated (0.9993): the data fix only a combination of them ...
    >>> round(r.regressed['m'], 3), round(r.regressed['cf'] * 1e6, 2)
    (2.0, 50.0)
    >>> r.ooip
    1427566.2426175773
    >>> r.regress_report['identifiable']
    False

Both parameters ended on their upper bounds and the fitted OOIP is below the 2.1 MMSTB already produced: four points cannot separate N, m and cf here, and the warnings say so. Fix m or cf from independent data and regress the other.

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
    57244701.589534335


.. note:: **Compressibility Guidance**

   The ``cw`` parameter controls how water compressibility is treated in the formation/water expansion term (Efw). Two approaches are common:

   - **Undersaturated (liquid-phase) compressibility**: Use ``cw_usat`` from ``brine.brine_props()``, which is the isothermal compressibility of brine at constant dissolved gas content. Appropriate when the material balance tracks free gas volumes separately, or when the water does not evolve dissolved gas over the depletion range.

   - **Saturated (pseudo) compressibility**: Use ``cw_sat`` from ``brine.brine_props()``, which includes the volume of differentially evolved gas. Appropriate when no separate water-gas tracking is performed. Dissolved gas in water can contribute compressibility as large as or larger than the liquid-phase value alone (Ramey, 1964).

   To compute ``cw`` from brine properties:

   .. code-block:: python

       >>> from pyrestoolbox import brine
       >>> props = brine.brine_props(p=4000, degf=220, wt=3)
       >>> cw_usat, cw_sat = props['cw']  # [undersaturated, saturated]

   Similarly, ``oil.oil_co(co_sat=True)`` returns ``[co_usat, co_sat]`` - the undersaturated value is the liquid-phase compressibility, while the saturated value includes gas evolution effects (Perrine's definition).


pyrestoolbox.matbal.aquifer_influx
==================================

.. code-block:: python

    aquifer_influx(t, p, B, td_per_t, red=None, M=7) -> dict

Cumulative aquifer influx by the Carter-Tracy (1960) method, the influx model behind ECLIPSE AQUCT, using the Van Everdingen-Hurst constant-terminal-rate pD from ``simtools.influence_tables`` (spline in ln tD for pD and its derivative):

We_n = We_{n-1} + (tD_n - tD_{n-1}) (B dp_n - We_{n-1} pD'_n) / (pD_n - tD_{n-1} pD'_n), with dp_n = p_i - p_n.

Units are consistent rather than fixed: We takes the units of B times pressure (rb/psi x psi = rb, or rm3/bar x bar = rm3), so no ``metric`` flag is needed. ``aquifer_constants`` builds B and tD per day for a radial aquifer from its properties in either unit system. Carter-Tracy approximates Van Everdingen-Hurst superposition and tracks it more closely as the steps shorten; 30-day steps or shorter are usually adequate. Feed the result to ``oil_matbal`` as ``We`` directly. ``gas_matbal`` takes ``We`` in rcf in field units, so multiply field-unit influx (rb) by 5.61458 first; metric (rm3) passes straight through.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - t
     - list or array
     - Times, strictly increasing; the first is the initial state (any unit, matching td_per_t)
   * - p
     - list or array
     - Pressure at the reservoir-aquifer boundary at each time (psia | barsa, matching B)
   * - B
     - float
     - Aquifer influx constant (rb/psi | rm3/bar)
   * - td_per_t
     - float
     - Dimensionless time per unit of t
   * - red
     - float
     - Aquifer to reservoir radius ratio re/ro (> 1). None (default) for infinite-acting
   * - M
     - int
     - Laplace inversion accuracy passed to ``simtools.influence_tables``. Default 7

.. list-table:: Returns (dict)
   :widths: 10 15 40
   :header-rows: 1

   * - Key
     - Type
     - Description
   * - 'We'
     - np.ndarray
     - Cumulative influx (units of B x pressure), zero at the first time
   * - 'tD'
     - np.ndarray
     - Dimensionless time at each t
   * - 'pD'
     - np.ndarray
     - Dimensionless pressure at each tD (0 at tD = 0)
   * - 'dpD'
     - np.ndarray
     - dpD/dtD at each tD (NaN at tD = 0)

Example: the 6-month-step case of Ahmed, Reservoir Engineering Handbook (4th ed.), Example 10-9 (infinite aquifer, B = 20.4 bbl/psi, tD = 0.9888 t days). The book prints 12,266, 42,546, 104,406 and 202,477 bbl; it rounds pD from the long-time approximation, and its second step carries an arithmetic slip (its own pD values give 44,546, not 42,546), which propagates. On its monthly-step table the two agree to 0.2-0.7% from month 6.

.. code-block:: python

    >>> r = matbal.aquifer_influx(t=[0, 182.5, 365, 547.5, 730], p=[2500, 2490, 2472, 2444, 2408],
    ...                           B=20.4, td_per_t=0.9888)
    >>> [round(x) for x in r['We']]
    [0, 12223, 44429, 106088, 203972]


pyrestoolbox.matbal.aquifer_constants
=====================================

.. code-block:: python

    aquifer_constants(phi, ct, h, ro, k, mu, theta=360, metric=False) -> dict

Radial aquifer influx constant and dimensionless time per day (Van Everdingen and Hurst, 1949): B = 2 pi phi ct h ro^2 (theta/360) and tD = k t / (phi mu ct ro^2), with the unit factors of the chosen system. In field units these are the familiar B = 1.119 phi ct h ro^2 (theta/360) rb/psi and tD per day = 0.0063283 k / (phi mu ct ro^2); in metric the Darcy constant is 0.0085270 (bar, m). The returned keys match ``aquifer_influx``, so the dict unpacks straight into it with t in days.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - phi
     - float
     - Aquifer porosity (fraction)
   * - ct
     - float
     - Aquifer total compressibility, water plus rock (1/psi | 1/barsa)
   * - h
     - float
     - Aquifer thickness (ft | m)
   * - ro
     - float
     - Reservoir (inner aquifer) radius (ft | m)
   * - k
     - float
     - Aquifer permeability (mD)
   * - mu
     - float
     - Aquifer water viscosity (cP)
   * - theta
     - float
     - Encroachment angle in degrees, 0 < theta <= 360. Default 360
   * - metric
     - bool
     - If True, inputs and B in Eclipse METRIC units. Default False

.. list-table:: Returns (dict)
   :widths: 10 15 40
   :header-rows: 1

   * - Key
     - Type
     - Description
   * - 'B'
     - float
     - Aquifer influx constant (rb/psi | rm3/bar)
   * - 'td_per_t'
     - float
     - Dimensionless time per day

Example: a half-circle edge aquifer, then Carter-Tracy influx over two years. The metric call is the same aquifer (B in rm3/bar; tD per day unchanged).

.. code-block:: python

    >>> from pyrestoolbox import matbal
    >>> c = matbal.aquifer_constants(phi=0.25, ct=7e-6, h=100, ro=3000, k=50, mu=0.5, theta=180)
    >>> round(c['B'], 2), round(c['td_per_t'], 6)
    (881.28, 0.04018)
    >>> r = matbal.aquifer_influx(t=[0, 182.5, 365, 547.5, 730], p=[2500, 2490, 2472, 2444, 2408], **c)
    >>> [round(x) for x in r['We']]
    [0, 42544, 149971, 349715, 659711]
    >>> m = matbal.aquifer_constants(phi=0.25, ct=1.015e-4, h=30.48, ro=914.4, k=50, mu=0.5, theta=180, metric=True)
    >>> round(m['B'], 2), round(m['td_per_t'], 6)
    (2031.62, 0.04019)

Class Objects
=============

.. list-table::
   :widths: 15 40
   :header-rows: 1

   * - Class
     - Description
   * - GasMatbalResult
     - Result from ``gas_matbal()``. See `Returns table above <#pyrestoolbox-matbal-gas-matbal>`__ for full attribute details
   * - OilMatbalResult
     - Result from ``oil_matbal()``. See `Returns table above <#pyrestoolbox-matbal-oil-matbal>`__ for full attribute details
