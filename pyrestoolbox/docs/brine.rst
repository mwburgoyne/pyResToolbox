===================================
Brine PVT with differing degrees of methane, CO2, or multicomponent gas saturation
===================================

Which one should I use?
----------------------

Three brine property models are available. They are **not** interchangeable
wrappers around one calculation: they carry three different *solubility* models,
and each is the better choice somewhere. The density and viscosity legs
downstream of the solubility are now shared, so the models agree on gas-free
brine density to 0.03 to 0.05% and on CH4-saturated viscosity to 0.33%; what
distinguishes them is how much gas they put into the brine.

.. list-table:: Choosing an entry point
   :widths: 22 20 48
   :header-rows: 1

   * - Equilibrium gas
     - Use
     - Why
   * - Pure CO2
     - ``CO2_Brine_Mixture``
     - Spycher-Pruess is more accurate than the S&W flash for pure CO2: 3.7% MARE against Yan (2011) where S&W scores 8.3%
   * - Any mixture, or any single gas other than pure CO2
     - ``SoreideWhitson``
     - The only route that handles mixtures, and the only one covering H2S, N2 and H2
   * - Methane only, oilfield units, following the textbook
     - ``brine_props``
     - Duan-Mao/McCain Ch. 4 solubility. Use ``SoreideWhitson`` instead if you want the S&W flash or any other gas
   * - Gas-free brine viscosity on its own
     - ``brine_viscosity``
     - Multi-salt, and the baseline the other three multiply

``brine.recommended_method()`` returns the same recommendation programmatically,
along with the solubility model it implies and the reason, so a script can record
which model produced a number:

.. code-block:: python

    >>> from pyrestoolbox import brine
    >>> brine.recommended_method(y_CO2=1.0)[0]
    'CO2_Brine_Mixture'
    >>> brine.recommended_method(y_CO2=0.5)[0]
    'SoreideWhitson'

What each returns
----------------------

**brine_props** — Methane-saturated brine using IAPWS-IF97 freshwater density with Spivey salt correction per McCain Petroleum Reservoir Fluid Properties pg 160. Includes effect of user specified salt concentration and degree of methane saturation.
Returns tuple of (Bw (rb/stb), Density (sg), viscosity (cP), Compressibility (1/psi), Rw GOR (scf/stb)).
**Changed in 3.7.4:** the returned viscosity now includes the dissolved-methane
correction. Before 3.7.4 it returned the gas-free viscosity even at
``ch4_sat=1.0``, understating the measured effect (Ostermann 1985: +3 to 6%) by
its full size. Expect the reported viscosity to rise by up to a few percent at
high methane saturation; density, Bw, Cw and Rsw are unchanged.

**CO2_Brine_Mixture** — CO2-saturated brine via Spycher-Pruess mutual solubility model. Returns a class object with calculated CO2 saturated brine property attributes. Retained deliberately rather than folded into ``SoreideWhitson``, because it is the more accurate route for pure CO2.

**SoreideWhitson** — Multicomponent gas-saturated brine via the Soreide-Whitson VLE model, using by default the refreshed BIP relationships of `Burgoyne & Nielsen (2026) <https://doi.org/10.1016/j.fluid.2026.114824>`_. Supports mixtures of C1, C2, C3, nC4, CO2, H2S, N2 and H2 in fresh or saline water.

Unit System Support
----------------------

All three brine models default to oilfield units (psia, degF, 1/psi, scf/stb) and accept ``metric=True`` to switch to Eclipse METRIC (barsa, degC, 1/bar, sm3/sm3):

- ``brine_props``: ``metric=False`` (default).
- ``CO2_Brine_Mixture``: ``metric=False`` (default).
- ``SoreideWhitson``: ``metric=False`` (default).

All three models compute freshwater density with IAPWS-IF97 Region 1, which is
valid up to the Region 1 limit of 100 MPa (14,503 psia). Supplying a higher
pressure raises ``ValueError``.

.. note::

    All "standard" volumes (Bw, Rs) use oilfield standard conditions (60 deg F, 14.696 psia) regardless of unit system.

.. list-table:: Method employed for different calculations (CO2_Brine_Mixture)
   :widths: 30 40
   :header-rows: 1

   * - Property
     - Calculation method
   * - Mutual Solubility between CO2 and Brine
     - Spycher & Pruess (2010), modified SRK Cubic EOS method
   * - Pure Water Density
     - IAPWS-IF97 Region 1 (international reference standard)
   * - Brine Salinity Correction
     - Spivey et al. (modified), per "Petroleum Reservoir Fluid Property Correlations", (McCain, Spivey & Lenn: Chapter 4)
   * - CO2 Corrected Brine Density
     - Molar volume of dissolved CO2 estimated with Garcia (2001) equation, used with xCO2 calculated from Spycher & Pruess, and CO2-free brine density to calculate insitu density
   * - Gas-free Brine viscosity
     - IAPWS-2008 water (Huber et al. 2009) x ion-additive Jones-Dole salt ratio (Appelo & Parkhurst, as implemented in PHREEQC) x Kestin (1981) measured pressure factor. **Changed in 3.7.3**; Mao-Duan (2009) remains available as ``route='mao_duan'``. See the ``brine_viscosity`` section below.
   * - CO2 Corrected Brine Viscosity
     - Calabrese et al. (2019) Eq. 25, a temperature-dependent increment fitted to 415 measured brine points. **Changed in 3.7.2**, superseding Islam-Carlson (2012), which carries no temperature dependence and overstates the effect by 3.5x at 366 K and 9.4x at 423 K.
     

pyrestoolbox.brine.brine_props
======================

.. code-block:: python

    brine_props(p, degf, wt=0, ch4_sat=0, metric=False) -> tuple

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float
     - Pressure (psia, or barsa if metric=True)
   * - degf
     - float
     - Temperature (deg F, or deg C if metric=True)
   * - wt
     - float
     - Salt weight% in the brine (0 - 100)
   * - ch4_sat
     - float
     - Degree of methane saturation (0 - 1). 0 = No Methane, 1 = 100% Methane saturated
   * - metric
     - bool
     - If True, treats input pressure & temperature as metric, otherwise treats as Field units. Default False

.. list-table:: Returns (tuple)
   :widths: 10 15 40
   :header-rows: 1

   * - Index
     - Type
     - Description
   * - [0]
     - float
     - Bw — water formation volume factor (rb/stb)
   * - [1]
     - float
     - Density (specific gravity relative to water)
   * - [2]
     - float
     - Viscosity (cP)
   * - [3]
     - list
     - Compressibility [undersaturated, saturated] (1/psi, or 1/barsa if metric=True). The undersaturated value (Cw[0]) is the isothermal compressibility of the brine at constant dissolved gas content. The saturated value (Cw[1]) is a pseudo-compressibility representing the average compressibility of the brine and differentially evolved gas system, accounting for both liquid compression and gas exsolution volume changes
   * - [4]
     - float
     - Rsw — solution gas-water ratio (scf/stb, or sm3/sm3 if metric=True)

Examples:

.. code-block:: python

    >>> from pyrestoolbox import brine
    >>> bw, lsg, visw, cw, rsw = brine.brine_props(p=160, degf=135, wt=1.5, ch4_sat=1.0)
    >>> print('Bw:', bw)
    >>> print('SGw:', lsg)
    >>> print('Visw:', visw)
    >>> print('Cw_usat:', cw[0])
    >>> print('Cw_sat:', cw[1])
    >>> print('Rsw:', rsw)
    Bw: 1.0152005799432318
    SGw: 0.9950108380379709
    Visw: 0.5047732519299403
    Cw_usat: 2.9696277255527504e-06
    Cw_sat: 0.0001539877228225709
    Rsw: 1.2567682353688225

.. note::

    When ``metric=True``, Cw values are returned in 1/barsa (instead of 1/psi) and Rsw in sm3/sm3 (instead of scf/stb).
    Bw (rb/stb), density (SG), and viscosity (cP) are unchanged by the unit system.

pyrestoolbox.brine.CO2_Brine_Mixture
======================

.. code-block:: python

    CO2_Brine_Mixture(pres, temp, ppm = 0, metric = False, cw_sat = False) -> class

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - pres
     - float
     - Pressure (Bar or psia)
   * - temp
     - float
     - Temperature (deg C or deg F)
   * - ppm
     - float
     - Parts per million (Wt) NaCl equivalent in brine (1% NaCl equivalent by Wt = 10,000 ppm)
   * - metric
     - Boolean
     - If True, treats input pressure & temperature as metric, otherwise treats as Field units.
   * - cw_sat
     - Boolean
     - If True, will also calculate saturated brine compressibility (doubles calculation time). Default is False.
     
     
.. list-table:: Returns (CO2_Brine_Mixture)
   :widths: 10 15 40
   :header-rows: 1

   * - Attribute
     - Type
     - Description
   * - .x
     - np.ndarray
     - CO2 and H2O in aqueous phase [xCO2, xH2O]
   * - .y
     - np.ndarray
     - CO2 and H2O in vapor phase [yCO2, yH2O]
   * - .xSalt
     - float
     - Mole fraction of NaCl in brine
   * - .rhoGas
     - float
     - CO2 rich gas density (gm/cm3)
   * - .bDen
     - list
     - Brine density [CO2 Saturated, Pure Brine, Freshwater] (gm/cm3)
   * - .bVis
     - list
     - Brine viscosity [CO2 Saturated, Pure Brine, Freshwater] (cP)
   * - .bVisblty
     - float
     - CO2 Saturated brine viscosibility (1/Bar or 1/psi)
   * - .bw
     - list
     - Brine formation volume factor [CO2 Saturated, Pure Brine, Freshwater] (rm3/sm3 or rb/stb)
   * - .Rs
     - float
     - CO2 Saturated Brine solution gas ratio, relative to standard conditions (sm3/sm3 or scf/stb)
   * - .Cf_usat
     - float
     - Brine undersaturated compressibility (1/Bar or 1/psi)
   * - .Cf_sat
     - float
     - Brine saturated pseudo-compressibility (1/Bar or 1/psi). Represents the average compressibility of the brine and differentially evolved gas system, accounting for both liquid compression and gas exsolution. Requires cw_sat input to be set True to calculate

                
Examples:

Usage example for 5000 psia x 275 deg F and 3% NaCl brine:

.. code-block:: python

    >>> from pyrestoolbox import brine
    >>> mix = brine.CO2_Brine_Mixture(pres = 5000, temp = 275, ppm = 30000, metric = False)
    >>> mix.bw  # Returns [CO2 Saturated, Pure Brine, Freshwater]
    [1.1091672843736888, 1.054302417027164, 1.0542033928155845]
    >>> mix.x  # Returns molar fractions in aqueous phase [xCO2, xH2O]
    array([0.02431225, 0.95743209])
    
Usage example for 175 Bara x 85 degC and 0% NaCl brine:

.. code-block:: python

    >>> mix = brine.CO2_Brine_Mixture(pres = 175, temp = 85, metric = True)
    >>> mix.Rs  # Returns sm3 dissolved CO2 / sm3 Brine
    24.743651168969475

pyrestoolbox.brine.make_pvtw_table
======================

.. note::

    In v3.0, the primary entry point for this function is ``simtools.make_pvtw_table()``.
    The ``brine.make_pvtw_table()`` function remains as a backward-compatible wrapper.

.. code-block:: python

    make_pvtw_table(pi, degf, wt=0, ch4_sat=0, pmin=500, pmax=10000, nrows=20, export=False, metric=False) -> dict

Generates a PVTW (water PVT) table over a pressure range using brine_props (IAPWS-IF97 freshwater with Spivey salt correction).
Optionally exports ECLIPSE PVTW keyword file and Excel spreadsheet.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - pi
     - float
     - Initial (reference) pressure (psia, or barsa if metric=True)
   * - degf
     - float
     - Temperature (deg F, or deg C if metric=True)
   * - wt
     - float
     - Salt weight% in the brine (0 - 100). Default 0
   * - ch4_sat
     - float
     - Degree of methane saturation (0 - 1). Default 0
   * - pmin
     - float
     - Minimum table pressure (psia, or barsa if metric=True). Default 500
   * - pmax
     - float
     - Maximum table pressure (psia, or barsa if metric=True). Default 10000
   * - nrows
     - int
     - Number of table rows. Default 20
   * - export
     - bool
     - If True, writes PVTW.INC and pvtw_table.xlsx. Default False
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units (barsa, deg C, sm3/sm3, 1/bar). Default False

.. list-table:: Returns (dict)
   :widths: 10 15 40
   :header-rows: 1

   * - Key
     - Type
     - Description
   * - table
     - DataFrame
     - Pressure, Bw, Density, Viscosity, Cw_usat, Cw_sat, Rsw
   * - pref
     - float
     - Reference pressure (psia, or barsa if metric=True)
   * - bw_ref
     - float
     - Bw at reference pressure (rb/stb)
   * - cw_ref
     - list
     - Compressibility [undersaturated, saturated] at reference pressure (1/psi, or 1/bar if metric=True)
   * - visw_ref
     - float
     - Viscosity at reference pressure (cP)
   * - rsw_ref
     - float
     - Rsw at reference pressure (scf/stb, or sm3/sm3 if metric=True)
   * - den_ref
     - float
     - Density (sg) at reference pressure

Examples:

.. code-block:: python

    >>> from pyrestoolbox import brine
    >>> result = brine.make_pvtw_table(pi=3000, degf=200, wt=0, ch4_sat=0)
    >>> result['bw_ref']
    1.027589195773527
    >>> result['cw_ref']
    [3.0887176266534516e-06, 3.0887176266534516e-06]
    >>> result['visw_ref']
    0.308131761431705

pyrestoolbox.brine.brine_viscosity
======================

::

    brine_viscosity(T_K, P_MPa, composition=None, salts=None, m=None, S=None, route=None, pressure_term=True) -> float

Gas-free brine viscosity in mPa s (equivalently cP). New in 3.7.3, and now the
baseline used internally by ``brine_props`` and ``SoreideWhitson``.

::

    mu = mu_water(T, P) x salt_ratio(T, ions) x pressure_factor(T, P, I)

**MULTI-SALT.** The salt ratio is additive over ions rather than fitted to one
salt, so any combination of the parameterised ions can be supplied. Exactly one
salinity input is required, and **a salinity given with no species named means
NaCl**, matching the density side:

.. code-block:: python

    from pyrestoolbox import brine

    brine.brine_viscosity(373.15, 30.0, m=3.0)                      # 3 molal NaCl
    brine.brine_viscosity(373.15, 30.0, S=0.1492)                    # weight fraction NaCl
    brine.brine_viscosity(373.15, 30.0, salts={'NaCl': 2.0, 'CaCl2': 1.0})
    brine.brine_viscosity(373.15, 30.0, composition={'Na+': 3.0, 'Cl-': 3.0})

    brine.salt_ratio(373.15, 30.0, salts={'KCl': 4.0})               # the ratio alone

``salts=`` accepts NaCl, KCl, CaCl2, MgCl2, CaBr2, Na2SO4, SrCl2 and BaCl2.
``composition=`` accepts the ions directly, which is the more general form.
An ion with no viscosity parameterisation raises rather than contributing zero.

.. list-table:: Inputs
   :widths: 12 12 46
   :header-rows: 1

   * - Input
     - Units
     - Description
   * - T_K
     - K
     - Temperature
   * - P_MPa
     - MPa
     - Pressure
   * - m
     - mol/kg
     - NaCl molality. One of m, S, salts or composition must be supplied
   * - S
     - fraction
     - NaCl weight fraction
   * - salts
     - dict
     - {salt name: molality}, e.g. {'NaCl': 2.0, 'CaCl2': 1.0}
   * - composition
     - dict
     - {ion: molality}, e.g. {'Na+': 3.0, 'Cl-': 3.0}
   * - route
     - str
     - 'jones_dole' (default, multi-salt) or 'mao_duan' (NaCl only, pre-3.7.3)
   * - pressure_term
     - bool
     - Apply Kestin's measured pressure factor. Default True

.. list-table:: Accuracy against measured viscosity (salt ratio only)
   :widths: 30 15 15
   :header-rows: 1

   * - Data
     - Mean
     - Max
   * - NaCl, Kestin (1981), 298-423 K, 0.1-35 MPa, 0.5-5 mol/kg
     - 0.300%
     - 0.78%
   * - KCl, Kestin (1981) companion tables, same box
     - 0.764%
     - 4.62%

**BEHAVIOUR CHANGE IN 3.7.3.** Existing NaCl viscosities move by 0.34% on
average and 1.36% at most, and on freshwater by +0.069%, which is the water leg
alone. This is a change of route rather than a correction of an error: on NaCl
the old and new salt ratios are statistically indistinguishable (0.443% against
0.415% mean against Kestin, both inside his stated +/-0.5%). What the new
baseline buys is the measured pressure dependence and salts other than sodium
chloride. Pass ``route='mao_duan'`` to reproduce pre-3.7.3 numbers exactly; that
route is NaCl only and raises on any other salt rather than silently treating it
as NaCl, which would cost 15.9% on average and 51.2% at worst on KCl.

**Three limits ship with this baseline and should be quoted with any result.**

1. **Ion mixing is not validated against measurement.** Every brine scored above
   is a single salt, because no measured viscosity of a mixed brine at known
   composition, temperature and pressure was found. The multi-salt capability
   reproduces the reference implementation to 0.0015% over 145 cases, which
   tests this implementation rather than the model. Describe mixed-brine results
   in those terms.
2. **Above 35 MPa the pressure factor is held, not extrapolated**, so the
   baseline is pressure-blind from 35 to 100 MPa. Kestin's own factor at 35 MPa
   is the largest pressure correction applied anywhere in the chain.
3. **The K+ worst case is 4.6%**, at 298 to 323 K and 4 to 5 mol/kg. Na+ is
   worst at 0.78% over the same box, so the ion-additive form is better on
   average and far better away from NaCl, not uniformly better.

No PHREEQC installation is required: the per-ion parameters are embedded and are
machine-verified against the source database.

pyrestoolbox.brine.SoreideWhitson
======================

.. code-block:: python

    SoreideWhitson(pres, temp, ppm=0, y_CO2=0, y_H2S=0, y_N2=0, y_H2=0, sg=0.65, metric=False, cw_sat=False, framework='proposed', salinity_method='gamma_phi', vphi_route='auto', *, p=None, degf=None, wt=None) -> class

Soreide-Whitson VLE model for multicomponent gas solubility in water/brine, using by
default the refreshed BIP relationships of `Burgoyne & Nielsen (2026) <https://doi.org/10.1016/j.fluid.2026.114824>`_,
with mass-balance density corrections and calibrated viscosity corrections. Supports gas mixtures containing any combination
of C1, C2, C3, nC4, CO2, H2S, N2 and H2 in fresh or saline water.

The hydrocarbon portion of the gas (1 - y_CO2 - y_H2S - y_N2 - y_H2) is automatically split among C1-C4
based on the gas specific gravity using constrained exponential decay to match the implied HC molecular weight.

.. list-table:: Method employed for different calculations (SoreideWhitson)
   :widths: 30 40
   :header-rows: 1

   * - Property
     - Calculation method
   * - Gas-Brine Equilibrium
     - Soreide-Whitson Peng-Robinson EOS VLE flash, refreshed BIPs of Burgoyne & Nielsen (2026) by default
   * - Pure Water Density
     - IAPWS-IF97 Region 1 (international reference standard)
   * - Brine Salinity Correction
     - Spivey et al. (modified), per "Petroleum Reservoir Fluid Property Correlations", (McCain, Spivey & Lenn: Chapter 4)
   * - Gas-Corrected Brine Density
     - Mass and volume balance (an identity, not a fitted mixing rule) with apparent molar volumes from the Soreide-Whitson modified PR route plus one volume shift per gas. **Changed in 3.7.2**; the Plyasunov (2019-2021) correlation is retained as ``vphi_route='plyasunov'``. A literature-anchored salinity shift is applied to V_phi from 3.7.3 (relative form from 3.7.4), disabled with ``salt_effect=False``.
   * - Gas-free Brine Viscosity
     - IAPWS-2008 water x ion-additive Jones-Dole salt ratio x Kestin measured pressure factor. **Changed in 3.7.3**; Mao-Duan (2009) remains available. See the ``brine_viscosity`` section below.
   * - Gas-Corrected Brine Viscosity
     - Per-gas multiplicative corrections: Calabrese et al. (2019) Eq. 25 for CO2 (**changed in 3.7.2**, superseding Islam-Carlson), an Arrhenius-times-Langmuir form refitted to all 23 Ostermann (SPE-14211, 1985) measurements for CH4, and Murphy & Gaines (1974) for H2S. C2H6 and N2 are measured nulls; H2, C3H8 and nC4H10 carry no correction because no data exists.

**Apparent molar volume route** (``vphi_route``, new in 3.7.2). ``'auto'`` is the
default and uses the Soreide-Whitson modified PR route with one dimensionless
volume shift per gas, falling back to the Plyasunov correlation for C3H8 and
nC4H10 and outside the PR route's validity box. ``'pr'`` and ``'plyasunov'``
force one route; ``'plyasunov'`` reproduces pre-3.7.2 results exactly. The PR
route needs one fitted number per gas against 35 coefficients per gas for the
correlation, matches or beats it on five of six gases against the calibration
densimetry, and reproduces an H2S temperature trend the correlation misses.
The delivered shifts are CH4 -0.109632, CO2 -0.037913, H2S -0.078975,
N2 -0.155288, H2 -0.177625, C2H6 -0.073142, C3H8 -0.112963 and
nC4H10 +0.110924.

Three ceilings should not be conflated: the arithmetic stops at 623.15 K, the
volume shift stops being fitted at 473.15 K, and **accuracy is claimed only to
about 450 K (300 to 350 degF)**, which is the only one of the three that is a
statement about the answer. N2, H2, C2H6, C3H8 and nC4H10 have no
temperature-resolved calibration data, so for those the temperature behaviour is
the equation of state unaided.

**Salinity shift on V_phi** (new in 3.7.3; relative form from 3.7.4). A
gas-generic dimensionless fraction is applied,
``V_eff = V_phi(T,P) * (1 + g(m))`` with
``g(m) = -1.7009 m/(1 + 0.090684 m)`` percent, giving -1.56% at 1 mol/kg and
-5.85% at 5 mol/kg (for CO2, -0.5 to -2.7 cm3/mol over the working range). It
is fixed entirely from Tiepel and Gubbins (1972) KCl dilatometry with no
parameter fitted to any brine-density data; the magnitude is known to about a
factor of two. Freshwater results are unchanged exactly; gas-saturated brine
densities move by up to a few hundredths of a percent. Pass
``salt_effect=False`` to ``brine.vphi_route.V_phi`` to recover freshwater
V_phi.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - pres
     - float
     - Pressure (Bar or psia)
   * - temp
     - float
     - Temperature (deg C or deg F)
   * - ppm
     - float
     - Parts per million (Wt) NaCl equivalent in brine. Default 0
   * - y_CO2
     - float
     - Mole fraction of CO2 in dry gas. Default 0
   * - y_H2S
     - float
     - Mole fraction of H2S in dry gas. Default 0
   * - y_N2
     - float
     - Mole fraction of N2 in dry gas. Default 0
   * - y_H2
     - float
     - Mole fraction of H2 in dry gas. Default 0
   * - sg
     - float
     - Gas specific gravity, used to estimate HC split among C1-C4. Default 0.65
   * - metric
     - Boolean
     - If True, treats input pressure & temperature as metric, otherwise treats as Field units. Default False
   * - cw_sat
     - Boolean
     - If True, will also calculate saturated brine compressibility. Default False
   * - framework
     - str
     - VLE framework used by the S&W engine. ``'proposed'`` (default, the `Burgoyne & Nielsen (2026) <https://doi.org/10.1016/j.fluid.2026.114824>`_ refreshed BIPs), ``'sw_original'`` (original 1992 published coefficients), or ``'dropin'`` (PR-EOS fitted with brine-aware water alpha). Affects kij and ks correlations.
   * - salinity_method
     - str
     - How salinity enters the flash. ``'gamma_phi'`` (default, Sechenov salting-out via activity coefficient), ``'embedded'`` (salinity inside kij — only compatible with ``'dropin'``/``'sw_original'``), ``'explicit'`` (brine treated as a component in the flash), ``'sechenov'`` and ``'auto'`` (both accepted aliases that normalise to ``gamma_phi``). ``framework='proposed'`` + ``salinity_method='embedded'`` emits a warning and falls back to ``gamma_phi``. Note that the Rust-accelerated flash runs only for ``framework='proposed'`` with ``salinity_method='gamma_phi'``; other combinations use the pure-Python flash.

.. list-table:: Returns (SoreideWhitson)
   :widths: 10 15 40
   :header-rows: 1

   * - Attribute
     - Type
     - Description
   * - .x
     - dict
     - Dissolved gas mole fractions in aqueous phase, e.g. {'CO2': 0.024, 'CH4': 0.0015}
   * - .x_total
     - float
     - Total dissolved gas mole fraction (sum of all x_i)
   * - .y
     - dict
     - Gas phase compositions (dry basis, normalized)
   * - .y_H2O
     - float
     - Water mole fraction in gas phase
   * - .water_content
     - dict
     - Water content of gas: {'y_H2O': ..., 'stb_mmscf': ..., 'lb_mmscf': ...}
   * - .bDen
     - list
     - Brine density [Gas Saturated, Gas-Free Brine, Freshwater] (gm/cm3)
   * - .bVis
     - list
     - Brine viscosity [Gas Saturated, Gas-Free Brine, Freshwater] (cP)
   * - .bVisblty
     - float
     - Gas-saturated brine viscosibility (1/Bar or 1/psi)
   * - .bw
     - list
     - Brine formation volume factor [Gas Saturated, Gas-Free Brine, Freshwater] (rm3/sm3 or rb/stb)
   * - .Rs
     - dict
     - Per-gas solution ratios, e.g. {'CO2': 15.2, 'CH4': 2.1} (sm3/sm3 or scf/stb)
   * - .Rs_total
     - float
     - Total solution gas ratio (sum of all per-gas Rs) (sm3/sm3 or scf/stb)
   * - .Cf_usat
     - float
     - Brine undersaturated compressibility (1/Bar or 1/psi)
   * - .Cf_sat
     - float
     - Brine saturated pseudo-compressibility (1/Bar or 1/psi). Represents the average compressibility of the brine and differentially evolved gas system, accounting for both liquid compression and gas exsolution. Requires cw_sat input to be set True to calculate
   * - .gas_comp
     - dict
     - Normalized gas composition used (including estimated HC split from SG)
   * - .MwBrine
     - float
     - Molecular weight of gas-free brine (g/mol)

Examples:

Pure CO2 case at 5000 psia x 275 deg F and 3% NaCl brine:

.. code-block:: python

    >>> from pyrestoolbox import brine
    >>> mix = brine.SoreideWhitson(pres=5000, temp=275, ppm=30000, y_CO2=1.0, metric=False)
    >>> mix.bDen  # Returns [Gas Saturated, Gas-Free Brine, Freshwater]
    [0.973398983496999, 0.968164592979362, 0.9476497407774847]
    >>> mix.Rs  # Returns per-gas Rs dict (scf/stb)
    {'CO2': 139.5790835259016}
    >>> mix.bw  # Returns [Gas Saturated, Gas-Free, Freshwater]
    [1.096613263262556, 1.0543023174291248, 1.0542033190822462]

Pure CH4 case (SG=0.554) at 5000 psia x 275 deg F and 3% NaCl brine:

.. code-block:: python

    >>> mix = brine.SoreideWhitson(pres=5000, temp=275, ppm=30000, y_CO2=0, sg=0.554, metric=False)
    >>> mix.Rs
    {'CH4': 21.21234560600256}
    >>> mix.bDen
    [0.9641838754552089, 0.968164592979362, 0.9476497407774847]

Mixed gas (10% CO2, 5% H2S, SG=0.7) at 200 Bar x 80 degC and 10,000 ppm NaCl:

.. code-block:: python

    >>> mix = brine.SoreideWhitson(pres=200, temp=80, ppm=10000, y_CO2=0.1, y_H2S=0.05, sg=0.7, metric=True)
    >>> mix.gas_comp  # Estimated gas composition including HC split
    {'CO2': 0.1, 'H2S': 0.05, 'CH4': 0.8133, 'C2H6': 0.0351, 'C3H8': 0.0015, 'nC4H10': 0.0001}
    >>> mix.Rs_total  # Total dissolved gas (sm3/sm3)
    8.5106741070893
    >>> mix.bDen
    [0.9855934589486185, 0.9871360082710434, 0.9804911502375318]

Pure CO2 fresh water at 175 Bar x 85 degC with saturated compressibility:

.. code-block:: python

    >>> mix = brine.SoreideWhitson(pres=175, temp=85, ppm=0, y_CO2=1.0, metric=True, cw_sat=True)
    >>> mix.Rs_total  # sm3 dissolved CO2 / sm3 Brine
    24.188037633302223
    >>> mix.Cf_sat
    0.00016012590421810821
    >>> mix.water_content
    {'y_H2O': 0.013982317814779299, 'stb_mmscf': 1.923030543083137, 'lb_mmscf': 673.2379475216799}

References:

- Soreide, I. and Whitson, C.H., "Peng-Robinson Predictions for Hydrocarbons, CO2, N2, and H2S with Pure Water and NaCl Brine", Fluid Phase Equilibria, 77, 217-240, 1992.
- Calabrese, C., McBride-Wright, M., Maitland, G.C. and Trusler, J.P.M. (2019), J. Chem. Eng. Data, 64, 3831-3847. (CO2 viscosity increment, and the CO2-brine density validation set)
- Ostermann, R.D. et al. (1985), "The Effect of Dissolved Gas on Reservoir Brine Viscosity", SPE 14211.
- Murphy, W.R. and Gaines, T.M. (1974), J. Chem. Eng. Data, 19(4), 359-362. (H2S)
- Huber, M.L. et al. (2009), "New International Formulation for the Viscosity of H2O", J. Phys. Chem. Ref. Data, 38(2), 101-125.
- Kestin, J., Khalifa, H.E. and Correia, R.J. (1981), J. Phys. Chem. Ref. Data, 10, 71. (measured NaCl viscosity, and the pressure factor)
- Appelo, C.A.J., Parkhurst, D.L. and Post, V.E.A. (2014), Geochim. Cosmochim. Acta, 125, 49-67. (ion-additive salt ratio)
- Hnedkovsky, L., Wood, R.H. and Majer, V. (1996), J. Chem. Thermodynamics, 28, 125-142. (the densimetry the volume shifts are fitted to)
- Plyasunov, A.V., "Values of the Apparent Molar Volumes...," Fluid Phase Equilibria, Parts I-IV, 2019-2021. (the fallback V_phi route)
- Garcia, J.E., "Density of Aqueous Solutions of CO2", LBNL Report 49023, 2001. (the CO2 V_phi cubic used by CO2_Brine_Mixture)
- Islam, A.W. and Carlson, E.S. (2012), "Viscosity Models and Effects of Dissolved CO2", Energy & Fuels, 26(8), 5330-5336. (superseded by Calabrese in 3.7.2; listed for provenance)

Density and Viscosity Calculation Details
----------------------

**Gas-corrected brine density** is a mass and volume balance. It is an identity rather than a
fitted mixing rule: mass is additive, and the apparent molar volume is *defined* so that the
solution volume equals the gas-free solvent volume plus the dissolved moles times V_phi.
Nothing in it is fitted or fittable, and it is gas-generic by construction.

For each dissolved species the apparent molar volume at infinite dilution V_phi(T,P) comes from
the Soreide-Whitson modified PR equation of state, using the exact relation
V_bar_2 = -(dP/dn2)/(dP/dV) at fixed temperature evaluated on the water-rich liquid root, plus
one dimensionless volume shift per gas. **This became the default in 3.7.2**; the Plyasunov
(2019-2021) A12-infinity model is retained as the fallback for C3H8 and nC4H10 and outside the
PR route's validity box, and is selectable with ``vphi_route='plyasunov'``. A literature-anchored
salinity shift is applied to V_phi from 3.7.3, as a relative fraction from 3.7.4.

For mixed dissolved gases, mole-fraction-weighted effective V_phi and MW are used. The weighting
is an exact algebraic identity; the approximation is the physics it rests on, namely that each
gas contributes its infinite-dilution volume with no solute-solute cross term, which is exact
only as the dissolved amount tends to zero:

.. code-block:: text

    V_phi_eff = sum(yi * V_phi_i) / sum(yi)
    MW_eff    = sum(yi * MW_i) / sum(yi)

    rho = (1 + x_gas * MW_eff / (Mw * x_w)) / (x_gas * V_phi_eff / (Mw * x_w) + 1/rho_brine)

**SOLVENT BASIS, which is easy to get wrong.** V_phi is defined against the whole gas-free
brine, so ``x_gas * V_phi`` adds to the volume of water *and salt together*. Where the dissolved
mole fraction is returned on a salt-free basis, as VLE models normally return it, the solvent
term ``Mw * x_w`` is the mass of water only and must be grossed up to brine mass by dividing by
(1 - S), with S the salt weight fraction. Pairing the water molecular weight with a brine
density instead omits the salt from both mass and volume and inflates the whole dissolved-gas
density effect by a factor approaching 1/(1 - S): 5.6% at 1 mol/kg NaCl, 11% at 2 and 28% at 5.
The two forms are identical in freshwater. This was corrected in 3.7.2.

**Which gases densify.** The sign is that of ``MW_gas - rho_brine * V_phi``, *not* of
MW/V_phi against 1. Of the gases in scope only CO2 is positive: at 298.15 K and 30 MPa it
carries 44.0 g/mol against a displaced 33.3 cm3/mol, a margin of +10.4 g/mol. **H2S lightens
brine**, by -0.9 g/mol at 298.15 K widening to -7.3 g/mol by 453 K, and above about 327 K it
lightens brine more than any other gas in the set. H2S is the near-miss case and at ambient
temperature its sign should be read as probable rather than established. CH4, N2, H2, C2H6,
C3H8 and nC4H10 lighten brine throughout.

**Gas-corrected brine viscosity** applies multiplicative per-gas corrections to the gas-free
baseline described in the ``brine_viscosity`` section above. Every correction is calibrated
directly to measured viscosity; none is scaled from the density change, which experiment
refutes, getting the sign wrong for three of five gases.

.. list-table:: Viscosity correction by gas
   :widths: 12 46 26
   :header-rows: 1

   * - Gas
     - Correction
     - Reference
   * - CO2
     - ln(mu/mu_brine) = e1 exp(-e2 (T/T0 - 1)) x, with e1 = 65.560, e2 = 2.468, T0 = 142 K
     - Calabrese et al. (2019) Eq. 25, 415 measured brine points. **Changed in 3.7.2**, superseding Islam-Carlson
   * - CH4
     - mu x (1 + A exp(B/T) x/(K + x)), with A = 1.71739196e-3, B = 1239.77535 K, K = 1.52860547e-3
     - Refitted here to all 23 Ostermann (SPE-14211, 1985) measurements rather than to his three summary values
   * - H2S
     - mu x (1 + 1.70 x)
     - Murphy-Gaines (1974), on their own solubility source. **About 27% uncertain**; five points, all below 309 K
   * - C2H6, N2
     - No correction
     - Measured nulls. N2's null bounds the effect below about 2.5% rather than establishing zero
   * - H2, C3H8, nC4H10
     - No correction
     - **No measurement exists.** The zero is a conservative default, not evidence; treat H2 as 0 to several percent

For gas mixtures the corrections are applied multiplicatively across all dissolved species.
That is a design choice rather than a validated one: no mixed-dissolved-gas viscosity data
exist, and the assumption plausibly overstates the combined effect at high total dissolved gas.

