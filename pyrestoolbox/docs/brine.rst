===================================
Brine PVT with differing degrees of methane, CO2, or multicomponent gas saturation
===================================

Three brine property models are available:

**brine_props** — Methane-saturated brine using IAPWS-IF97 freshwater density with Spivey salt correction per McCain Petroleum Reservoir Fluid Properties pg 160. Includes effect of user specified salt concentration and degree of methane saturation.
Returns tuple of (Bw (rb/stb), Density (sg), viscosity (cP), Compressibility (1/psi), Rw GOR (scf/stb))

**CO2_Brine_Mixture** — CO2-saturated brine via Spycher-Pruess mutual solubility model. Returns a class object with calculated CO2 saturated brine property attributes.

**SoreideWhitson** — Multicomponent gas-saturated brine via Soreide-Whitson (1992) VLE model. Supports mixtures of C1, C2, C3, nC4, CO2, H2S, N2 and H2 in fresh or saline water.

Unit System Support
----------------------

All three brine models support both oilfield and metric unit systems via a ``metric`` parameter:

- ``brine_props``: ``metric=False`` (default). When ``metric=True``, pressure is in barsa, temperature is in degC, compressibility is in 1/barsa, and Rsw is in sm3/sm3.
- ``CO2_Brine_Mixture``: ``metric=True`` (default). When ``metric=False``, pressure is in psia and temperature is in degF.
- ``SoreideWhitson``: ``metric=True`` (default). When ``metric=False``, pressure is in psia and temperature is in degF.

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
   * - Pure Brine viscosity
     - Mao-Duan (2009).
   * - CO2 Corrected Brine Viscosity
     - "Viscosity Models and Effects of Dissolved CO2", Islam-Carlson (2012) to adjust the pure brine viscosity for xCO2 calculated from Spycher & Pruess.
     

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

Examples:

.. code-block:: python

    >>> from pyrestoolbox import brine
    >>> bw, lsg, visw, cw, rsw = brine.brine_props(p=160, degf=135, wt=1.5, ch4_sat=1.0)
    >>> print('Bw:', bw)
    >>> print('SGw:', lsg)
    >>> print('Visw:', visw)
    >>> print('Cw:', cw)
    >>> print('Rsw:', rsw)
    Bw: 1.0152007040056148
    SGw: 0.9950108179684295
    Visw: 0.4994004662758671
    Cw: 0.0001539690974662865
    Rsw: 1.2540982731813703

.. note::

    When ``metric=True``, Cw is returned in 1/barsa (instead of 1/psi) and Rsw in sm3/sm3 (instead of scf/stb).
    Bw (rb/stb), density (SG), and viscosity (cP) are unchanged by the unit system.

pyrestoolbox.brine.CO2_Brine_Mixture
======================

.. code-block:: python

    CO2_Brine_Mixture(pres, temp, ppm = 0, metric = True) -> class

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
     
     
.. list-table:: Results
   :widths: 10 15 40
   :header-rows: 1

   * - Class Attribute
     - Unit
     - Description
   * - .x
     - Mole fractions
     - CO2 and H2O in aqueous phase [xCO2, xH2O]
   * - .y
     - Mole fractions
     - CO2 and H2O in vapor phase [yCO2, yH2O]
   * - .xSalt
     - Mole Fraction
     - Mole fraction of NaCl in brine
   * - .rhoGas
     - (gm/cm3)
     - CO2 rich gas density.
   * - .bDen
     - (gm/cm3)
     - Brine density [CO2 Saturated, Pure Brine, Freshwater]
   * - .bVis
     - (cP)
     - Brine viscosity [CO2 Saturated, Pure Brine, Freshwater]
   * - .bVisblty
     - (1/Bar or 1/psi)
     - CO2 Saturated brine viscosibility.
   * - .bw
     - (rm3/smr or rb/stb)
     - Brine formation volume factor  [CO2 Saturated, Pure Brine, Freshwater]
   * - .Rs
     - (sm3/sm3 or scf/stb)
     - CO2 Saturated Brine solution gas ratio, relative to standard conditions
   * - .Cf_usat
     - (1/Bar or 1/psi)
     - Brine undersaturated compressibility 
   * - .Cf_ssat
     - (1/Bar or 1/psi)
     - Brine saturated compressibility. Requires cw_sat input to be set True to calculate

                
Examples:

Usage example for 5000 psia x 275 deg F and 3% NaCl brine:

.. code-block:: python

    >>> from pyrestoolbox import brine
    >>> mix = brine.CO2_Brine_Mixture(pres = 5000, temp = 275, ppm = 30000, metric = False)
    >>> mix.bw  # Returns [CO2 Saturated, Pure Brine, Freshwater]
    [1.108578337107381, 1.054302417027164, 1.0542033928155845]
    >>> mix.x  # Returns molar fractions in aqueous phase [xCO2, xH2O]
    array([0.02431225, 0.95743209])
    
Usage example for 175 Bara x 85 degC and 0% NaCl brine:

.. code-block:: python

    >>> mix = brine.CO2_Brine_Mixture(pres = 175, temp = 85)
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

.. list-table:: Return dict keys
   :widths: 10 15 40
   :header-rows: 1

   * - Key
     - Type
     - Description
   * - table
     - DataFrame
     - Pressure, Bw, Density, Viscosity, Cw, Rsw
   * - pref
     - float
     - Reference pressure (psia, or barsa if metric=True)
   * - bw_ref
     - float
     - Bw at reference pressure (rb/stb)
   * - cw_ref
     - float
     - Compressibility at reference pressure (1/psi, or 1/bar if metric=True)
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
    3.0887176266534516e-06
    >>> result['visw_ref']
    0.3083544960904146

pyrestoolbox.brine.SoreideWhitson
======================

.. code-block:: python

    SoreideWhitson(pres, temp, ppm=0, y_CO2=0, y_H2S=0, y_N2=0, y_H2=0, sg=0.65, metric=True, cw_sat=False) -> class

Soreide-Whitson (1992) VLE model for multicomponent gas solubility in water/brine, with Garcia/Plyasunov
density corrections and calibrated viscosity corrections. Supports gas mixtures containing any combination
of C1, C2, C3, nC4, CO2, H2S, N2 and H2 in fresh or saline water.

The hydrocarbon portion of the gas (1 - y_CO2 - y_H2S - y_N2 - y_H2) is automatically split among C1-C4
based on the gas specific gravity using constrained exponential decay to match the implied HC molecular weight.

.. list-table:: Method employed for different calculations (SoreideWhitson)
   :widths: 30 40
   :header-rows: 1

   * - Property
     - Calculation method
   * - Gas-Brine Equilibrium
     - Soreide-Whitson (1992), Peng-Robinson EOS VLE flash
   * - Pure Water Density
     - IAPWS-IF97 Region 1 (international reference standard)
   * - Brine Salinity Correction
     - Spivey et al. (modified), per "Petroleum Reservoir Fluid Property Correlations", (McCain, Spivey & Lenn: Chapter 4)
   * - Gas-Corrected Brine Density
     - Garcia (2001) Eq. 18 mixing rule with Plyasunov (2019-2021) apparent molar volumes for each dissolved gas
   * - Pure Brine Viscosity
     - Mao-Duan (2009)
   * - Gas-Corrected Brine Viscosity
     - Per-gas multiplicative corrections: Islam-Carlson (2012) for CO2, Murphy-Gaines (1974) for H2S, Ostermann (SPE-14211, 1985) for CH4

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
     - If True, treats input pressure & temperature as metric, otherwise treats as Field units. Default True
   * - cw_sat
     - Boolean
     - If True, will also calculate saturated brine compressibility. Default False

.. list-table:: Results
   :widths: 10 15 40
   :header-rows: 1

   * - Class Attribute
     - Unit
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
     - (gm/cm3)
     - Brine density [Gas Saturated, Gas-Free Brine, Freshwater]
   * - .bVis
     - (cP)
     - Brine viscosity [Gas Saturated, Gas-Free Brine, Freshwater]
   * - .bVisblty
     - (1/Bar or 1/psi)
     - Gas-saturated brine viscosibility
   * - .bw
     - (rm3/sm3 or rb/stb)
     - Brine formation volume factor [Gas Saturated, Gas-Free Brine, Freshwater]
   * - .Rs
     - dict (sm3/sm3 or scf/stb)
     - Per-gas solution ratios, e.g. {'CO2': 15.2, 'CH4': 2.1}
   * - .Rs_total
     - (sm3/sm3 or scf/stb)
     - Total solution gas ratio (sum of all per-gas Rs)
   * - .Cf_usat
     - (1/Bar or 1/psi)
     - Brine undersaturated compressibility
   * - .Cf_sat
     - (1/Bar or 1/psi)
     - Brine saturated compressibility. Requires cw_sat input to be set True to calculate
   * - .gas_comp
     - dict
     - Normalized gas composition used (including estimated HC split from SG)
   * - .MwBrine
     - (g/mol)
     - Molecular weight of gas-free brine

Examples:

Pure CO2 case at 5000 psia x 275 deg F and 3% NaCl brine:

.. code-block:: python

    >>> from pyrestoolbox import brine
    >>> mix = brine.SoreideWhitson(pres=5000, temp=275, ppm=30000, y_CO2=1.0, metric=False)
    >>> mix.bDen  # Returns [Gas Saturated, Gas-Free Brine, Freshwater]
    [0.9733769457162755, 0.968164592979362, 0.9476497407774847]
    >>> mix.Rs  # Returns per-gas Rs dict (scf/stb)
    {'CO2': 140.90858294709142}
    >>> mix.bw  # Returns [Gas Saturated, Gas-Free, Freshwater]
    [1.0968991160573776, 1.0543023174291248, 1.0542033190822462]

Pure CH4 case (SG=0.554) at 5000 psia x 275 deg F and 3% NaCl brine:

.. code-block:: python

    >>> mix = brine.SoreideWhitson(pres=5000, temp=275, ppm=30000, y_CO2=0, sg=0.554, metric=False)
    >>> mix.Rs
    {'CH4': 21.414423540331008}
    >>> mix.bDen
    [0.9641137202631425, 0.968164592979362, 0.9476497407774847]

Mixed gas (10% CO2, 5% H2S, SG=0.7) at 200 Bar x 80 degC and 10,000 ppm NaCl:

.. code-block:: python

    >>> mix = brine.SoreideWhitson(pres=200, temp=80, ppm=10000, y_CO2=0.1, y_H2S=0.05, sg=0.7, metric=True)
    >>> mix.gas_comp  # Estimated gas composition including HC split
    {'CO2': 0.1, 'H2S': 0.05, 'CH4': 0.8133, 'C2H6': 0.0351, 'C3H8': 0.0015, 'nC4H10': 0.0001}
    >>> mix.Rs_total  # Total dissolved gas (sm3/sm3)
    6.32875877743837
    >>> mix.bDen
    [0.9854845215724698, 0.9871360082710434, 0.9804911502375318]

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
- Garcia, J.E., "Density of Aqueous Solutions of CO2", LBNL Report 49023, 2001.
- Plyasunov, A.V., "Values of the Apparent Molar Volumes...," Fluid Phase Equilibria, Parts I-IV, 2019-2021.
- Islam, A.W. and Carlson, E.S. (2012), "Viscosity Models and Effects of Dissolved CO2", Energy & Fuels, 26(8), 5330-5336.
- Ostermann, R.D. et al. (1985), "The Effect of Dissolved Gas on Reservoir Brine Viscosity", SPE 14211.

Density and Viscosity Calculation Details
----------------------

**Gas-corrected brine density** uses the Garcia (2001) Eq. 18 mixing rule, which is thermodynamically generic
and applicable to any dissolved gas or mixture of dissolved gases. For each dissolved gas species, the apparent
molar volume at infinite dilution V_phi(T,P) is computed using the Plyasunov (2019-2021) A12-infinity model
with IAPWS-IF97 Region 1 pure water properties. This provides rigorous T- and P-dependent V_phi for all 8
supported gas species (H2, N2, CH4, CO2, C2H6, C3H8, nC4H10, H2S).

For mixed dissolved gases, mole-fraction-weighted effective V_phi and MW are used:

.. code-block:: text

    V_phi_eff = sum(yi * V_phi_i) / sum(yi)
    MW_eff    = sum(yi * MW_i) / sum(yi)

    rho = (1 + x_gas * MW_eff / (MW_brine * x_brine)) / (x_gas * V_phi_eff / (MW_brine * x_brine) + 1/rho_brine)

where x_gas is the total dissolved gas mole fraction, x_brine = 1 - x_gas, and rho_brine is the gas-free
brine density (IAPWS-IF97 freshwater with Spivey salt correction).

Gases with MW/V_phi > 1 (CO2, H2S) increase brine density; gases with MW/V_phi < 1 (CH4, N2, H2, C2H6, C3H8)
decrease it.

**Gas-corrected brine viscosity** applies multiplicative per-gas corrections to the gas-free Mao-Duan (2009)
brine viscosity:

.. list-table:: Viscosity correction by gas
   :widths: 15 40 20
   :header-rows: 1

   * - Gas
     - Correction
     - Reference
   * - CO2
     - mu x (1 + 4.65 x xCO2^1.0134)
     - Islam-Carlson (2012)
   * - H2S
     - mu x (1 + 1.5 x xH2S^1.0134)
     - Murphy-Gaines (1974)
   * - CH4
     - mu x max(1.109 - 5.98e-4*T + 1.0933e-6*T^2, 1.0), T in degF
     - Ostermann (SPE-14211, 1985)
   * - Other
     - No correction (factor = 1.0)
     - Negligible effect at typical dissolved concentrations

For gas mixtures, the corrections are applied multiplicatively across all dissolved species.

