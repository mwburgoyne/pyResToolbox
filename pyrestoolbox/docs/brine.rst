===================================
Brine PVT with differing degrees of methane or CO2 saturation
===================================

For Brine:Methane (brine_props) - Calculates Brine properties from modified Spivey Correlation per McCain Petroleum Reservoir Fluid Properties pg 160. Includes effect of user specified salt concentration and degree of methane saturation.
Returns tuple of (Bw (rb/stb), Density (sg), viscosity (cP), Compressibility (1/psi), Rw GOR (scf/stb))

For Brine:CO2 - Returns a class object with calculated CO2 saturated brine property attributes with the following methods

.. list-table:: Method employed for different calculations
   :widths: 30 40
   :header-rows: 1

   * - Property
     - Calculation method
   * - Mutual Solubility between CO2 and Brine
     - Spycher & Pruess (2010), modified SRK Cubic EOS method
   * - Pure Brine Density
     - Spivey et al. (modified), per "Petroleum Reservoir Fluid Property Correlations", (McCain, Spivey & Lenn: Chapter 4)
   * - CO2 Corrected Brine Density
     - Molar volume of dissolved CO2 estimated with Garcia (2001) equation, used with xCO2 calculated from Spycher & Pruess, and CO2-free brine density from Spivey et al to calculate insitu density
   * - Pure Brine viscosity
     - Mao-Duan (2009).
   * - CO2 Corrected Brine Viscosity
     - "Viscosity Models and Effects of Dissolved CO2", Islam-Carlson (2012) to adjust the pure brine viscosity for xCO2 calculated from Spycher & Pruess.  
     

pyrestoolbox.brine.brine_props
======================

.. code-block:: python

    brine_props(p, degf, wt, ch4_sat) -> tuple

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float
     - Pressure (psia)
   * - degf
     - float
     - Temperature (deg F)
   * - wt
     - float
     - Salt weight% in the brine (0 - 100)
   * - ch4_sat
     - float
     - Degree of methane saturation (0 - 1). 0 = No Methane, 1 = 100% Methane saturated

Examples:

.. code-block:: python

    >>> from pyrestoolbox import brine
    >>> bw, lsg, visw, cw, rsw = brine.brine_props(p=160, degf=135, wt=1.5, ch4_sat=1.0)
    >>> print('Bw:', bw)
    >>> print('SGw:', lsg)
    >>> print('Visw:', visw)
    >>> print('Cw:', cw)
    >>> print('Rsw:', rsw)
    Bw: 1.0151710978322923
    SGw: 0.9950036658248123
    Visw: 0.4993957925685796
    Cw: 0.00015465633691558178
    Rsw: 1.2549011339427625

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
    [1.1085795290443725, 1.0543051245909865, 1.0542061001251017]
    >>> mix.x  # Returns molar fractions in aqueous phase [xCO2, xH2O]
    array([0.02431245, 0.95743175])
    
Usage example for 175 Bara x 85 degC and 0% NaCl brine:

.. code-block:: python

    >>> mix = brine.CO2_Brine_Mixture(pres = 175, temp = 85)
    >>> mix.Rs  # Returns sm3 dissolved CO2 / sm3 Brine
    24.742923469934272   

