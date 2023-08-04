===================================
CO2 Saturated Brine PVT
===================================

Returns a class object with calculated CO2 saturated brine property attributes with the following methods

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
     
pyrestoolbox.CO2_Brine_Mixture
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
     - Parts per million (Wt) NaCl in the brine (1% by Wt = 10,000 ppm)
   * - metric
     - Boolean
     - If true, treats input pressure & temperature as metric, otherwise treats as Field units.
     
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
     - Brine density [CO2 Saturated, Pure Brine].
   * - .bVis
     - (cP)
     - Brine viscosity [CO2 Saturated, Pure Brine]
   * - .bVisblty
     - (1/Bar or 1/psi)
     - CO2 Saturated brine viscosibility.
   * - .bw
     - (rm3/smr or rb/stb)
     - Brine formation volume factor  [CO2 Saturated, Pure Brine]
   * - .Rs
     - (sm3/sm3 or scf/stb)
     - CO2 Saturated Brine solution gas ratio
   * - .Cf_usat
     - (1/Bar or 1/psi)
     - Brine undersaturated compressibility 

.. list-table:: Optional class functions
   :widths: 25 40
   :header-rows: 1

   * - Function
     - Description
   * - .ezrokhi(lower_degC, upper_degC)
     - Calculate Ezrokhi coefficients for effects of dissolved CO2 on brine density and viscosity

Calculates and populates the following additional attributes;

.EzrokhiDenA : List of A_CO2's for equation Ai(T) = A[0] + A[1] * degC + A[2] * degC**2, for Ezrokhi density adjustment

.EzrokhiVisB : List of B_CO2's for equation Bi(T) = B[0] + B[1] * degC + B[2] * degC**2, for Ezrokhi viscosity adjustment
                
Examples:

Usage example for 5000 psia x 275 deg F and 3% NaCl brine:

.. code-block:: python

    >>> from pyrestoolbox import pyrestoolbox as rtb
    >>> mix = rtb.CO2_Brine_Mixture(pres = 5000, temp = 275, ppm = 30000, metric = False)
    >>> mix.bw  # Returns [CO2 Saturated Brine Bw, Pure Brine Bw]
    [1.1085535365828796, 1.0543051557116332]
    >>> mix.x  # Returns molar fractions in aqueous phase [xCO2, xH2O]
    array([0.02431431, 0.96633611])
    
Usage example for 175 Bara x 85 degC and 0% NaCl brine:

.. code-block:: python

    >>> mix = rtb.CO2_Brine_Mixture(pres = 175, temp = 85)
    >>> mix.Rs  # Returns sm3 dissolved CO2 / sm3 Brine
    24.502168045494223   
    >>> mix.ezrokhi(lower_degC = 70, upper_degC = 140)
    >>> mix.EzrokhiDenA, mix.EzrokhiVisB
    ([0.05812057460992598, 0.0007985131535236179, -5.863680161690182e-06], [0.7819371846535907, 0, 0])
