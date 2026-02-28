===================================
Gas PVT & Flow
===================================

Calculation Methods and Class Objects
=====================================
pyResToolBox uses class objects to track calculation options through the functions. Class objects can be set via strings or explicitly via object options

.. list-table:: Method Variables & Class Objects
   :widths: 10 15 40
   :header-rows: 1

   * - Class Variable
     - Class Object 
     - Class Description & Options
   * - zmethod
     - z_method
     - Method for calculating gas Z-Factor. Defaults to 'DAK'. 
       Options are:
        + 'DAK': Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al. - Slowest, Most Accurate
        + 'HY': Hall & Yarborough (1973) - Second Fastest
        + 'WYW': Wang, Ye & Wu (2021) - Fastest, Least Accurate
        + 'BUR'/'BNS': Tuned five component Peng Robinson EOS model, Burgoyne, Nielsen & Stanko (2025), `SPE-229932-MS <https://doi.org/10.2118/229932-MS>`_ - Fast, reliable, and able to handle wide range of mixture types including CO2, H2S, N2 and H2 at concentrations up to pure inerts. More information about the method can be found `here <https://github.com/mwburgoyne/5_Component_PengRobinson_Z-Factor>`_
   * - cmethod
     - c_method
     - Method for calculating gas critical properties. Defaults to 'PMC' 
       Options are:
        + 'SUT': Sutton with Wichert & Aziz non-hydrocarbon corrections
        + 'PMC': Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
        + 'BUR': Correlation tuned to return the critical properties of the pure hydrocarbon component of a mixture, required for the 'BUR' (tuned Peng Robinson) Z-Factor method. More information about the method can be found `here <https://github.com/mwburgoyne/5_Component_PengRobinson_Z-Factor>`_ 


Users can specify which calculation method to use either by passing an option string, or a class object to any given function. The implementation of class objects should make it easier to program in an IDE that supports type hinting

Examples:

Calculating gas Z-Factor of pure methane using DAK and PMC for critical properties

.. code-block:: python

    >>> from pyrestoolbox import gas
    >>> gas.gas_z(p=2350, sg=0.68, degf = 180, zmethod='DAK', cmethod='PMC')
    0.8785399927100872

.. note::

   **BNS method strongly recommended for high-inert or hydrogen-bearing gases.**
   When modelling gases with significant CO2 (>10%), H2S, N2, or any H2 content, the ``'BNS'`` Z-factor and critical property methods (tuned 5-component Peng-Robinson EOS) should be used. Standard correlations (DAK/HY/WYW with PMC/SUT) were developed for sweet natural gas and become unreliable at high impurity concentrations. The BNS method handles the full range from pure hydrocarbon to 100% inerts. When ``h2 > 0``, BNS is auto-selected. For pure or high-concentration CO2/H2S/N2 gases, explicitly set ``zmethod='BNS', cmethod='BNS'``.

Calculating gas Z-Factor of pure CO2

.. code-block:: python

    >>> gas.gas_z(p=2350, sg=0.68, degf = 180, co2 = 1.0, zmethod='BUR', cmethod='BUR')
    0.5258309021348752
    
Calculating gas sg, and then gas Z-Factor of a mixture of 5% CO2, 10% H2S, 0% N2 and 20% H2 (remainder natural gas with MW = 19).

.. code-block:: python

    >>> gsg = gas.gas_sg(hc_mw = 19.0, co2 = 0.05, h2s = 0.10, n2 = 0, h2 = 0.20)
    >>> gas.gas_z(p=2350, sg=gsg, degf = 180, co2 = 0.05, h2s = 0.10, n2 = 0, h2 = 0.20, zmethod='BUR', cmethod='BUR')
    0.9048153036714465


Unit System Support
===================
All gas module public functions accept an optional ``metric=False`` parameter. When ``metric=True``, inputs and outputs use Eclipse METRIC units (barsa, deg C, m, kg/m3, etc.) instead of oilfield units (psia, deg F, ft, lb/cuft, etc.). See individual function parameter tables for specific unit mappings.

.. note::

   **Standard conditions:** All "standard" volumes (scf, sm3, Bg, Rs, etc.) are calculated using oilfield standard conditions (60 deg F, 14.696 psia) regardless of the ``metric`` setting. This produces a ~0.06% systematic offset vs true Eclipse METRIC standard conditions (15 deg C, 1.01325 barsa) for gas volumes.


Function List
=============

.. list-table:: Gas Functions
   :widths: 15 40
   :header-rows: 1

   * - Task
     - Function
   * - Gas Tc & Pc Calculation
     - `pyrestoolbox.gas.gas_tc_pc`_  
   * - Gas Z-Factor Calculation
     - `pyrestoolbox.gas.gas_z`_
   * - Gas Viscosity or Viscosity * Z product
     - `pyrestoolbox.gas.gas_ug`_
   * - Gas Compressibility
     - `pyrestoolbox.gas.gas_cg`_
   * - Gas Formation Volume Factor
     - `pyrestoolbox.gas.gas_bg`_  
   * - Gas Density
     - `pyrestoolbox.gas.gas_den`_  
   * - Gas Water of Condensation
     - `pyrestoolbox.gas.gas_water_content`_
   * - Convert P/Z to P
     - `pyrestoolbox.gas.gas_ponz2p`_
   * - Convert Gas Gradient to SG
     - `pyrestoolbox.gas.gas_grad2sg`_
   * - Delta Pseudopressure
     - `pyrestoolbox.gas.gas_dmp`_
   * - Gas Condensate FWS SG
     - `pyrestoolbox.gas.gas_fws_sg`_
   * - Gas Flow Rate Radial
     - `pyrestoolbox.gas.gas_rate_radial`_  
   * - Gas Mixture SG
     - `pyrestoolbox.gas.gas_sg`_
   * - Gas Flow Rate Linear
     - `pyrestoolbox.gas.gas_rate_linear`_
  

pyrestoolbox.gas.gas_tc_pc
======================

.. code-block:: python

    gas_tc_pc(sg, co2 = 0, h2s = 0, n2 = 0, h2 = 0, cmethod = 'PMC', tc = 0, pc = 0, metric = False) -> tuple

Returns a tuple of critical temperature (deg R, or K if metric=True) and critical pressure (psia, or barsa if metric=True) for hydrocarbon gas. If one or both of the tc and pc parameters are set to be non-zero, then this function will return that unchanged value for the corresponding critical parameter.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - sg
     - float
     - Gas SG relative to air   
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined 
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to zero if undefined. If positive fraction, cmethod will override to 'BUR'
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - tc
     - float
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

Examples:

.. code-block:: python

    >>> gas.gas_tc_pc(sg=0.7, co2 = 0.15)
    (363.9387708314338, 738.3190067714969)

    >>> gas.gas_tc_pc(sg=0.7, co2 = 0.15, tc=365, cmethod='SUT')
    (365, 709.2356299485114)

Using metric units (returns Tc in K, Pc in barsa):

.. code-block:: python

    >>> gas.gas_tc_pc(sg=0.7, co2 = 0.15, metric=True)
    (202.18820601746322, 50.90644977899459)

pyrestoolbox.gas.gas_z
==================

.. code-block:: python

    gas_z(p, sg, degf, zmethod='DAK', cmethod='PMC', co2 = 0, h2s = 0, n2 = 0, h2 = 0, tc = 0, pc = 0, metric = False) -> float or np.array

Returns gas Z-factor (either float or Numpy array depending upon type of p specified) using specified method. 
A float or list / array can be used for p, returning corresponding 1-D array of Z-Factors. The cmethod will be used to calculate critical gas parameters unless tc and/or pc are explicitly set to be non-zero. This option enables users to use precalculate gas critical properties and so avoid repeated duplicated critical property calculations when compute time is an issue


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float, list or np.array
     - Gas pressure (psia, or barsa if metric=True)
   * - sg
     - float
     - Gas SG relative to air
   * - degf
     - float
     - Reservoir Temperature (deg F, or deg C if metric=True)
   * - zmethod
     - string or z_method
     - Method for calculating gas Z-factor. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to zero if undefined. zmethod and cmethod get overriden to 'BUR' if positive fraction.
   * - tc
     - float
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

Examples:

.. code-block:: python

    >>> gas.gas_z(p=1000, sg=0.75, degf=160, n2 = 0.02, co2 = 0.17)
    0.9138558878125714

    >>> gas.gas_z(p=1000, sg=0.75, degf=160, n2 = 0.02, co2 = 0.17, zmethod='HY')
    0.9142136711443208

    >>> gas.gas_z(p=[1000, 2000], sg=0.75, degf=160, cmethod='SUT', n2 = 0.02, co2 = 0.17)
    array([0.91900003, 0.87160514])

pyrestoolbox.gas.gas_ug
===================

.. code-block:: python

    gas_ug(p, sg, degf, zmethod ='DAK', cmethod = 'PMC', co2 = 0, h2s = 0, n2 = 0, h2 = 0, tc = 0, pc = 0, zee = 0, ugz = False, metric = False) -> float or np.array

Returns gas viscosity (cP) using Lee, Gonzalez & Eakin (1966) correlation unless the 'BUR' method for Z-Factor is selected in which case a tuned LBC method is used. 
A float or list / array can be used for p, returning corresponding 1-D array of gas viscosities. The cmethod will be used to calculate critical gas parameters unless tc and/or pc are explicitly set to be non-zero. This option enables users to use pre-calculate gas critical properties and so avoid repeated duplicated critical property calculations when compute time is an issue
Furnishing a positive value for zee means it will be used instead of calculating Z, and setting ugz to True will return the viscosity * Z product instead of viscosity alone


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float, list or np.array
     - Gas pressure (psia, or barsa if metric=True)
   * - sg
     - float
     - Gas SG relative to air
   * - degf
     - float
     - Reservoir Temperature (deg F, or deg C if metric=True)
   * - zmethod
     - string or z_method
     - Method for calculating gas Z-factor. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to zero if undefined. Overrides methods to 'BUR' if positive fraction.
   * - tc
     - float
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified
   * - zee
     - float, list or np.array
     - Z-Factor value, or list of values of same length as P, in case recalculation of Z-Factors is not needed. If undefined, will trigger Z-Factor calculation.
   * - ugz
     - boolean
     - Boolean flag that if True returns ug * Z instead of ug
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

Examples:

.. code-block:: python

    >>> gas.gas_ug(p=1000, sg=0.75, degf=180, zmethod ='HY', cmethod = 'SUT')
    0.014118890100250796

    >>> gas.gas_ug(p=1000, sg=0.75, degf=180)
    0.014110092961853301
    
    
pyrestoolbox.gas.gas_cg
===================

.. code-block:: python

    gas_cg(p, sg, degf, co2 = 0, h2s = 0, n2 = 0, h2 = 0, tc = 0, pc = 0, zmethod = 'DAK', cmethod ='PMC', metric = False) -> float or np.array

Returns gas compressibility (1/psi, or 1/barsa if metric=True).
A float or list / array can be used for p, returning corresponding 1-D array of gas compressibility's. The cmethod will be used to calculate critical gas parameters unless tc and/or pc are explicitly set to be non-zero. This option enables users to use precalculated gas critical properties and so avoid repeated duplicated critical property calculations when compute time is an issue


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float, list or np.array
     - Gas pressure (psia, or barsa if metric=True)
   * - sg
     - float
     - Gas SG relative to air
   * - degf
     - float
     - Reservoir Temperature (deg F, or deg C if metric=True)
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to zero if undefined. If positive fraction, cmethod will override to 'BUR'
   * - tc
     - float
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

Examples:

.. code-block:: python

    >>> gas.gas_cg(p=2000, sg=0.68, degf=120, co2=0.05)
    0.0005374854430839333

    >>> gas.gas_cg(p=np.array([1000,2000]), sg=0.68, degf=120, co2=0.05)
    array([0.00110369, 0.00053749])
    

pyrestoolbox.gas.gas_bg
===================

.. code-block:: python

    gas_bg(p, sg, degf, zmethod='DAK', cmethod = 'PMC', co2 = 0, h2s = 0, n2 = 0, h2 = 0, tc = 0, pc = 0, metric = False) -> float or np.array

Returns gas formation volume factor (rcf/scf, or rm3/sm3 if metric=True).
A float or list / array can be used for p, returning corresponding 1-D array of gas FVF's. The cmethod will be used to calculate critical gas parameters unless tc and/or pc are explicitly set to be non-zero. This option enables users to use precalculate gas critical properties and so avoid repeated duplicated critical property calculations when compute time is an issue.


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float, list or np.array
     - Gas pressure (psia, or barsa if metric=True)
   * - sg
     - float
     - Gas SG relative to air
   * - degf
     - float
     - Reservoir Temperature (deg F, or deg C if metric=True)
   * - zmethod
     - string or z_method
     - Method for calculating gas Z-factor. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to zero if undefined. If positive fraction, cmethod will override to 'BUR'
   * - tc
     - float
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

Examples:

.. code-block:: python

    >>> gas.gas_bg (p=3000, sg=0.78, degf=240)
    0.005927563975073749

    >>> 1 / gas.gas_bg (p=[3000, 5000], sg=0.78, degf=240)
    array([168.70336688, 249.54573283])

pyrestoolbox.gas.gas_den
=====================

.. code-block:: python

    gas_den(p, sg, degf, zmethod ='DAK', cmethod ='PMC', co2 = 0, h2s = 0, n2 = 0, h2 = 0, tc = 0, pc = 0, metric = False) -> float or np.array

Returns gas density (lb/cuft, or kg/m3 if metric=True).
A float or list / array can be used for p, returning corresponding 1-D array of gas densities. The cmethod will be used to calculate critical gas parameters unless tc and/or pc are explicitly set to be non-zero. This option enables users to use precalculate gas critical properties and so avoid repeated duplicated critical property calculations when compute time is an issue


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float, list or np.array
     - Gas pressure (psia, or barsa if metric=True)
   * - sg
     - float
     - Gas SG relative to air
   * - degf
     - float
     - Reservoir Temperature (deg F, or deg C if metric=True)
   * - zmethod
     - string or z_method
     - Method for calculating gas Z-factor. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to zero if undefined. If positive fraction, cmethod will override to 'BUR'
   * - tc
     - float
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

Examples:

.. code-block:: python

    >>> gas.gas_den (p=2000, sg=0.75, degf=150, zmethod ='HY', cmethod ='SUT', n2 = 0.02, co2 = 0.15, h2s = 0.02)
    7.736656004563576
    

pyrestoolbox.gas.gas_water_content
==============================

.. code-block:: python

    gas_water_content(p, degf, salinity=0, metric = False) -> float

Returns saturated volume of water vapor in natural gas (stb/mmscf). From 'PVT and Phase Behaviour Of Petroleum Reservoir Fluids' by Ali Danesh, with salinity correction.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float
     - Gas pressure (psia, or barsa if metric=True)
   * - degf
     - float
     - Reservoir Temperature (deg F, or deg C if metric=True)
   * - salinity
     - float
     - Water salinity (wt% NaCl). Defaults to 0 (freshwater). Higher salinity reduces water content
   * - metric
     - bool
     - If True, pressure in barsa, temperature in deg C. Defaults to False

Examples:

.. code-block:: python

    >>> gas.gas_water_content(p=1500, degf=165)
    0.6474226409378979

    >>> gas.gas_water_content(p=1500, degf=165, salinity=5)
    0.628635730743162

pyrestoolbox.gas.gas_ponz2p
=======================

.. code-block:: python

    gas_ponz2p(poverz, sg, degf, zmethod='DAK', cmethod='PMC', co2 = 0, h2s = 0, n2 = 0, h2 = 0, tc = 0, pc = 0, rtol = 1E-7, metric = False) -> float or np.array

Returns gas pressure (psia, or barsa if metric=True) corresponding to a value of P/Z, iteratively solving with specified zmethod via bisection.
A float or list / array can be used for poverz, returning corresponding 1-D array of pressures. The cmethod will be used to calculate critical gas parameters unless tc and/or pc are explicitly set to be non-zero. This option enables users to use precalculate gas critical properties and so avoid repeated duplicated critical property calculations when compute time is an issue


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - poverz
     - float, list or np.array
     - Gas pressure / Z-factor (psia, or barsa if metric=True)
   * - sg
     - float
     - Gas SG relative to air
   * - degf
     - float
     - Reservoir Temperature (deg F, or deg C if metric=True)
   * - zmethod
     - string or z_method
     - Method for calculating gas Z-factor. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to zero if undefined. If positive fraction, cmethod will override to 'BUR'
   * - tc
     - float
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified
   * - rtol
     - float
     - relative solution tolerance as compared with abs([User P/Z - Calculated P/Z] / [User P/Z])
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

Examples:

.. code-block:: python

    >>> gas.gas_ponz2p(poverz=2500, sg=0.75, degf=165)
    2081.5489292144775

    >>> gas.gas_ponz2p(poverz=[2500,5000], sg=0.75, degf=165)
    array([2081.54892921, 4856.97983205])
    
pyrestoolbox.gas.gas_grad2sg
========================

.. code-block:: python

    gas_grad2sg( grad, p, degf, zmethod='DAK', cmethod='PMC', co2 = 0, h2s = 0, n2 = 0, h2 = 0, tc = 0, pc = 0, rtol = 1E-7, metric = False) -> float

Returns gas specific gravity consistent with observed gas gradient. Calculated through iterative solution method. Will fail if gas SG is below 0.55, or greater than 1.75

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - grad
     - float
     - Observed gas gradient (psi/ft, or bar/m if metric=True)
   * - p
     - float, list or np.array
     - Pressure at observation (psia, or barsa if metric=True)
   * - degf
     - float
     - Reservoir Temperature (deg F, or deg C if metric=True)
   * - zmethod
     - string or z_method
     - Method for calculating gas Z-factor. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to zero if undefined. If positive fraction, cmethod will override to 'BUR'
   * - tc
     - float
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified
   * - rtol
     - float
     - relative solution tolerance as compared with abs([User grad - Calculated grad] / [User grad])
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

Examples:

.. code-block:: python

    >>> gas.gas_grad2sg(grad=0.0657, p=2500, degf=175)
    0.7495803994806547
    

pyrestoolbox.gas.gas_dmp
=====================

.. code-block:: python

    gas_dmp(p1, p2, degf, sg, zmethod='DAK', cmethod = 'PMC', co2 = 0, h2s = 0, n2 = 0, h2 = 0, tc = 0, pc = 0, metric = False) -> float

Returns gas pseudo-pressure integral (psi2/cP, or bar2/cP if metric=True) between two pressure points. Will return a positive value if p1 < p2, and a negative value if p1 > p2.
Integrates the equation: m(p) = 2 * p / (ug * z)

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p1
     - float, list or np.array
     - First gas pressure (psia, or barsa if metric=True)
   * - p2
     - float, list or np.array
     - Second gas pressure (psia, or barsa if metric=True)
   * - sg
     - float
     - Gas SG relative to air.
   * - degf
     - float
     - Reservoir Temperature (deg F, or deg C if metric=True)
   * - zmethod
     - string or z_method
     - Method for calculating gas Z-factor. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to zero if undefined. If positive fraction, cmethod will override to 'BUR'
   * - tc
     - float
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

Examples:

.. code-block:: python

    >>> gas.gas_dmp(p1=1000, p2=2000, degf=185, sg=0.78, zmethod='HY', cmethod = 'SUT', n2 = 0.05, co2 = 0.1, h2s = 0.02)
    213690308.9907268

    >>> gas.gas_dmp(p1=2000, p2=1000, degf=185, sg=0.78, tc = 371, pc = 682)
    -213713909.36339885
        
pyrestoolbox.gas.gas_fws_sg
=======================

.. code-block:: python

    gas_fws_sg(sg_g, cgr, api_st, metric = False) -> float

Estimates specific gravity of full-wellstream (FWS) gas from gas-condensate well. Calculates from weighted average surface gas SG, CGR and API. Uses Standing correlation to estimate condensate MW from API.
Returns SG of FWS gas

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - sg_g
     - float
     - Specific gravity of weighted average surface gas (relative to air)
   * - cgr
     - float
     - Condensate gas ratio (stb/mmscf, or sm3/sm3 if metric=True)
   * - api_st
     - float
     - Density of stock tank liquid (API)
   * - metric
     - bool
     - If True, CGR input in sm3/sm3. Defaults to False

Examples:

.. code-block:: python

    >>> gas.gas_fws_sg(sg_g=0.855, cgr=30, api_st=53)
    0.937116010334538
    
    
pyrestoolbox.gas.gas_rate_radial
======================

.. code-block:: python

    gas_rate_radial(k, h, pr, pwf, r_w, r_ext, degf, zmethod='DAK, cmethod='PMC', S = 0, D = 0, sg = 0.75, n2 = 0, co2 = 0, h2s = 0, tc  = 0, pc = 0, metric = False) -> float or np.array

Returns gas rate (Mscf/day, or sm3/d if metric=True) for radial flow using Darcy pseudo steady state equation & gas pseudopressure.
Arrays can be used for any one of k, h, pr or pwf, returning corresponding 1-D array of rates. Using more than one input array -- while not prohibited -- will not return expected results.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - k
     - float, list or np.array
     - Effective permeability to gas flow (mD)
   * - h
     - float, list or np.array
     - Net height for flow (ft, or m if metric=True)
   * - pr
     - float, list or np.array
     - Reservoir pressure (psia, or barsa if metric=True)
   * - pwf
     - float, list or np.array
     - BHFP (psia, or barsa if metric=True)
   * - r_w
     - float
     - Wellbore Radius (ft, or m if metric=True)
   * - r_ext
     - float
     - External Reservoir Radius (ft, or m if metric=True)
   * - degf
     - float
     - Reservoir Temperature (deg F, or deg C if metric=True)
   * - zmethod
     - string or z_method
     - Method for calculating gas Z-factor. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - tc
     - float
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - S
     - float
     - Skin. Defaults to zero if undefined
   * - D
     - float
     - Non Darcy Skin Factor (day/Mscf, or day/sm3 if metric=True). Defaults to zero if undefined
   * - sg
     - float
     - Gas SG relative to air, Defaults to 0.75 if undefined
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

Examples:

.. code-block:: python

    >>> gas.gas_rate_radial(k=5, h=50, pr=2000, pwf=750, r_w=0.3, r_ext=1500, degf=180, sg = 0.75, D = 0.01, S=5)
    2078.9101970773477

    >>> gas.gas_rate_radial(k=1, h=50, pr=[2000,1000], pwf=750, r_w=0.3, r_ext=1500, degf=180, sg = 0.75, D = 0.01, S=5)
    array([704.29202227, 135.05317439])
    

pyrestoolbox.gas.gas_sg
======================

.. code-block:: python

    gas_sg(hc_mw, co2, h2s, n2, h2) -> float

Returns specific gravity of a gas mixture (relative to air) given the hydrocarbon molecular weight and molar fractions of non-hydrocarbon components. Useful for converting a known gas composition into the SG input required by other gas functions.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - hc_mw
     - float
     - Molecular weight of the hydrocarbon portion of the gas (g/gmol). Typical range 16 (pure C1) to ~30+ (wet gas)
   * - co2
     - float
     - Molar fraction of CO2 in the total gas mixture
   * - h2s
     - float
     - Molar fraction of H2S in the total gas mixture
   * - n2
     - float
     - Molar fraction of Nitrogen in the total gas mixture
   * - h2
     - float
     - Molar fraction of Hydrogen in the total gas mixture

Examples:

.. code-block:: python

    >>> gas.gas_sg(hc_mw=19.0, co2=0.05, h2s=0.10, n2=0, h2=0.20)
    0.6338246461857093

    >>> gas.gas_sg(hc_mw=16.043, co2=0, h2s=0, n2=0, h2=0)
    0.5537797721781152


pyrestoolbox.gas.gas_rate_linear
======================

.. code-block:: python

    gas_rate_linear(k, pr, pwf, area, length, degf, zmethod='DAK, cmethod='PMC', sg = 0.75, n2 = 0, co2 = 0, h2s = 0, tc  = 0, pc = 0, metric = False) -> float or np.array

Returns gas rate (Mscf/day, or sm3/d if metric=True) for linear flow using Darcy steady state equation & gas pseudopressure.
Arrays can be used for any one of k, pr, pwf or area, returning corresponding 1-D array of rates. Using more than one input array -- while not prohibited -- will not return expected results.


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - k
     - float, list or np.array
     - Effective permeability to gas flow (mD)
   * - pr
     - float, list or np.array
     - Reservoir pressure (psia, or barsa if metric=True)
   * - pwf
     - float, list or np.array
     - BHFP (psia, or barsa if metric=True)
   * - area
     - float
     - Net cross-sectional area perpendicular to direction of flow (ft2, or m2 if metric=True)
   * - length
     - float
     - Linear distance of fluid flow (ft, or m if metric=True)
   * - degf
     - float
     - Reservoir Temperature (deg F, or deg C if metric=True)
   * - zmethod
     - string or z_method
     - Method for calculating gas Z-factor. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - tc
     - float
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - sg
     - float
     - Gas SG relative to air, Defaults to 0.75 if undefined
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

Examples:

.. code-block:: python

    >>> gas.gas_rate_linear(k=0.1, area=50, length=200, pr=2000, pwf=250, degf=180, sg = 0.8)
    8.202200317597859

    >>> gas.gas_rate_linear(k=0.1, area=50, length=200, pr=[2000, 1000, 500], pwf=250, degf=180, sg = 0.8)
    array([8.20220032, 2.10691337, 0.42685002])
    
