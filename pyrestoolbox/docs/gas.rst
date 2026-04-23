===================================
Gas PVT & Flow
===================================

Gas property calculations for hydrocarbon mixtures with optional impurity components (CO2, H2S, N2, H2). Includes Z-factor, viscosity, density, compressibility, formation volume factor, pseudopressure, gas flow rates (radial and linear), water content, and hydrate formation prediction with inhibitor dosing.

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
   * - hydmethod
     - hyd_method
     - Method for calculating gas hydrate formation temperature. Defaults to 'TOWLER'.
       Options are:
        + 'TOWLER': Towler & Mokhatab (2005) — T-explicit, gas gravity based. Hydrocarbon Processing 84, pp 61-62
        + 'MOTIEE': Motiee (1991) — T-explicit, gas gravity based. Hydrocarbon Processing 70, pp 98-99
   * - inhibitor_type
     - inhibitor
     - Thermodynamic hydrate inhibitor type. Optional (None = no inhibitor).
       Options are:
        + 'MEOH': Methanol
        + 'MEG': Monoethylene Glycol
        + 'DEG': Diethylene Glycol
        + 'TEG': Triethylene Glycol
        + 'ETOH': Ethanol


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

.. note::

   **BNS coupling and user-supplied Tc/Pc.**
   ``zmethod`` and ``cmethod`` are coupled for the BNS framework: if either is BNS, both are forced to BNS (a ``UserWarning`` is emitted if this overrules a user-chosen non-BNS method). Non-BNS methods are not coupled.

   User-supplied ``tc`` and ``pc`` are always respected, but their meaning depends on the critical-property method:

   - **SUT / PMC**: ``tc`` and ``pc`` are the *mixture* pseudo-critical values. They replace the correlation output entirely.
   - **BNS**: ``tc`` and ``pc`` are the *inert-free hydrocarbon* pseudo-critical values. They replace only the hydrocarbon pseudo-component in the 5-component EOS; inert Tc/Pc (CO2, H2S, N2, H2) remain the BNS internal per-component constants.

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
   * - Gas Hydrate Prediction
     - `pyrestoolbox.gas.gas_hydrate`_
   * - Forchheimer HVF Coefficient (β)
     - `pyrestoolbox.gas.gas_hvf_beta`_
   * - Non-Darcy (HVF) Skin
     - `pyrestoolbox.gas.gas_non_darcy_skin`_
   * - Partial-Penetration Pseudoskin
     - `pyrestoolbox.gas.gas_partial_penetration_skin`_
   * - Gas PVT Wrapper
     - `pyrestoolbox.gas.GasPVT`_


pyrestoolbox.gas.gas_tc_pc
======================

.. code-block:: python

    gas_tc_pc(sg, co2 = 0, h2s = 0, n2 = 0, h2 = 0, cmethod = 'PMC', tc = 0, pc = 0, metric = False) -> tuple

Returns a tuple of critical temperature (deg R, or K if metric=True) and critical pressure (psia, or barsa if metric=True). If one or both of the ``tc`` and ``pc`` parameters are set to be non-zero, then this function will return that unchanged value for the corresponding critical parameter.

For ``cmethod='SUT'`` and ``cmethod='PMC'``, the returned values describe *mixture* pseudo-critical properties (full gas, inerts included).

For ``cmethod='BNS'``, the returned values describe the *inert-free hydrocarbon* pseudo-critical properties only. The supplied inert fractions (``co2``, ``h2s``, ``n2``, ``h2``) are used solely to back out the inert-free hydrocarbon specific gravity from the overall mixture SG; the BNS pseudo-critical correlation is then applied to that hydrocarbon SG. Inert species (CO2, H2S, N2, H2) carry their own per-component Tc/Pc internally within the BNS 5-component PR-EOS, and are *not* included in the value returned here.

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
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

.. list-table:: Returns (tuple)
   :widths: 10 15 40
   :header-rows: 1

   * - Index
     - Type
     - Description
   * - [0]
     - float
     - Critical temperature (deg R, or K if metric=True)
   * - [1]
     - float
     - Critical pressure (psia, or barsa if metric=True)

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
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float or np.ndarray
     - Gas Z-factor (dimensionless). Returns same type as input p

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
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
   * - zee
     - float, list or np.array
     - Z-Factor value, or list of values of same length as P, in case recalculation of Z-Factors is not needed. If undefined, will trigger Z-Factor calculation.
   * - ugz
     - boolean
     - Boolean flag that if True returns ug * Z instead of ug
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float or np.ndarray
     - Gas viscosity (cP), or viscosity * Z product if ugz=True. Returns same type as input p

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
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float or np.ndarray
     - Gas compressibility (1/psi, or 1/barsa if metric=True). Returns same type as input p

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
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float or np.ndarray
     - Gas formation volume factor (rcf/scf, or rm3/sm3 if metric=True). Returns same type as input p

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
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float or np.ndarray
     - Gas density (lb/cuft, or kg/m3 if metric=True). Returns same type as input p

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

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float
     - Saturated water vapor content (stb/mmscf)

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
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
   * - rtol
     - float
     - relative solution tolerance as compared with abs([User P/Z - Calculated P/Z] / [User P/Z])
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float or np.ndarray
     - Gas pressure (psia, or barsa if metric=True). Returns same type as input poverz

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

Returns gas specific gravity consistent with observed gas gradient. Calculated through iterative solution method. Bisection bounds span pure H2 (SG ~0.070) to 3.0 to accommodate H2-blend and CO2-rich compositions; solutions outside this range will fail.

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
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
   * - rtol
     - float
     - relative solution tolerance as compared with abs([User grad - Calculated grad] / [User grad])
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float
     - Gas specific gravity (relative to air)

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
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float
     - Pseudo-pressure integral (psi2/cP, or bar2/cP if metric=True)

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

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float
     - Full wellstream gas specific gravity (relative to air)

Examples:

.. code-block:: python

    >>> gas.gas_fws_sg(sg_g=0.855, cgr=30, api_st=53)
    0.937116010334538
    
    
pyrestoolbox.gas.gas_rate_radial
======================

.. code-block:: python

    gas_rate_radial(k, h, pr, pwf, r_w, r_ext, degf, zmethod='DAK', cmethod='PMC', S=0, D=0, sg=0.75, co2=0, h2s=0, n2=0, h2=0, tc=0, pc=0, gas_pvt=None, metric=False) -> float or np.array

Returns gas rate (Mscf/day, or sm3/d if metric=True) for radial flow using Darcy pseudo steady state equation & gas pseudopressure.
Arrays can be used for any one of k, h, pr or pwf, returning corresponding 1-D array of rates. Using more than one input array -- while not prohibited -- will not return expected results.

A ``gas_pvt`` object can be provided instead of individual gas composition parameters. When ``gas_pvt`` is used, sg, co2, h2s, n2, h2, zmethod, cmethod, tc, and pc are extracted from the object.

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
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
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
   * - gas_pvt
     - GasPVT
     - GasPVT object. If provided, sg/co2/h2s/n2/h2/zmethod/cmethod/tc/pc are extracted from it
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float or np.ndarray
     - Gas rate (Mscf/day, or sm3/d if metric=True). Returns same type as the array input

Examples:

.. code-block:: python

    >>> gas.gas_rate_radial(k=5, h=50, pr=2000, pwf=750, r_w=0.3, r_ext=1500, degf=180, sg = 0.75, D = 0.01, S=5)
    2078.9101970773477

    >>> gas.gas_rate_radial(k=1, h=50, pr=[2000,1000], pwf=750, r_w=0.3, r_ext=1500, degf=180, sg = 0.75, D = 0.01, S=5)
    array([704.29202227, 135.05317439])

Using a GasPVT object:

.. code-block:: python

    >>> gpvt = gas.GasPVT(sg=0.75, co2=0.05)
    >>> gas.gas_rate_radial(k=5, h=50, pr=2000, pwf=750, r_w=0.3, r_ext=1500, degf=180, gas_pvt=gpvt, S=5, D=0.01)
    2072.675775394653
    

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

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float
     - Gas mixture specific gravity (relative to air)

Examples:

.. code-block:: python

    >>> gas.gas_sg(hc_mw=19.0, co2=0.05, h2s=0.10, n2=0, h2=0.20)
    0.6338246461857093

    >>> gas.gas_sg(hc_mw=16.043, co2=0, h2s=0, n2=0, h2=0)
    0.5537797721781152


pyrestoolbox.gas.gas_rate_linear
======================

.. code-block:: python

    gas_rate_linear(k, pr, pwf, area, length, degf, zmethod='DAK', cmethod='PMC', sg=0.75, co2=0, h2s=0, n2=0, h2=0, tc=0, pc=0, gas_pvt=None, metric=False) -> float or np.array

Returns gas rate (Mscf/day, or sm3/d if metric=True) for linear flow using Darcy steady state equation & gas pseudopressure.
Arrays can be used for any one of k, pr, pwf or area, returning corresponding 1-D array of rates. Using more than one input array -- while not prohibited -- will not return expected results.

A ``gas_pvt`` object can be provided instead of individual gas composition parameters.

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
     - Critical gas temperature (deg R, or K if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Tc (inert Tc stay at BNS internal constants)
   * - pc
     - float
     - Critical gas pressure (psia, or barsa if metric=True). Uses cmethod correlation if not specified. For BNS, overrides only the hydrocarbon pseudo-component Pc (inert Pc stay at BNS internal constants)
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
   * - gas_pvt
     - GasPVT
     - GasPVT object. If provided, sg/co2/h2s/n2/h2/zmethod/cmethod/tc/pc are extracted from it
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units. Defaults to False

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float or np.ndarray
     - Gas rate (Mscf/day, or sm3/d if metric=True). Returns same type as the array input

Examples:

.. code-block:: python

    >>> gas.gas_rate_linear(k=0.1, area=50, length=200, pr=2000, pwf=250, degf=180, sg = 0.8)
    8.202200317597859

    >>> gas.gas_rate_linear(k=0.1, area=50, length=200, pr=[2000, 1000, 500], pwf=250, degf=180, sg = 0.8)
    array([8.20220032, 2.10691337, 0.42685002])

Using a GasPVT object:

.. code-block:: python

    >>> gpvt = gas.GasPVT(sg=0.8)
    >>> gas.gas_rate_linear(k=0.1, area=50, length=200, pr=2000, pwf=250, degf=180, gas_pvt=gpvt)
    8.202199932859799


pyrestoolbox.gas.gas_hydrate
=============================

.. code-block:: python

    gas_hydrate(p, degf, sg, hydmethod='TOWLER', inhibitor_type=None, inhibitor_wt_pct=0,
                co2=0, h2s=0, n2=0, h2=0, p_res=None, degf_res=None,
                additional_water=0, metric=False) -> HydrateResult

Returns a ``HydrateResult`` dataclass with gas hydrate formation temperature (HFT), hydrate formation pressure (HFP), subcooling, hydrate window assessment, thermodynamic inhibitor calculations, a full water balance between reservoir and operating conditions, and inhibitor injection rates.

Two HFT correlations are available: Motiee (1991) and Towler & Mokhatab (2005). Hydrate formation pressure is computed by bisection inversion of the HFT correlation. Inhibitor temperature depression uses the Østergaard et al. (2005) cubic polynomial, and required inhibitor concentration is computed by Newton-Raphson inversion of the same polynomial. The required concentration is capped at the physical maximum for each inhibitor type.

**Water balance.** The gas leaves the reservoir saturated with vaporized water at reservoir P,T (``p_res``, ``degf_res``). At the operating point (lower P,T), the gas can hold less water vapor — the excess condenses as liquid. The function computes vaporized water at both conditions and reports the condensed amount. Any free liquid water entrained from the reservoir (``additional_water``) is added to the condensed water to give the total liquid water at the operating point. When gas composition is provided (``co2``/``h2s``/``n2``/``h2``), the SoreideWhitson VLE model is used; otherwise the Danesh correlation. If ``p_res``/``degf_res`` are not provided, the operating ``p``/``degf`` are used for both (no condensation).

**Inhibitor injection rate** is calculated from the **total liquid water** at the operating point (condensed + free water) and the required inhibitor concentration. Only liquid water needs treatment — vaporized water does not form hydrates. The ``inhibitor_wt_pct`` parameter represents the concentration of inhibitor in the aqueous phase (water + inhibitor mixture), not the mass fraction of the total stream. Injection rates are reported in lb/MMscf and gal/MMscf (oilfield) or kg/sm3 and L/sm3 (metric).

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float
     - Operating pressure at hydrate assessment point, e.g. wellhead (psia | barsa)
   * - degf
     - float
     - Operating temperature at hydrate assessment point (deg F | deg C)
   * - sg
     - float
     - Gas specific gravity (air = 1.0)
   * - hydmethod
     - string or hyd_method
     - Hydrate formation correlation. 'TOWLER' (default) or 'MOTIEE'. `Calculation Methods and Class Objects`_.
   * - inhibitor_type
     - string or inhibitor
     - Thermodynamic hydrate inhibitor. 'MEOH', 'MEG', 'DEG', 'TEG', 'ETOH', or None (default = no inhibitor). `Calculation Methods and Class Objects`_.
   * - inhibitor_wt_pct
     - float
     - Weight percent of inhibitor in aqueous phase (0-100). Defaults to 0
   * - co2
     - float
     - CO2 mole fraction (0-1). Enables composition-aware water content via SoreideWhitson. Defaults to 0
   * - h2s
     - float
     - H2S mole fraction (0-1). Defaults to 0
   * - n2
     - float
     - N2 mole fraction (0-1). Defaults to 0
   * - h2
     - float
     - H2 mole fraction (0-1). Defaults to 0
   * - p_res
     - float or None
     - Reservoir pressure where gas equilibrated with water (psia | barsa). Controls water content. If None, uses ``p``. Defaults to None
   * - degf_res
     - float or None
     - Reservoir temperature where gas equilibrated with water (deg F | deg C). Controls water content. If None, uses ``degf``. Defaults to None
   * - additional_water
     - float
     - Free liquid water entrained in the gas stream from the reservoir, e.g. mobile formation water (stb/MMscf | sm3/sm3). Added to condensed water for inhibitor dosing. Defaults to 0
   * - metric
     - bool
     - If True, inputs/outputs use Eclipse METRIC units (barsa, deg C). Defaults to False

.. list-table:: Returns (HydrateResult)
   :widths: 10 15 40
   :header-rows: 1

   * - Attribute
     - Type
     - Description
   * - hft
     - float
     - Hydrate formation temperature at operating pressure (deg F | deg C)
   * - hfp
     - float
     - Hydrate formation pressure at operating temperature (psia | barsa)
   * - subcooling
     - float
     - HFT minus operating temperature (deg F | deg C delta). Positive = in hydrate window
   * - in_hydrate_window
     - bool
     - True if operating temperature is below HFT
   * - inhibited_hft
     - float
     - HFT after inhibitor depression (deg F | deg C), or NaN if no inhibitor specified
   * - inhibitor_depression
     - float
     - Temperature depression from inhibitor (deg F | deg C delta), or 0
   * - required_inhibitor_wt_pct
     - float
     - Wt% inhibitor in aqueous phase needed to bring HFT below operating temperature, capped at physical maximum for selected inhibitor. 0 if no inhibitor or outside hydrate window
   * - max_inhibitor_wt_pct
     - float
     - Maximum valid wt% for selected inhibitor type (MEOH: 25, MEG: 70, DEG: 70, TEG: 50, ETOH: 30). 0 if no inhibitor specified
   * - inhibitor_underdosed
     - bool
     - True if ``required_inhibitor_wt_pct`` exceeds ``max_inhibitor_wt_pct`` for the selected inhibitor type — meaning this inhibitor **cannot** provide sufficient depression even at its physical maximum concentration. This does NOT indicate whether the *applied* ``inhibitor_wt_pct`` is sufficient; compare ``inhibited_hft`` to the operating temperature to determine if the applied dose provides protection
   * - water_vaporized_res
     - float
     - Vaporized water content at reservoir P,T (stb/MMscf | sm3/sm3). This is the water the gas picked up in the reservoir
   * - water_vaporized_op
     - float
     - Vaporized water content at operating P,T (stb/MMscf | sm3/sm3). This is the water the gas can still hold at the assessment point
   * - water_condensed
     - float
     - Water that condensed from vapor between reservoir and operating conditions (stb/MMscf | sm3/sm3). Equals max(water_vaporized_res - water_vaporized_op, 0)
   * - free_water
     - float
     - Free liquid water entrained from the reservoir (stb/MMscf | sm3/sm3). Equals the ``additional_water`` input
   * - total_liquid_water
     - float
     - Total liquid water at operating point: condensed + free water (stb/MMscf | sm3/sm3). This is the water that needs inhibitor treatment
   * - inhibitor_mass_rate
     - float
     - Required inhibitor mass injection rate based on total liquid water (lb/MMscf | kg/sm3). 0 if outside hydrate window or no inhibitor
   * - inhibitor_vol_rate
     - float
     - Required inhibitor volume injection rate based on total liquid water (gal/MMscf | L/sm3). 0 if outside hydrate window or no inhibitor

Examples:

.. code-block:: python

    >>> from pyrestoolbox import gas
    >>> r = gas.gas_hydrate(p=1000, degf=60, sg=0.65)
    >>> r.hft
    62.918902535978695
    >>> r.hfp
    814.0829443359373
    >>> r.subcooling
    2.9189025359786953
    >>> r.in_hydrate_window
    True

Using Motiee method:

.. code-block:: python

    >>> r = gas.gas_hydrate(p=1000, degf=60, sg=0.65, hydmethod='MOTIEE')
    >>> r.hft
    96.20360674625366

With MEG inhibitor (25 wt%):

.. code-block:: python

    >>> r = gas.gas_hydrate(p=2000, degf=80, sg=0.7, hydmethod='MOTIEE', inhibitor_type='MEG', inhibitor_wt_pct=25)
    >>> r.inhibited_hft
    97.97592063045367
    >>> r.inhibitor_depression
    11.0109375
    >>> r.required_inhibitor_wt_pct
    58.27975454526845

Using metric units (barsa, deg C):

.. code-block:: python

    >>> r = gas.gas_hydrate(p=100, degf=20, sg=0.7, metric=True)
    >>> r.hft
    21.017623569244456
    >>> r.hfp
    87.78593833339396

MEOH inhibitor with capping and injection rate (reservoir P,T specified, MEOH max = 25 wt%):

.. code-block:: python

    >>> r = gas.gas_hydrate(p=2000, degf=60, sg=0.7, hydmethod='MOTIEE', inhibitor_type='MEOH',
    ...                      p_res=4000, degf_res=250)
    >>> r.inhibitor_underdosed
    True
    >>> r.required_inhibitor_wt_pct
    25.0
    >>> r.max_inhibitor_wt_pct
    25.0
    >>> r.water_condensed
    1.6171304990445141
    >>> r.inhibitor_mass_rate
    188.77303358846294
    >>> r.inhibitor_vol_rate
    28.596710738497325

With CO2 composition and reservoir P,T (SoreideWhitson water content):

.. code-block:: python

    >>> r = gas.gas_hydrate(p=1000, degf=60, sg=0.65, hydmethod='MOTIEE', inhibitor_type='MEG', co2=0.05,
    ...                      p_res=3000, degf_res=200)
    >>> r.water_vaporized_res
    0.9105104447421333
    >>> r.water_condensed
    0.8606277674562288
    >>> r.required_inhibitor_wt_pct
    68.1447380066354
    >>> r.inhibitor_underdosed
    False

With reservoir P,T (gas equilibrated at reservoir, hydrate assessment at wellhead):

.. code-block:: python

    >>> r = gas.gas_hydrate(p=500, degf=40, sg=0.7, hydmethod='MOTIEE', inhibitor_type='MEG',
    ...                      p_res=4000, degf_res=250)
    >>> r.water_vaporized_res
    1.651022101177945
    >>> r.inhibitor_mass_rate
    1314.312441059706
    >>> r.inhibitor_underdosed
    True


pyrestoolbox.gas.gas_hvf_beta
=============================

.. code-block:: python

    gas_hvf_beta(k, method='FK', phi=0.0, metric=False) -> float

Returns the Forchheimer high-velocity-flow (inertial) coefficient β, used in non-Darcy flow calculations. Three correlations are available. All return β in 1/ft (or 1/m when ``metric=True``) for permeability ``k`` in md.

.. list-table:: β correlations
   :widths: 8 25 20
   :header-rows: 1

   * - method
     - Correlation
     - Source
   * - ``'FK'``
     - β = 2.172e10 · k⁻¹·²⁰¹
     - Log-log fit of Firoozabadi & Katz (1979) JPT Feb, pp.211-216 (SPE-6827) consolidated-rock β(k) chart. Default.
   * - ``'JONES'``
     - β = 6.15e10 · k⁻¹·⁵⁵
     - Jones, S.C. (1987) SPE-16949
   * - ``'TCK'``
     - β = 1.88e10 / (k¹·⁴⁷ · φ⁰·⁵³)
     - Tek, Coats, Katz (1962) JPT Jul, pp.799-806. Requires ``phi``.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - k
     - float
     - Permeability (md). To evaluate β at a damaged-zone permeability, pass ``k' = k * krg``.
   * - method
     - string
     - ``'FK'`` (default), ``'JONES'``, or ``'TCK'``.
   * - phi
     - float
     - Porosity fraction. Required for ``method='TCK'``.
   * - metric
     - bool
     - If True, returns β in 1/m; else 1/ft. Defaults to False.

Examples:

.. code-block:: python

    >>> gas.gas_hvf_beta(100.0)                           # FK, 100 md
    86071589.04028143
    >>> gas.gas_hvf_beta(100.0, method='JONES')
    48851186.4355433
    >>> gas.gas_hvf_beta(100.0, method='TCK', phi=0.25)
    45003847.531218514
    >>> gas.gas_hvf_beta(100.0, metric=True)              # 1/m
    282387103.1505296


pyrestoolbox.gas.gas_non_darcy_skin
===================================

.. code-block:: python

    gas_non_darcy_skin(qg, k, h_perf, rw, mug, sg, krg=1.0,
                       beta_method='FK', phi=0.0, metric=False) -> dict

Returns the rate-dependent (high-velocity-flow) pseudoskin for a gas well. Computes β via ``gas_hvf_beta``, the non-Darcy coefficient D following Jones (1987) SPE-16949 / Odeh, Moreland, Schueler (1975) JPT Dec pp.1501-1504, and the skin S\ :sub:`hvf` = D · q\ :sub:`g`.

Formula (field units):

    D [day/MSCF] = 2.222×10⁻¹⁵ · β · γ\ :sub:`g` · k / (μ\ :sub:`g` · h\ :sub:`perf` · r\ :sub:`w`)

where β is in 1/ft, k in md, μ\ :sub:`g` in cP, h\ :sub:`perf` and r\ :sub:`w` in ft.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - qg
     - float
     - Gas rate (MSCF/D | sm3/D).
   * - k
     - float
     - Absolute permeability (md).
   * - h_perf
     - float
     - Perforated interval thickness (ft | m).
   * - rw
     - float
     - Wellbore radius (ft | m).
   * - mug
     - float
     - Gas viscosity at reservoir conditions (cP). Obtain from ``gas_ug`` or ``GasPVT.viscosity()``.
   * - sg
     - float
     - Gas specific gravity relative to air.
   * - krg
     - float
     - Gas relative permeability at critical oil saturation. If < 1.0, β is evaluated at the damaged-zone permeability ``k' = k * krg``, producing a more pessimistic S\ :sub:`hvf`. Defaults to 1.0 (undamaged).
   * - beta_method
     - string
     - ``'FK'`` | ``'JONES'`` | ``'TCK'`` (see `pyrestoolbox.gas.gas_hvf_beta`_). Defaults to ``'FK'``.
   * - phi
     - float
     - Porosity fraction. Required for ``beta_method='TCK'``.
   * - metric
     - bool
     - If True, inputs in (sm3/D, m, m) and D returned in day/sm3; else (MSCF/D, ft, ft) and D in day/MSCF. Defaults to False.

Returns a ``dict`` with keys:

.. list-table:: Returns
   :widths: 10 40
   :header-rows: 1

   * - Key
     - Description
   * - ``'beta'``
     - Forchheimer coefficient (1/ft | 1/m)
   * - ``'D'``
     - Non-Darcy coefficient (day/MSCF | day/sm3). Suitable for passing directly to ``gas_rate_radial(D=...)``.
   * - ``'S_hvf'``
     - Rate-dependent skin (dimensionless).

Example (field units):

.. code-block:: python

    >>> r = gas.gas_non_darcy_skin(qg=10000, k=100, h_perf=100, rw=0.33,
    ...                             mug=0.025, sg=0.7)
    >>> round(r['beta'], 0)
    86071589.0
    >>> round(r['D'] * 1e5, 4)   # D in units of 1e-5 day/MSCF
    1.6227
    >>> round(r['S_hvf'], 4)
    0.1623

Using damaged-zone β (krg < 1) and the TCK correlation:

.. code-block:: python

    >>> r = gas.gas_non_darcy_skin(qg=10000, k=100, h_perf=100, rw=0.33,
    ...                             mug=0.025, sg=0.7, krg=0.7,
    ...                             beta_method='TCK', phi=0.22)
    >>> round(r['S_hvf'], 4)
    0.1534

Metric units (sm3/D, m) — same well, same skin:

.. code-block:: python

    >>> r = gas.gas_non_darcy_skin(qg=283168.5, k=100, h_perf=30.48, rw=0.1006,
    ...                             mug=0.025, sg=0.7, metric=True)
    >>> round(r['S_hvf'], 4)
    0.1622


pyrestoolbox.gas.gas_partial_penetration_skin
=============================================

.. code-block:: python

    gas_partial_penetration_skin(htot, htop, hbot, rw, kh_kv=10.0) -> float

Returns the partial-penetration pseudoskin S\ :sub:`p` for a vertical well in a radial anisotropic reservoir, using the analytical series solution of Streltsova-Adams, T.D. (1979) "Pressure Drawdown in a Well with Limited Flow Entry," SPE J. Nov 1979, pp.1469-1476 (SPE-7486).

The series is summed directly (up to 20,000 terms, vectorised) with Bessel K₀ evaluated via ``scipy.special.k0``. A warning is emitted when the tail of the summation suggests the series has not fully converged (typical for very thin wellbores, i.e. r\ :sub:`w`/h\ :sub:`tot` ≪ 10⁻³ combined with small k\ :sub:`v`).

**All length inputs must share a consistent unit** (ft or m — only ratios enter the formula), so there is no separate ``metric`` flag on this function.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - htot
     - float
     - Total formation thickness (no-flow to no-flow).
   * - htop
     - float
     - Distance from formation top to top of perforated interval.
   * - hbot
     - float
     - Distance from formation top to bottom of perforated interval.
   * - rw
     - float
     - Wellbore radius.
   * - kh_kv
     - float
     - Horizontal-to-vertical permeability anisotropy ratio (k\ :sub:`h`/k\ :sub:`v`). Defaults to 10.0. Set to 0 to treat vertical permeability as negligible (returns S\ :sub:`p` = 0).

Example (field units, ft):

.. code-block:: python

    >>> gas.gas_partial_penetration_skin(htot=164.04, htop=32.81, hbot=147.64,
    ...                                   rw=0.3543, kh_kv=8.0)
    2.155398050852645

Metric call using meters — same well, dimensionless result is essentially identical:

.. code-block:: python

    >>> gas.gas_partial_penetration_skin(htot=50.0, htop=10.0, hbot=45.0,
    ...                                   rw=0.108, kh_kv=8.0)
    2.155481748242132

Fully perforated interval (S\ :sub:`p` → 0):

.. code-block:: python

    >>> round(gas.gas_partial_penetration_skin(htot=100, htop=0, hbot=100, rw=0.33), 10)
    0.0


pyrestoolbox.gas.GasPVT
========================

.. code-block:: python

    GasPVT(sg=0.75, co2=0, h2s=0, n2=0, h2=0, zmethod='DAK', cmethod='PMC', metric=False)

Stores gas composition and method choices. Pre-computes critical temperature and pressure via ``gas_tc_pc()`` so they are not recalculated per call. Automatically selects BNS zmethod/cmethod when h2 > 0. Can be passed directly to ``fbhp()`` and ``operating_point()`` for VLP calculations, and to ``gas_rate_radial()`` and ``gas_rate_linear()`` for IPR rate calculations.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - sg
     - float
     - Gas SG relative to air. Defaults to 0.75
   * - co2
     - float
     - Molar fraction of CO2. Defaults to 0
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to 0
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to 0
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to 0. If positive, overrides zmethod and cmethod to 'BNS'
   * - zmethod
     - string or z_method
     - Method for Z-factor calculation. Defaults to 'DAK'
   * - cmethod
     - string or c_method
     - Method for critical property calculation. Defaults to 'PMC'
   * - metric
     - bool
     - If True, methods accept/return Eclipse METRIC units (barsa, deg C, kg/m3, rm3/sm3). Defaults to False

.. list-table:: Methods
   :widths: 15 40
   :header-rows: 1

   * - Method
     - Description
   * - ``z(p, degf)``
     - Returns Z-factor at pressure p (psia, or barsa if metric=True) and temperature degf (deg F, or deg C if metric=True)
   * - ``viscosity(p, degf)``
     - Returns gas viscosity (cP)
   * - ``density(p, degf)``
     - Returns gas density (lb/cuft, or kg/m3 if metric=True)
   * - ``bg(p, degf)``
     - Returns gas formation volume factor (rcf/scf, or rm3/sm3 if metric=True)

Examples:

.. code-block:: python

    >>> from pyrestoolbox import gas
    >>> gpvt = gas.GasPVT(sg=0.65, co2=0.1)
    >>> gpvt.z(2000, 180)
    0.9026719828498643
    >>> gpvt.viscosity(2000, 180)
    0.016666761192334678
    >>> gpvt.density(2000, 180)
    6.077743278791424
    >>> gpvt.bg(2000, 180)
    0.008164459661048012

Using metric units (pressure in barsa, temperature in deg C):

.. code-block:: python

    >>> gpvt_m = gas.GasPVT(sg=0.65, co2=0.1, metric=True)
    >>> gpvt_m.z(137.9, 82.2)
    0.9026433033953588
    >>> gpvt_m.density(137.9, 82.2)
    97.36871728783241

``GasPVT`` also exposes two skin helpers that reuse stored PVT state for viscosity:

.. code-block:: python

    >>> gpvt = gas.GasPVT(sg=0.70)
    >>> r = gpvt.non_darcy_skin(qg=10000, p=3000, degf=200,
    ...                          k=100, h_perf=100, rw=0.33, krg=0.7)
    >>> round(r['S_hvf'], 4)
    0.3044
    >>> round(gpvt.partial_penetration_skin(htot=50, htop=10, hbot=45,
    ...                                      rw=0.108, kh_kv=8), 4)
    2.1555

