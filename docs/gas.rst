===================================
Gas PVT
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
        + 'LIN': An explicit linearized form that is faster than DAK or HY `(2015) <https://link.springer.com/article/10.1007/s13202-015-0209-3>`_
        + 'DAK': Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
        + 'HY': Hall & Yarborough (1973)
   * - cmethod
     - c_method
     - Method for calculating gas critical properties. Defaults to 'PMC' 
       Options are:
        + 'SUT': Sutton with Wichert & Aziz non-hydrocarbon corrections
        + 'PMC': Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.


Users can specify which calculation method to use either by passing an option string, or a class object to any given function. The implementation of class objects should make it easier to program in an IDE that supports type hinting

Examples:

Calculating bubble point pressure with Standing correlation via option string, and then via class object

.. code-block:: python

    >>> rtb.oil_pbub(api=43, degf=185, rsb=2350, sg_g =0.72, pbmethod ='STAN')
    6406.067846808766
    
    >>> rtb.oil_pbub(api=43, degf=185, rsb=2350, sg_g =0.72, pbmethod = rtb.pb_method.STAN)
    6406.067846808766


Function List
=============

+-------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| Gas PVT                 | - Gas Tc & Pc Calculation: `pyrestoolbox.gas_tc_pc`_                                                                            |
|                         | - Gas Z-Factor Calculation: `gas_z(...) <./docs/api.html#pyrestoolbox.gas_z>`_                                                  |
|                         | - Gas Viscosity: `gas_ug(...) <./docs/api.html#pyrestoolbox.gas_ug>`_                                                           |
|                         | - Gas Viscosity * Z: `gas_ugz(...) <./docs/api.html#pyrestoolbox.gas_ugz>`_                                                     |
|                         | - Gas Compressibility: `gas_cg(...) <./docs/api.html#pyrestoolbox.gas_cg>`_                                                     |
|                         | - Gas Formation Volume Factor: `gas_bg(...) <./docs/api.html#pyrestoolbox.gas_bg>`_                                             |   
|                         | - Gas Density: `gas_den(...) <./docs/api.html#pyrestoolbox.gas_den>`_                                                           |
|                         | - Gas Water of Condensation: `gas_water_content(...) <./docs/api.html#pyrestoolbox.gas_water_content>`_                         |                       
|                         | - Convert P/Z to P: `gas_ponz2p(...) <./docs/api.html#pyrestoolbox.gas_ponz2p>`_                                                |
|                         | - Convert Gas Gradient to SG: `gas_grad2sg(...) <./docs/api.html#pyrestoolbox.gas_grad2sg>`_                                    |            
|                         | - Delta Pseudopressure: `gas_dmp(...) <./docs/api.html#pyrestoolbox.gas_dmp>`_                                                  |
|                         | - Gas Condensate FWS SG: `gas_fws_sg(...) <./docs/api.html#pyrestoolbox.gas_fws_sg>`_                                           |
+-------------------------+---------------------------------------------------------------------------------------------------------------------------------+


pyrestoolbox.gas_tc_pc
======================

.. code-block:: python

    gas_tc_pc(sg, n2 = 0, co2 = 0, h2s = 0, cmethod = 'PMC', tc = 0, pc = 0) -> tuple

Returns a tuple of critical temperature (deg R) and critical pressure (psia) for hydrocarbon gas. If one or both of the tc and pc parameters are set to be non-zero, then this function will return that unchanged value for the corresponding critical parameter.

.. list-table:: Method Variables & Class Objects
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - sg
     - float
     - Gas SG relative to air  
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined  
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - cmethod
     - String or c_method `class object <Calculation Methods and Class Objects>`_
     - Method for calculating gas critical parameters  
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified  

Examples:

.. code-block:: python

    >>> rtb.gas_tc_pc(sg=0.7, co2 = 0.15)
    (363.9387708314338, 738.3190067714969)
    
    >>> rtb.gas_tc_pc(sg=0.7, co2 = 0.15, tc=365, cmethod='SUT')
    (365, 709.2389730048743)

