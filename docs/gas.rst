===================================
Gas PVT
===================================

Function List
=============

+-------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| Gas PVT                 | - Gas Tc & Pc Calculation: `gas_tc_pc(...) <./docs/api.html#pyrestoolbox.gas_tc_pc>`_                                           |
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
===============

.. code-block:: python

    gas_tc_pc(sg, n2 = 0, co2 = 0, h2s = 0, cmethod = 'PMC', tc = 0, pc = 0) -> tuple

Returns a tuple of critical temperature (deg R) and critical pressure (psia) for hydrocarbon gas. If one or both of the tc and pc parameters are set to be non-zero, then this function will simply return that value for the corresponding critical parameter.

+---------------------------------------------------+---------------------------------------------------------------------------------------------+
| Parameter     |  Type                             | Description                                                                                 |
+---------------+-----------------------------------+---------------------------------------------------------------------------------------------+
| sg            | float                             | Gas SG relative to air                                                                      |
+---------------+-----------------------------------+---------------------------------------------------------------------------------------------+
| n2            | float                             | Molar fraction of Nitrogen. Defaults to zero if undefined                                   |
+---------------+-----------------------------------+---------------------------------------------------------------------------------------------+
| co2           | float                             | Molar fraction of CO2. Defaults to zero if undefined                                        |
+---------------+-----------------------------------+---------------------------------------------------------------------------------------------+
| h2s           | float                             | Molar fraction of H2S. Defaults to zero if undefined                                        |
+---------------+-----------------------------------+---------------------------------------------------------------------------------------------+
| cmethod       | String or c_method class object   | Method for calculating gas critical parameters                                              |
+---------------+-----------------------------------+---------------------------------------------------------------------------------------------+
| tc            | float                             | Critical gas temperature (deg R). Uses cmethod correlation if not specified                 |
+---------------+-----------------------------------+---------------------------------------------------------------------------------------------+
| n2            | float                             | Critical gas pressure (psia). Uses cmethod correlation if not specified                     |
+---------------+-----------------------------------+---------------------------------------------------------------------------------------------+

Examples:

.. code-block:: python

    >>> rtb.gas_tc_pc(sg=0.7, co2 = 0.15)
    (363.9387708314338, 738.3190067714969)
    
    >>> rtb.gas_tc_pc(sg=0.7, co2 = 0.15, tc=365, cmethod='SUT')
    (365, 709.2389730048743)

Returns a tuple of critical temperature (deg R) and critical pressure (psia) for hydrocarbon 