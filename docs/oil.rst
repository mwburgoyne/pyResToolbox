===================================
Oil PVT
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
   * - pbmethod
     - pb_method
     - Method for calculating bubble point pressure of oil. Defaults to 'VELAR'. 
       Choices are:
        + 'STAN': Standing Correlation (1947)
        + 'VALMC': Valko-McCain Correlation (2003)
        + 'VELAR': Velarde, Blasingame 
   * - rsmethod
     - rs_method
     - Method for calculating solution gas-oil ratio. Defaults to 'VELAR'
       Options are:
        + 'VELAR': Velarde, Blasingame & McCain (1999)
        + 'STAN': Standing Correlation (1947), using form from 
        + 'VASBG': Vasquez & Beggs Correlation (1984)
   * - comethod
     - co_method
     - Method for calculating oil compressibility. Defaults to 'EXPLT'.  
       Options are:
        + 'EXPLT': Explicit calculation with numerical derivatives assessed and calculation via co = -1/bo*(dbodp - bg*drsdp/5.61458). 
   * - denomethod
     - deno_method
     - Method for calculating oil density. Defaults to 'SWMH'
       Options are:
        + 'SWMH': Standing, White, McCain-Hill (1995)
   * - bomethod
     - bo_method
     - Method for calculating oil formation volume factor. Defaults to 'MCAIN'
       Options are:
        + 'STAN': Standing Correlation
        + 'MCAIN': McCain approach, calculating from densities – Default

Users can specify which calculation method to use either by passing an option string, or a class object to any given function. The implementation of class objects should make it easier to program in an IDE that supports type hinting

Examples:

Calculating bubble point pressure with Standing correlation via option string, and then via class object

.. code-block:: python

    >>> rtb.oil_pbub(api=43, degf=185, rsb=2350, sg_g = 0.72, pbmethod ='STAN')
    6406.067846808766
    
    >>> rtb.oil_pbub(api=43, degf=185, rsb=2350, sg_g = 0.72, pbmethod = rtb.pb_method.STAN)
    6406.067846808766


Function List
=============

.. list-table:: Gas Functions
   :widths: 15 40
   :header-rows: 1

   * - Task
     - Function
   * - Oil Density from MW 
     - `pyrestoolbox.oil_ja_sg`_  
   * - Oil Critical Properties with Twu
     - `pyrestoolbox.oil_twu_props`_
   * - Incrememtal GOR post Separation
     - `pyrestoolbox.oil_rs_st`_
   * - Oil Bubble Point Pressure
     - `pyrestoolbox.oil_pbub`_
   * - Oil GOR at Pb
     - `pyrestoolbox.oil_rs_bub`_
   * - Oil GOR at P
     - `pyrestoolbox.oil_rs`_  
   * - Oil Compressibility
     - `pyrestoolbox.oil_co`_  
   * - Oil Density
     - `pyrestoolbox.oil_deno`_
   * - Oil Formation Volume Factor
     - `pyrestoolbox.oil_bo`_
   * - Oil Viscosity
     - `pyrestoolbox.oil_viso`_
   * - Generate Black Oil Table data
     - `pyrestoolbox.make_bot_og`_
   * - Estimate soln gas SG from oil
     - `pyrestoolbox.sg_evolved_gas`_
   * - Estimate SG of gas post separator
     - `pyrestoolbox.sg_st_gas`_
   * - Calculate weighted average surface gas SG
     - `pyrestoolbox.sgg_wt_avg`_  

pyrestoolbox.oil_ja_sg
======================

.. code-block:: python

    oil_ja_sg(mw, ja) -> float

Returns liquid hydrocarbon specific gravity (SG relative to water) at stock tank conditions using Jacoby Aromaticity Factor relationship 
.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - mw
     - float
     - Molecular weight of the liquid (g/gmole or lb/lb.mol)   
   * - ja
     - float
     - Jacoby aromaticity factor, vary between 0 (Paraffins) - 1 (Aromatic).

Examples:

.. code-block:: python

    >>> rtb.ja_sg(mw=150, ja=0.5)
    0.8583666666666667

pyrestoolbox.oil_twu_props
==================

.. code-block:: python

    oil_twu_props(mw, ja = 0, sg = 0, damp = 1) -> tuple

Returns tuple of liquid critical parameters - sg, tb (R), tc (R), pc (psia), vc (ft3/lbmol) - using correlation method from Twu (1984). Modified with damping factor proposed by A. Zick between 0 (paraffin) and 1 (original Twu)
If sg is left as default, the Jacoby relationship shall be used to estimate specific gravity


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - mw
     - float
     - Molecular weight of the liquid (g/gmole or lb/lb.mol)   
   * - ja
     - float
     - Jacoby aromaticity factor, vary between 0 (Paraffins) - 1 (Aromatic).
   * - sg
     - float
     - Specific gravity of the liquid (fraction relative to water density). Will use Jacoby method to estimate sg from mw if undefined
   * - damp
     - float
     - Damping factor proposed by A. Zick, modifying Eq 5.78 in the Whitson Monograph to permit some flexibility in the degree of parafinicity. Varies between 0 (paraffin) and 1 (original Twu). Defaults to 1
     
Examples:

.. code-block:: python

    >>> rtb.twu_liq_props(mw=225, ja = 0.5)
    (0.8954444444444445,
    1068.3961103813851,
    1422.4620493584146,
    264.23402773211745,
    13.498328588856445)

    
pyrestoolbox.oil_rs_st
===================

.. code-block:: python

    oil_rs_st(psp, degf_sp, api) -> float

Estimates incremental gas evolved from separator liquid as it equilibrates to stock tank conditions (scf/stb). Correlation reproduced from Valko McCain 2003 paper Eq 3-2
Rsb = Rsp + Rst (Solution GOR at bubble point = Separator GOR + Stock Tank GOR). 
In absence of separator properties, a simple linear relationship with Rsp could be used instead.
rs_st = 0.1618 * Separator GOR (Adapted from Eq 3-4 in Valko McCain 2003 paper)


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - psp
     - float 
     - Separator pressure (psia). 
   * - degf_sp
     - float
     - Separator temperature (deg F).
   * - api
     - float
     - Density of stock tank liquid (API)
     
Examples:

.. code-block:: python

    >>> rtb.oil_rs_st(psp=114.7, degf_sp=80, api=38)
    4.176458005559282
    
pyrestoolbox.pyrestoolbox.oil_pbub
====================

.. code-block:: python

    pbub(api, degf, rsb, sg_g =0, sg_sp =0, pbmethod ='VALMC') -> float

Returns bubble point pressure (psia) calculated with different correlations. 
At least one of sg_g and sg_sp must be supplied. This function will make simple assumption to estimate missing gas sg if only one is provided.


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - api
     - float 
     - Density of stock tank liquid (API)
   * - degf
     - float
     - Oil temperature (deg F).   
   * - sg_g
     - float
     - Weighted average specific gravity of surface gas, inclusive of gas evolved after separation (relative to air). 
   * - sg_sp
     - float
     - Specific gravity of separator gas (relative to air)
   * - pbmethod
     - str or pb_method class object
     - The method of Pb calculation to be employed. See Class Objects section for details.

Examples:

.. code-block:: python

    >>> rtb.oil_pbub(api=43, degf=185, rsb=2350, sg_g =0.72)
    5179.51086900132
    
    >>> rtb.oil_pbub(api=43, degf=185, rsb=2350, sg_sp = 0.72, pbmethod ='STAN')
    6390.285894698239
    
    
pyrestoolbox.oil_rs_bub
===================

.. code-block:: python

    oil_rs_bub(api, degf, pb, sg_g =0, sg_sp =0, pbmethod ='VALMC', rsmethod='VELAR') -> float

Returns solution GOR (scf/stb) at bubble point pressure. Uses the inverse of the Bubble point pressure correlations, with the same method families. Note: At low pressures, the VALMC method will fail (generally when rsb < 10 scf/stb). The VALMC method will revert to the VELAR method in these cases
At least one of sg_g and sg_sp must be supplied. This function will make simple assumption to estimate missing gas sg if only one is provided.


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - api
     - Density of stock tank liquid (API)
   * - degf
     - float
     - Oil Temperature (deg F)
   * - sg_g
     - float
     - Weighted average specific gravity of surface gas, inclusive of gas evolved after separation (relative to air).   
   * - sg_sp
     - float
     - Specific gravity of separator gas (relative to air).  
   * - pbmethod
     - string or pb_method
     - The method of Pb calculation to be employed. `Calculation Methods and Class Objects`_.
   * - rsmethod
     - string or rs_method
     - The method of Rs calculation to be employed. `Calculation Methods and Class Objects`_.

Examples:

.. code-block:: python

    >>> rtb.oil_rs_bub(api=43, degf=185, pb=5179.5, sg_sp = 0.72)
    2370.1268120547115
    

pyrestoolbox.oil_rs
===================

.. code-block:: python

    oil_rs(api, degf, sg_sp, p, pb =0, rsb =0, rsmethod='VELAR', pbmethod='VALMC') -> float

Returns solution gas oil ratio (scf/stb) calculated with different correlations. Either pb, rsb or both need to be specified. If one is missing, the other will be calculated from correlation

   * - Parameter
     - Type
     - Description
   * - api
     - Density of stock tank liquid (API)
   * - degf
     - float
     - Oil Temperature (deg F)
   * - sg_g
     - float
     - Weighted average specific gravity of surface gas, inclusive of gas evolved after separation (relative to air).   
   * - p
     - float
     - Pressure (psia). 
   * - pb
     - float
     - Original bubble point pressure (psia) 
   * - rsb
     - float
     - Original solution GOR at original bubble point pressure (scf/stb)
   * - rsmethod
     - string or rs_method
     - The method of Rs calculation to be employed. `Calculation Methods and Class Objects`_.
   * - pbmethod
     - string or pb_method
     - The method of Pb calculation to be employed. `Calculation Methods and Class Objects`_.

Examples:

.. code-block:: python

    >>> rtb.oil_rs(api = 43, degf = 185, sg_sp =0 .72, p = 3000, pb = 5179.5, rsb = 2370)
    1017.9424240354475
    
    >>> rtb.oil_rs(api=43, degf=185, sg_sp=0.72, p=3000, rsb =2370)
    1010.0669446829819
    
    >>> rtb.oil_rs(api=43, degf=185, sg_sp=0.72, p=3000, pb =5180)
    1018.0168241109982
    
    >>> rtb.oil_rs(api=43, degf=185, sg_sp=0.72, p=3000, pb =5180, rsmethod ='STAN')
    949.754509735243

pyrestoolbox.oil_co
=====================

.. code-block:: python

    oil_co(p, api,  degf, sg_sp =0, sg_g =0, pb =0, rsb =0, comethod='EXPLT', zmethod='DAK', rsmethod='VELAR', cmethod='PMC', denomethod='SWMH', bomethod='MCAIN', pbmethod='VALMC') -> float

Returns oil compressibility (1/psi) calculated with Co = -1/Bo *[dBodp - Bg*dRsdp], using correlation values and their numerical derivatives. 
At least one of sg_g and sg_sp must be supplied. This function will make simple assumption to estimate missing gas sg if only one is provided.
Either pb, rsb or both need to be specified. If one is missing, the other will be calculated from correlation 


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float
     - Pressure (psia).
   * - api
     - Density of stock tank liquid (API)
   * - degf
     - float
     - Oil Temperature (deg F)
   * - sg_sp
     - float
     - Separator gas gravity (relative to air). 
   * - sg_g
     - float
     - Weighted average specific gravity of surface gas, inclusive of gas evolved after separation (relative to air).   
   * - pb
     - float
     - Original bubble point pressure (psia) 
   * - rsb
     - float
     - Original solution GOR at original bubble point pressure (scf/stb)
   * - comethod
     - string or co_method
     - The method of Compressibility calculation to be employed. `Calculation Methods and Class Objects`_.
   * - zmethod
     - string or z_method
     - The method of gas z-factor calculation to be employed. `Calculation Methods and Class Objects`_.
   * - rsmethod
     - string or rs_method
     - The method of Rs calculation to be employed. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - The method of critical gas property calculation to be employed. `Calculation Methods and Class Objects`_.
   * - denomethod
     - string or deno_method
     - The method of live oil density  calculation to be employed. `Calculation Methods and Class Objects`_.
   * - bomethod
     - string or bo_method
     - The method of Bo calculation to be employed. `Calculation Methods and Class Objects`_.
   * - pbmethod
     - string or pb_method
     - The method of Rs calculation to be employed. `Calculation Methods and Class Objects`_.


Examples:

.. code-block:: python

    >>> rtb.oil_co(p = 4500, api = 47, degf = 180, sg_sp = 0.72, rsb = 2750)
    8.807199545797315e-05
    
    >>> rtb.oil_co(p=2000, api=47, degf=180, sg_sp =0.72, rsb =2750, pb=4945)
    0.00023290195865185949
    

pyrestoolbox.pyrestoolbox.oil_deno
==============================

.. code-block:: python

    oil_deno(p, degf, rs, rsb, sg_g = 0, sg_sp = 0, pb = 1e6, sg_o =0, api =0, denomethod='SWMH') -> float

Returns live oil density (lb/cuft). 
At least one of sg_g and sg_sp must be supplied. This function will make simple assumption to estimate missing gas sg if only one is provided.
At least one of sg_o and api must be supplied. One will be calculated from the other if only one supplied. If both specified, api will be used.
pb only needs to be set when pressures are above pb. For saturated oil, this can be left as default

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float
     - Pressure (psia).
   * - rs
     - float
     - Solution GOR at pressure of interest (scf/stb).
   * - rsb
     - float
     - Original solution GOR at original bubble point pressure (scf/stb)
   * - sg_g
     - float
     - Weighted average specific gravity of surface gas, inclusive of gas evolved after separation (relative to air).   
   * - sg_sp
     - float
     - Separator gas gravity (relative to air). 
   * - pb
     - float
     - Original bubble point pressure (psia) 
   * - sg_sto
     - float
     - Specific gravity of stock tank liquid (rel water). Will calculate from api if not specified
   * - api
     - Density of stock tank liquid (API). Will calculate from sg_sto if not specified
   * - denomethod
     - string or deno_method
     - The method of live oil density  calculation to be employed. `Calculation Methods and Class Objects`_.


Examples:

.. code-block:: python

    >>> rtb.oil_deno(p=2000, degf=165, rs=1000, rsb=2000, sg_g = 0.72, api =38)
    40.95314479616728 

pyrestoolbox.oil_bo
=======================

.. code-block:: python

    oil_bo(p, pb, degf, rs, rsb, sg_o, sg_g =0, sg_sp =0, bomethod='MCAIN', denomethod='SWMH') -> float

Returns oil formation volume factor calculated with different correlations
At least one of sg_g and sg_sp must be supplied. This function will make simple assumption to estimate missing gas sg if only one is provided.


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float
     - Pressure (psia).
   * - pb
     - float
     - Original bubble point pressure (psia) 
   * - degf
     - float
     - Oil temperature (deg F).
   * - rs
     - float
     - Solution GOR at pressure of interest (scf/stb).
   * - rsb
     - float
     - Original solution GOR at original bubble point pressure (scf/stb)
   * - sg_o
     - float
     - Specific gravity of live oil (rel water).
   * - sg_g
     - float
     - Weighted average specific gravity of surface gas, inclusive of gas evolved after separation (relative to air).   
   * - sg_sp
     - float
     - Separator gas gravity (relative to air). 
   * - bomethod
     - string or bo_method
     - The method of oil FVF calculation to be employed. See Class Objects section for details.
   * - denomethod
     - string or deno_method
     - The method of live oil density  calculation to be employed. `Calculation Methods and Class Objects`_.

Examples:

.. code-block:: python

    >>> rtb.oil_bo(p=2000, pb=3000, degf=165, rs=1000, sg_o=0.8, sg_g =0.68)
    1.5038177989551806   
    
    >>> rtb.oil_bo(p=2000, pb=3000, degf=165, rs=1000, sg_o=0.8, sg_g =0.68, bomethod='STAN')
    1.5393786735904431
    
pyrestoolbox.oil_viso
========================

.. code-block:: python

    gas_grad2sg( grad, p, degf, zmethod='DAK', cmethod='PMC', n2 = 0, co2 = 0, h2s = 0, tc = 0, pc = 0, rtol = 1E-7) -> float

Returns gas specific gravity consistent with observed gas gradient. Calculated through iterative solution method. Will fail if gas SG is below 0.55, or greater than 1.75

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - grad
     - float
     - Observed gas gradient (psi/ft)
   * - p
     - float, list or np.array 
     - Pressure at observation (psia)
   * - degf
     - float
     - Reservoir Temperature (deg F)
   * - zmethod
     - string or z_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined  
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified  

Examples:

.. code-block:: python

    >>> rtb.gas_grad2sg(grad=0.0657, p=2500, degf=175)
    0.7500786632299423   
    

pyrestoolbox.gas_dmp
=====================

.. code-block:: python

    gas_dmp(p1, p2, degf, sg, zmethod='DAK', cmethod = 'PMC', n2 = 0, co2 = 0, h2s = 0, tc = 0, pc = 0) -> float

Returns gas pseudo-pressure integral between two pressure points. Will return a positive value if p1 < p2, and a negative value if p1 > p2. 
Integrates the equation: m(p) = 2 * p / (ug * z) 

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p1
     - float, list or np.array 
     - First gas pressure (psia)
   * - p2
     - float, list or np.array 
     - Second gas pressure (psia)
   * - sg
     - float
     - Gas SG relative to air.
   * - degf
     - float
     - Reservoir Temperature (deg F)
   * - zmethod
     - string or z_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined  
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified  

Examples:

.. code-block:: python

    >>> rtb.gas_dmp(p1=1000, p2=2000, degf=185, sg=0.78, zmethod='HY', cmethod = 'SUT', n2 = 0.05, co2 = 0.1, h2s = 0.02)
    3690873383.43637  
    
    >>> rtb.gas_dmp(p1=2000, p2=1000, degf=185, sg=0.78, tc = 371, pc = 682)
    -3691052075.812854
        
pyrestoolbox.gas_fws_sg
=======================

.. code-block:: python

    gas_fws_sg(sg_g, cgr, api_st) -> float

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
     - Condensate gas ratio (stb/mmscf). 
   * - api_st
     - float
     - Density of stock tank liquid (API)  

Examples:

.. code-block:: python

    >>> rtb.gas_fws_sg(sg_g=0.855, cgr=30, api_st=53)
    0.9371015922844881
