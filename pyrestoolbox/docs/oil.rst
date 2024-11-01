===================================
Oil PVT & Flow Rates
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
        + 'BUR': Fast, can handle 100% inerts and Hydrogen. Tuned 5 component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024)
   * - cmethod
     - c_method
     - Method for calculating gas critical properties. Defaults to 'PMC' 
       Options are:
        + 'SUT': Sutton with Wichert & Aziz non-hydrocarbon corrections
        + 'PMC': Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
        + 'BUR': Burgoyne method (2024). If h2 > 0, or the 'BUR' method is used for Z-Factor then 'BUR' will automatically be used
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
        + 'STAN': Standing Correlation (1947)
        + 'VALMC': Valko-McCain Correlation (2003) - Only for oil_rs_bub (Rs at Pb)
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

    >>> from pyrestoolbox import oil
    >>> oil.oil_pbub(api=43, degf=185, rsb=2350, sg_g = 0.72, pbmethod ='STAN')
    6406.067846808766
    
    >>> oil.oil_pbub(api=43, degf=185, rsb=2350, sg_g = 0.72, pbmethod = oil.pb_method.STAN)
    6406.067846808766


Function List
=============

.. list-table:: Gas Functions
   :widths: 15 40
   :header-rows: 1

   * - Task
     - Function
   * - Oil Density from MW 
     - `pyrestoolbox.oil.oil_ja_sg`_  
   * - Oil Critical Properties with Twu
     - `pyrestoolbox.oil.oil_twu_props`_
   * - Incrememtal GOR post Separation
     - `pyrestoolbox.oil.oil_rs_st`_
   * - Oil Bubble Point Pressure
     - `pyrestoolbox.oil.oil_pbub`_
   * - Oil GOR at Pb
     - `pyrestoolbox.oil.oil_rs_bub`_
   * - Oil GOR at P
     - `pyrestoolbox.oil.oil_rs`_  
   * - Oil Compressibility
     - `pyrestoolbox.oil.oil_co`_  
   * - Oil Density
     - `pyrestoolbox.oil.oil_deno`_
   * - Oil Formation Volume Factor
     - `pyrestoolbox.oil.oil_bo`_
   * - Oil Viscosity
     - `pyrestoolbox.oil.oil_viso`_
   * - Generate Black Oil Table data
     - `pyrestoolbox.oil.make_bot_og`_
   * - Estimate soln gas SG from oil
     - `pyrestoolbox.oil.sg_evolved_gas`_
   * - Estimate SG of gas post separator
     - `pyrestoolbox.oil.sg_st_gas`_
   * - Weighted average surface gas SG
     - `pyrestoolbox.oil.sgg_wt_avg`_  
   * - Oil API from SG
     - `pyrestoolbox.oil.oil_api`_  
   * - Oil SG from API
     - `pyrestoolbox.oil.oil_sg`_  
   * - Oil Flow Rate Radial
     - `pyrestoolbox.oil.oil_rate_radial`_
   * - Oil Flow Rate Linear
     - `pyrestoolbox.oil.oil_rate_linear`_

pyrestoolbox.oil.oil_ja_sg
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

    >>> oil.ja_sg(mw=150, ja=0.5)
    0.8583666666666667

pyrestoolbox.oil.oil_twu_props
==================

.. code-block:: python

    oil_twu_props(mw, ja = 0, sg = 0, damp = 1) -> tuple

Returns tuple of liquid critical parameters - sg, tb (R), tc (R), pc (psia), vc (ft3/lbmol) - using correlation method from Twu (1984). Modified with damping factor proposed by A. Zick between 0 (paraffin) and 1 (original Twu). 
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

    >>> oil.oil_twu_props(mw=225, ja = 0.5)
    (0.8954444444444445,
    1068.3961103813851,
    1422.4620493584146,
    264.23402773211745,
    13.498328588856445)

    
pyrestoolbox.oil.oil_rs_st
===================

.. code-block:: python

    oil_rs_st(psp, degf_sp, api) -> float

Estimates incremental gas evolved from separator liquid as it equilibrates to stock tank conditions (scf/stb). Correlation reproduced from Valko McCain 2003 paper Eq 3-2. 

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

    >>> oil.oil_rs_st(psp=114.7, degf_sp=80, api=38)
    4.176458005559282
    
pyrestoolbox.oil.oil_pbub
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

    >>> oil.oil_pbub(api=43, degf=185, rsb=2350, sg_g =0.72)
    5179.51086900132
    
    >>> oil.oil_pbub(api=43, degf=185, rsb=2350, sg_sp = 0.72, pbmethod ='STAN')
    6390.285894698239
    
    
pyrestoolbox.oil.oil_rs_bub
===================

.. code-block:: python

    oil_rs_bub(api, degf, pb, sg_g =0, sg_sp =0, pbmethod ='VALMC', rsmethod='VELAR') -> float

Returns solution GOR (scf/stb) at bubble point pressure. Uses the inverse of the Bubble point pressure correlations, with the same method families. Note: At low pressures, the VALMC method will fail (generally when rsb < 10 scf/stb). The VALMC method will revert to the Standing method in these cases. 
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

    >>> oil.oil_rs_bub(api=43, degf=185, pb=5179.5, sg_sp = 0.72)
    2370.1268120547115
    

pyrestoolbox.oil.oil_rs
===================

.. code-block:: python

    oil_rs(api, degf, sg_sp, p, pb =0, rsb =0, rsmethod='VELAR', pbmethod='VALMC') -> float

Returns solution gas oil ratio (scf/stb) calculated with different correlations. Either pb, rsb or both need to be specified. If one is missing, the other will be calculated from correlation

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

    >>> oil.oil_rs(api = 43, degf = 185, sg_sp =0 .72, p = 3000, pb = 5179.5, rsb = 2370)
    1017.9424240354475
    
    >>> oil.oil_rs(api=43, degf=185, sg_sp=0.72, p=3000, rsb =2370)
    1010.0669446829819
    
    >>> oil.oil_rs(api=43, degf=185, sg_sp=0.72, p=3000, pb =5180)
    1018.0168241109982
    
    >>> oil.oil_rs(api=43, degf=185, sg_sp=0.72, p=3000, pb =5180, rsmethod ='STAN')
    949.754509735243

pyrestoolbox.oil.oil_co
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

    >>> oil.oil_co(p = 4500, api = 47, degf = 180, sg_sp = 0.72, rsb = 2750)
    8.807199545797315e-05
    
    >>> oil.oil_co(p=2000, api=47, degf=180, sg_sp =0.72, rsb =2750, pb=4945)
    0.00023290195865185949
    

pyrestoolbox.oil.oil_deno
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

    >>> oil.oil_deno(p=2000, degf=165, rs=1000, rsb=2000, sg_g = 0.72, api =38)
    40.95314479616728 

pyrestoolbox.oil.oil_bo
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
     - The method of oil FVF calculation to be employed. `Calculation Methods and Class Objects`_.
   * - denomethod
     - string or deno_method
     - The method of live oil density  calculation to be employed. `Calculation Methods and Class Objects`_.

Examples:

.. code-block:: python

    >>> oil.oil_bo(p=2000, pb=3000, degf=165, rs=1000, sg_o=0.8, sg_g =0.68)
    1.5038177989551806   
    
    >>> oil.oil_bo(p=2000, pb=3000, degf=165, rs=1000, sg_o=0.8, sg_g =0.68, bomethod='STAN')
    1.5393786735904431
    
pyrestoolbox.oil.oil_viso
========================

.. code-block:: python

    oil_viso(p, api, degf, pb, rs) -> float

Returns Oil Viscosity with Beggs-Robinson (1975) correlation at saturated pressures and Petrosky-Farshad (1995) at undersaturated pressures

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float
     - Pressure at observation (psia)
   * - api
     - float 
     - Stock tank oil density (degrees API)
   * - degf
     - float
     - Oil Temperature (deg F)
   * - pb
     - float
     - Original bubble point pressure of the oil (psia)
   * - rs
     - float
     - Solution GOR at pressure of interest (scf/stb).

Examples:

.. code-block:: python

    >>> oil.oil_viso(p=2000, api=38, degf=165, pb=3500, rs=1000)
    0.416858469042502
    

pyrestoolbox.oil.make_bot_og
=====================

.. code-block:: python

    make_bot_og(pi, api, degf, sg_g, pmax, pb =0, rsb =0, pmin =14.7, nrows = 20, wt =0, ch4_sat =0, comethod='EXPLT', zmethod='DAK', rsmethod='VELAR', cmethod'PMC', denomethod='SWMH', bomethod='MCAIN', pbmethod='VALMC', export=False) -> tuple

Creates data required for Oil-Gas-Water black oil tables. Returns dictionary of results, with index:
 - bot: Pandas table of blackoil data (for PVTO == False), or Saturated properties to pmax (if PVTO == True)
 - deno: ST Oil Density (lb/cuft)
 - deng: ST Gas Density (lb/cuft)
 - denw: Water Density at Pi (lb/cuft), 
 - cw: Water Compressibility at Pi (1/psi)
 - uw: Water Viscosity at Pi (cP))
 - pb: Bubble point pressure either calculated (if only Rsb provided), or supplied by user
 - rsb: Solution GOR at Pb either calculated (if only Pb provided), or supplied by user
 - rsb_scale: The scaling factor that was needed to match user supplied Pb and Rsb
 - usat: a list of understaurated values (if PVTO == True) [usat_p, usat_bo, usat_uo]. This will be empty if PVTO == False

If user species Pb or Rsb only, the corresponding property will be calculated. If both Pb and Rsb are specified, then Pb calculations will be adjusted to honor both


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - pi
     - float
     - Initial Pressure (psia).
   * - api
     - float
     - Density of stock tank liquid (API)
   * - degf
     - float
     - Oil Temperature (deg F)
   * - sg_g
     - float
     - Weighted average specific gravity of surface gas, inclusive of gas evolved after separation (relative to air).   
   * - pmax
     - float
     - Maximum pressure to calculate properties to.  
   * - pb
     - float
     - Original bubble point pressure (psia). Calculated from rsb if not specified.
   * - rsb
     - float
     - Original solution GOR at original bubble point pressure (scf/stb). Calculated from pb if not specified.
   * - pmin
     - float
     - Minimum pressure to evaluate pressures down to. Default = 25 psia
   * - nrows
     - int
     - Number of rows of table data to return
   * - wt
     - float
     - Salt wt% in brine (0-100).
   * - ch4_sat
     - float
     - Degree of methane saturation in the brine (0 - 1)
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
   * - export
     - bool
     - Boolean flag that controls whether to export full table to excel, and export separate PVDG, PVDO (and PVTO if requested) include files. Default is False.
   * - pvto
     - bool
     - Boolean flag that controls whether the pvto live oil Eclipse format will be generated. 
         - extends saturated pressures up to maximum pressure
         - generates undersaturated oil propeties at each pressure step
         - writes out pvto include file if export == True
         
Examples:

.. code-block:: python

    >>> results = oil.make_bot_og(pvto=False, pi=4000, api=38, degf=175, sg_g=0.68, pmax=5500, pb=4500, nrows=10, export=True)
    >>> df, st_deno, st_deng, res_denw, res_cw, visw, pb, rsb, rsb_frac, usat = results['bot'], results['deno'], results['deng'], results['denw'], results['cw'], results['uw'], results['pb'], results['rsb'], results['rsb_scale'], results['usat']
    >>> df
.. image:: https://github.com/mwburgoyne/pyResToolbox/blob/main/docs/img/bot_img.png
    :alt: Black Oil Table DataFrame

pyrestoolbox.oil.sg_evolved_gas
==============================

.. code-block:: python

    sg_evolved_gas(p, degf, rsb, api, sg_sp) -> float

Returns estimated specific gravity of gas evolved from insitu-oil due to depressurization below Pb. Uses McCain & Hill Correlation (1995, SPE 30773) 

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float
     - Pressure (psia).
   * - degf
     - float
     - Oil Temperature (deg F).
   * - rsb
     - float
     - Original solution GOR at original bubble point pressure (scf/stb)
   * - api
     - float
     - Density of stock tank liquid (API). Will calculate from sg_sto if not specified
   * - sg_sp
     - float
     - Separator gas gravity (relative to air). 


Examples:

.. code-block:: python

    >>> oil.sg_st_gas(114.7, rsp=1500, api=42, sg_sp=0.72, degf_sp=80)
    1.1923932340625523 
    

pyrestoolbox.oil.sg_st_gas
=======================

.. code-block:: python

    sg_st_gas(psp, rsp, api, sg_sp, degf_sp) -> float

Estimates specific gravity of gas evolving from liquid exiting the separator. Returns sg_st (Stock Tank gas specific gravity relative to air). Correlation reproduced from Valko McCain 2003 paper Eq 4-2

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - psp
     - float
     - Separator pressure (psia). 
   * - rsp
     - float
     - Separator GOR (separator scf / stb). 
   * - api
     - float
     - Density of stock tank liquid (API)
   * - degf_sp
     - float
     - Separator temperature (deg F). 
     
Examples:

.. code-block:: python

    >>> oil.sg_st_gas(114.7, rsp=1500, api=42, sg_sp=0.72, degf_sp=80)
    1.1923932340625523


        
pyrestoolbox.oil.sgg_wt_avg
=======================

.. code-block:: python

    sgg_wt_avg(sg_sp, rsp, sg_st, rst) -> float

Calculates weighted average specific gravity of surface gas (sg_g) from separator and stock tank gas properties. Returns sg_g (Weighted average surface gas SG relative to air). From McCain Correlations book, Eq 3.4

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - sg_sp
     - float
     - Separator pressure (psia). 
   * - rsp
     - float
     - Separator GOR (separator scf / stb). 
   * - sg_st
     - float
     - Specific gravity of incremental gas evolved from separator liquid as it equilibrates to stock tank conditions (relative to air)  
   * - rst
     - float
     - Incremental gas evolved from separator liquid as it equilibrates to stock tank conditions (scf/stb). 
     
Examples:

.. code-block:: python

    >>> oil.sgg_wt_avg (sg_sp=0.72, rsp=1000, sg_st=1.1, rst=5)
    0.7218905472636816


pyrestoolbox.oil.oil_api
=======================

.. code-block:: python

    oil_api(sg_value) -> float

Returns oil API given specific gravity value of oil

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - sg_value
     - float
     - Specific gravity (relative to water)
     
Examples:

.. code-block:: python

    >>> oil.oil_api (sg_value=0.82)
    41.0609756097561


pyrestoolbox.oil.oil_sg
=======================

.. code-block:: python

    oil_sg(api_value) -> float

Returns oil specific gravity given API value of oil

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - api_value
     - float
     - Oil API density
     
Examples:

.. code-block:: python

    >>> oil.oil_sg(api_value=45)
    0.8016997167138811
    
pyrestoolbox.oil.oil_rate_radial
======================

.. code-block:: python

    oil_rate_radial(k, h, pr, pwf, r_w, r_ext, uo, bo, S = 0, vogel = False, pb = 0) -> float or np.array

Returns liquid rate (stb/day) for radial flow using Darcy pseudo steady state equation with optional Vogel correction.
Arrays can be used for any one of k, h, pr or pwf, returning corresponding 1-D array of rates. Using more than one input array – while not prohibited - will not return expected results 

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
     - Net height for flow (ft).
   * - pr
     - float, list or np.array
     - Reservoir pressure (psia)
   * - pwf
     - float, list or np.array
     - BHFP (psia).
   * - r_w
     - float
     - Wellbore Radius (ft).
   * - r_ext
     - float
     - External Reservoir Radius (ft).
   * - uo
     - float
     - Liquid viscosity (cP). 
   * - bo
     - float
     - Liquid formation volume factor (rb/stb)
   * - S
     - float
     - Skin. Defaults to zero if undefined
   * - vogel
     - bool
     - Boolean flag indicating whether to use vogel Pb correction. Defaults to False
   * - pb
     - float
     - Bubble point pressure. Used only when Vogel correction is invoked
     
Examples:

.. code-block:: python

    >>> oil.oil_rate_radial(k=20, h=20, pr=1500, pwf=250, r_w=0.3, r_ext=1500, uo=0.8, bo=1.4, vogel=True, pb=1800)
    213.8147848023242
    
    >>> oil.oil_rate_radial(k=20, h=20, pr=[1500, 2000], pwf=250, r_w=0.3, r_ext=1500, uo=0.8, bo=1.4, vogel=True, pb=1800)
    array([213.8147848 , 376.58731835])
    
pyrestoolbox.oil.oil_rate_linear
======================

.. code-block:: python

    oil_rate_linear(k, pr, pwf, area, length, uo, bo, vogel = False, pb = 0) -> float or np.array

Returns liquid rate (stb/day) for linear flow using Darcy steady state equation with optional Vogel correction.
Arrays can be used for any one of k, pr, pwf or area, returning corresponding 1-D array of rates. Using more than one input array – while not prohibited - will not return expected results 

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
     - Reservoir pressure (psia)
   * - pwf
     - float, list or np.array
     - BHFP (psia).
   * - area
     - float, list or np.array
     - Net cross-sectional area perpendicular to direction of flow (ft2)
   * - length
     - float
     - Linear distance of fluid flow (ft)
   * - bo
     - float
     - Liquid formation volume factor (rb/stb)
   * - vogel
     - bool
     - Boolean flag indicating whether to use vogel Pb correction. Defaults to False
   * - pb
     - float
     - Bubble point pressure. Used only when Vogel correction is invoked
     
Examples:

.. code-block:: python

    >>> oil.oil_rate_linear(k=0.1, area=15000, pr=3000, pwf=500, length=500, uo=0.4, bo=1.5)
    14.08521246363274
    
    >>> oil.oil_rate_linear(k=[0.1, 1, 5, 10], area=15000, pr=3000, pwf=500, length=500, uo=0.4, bo=1.5)
    array([  14.08521246,  140.85212464,  704.26062318, 1408.52124636])