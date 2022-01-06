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
    
pyrestoolbox.oil_pbub
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
    

pyrestoolbox.oil_deno
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

    >>> rtb.oil_viso(p=2000, api=38, degf=165, pb=3500, rs=1000)
    0.416858469042502
    

pyrestoolbox.make_bot_og
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
     - Boolean flag that controls whether to export full table to excel, and separate PVDG and PVDO include files. Default is False.
   * - pvto
     - bool
     - Boolean flag that controls whether the pvto live oil Eclipse format will be generated. 
         - extends bubble point line up to maximum pressure
         - generates undersaturated oil propeties
         - writes out pvto include file
         
Examples:

.. code-block:: python

    >>> df, st_deno, st_deng, res_denw, res_cw, visw = rtb.make_bot(pi=4000, api=38, degf=175, sg_g=0.68, pmax=5000, pb=3900, rsb=2300, nrows=10)
    >>> df
        	Pressure (psia)	Rs (scf/stb)	Bo (rb/stb)	uo (cP)	Gas Z (v/v)	Bg (rb/mscf	ug (cP)	Bw (rb/stb)	uw (cP)
    0	14.7	0	1.062457	1.83067	0.998651	2.17E-04	0.012685	1.027687	0.356003
    1	726.8857	520.4289	1.289919	0.562439	0.938669	4.13E-06	0.013612	1.025392	0.357469
    2	1439.071	855.4766	1.440451	0.429764	0.895074	1.99E-06	0.015214	1.023133	0.358924
    3	2151.257	1199.924	1.598154	0.355871	0.875372	1.30E-06	0.017376	1.02091	0.360368
    4	2863.443	1595.424	1.782423	0.303031	0.881091	9.84E-07	0.019924	1.018722	0.361802
    5	3575.629	2068.841	2.005794	0.261505	0.908237	8.12E-07	0.02262	1.016568	0.363225
    6	3900	2316.211	2.122961	0.245254	0.925664	7.59E-07	0.02385	1.015597	0.36387
    7	4000	2300	2.104396	0.274803	0.931556	7.45E-07	0.024226	1.015299	0.364069
    8	4287.814	2300	2.076106	0.357028	0.949547	7.08E-07	0.025303	1.014446	0.364638
    9	5000	2300	2.01238	0.560491	0.997923	6.38E-07	0.027942	1.012356	0.366042

    

pyrestoolbox.sg_evolved_gas
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

    >>> rtb.sg_st_gas(114.7, rsp=1500, api=42, sg_sp=0.72, degf_sp=80)
    1.1923932340625523 
    

pyrestoolbox.sg_st_gas
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

    >>> rtb.sg_st_gas(114.7, rsp=1500, api=42, sg_sp=0.72, degf_sp=80)
    1.1923932340625523


        
pyrestoolbox.sgg_wt_avg
=======================

.. code-block:: python

    sgg_wt_avg (sg_sp, rsp, sg_st, rst) -> float

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

    >>> rtb.sgg_wt_avg (sg_sp=0.72, rsp=1000, sg_st=1.1, rst=5)
    0.7218905472636816
