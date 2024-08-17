===================================
``pyrestoolbox``
===================================

-----------------------------
A collection of Reservoir Engineering Utilities
-----------------------------

This set of functions focuses on those that the author uses often while crafting programming solutions. 
These are the scripts that are often copy/pasted from previous work - sometimes slightly modified - resulting in a trail of slightly different versions over the years. Some attempt has been made here to make this implementation flexible enough such that it can be relied on as-is going forward.

Note: Version 2.x now refactors functions into different modules, requiring seperate imports

Includes functions to perform simple calculations including;

- Inflow for oil and gas
- PVT Calculations for oil
- PVT calculation for gas, including up to 100% inerts for CO2, H2S, N2 and H2
- Return critical parameters for typical components
- Creation of Black Oil Table information
- Creation of layered permeability distribution consistent with a Lorenze heterogeneity factor
- Extract problem cells information from Intesect (IX) print files
- Generation of AQUTAB include file influence functions for use in ECLIPSE
- Creation of Corey and LET relative permeability tables in Eclipse format
- Calculation of Methane and CO2 saturated brine properties

`Changelist <https://github.com/mwburgoyne/pyResToolbox/blob/main/docs/changelist.rst>`_ 

Upgrade previous installations with

.. code-block:: shell

    pip install pyrestoolbox --upgrade


Module List
=============

+----------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| `gas <https://github.com/mwburgoyne/pyResToolbox/tree/main/pyrestoolbox/docs/gas.rst>`_            | - Gas Tc & Pc Calculation                                       |
|                                                                                                    | - Gas Z-Factor Calculation                                      |
|                                                                                                    | - Gas Viscosity                                                 |
|                                                                                                    | - Gas Viscosity * Z                                             |
|                                                                                                    | - Gas Compressibility                                           |
|                                                                                                    | - Gas Formation Volume Factor                                   |
|                                                                                                    | - Gas Density                                                   |
|                                                                                                    | - Gas Water of Condensation                                     |
|                                                                                                    | - Convert P/Z to P                                              |
|                                                                                                    | - Convert Gas Gradient to SG                                    |
|                                                                                                    | - Delta Pseudopressure                                          |
|                                                                                                    | - Gas Condensate FWS SG                                         |
|                                                                                                    | - Gas Flow Rate Radial                                          |
|                                                                                                    | - Gas Flow Rate Linear                                          |
+----------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| `oil <https://github.com/mwburgoyne/pyResToolbox/tree/main/pyrestoolbox/docs/oil.rst>`_            | - Oil Density from MW                                           |
|                                                                                                    | - Oil Critical Properties with Twu                              |
|                                                                                                    | - Incrememtal GOR post Separation                               |
|                                                                                                    | - Oil Bubble Point Pressure                                     |
|                                                                                                    | - Oil GOR at Pb                                                 |
|                                                                                                    | - Oil GOR at P                                                  |
|                                                                                                    | - Oil Compressibility                                           |
|                                                                                                    | - Oil Density                                                   |
|                                                                                                    | - Oil Formation Volume Factor                                   |
|                                                                                                    | - Oil Viscosity                                                 |
|                                                                                                    | - Generate Black Oil Table data                                 |
|                                                                                                    | - Estimate soln gas SG from oil                                 |
|                                                                                                    | - Estimate SG of gas post separator                             |
|                                                                                                    | - Calculate weighted average surface gas SG                     |
|                                                                                                    | - Oil API to SG                                                 |
|                                                                                                    | - Oil SG to API                                                 |
|                                                                                                    | - Oil Flow Rate Radial                                          |
|                                                                                                    | - Oil Flow Rate Linear                                          |
+----------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| `library <https://github.com/mwburgoyne/pyResToolbox/tree/main/pyrestoolbox/docs/library.rst>`_    | - Return critical parameters for typical single components      |
+----------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| `brine <https://github.com/mwburgoyne/pyResToolbox/tree/main/pyrestoolbox/docs/brine.rst>`_        | - Calculate suite of brine properties with variable methane     |
|                                                                                                    | - Calculate suite of CO2 saturated brine properties             |
+----------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| `layer <https://github.com/mwburgoyne/pyResToolbox/tree/main/pyrestoolbox/docs/layer.rst>`_        | - Lorenz coefficient from Beta value                            |
|                                                                                                    | - Lorenz coefficient from flow fraction                         |
|                                                                                                    | - Lorenz coefficient to flow fraction                           |
|                                                                                                    | - Lorenz coefficient to permeability array                      |
+----------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| `simtools <https://github.com/mwburgoyne/pyResToolbox/tree/main/pyrestoolbox/docs/simtools.rst>`_  | - Summarize IX convergence errors from PRT file                 |
|                                                                                                    | - Create Aquifer Influence Functions                            |
|                                                                                                    | - Perform recursive ECL or IX deck zip/check for INCLUDE files  |
|                                                                                                    | - Solve Rachford Rice for user specified feed Zis and Ki's      |
|                                                                                                    | - Create sets of rel perm tables                                |
+----------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+


Getting Started
===============

Install the library with  `pip <https://pip.pypa.io/en/stable/>`_:

.. code-block:: shell

    pip install pyrestoolbox


Import library into your project and start using. 

A simple example below of estimating oil bubble point pressure.

.. code-block:: python

    >>> from pyrestoolbox import oil
    >>> oil.oil_pbub(api=43, degf=185, rsb=2350, sg_g =0.72, pbmethod ='VALMC')
    5179.51086900132
    
A set of Gas-Oil relative permeability curves with the LET method

.. code-block:: python

    >>> import matplotlib.pyplot as plt
    >>> from pyrestoolbox import simtools 
    >>> df = simtools.rel_perm(rows=25, krtable='SGOF', krfamily='LET', kromax =1, krgmax =1, swc =0.2, sorg =0.15, Lo=2.5, Eo = 1.25, To = 1.75, Lg = 1.2, Eg = 1.5, Tg = 2.0)
    >>> plt.plot(df['Sg'], df['Krgo'], c = 'r', label='Gas')
    >>> plt.plot(df['Sg'], df['Krog'], c = 'g', label='Oil')
    >>> plt.title('SGOF Gas Oil LET Relative Permeability Curves')
    >>> plt.xlabel('Sg')
    >>> plt.ylabel('Kr')
    >>> plt.legend()
    >>> plt.grid('both')
    >>> plt.plot()

.. image:: https://github.com/mwburgoyne/pyResToolbox/blob/main/pyrestoolbox/docs/img/sgof.png
    :alt: SGOF Relative Permeability Curves

Or a set of Water-Oil relative permeability curves with the Corey method

.. code-block:: python

    >>> df = simtools.rel_perm(rows=25, krtable='SWOF', kromax =1, krwmax =0.25, swc =0.15, swcr = 0.2, sorw =0.15, no=2.5, nw=1.5)
    >>> plt.plot(df['Sw'], df['Krow'], c = 'g', label='Oil')
    >>> plt.plot(df['Sw'], df['Krwo'], c = 'b', label='Water')
    >>> plt.title('SWOF Water Oil Corey Relative Permeability Curves')
    >>> plt.xlabel('Sw')
    >>> plt.ylabel('Kr')
    >>> plt.legend()
    >>> plt.grid('both')
    >>> plt.plot()
    
.. image:: https://github.com/mwburgoyne/pyResToolbox/blob/main/pyrestoolbox/docs/img/swof.png
    :alt: SWOF Relative Permeability Curves

A set of dimensionless pressures for the constant terminal rate Van Everdingin & Hurst aquifer, along with an AQUTAB.INC export for use in ECLIPSE.

.. code-block:: python

    >>> ReDs = [1.5, 2, 3, 5, 10, 25, 1000]
    >>> tds, pds = simtools.influence_tables(ReDs=ReDs, export=True)
    >>> 
    >>> for p, pd in enumerate(pds):
    >>>     plt.plot(tds, pd, label = str(ReDs[p]))
    >>>     
    >>> plt.xscale('log')
    >>> plt.yscale('log')
    >>> plt.legend(loc='upper left')
    >>> plt.grid(which='both')
    >>> plt.xlabel('Dimensionless Time (tD)')
    >>> plt.ylabel('Dimensionless Pressure Drop (PD)')
    >>> plt.title('Constant Terminal Rate Solution')
    >>> plt.show()
    
.. image:: https://github.com/mwburgoyne/pyResToolbox/blob/main/pyrestoolbox/docs/img/influence.png
    :alt: Constant Terminal Rate influence tables

Or creating black oil table information for oil

.. code-block:: python

    >>> results = oil.make_bot_og(pi=4000, api=38, degf=175, sg_g=0.68, pmax=5000, pb=3900, rsb=2300, nrows=50)
    >>> df, st_deno, st_deng, res_denw, res_cw, visw, pb, rsb, rsb_frac, usat = results['bot'], results['deno'], results['deng'], results['denw'], results['cw'], results['uw'], results['pb'], results['rsb'], results['rsb_scale'], results['usat']
    >>> 
    >>> print('Stock Tank Oil Density:', st_deno, 'lb/cuft')
    >>> print('Stock Tank Gas Density:', st_deng, 'lb/cuft')
    >>> print('Reservoir Water Density:', res_denw, 'lb/cuft')
    >>> print('Reservoir Water Compressibility:', res_cw, '1/psi')
    >>> print('Reservoir Water Viscosity:', visw,'cP')
    >>> 
    >>> fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10,10))
    >>> ax1.plot(df['Pressure (psia)'], df['Rs (mscf/stb)'])
    >>> ax2.plot(df['Pressure (psia)'], df['Bo (rb/stb)'])
    >>> ax3.plot(df['Pressure (psia)'], df['uo (cP)'])
    >>> ax4.semilogy(df['Pressure (psia)'], df['Co (1/psi)'])
    >>> 
    >>> fig.suptitle('Black Oil Properties')
    >>> ax1.set_title("Rs vs P")
    >>> ax1.set_ylabel('Rs (mscf/stb)')
    >>> ax1.set_xlabel('Pressure (psia)')
    >>> ax1.grid('both')
    >>> 
    >>> ax2.set_title("Bo vs P")
    >>> ax2.set_ylabel('Bo (rb/stb)')
    >>> ax2.set_xlabel('Pressure (psia)')
    >>> ax2.grid('both')
    >>> 
    >>> ax3.set_title("Viso vs P")
    >>> ax3.set_xlabel('Pressure (psia)')
    >>> ax3.set_ylabel('Viscosity (cP)')
    >>> ax3.grid('both')
    >>> 
    >>> ax4.set_title("Co vs P")
    >>> ax4.set_ylabel('Co (1/psi)')
    >>> ax4.set_xlabel('Pressure (psia)')
    >>> ax4.grid('both')
    >>> 
    >>> plt.tight_layout()
    >>> plt.show()
    Iteratively solving for Rsb fraction to use in order to harmonize user specified Pb and Rsb
    
    Stock Tank Oil Density: 52.09203539823009 lb/cuft
    Stock Tank Gas Density: 0.052046870460837856 lb/cuft
    Reservoir Water Density: 61.40223160167964 lb/cuft
    Reservoir Water Compressibility: 2.930237693350768e-06 1/psi
    Reservoir Water Viscosity: 0.3640686136171888 cP

.. image:: https://github.com/mwburgoyne/pyResToolbox/blob/main/pyrestoolbox/docs/img/bot.png
    :alt: Black Oil Properties
    
And gas

.. code-block:: python

    >>> fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10,10))
    >>> ax1.semilogy(df['Pressure (psia)'], df['Bg (rb/mscf'])
    >>> ax2.plot(df['Pressure (psia)'], df['ug (cP)'])
    >>> ax3.plot(df['Pressure (psia)'], df['Gas Z (v/v)'])
    >>> ax4.semilogy(df['Pressure (psia)'], df['Cg (1/psi)'])
    >>> ...
    >>> plt.show()

.. image:: https://github.com/mwburgoyne/pyResToolbox/blob/main/pyrestoolbox/docs/img/dry_gas.png
    :alt: Dry Gas Properties
    
With ability to generate Live Oil PVTO style table data as well

.. code-block:: python

    >>> pb = 4500
    >>> results = oil.make_bot_og(pvto=True, pi=4000, api=38, degf=175, sg_g=0.68, pmax=5500, pb=pb, nrows=25, export=True)
    >>> df, st_deno, st_deng, res_denw, res_cw, visw, pb, rsb, rsb_frac, usat = results['bot'], results['deno'], results['deng'], results['denw'], results['cw'], results['uw'], results['pb'], results['rsb'], results['rsb_scale'], results['usat']
    >>> 
    >>> if len(usat) == 0:
    >>>     usat_flag = False
    >>> else:
    >>>     usat_flag=True
    >>>     usat_p, usat_bo, usat_uo = usat 
    >>> 
    >>> try:
    >>>     pb_idx = df['Pressure (psia)'].tolist().index(pb)
    >>>     bob = df['Bo (rb/stb)'].iloc[pb_idx]
    >>>     rsb = df['Rs (mscf/stb)'].iloc[pb_idx]
    >>>     uob = df['uo (cP)'].iloc[pb_idx]
    >>>     cob = df['Co (1/psi)'].iloc[pb_idx]
    >>>     no_pb = False
    >>> except:
    >>>     print('Pb was > Pmax')
    >>>     no_pb = True
    >>> 
    >>> print('Pb (psia):', pb)
    >>> print('Bob (rb/stb):', bob)
    >>> print('Rsb (mscf/stb):', rsb)
    >>> print('Rsb Scaling Required:', rsb_frac)
    >>> print('Visob (cP):', uob)
    >>> print('Cob (1/psi):', cob,'\n')
    >>> print('Stock Tank Oil Density:', st_deno, 'lb/cuft')
    >>> print('Stock Tank Gas Density:', st_deng, 'lb/cuft')
    >>> print('Reservoir Water Density:', res_denw, 'lb/cuft')
    >>> print('Reservoir Water Compressibility:', res_cw, '1/psi')
    >>> print('Reservoir Water Viscosity:', visw,'cP')
    >>> 
    >>> fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10,10))
    >>> ax1.plot(df['Pressure (psia)'], df['Rs (mscf/stb)'])
    >>> ax2.plot(df['Pressure (psia)'], df['Bo (rb/stb)'])
    >>> ax3.plot(df['Pressure (psia)'], df['uo (cP)'])
    >>> ax4.semilogy(df['Pressure (psia)'], df['Co (1/psi)'])
    >>> 
    >>> ax1.plot([pb], [rsb], 'o', c='r')
    >>> ax2.plot([pb], [bob], 'o', c='r')
    >>> ax3.plot([pb], [uob], 'o', c='r')
    >>> ax4.plot([pb], [cob], 'o', c='r')
    >>> 
    >>> if usat_flag:
    >>>     if no_pb == False:
    >>>         for i in range(len(usat_bo)):
    >>>             ax2.plot(usat_p[i], usat_bo[i], c='k')
    >>>             ax3.plot(usat_p[i], usat_uo[i], c='k')
    >>> 
    >>> fig.suptitle('Black Oil Properties')
    >>> ..
    >>> ..
    >>> plt.show()
    Pb (psia): 4500
    Bob (rb/stb): 1.6072798403441817
    Rsb (mscf/stb): 1.2863705330979234
    Rsb Scaling Required: 0.9713981737449556
    Visob (cP): 0.3422139569449832
    Cob (1/psi): 5.711273668114706e-05 
    
    Stock Tank Oil Density: 52.05522123893805 lb/cuft
    Stock Tank Gas Density: 0.052025361717109773 lb/cuft
    Reservoir Water Density: 61.40223160167964 lb/cuft
    Reservoir Water Compressibility: 2.930237693350768e-06 1/psi
    Reservoir Water Viscosity: 0.3640686136171888 cP
    
.. image:: https://github.com/mwburgoyne/pyResToolbox/blob/main/pyrestoolbox/docs/img/bot_PVTO.png
    :alt: Live Oil Properties


Development
===========
``pyrestoolbox`` is maintained by Mark W. Burgoyne (`<https://github.com/mwburgoyne>`_).
