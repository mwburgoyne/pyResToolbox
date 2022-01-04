===================================
``pyrestoolbox``
===================================

-----------------------------
A collection of Reservoir Engineering Utilities
-----------------------------

This set of functions focuses on those that the author uses often while creating programmatic solutions. These are the scripts that are often copy/pasted from previous work - sometimes slightly modified - resulting in a trail of slightly different versions over the years. Some attempt has been made here to make this implementation flexible enough such that it can be relied as-is going forward.

Includes functions to perform simple calculations including;
- Inflow for oil and gas
- PVT Calculations for oil
- PVT calculation for gas
- Creation of Black Oil Table information
- Creation of layered permeability distribution consistent with a Lorenze heterogeneity factor
- Extract problem cells information from Intesect (IX) print files
- Creation of Corey and LET relative permeability tables in Eclipse format

This is the initial public release, with improvements and additions expected over time. Apologies that it is only in oilfield units, with no current plans to add multi-unit support.

The current function list is as follows

+----------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| Inflow                     | `gas_rate_radial(...) </docs/api.html#pyrestoolbox.gas_rate_radial>`_,                                                          |
|                            | `gas_rate_linear(...) </docs/api.html#pyrestoolbox.gas_rate_linear>`_,                                                          |
|                            | `oil_rate_radial(...) </docs/api.html#pyrestoolbox.pyrestoolbox.oil_rate_radial>`_,                                             |
|                            | `oil_rate_linear(...) </docs/api.html#pyrestoolbox.pyrestoolbox.oil_rate_radial>`_,                                             |
+----------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| Gas PVT                    | `gas_tc_pc(...) </docs/api.html#pyrestoolbox.gas_tc_pc>`_,                                                                      |
|                            | `gas_z(...) </docs/api.html#pyrestoolbox.gas_z>`_,                                                                              |
|                            | `gas_ug(...) </docs/api.html#pyrestoolbox.gas_ug>`_,                                                                            |       
|                            | `gas_ugz(...) </docs/api.html#pyrestoolbox.gas_ugz>`_,                                                                          |         
|                            | `gas_cg(...) </docs/api.html#pyrestoolbox.gas_cg>`_,                                                                            |       
|                            | `gas_bg(...) </docs/api.html#pyrestoolbox.gas_bg>`_,                                                                            |       
|                            | `gas_den(...) </docs/api.html#pyrestoolbox.gas_den>`_,                                                                          |         
|                            | `gas_water_content(...) </docs/api.html#pyrestoolbox.gas_water_content>`_,                                                      |                             
|                            | `gas_ponz2p(...) </docs/api.html#pyrestoolbox.gas_ponz2p>`_,                                                                    |               
|                            | `gas_grad2sg(...) </docs/api.html#pyrestoolbox.gas_grad2sg>`_,                                                                  |                 
|                            | `gas_dmp(...) </docs/api.html#pyrestoolbox.gas_dmp>`_,                                                                          |
|                            | `gas_dsg2rv(...) </docs/api.html#pyrestoolbox.gas_dsg2rv>`_,                                                                    |
+----------------------------+---------------------------------------------------------------------------------------------------------------------------------+  
| Oil PVT                    | `oil_ja_sg(...) </docs/api.html#pyrestoolbox.oil_ja_sg>`_,                                                                      |
|                            | `oil_twu_props(...) </docs/api.html#pyrestoolbox.oil_twu_props>`_,                                                              |
|                            | `oil_ja_sg(...) </docs/api.html#pyrestoolbox.oil_ja_sg>`_,                                                                      |
|                            | `oil_rs_st(...) </docs/api.html#pyrestoolbox.oil_rs_st>`_,                                                                      |
|                            | `oil_pbub(...) </docs/api.html#pyrestoolbox.oil_pbub>`_,                                                                        |
|                            | `oil_rs_bub(...) </docs/api.html#oil_rs_bub>`_,                                                                                 |
|                            | `oil_rs(...) </docs/api.html#pyrestoolbox.oil_rs>`_,                                                                            |
|                            | `oil_co(...) </docs/api.html#pyrestoolbox.oil_co>`_,                                                                            |
|                            | `oil_deno(...) </docs/api.html#pyrestoolbox.oil_deno>`_,                                                                        |
|                            | `oil_bo(...) </docs/api.html#pyrestoolbox.oil_bo>`_,                                                                            |
|                            | `oil_viso(...) </docs/api.html#pyrestoolbox.oil_viso>`_,                                                                        |
|                            | `make_bot(...) </docs/api.html#pyrestoolbox.make_bot>`_,                                                                        |
|                            | `sg_evolved_gas(...) </docs/api.html#pyrestoolbox.sg_evolved_gas>`_,                                                            |
|                            | `sg_st_gas(...) </docs/api.html#pyrestoolbox.sg_st_gas>`_,                                                                      |
|                            | `sgg_wt_avg(...) </docs/api.html#pyrestoolbox.sgg_wt_avg>`_,                                                                    |
+----------------------------+---------------------------------------------------------------------------------------------------------------------------------+  
| Water PVT                  | `brine_props(...) </docs/api.html#pyrestoolbox.brine_props>`_,                                                                  |
+----------------------------+---------------------------------------------------------------------------------------------------------------------------------+  
| Permeability Layering      | `lorenzfromb(...) </docs/api.html#pyrestoolbox.lorenzfromb>`_,                                                                  |
|                            | `lorenz_from_flow_fraction(...) </docs/api.html#pyrestoolbox.lorenz_from_flow_fraction>`_,                                      |
|                            | `lorenz_2_flow_frac(...) </docs/api.html#pyrestoolbox.lorenz_2_flow_frac>`_,                                                    |
|                            | `lorenz_2_layers(...) </docs/api.html#pyrestoolbox.lorenz_2_layers>`_,                                                          |        
+----------------------------+---------------------------------------------------------------------------------------------------------------------------------+  
| Simulation Helpers         | `ix_extract_problem_cells(...) </docs/api.html#pyrestoolbox.ix_extract_problem_cells>`_                                         |
+----------------------------+---------------------------------------------------------------------------------------------------------------------------------+  
| Relative Permeability      | `rel_perm(...) </docs/api.html#pyrestoolbox.rel_perm>`_,                                                                           |
+----------------------------+---------------------------------------------------------------------------------------------------------------------------------+


Getting Started
===============

Install the library with `pip <https://pip.pypa.io/en/stable/>`_:

.. code-block:: shell

    pip install pyrestoolbox


Import library into your project and start using. 

A simple example below of estimating oil bubble point pressure.

.. code-block:: python

    >>> import restoolbox as rtb
    >>> rtb.oil_pbub(api=43, degf=185, rsb=2350, sg_g =0.72, pbmethod ='VALMC')
    5179.51086900132


Or creating black oil table information

.. code-block:: python

    >>> import matplotlib.pyplot as plt
    >>> df, st_deno, st_deng, res_denw, res_cw, visw = rtb.make_bot(pi=4000, api=38, degf=175, sg_g=0.68, pmax=5000, pb=3900, rsb=2300, nrows=50)
    >>> print('Stock Tank Oil Density:', st_deno, 'lb/cuft')
    >>> print('Stock Tank Gas Density:', st_deng, 'lb/cuft')
    >>> print('Reservoir Water Density:', res_denw, 'lb/cuft')
    >>> print('Reservoir Water Compressibility:', res_cw, '1/psi')
    >>> print('Reservoir Water Viscosity:', visw,'cP')

    >>> fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10,10))
    >>> ax1.plot(df['Pressure (psia)'], df['Rs (scf/stb)'])
    >>> ax2.plot(df['Pressure (psia)'], df['Bo (rb/stb)'])
    >>> ax3.plot(df['Pressure (psia)'], df['uo (cP)'])
    >>> ax4.semilogy(df['Pressure (psia)'], df['Co (1/psi)'])
    >>> ...
    >>> plt.show()
    Stock Tank Oil Density: 52.05522123893805 lb/cuft
	Stock Tank Gas Density: 0.052025361717109773 lb/cuft
	Reservoir Water Density: 61.40223160167964 lb/cuft
	Reservoir Water Compressibility: 2.930237693350768e-06 1/psi
	Reservoir Water Viscosity: 0.3640686136171888 cP

.. image:: https://github.com/vinomarkus/pyResToolbox/blob/main/docs/img/bot.png
    :alt: Black Oil Properties

.. code-block:: python
    >>> fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10,10))
    >>> ax1.semilogy(df['Pressure (psia)'], df['Bg (rb/mscf'])
    >>> ax2.plot(df['Pressure (psia)'], df['ug (cP)'])
    >>> ax3.plot(df['Pressure (psia)'], df['Gas Z (v/v)'])
    >>> ax4.semilogy(df['Pressure (psia)'], df['Cg (1/psi)'])
    >>> ...
    >>> plt.show()

.. image:: https://github.com/vinomarkus/pyResToolbox/blob/main/docs/img/dry_gas.png
    :alt: Dry Gas Properties
    
A set of Gas-Oil relative permeability curves with LET method

.. code-block:: python

    >>> df = rtb.rel_perm(rows=25, krtable='SGOF', krfamily='LET', kromax =1, krgmax =1, swc =0.2, sorg =0.15, Lo=2.5, Eo = 1.25, To = 1.75, Lg = 1.2, Eg = 1.5, Tg = 2.0)
    >>> plt.plot(df['Sg'], df['Krgo'], c = 'r', label='Gas')
    >>> plt.plot(df['Sg'], df['Krog'], c = 'g', label='Oil')
    >>> plt.title('SGOF Gas Oil LET Relative Permeability Curves')
    >>> plt.xlabel('Sg')
    >>> plt.ylabel('Kr')
    >>> plt.legend()
    >>> plt.grid('both')
    >>> plt.plot()

.. image:: https://github.com/vinomarkus/pyResToolbox/blob/main/docs/img/sgof.png
    :alt: SGOF Relative Permeability Curves

Or a set of Water-Oil curves with Corey method
.. code-block:: python
    >>> df = rtb.rel_perm(rows=25, krtable='SWOF', kromax =1, krwmax =0.25, swc =0.15, swcr = 0.2, sorw =0.15, no=2.5, nw=1.5)
    >>> plt.plot(df['Sw'], df['Krow'], c = 'g', label='Oil')
    >>> plt.plot(df['Sw'], df['Krwo'], c = 'b', label='Water')
    >>> plt.title('SWOF Water Oil Corey Relative Permeability Curves')
    >>> plt.xlabel('Sw')
    >>> plt.ylabel('Kr')
    >>> plt.legend()
    >>> plt.grid('both')
    >>> plt.plot()
    
.. image:: https://github.com/vinomarkus/pyResToolbox/blob/main/docs/img/swof.png
    :alt: SWOF Relative Permeability Curves

See the `API documentation </docs/api.html>`_ for a complete listing and usage examples.


Development
===========
``pyrestoolbox`` is maintained by Mark W. Burgoyne (`<https://github.com/vinomarkus>`_).
