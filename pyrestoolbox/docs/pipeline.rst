=============
Gas Pipelines
=============

Steady-state, isothermal gas flow in a single pipeline segment, following Menon (2005), Gas Pipeline Hydraulics, CRC Press, Chapter 2 (doi 10.1201/9781420038224). The General Flow equation is offered with Colebrook-White or AGA NB-13 transmission factors, alongside the empirical Weymouth, Panhandle A, Panhandle B and IGT equations, all with the elevation correction. Z and viscosity come from the gas module at the average pipe pressure.


pyrestoolbox.pipeline.gas_pipeline
==================================

.. code-block:: python

    gas_pipeline(d, length, degf, sg=0.65, p1=None, p2=None, qg=None,
                 method='COLEBROOK', efficiency=1.0, roughness=None,
                 drag_factor=0.96, dz=0.0, co2=0, h2s=0, n2=0, h2=0,
                 zmethod='DAK', cmethod='PMC', metric=False) -> dict

Solves for whichever one of upstream pressure ``p1``, downstream pressure ``p2`` and rate ``qg`` is not given. The rate is

.. code-block:: text

    General Flow:  Q = 38.77 E F (Tb/Pb) [(P1^2 - e^s P2^2) / (G Tf Le Z)]^0.5 D^2.5
    Le = L (e^s - 1) / s,   s = 0.0375 G (H2 - H1) / (Tf Z)

with Q in scf/d, D in inches, L in miles, Tf in deg R and F = 2 / sqrt(f_Darcy). ``'COLEBROOK'`` takes f from the Colebrook-White equation (64/Re below Re 2000). ``'AGA'`` takes F as the smaller of the rough-pipe value 4 log10(3.7 D / e) and the partially turbulent value 4 Df log10(Re / (1.4125 Ft)), with Ft the Von Karman smooth-pipe factor. The empirical equations replace F, the exponents and the constant:

.. list-table:: Empirical equations (Menon, Eqs. 2.30, 2.34, 2.36, 2.38)
   :widths: 15 40
   :header-rows: 1

   * - Method
     - Rate
   * - ``'WEYMOUTH'``
     - Q = 433.5 E (Tb/Pb) [(P1^2 - e^s P2^2) / (G Tf Le Z)]^0.5 D^2.667
   * - ``'PANHANDLE_A'``
     - Q = 435.87 E (Tb/Pb)^1.0788 [(P1^2 - e^s P2^2) / (G^0.8539 Tf Le Z)]^0.5394 D^2.6182
   * - ``'PANHANDLE_B'``
     - Q = 737 E (Tb/Pb)^1.02 [(P1^2 - e^s P2^2) / (G^0.961 Tf Le Z)]^0.51 D^2.53
   * - ``'IGT'``
     - Q = 136.9 E (Tb/Pb) [(P1^2 - e^s P2^2) / (G^0.8 Tf Le mu^0.2)]^0.555 D^2.667, mu in lb/ft-s

Z and viscosity are evaluated at the flowing temperature and the average pressure 2/3 (p1 + p2 - p1 p2 / (p1 + p2)). Standard conditions are the library's 14.696 psia and 60 deg F; Menon's examples use 14.7 psia, about 0.03% different in rate. The segment has one temperature and one average pressure, so a long line with a large pressure drop should be split into segments and the calls chained (outlet of one is the inlet of the next). The tests reproduce Menon's Examples 8, 9, 11 and 12 to the printed digits. His Example 7 prints 380.99 MMscf/d, but that arithmetic used Tb = 540 deg R; with 520 deg R the equation gives 366.88.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - d
     - float
     - Pipe inside diameter (in | mm)
   * - length
     - float
     - Segment length (ft | m)
   * - degf
     - float
     - Average flowing gas temperature (deg F | deg C)
   * - sg
     - float
     - Gas SG relative to air. Default 0.65
   * - p1, p2
     - float
     - Upstream and downstream pressure (psia | barsa). Give two of p1, p2 and qg
   * - qg
     - float
     - Gas rate (Mscf/d | sm3/d)
   * - method
     - str
     - 'COLEBROOK' (default), 'AGA', 'WEYMOUTH', 'PANHANDLE_A', 'PANHANDLE_B' or 'IGT'
   * - efficiency
     - float
     - Pipeline efficiency E, multiplies the rate. Default 1.0 (Panhandle and Weymouth are usually run at 0.88-0.95)
   * - roughness
     - float
     - Absolute roughness (in | mm), COLEBROOK and AGA. Default 0.0007 in (0.01778 mm), bare steel
   * - drag_factor
     - float
     - AGA drag factor Df. Default 0.96 (average bend index, bare steel); 0.90-0.99 per Menon Table 1
   * - dz
     - float
     - Outlet elevation minus inlet elevation (ft | m). Default 0
   * - co2, h2s, n2, h2
     - float
     - Mole fractions in the gas
   * - zmethod, cmethod
     - str
     - Z-factor and critical-property methods, as for ``gas.gas_z``
   * - metric
     - bool
     - If True, inputs and outputs in Eclipse METRIC units

.. list-table:: Returns (dict)
   :widths: 10 15 40
   :header-rows: 1

   * - Key
     - Type
     - Description
   * - 'p1', 'p2'
     - float
     - Upstream and downstream pressure (psia | barsa)
   * - 'qg'
     - float
     - Gas rate (Mscf/d | sm3/d)
   * - 'pavg', 'z'
     - float
     - Average pressure (psia | barsa) and Z at it
   * - 'reynolds'
     - float
     - Reynolds number
   * - 'transmission_factor', 'friction_factor'
     - float
     - F the rate implies in the General Flow equation (for the empirical methods, their equivalent F) and the Darcy f = 4 / F^2
   * - 's'
     - float
     - Elevation parameter

Examples (NPS 18 line, 17.5 in ID, 20 miles, 0.6 SG gas at 70 deg F):

.. code-block:: python

    >>> from pyrestoolbox import pipeline
    >>> r = pipeline.gas_pipeline(d=17.5, length=105600, degf=70, sg=0.6, p1=1014.7, qg=150000)
    >>> round(r['p2'], 2), round(r['z'], 4), round(r['transmission_factor'], 2)
    (960.07, 0.879, 19.44)
    >>> r = pipeline.gas_pipeline(d=17.5, length=105600, degf=70, sg=0.6, p1=1014.7, qg=150000,
    ...                           method='PANHANDLE_B', efficiency=0.95)
    >>> round(r['p2'], 2)
    970.33
    >>> round(pipeline.gas_pipeline(d=17.5, length=105600, degf=70, sg=0.6, p1=1014.7, p2=900)['qg'])
    214949
    >>> r = pipeline.gas_pipeline(d=17.5, length=105600, degf=70, sg=0.6, p2=900, qg=150000, dz=1500)
    >>> round(r['p1'], 2), round(r['s'], 5)
    (991.45, 0.07212)
