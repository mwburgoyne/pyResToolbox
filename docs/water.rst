===================================
Brine PVT
===================================

Calculates Brine properties from modified Spivey Correlation per McCain Petroleum Reservoir Fluid Properties pg 160. Includes effect of user specified salt concentration and degree of methane satuartion.

Returns tuple of (Bw (rb/stb), Density (sg), viscosity (cP), Compressibility (1/psi), Rw GOR (scf/stb))

pyrestoolbox.brine_props
======================

.. code-block:: python

    brine_props(p, degf, wt, ch4_sat) -> tuple

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float
     - Pressure (psia)
   * - degf
     - float
     - Temperature (deg F)
   * - wt
     - float
     - Salt weight% in the brine (0-100)
   * - ch4_sat
     - float
     - Degree of methane saturation (0 - 1). 0 = No Methane, 1 = 100% Methane saturated

Examples:

.. code-block:: python

    >>> from pyrestoolbox import pyrestoolbox as rtb
    >>> bw, lden, visw, cw, rsw = rtb.brine_props(p=160, degf=135, wt=1.5, ch4_sat=1.0)
    >>> print('Bw:', bw)
    >>> print('Denw:', lden)
    >>> print('Visw:', visw)
    >>> print('Cw:', cw)
    >>> print('Rsw:', rsw)
    Bw: 1.0151710978322923
    Denw: 0.9950036658248123
    Visw: 0.4993957925685796
    Cw: 0.00015465633691558178
    Rsw: 1.2549011339427625
