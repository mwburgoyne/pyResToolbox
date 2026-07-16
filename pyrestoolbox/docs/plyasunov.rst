===================================
Plyasunov Dissolved Gas Volumes
===================================

Infinite-dilution apparent molar volumes of gases dissolved in water, after
Plyasunov (Fluid Phase Equilibria, 2019-2021, Parts I-IV), with pure water
properties from IAPWS-IF97 Region 1. This is a supporting module: the brine
module uses it to correct brine density for dissolved gas
(``SoreideWhitson``), and most users will consume it indirectly through
``brine``. It is documented here because the functions are stable, physically
meaningful on their own, and useful wherever a partial molar volume of a
dissolved gas is needed.

Supported gases (case-insensitive strings): ``H2``, ``N2``, ``CH4``, ``CO2``,
``C2H6``, ``C3H8``, ``NC4H10`` (aliases ``N-C4H10``, ``NC4``) and ``H2S``.
An unknown gas raises ``ValueError`` listing valid options.

Unlike the rest of pyResToolbox, this module works in SI units throughout:
temperature in Kelvin, pressure in MPa, volumes in cm3/mol. Validity follows
IAPWS-IF97 Region 1: 273.15 K to 623.15 K, water saturation pressure to
100 MPa.

Import pattern:

.. code-block:: python

    >>> from pyrestoolbox import plyasunov

Function List
=============

.. list-table:: Plyasunov Functions
   :widths: 15 40
   :header-rows: 1

   * - Task
     - Function
   * - Apparent molar volume at infinite dilution
     - `pyrestoolbox.plyasunov.V_phi`_
   * - Infinite-dilution partial molar volume (same quantity)
     - `pyrestoolbox.plyasunov.V2_inf`_
   * - Cross second virial coefficient gas-water
     - `pyrestoolbox.plyasunov.B12`_
   * - Dimensionless A12_inf group
     - `pyrestoolbox.plyasunov.A12_inf`_
   * - Gas molecular weight lookup
     - `pyrestoolbox.plyasunov.gas_mw`_

pyrestoolbox.plyasunov.V_phi
============================

.. code-block:: python

    V_phi(gas, T, P) -> float

Returns the apparent molar volume of the dissolved gas at infinite dilution
(cm3/mol). Alias for ``V2_inf`` - at infinite dilution the apparent and
partial molar volumes coincide. Results are LRU-cached (256 entries) on
``(gas, T, P)``.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - gas
     - str
     - Gas name, case-insensitive. One of H2, N2, CH4, CO2, C2H6, C3H8, NC4H10, H2S
   * - T
     - float
     - Temperature (K)
   * - P
     - float
     - Pressure (MPa)

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float
     - Apparent molar volume at infinite dilution (cm3/mol)

Examples:

.. code-block:: python

    >>> plyasunov.V_phi('CO2', 373.15, 20.0)
    37.60783904679587

    >>> plyasunov.V_phi('CH4', 350.0, 10.0)
    38.77012223993382

    >>> plyasunov.V_phi('H2', 298.15, 10.0)
    25.966076587141185

pyrestoolbox.plyasunov.V2_inf
=============================

.. code-block:: python

    V2_inf(gas, T, P) -> float

Returns the infinite-dilution partial molar volume (cm3/mol) via
``V2_inf = A12_inf * kappa_T * R * T``, with the isothermal compressibility
``kappa_T`` and density of pure water from IAPWS-IF97 Region 1. Inputs and
returns as for ``V_phi``.

.. code-block:: python

    >>> plyasunov.V2_inf('CO2', 373.15, 20.0)
    37.60783904679587

pyrestoolbox.plyasunov.B12
==========================

.. code-block:: python

    B12(gas, T) -> float

Returns the cross second virial coefficient for the gas-water pair
(cm3/mol). Sources per gas: H2/N2/CH4 from Plyasunov Part IV power series,
CO2 from Hellmann (2019), C2H6/C3H8/NC4H10 from square-well group
contributions, H2S from square-well parameters.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - gas
     - str
     - Gas name, case-insensitive. One of H2, N2, CH4, CO2, C2H6, C3H8, NC4H10, H2S
   * - T
     - float
     - Temperature (K)

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float
     - Cross second virial coefficient B12 (cm3/mol)

Examples:

.. code-block:: python

    >>> plyasunov.B12('CO2', 373.15)
    -98.40309129210485

    >>> plyasunov.B12('CH4', 298.15)
    -55.12819074948947

pyrestoolbox.plyasunov.A12_inf
==============================

.. code-block:: python

    A12_inf(gas, T, rho1_star) -> float

Returns the dimensionless group ``A12_inf = V2_inf / (kappa_T * R * T)``
evaluated from the density polynomial (Plyasunov Eq. 32). Mostly of interest
when the water density is prescribed independently of IAPWS-IF97.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - gas
     - str
     - Gas name, case-insensitive. One of H2, N2, CH4, CO2, C2H6, C3H8, NC4H10, H2S
   * - T
     - float
     - Temperature (K)
   * - rho1_star
     - float
     - Pure water density (kg/m3)

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float
     - A12_inf (dimensionless)

Examples:

.. code-block:: python

    >>> plyasunov.A12_inf('CO2', 373.15, 958.35)
    25.353028387016067

pyrestoolbox.plyasunov.gas_mw
=============================

.. code-block:: python

    gas_mw(gas) -> float

Returns the molecular weight of the gas (g/mol) from the internal ``MW_GAS``
dictionary. Case-insensitive; raises ``ValueError`` for an unknown gas.

Examples:

.. code-block:: python

    >>> plyasunov.gas_mw('CO2')
    44.0095

    >>> plyasunov.gas_mw('nc4h10')
    58.122
