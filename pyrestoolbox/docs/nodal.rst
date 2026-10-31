====================
Nodal Analysis & VLP
====================

Vertical Lift Performance (VLP) calculations, Inflow Performance Relationships (IPR), and nodal operating point analysis for gas and oil wells. The typical workflow is: define a wellbore geometry using ``Completion`` (with optional ``WellSegment`` objects for deviated wells), define reservoir properties using ``Reservoir``, then find the operating point where VLP and IPR curves intersect using ``operating_point()``. Individual VLP and IPR curves can also be generated independently.


Getting Started
===============

A minimal example finding the operating point for a gas well:

.. code-block:: python

    >>> from pyrestoolbox import nodal
    >>> c = nodal.Completion(tid=2.441, length=10000, tht=100, bht=200)
    >>> r = nodal.Reservoir(pr=3000, degf=200, k=10, h=50, re=1500, rw=0.35, S=2, D=0.001)
    >>> result = nodal.operating_point(thp=500, completion=c, reservoir=r, vlpmethod='HB', well_type='gas', gsg=0.65)
    >>> round(result['rate'], 1)
    10578.1
    >>> round(result['bhp'], 1)
    1583.4


Function List
=============

.. list-table:: Nodal Functions
   :widths: 15 40
   :header-rows: 1

   * - Task
     - Function
   * - Wellbore Segment
     - `pyrestoolbox.nodal.WellSegment`_
   * - Wellbore Completion
     - `pyrestoolbox.nodal.Completion`_
   * - Reservoir Description
     - `pyrestoolbox.nodal.Reservoir`_
   * - Flowing BHP
     - `pyrestoolbox.nodal.fbhp`_
   * - Tubing Head Pressure (inverse)
     - `pyrestoolbox.nodal.fthp`_
   * - VLP Outflow Curve
     - `pyrestoolbox.nodal.outflow_curve`_
   * - IPR Inflow Curve
     - `pyrestoolbox.nodal.ipr_curve`_
   * - Nodal Operating Point
     - `pyrestoolbox.nodal.operating_point`_
   * - Liquid-Loading Critical Gas Rate
     - `pyrestoolbox.nodal.critical_gas_rate`_
   * - Erosional Velocity (API RP 14E)
     - `pyrestoolbox.nodal.erosional_velocity`_
   * - Wellhead Choke (Gilbert family)
     - `pyrestoolbox.nodal.choke_flow`_
   * - Wellhead Choke (Sachdeva, critical and subcritical)
     - `pyrestoolbox.nodal.choke_flow_sachdeva`_


pyrestoolbox.nodal.WellSegment
==============================

.. code-block:: python

    WellSegment(md, id, deviation=0, roughness=None, metric=False)

Single wellbore segment with uniform geometry and deviation. Used to build multi-segment completions for deviated and horizontal wells.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - md
     - float
     - Measured depth of this segment (ft, or m if metric=True)
   * - id
     - float
     - Internal diameter (inches, or mm if metric=True)
   * - deviation
     - float
     - Deviation from vertical (degrees). 0 = vertical, 90 = horizontal. Defaults to 0
   * - roughness
     - float
     - Pipe roughness (inches, or mm if metric=True). Defaults to 0.0006 in / 0.01524 mm
   * - metric
     - bool
     - If True, interpret inputs in Eclipse METRIC units. Defaults to False

.. list-table:: Properties
   :widths: 15 40
   :header-rows: 1

   * - Property
     - Description
   * - ``tvd``
     - True vertical depth contribution of this segment (ft). Equals ``md * cos(deviation)``
   * - ``theta``
     - Angle from horizontal (radians) for multiphase correlations. Equals ``pi/2 - radians(deviation)``

Examples:

.. code-block:: python

    >>> from pyrestoolbox import nodal
    >>> seg = nodal.WellSegment(md=10000, id=2.441, deviation=0)
    >>> round(seg.tvd, 1)
    10000.0
    >>> seg = nodal.WellSegment(md=10000, id=2.441, deviation=45)
    >>> round(seg.tvd, 1)
    7071.1


pyrestoolbox.nodal.Completion
=============================

.. code-block:: python

    Completion(tid, length, tht, bht, rough=None, cid=0, crough=None, mpd=0, metric=False)
    Completion(segments=[...], tht=..., bht=..., metric=False)

Wellbore completion description for VLP calculations. Can be constructed in two ways:

**Legacy mode** (positional arguments): Defines a simple vertical wellbore with optional casing section below the tubing shoe (when ``mpd > length`` and ``cid > 0``).

**Segment mode** (``segments`` keyword): Accepts a list of ``WellSegment`` objects for arbitrary multi-segment wellbore definitions with deviation support. Temperature is interpolated linearly over true vertical depth (the geothermal gradient acts on TVD, so a deviated segment spans a smaller temperature interval than a vertical one of the same MD). In segment mode, ``metric`` only converts ``tht`` and ``bht``; segment dimensions are handled by each ``WellSegment``'s own ``metric`` flag.

.. list-table:: Legacy Mode Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - tid
     - float
     - Tubing ID (inches, or mm if metric=True)
   * - length
     - float
     - Tubing length from wellhead to tubing shoe (ft, or m if metric=True)
   * - tht
     - float
     - Tubing head (wellhead) temperature (deg F, or deg C if metric=True)
   * - bht
     - float
     - Bottom hole temperature (deg F, or deg C if metric=True)
   * - rough
     - float
     - Tubing roughness (inches, or mm if metric=True). Defaults to 0.0006 in / 0.01524 mm
   * - cid
     - float
     - Casing ID below tubing shoe (inches, or mm if metric=True). Defaults to 0 (no casing section)
   * - crough
     - float
     - Casing roughness (inches, or mm if metric=True). Defaults to 0.0006 in / 0.01524 mm
   * - mpd
     - float
     - Mid-perforation depth (ft, or m if metric=True). Defaults to length (no casing section)
   * - metric
     - bool
     - If True, interpret inputs in Eclipse METRIC units. Defaults to False

.. list-table:: Segment Mode Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - segments
     - list
     - List of WellSegment objects defining the wellbore from surface to bottom
   * - tht
     - float
     - Tubing head (wellhead) temperature (deg F, or deg C if metric=True)
   * - bht
     - float
     - Bottom hole temperature (deg F, or deg C if metric=True)
   * - metric
     - bool
     - If True, convert tht/bht from deg C. Defaults to False

.. list-table:: Properties
   :widths: 15 40
   :header-rows: 1

   * - Property
     - Description
   * - ``segments``
     - List of WellSegment objects defining the wellbore
   * - ``total_md``
     - Total measured depth across all segments (ft)
   * - ``total_tvd``
     - Total true vertical depth across all segments (ft)
   * - ``has_casing_section``
     - True if cid > 0 and mpd > length (legacy mode)
   * - ``casing_length``
     - Length of casing section (ft). Returns 0 if no casing section
   * - ``tubing_end_temperature``
     - Temperature at tubing shoe (deg F). Interpolated linearly between tht and bht

Examples:

.. code-block:: python

    >>> from pyrestoolbox import nodal
    >>> c = nodal.Completion(tid=2.441, length=10000, tht=100, bht=200)
    >>> c.has_casing_section
    False

Two-section completion with casing below tubing shoe:

.. code-block:: python

    >>> c2 = nodal.Completion(tid=2.441, length=9000, tht=100, bht=200, cid=6.184, mpd=10000)
    >>> c2.has_casing_section
    True
    >>> c2.casing_length
    1000

Multi-segment deviated well using WellSegment:

.. code-block:: python

    >>> segs = [nodal.WellSegment(md=5000, id=2.441, deviation=0), nodal.WellSegment(md=3000, id=2.441, deviation=45), nodal.WellSegment(md=2000, id=2.441, deviation=60)]
    >>> c3 = nodal.Completion(segments=segs, tht=100, bht=250)
    >>> c3.total_md
    10000
    >>> round(c3.total_tvd, 1)
    8121.3


Completion.geometry_at_md
-------------------------

.. code-block:: python

    completion.geometry_at_md(md)

Returns wellbore geometry at a given measured depth as a dictionary with keys ``'md'``, ``'tvd'``, ``'id'``, ``'deviation'``, and ``'roughness'``. All values are returned in the same unit system used to construct the Completion (oilfield or metric). Raises ``ValueError`` if ``md`` is outside the range ``[0, total_md]``.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - md
     - float
     - Measured depth from surface (ft, or m if metric=True)

.. list-table:: Returns (dict)
   :widths: 10 40
   :header-rows: 1

   * - Key
     - Description
   * - ``md``
     - Queried measured depth (ft | m)
   * - ``tvd``
     - True vertical depth at this MD (ft | m)
   * - ``id``
     - Internal diameter at this MD (inches | mm)
   * - ``deviation``
     - Deviation from vertical at this MD (degrees)
   * - ``roughness``
     - Pipe roughness at this MD (inches | mm)

Examples:

Query geometry within a multi-segment deviated well:

.. code-block:: python

    >>> segs = [nodal.WellSegment(md=5000, id=2.441, deviation=0), nodal.WellSegment(md=3000, id=2.441, deviation=45), nodal.WellSegment(md=2000, id=4.0, deviation=60)]
    >>> c = nodal.Completion(segments=segs, tht=100, bht=250)
    >>> g = c.geometry_at_md(6500)
    >>> round(g['tvd'], 1)
    6060.7
    >>> g['id']
    2.441
    >>> g['deviation']
    45

Query geometry in a legacy completion with casing section:

.. code-block:: python

    >>> c2 = nodal.Completion(tid=2.441, length=9000, tht=100, bht=200, cid=6.184, mpd=10000)
    >>> g2 = c2.geometry_at_md(9500)
    >>> g2['id']
    6.184
    >>> g2['tvd']
    9500.0


Completion.profile
------------------

.. code-block:: python

    completion.profile()

Returns a pandas DataFrame showing the wellbore profile at all segment boundaries. Crossover points (where ID or deviation changes) appear as two rows at the same MD/TVD with different properties, showing the transition from one segment to the next. Units match the Completion construction unit system.

.. list-table:: Output Columns
   :widths: 10 40
   :header-rows: 1

   * - Column
     - Description
   * - ``MD``
     - Measured depth (ft | m)
   * - ``TVD``
     - True vertical depth (ft | m)
   * - ``Deviation``
     - Deviation from vertical (degrees)
   * - ``ID``
     - Internal diameter (inches | mm)
   * - ``Roughness``
     - Pipe roughness (inches | mm)

Examples:

Profile of a multi-segment deviated well:

.. code-block:: python

    >>> segs = [nodal.WellSegment(md=5000, id=2.441, deviation=0), nodal.WellSegment(md=3000, id=2.441, deviation=45), nodal.WellSegment(md=2000, id=4.0, deviation=60)]
    >>> c = nodal.Completion(segments=segs, tht=100, bht=250)
    >>> df = c.profile()
    >>> len(df)
    6
    >>> list(df.columns)
    ['MD', 'TVD', 'Deviation', 'ID', 'Roughness']

Profile of a legacy completion with casing section:

.. code-block:: python

    >>> c2 = nodal.Completion(tid=2.441, length=9000, tht=100, bht=200, cid=6.184, mpd=10000)
    >>> df2 = c2.profile()
    >>> len(df2)
    4
    >>> df2['ID'].iloc[1]
    2.441
    >>> df2['ID'].iloc[2]
    6.184


pyrestoolbox.nodal.Reservoir
============================

.. code-block:: python

    Reservoir(pr, degf, k, h, re, rw, S=0, D=0, metric=False)

Reservoir description for IPR calculations.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - pr
     - float
     - Reservoir pressure (psia, or barsa if metric=True)
   * - degf
     - float
     - Reservoir temperature (deg F, or deg C if metric=True)
   * - k
     - float
     - Permeability (mD)
   * - h
     - float
     - Net pay thickness (ft, or m if metric=True)
   * - re
     - float
     - Drainage radius (ft, or m if metric=True)
   * - rw
     - float
     - Wellbore radius (ft, or m if metric=True)
   * - S
     - float
     - Skin factor. Defaults to 0
   * - D
     - float
     - Non-Darcy coefficient (day/mscf for gas, or day/sm3 if metric=True). Defaults to 0
   * - metric
     - bool
     - If True, interpret inputs in Eclipse METRIC units. Defaults to False

Examples:

.. code-block:: python

    >>> r = nodal.Reservoir(pr=3000, degf=200, k=10, h=50, re=1500, rw=0.35, S=2, D=0.001)
    >>> r.pr
    3000


pyrestoolbox.nodal.fbhp
=======================

.. code-block:: python

    fbhp(thp, completion, vlpmethod='WG', well_type='gas', gas_pvt=None, oil_pvt=None, qg_mmscfd=None, cgr=0, qw_bwpd=0, oil_vis=1.0, api=45, pr=0, qt_stbpd=0, gor=0, wc=0, wsg=1.07, injection=False, gsg=0.65, pb=0, rsb=0, sgsp=0.65, metric=False, return_profile=False, qg_mscfd=0) -> float

Returns flowing bottom hole pressure (psia, or barsa if metric=True) using the specified VLP correlation. Supports both gas and oil wells. The Completion object can define vertical wells (legacy mode) or multi-segment deviated wells using WellSegment objects. The ``sin(theta)`` multiplier on hydrostatic gradient enables proper deviated and horizontal well support.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - thp
     - float
     - Tubing head pressure (psia, or barsa if metric=True)
   * - completion
     - Completion
     - Completion object describing the wellbore
   * - vlpmethod
     - string or vlp_method
     - VLP method. `Calculation Methods and Class Objects`_
   * - well_type
     - string
     - 'gas' or 'oil'
   * - gas_pvt
     - GasPVT
     - Gas PVT object. If provided, the gas-well VLP march evaluates Z and gas viscosity with its methods and composition (as ``gas_pvt.z()`` and ``gas_pvt.viscosity()``: DAK with PMC by default, or HY or BNS as configured, impurities and H2/He included), and IPR calculations use the same. Without it the march uses Hall-Yarborough Z on Sutton pseudo-criticals of ``gsg`` with Lee-Gonzalez-Eakin viscosity
   * - oil_pvt
     - OilPVT
     - Oil PVT object. If provided for oil wells, extracts api, sgsp, pb, rsb
   * - qg_mscfd
     - float
     - Gas rate (Mscf/d, or sm3/d if metric=True). Gas wells only. Gas rates are Mscf/d throughout the library from 3.8.0
   * - qg_mmscfd
     - float
     - Deprecated alias for the gas rate in MMscf/d (sm3/d if metric=True); converted to ``qg_mscfd`` with a ``DeprecationWarning``. Do not pass both
   * - cgr
     - float
     - Condensate-gas ratio (STB/MMscf, or sm3/sm3 if metric=True). Gas wells only
   * - qw_bwpd
     - float
     - Water rate (STB/d, or sm3/d if metric=True). Gas wells only
   * - oil_vis
     - float
     - Oil/condensate viscosity (cP). Defaults to 1.0. Gas wells only
   * - api
     - float
     - Condensate or oil API gravity. Defaults to 45
   * - pr
     - float
     - Reservoir pressure for condensate dropout model (psia, or barsa if metric=True). 0 disables. Gas wells only
   * - qt_stbpd
     - float
     - Total liquid rate (STB/d, or sm3/d if metric=True). Oil wells only
   * - gor
     - float
     - Gas-oil ratio (scf/STB, or sm3/sm3 if metric=True). Oil wells only
   * - wc
     - float
     - Water cut (fraction 0-1). Oil wells only
   * - wsg
     - float
     - Water specific gravity. Defaults to 1.07. Also sets the water-phase viscosity: converted to NaCl salinity with McCain's gamma_w = 1 + 0.695e-6 x ppm (clamped to 0-260,000 ppm; wsg <= 1 is fresh water), then the gas-free brine viscosity ``brine_props`` reports (IAPWS-2008 water, Jones-Dole salt ratio, Kestin pressure factor) at each segment's pressure and temperature
   * - injection
     - bool
     - True for injection wells. Defaults to False
   * - gsg
     - float
     - Gas specific gravity relative to air. Defaults to 0.65
   * - pb
     - float
     - Bubble point pressure (psia, or barsa if metric=True). Required for oil wells
   * - rsb
     - float
     - Solution GOR at Pb (scf/STB, or sm3/sm3 if metric=True). Required for oil wells
   * - sgsp
     - float
     - Separator gas specific gravity. Defaults to 0.65
   * - metric
     - bool
     - If True, interpret inputs and return output in Eclipse METRIC units. Defaults to False

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float
     - Flowing bottom hole pressure (psia, or barsa if metric=True)

Examples:

Gas well:

.. code-block:: python

    >>> from pyrestoolbox import nodal
    >>> c = nodal.Completion(tid=2.441, length=10000, tht=100, bht=200)
    >>> nodal.fbhp(thp=500, completion=c, vlpmethod='HB', well_type='gas', qg_mscfd=5000, gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    961.594453485463

Oil well:

.. code-block:: python

    >>> c = nodal.Completion(tid=2.441, length=8000, tht=100, bht=180)
    >>> nodal.fbhp(thp=200, completion=c, vlpmethod='HB', well_type='oil', qt_stbpd=2000, gor=800, wc=0.3, gsg=0.65, pb=2500, rsb=500, sgsp=0.65, api=35)
    1901.8236509686349

Oil well using OilPVT object:

.. code-block:: python

    >>> from pyrestoolbox import oil
    >>> opvt = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    >>> nodal.fbhp(thp=200, completion=c, vlpmethod='HB', well_type='oil', oil_pvt=opvt, qt_stbpd=2000, gor=800, wc=0.3, gsg=0.65)
    1902.8150784867205

Deviated well using WellSegment:

.. code-block:: python

    >>> segs = [nodal.WellSegment(md=5000, id=2.441, deviation=0), nodal.WellSegment(md=5000, id=2.441, deviation=45)]
    >>> c_dev = nodal.Completion(segments=segs, tht=100, bht=200)
    >>> nodal.fbhp(thp=500, completion=c_dev, vlpmethod='HB', well_type='gas', qg_mscfd=5000, gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    932.8973511781188

.. note::

   Not all VLP methods are equally suitable for deviated wells. For wells with deviation > 30 degrees, prefer **BB** or **WG** which have explicit inclination modelling. See `VLP Method Suitability for Deviated and Horizontal Wells`_ for detailed guidance.

Pressure traverse output
------------------------

Pass ``return_profile=True`` to ``fbhp`` to get the wellbore pressure profile at each segment boundary (in addition to the bottom-hole pressure):

.. code-block:: python

    >>> segs = [nodal.WellSegment(md=5000, id=2.441, deviation=0),
    ...         nodal.WellSegment(md=5000, id=2.441, deviation=45)]
    >>> c = nodal.Completion(segments=segs, tht=100, bht=200)
    >>> profile = nodal.fbhp(thp=500, completion=c, vlpmethod='HB', well_type='gas',
    ...                      qg_mscfd=5000, gsg=0.65, cgr=10, qw_bwpd=10, api=45,
    ...                      oil_vis=1.0, return_profile=True)
    >>> profile['md']         # cumulative MD at each segment boundary
    array([    0.,  5000., 10000.])
    >>> profile['tvd']        # cumulative TVD (respects deviation)
    array([   0.        , 5000.        , 8535.53390593])
    >>> profile['p']          # pressure at each boundary (psia; barsa if metric)
    # array([500., ..., bhp_final])
    >>> profile['bhp']        # scalar BHP, same as return_profile=False


pyrestoolbox.nodal.fthp
=======================

.. code-block:: python

    fthp(bhp, completion, vlpmethod='WG', well_type='gas', ..., thp_min=None, thp_max=None, tol=1e-3) -> float

Inverse of ``fbhp``. Solves for the tubing head pressure that produces the target flowing bottom-hole pressure, for the given flow parameters and VLP correlation. Internally wraps ``bisect_solve`` around ``fbhp`` between ``thp_min`` and ``thp_max``. ``thp_min=None`` resolves to 14.7 psia (1.01325 barsa if metric=True); ``thp_max=None`` resolves to 20,000 psia (1378.95 barsa if metric=True). ``tol`` is the absolute THP tolerance in psia (barsa if metric=True).

All flow / well / unit kwargs accepted by ``fbhp`` (``qg_mscfd``, ``qt_stbpd``, ``gor``, ``wc``, ``api``, ``pb``, ``rsb``, ``metric`` ...) are accepted here identically. ``metric=True`` switches every pressure, rate and ratio input and the return to Eclipse METRIC units, exactly as in ``fbhp``.

.. code-block:: python

    >>> c = nodal.Completion(tid=2.441, length=8000, tht=100, bht=180)
    >>> thp_target = 800.0
    >>> bhp = nodal.fbhp(thp=thp_target, completion=c, vlpmethod='HB', well_type='gas',
    ...                  qg_mscfd=3000, gsg=0.65, cgr=0, qw_bwpd=0, api=45, oil_vis=1.0)
    >>> thp_recovered = nodal.fthp(bhp=bhp, completion=c, vlpmethod='HB', well_type='gas',
    ...                             qg_mscfd=3000, gsg=0.65, cgr=0, qw_bwpd=0, api=45, oil_vis=1.0)
    >>> abs(thp_recovered - thp_target) < 1.0
    True

Useful for: wellhead choke sizing (target BHP above the reservoir is known, need surface pressure), back-calculating THP from downhole gauge data, or matching shut-in surface pressures against a downhole model.


pyrestoolbox.nodal.outflow_curve
================================

.. code-block:: python

    outflow_curve(thp, completion, vlpmethod='WG', well_type='gas', gas_pvt=None, oil_pvt=None, rates=None, n_points=20, max_rate=None, cgr=0, qw_bwpd=0, oil_vis=1.0, api=45, pr=0, gor=0, wc=0, wsg=1.07, injection=False, gsg=0.65, pb=0, rsb=0, sgsp=0.65, metric=False, n_rates=None) -> dict

Returns VLP outflow curve as a dictionary with keys ``'rate'`` (also aliased as ``'rates'``) and ``'bhp'``. ``n_rates`` is a deprecated alias for ``n_points`` (takes precedence if both are given). Evaluates ``fbhp()`` at each rate point. Rates are Mscf/d for gas wells (sm3/d if metric=True), STB/d for oil wells (sm3/d if metric=True). BHP is in psia (barsa if metric=True).

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - thp
     - float
     - Tubing head pressure (psia, or barsa if metric=True)
   * - completion
     - Completion
     - Completion object
   * - vlpmethod
     - string or vlp_method
     - VLP method. `Calculation Methods and Class Objects`_
   * - well_type
     - string
     - 'gas' or 'oil'
   * - rates
     - list
     - Explicit list of rates to evaluate (same units as output rates). If None, auto-generated
   * - n_points
     - int
     - Number of rate points if rates is None. Defaults to 20. ``n_rates`` is kept as a deprecated alias for backward compatibility.
   * - max_rate
     - float
     - Maximum rate for auto-generation. Defaults to 50,000 Mscf/d (gas) or 10,000 STB/d (oil)
   * - metric
     - bool
     - If True, interpret inputs and return output in Eclipse METRIC units. Defaults to False
   * - (other)
     -
     - Same parameters as ``fbhp()``

.. list-table:: Returns (dict)
   :widths: 10 15 40
   :header-rows: 1

   * - Key
     - Type
     - Description
   * - 'rate' (alias 'rates')
     - list
     - Rate values (Mscf/d for gas, STB/d for oil; sm3/d if metric). Gas was MMscf/d before 3.8.0
   * - 'bhp'
     - list
     - Flowing BHP at each rate (psia, or barsa if metric=True)

Examples:

.. code-block:: python

    >>> c = nodal.Completion(tid=2.441, length=10000, tht=100, bht=200)
    >>> result = nodal.outflow_curve(thp=500, completion=c, vlpmethod='HB', well_type='gas', rates=[2000, 5000, 10000, 15000, 20000], gsg=0.65)
    >>> result['rates']
    [2000, 5000, 10000, 15000, 20000]
    >>> [round(b, 1) for b in result['bhp']]
    [681.6, 935.2, 1512.3, 2138.0, 2777.1]


pyrestoolbox.nodal.ipr_curve
============================

.. code-block:: python

    ipr_curve(reservoir, well_type='gas', gas_pvt=None, oil_pvt=None, n_points=20, min_pwf=None, wc=0, wsg=1.07, bo=1.2, uo=1.0, gsg=0.65, metric=False) -> dict

Returns IPR (Inflow Performance Relationship) curve as a dictionary with keys ``'pwf'`` and ``'rate'``. Pressures are in psia (barsa if metric=True), rates in Mscf/d for gas (sm3/d if metric=True) or STB/d for oil (sm3/d if metric=True).

For gas wells, uses pseudopressure deliverability via ``gas.gas_rate_radial()``. Returns rates in Mscf/d.

For oil wells with OilPVT: uses Darcy above Pb, Vogel below Pb. Without OilPVT: simple Darcy. Returns rates in STB/d.

.. note::

   Gas rates are Mscf/d everywhere from 3.8.0: ``ipr_curve()``, ``outflow_curve()``, ``operating_point()`` (its ``rate``, ``vlp`` and ``ipr``) and ``fbhp``/``fthp`` (``qg_mscfd``) share one unit, so IPR and VLP curves plot on the same axis directly. Before 3.8.0 the VLP side used MMscf/d; ``qg_mmscfd`` remains as a deprecated alias.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - reservoir
     - Reservoir
     - Reservoir object
   * - well_type
     - string
     - 'gas', 'oil', or 'water'
   * - gas_pvt
     - GasPVT
     - Gas PVT object. Used for gas well composition (sg, co2, h2s, n2, h2, he)
   * - oil_pvt
     - OilPVT
     - Oil PVT object. Used for oil well Pb, Rs, Bo, viscosity
   * - n_points
     - int
     - Number of pressure points. Defaults to 20
   * - min_pwf
     - float
     - Minimum flowing BHP (psia, or barsa if metric=True). Defaults to 14.7 psia / 1.01325 barsa
   * - wc
     - float
     - Water cut (fraction 0-1). For oil wells
   * - wsg
     - float
     - Water specific gravity. Defaults to 1.07
   * - bo
     - float
     - Oil FVF (rb/STB). Used if oil_pvt not provided. Defaults to 1.2
   * - uo
     - float
     - Oil viscosity (cP). Used if oil_pvt not provided. Defaults to 1.0
   * - gsg
     - float
     - Gas specific gravity. Used if gas_pvt not provided. Defaults to 0.65
   * - metric
     - bool
     - If True, interpret inputs and return output in Eclipse METRIC units. Defaults to False

.. list-table:: Returns (dict)
   :widths: 10 15 40
   :header-rows: 1

   * - Key
     - Type
     - Description
   * - 'pwf'
     - list
     - Flowing BHP values (psia, or barsa if metric=True)
   * - 'rate'
     - list
     - Rate values (Mscf/d for gas, STB/d for oil; sm3/d if metric)

Examples:

.. code-block:: python

    >>> r = nodal.Reservoir(pr=3000, degf=200, k=10, h=50, re=1500, rw=0.35, S=2, D=0.001)
    >>> ipr = nodal.ipr_curve(r, well_type='gas', gsg=0.65, n_points=5)
    >>> [round(p, 1) for p in ipr['pwf']]
    [14.7, 761.0, 1507.4, 2253.7, 3000.0]
    >>> [round(q, 1) for q in ipr['rate']]
    [13458.9, 12814.7, 10863.4, 7292.7, 0.0]


pyrestoolbox.nodal.operating_point
==================================

.. code-block:: python

    operating_point(thp, completion, reservoir, vlpmethod='WG', well_type='gas', gas_pvt=None, oil_pvt=None, cgr=0, qw_bwpd=0, oil_vis=1.0, api=45, gor=0, wc=0, wsg=1.07, injection=False, gsg=0.65, pb=0, rsb=0, sgsp=0.65, bo=1.2, uo=1.0, n_points=25, metric=False) -> dict

Finds the operating point where VLP outflow curve intersects the IPR inflow curve via bisection.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - thp
     - float
     - Tubing head pressure (psia, or barsa if metric=True)
   * - completion
     - Completion
     - Completion object
   * - reservoir
     - Reservoir
     - Reservoir object
   * - vlpmethod
     - string or vlp_method
     - VLP method. `Calculation Methods and Class Objects`_
   * - well_type
     - string
     - 'gas' or 'oil'
   * - n_points
     - int
     - Number of points for VLP and IPR curves. Defaults to 25
   * - metric
     - bool
     - If True, all inputs and outputs use Eclipse METRIC units. Defaults to False
   * - (other)
     -
     - Same parameters as ``fbhp()`` and ``ipr_curve()``

.. list-table:: Returns (dict)
   :widths: 10 15 40
   :header-rows: 1

   * - Key
     - Type
     - Description
   * - 'rate'
     - float
     - Operating rate (Mscf/d for gas, STB/d for oil; sm3/d if metric). Gas was MMscf/d before 3.8.0
   * - 'bhp'
     - float
     - Operating flowing BHP (psia, or barsa if metric=True)
   * - 'vlp'
     - dict
     - VLP curve: {'rates': [...], 'bhp': [...]}
   * - 'ipr'
     - dict
     - IPR curve: {'pwf': [...], 'rate': [...]}
   * - 'converged'
     - bool
     - True if the VLP/IPR intersection solve converged

Examples:

Gas well operating point:

.. code-block:: python

    >>> c = nodal.Completion(tid=2.441, length=10000, tht=100, bht=200)
    >>> r = nodal.Reservoir(pr=3000, degf=200, k=10, h=50, re=1500, rw=0.35, S=2, D=0.001)
    >>> result = nodal.operating_point(thp=500, completion=c, reservoir=r, vlpmethod='HB', well_type='gas', gsg=0.65)
    >>> round(result['rate'], 1)
    10578.1
    >>> round(result['bhp'], 1)
    1583.4

Oil well operating point:

.. code-block:: python

    >>> from pyrestoolbox import oil
    >>> c = nodal.Completion(tid=2.441, length=8000, tht=100, bht=180)
    >>> r = nodal.Reservoir(pr=3000, degf=180, k=50, h=30, re=1000, rw=0.35)
    >>> opvt = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    >>> result = nodal.operating_point(thp=200, completion=c, reservoir=r, vlpmethod='HB', well_type='oil', oil_pvt=opvt, gor=800, wc=0.3, gsg=0.65)
    >>> round(result['rate'], 1)   # total liquid, STB/d (oil 1643.9 at 30% water cut)
    2348.5
    >>> round(result['bhp'], 1)
    2013.1


pyrestoolbox.nodal.critical_gas_rate
====================================

.. code-block:: python

    critical_gas_rate(p, degf, sg, tid, method='TURNER', liquid='water', rho_liq=None, ift=None, co2=0, h2s=0, n2=0, h2=0, zmethod='DAK', cmethod='PMC', metric=False) -> dict

Minimum gas rate to lift liquid droplets out of the tubing (onset of liquid loading), from the droplet model v = K [ift (rho_liq - rho_g)]^0.25 / rho_g^0.5 (ft/s, dyne/cm, lbm/ft3). ``'TURNER'`` uses K = 1.912, Turner, Hubbard & Dukler (1969) with their 20% uplift; ``'COLEMAN'`` uses the unadjusted K = 1.593, Coleman et al. (1991), who found the uplift unnecessary at wellhead pressures below about 500 psia. The liquid presets are Turner's: water 67 lbm/ft3 and 60 dyne/cm, condensate 45 lbm/ft3 and 20 dyne/cm. Evaluate at the wellhead and at the bottom of the tubing; the larger rate controls. ``p`` may be a scalar or array; outputs follow it.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float, list or array
     - Flowing pressure (psia | barsa)
   * - degf
     - float
     - Flowing temperature (deg F | deg C)
   * - sg
     - float
     - Gas SG relative to air
   * - tid
     - float
     - Tubing inside diameter (in | mm)
   * - method
     - string
     - 'TURNER' (default) or 'COLEMAN'
   * - liquid
     - string
     - 'water' (default) or 'condensate'; sets the default liquid density and interfacial tension
   * - rho_liq
     - float
     - Liquid density, overrides the preset (lbm/ft3 | kg/m3)
   * - ift
     - float
     - Gas-liquid interfacial tension, overrides the preset (dyne/cm | mN/m)
   * - co2, h2s, n2, h2
     - float
     - Mole fractions in the gas. Default 0
   * - zmethod, cmethod
     - string
     - Z-factor and critical-property methods, as for ``gas.gas_z``. Defaults 'DAK', 'PMC'
   * - metric
     - bool
     - If True, inputs and outputs in Eclipse METRIC units. Defaults to False

.. list-table:: Returns (dict)
   :widths: 10 15 40
   :header-rows: 1

   * - Key
     - Type
     - Description
   * - 'v_crit'
     - float or array
     - Critical velocity (ft/s | m/s)
   * - 'q_crit'
     - float or array
     - Critical gas rate (Mscf/d | sm3/d)
   * - 'z'
     - float or array
     - Gas Z-factor
   * - 'rho_gas'
     - float or array
     - Gas density (lbm/ft3 | kg/m3)

Examples:

.. code-block:: python

    >>> from pyrestoolbox import nodal
    >>> r = nodal.critical_gas_rate([150, 1200], degf=120, sg=0.65, tid=2.441)
    >>> [round(q, 1) for q in r['q_crit']]
    [584.4, 1722.4]
    >>> round(nodal.critical_gas_rate(150, 120, 0.65, 2.441, method='COLEMAN', liquid='condensate')['q_crit'], 1)
    334.7


pyrestoolbox.nodal.erosional_velocity
=====================================

.. code-block:: python

    erosional_velocity(p, degf, sg, glr=None, sg_liq=1.0, c=100.0, tid=None, co2=0, h2s=0, n2=0, h2=0, zmethod='DAK', cmethod='PMC', metric=False) -> dict

API RP 14E erosional velocity, Ve = C / sqrt(rho_m). With ``glr`` given, rho_m is the RP 14E no-slip mixture density, (12409 S_l P + 2.7 R S_g P) / (198.7 P + R T Z); with ``glr=None`` (dry gas) it is the gas density. ``c`` is always in field units, (lbm/ft s^2)^0.5: RP 14E gives 100 for continuous and 125 for intermittent service, and 150-200 is used for solids-free, non-corrosive service. It is a screening limit, not an erosion-rate prediction. When ``tid`` is given the maximum rate at Ve is also returned: liquid rate from the RP 14E minimum-area relation, A = (9.35 + Z R T / (21.25 P)) / Ve in2 per 1000 bbl/d, or gas rate for dry gas.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float, list or array
     - Flowing pressure (psia | barsa)
   * - degf
     - float
     - Flowing temperature (deg F | deg C)
   * - sg
     - float
     - Gas SG relative to air
   * - glr
     - float
     - Gas-liquid ratio at standard conditions (scf/STB | sm3/sm3). None (default) for dry gas
   * - sg_liq
     - float
     - Liquid SG relative to water. Defaults to 1.0
   * - c
     - float
     - RP 14E C-factor (field units in both unit systems). Defaults to 100
   * - tid
     - float
     - Optional pipe inside diameter (in | mm)
   * - co2, h2s, n2, h2
     - float
     - Mole fractions in the gas. Default 0
   * - zmethod, cmethod
     - string
     - Z-factor and critical-property methods, as for ``gas.gas_z``
   * - metric
     - bool
     - If True, inputs and outputs in Eclipse METRIC units. Defaults to False

.. list-table:: Returns (dict)
   :widths: 10 15 40
   :header-rows: 1

   * - Key
     - Type
     - Description
   * - 've'
     - float or array
     - Erosional velocity (ft/s | m/s)
   * - 'rho_mix'
     - float or array
     - Mixture (or gas) density (lbm/ft3 | kg/m3)
   * - 'z'
     - float or array
     - Gas Z-factor
   * - 'q_liq_max'
     - float or array
     - Liquid rate at Ve (STB/d | sm3/d); only when tid and glr are given
   * - 'q_gas_max'
     - float or array
     - Gas rate at Ve (Mscf/d | sm3/d); only when tid is given and glr is None

Examples:

.. code-block:: python

    >>> r = nodal.erosional_velocity(1000, degf=140, sg=0.7, glr=1500, tid=2.441)
    >>> round(r['ve'], 2), round(r['q_liq_max'], 1)
    (25.64, 2544.8)
    >>> round(nodal.erosional_velocity(3000, degf=180, sg=0.7, tid=3.958)['q_gas_max'], 1)
    43830.4


pyrestoolbox.nodal.choke_flow
=============================

.. code-block:: python

    choke_flow(glr, p_up=None, q_liq=None, size=None, method='GILBERT', p_down=None, metric=False) -> dict

Gilbert-family critical-flow choke correlations, p_up = A q GLR^B / S^C with p_up in psig, q gross liquid in STB/d, GLR in scf/STB and S the bean size in 64ths of an inch. Give exactly two of ``p_up``, ``q_liq`` and ``size``; the third is solved. These are empirical fits to critical (choked) flow of oil wells and take no account of downstream pressure. When ``p_down`` is given and p_down / p_up exceeds 1/1.7 (Gilbert's criterion, absolute pressures), the flow is subcritical, a warning is issued and ``'critical'`` is False.

.. list-table:: Coefficients (GLR in scf/STB)
   :widths: 15 10 10 10
   :header-rows: 1

   * - Method
     - A
     - B
     - C
   * - 'GILBERT' (1954)
     - 10.00
     - 0.546
     - 1.89
   * - 'ROS' (1960)
     - 17.40
     - 0.500
     - 2.00
   * - 'BAXENDELL' (1958)
     - 9.56
     - 0.546
     - 1.93
   * - 'ACHONG' (1961)
     - 3.82
     - 0.650
     - 1.88

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - glr
     - float
     - Gross gas-liquid ratio (scf/STB | sm3/sm3)
   * - p_up
     - float
     - Upstream (wellhead) pressure (psia | barsa)
   * - q_liq
     - float
     - Gross liquid rate (STB/d | sm3/d)
   * - size
     - float
     - Bean size in 64ths of an inch (both unit systems)
   * - method
     - string
     - 'GILBERT' (default), 'ROS', 'BAXENDELL' or 'ACHONG'
   * - p_down
     - float
     - Optional downstream pressure for the critical-flow check (psia | barsa)
   * - metric
     - bool
     - If True, inputs and outputs in Eclipse METRIC units. Defaults to False

.. list-table:: Returns (dict)
   :widths: 10 15 40
   :header-rows: 1

   * - Key
     - Type
     - Description
   * - 'p_up'
     - float
     - Upstream pressure (psia | barsa)
   * - 'q_liq'
     - float
     - Gross liquid rate (STB/d | sm3/d)
   * - 'size'
     - float
     - Bean size (64ths of an inch)
   * - 'critical'
     - bool or None
     - Critical-flow check result; None when p_down is not given

Examples:

.. code-block:: python

    >>> round(nodal.choke_flow(glr=1200, q_liq=800, size=24)['p_up'], 1)
    960.3
    >>> round(nodal.choke_flow(glr=1200, p_up=450, q_liq=800, method='ACHONG')['size'], 2)
    32.72


pyrestoolbox.nodal.choke_flow_sachdeva
======================================

.. code-block:: python

    choke_flow_sachdeva(p_down, degf, glr, p_up=None, q_liq=None, size=None, sg=0.65,
                        sg_liq=1.0, cd=0.75, k=1.3, cl=None, co2=0, h2s=0, n2=0, h2=0,
                        zmethod='DAK', cmethod='PMC', metric=False) -> dict

Mechanistic two-phase choke model of Sachdeva, Schmidt, Brill & Blais (1986), Two-Phase Flow Through Chokes, SPE-15657-MS (doi 10.2118/15657-MS), valid in both critical and subcritical flow. Give exactly two of ``p_up``, ``q_liq`` and ``size``; the third is solved. Unlike the Gilbert family it sees the downstream pressure and the fluid properties.

The critical pressure ratio y_c solves the paper's Eq. 2 (iterated, since the throat gas volume depends on y), with the polytropic exponent n = 1 + x (Cp - Cv) / (x Cv + (1 - x) CL). When p_down / p_up <= y_c the flow is critical and the throat ratio is y_c; otherwise it is subcritical at the actual ratio. The mass flux is

.. code-block:: text

    G2 = CD {2 gc 144 P1 rho_m2^2 [(1 - x)(1 - y) VL + x k / (k - 1) (VG1 - y VG2)]}^0.5
    VG2 = VG1 y^(-1/k),   1 / rho_m2 = x VG2 + (1 - x) VL

with x the free-gas mass quality, V specific volumes (ft3/lbm) and P1 in psia. As in the paper, only free gas is considered (``glr`` is the free-gas ratio and ``q_liq`` the liquid rate at upstream conditions), and the gas expands adiabatically with an ideal-gas exponent ``k``. The y_c(x) curves of the paper's Fig. 5 are reproduced to within 0.025 (the reading resolution) across its five cases; in the pure-gas limit y_c is (2 / (k + 1))^(k / (k - 1)). Note that y_c from Eq. 2 is not exactly the ratio that maximises the Eq. 4 flux (the boundary uses the polytropic n, the flux the adiabatic k); the model is implemented as published.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p_down
     - float
     - Downstream pressure (psia | barsa)
   * - degf
     - float
     - Upstream temperature (deg F | deg C)
   * - glr
     - float
     - Free gas-liquid ratio (scf/STB | sm3/sm3), > 0
   * - p_up, q_liq, size
     - float
     - Upstream pressure (psia | barsa), liquid rate (STB/d | sm3/d), choke diameter in 64ths of an inch (both unit systems). Give two
   * - sg, sg_liq
     - float
     - Gas SG (air = 1) and liquid SG (water = 1). Defaults 0.65 and 1.0
   * - cd
     - float
     - Discharge coefficient: 0.75 (default) for a choke in a housing with an upstream elbow, 0.85 for undisturbed upstream flow
   * - k
     - float
     - Gas heat-capacity ratio Cp/Cv. Default 1.3
   * - cl
     - float
     - Liquid specific heat (Btu/(lbm degF) | kJ/(kg K)). Default 0.8 Btu/(lbm degF), the paper's base case; water about 1.0, crude about 0.5
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
   :widths: 10 40
   :header-rows: 1

   * - Key
     - Description
   * - 'p_up', 'q_liq', 'size'
     - Upstream pressure (psia | barsa), liquid rate (STB/d | sm3/d) and choke size (64ths)
   * - 'q_gas'
     - Free-gas rate (Mscf/d | sm3/d)
   * - 'critical'
     - True when p_down / p_up <= y_c
   * - 'y_crit', 'y'
     - Critical and actual downstream/upstream pressure ratios
   * - 'x'
     - Free-gas mass quality

Examples:

.. code-block:: python

    >>> r = nodal.choke_flow_sachdeva(p_down=100, degf=120, glr=800, p_up=400, size=32, sg=0.7, sg_liq=0.85)
    >>> round(r['q_liq'], 1), r['critical'], round(r['y_crit'], 4)
    (742.4, True, 0.5929)
    >>> round(nodal.choke_flow_sachdeva(p_down=350, degf=120, glr=800, p_up=400, size=32, sg=0.7, sg_liq=0.85)['q_liq'], 1)
    503.6
    >>> round(nodal.choke_flow_sachdeva(p_down=100, degf=120, glr=800, q_liq=1000, size=32, sg=0.7, sg_liq=0.85)['p_up'], 2)
    540.78


Calculation Methods and Class Objects
=====================================
pyResToolBox uses class objects to track calculation options through the functions. Class objects can be set via strings or explicitly via object options

.. list-table:: Method Variables & Class Objects
   :widths: 10 15 40
   :header-rows: 1

   * - Class Variable
     - Class Object
     - Class Description & Options
   * - vlpmethod
     - vlp_method
     - Multiphase flow correlation for VLP calculations. Defaults to 'WG'.
       Options are:

       + 'HB': Hagedorn-Brown (1965) with Orkiszewski bubble flow correction
       + 'WG': Woldesemayat-Ghajar (2007) drift-flux model
       + 'GRAY': Gray (1978) effective roughness with acceleration term
       + 'BB': Beggs & Brill (1973), with the 0.924 uphill holdup factor usually credited to Payne et al. (1979); injection uses Beggs & Brill's downhill inclination correction

Users can specify which calculation method to use either by passing an option string, or a class object to any given function.

Examples:

.. code-block:: python

    >>> from pyrestoolbox import nodal
    >>> c = nodal.Completion(tid=2.441, length=10000, tht=100, bht=200)
    >>> nodal.fbhp(thp=500, completion=c, vlpmethod='HB', well_type='gas', qg_mscfd=5000, gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    961.594453485463

Comparing all four VLP methods for the same gas well:

.. code-block:: python

    >>> nodal.fbhp(thp=500, completion=c, vlpmethod='HB', well_type='gas', qg_mscfd=5000, gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    961.594453485463
    >>> nodal.fbhp(thp=500, completion=c, vlpmethod='WG', well_type='gas', qg_mscfd=5000, gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    1127.2904240222729
    >>> nodal.fbhp(thp=500, completion=c, vlpmethod='GRAY', well_type='gas', qg_mscfd=5000, gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    1066.755704999179
    >>> nodal.fbhp(thp=500, completion=c, vlpmethod='BB', well_type='gas', qg_mscfd=5000, gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    1224.730727621603


VLP Method Suitability for Deviated and Horizontal Wells
========================================================

.. warning::

   The four VLP correlations implemented here were originally developed and validated against vertical or near-vertical well data. When applied to high-angle or horizontal wellbores, their accuracy varies significantly. Users should understand these limitations before relying on results for deviated completions.

The deviation support in pyResToolbox applies a ``sin(theta)`` correction to the hydrostatic gradient term in each correlation and, where applicable, passes the inclination angle through to holdup and flow-pattern sub-models. This is the standard first-order treatment for inclined pipe, but it does not address all the physics that change at high angles.

.. list-table:: Method Suitability by Inclination
   :widths: 12 12 40
   :header-rows: 1

   * - Method
     - Deviation Range
     - Notes
   * - **BB** (Beggs & Brill)
     - 0 -- 90 deg
     - The only correlation in this set originally developed for all inclinations. The inclination correction (Beggs and Brill 1973) was fitted to data from -90 to +90 degrees, with one coefficient set per flow pattern uphill and a single set for all downhill flow; injection wells (fluid flowing down the tubing) use the downhill set. Preferred first choice for deviated and horizontal wells.
   * - **WG** (Woldesemayat-Ghajar)
     - 0 -- 90 deg
     - Drift-flux formulation with an inclination-dependent drift velocity term. The void-fraction model explicitly includes ``sin(theta)`` and ``cos(theta)`` terms derived from a broad multi-angle dataset. Suitable across the full inclination range.
   * - **HB** (Hagedorn-Brown)
     - 0 -- ~30 deg
     - Developed from vertical well data only. The holdup correlation charts were fitted to vertical upflow experiments. At moderate deviations the ``sin(theta)`` hydrostatic correction provides a reasonable approximation, but the holdup model itself has no inclination dependence. Use with caution beyond ~30 degrees from vertical; results become increasingly unreliable at high angles.
   * - **GRAY** (Gray)
     - 0 -- ~30 deg
     - Developed for vertical gas-condensate wells. The effective roughness and liquid holdup models assume vertical counter-current film flow. Like HB, only the hydrostatic term is corrected for inclination. Not recommended for wells with significant deviation.

**Recommendations for high-angle wells:**

- For wells with deviation > 30 degrees, prefer **BB** or **WG**. Both have explicit inclination modelling in their holdup / void-fraction formulations.
- For near-horizontal wells (deviation > 70 degrees), **BB** is the most widely validated choice. Flow-pattern transitions change fundamentally in horizontal pipe (stratified flow becomes dominant), and the Beggs & Brill flow-pattern map was designed to capture this.
- For vertical and near-vertical gas and gas-condensate wells (deviation < 30 degrees), all four methods are appropriate. HB and Gray were specifically developed for this regime and may give more accurate results than BB. Gray (API 14B) is a gas-well method: it is accepted for ``well_type='oil'`` for completeness but has no basis for oil or water-dominated flow, so prefer HB, BB or WG for oil wells.
- When in doubt, run multiple methods and compare. Large disagreement between BB/WG and HB/Gray at high angles is expected and indicates that the inclination-insensitive correlations are outside their validation envelope.
- These correlations should not be used as a substitute for a calibrated mechanistic multiphase model when accurate pressure predictions are required for field development decisions in deviated or horizontal wells.


Unit System Support
===================

All nodal classes and functions accept a ``metric=False`` parameter. When set to ``True``, inputs and outputs use Eclipse METRIC units:

.. list-table:: Unit Convention
   :widths: 15 20 20
   :header-rows: 1

   * - Quantity
     - Field (metric=False)
     - Metric (metric=True)
   * - Pressure
     - psia
     - barsa
   * - Temperature
     - deg F
     - deg C
   * - Length / Depth
     - ft
     - m
   * - Pipe diameter
     - inches
     - mm
   * - Pipe roughness
     - inches
     - mm
   * - Gas rate
     - Mscf/d
     - sm3/d
   * - Liquid rate
     - STB/d
     - sm3/d
   * - CGR
     - STB/MMscf
     - sm3/sm3
   * - GOR
     - scf/STB
     - sm3/sm3
   * - Non-Darcy D
     - day/mscf
     - day/sm3

The data classes (``WellSegment``, ``Completion``, ``Reservoir``) convert metric inputs to oilfield units internally, so all calculation engines operate in field units regardless of the user-facing unit system. Standard conditions in metric mode are 1.01325 barsa and 15.0 deg C.


PVT Wrapper Objects
===================

``GasPVT`` and ``OilPVT`` are convenience classes defined in the ``gas`` and ``oil`` modules respectively. See the `gas module documentation <gas.rst>`_ and `oil module documentation <oil.rst>`_ for full constructor signatures, parameter tables, and examples.

They can be passed directly to ``fbhp()``, ``outflow_curve()``, ``ipr_curve()``, and ``operating_point()`` for property evaluation. They can also be passed to ``gas_rate_radial()``, ``gas_rate_linear()``, ``oil_rate_radial()``, and ``oil_rate_linear()`` for IPR rate calculations.
