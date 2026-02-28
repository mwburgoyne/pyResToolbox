===================================
Nodal Analysis & VLP
===================================

Vertical Lift Performance (VLP) calculations, Inflow Performance Relationships (IPR), and nodal operating point analysis for gas and oil wells.

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
        + 'BB': Beggs & Brill (1973) with Payne et al. correction

Users can specify which calculation method to use either by passing an option string, or a class object to any given function.

Examples:

.. code-block:: python

    >>> from pyrestoolbox import nodal
    >>> c = nodal.Completion(tid=2.441, length=10000, tht=100, bht=200)
    >>> nodal.fbhp(thp=500, completion=c, vlpmethod='HB', well_type='gas', qg_mmscfd=5.0, gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    952.6868477414688

Comparing all four VLP methods for the same gas well:

.. code-block:: python

    >>> nodal.fbhp(thp=500, completion=c, vlpmethod='HB', well_type='gas', qg_mmscfd=5.0, gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    952.6868477414688
    >>> nodal.fbhp(thp=500, completion=c, vlpmethod='WG', well_type='gas', qg_mmscfd=5.0, gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    1172.8626065704736
    >>> nodal.fbhp(thp=500, completion=c, vlpmethod='GRAY', well_type='gas', qg_mmscfd=5.0, gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    1940.034905804135
    >>> nodal.fbhp(thp=500, completion=c, vlpmethod='BB', well_type='gas', qg_mmscfd=5.0, gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    1213.2399909265969


VLP Method Suitability for Deviated and Horizontal Wells
=========================================================

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
     - The only correlation in this set originally developed for all inclinations. The flow-pattern map and inclination correction factor (Payne et al. 1979) were fitted to data covering the full range from vertical to horizontal. Preferred first choice for deviated and horizontal wells.
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
- For vertical and near-vertical wells (deviation < 30 degrees), all four methods are appropriate. HB and Gray were specifically developed for this regime and may give more accurate results than BB for gas and gas-condensate wells.
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
     - MMscf/d
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

GasPVT and OilPVT are convenience classes that store fluid characterization parameters and method choices, pre-compute critical properties, and expose simple property methods. They are defined in the ``gas`` and ``oil`` modules respectively.

pyrestoolbox.gas.GasPVT
------------------------

.. code-block:: python

    GasPVT(sg=0.75, co2=0, h2s=0, n2=0, h2=0, zmethod='DAK', cmethod='PMC', metric=False)

Stores gas composition and method choices. Pre-computes critical temperature and pressure via ``gas_tc_pc()`` so they are not recalculated per call. Automatically selects BNS zmethod/cmethod when h2 > 0.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - sg
     - float
     - Gas SG relative to air. Defaults to 0.75
   * - co2
     - float
     - Molar fraction of CO2. Defaults to 0
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to 0
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to 0
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to 0. If positive, overrides zmethod and cmethod to 'BNS'
   * - zmethod
     - string or z_method
     - Method for Z-factor calculation. Defaults to 'DAK'
   * - cmethod
     - string or c_method
     - Method for critical property calculation. Defaults to 'PMC'
   * - metric
     - bool
     - If True, methods accept/return Eclipse METRIC units (barsa, deg C, kg/m3, rm3/sm3). Defaults to False

.. list-table:: Methods
   :widths: 15 40
   :header-rows: 1

   * - Method
     - Description
   * - ``z(p, degf)``
     - Returns Z-factor at pressure p (psia, or barsa if metric=True) and temperature degf (deg F, or deg C if metric=True)
   * - ``viscosity(p, degf)``
     - Returns gas viscosity (cP)
   * - ``density(p, degf)``
     - Returns gas density (lb/cuft, or kg/m3 if metric=True)
   * - ``bg(p, degf)``
     - Returns gas formation volume factor (rcf/scf, or rm3/sm3 if metric=True)

Examples:

.. code-block:: python

    >>> from pyrestoolbox import gas
    >>> gpvt = gas.GasPVT(sg=0.65, co2=0.1)
    >>> gpvt.z(2000, 180)
    0.9026719828498643
    >>> gpvt.viscosity(2000, 180)
    0.016666761192334678
    >>> gpvt.density(2000, 180)
    6.077743278791424
    >>> gpvt.bg(2000, 180)
    0.008164459661048012

Using metric units (pressure in barsa, temperature in deg C):

.. code-block:: python

    >>> gpvt_m = gas.GasPVT(sg=0.65, co2=0.1, metric=True)
    >>> gpvt_m.z(137.9, 82.2)
    0.9026433033953588
    >>> gpvt_m.density(137.9, 82.2)
    97.36871728783241


pyrestoolbox.oil.OilPVT
------------------------

.. code-block:: python

    OilPVT(api, sg_sp, pb, rsb, sg_g=0, rsmethod='VELAR', pbmethod='VALMC', bomethod='MCAIN', metric=False)

Stores oil characterization parameters and method choices. Computes ``sg_o`` from API in constructor. Can be passed directly to ``fbhp()`` and ``operating_point()`` for oil wells.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - api
     - float
     - Stock tank oil density (deg API)
   * - sg_sp
     - float
     - Separator gas specific gravity (relative to air)
   * - pb
     - float
     - Bubble point pressure (psia, or barsa if metric=True)
   * - rsb
     - float
     - Solution GOR at Pb (scf/STB, or sm3/sm3 if metric=True)
   * - sg_g
     - float
     - Weighted average surface gas SG. Estimated from sg_sp if not provided
   * - rsmethod
     - string or rs_method
     - Method for Rs calculation. Defaults to 'VELAR'
   * - pbmethod
     - string or pb_method
     - Method for Pb calculation. Defaults to 'VALMC'
   * - bomethod
     - string or bo_method
     - Method for Bo calculation. Defaults to 'MCAIN'
   * - metric
     - bool
     - If True, constructor inputs (pb, rsb) and method inputs/outputs use Eclipse METRIC units. Defaults to False

.. list-table:: Methods
   :widths: 15 40
   :header-rows: 1

   * - Method
     - Description
   * - ``rs(p, degf)``
     - Returns solution GOR (scf/STB, or sm3/sm3 if metric=True) at pressure p (psia, or barsa if metric=True) and temperature degf (deg F, or deg C if metric=True)
   * - ``bo(p, degf, rs=None)``
     - Returns oil FVF (rb/STB, or rm3/sm3 if metric=True). Optionally pass pre-calculated rs to avoid redundant calculation
   * - ``density(p, degf, rs=None)``
     - Returns live oil density (lb/cuft, or kg/m3 if metric=True)
   * - ``viscosity(p, degf, rs=None)``
     - Returns oil viscosity (cP)

Examples:

.. code-block:: python

    >>> from pyrestoolbox import oil
    >>> opvt = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    >>> opvt.rs(2000, 180)
    403.58333168879415
    >>> opvt.bo(2000, 180)
    1.22370082673546
    >>> opvt.density(2000, 180)
    46.23700811760461
    >>> opvt.viscosity(2000, 180)
    0.7187504436478858

Using metric units (pb in barsa, rsb in sm3/sm3):

.. code-block:: python

    >>> opvt_m = oil.OilPVT(api=35, sg_sp=0.65, pb=172.4, rsb=89, metric=True)
    >>> opvt_m.rs(137.9, 82.2)
    71.82727018664512
    >>> opvt_m.density(137.9, 82.2)
    740.7086089268661


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
   * - VLP Outflow Curve
     - `pyrestoolbox.nodal.outflow_curve`_
   * - IPR Inflow Curve
     - `pyrestoolbox.nodal.ipr_curve`_
   * - Nodal Operating Point
     - `pyrestoolbox.nodal.operating_point`_


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

**Segment mode** (``segments`` keyword): Accepts a list of ``WellSegment`` objects for arbitrary multi-segment wellbore definitions with deviation support. Temperature is interpolated linearly over the total measured depth. In segment mode, ``metric`` only converts ``tht`` and ``bht``; segment dimensions are handled by each ``WellSegment``'s own ``metric`` flag.

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
--------------------------

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
-------------------

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
========================

.. code-block:: python

    fbhp(thp, completion, vlpmethod='WG', well_type='gas', gas_pvt=None, oil_pvt=None, qg_mmscfd=0, cgr=0, qw_bwpd=0, oil_vis=1.0, api=45, pr=0, qt_stbpd=0, gor=0, wc=0, wsg=1.07, injection=False, gsg=0.65, pb=0, rsb=0, sgsp=0.65, metric=False) -> float

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
     - Gas PVT object. If provided, extracts gas composition (sg, co2, h2s, n2, h2) and method selections for IPR calculations
   * - oil_pvt
     - OilPVT
     - Oil PVT object. If provided for oil wells, extracts api, sgsp, pb, rsb
   * - qg_mmscfd
     - float
     - Gas rate (MMscf/d, or sm3/d if metric=True). Gas wells only
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
     - Water specific gravity. Defaults to 1.07
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

Examples:

Gas well:

.. code-block:: python

    >>> from pyrestoolbox import nodal
    >>> c = nodal.Completion(tid=2.441, length=10000, tht=100, bht=200)
    >>> nodal.fbhp(thp=500, completion=c, vlpmethod='HB', well_type='gas', qg_mmscfd=5.0, gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    952.6868477414688

Oil well:

.. code-block:: python

    >>> c = nodal.Completion(tid=2.441, length=8000, tht=100, bht=180)
    >>> nodal.fbhp(thp=200, completion=c, vlpmethod='HB', well_type='oil', qt_stbpd=2000, gor=800, wc=0.3, gsg=0.65, pb=2500, rsb=500, sgsp=0.65, api=35)
    2256.2340921828286

Oil well using OilPVT object:

.. code-block:: python

    >>> from pyrestoolbox import oil
    >>> opvt = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    >>> nodal.fbhp(thp=200, completion=c, vlpmethod='HB', well_type='oil', oil_pvt=opvt, qt_stbpd=2000, gor=800, wc=0.3, gsg=0.65)
    2258.2220198140935

Deviated well using WellSegment:

.. code-block:: python

    >>> segs = [nodal.WellSegment(md=5000, id=2.441, deviation=0), nodal.WellSegment(md=5000, id=2.441, deviation=45)]
    >>> c_dev = nodal.Completion(segments=segs, tht=100, bht=200)
    >>> nodal.fbhp(thp=500, completion=c_dev, vlpmethod='HB', well_type='gas', qg_mmscfd=5.0, gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    923.092017723091

.. note::

   Not all VLP methods are equally suitable for deviated wells. For wells with deviation > 30 degrees, prefer **BB** or **WG** which have explicit inclination modelling. See `VLP Method Suitability for Deviated and Horizontal Wells`_ for detailed guidance.


pyrestoolbox.nodal.outflow_curve
================================

.. code-block:: python

    outflow_curve(thp, completion, vlpmethod='WG', well_type='gas', gas_pvt=None, oil_pvt=None, rates=None, n_rates=20, max_rate=None, cgr=0, qw_bwpd=0, oil_vis=1.0, api=45, pr=0, gor=0, wc=0, wsg=1.07, injection=False, gsg=0.65, pb=0, rsb=0, sgsp=0.65, metric=False) -> dict

Returns VLP outflow curve as a dictionary with keys ``'rates'`` and ``'bhp'``. Evaluates ``fbhp()`` at each rate point. Rates are MMscf/d for gas wells (sm3/d if metric=True), STB/d for oil wells (sm3/d if metric=True). BHP is in psia (barsa if metric=True).

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
   * - n_rates
     - int
     - Number of rate points if rates is None. Defaults to 20
   * - max_rate
     - float
     - Maximum rate for auto-generation. Defaults to 50 MMscf/d (gas) or 10000 STB/d (oil)
   * - metric
     - bool
     - If True, interpret inputs and return output in Eclipse METRIC units. Defaults to False
   * - (other)
     -
     - Same parameters as ``fbhp()``

Examples:

.. code-block:: python

    >>> c = nodal.Completion(tid=2.441, length=10000, tht=100, bht=200)
    >>> result = nodal.outflow_curve(thp=500, completion=c, vlpmethod='HB', well_type='gas', rates=[2.0, 5.0, 10.0, 15.0, 20.0], gsg=0.65)
    >>> result['rates']
    [2.0, 5.0, 10.0, 15.0, 20.0]
    >>> [round(b, 1) for b in result['bhp']]
    [676.8, 925.7, 1498.7, 2121.5, 2757.8]


pyrestoolbox.nodal.ipr_curve
=============================

.. code-block:: python

    ipr_curve(reservoir, well_type='gas', gas_pvt=None, oil_pvt=None, n_points=20, min_pwf=None, wc=0, wsg=1.07, bo=1.2, uo=1.0, gsg=0.65, metric=False) -> dict

Returns IPR (Inflow Performance Relationship) curve as a dictionary with keys ``'pwf'`` and ``'rate'``. Pressures are in psia (barsa if metric=True), rates in Mscf/d for gas (sm3/d if metric=True) or STB/d for oil (sm3/d if metric=True).

For gas wells, uses pseudopressure deliverability via ``gas.gas_rate_radial()``. Returns rates in Mscf/d.

For oil wells with OilPVT: uses Darcy above Pb, Vogel below Pb. Without OilPVT: simple Darcy. Returns rates in STB/d.

.. warning::

   **Gas rate unit mismatch between IPR and VLP:** Gas IPR rates from ``ipr_curve()`` are returned in **Mscf/d**, while VLP rates from ``outflow_curve()`` and the operating rate from ``operating_point()`` are in **MMscf/d**. When plotting IPR and VLP curves on the same axis, divide IPR rates by 1000. The ``operating_point()`` function handles this conversion internally.

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
     - Gas PVT object. Used for gas well composition (sg, co2, h2s, n2, h2)
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

Examples:

.. code-block:: python

    >>> r = nodal.Reservoir(pr=3000, degf=200, k=10, h=50, re=1500, rw=0.35, S=2, D=0.001)
    >>> ipr = nodal.ipr_curve(r, well_type='gas', gsg=0.65, n_points=5)
    >>> [round(p, 1) for p in ipr['pwf']]
    [14.7, 761.0, 1507.4, 2253.7, 3000.0]
    >>> [round(q, 1) for q in ipr['rate']]
    [13456.5, 12812.2, 10861.0, 7290.7, 0.0]


pyrestoolbox.nodal.operating_point
===================================

.. code-block:: python

    operating_point(thp, completion, reservoir, vlpmethod='WG', well_type='gas', gas_pvt=None, oil_pvt=None, cgr=0, qw_bwpd=0, oil_vis=1.0, api=45, gor=0, wc=0, wsg=1.07, gsg=0.65, pb=0, rsb=0, sgsp=0.65, bo=1.2, uo=1.0, n_points=25, metric=False) -> dict

Finds the operating point where VLP outflow curve intersects the IPR inflow curve via bisection. Returns a dictionary with keys:

- ``'rate'``: Operating rate (MMscf/d for gas / sm3/d if metric, STB/d for oil / sm3/d if metric)
- ``'bhp'``: Operating flowing BHP (psia, or barsa if metric=True)
- ``'vlp'``: VLP curve dict ``{'rates': [...], 'bhp': [...]}``
- ``'ipr'``: IPR curve dict ``{'pwf': [...], 'rate': [...]}``

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

Examples:

Gas well operating point:

.. code-block:: python

    >>> c = nodal.Completion(tid=2.441, length=10000, tht=100, bht=200)
    >>> r = nodal.Reservoir(pr=3000, degf=200, k=10, h=50, re=1500, rw=0.35, S=2, D=0.001)
    >>> result = nodal.operating_point(thp=500, completion=c, reservoir=r, vlpmethod='HB', well_type='gas', gsg=0.65)
    >>> round(result['rate'], 2)
    10.61
    >>> round(result['bhp'], 1)
    1573.7

Oil well operating point:

.. code-block:: python

    >>> from pyrestoolbox import oil
    >>> c = nodal.Completion(tid=2.441, length=8000, tht=100, bht=180)
    >>> r = nodal.Reservoir(pr=3000, degf=180, k=50, h=30, re=1000, rw=0.35)
    >>> opvt = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    >>> result = nodal.operating_point(thp=200, completion=c, reservoir=r, vlpmethod='HB', well_type='oil', oil_pvt=opvt, gor=800, wc=0.3, gsg=0.65)
    >>> round(result['rate'], 1)
    1412.6
    >>> round(result['bhp'], 1)
    2193.0
