=========================================
pyResToolbox Documentation Index
=========================================

This directory ships inside the installed package. Locate it at runtime with
``pyrestoolbox.docs_dir()``. Each RST file documents one module, with a
function list, per-function signatures, input/return tables and worked
examples whose outputs are verified by the test suite.

Not sure which correlation method to use? Call
``pyrestoolbox.recommend.recommend_methods(...)`` first - it returns
structured recommendations with rationale.

Module documentation
====================

===================  ==========================================================================
File                 Content
===================  ==========================================================================
gas.rst              Gas PVT (Z-factor, viscosity, Bg, Cg, density, pseudopressure), flow
                     rates, pseudoskins, hydrates, GasPVT class
oil.rst              Oil PVT (Pb, Rs, Bo, density, viscosity, compressibility), flow rates,
                     harmonization, black oil tables, OilPVT class
brine.rst            Brine properties: CH4-saturated (brine_props), CO2-saturated
                     (CO2_Brine_Mixture), multicomponent gas-saturated (SoreideWhitson)
plyasunov.rst        Plyasunov infinite-dilution partial molar volume model (supporting
                     module used by brine density calculations)
nodal.rst            VLP/IPR nodal analysis: fbhp, fthp, outflow_curve, ipr_curve,
                     operating_point, Completion/Reservoir/WellSegment classes
dca.rst              Decline curve analysis: Arps, Duong, modified hyperbolic, two-segment
                     hyperbolic, fitting, forecasting, EUR
matbal.rst           Material balance: P/Z gas and Havlena-Odeh oil
simtools.rst         Simulation helpers: rel-perm tables, aquifer influence, VFP tables,
                     PVT tables, deck checking, IX PRT parsing
layer.rst            Lorenz heterogeneity and permeability layering
library.rst          Component critical property database
recommend.rst        Correlation method recommendation engine
sensitivity.rst      Parameter sweeps and tornado analysis
changelist.rst       Version history
===================  ==========================================================================

Example notebooks
=================

===========================  ====================================================
File                         Content
===========================  ====================================================
examples.ipynb               All documented examples plus BNS/BUR demonstrations
nodal_examples.ipynb         Multi-segment nodal framework examples
nodal_hydrate_demo.ipynb     VLP, VFPPROD, nodal solutions and hydrate analysis
===========================  ====================================================

Conventions
===========

- Import pattern: ``from pyrestoolbox import <module>`` then ``<module>.function(...)``
- Oilfield units by default (psia, deg F, ft, mD, cP); ``metric=True`` switches
  to Eclipse METRIC (barsa, deg C, m). Standard volumes at 60 deg F, 14.696 psia
- Gas and brine functions accept scalars, lists or numpy arrays; the oil module
  is scalar-only
- Correlation methods are selected with strings (or Enum objects); an invalid
  string raises ``ValueError`` listing the valid options
