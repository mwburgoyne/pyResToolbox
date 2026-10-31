================================
pyResToolbox Documentation Index
================================

This directory ships inside the installed package. Locate it at runtime with
``pyrestoolbox.docs_dir()``. Each RST file documents one module, with a
function list, per-function signatures, input/return tables and worked
examples whose outputs are verified by the test suite.

Not sure which correlation method to use? Call
``pyrestoolbox.recommend.recommend_methods(...)`` first - it returns
structured recommendations with rationale.

Find a function by task
=======================

Each entry is ``module.function``; the module's RST file documents it.

====================================================  ==========================================================
Task                                                  Function
====================================================  ==========================================================
Gas Z-factor, viscosity, density, Bg, Cg              gas.gas_z, gas.gas_ug, gas.gas_den, gas.gas_bg, gas.gas_cg
Gas critical properties, mixture SG                   gas.gas_tc_pc, gas.gas_sg
Gas enthalpy, Cp, Cv, Joule-Thomson                   gas.gas_thermal
Gas pseudopressure                                    gas.gas_dmp
Gas SG from gradient; p from p/z                      gas.gas_grad2sg, gas.gas_ponz2p
Condensate full-wellstream SG                         gas.gas_fws_sg
Water content of gas                                  gas.gas_water_content
Gas well rate (radial, linear, from pseudopressure)   gas.gas_rate_radial, gas.gas_rate_linear, gas.darcy_gas
Gas deliverability test (C,n or LIT, AOF)             gas.gas_deliverability
Non-Darcy and partial-penetration skin                gas.gas_hvf_beta, gas.gas_non_darcy_skin, gas.gas_partial_penetration_skin
Hydrate formation and inhibitor dosing                gas.gas_hydrate
Oil bubble point, Rs, Bo, density, viscosity, Co      oil.oil_pbub, oil.oil_rs, oil.oil_bo, oil.oil_deno, oil.oil_viso, oil.oil_co
Tune oil correlations to lab Pb, Rsb, viscosity       oil.oil_harmonize, oil.oil_harmonize_pb_rsb
Oil well rate (radial, linear)                        oil.oil_rate_radial, oil.oil_rate_linear
Gas gravities (evolved, stock tank, weighted)         oil.sg_evolved_gas, oil.sg_st_gas, oil.sgg_wt_avg
Black oil tables (PVTO/PVDG)                          oil.make_bot_og
Brine density, viscosity, Bw, gas solubility          brine.brine_props, brine.CO2_Brine_Mixture, brine.SoreideWhitson
Water PVT table (PVTW)                                brine.make_pvtw_table
Flowing BHP or THP (VLP)                              nodal.fbhp, nodal.fthp
VLP and IPR curves, operating point                   nodal.outflow_curve, nodal.ipr_curve, nodal.operating_point
Liquid loading (Turner, Coleman critical rate)        nodal.critical_gas_rate
Erosional velocity (API RP 14E)                       nodal.erosional_velocity
Wellhead choke (Gilbert, Ros, Baxendell, Achong)      nodal.choke_flow
Wellhead choke, critical and subcritical (Sachdeva)   nodal.choke_flow_sachdeva
Gas pipeline pressure drop or capacity                pipeline.gas_pipeline
VFP tables for simulation                             simtools.make_vfpprod, simtools.make_vfpinj
Decline fit, forecast, EUR                            dca.fit_decline, dca.fit_decline_cum, dca.forecast, dca.eur
Decline b-factor and EUR uncertainty                  dca.decline_b_range
GOR / WOR trend fit                                   dca.fit_ratio, dca.ratio_forecast
Gas material balance (p/z), oil Havlena-Odeh          matbal.gas_matbal, matbal.oil_matbal
Aquifer influx history (Carter-Tracy)                 matbal.aquifer_influx
Aquifer B and tD constants (field or metric)          matbal.aquifer_constants
Aquifer influence tables (Van Everdingen-Hurst)       simtools.influence_tables
Relative permeability tables and fits                 simtools.rel_perm_table, simtools.fit_rel_perm, simtools.fit_rel_perm_best
Rachford-Rice flash                                   simtools.rr_solver
Check a simulation deck; IX problem cells             simtools.zip_check_sim_deck, simtools.ix_extract_problem_cells
Waterflood connectivity, CRM fit and forecast         crm.crm_fit, crm.crm_oil_cut_fit
Layer permeabilities from a Lorenz coefficient        layer.lorenz_2_layers, layer.lorenzfromb, layer.lorenz_from_flow_fraction
Component critical properties                         library.prop
Which correlation to use                              recommend.recommend_methods
Sensitivity sweep, tornado chart                      sensitivity.sweep, sensitivity.tornado
====================================================  ==========================================================

Module documentation
====================

===================  ==========================================================================
File                 Content
===================  ==========================================================================
gas.rst              Gas PVT (Z-factor, viscosity, Bg, Cg, density, pseudopressure),
                     caloric properties (H, Cp, Cv, JT), flow rates, deliverability
                     tests, pseudoskins, hydrates, GasPVT class
oil.rst              Oil PVT (Pb, Rs, Bo, density, viscosity, compressibility), flow rates,
                     harmonization, black oil tables, OilPVT class
brine.rst            Brine properties: CH4-saturated (brine_props), CO2-saturated
                     (CO2_Brine_Mixture), multicomponent gas-saturated (SoreideWhitson)
plyasunov.rst        Plyasunov infinite-dilution partial molar volume model (supporting
                     module used by brine density calculations)
nodal.rst            VLP/IPR nodal analysis: fbhp, fthp, outflow_curve, ipr_curve,
                     operating_point, Completion/Reservoir/WellSegment classes;
                     liquid loading, erosional velocity, wellhead chokes
pipeline.rst         Gas pipeline flow: General Flow (Colebrook, AGA), Weymouth,
                     Panhandle A/B, IGT, elevation correction
crm.rst              Capacitance-resistance models (CRMP, CRMT) for waterfloods,
                     power-law oil cut
dca.rst              Decline curve analysis: Arps, Duong, modified hyperbolic, two-segment
                     hyperbolic, fitting, forecasting, EUR
matbal.rst           Material balance: P/Z gas and Havlena-Odeh oil, Carter-Tracy
                     aquifer influx
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
examples.ipynb               Worked examples across the modules plus BNS demonstrations
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
