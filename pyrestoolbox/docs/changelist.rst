Changelist in 3.0.2:

- **gas_hydrate()**: New function for gas hydrate formation prediction and thermodynamic inhibitor calculations. Returns a ``HydrateResult`` dataclass with hydrate formation temperature (HFT), hydrate formation pressure (HFP), subcooling, hydrate window assessment, inhibitor temperature depression, and required inhibitor concentration. Two HFT correlations: Motiee (1991) and Towler & Mokhatab (2005). Inhibitor depression via Østergaard et al. (2005). Supports five inhibitor types (MeOH, MEG, DEG, TEG, EtOH). Full Eclipse METRIC unit support.
- **gas_hydrate() water content**: Computes equilibrium water content of the gas. Uses the SoreideWhitson VLE model when gas composition (``co2``, ``h2s``, ``n2``, ``h2``) is provided, otherwise uses the Danesh correlation. Separate ``p_res``/``degf_res`` parameters allow specifying reservoir conditions for water content (where gas equilibrated with water) independently of the hydrate assessment point.
- **gas_hydrate() inhibitor capping**: ``required_inhibitor_wt_pct`` is now capped at the physical maximum for each inhibitor type (MEOH: 25%, MEG: 70%, DEG: 70%, TEG: 50%, ETOH: 30%). New ``max_inhibitor_wt_pct`` and ``inhibitor_underdosed`` fields indicate whether the required concentration exceeds the achievable maximum.
- **gas_hydrate() water balance**: Full water balance between reservoir and operating conditions. Reports vaporized water at both reservoir P,T (``water_vaporized_res``) and operating P,T (``water_vaporized_op``), condensed water (``water_condensed``), free water (``free_water``), and total liquid water (``total_liquid_water``). Only liquid water (condensed + free) needs inhibitor treatment.
- **gas_hydrate() injection rate**: New ``inhibitor_mass_rate`` (lb/MMscf | kg/sm3) and ``inhibitor_vol_rate`` (gal/MMscf | L/sm3) fields compute the required inhibitor injection rate from total liquid water (condensed + ``additional_water``) and required concentration.
- **hyd_method enum**: New enum for hydrate formation correlation selection (MOTIEE, TOWLER).
- **inhibitor enum**: New enum for thermodynamic hydrate inhibitor selection (MEOH, MEG, DEG, TEG, ETOH).
- New unit conversion constants in ``constants``: ``LB_TO_KG``, ``GAL_TO_LITER``, ``LB_PER_MMSCF_TO_KG_PER_SM3``, ``GAL_PER_MMSCF_TO_L_PER_SM3`` and their inverses.
- **oil_co() bugfix**: Fixed three bugs in oil compressibility numerical derivative. (1) Derivative stencil could cross the bubble point, mixing saturated and undersaturated physics and producing negative Co. (2) Near-Pb one-sided derivative used half the step span of the symmetric derivative without normalization, causing a 2x discontinuity. (3) The derivative formula omitted division by the pressure step size, making results proportional to dp rather than true 1/psi derivatives. All three bugs are now fixed; Co is positive at all pressures with smooth transitions across Pb.
- **_build_bot_tables() bugfix**: Fixed incorrect ``rsb`` parameter passed to ``oil_co()`` during black oil table generation. Was passing current pressure Rs (``rss[-1]``) instead of the actual bubble-point Rs, causing oil_co to compute derivatives around the wrong bubble point at sub-Pb pressures.
- Fixed shadowed duplicate ``test_doc_make_pvtw_table`` in test suite (renamed to ``test_doc_make_pvtw_table_keys``).
- Fixed ``MANIFEST.in`` to reference ``README.rst`` instead of non-existent ``README.md``.
- 416 validation tests (up from 375 in 3.0.1).


Changelist in 3.0.1:

- **OilPVT auto-harmonization**: ``rsb`` is now optional (default 0). When ``degf`` is provided (> 0), the constructor calls ``oil_harmonize()`` internally to resolve consistent Pb, Rsb, rsb_frac, and vis_frac from a single call. Accepts ``uo_target`` and ``p_uo`` for viscosity tuning. ``OilPVT.from_harmonize()`` remains as a deprecated thin wrapper.
- **oil_harmonize()**: New function (replaces ``oil_harmonize_pb_rsb``) that resolves consistent Pb, Rsb, rsb_frac, and ``vis_frac`` viscosity scaling factor from user inputs. Accepts ``uo_target`` and ``p_uo`` parameters to compute vis_frac = uo_target / uo_corr. The deprecated ``oil_harmonize_pb_rsb()`` wrapper remains for backward compatibility, returning the original 3-tuple.
- **OilPVT.vis_frac and OilPVT.rsb_frac**: New ``vis_frac`` and ``rsb_frac`` parameters on ``OilPVT`` constructor (both default 1.0). All ``viscosity()`` outputs are multiplied by vis_frac, and ``rs()`` applies rsb_frac scaling. Both factors flow through to VLP segment calculations (``fbhp()``, ``outflow_curve()``, ``operating_point()``) and ``make_bot_og()`` BOT generation.
- **oil_rate_radial / oil_rate_linear**: Accept ``oil_pvt`` and ``degf`` parameters. When an ``OilPVT`` object is provided, uo, bo, and pb are extracted automatically, and Vogel correction is enabled.
- **gas_rate_radial / gas_rate_linear**: Accept ``gas_pvt`` parameter. When a ``GasPVT`` object is provided, sg, composition, method choices, and pre-computed Tc/Pc are extracted automatically.
- **make_bot_og()**: New ``vis_frac`` parameter (default 1.0) that scales all oil viscosity values in generated black oil tables (PVDO, PVTO). Results dict now includes ``vis_frac`` key.
- **Completion.geometry_at_md()**: New method returning wellbore geometry (TVD, ID, deviation, roughness) at any measured depth along the completion. Supports both oilfield and metric unit systems.
- **Completion.profile()**: New method returning a pandas DataFrame of the wellbore profile at all segment boundaries, including crossover rows where geometry changes. Columns: MD, TVD, Deviation, ID, Roughness.
- Added explicit ``__all__`` exports and module docstrings to all modules for improved discoverability via ``dir()`` and ``help()``.
- ``OilPVT`` docs moved to ``oil.rst``; ``GasPVT`` docs moved to ``gas.rst``. Cross-references remain in ``nodal.rst``.
- 375 validation tests (up from 318 in 3.0.0).


Changelist in 3.0.0:

- **Nodal Analysis module**: Added VLP (Vertical Lift Performance), IPR (Inflow Performance Relationship), and operating point analysis with four multiphase flow correlations (Hagedorn-Brown, Woldesemayat-Ghajar, Gray, Beggs & Brill). Supports gas and oil wells with multi-segment deviated/horizontal wellbores via the new ``WellSegment`` and ``Completion`` classes. PVT convenience wrappers ``GasPVT`` and ``OilPVT`` for use with nodal and VFP functions.
- **VFP Table Generation**: New ``make_vfpprod()`` and ``make_vfpinj()`` functions in simtools for generating Eclipse VFPPROD and VFPINJ lift curve tables directly from wellbore geometry and fluid properties using the nodal VLP correlations.
- **Jerauld Relative Permeability Model**: Added the Jerauld (Arco) two-parameter kr model (``krfamily='JER'``) as a third option alongside Corey and LET in ``rel_perm_table()``.
- **Relative Permeability Curve Fitting**: New ``fit_rel_perm()`` function for fitting Corey, LET, or Jerauld models to measured kr data using least-squares optimization. ``fit_rel_perm_best()`` tries all three models and returns the best fit.
- **LET Physicality Check**: New ``is_let_physical()`` function to verify monotonicity and concavity of LET curves.
- **Simulation Workflow Consolidation**: ``make_bot_og()`` (black oil tables) and ``make_pvtw_table()`` (water PVT tables) now accessible via the ``simtools`` module, consolidating simulation-oriented functions in one place. Original locations (``oil.make_bot_og``, ``brine.make_pvtw_table``) remain as backward-compatible wrappers.
- **Multicomponent Gas-Saturated Brine**: New ``SoreideWhitson`` class for multi-gas brine properties using the Soreide-Whitson VLE framework with Garcia/Plyasunov density corrections and calibrated viscosity corrections. Supports CH4, C2-C4, CO2, H2S, N2, and H2.
- **IAPWS-IF97 Freshwater Density**: All brine models now use IAPWS-IF97 Region 1 for freshwater density base, improving accuracy across temperature and pressure ranges.
- **BNS Z-Factor Improvements**: Fugacity-based root selection for sub-critical conditions, vectorized Halley cubic solver, and tuned LBC viscosity model for hydrogen-containing mixtures.
- **Pseudopressure Performance**: Batch Gauss-Legendre quadrature replacing scipy integration for gas pseudopressure calculations, eliminating scipy dependency from the gas module.
- **Eclipse METRIC unit support**: All public PVT, flow rate, and simulation table functions now accept ``metric=False`` parameter. When ``metric=True``, inputs and outputs use Eclipse METRIC units (barsa, deg C, m, sm3/d, sm3/sm3, kg/m3, 1/bar). Applies to all gas, oil, brine, nodal, and simtools functions including ``GasPVT``, ``OilPVT``, ``WellSegment``, ``Completion``, ``Reservoir``, VFP table generation, and black oil table generation. Standard volumes always reference oilfield standard conditions (60 deg F, 14.696 psia).
- Numerous bugfixes, test suite expansion (318 tests), and code hardening across all modules.


Changelist in 2.2:

- Bugfixes.


Changelist in 2.1.3:

- Updated viscosity parameters for BUR method.


Changelist in 2.1.2:

- Fixed bug in implementation of Velarde, Blasingame & McCain Oil Rs calculation.


Changelist in 2.1.0:

- Fixed variable Typing issue that caused problems with Python 3.9 and older.
- Added reference to the Burgoyne ('BUR') methods for gas Z-Factor and critical property correlation


Changelist in 2.0.0:

- Modified the new Z-Factor method, 'BUR', now a tuned five component Peng Robinson method that is fast and stable and able to handle up to 100% of CO2, H2S, N2 or H2 as well as natural gas. Viscosities are calculated with a tuned LBC model.
- Refactored all code to split into modules for ease of future maintenance

Changelist in 1.4.4:

- Added in new Z-Factor method, 'BUR', which is a tuned five component Peng Robinson method that is fast and stable 

Changelist in 1.4.2:

- Corrected CO2 solubility calculations when two roots in CO2 liquid phase

Changelist in 1.4.1:

- Added calculation of Ezrokhi coefficients for brine density and viscosity with dissolved CO2

Changelist in 1.4.0:

- Introduced CO2 saturated brine calculations using Spycher & Pruess modified SRK EOS method
- Rectified an error introduced in Gas Z-Factor calculations due to errant indentation

Changelist in 1.3.9:

- Tweaks to speed DAK and Hall & Yarborough Z-Factor calculations

Changelist in 1.3.8:

- Fix bug in Hall & Yarborough Z-Factor algorithm

Changelist in 1.3.5:

- Fix bug in ECL deck zip/check recursion


Changelist in 1.3.4:

- Extend ECL deck zip/check function to handle IX formatted decks, and support zipping multiple decks at once.


Changelist in 1.3.2:

- Added robust Rachford Rice solver in Simulation Helpers
- Moved relative permeability functions and simulation helpers to seperate .simtools module


Changelist in v1.2.0:

- Added Component Critical Property Library


Changelist in v1.1.4:

- Attempting to fix reported issue with required dependencies not installing correctly


Changelist in v1.1:

- Fix API to SG calculation (141.4 vs 141.5)
- Added lower limit to first rho_po estimate for Oil Density with McCain method to avoid negative values with high Rs
- Added oil_sg and oil_api functions
- Modified HY Z-Factor solve algorithm to improve robustness
- Modified DAK Z-Factor solve algorithm to improve robustness
- Added Gas Z-Factor correlation from Wang, Ye & Wu (2021)
- Removed 'LIN' Z-Factor method due to significant errors above 12,000 psi. Use WYW method instead if speed needed.