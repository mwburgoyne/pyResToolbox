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
- Numerous bugfixes, test suite expansion (313 tests), and code hardening across all modules.


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