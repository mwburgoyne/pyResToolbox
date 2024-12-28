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