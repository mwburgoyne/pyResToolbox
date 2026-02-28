"""
pyrestoolbox
===================================

-----------------------------------------------
A collection of Reservoir Engineering Utilities
-----------------------------------------------

This set of functions focuses on those that the author uses often while crafting programming solutions.
These are the scripts that are often copy/pasted from previous work - sometimes slightly modified - resulting
in a trail of slightly different versions over the years. Some attempt has been made here to make this
implementation flexible enough such that it can be relied on as-is going forward.

Note: Version 2.x refactors functions into different modules, requiring separate imports.

Usage: ``from pyrestoolbox import <module>`` then ``<module>.function(...)``.
All inputs use oilfield units (psia, degF, ft, mD, cP) unless otherwise noted.

Modules
-------
- **gas** — Gas PVT (Z-factor, viscosity, Bg, Cg, density, pseudopressure) and flow rates.
  Four Z-factor methods: DAK, HY, WYW, BNS. For high CO2/H2S/N2 or any H2, use BNS.
  Includes GasPVT convenience class.
- **oil** — Oil PVT (Pb, Rs, Bo, density, viscosity, compressibility) and flow rates.
  Black oil table generation (PVTO/PVDO/PVDG). Includes OilPVT convenience class.
- **brine** — Brine properties with three models: methane-saturated (brine_props),
  CO2-saturated (CO2_Brine_Mixture), multicomponent gas-saturated (SoreideWhitson).
  Supports C1-C4, CO2, H2S, N2, H2 in fresh or saline water.
- **nodal** — VLP/IPR nodal analysis with four multiphase correlations (HB, WG, GRAY, BB).
  Multi-segment deviated/horizontal well support. Operating point solver.
- **layer** — Lorenz heterogeneity framework and permeability distribution generation.
- **library** — Component critical property database with EOS-model-specific parameters.
- **simtools** — Simulation helpers: IX PRT parsing, rel-perm tables (Corey/LET),
  Van Everdingen-Hurst aquifer influence, Rachford-Rice flash, deck file checking.
"""

submodules = [
    'brine',
    'classes',
    'constants',
    'gas',
    'layer',
    'library',
    'nodal',
    'oil',
    'shared_fns',
    'simtools',
    'validate'
]

__all__ = submodules 

import importlib

def __dir__():
    return __all__


def __getattr__(name):
    if name in submodules:
        return importlib.import_module(f'pyrestoolbox.{name}')
    else:
        try:
            return globals()[name]
        except KeyError:
            raise AttributeError(
                f"Module 'pyrestoolbox' has no attribute '{name}'"
            )