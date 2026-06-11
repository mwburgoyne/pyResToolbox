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
- **dca** — Decline curve analysis: Arps (exponential/hyperbolic/harmonic), Duong,
  model fitting, forecasting, and EUR calculations.
- **matbal** — Material balance: P/Z gas material balance and Havlena-Odeh oil
  material balance for OGIP/OOIP estimation.
- **recommend** — Method recommendation engine for selecting appropriate correlations
  based on fluid composition, API gravity, and well deviation.
- **sensitivity** — Sensitivity analysis framework: parameter sweeps and tornado charts.
"""

submodules = [
    'brine',
    'classes',
    'constants',
    'dca',
    'gas',
    'layer',
    'library',
    'matbal',
    'nodal',
    'oil',
    'recommend',
    'sensitivity',
    'shared_fns',
    'simtools',
    'validate'
]

__all__ = submodules 

import importlib

def __dir__():
    return __all__


def _get_version():
    """Resolve the package version lazily.

    A repo checkout has pyproject.toml next to the package directory; read it
    so the version matches the imported code even when a different release is
    installed in site-packages. Installed packages have no adjacent
    pyproject.toml and use distribution metadata instead.
    """
    import os
    import re
    pyproject = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'pyproject.toml',
    )
    try:
        with open(pyproject, encoding='utf-8') as f:
            text = f.read()
        if re.search(r'^name\s*=\s*"pyrestoolbox"', text, re.M):
            match = re.search(r'^version\s*=\s*"([^"]+)"', text, re.M)
            if match:
                return match.group(1)
    except OSError:
        pass
    try:
        from importlib.metadata import version, PackageNotFoundError
        try:
            return version('pyrestoolbox')
        except PackageNotFoundError:
            pass
    except ImportError:  # Python < 3.8 has no importlib.metadata
        pass
    return 'unknown'


def __getattr__(name):
    if name in submodules:
        return importlib.import_module(f'pyrestoolbox.{name}')
    if name == '__version__':
        return _get_version()
    raise AttributeError(
        f"Module 'pyrestoolbox' has no attribute '{name}'"
    )