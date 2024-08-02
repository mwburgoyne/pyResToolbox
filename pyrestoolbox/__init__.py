"""
pyrestoolbox
===================================

-----------------------------------------------
A collection of Reservoir Engineering Utilities
-----------------------------------------------

This set of functions focuses on those that the author uses often while crafting programming solutions. 
These are the scripts that are often copy/pasted from previous work - sometimes slightly modified - resulting in a trail of slightly different versions over the years. Some attempt has been made here to make this implementation flexible enough such that it can be relied on as-is going forward.

Note: Version 2.x now refactors functions into different modules, requiring seperate imports

Includes functions to perform simple calculations including;

- Inflow for oil and gas
- PVT Calculations for oil
- PVT calculation for gas
- Return critical parameters for typical components
- Creation of Black Oil Table information
- Creation of layered permeability distribution consistent with a Lorenze heterogeneity factor
- Extract problem cells information from Intesect (IX) print files
- Generation of AQUTAB include file influence functions for use in ECLIPSE
- Creation of Corey and LET relative permeability tables in Eclipse format
- Calculation of Methane and CO2 saturated brine properties


"""

submodules = [
    'brine',
    'classes',
    'constants',
    'gas',
    'layer',
    'library',
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