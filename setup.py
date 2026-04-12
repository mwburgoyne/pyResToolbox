#!/usr/bin/python3
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import os, re

ROOT = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(ROOT, 'README.rst'), 'r', encoding='utf-8') as f:
    long_description = f.read()

# Single-source version from pyproject.toml
with open(os.path.join(ROOT, 'pyproject.toml'), 'r', encoding='utf-8') as f:
    _version = re.search(r'^version\s*=\s*"([^"]+)"', f.read(), re.MULTILINE).group(1)

setup(
    name='pyrestoolbox',
    include_package_data=True,
    version=_version,
    packages=find_packages(exclude=['pyrestoolbox.tests', 'pyrestoolbox.tests.*']),
    description='pyResToolbox - A collection of Reservoir Engineering Utilities',
    license="GNU General Public License v3 or later (GPLv3+)",
    long_description=long_description,
    long_description_content_type='text/x-rst',
    author='Mark W. Burgoyne',
    author_email='mark.w.burgoyne@gmail.com',
    url='https://github.com/mwburgoyne/pyResToolbox',
    keywords=['restoolbox', 'petroleum', 'reservoir'],
    classifiers=[],
    install_requires=[
        'numpy',
        'scipy',
        'pandas',
        'tabulate',
        'ilt-inversion',
        'mpmath',
        'gmpy2',
        'python-flint',
        'openpyxl'
    ]
)