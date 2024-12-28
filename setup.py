#!/usr/bin/python3
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import os

ROOT = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(ROOT, 'README.md'), 'r', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='pyrestoolbox',
    include_package_data=True,
    version='2.1.3',  # Ideally should be same as your GitHub release tag version
    packages=find_packages(),
    description='pyResToolbox - A collection of Reservoir Engineering Utilities',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Mark W. Burgoyne',
    author_email='mark.w.burgoyne@gmail.com',
    url='https://github.com/mwburgoyne/pyResToolbox',
    keywords=['restoolbox', 'petroleum', 'reservoir'],
    classifiers=[],
    install_requires=[
        'requests',
        'numpy',
        'scipy',
        'pandas',
        'tabulate',
        'gwr_inversion',
        'mpmath',
        'openpyxl',
        'setuptools'
    ]
)