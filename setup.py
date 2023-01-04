#!/usr/bin/python3
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
#import markdown
import importlib
import collections
import os
import pkg_resources

ROOT = os.path.abspath(os.path.dirname(__file__))
f = open('C:\\Users\\vinom\OneDrive - Santos\\Work in progress\\ML Learning tidbits\\0-ResEng\\pyResToolbox\\source_dir\\README.md', 'r').read()
long_description= f #markdown.markdown( f.read() )

setup(
    name = 'pyrestoolbox',
    packages = find_packages(),
    include_package_data=True,
    version = '1.3.8',  # Ideally should be same as your GitHub release tag varsion
    description = 'pyResToolbox - A collection of Reservoir Engineering Utilities',
    long_description= long_description,
    long_description_content_type = 'text/markdown',
    author = 'Mark W. Burgoyne',
    author_email = 'mark.w.burgoyne@gmail.com',
    url = 'https://github.com/mwburgoyne/pyResToolbox',
    keywords = ['restoolbox', 'petroleum', 'reservoir'],
    classifiers = [],
    install_requires=[
        'requests',
        'numpy',
        'scipy',
        'pandas',
        'tabulate',
        'gwr_inversion',
        'mpmath',
        'openpyxl'
    ]
)