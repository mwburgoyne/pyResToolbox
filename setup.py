# -*- coding: utf-8 -*-

# Learn more: https://github.com/kennethreitz/setup.py

from setuptools import setup, find_packages


with open('README.rst') as f:
    readme = f.read()

with open('license.txt') as f:
    license = f.read()

setup(
    name='pyrestoolbox',
    version='0.0.1',
    description='pyResToolbox - A collection of Reservoir Engineering Utilities',
    long_description=readme,
    author='Mark Burgoyne',
    author_email='mark.w.burgoyne@gmail.com',
    url='https://github.com/vinomarkus/pyResToolbox',
    license=license,
    packages=find_packages(exclude=('tests', 'docs'))
)