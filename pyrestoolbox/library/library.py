#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
    pyResToolbox - A collection of Reservoir Engineering Utilities
              Copyright (C) 2022, Mark Burgoyne

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    The GNU General Public License can be found in the LICENSE directory,
    and at  <https://www.gnu.org/licenses/>.

          Contact author at mark.w.burgoyne@gmail.com
"""

import sys
from collections import Counter
import glob
from enum import Enum
import pkg_resources

import numpy as np
import numpy.typing as npt
import pandas as pd
from tabulate import tabulate


class component_library:
    def __init__(self, model='PR79'):
        path = 'component_library.xlsx'
        filepath = pkg_resources.resource_filename(__name__, path)
        self.df = pd.read_excel(filepath, engine="openpyxl")
        self.model = model
        self.all_cols = ['Name', 'MW', 'Tc_R', 'Pc_psia',
                         'Zc', 'Pchor', 'Vc_cuft_per_lbmol']
        self.model_cols = ['Acentric', 'VTran', 'Tb_F', 'SpGr']
        self.all_dics = {}
        self.model_dics = {}
        self.components = self.df['Component'].tolist()
        self.names = self.df['Name'].tolist()
        self.property_list = self.all_cols + self.model_cols
        # Create dictionaries for all the model agnostic properties
        for col in self.all_cols:
            self.all_dics[col.upper()] = dict(
                zip(self.df['Component'], self.df[col]))
        # And then for all the model specific properties
        self.models = ['PR79', 'PR77', 'SRK', 'RK']
        for model in self.models:
            model_dic = {}
            for col in self.model_cols:
                model_dic[col.upper()] = dict(
                    zip(self.df['Component'], self.df[model+'-'+col]))
            self.model_dics[model] = model_dic

    def prop(self, comp, prop, model='PR79'):
        comp = comp.upper()
        if comp not in self.components:
            return 'Component not in Library. Choose from Name, MW, Tc_R, Pc_psia, Pchor, Zc, Vc_cuft_per_lbmol, Acentric, VTran, Tb_F, SpGr'
        prop = prop.upper()
        props = [x.upper() for x in self.property_list]
        if model.upper() not in self.models:
            return 'Incorrect Model Name'
        dic = self.model_dics[model.upper()]
        if prop == 'ALL':
            return [self.all_dics[p.upper()][comp] for p in self.all_cols]+[dic[p.upper()][comp] for p in self.model_cols]
        props = [x.upper() for x in self.property_list]
        if prop not in props:
            return 'Property not in Library'
        if prop in [x.upper() for x in self.all_cols]:
            return self.all_dics[prop][comp]
        if prop in [x.upper() for x in self.model_cols]:
            return dic[prop][comp]
        return 'Component or Property not in library'

comp_library = component_library()