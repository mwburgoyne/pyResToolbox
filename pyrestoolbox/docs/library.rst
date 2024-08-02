===================================
Component Critical Properties Library
===================================

Returns component Critical properties for common species

Function List
=============

.. list-table:: Component Library Functions
   :widths: 15 40
   :header-rows: 1

   * - Task
     - Function
   * - Return Critical Property 
     - `pyrestoolbox.library.prop`_  
   * - List of all components in Library
     - `pyrestoolbox.library.components`_
   * - List of all long-form names of components in Library
     - `pyrestoolbox.library.names`_
   * - Returns a list of all critical properties available for components in the Library
     - `pyrestoolbox.library.property_list`_
   * - Returns a list of all valid EOS models in the Library 
     - `pyrestoolbox.library.models`_
   * - Returns a DataFrame of the Library 
     - `pyrestoolbox.library.df`_
     
     

pyrestoolbox.library.prop
======================

.. code-block:: python

    pyrestoolbox.library.prop(comp, prop, model = 'PR79') -> float:

Returns specified critical property of a component    

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - comp
     - string
     - The component for which critical properties are to be returned (use rtb.comp_library.components to get available list)
   * - prop
     - str
     - The critical property to be returned. Selections are 'MW', 'Tc_R', 'Pc_psia', 'Visc_Zc', 'Pchor', 'Vc_cuft_per_lbmol', 'Acentric', 'VTran', 'Tb_F', 'SpGr'.
   * - model
     - str
     - The EOS model. Selections are 'PR79', 'PR77', 'SRK', 'RK', with default 'PR79'. Used for EOS specific critical properties for Acentric, VTran, Tb_F and SpGr.


Examples:

.. code-block:: python

    >>> from pyrestoolbox import library
    >>> library.prop(comp = 'CH4', prop = 'Pc_psia')
    667.029

.. code-block:: python

    >>> library.prop(comp = 'C3', prop = 'VTran'), rtb.comp_library.prop(comp = 'C3', prop = 'VTran', model='SRK')
    (-0.06381, 0.09075)


pyrestoolbox.library.components
======================

.. code-block:: python

    pyrestoolbox.library.components -> list:

Returns a list of all available components in the library    


Example:

.. code-block:: python

    >>> print(library.components)
    ['HE', 'HELIUM', 'NE', 'NEON', 'AR', 'ARGON', 'KR', 'KRYPTON', 'XE', 'XEON', 'RN', 'RADON', 'H2', 'HYDROGEN', 'N2', 'NITROGEN', 'CO', 'O2', 'OXYGEN', 'NO', 'N2O', 'CO2', 'H2S', 'NH3', 'AMMONIA', 'SO2', 'NO2', 'N2O4', 'H2O', 'WATER', 'C1', 'CH4', 'METHANE', 'C2', 'C2H6', 'ETHANE', 'C3', 'C3H8', 'PROPANE', 'C-C3', 'CYCLO-C3', 'C-PROPANE', 'CYCLOPROP', 'I-C4', 'ISO-C4', 'I-BUTANE', 'ISOBUTANE', 'N-C4', 'N-BUTANE', 'BUTANE', 'NEO-C5', 'NEOPENTAN', 'C-C4', 'CYCLO-C4', 'C-BUTANE', 'CYCLOBUTA', 'I-C5', 'ISO-C5', 'I-PENTANE', 'ISOPENTAN', 'N-C5', 'N-PENTANE', 'PENTANE', 'C-C5', 'CYCLO-C5', 'C-PENTANE', 'CYCLOPENT', '22DM-C4', '22DM-BUTA', '23DM-C4', '23DM-BUTA', '2M-C5', '2M-PENTAN', '3M-C5', '3M-PENTAN', 'N-C6', 'N-HEXANE', 'HEXANE', 'MC-C5', 'MC-PENTAN', '22DM-C5', '22DM-PENT', 'BENZENE', '24DM-C5', '24DM-PENT', 'C-C6', 'CYCLO-C6', 'C-HEXANE', 'CYCLOHEXA', '223TM-C4', '223TM-BUT', '33DM-C5', '33DM-PENT', '23DM-C5', '23DM-PENT', '2M-C6', '2M-HEXANE', '3M-C6', '3M-HEXANE', '3E-C5', '3E-PENTAN', 'N-C7', 'N-HEPTANE', 'HEPTANE', 'MC-C6', 'MC-HEXANE', 'EC-C5', 'EC-PENTAN', 'TOLUENE', 'C-C7', 'CYCLO-C7', 'C-HEPTANE', 'CYCLOHEPT', 'N-C8', 'N-OCTANE', 'OCTANE', 'E-BENZENE', 'P-XYLENE', 'M-XYLENE', 'O-XYLENE', 'N-C9', 'N-NONANE', 'NONANE', 'C-C8', 'CYCLO-C8', 'C-OCTANE', 'CYCLOOCTA', 'CUMENE', 'I-C3-BENZ', '1ME-BENZE', 'P-BENZENE', '1E4M-BENZ', '135TM-BEN', '124TM-BEN', 'N-C10', 'N-DECANE', 'DECANE', '123TM-BEN', 'N-C11', 'N-UNDECAN', 'UNDECANE', 'N-C12', 'N-DODECAN', 'DODECANE', 'NAPTHALEN', 'N-C13', 'N-TRIDECA', 'TRIDECANE', '2M-NAPTHA', '1M-NAPTHA', 'N-C14', 'N-TETRADE', 'TETRADECA', 'DPH-C1', 'DPH-METHA', 'N-C15', 'N-PENTADE', 'PENTADECA', 'N-C16', 'N-HEXADEC', 'HEXADECAN', 'N-C17', 'N-HEPTADE', 'HEPTADECA', 'N-C18', 'N-OCTADEC', 'OCTADECAN', 'N-C19', 'N-NONADEC', 'NONADECAN', '12DPH-BEN', 'PHENANTHR', 'ANTHRACEN', 'N-C20', 'N-EICOSAN', 'EICOSANE', 'N-C21', 'N-HENEICO', 'HENEICOSA', '13DPH-BEN', 'N-C22', 'N-DOCOSAN', 'DOCOSANE', '14DPH-BEN', 'N-C23', 'N-TRICOSA', 'TRICOSANE', 'N-C24', 'N-TETRACO', 'TETRACOSA', 'N-C25', 'N-C26', 'N-C27', 'N-C28', 'N-C29', 'N-C30', 'N-C31', 'N-C32', 'N-C33', 'N-C34']



pyrestoolbox.library.names
======================

.. code-block:: python

    pyrestoolbox.library.names -> list:

Returns a list of long-form names of all components available in the Library   


Example:

.. code-block:: python

    >>> print(library.names)
    ['Helium', 'Helium', 'Neon', 'Neon', 'Argon', 'Argon', 'Krypton', 'Krypton', 'Xenon', 'Xenon', 'Radon', 'Radon', 'Hydrogen', 'Hydrogen', 'Nitrogen', 'Nitrogen', 'Carbon Monoxide', 'Oxygen', 'Oxygen', 'Nitric Oxide', 'Nitrous Oxide', 'Carbon Dioxide', 'Hydrogen Sulfide', 'Ammonia', 'Ammonia', 'Sulfur Dioxide', 'Nitrogen Dioxide', 'Nitrogen Tetroxide', 'Water', 'Water', 'Methane', 'Methane', 'Methane', 'Ethane', 'Ethane', 'Ethane', 'Propane', 'Propane', 'Propane', 'Cyclopropane', 'Cyclopropane', 'Cyclopropane', 'Cyclopropane', 'Isobutane', 'Isobutane', 'Isobutane', 'Isobutane', 'Butane', 'Butane', 'Butane', 'Neopentane', 'Neopentane', 'Cyclobutane', 'Cyclobutane', 'Cyclobutane', 'Cyclobutane', 'Isopentane', 'Isopentane', 'Isopentane', 'Isopentane', 'Pentane', 'Pentane', 'Pentane', 'Cyclopentane', 'Cyclopentane', 'Cyclopentane', 'Cyclopentane', '2,2-Dimethylbutane', '2,2-Dimethylbutane', '2,3-Dimethylbutane', '2,3-Dimethylbutane', '2-Methylpentane', '2-Methylpentane', '3-Methylpentane', '3-Methylpentane', 'Hexane', 'Hexane', 'Hexane', 'Methylcyclopentane', 'Methylcyclopentane', '2,2-Dimethylpentane', '2,2-Dimethylpentane', 'Benzene', '2,4-Dimethylpentane', '2,4-Dimethylpentane', 'Cyclohexane', 'Cyclohexane', 'Cyclohexane', 'Cyclohexane', '2,2,3-Trimethylbutane', '2,2,3-Trimethylbutane', '3,3-Dimethylpentane', '3,3-Dimethylpentane', '2,3-Dimethylpentane', '2,3-Dimethylpentane', '2-Methylhexane', '2-Methylhexane', '3-Methylhexane', '3-Methylhexane', '3-Ethylpentane', '3-Ethylpentane', 'Heptane', 'Heptane', 'Heptane', 'Methylcyclohexane', 'Methylcyclohexane', 'Ethylcyclopentane', 'Ethylcyclopentane', 'Toluene', 'Cycloheptane', 'Cycloheptane', 'Cycloheptane', 'Cycloheptane', 'Octane', 'Octane', 'Octane', 'Ethylbenzene', 'p-Xylene', 'm-Xylene', 'o-Xylene', 'Nonane', 'Nonane', 'Nonane', 'Cyclooctane', 'Cyclooctane', 'Cyclooctane', 'Cyclooctane', 'Cumene', 'Cumene', 'Cumene', 'Propylbenzene', '1-Ethyl-4-methylbenzene', '1,3,5-Trimethylbenzene', '1,2,4-Trimethylbenzene', 'Decane', 'Decane', 'Decane', '1,2,3-Trimethylbenzene', 'Undecane', 'Undecane', 'Undecane', 'Dodecane', 'Dodecane', 'Dodecane', 'Napthalene', 'Tridecane', 'Tridecane', 'Tridecane', '2-Methylnapthalene', '1-Methylnapthalene', 'Tetradecane', 'Tetradecane', 'Tetradecane', 'Diphenylmethane', 'Diphenylmethane', 'Pentadecane', 'Pentadecane', 'Pentadecane', 'Hexadecane', 'Hexadecane', 'Hexadecane', 'Heptadecane', 'Heptadecane', 'Heptadecane', 'Octadecane', 'Octadecane', 'Octadecane', 'Nonadecane', 'Nonadecane', 'Nonadecane', '1,2-Diphenylbenzene', 'Phenanthrene', 'Anthracene', 'Eicosane', 'Eicosane', 'Eicosane', 'Heneicosane', 'Heneicosane', 'Heneicosane', '1,3-Diphenylbenzene', 'Docosane', 'Docosane', 'Docosane', '1,4-Diphenylbenzene', 'Tricosane', 'Tricosane', 'Tricosane', 'Tetracosane', 'Tetracosane', 'Tetracosane', 'N-C25', 'N-C26', 'N-C27', 'N-C28', 'N-C29', 'N-C30', 'N-C31', 'N-C32', 'N-C33', 'N-C34']


pyrestoolbox.library.property_list
======================

.. code-block:: python

    pyrestoolbox.library.property_list -> list:

Returns a list of all critical properties available for components in the Library   


Example:

.. code-block:: python

    >>> print(library.property_list)
    ['Name', 'MW', 'Tc_R', 'Pc_psia', 'Visc_Zc', 'Pchor', 'Vc_cuft_per_lbmol', 'Acentric', 'VTran', 'Tb_F', 'SpGr']


pyrestoolbox.library.models
======================

.. code-block:: python

    pyrestoolbox.library.models -> list:

Returns a list of all valid EOS models in the Library   


Example:

.. code-block:: python

    >>> print(library.models)
    ['PR79', 'PR77', 'SRK', 'RK']


pyrestoolbox.library.df
======================

.. code-block:: python

    pyrestoolbox.library.df -> pandas.DataFrame:

Returns a dataframe of the Library data   


Example:

.. code-block:: python

    >>> library.df

.. image:: https://github.com/mwburgoyne/pyResToolbox/blob/main/docs/img/properties_df.png
    :alt: DataFrame of Component Library data