===================================
Inflow
===================================

Calculation Methods and Class Objects
=====================================
pyResToolBox uses class objects to track calculation options through the functions. Class objects can be set via strings or explicitly via object options

.. list-table:: Method Variables & Class Objects
   :widths: 10 15 40
   :header-rows: 1

   * - Class Variable
     - Class Object 
     - Class Description & Options
   * - zmethod
     - z_method
     - Method for calculating gas Z-Factor. Defaults to 'DAK'. 
       Options are:
        + 'DAK': Dranchuk & Abou-Kassem (1975) using from Equations 2.7-2.8 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al. - Slowest, Most Accurate
        + 'HY': Hall & Yarborough (1973) - Second Fastest
        + 'WYW': Wang, Ye & Wu (2021) - Fastest, Least Accurate
   * - cmethod
     - c_method
     - Method for calculating gas critical properties. Defaults to 'PMC' 
       Options are:
        + 'SUT': Sutton with Wichert & Aziz non-hydrocarbon corrections
        + 'PMC': Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.


Users can specify which calculation method to use either by passing an option string, or a class object to any given function. The implementation of class objects should make it easier to program in an IDE that supports type hinting


Function List
=============

.. list-table:: Inflow Functions
   :widths: 15 40
   :header-rows: 1

   * - Task
     - Function
   * - Gas Flow Rate Radial
     - `pyrestoolbox.gas_rate_radial`_  
   * - Gas Flow Rate Linear 
     - `pyrestoolbox.gas_rate_linear`_
   * - Oil Flow Rate Radial
     - `pyrestoolbox.oil_rate_radial`_
   * - Oil Flow Rate Linear
     - `pyrestoolbox.oil_rate_linear`_
  

pyrestoolbox.gas_rate_radial
======================

.. code-block:: python

    gas_rate_radial(k, h, pr, pwf, r_w, r_ext, degf, zmethod='DAK, cmethod='PMC', S = 0, D = 0, sg = 0.75, n2 = 0, co2 = 0, h2s = 0, tc  = 0, pc = 0) -> float or np.array

Returns gas rate (mscf/day) for radial flow using Darcy pseudo steady state equation & gas pseudopressure. 
Arrays can be used for any one of k, h, pr or pwf, returning corresponding 1-D array of rates. Using more than one input array – while not prohibited - will not return expected results 

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - k
     - float, list or np.array
     - Effective permeability to gas flow (mD)
   * - h
     - float, list or np.array
     - Net height for flow (ft).
   * - pr
     - float, list or np.array
     - Reservoir pressure (psia)
   * - pwf
     - float, list or np.array
     - BHFP (psia).
   * - r_w
     - float
     - Wellbore Radius (ft).
   * - r_ext
     - float
     - External Reservoir Radius (ft).
   * - degf
     - float
     - Reservoir Temperature (deg F). 
   * - zmethod
     - string or z_method
     - Method for calculating gas Z-factor. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified 
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined  
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - S
     - float
     - Skin. Defaults to zero if undefined
   * - D
     - float
     - Non Darcy Skin Factor (day/mscf). Defaults to zero if undefined
   * - sg
     - float
     - Gas SG relative to air, Defaults to 0.75 if undefined
     
Examples:

.. code-block:: python

    >>> from pyrestoolbox import pyrestoolbox as rtb
    >>> rtb.gas_rate_radial(k=5, h=50, pr=2000, pwf=750, r_w=0.3, r_ext=1500, degf=180, sg = 0.75, D = 0.01, S=5)
    10269.669190157822
    
    >>> rtb.gas_rate_radial(k=1, h=50, pr=[2000,1000], pwf=750, r_w=0.3, r_ext=1500, degf=180, sg = 0.75, D = 0.01, S=5)
    array([4273.15956785, 1177.38697977])
    

pyrestoolbox.gas_rate_linear
======================

.. code-block:: python

    gas_rate_linear(k, pr, pwf, area, length, degf, zmethod='DAK, cmethod='PMC', sg = 0.75, n2 = 0, co2 = 0, h2s = 0, tc  = 0, pc = 0) -> float or np.array

Returns gas rate (mscf/day) for linear flow using Darcy steady state equation & gas pseudopressure. 
Arrays can be used for any one of k, pr, pwf or area, returning corresponding 1-D array of rates. Using more than one input array – while not prohibited - will not return expected results 


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - k
     - float, list or np.array
     - Effective permeability to gas flow (mD)
   * - pr
     - float, list or np.array
     - Reservoir pressure (psia)
   * - pwf
     - float, list or np.array
     - BHFP (psia).
   * - area
     - float
     - Net cross-sectional area perpendicular to direction of flow (ft2)
   * - length
     - float
     - Linear distance of fluid flow (ft).
   * - degf
     - float
     - Reservoir Temperature (deg F). 
   * - zmethod
     - string or z_method
     - Method for calculating gas Z-factor. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified 
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined  
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - sg
     - float
     - Gas SG relative to air, Defaults to 0.75 if undefined
     
Examples:

.. code-block:: python

    >>> rtb.gas_rate_linear(k=0.1, area=50, length=200, pr=2000, pwf=250, degf=180, sg = 0.8)
    21.87803915816601
    
    >>> rtb.gas_rate_linear(k=0.1, area=50, length=200, pr=[2000, 1000, 500], pwf=250, degf=180, sg = 0.8)
    array([21.87803916,  4.89593662,  0.94342881])
    
pyrestoolbox.oil_rate_radial
======================

.. code-block:: python

    oil_rate_radial(k, h, pr, pwf, r_w, r_ext, uo, bo, S = 0, vogel = False, pb = 0) -> float or np.array

Returns liquid rate (stb/day) for radial flow using Darcy pseudo steady state equation with optional Vogel correction.
Arrays can be used for any one of k, h, pr or pwf, returning corresponding 1-D array of rates. Using more than one input array – while not prohibited - will not return expected results 

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - k
     - float, list or np.array
     - Effective permeability to gas flow (mD)
   * - h
     - float, list or np.array
     - Net height for flow (ft).
   * - pr
     - float, list or np.array
     - Reservoir pressure (psia)
   * - pwf
     - float, list or np.array
     - BHFP (psia).
   * - r_w
     - float
     - Wellbore Radius (ft).
   * - r_ext
     - float
     - External Reservoir Radius (ft).
   * - uo
     - float
     - Liquid viscosity (cP). 
   * - bo
     - float
     - Liquid formation volume factor (rb/stb)
   * - S
     - float
     - Skin. Defaults to zero if undefined
   * - vogel
     - bool
     - Boolean flag indicating whether to use vogel Pb correction. Defaults to False
   * - pb
     - float
     - Bubble point pressure. Used only when Vogel correction is invoked
     
Examples:

.. code-block:: python

    >>> rtb.oil_rate_radial(k=20, h=20, pr=1500, pwf=250, r_w=0.3, r_ext=1500, uo=0.8, bo=1.4, vogel=True, pb=1800)
    213.8147848023242
    
    >>> rtb.oil_rate_radial(k=20, h=20, pr=[1500, 2000], pwf=250, r_w=0.3, r_ext=1500, uo=0.8, bo=1.4, vogel=True, pb=1800)
    array([213.8147848 , 376.58731835])
    
pyrestoolbox.oil_rate_linear
======================

.. code-block:: python

    oil_rate_linear(k, pr, pwf, area, length, uo, bo, vogel = False, pb = 0) -> float or np.array

Returns liquid rate (stb/day) for linear flow using Darcy steady state equation with optional Vogel correction.
Arrays can be used for any one of k, pr, pwf or area, returning corresponding 1-D array of rates. Using more than one input array – while not prohibited - will not return expected results 

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - k
     - float, list or np.array
     - Effective permeability to gas flow (mD)
   * - pr
     - float, list or np.array
     - Reservoir pressure (psia)
   * - pwf
     - float, list or np.array
     - BHFP (psia).
   * - area
     - float, list or np.array
     - Net cross-sectional area perpendicular to direction of flow (ft2)
   * - length
     - float
     - Linear distance of fluid flow (ft)
   * - bo
     - float
     - Liquid formation volume factor (rb/stb)
   * - vogel
     - bool
     - Boolean flag indicating whether to use vogel Pb correction. Defaults to False
   * - pb
     - float
     - Bubble point pressure. Used only when Vogel correction is invoked
     
Examples:

.. code-block:: python

    >>> rtb.oil_rate_linear(k=0.1, area=15000, pr=3000, pwf=500, length=500, uo=0.4, bo=1.5)
    14.08521246363274
    
    >>> rtb.oil_rate_linear(k=[0.1, 1, 5, 10], area=15000, pr=3000, pwf=500, length=500, uo=0.4, bo=1.5)
    array([  14.08521246,  140.85212464,  704.26062318, 1408.52124636])