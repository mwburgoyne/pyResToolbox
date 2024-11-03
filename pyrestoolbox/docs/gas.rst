===================================
Gas PVT & Flow
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
        + 'BUR': Tuned five component Peng Robinson EOS model (Unpublished, created by M. Burgoyne 2024) - Fast, reliable, and able to handle wide range of mixture types including CO2, H2S, N2 and H2 at concentrations up to pure inerts. More information about the method can be found `here <https://github.com/mwburgoyne/5_Component_PengRobinson_Z-Factor>`_  
   * - cmethod
     - c_method
     - Method for calculating gas critical properties. Defaults to 'PMC' 
       Options are:
        + 'SUT': Sutton with Wichert & Aziz non-hydrocarbon corrections
        + 'PMC': Piper, McCain & Corredor (1999) correlation, using equations 2.4 - 2.6 from 'Petroleum Reservoir Fluid Property Correlations' by W. McCain et al.
        + 'BUR': Correlation tuned to return the critical properties of the pure hydrocarbon component of a mixture, required for the 'BUR' (tuned Peng Robinson) Z-Factor method. More information about the method can be found `here <https://github.com/mwburgoyne/5_Component_PengRobinson_Z-Factor>`_ 


Users can specify which calculation method to use either by passing an option string, or a class object to any given function. The implementation of class objects should make it easier to program in an IDE that supports type hinting

Examples:

Calculating gas Z-Factor of pure methane using DAK and PMC for critical properties

.. code-block:: python

    >>> from pyrestoolbox import gas
    >>> gas.gas_z(p=2350, sg=0.68, degf = 180, zmethod='DAK', cmethod='PMC')
    0.8785390750839772
    
Calculating gas Z-Factor of pure CO2

.. code-block:: python

    >>> gas.gas_z(p=2350, sg=0.68, degf = 180, co2 = 1.0, zmethod='BUR', cmethod='BUR')
    0.5258204505475871
    
Calculating gas sg, and then gas Z-Factor of a mixture of 5% CO2, 10% H2S, 0% N2 and 20% H2 (remainder natural gas with MW = 19).

.. code-block:: python

    >>> gsg = gas.gas_sg(hc_mw = 19.0, co2 = 0.05, h2s = 0.10, n2 = 0, h2 = 0.20)
    >>> gas.gas_z(p=2350, sg=gsg, degf = 180, co2 = 0.05, h2s = 0.10, n2 = 0, h2 = 0.20, zmethod='BUR', cmethod='BUR')
    0.9182754436816966


Function List
=============

.. list-table:: Gas Functions
   :widths: 15 40
   :header-rows: 1

   * - Task
     - Function
   * - Gas Tc & Pc Calculation
     - `pyrestoolbox.gas.gas_tc_pc`_  
   * - Gas Z-Factor Calculation
     - `pyrestoolbox.gas.gas_z`_
   * - Gas Viscosity or Viscosity * Z product
     - `pyrestoolbox.gas.gas_ug`_
   * - Gas Compressibility
     - `pyrestoolbox.gas.gas_cg`_
   * - Gas Formation Volume Factor
     - `pyrestoolbox.gas.gas_bg`_  
   * - Gas Density
     - `pyrestoolbox.gas.gas_den`_  
   * - Gas Water of Condensation
     - `pyrestoolbox.gas.gas_water_content`_
   * - Convert P/Z to P
     - `pyrestoolbox.gas.gas_ponz2p`_
   * - Convert Gas Gradient to SG
     - `pyrestoolbox.gas.gas_grad2sg`_
   * - Delta Pseudopressure
     - `pyrestoolbox.gas.gas_dmp`_
   * - Gas Condensate FWS SG
     - `pyrestoolbox.gas.gas_fws_sg`_
   * - Gas Flow Rate Radial
     - `pyrestoolbox.gas.gas_rate_radial`_  
   * - Gas Flow Rate Linear 
     - `pyrestoolbox.gas.gas_rate_linear`_
  

pyrestoolbox.gas.gas_tc_pc
======================

.. code-block:: python

    gas_tc_pc(sg, co2 = 0, h2s = 0, n2 = 0, h2 = 0, cmethod = 'PMC', tc = 0, pc = 0) -> tuple

Returns a tuple of critical temperature (deg R) and critical pressure (psia) for hydrocarbon gas. If one or both of the tc and pc parameters are set to be non-zero, then this function will return that unchanged value for the corresponding critical parameter.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - sg
     - float
     - Gas SG relative to air   
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined 
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to zero if undefined. If positive fraction, cmethod will override to 'BUR'
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified  

Examples:

.. code-block:: python

    >>> gas.gas_tc_pc(sg=0.7, co2 = 0.15)
    (363.9387708314338, 738.3190067714969)
    
    >>> gas.gas_tc_pc(sg=0.7, co2 = 0.15, tc=365, cmethod='SUT')
    (365, 709.2389730048743)

pyrestoolbox.gas.gas_z
==================

.. code-block:: python

    gas_z(p, sg, degf, zmethod='DAK', cmethod='PMC', co2 = 0, h2s = 0, n2 = 0, h2 = 0, tc = 0, pc = 0) -> float or np.array

Returns gas Z-factor (either float or Numpy array depending upon type of p specified) using specified method. 
A float or list / array can be used for p, returning corresponding 1-D array of Z-Factors. The cmethod will be used to calculate critical gas parameters unless tc and/or pc are explicitly set to be non-zero. This option enables users to use precalculate gas critical properties and so avoid repeated duplicated critical property calculations when compute time is an issue


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float, list or np.array 
     - Gas pressure (psia)
   * - sg
     - float
     - Gas SG relative to air  
   * - degf
     - float
     - Reservoir Temperature (deg F)
   * - zmethod
     - string or z_method
     - Method for calculating gas Z-factor. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_. 
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined 
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to zero if undefined. zmethod and cmethod get overriden to 'BUR' if positive fraction.
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified  

Examples:

.. code-block:: python

    >>> gas.gas_z(p=1000, sg=0.75, degf=160, n2 = 0.02, co2 = 0.17)
    0.9140707840075585
    
    >>> gas.gas_z(p=1000, sg=0.75, degf=160, n2 = 0.02, co2 = 0.17, zmethod='HY')
    0.9141985223206194
    
    >>> gas.gas_z(p=[1000, 2000], sg=0.75, degf=160, cmethod='SUT', n2 = 0.02, co2 = 0.17)
    array([0.91920553, 0.87196032])
 
 def gas_ug(
    p: npt.ArrayLike,
    sg: float,
    degf: float,
    zmethod: z_method = z_method.DAK,
    cmethod: c_method = c_method.PMC,
    co2: float = 0,
    h2s: float = 0,
    n2: float = 0,
    h2: float = 0,
    tc: float = 0,
    pc: float = 0,
    zee: float = 0,
    ugz = False
    
pyrestoolbox.gas.gas_ug
===================

.. code-block:: python

    gas_ug(p, sg, degf, zmethod ='DAK', cmethod = 'PMC', co2 = 0, h2s = 0, n2 = 0, h2 = 0, tc = 0, pc = 0, zee = 0, ugz = False) -> float or np.array

Returns gas viscosity (cP) using Lee, Gonzalez & Eakin (1966) correlation unless the 'BUR' method for Z-Factor is selected in which case a tuned LBC method is used. 
A float or list / array can be used for p, returning corresponding 1-D array of gas viscosities. The cmethod will be used to calculate critical gas parameters unless tc and/or pc are explicitly set to be non-zero. This option enables users to use pre-calculate gas critical properties and so avoid repeated duplicated critical property calculations when compute time is an issue
Furnishing a positive value for zee means it will be used instead of calculating Z, and setting ugz to True will return the viscosity * Z product instead of viscosity alone


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float, list or np.array 
     - Gas pressure (psia)
   * - sg
     - float
     - Gas SG relative to air  
   * - degf
     - float
     - Reservoir Temperature (deg F)
   * - zmethod
     - string or z_method
     - Method for calculating gas Z-factor. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_. 
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined 
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to zero if undefined. Overrides methods to 'BUR' if positive fraction.
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified  
   * - zee
     - float, list or np.array 
     - Z-Factor value, or list of values of same length as P, in case recalculation of Z-Factors is not needed. If undefined, will trigger Z-Factor calculation.
   * - ugz
     - boolean
     - Boolean flag that if True returns ug * Z instead of ug 

Examples:

.. code-block:: python

    >>> gas.gas_ug(p=1000, sg=0.75, degf=180, zmethod ='HY', cmethod = 'SUT')
    0.0141231843661131
    
    >>> gas.gas_ug(p=1000, sg=0.75, degf=180)
    0.014114198868648963
    
    
pyrestoolbox.gas.gas_cg
===================

.. code-block:: python

    gas_cg(p, sg, degf, co2 = 0, h2s = 0, n2 = 0, h2 = 0, tc = 0, pc = 0, zmethod = 'DAK', cmethod ='PMC') -> float or np.array

Returns gas compressibility (1/psi) 
A float or list / array can be used for p, returning corresponding 1-D array of gas compressibility's. The cmethod will be used to calculate critical gas parameters unless tc and/or pc are explicitly set to be non-zero. This option enables users to use precalculated gas critical properties and so avoid repeated duplicated critical property calculations when compute time is an issue


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float, list or np.array 
     - Gas pressure (psia)
   * - sg
     - float
     - Gas SG relative to air  
   * - degf
     - float
     - Reservoir Temperature (deg F)
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined 
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to zero if undefined. If positive fraction, cmethod will override to 'BUR'
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified  

Examples:

.. code-block:: python

    >>> gas.gas_cg(p=2000, sg=0.68, degf=120, co2=0.05)
    0.0005375634134905346
    
    >>> gas.gas_cg(p=np.array([1000,2000]), sg=0.68, degf=120, co2=0.05)
    array([0.0011039 , 0.00053756])
    

pyrestoolbox.gas.gas_bg
===================

.. code-block:: python

    gas_bg(p, sg, degf, zmethod='DAK', cmethod = 'PMC', co2 = 0, h2s = 0, n2 = 0, h2 = 0, tc = 0, pc = 0) -> float or np.array

Returns gas formation volume factor (rcf/scf). 
A float or list / array can be used for p, returning corresponding 1-D array of gas FVF's. The cmethod will be used to calculate critical gas parameters unless tc and/or pc are explicitly set to be non-zero. This option enables users to use precalculate gas critical properties and so avoid repeated duplicated critical property calculations when compute time is an issue.


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float, list or np.array 
     - Gas pressure (psia)
   * - sg
     - float
     - Gas SG relative to air  
   * - degf
     - float
     - Reservoir Temperature (deg F)
   * - zmethod
     - string or z_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined 
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to zero if undefined. If positive fraction, cmethod will override to 'BUR'
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified  

Examples:

.. code-block:: python

    >>> gas.gas_bg (p=3000, sg=0.78, degf=240)
    0.005930983977679231
    
    >>> 1 / gas.gas_bg (p=[3000, 5000], sg=0.78, degf=240)
    array([168.60608691, 249.6801909 ])

pyrestoolbox.gas.gas_den
=====================

.. code-block:: python

    gas_den(p, sg, degf, zmethod ='DAK', cmethod ='PMC', co2 = 0, h2s = 0, n2 = 0, h2 = 0, tc = 0, pc = 0) -> float or np.array

Returns gas density (lb/cuft) 
A float or list / array can be used for p, returning corresponding 1-D array of gas densities. The cmethod will be used to calculate critical gas parameters unless tc and/or pc are explicitly set to be non-zero. This option enables users to use precalculate gas critical properties and so avoid repeated duplicated critical property calculations when compute time is an issue


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float, list or np.array 
     - Gas pressure (psia)
   * - sg
     - float
     - Gas SG relative to air  
   * - degf
     - float
     - Reservoir Temperature (deg F)
   * - zmethod
     - string or z_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined 
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to zero if undefined. If positive fraction, cmethod will override to 'BUR'
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified  

Examples:

.. code-block:: python

    >>> gas.gas_den (p=2000, sg=0.75, degf=150, zmethod ='HY', cmethod ='SUT', n2 = 0.02, co2 = 0.15, h2s = 0.02)
    7.728991860473501
    

pyrestoolbox.gas.gas_water_content
==============================

.. code-block:: python

    gas_water_content(p, degf) -> float

Returns saturated volume of water vapor in natural gas (stb/mmscf). From 'PVT and Phase Behaviour Of Petroleum Reservoir Fluids' by Ali Danesh.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p
     - float 
     - Gas pressure (psia)
   * - degf
     - float
     - Reservoir Temperature (deg F)

Examples:

.. code-block:: python

    >>> gas.gas_water_content(p=1500, degf=165)
    0.6521546577394491  

pyrestoolbox.gas.gas_ponz2p
=======================

.. code-block:: python

    gas_ponz2p(poverz, sg, degf, zmethod='DAK', cmethod='PMC', co2 = 0, h2s = 0, n2 = 0, h2 = 0, tc = 0, pc = 0, rtol = 1E-7) -> float or np.array

Returns gas pressure corresponding to a value of P/Z, iteratively solving with specified zmethod via bisection.
A float or list / array can be used for poverz, returning corresponding 1-D array of pressures. The cmethod will be used to calculate critical gas parameters unless tc and/or pc are explicitly set to be non-zero. This option enables users to use precalculate gas critical properties and so avoid repeated duplicated critical property calculations when compute time is an issue


.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - poverz
     - float, list or np.array 
     - Gas pressure / Z-factor (psia)
   * - sg
     - float
     - Gas SG relative to air  
   * - degf
     - float
     - Reservoir Temperature (deg F)
   * - zmethod
     - string or z_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined 
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to zero if undefined. If positive fraction, cmethod will override to 'BUR'
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified 
   * - rtol
     - float
     - relative solution tolerance as compared with abs([User P/Z - Calculated P/Z] / [User P/Z])

Examples:

.. code-block:: python

    >>> gas.gas_ponz2p(poverz=2500, sg=0.75, degf=165)
    2082.5648307800293   
    
    >>> gas.gas_ponz2p(poverz=[2500,5000], sg=0.75, degf=165)
    array([2082.56483078, 4890.62070847])
    
pyrestoolbox.gas.gas_grad2sg
========================

.. code-block:: python

    gas_grad2sg( grad, p, degf, zmethod='DAK', cmethod='PMC', co2 = 0, h2s = 0, n2 = 0, h2 = 0, tc = 0, pc = 0, rtol = 1E-7) -> float

Returns gas specific gravity consistent with observed gas gradient. Calculated through iterative solution method. Will fail if gas SG is below 0.55, or greater than 1.75

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - grad
     - float
     - Observed gas gradient (psi/ft)
   * - p
     - float, list or np.array 
     - Pressure at observation (psia)
   * - degf
     - float
     - Reservoir Temperature (deg F)
   * - zmethod
     - string or z_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined 
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to zero if undefined. If positive fraction, cmethod will override to 'BUR'
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified  
   * - rtol
     - float
     - relative solution tolerance as compared with abs([User grad - Calculated grad] / [User grad])
     
Examples:

.. code-block:: python

    >>> gas.gas_grad2sg(grad=0.0657, p=2500, degf=175)
    0.7500786632299423   
    

pyrestoolbox.gas.gas_dmp
=====================

.. code-block:: python

    gas_dmp(p1, p2, degf, sg, zmethod='DAK', cmethod = 'PMC', co2 = 0, h2s = 0, n2 = 0, h2 = 0, tc = 0, pc = 0) -> float

Returns gas pseudo-pressure integral between two pressure points. Will return a positive value if p1 < p2, and a negative value if p1 > p2. 
Integrates the equation: m(p) = 2 * p / (ug * z) 

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - p1
     - float, list or np.array 
     - First gas pressure (psia)
   * - p2
     - float, list or np.array 
     - Second gas pressure (psia)
   * - sg
     - float
     - Gas SG relative to air.
   * - degf
     - float
     - Reservoir Temperature (deg F)
   * - zmethod
     - string or z_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - cmethod
     - string or c_method
     - Method for calculating gas critical parameters. `Calculation Methods and Class Objects`_.
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined 
   * - h2
     - float
     - Molar fraction of Hydrogen. Defaults to zero if undefined. If positive fraction, cmethod will override to 'BUR'
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified  

Examples:

.. code-block:: python

    >>> gas.gas_dmp(p1=1000, p2=2000, degf=185, sg=0.78, zmethod='HY', cmethod = 'SUT', n2 = 0.05, co2 = 0.1, h2s = 0.02)
    3690873383.43637  
    
    >>> gas.gas_dmp(p1=2000, p2=1000, degf=185, sg=0.78, tc = 371, pc = 682)
    -3691052075.812854
        
pyrestoolbox.gas.gas_fws_sg
=======================

.. code-block:: python

    gas_fws_sg(sg_g, cgr, api_st) -> float

Estimates specific gravity of full-wellstream (FWS) gas from gas-condensate well. Calculates from weighted average surface gas SG, CGR and API. Uses Standing correlation to estimate condensate MW from API.
Returns SG of FWS gas 

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - sg_g
     - float
     - Specific gravity of weighted average surface gas (relative to air) 
   * - cgr
     - float
     - Condensate gas ratio (stb/mmscf). 
   * - api_st
     - float
     - Density of stock tank liquid (API)  

Examples:

.. code-block:: python

    >>> gas.gas_fws_sg(sg_g=0.855, cgr=30, api_st=53)
    0.9371015922844881
    
    
pyrestoolbox.gas.gas_rate_radial
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
    >>> gas.gas_rate_radial(k=5, h=50, pr=2000, pwf=750, r_w=0.3, r_ext=1500, degf=180, sg = 0.75, D = 0.01, S=5)
    10269.669190157822
    
    >>> gas.gas_rate_radial(k=1, h=50, pr=[2000,1000], pwf=750, r_w=0.3, r_ext=1500, degf=180, sg = 0.75, D = 0.01, S=5)
    array([4273.15956785, 1177.38697977])
    

pyrestoolbox.gas.gas_rate_linear
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

    >>> gas.gas_rate_linear(k=0.1, area=50, length=200, pr=2000, pwf=250, degf=180, sg = 0.8)
    21.87803915816601
    
    >>> gas.gas_rate_linear(k=0.1, area=50, length=200, pr=[2000, 1000, 500], pwf=250, degf=180, sg = 0.8)
    array([21.87803916,  4.89593662,  0.94342881])
    
