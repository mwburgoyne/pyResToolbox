===================================
Gas PVT
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

Examples:

Calculating bubble point pressure with Standing correlation via option string, and then via class object

.. code-block:: python

    >>> rtb.oil_pbub(api=43, degf=185, rsb=2350, sg_g =0.72, pbmethod ='STAN')
    6406.067846808766
    
    >>> rtb.oil_pbub(api=43, degf=185, rsb=2350, sg_g =0.72, pbmethod = rtb.pb_method.STAN)
    6406.067846808766


Function List
=============

.. list-table:: Gas Functions
   :widths: 15 40
   :header-rows: 1

   * - Task
     - Function
   * - Gas Tc & Pc Calculation
     - `pyrestoolbox.gas_tc_pc`_  
   * - Gas Z-Factor Calculation
     - `pyrestoolbox.gas_z`_
   * - Gas Viscosity
     - `pyrestoolbox.gas_ug`_
   * - Gas Viscosity * Z
     - `pyrestoolbox.gas_ugz`_
   * - Gas Compressibility
     - `pyrestoolbox.gas_cg`_
   * - Gas Formation Volume Factor
     - `pyrestoolbox.gas_bg`_  
   * - Gas Density
     - `pyrestoolbox.gas_den`_  
   * - Gas Water of Condensation
     - `pyrestoolbox.gas_water_content`_
   * - Convert P/Z to P
     - `pyrestoolbox.gas_ponz2p`_
   * - Convert Gas Gradient to SG
     - `pyrestoolbox.gas_grad2sg`_
   * - Delta Pseudopressure
     - `pyrestoolbox.gas_dmp`_
   * - Gas Condensate FWS SG
     - `pyrestoolbox.gas_fws_sg`_
  

pyrestoolbox.gas_tc_pc
======================

.. code-block:: python

    gas_tc_pc(sg, n2 = 0, co2 = 0, h2s = 0, cmethod = 'PMC', tc = 0, pc = 0) -> tuple

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
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined  
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
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

    >>> rtb.gas_tc_pc(sg=0.7, co2 = 0.15)
    (363.9387708314338, 738.3190067714969)
    
    >>> rtb.gas_tc_pc(sg=0.7, co2 = 0.15, tc=365, cmethod='SUT')
    (365, 709.2389730048743)

pyrestoolbox.gas_z
==================

.. code-block:: python

    gas_z(p, sg, degf, zmethod='DAK', cmethod='PMC', n2 = 0, co2 = 0, h2s = 0, tc = 0, pc = 0) -> float or np.array

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
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined  
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified  

Examples:

.. code-block:: python

    >>> rtb.gas_z(p=1000, sg=0.75, degf=160, n2 = 0.02, co2 = 0.17)
    0.9140707840075585
    
    >>> rtb.gas_z(p=1000, sg=0.75, degf=160, n2 = 0.02, co2 = 0.17, zmethod='LIN')
    0.9131105248098116
    
    >>> rtb.gas_z(p=[1000, 2000], sg=0.75, degf=160, cmethod='SUT', n2 = 0.02, co2 = 0.17)
    array([0.91920553, 0.87196032])
    
pyrestoolbox.gas_ug
===================

.. code-block:: python

    gas_ug(p, sg, degf, zmethod ='DAK', cmethod = 'PMC', n2 = 0, co2 = 0, h2s = 0, tc = 0, pc = 0) -> float or np.array

Returns gas viscosity (cP) using Lee, Gonzalez & Eakin (1966) correlation. 
A float or list / array can be used for p, returning corresponding 1-D array of gas viscosities. The cmethod will be used to calculate critical gas parameters unless tc and/or pc are explicitly set to be non-zero. This option enables users to use pre-calculate gas critical properties and so avoid repeated duplicated critical property calculations when compute time is an issue


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
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined  
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified  

Examples:

.. code-block:: python

    >>> rtb.gas_ug(p=1000, sg=0.75, degf=180, zmethod ='HY', cmethod = 'SUT')
    0.0141231843661131
    
    >>> rtb.gas_ug(p=1000, sg=0.75, degf=180)
    0.014114198868648963
    
pyrestoolbox.gas_ugz
====================

.. code-block:: python

    gas_ugz(p, sg, degf, zee) -> float or np.array

Returns gas viscosity*Z-factor product (cP) using Lee, Gonzalez & Eakin (1966) correlation, utilizing a precaculated Z-factor
A float or list / array can be used for p and zee, returning a 1-D array of gas viscosity*Z-factor products. 
Using the gas_ugz function instead of the product of the gas_ug and gas_z functions removes duplications in calculating the z-factor as well as the critical properties


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
   * - zee
     - float, list or numpy array
     - Gas Z-factor(s)

Examples:

.. code-block:: python

    >>> rtb.gas_ugz(p=[1000,2000], sg=0.75, degf=140, zee=[0.9,1.0])
    array([0.01219254, 0.01600964])
    
    >>> rtb.gas_ugz(p=1000, sg=0.75, degf=140, zee=0.9)
    0.012192537840814146
    
    
pyrestoolbox.gas_cg
===================

.. code-block:: python

    gas_cg(p, sg, degf, n2 = 0, co2 = 0, h2s = 0, tc = 0, pc = 0, cmethod ='PMC') -> float or np.array

Returns gas compressibility (1/psi) using the 'DAK' Dranchuk & Abou-Kassem (1975) Z-Factor & Critical property correlation values if tc and/or pc not explicitly specified
A float or list / array can be used for p, returning corresponding 1-D array of gas compressibility's. The cmethod will be used to calculate critical gas parameters unless tc and/or pc are explicitly set to be non-zero. This option enables users to use precalculate gas critical properties and so avoid repeated duplicated critical property calculations when compute time is an issue


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
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined  
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified  

Examples:

.. code-block:: python

    >>> rtb.gas_cg(p=2000, sg=0.68, degf=120, co2=0.05)
    0.0005375634134905346
    
    >>> rtb.gas_cg(p=np.array([1000,2000]), sg=0.68, degf=120, co2=0.05)
    array([0.0011039 , 0.00053756])
    

pyrestoolbox.gas_bg
===================

.. code-block:: python

    gas_bg(p, sg, degf, zmethod='DAK', cmethod = 'PMC', n2 = 0, co2 = 0, h2s = 0, tc = 0, pc = 0) -> float or np.array

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
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined  
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified  

Examples:

.. code-block:: python

    >>> rtb.gas_bg (p=3000, sg=0.78, degf=240)
    0.005930983977679231
    
    >>> 1 / rtb.gas_bg (p=[3000, 5000], sg=0.78, degf=240)
    array([168.60608691, 249.6801909 ])

pyrestoolbox.gas_den
=====================

.. code-block:: python

    gas_den(p, sg, degf, zmethod ='DAK', cmethod ='PMC', n2 = 0, co2 = 0, h2s = 0, tc = 0, pc = 0) -> float or np.array

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
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined  
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified  

Examples:

.. code-block:: python

    >>> rtb.gas_den (p=2000, sg=0.75, degf=150, zmethod ='HY', cmethod ='SUT', n2 = 0.02, co2 = 0.15, h2s = 0.02)
    7.728991860473501
    

pyrestoolbox.gas_water_content
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

    >>> rtb.gas_water_content(p=1500, degf=165)
    0.6521546577394491  

pyrestoolbox.gas_ponz2p
=======================

.. code-block:: python

    gas_ponz2p(poverz, sg, degf, zmethod='DAK', cmethod='PMC', n2 = 0, co2 = 0, h2s = 0, tc = 0, pc = 0, rtol = 1E-7) -> float or np.array

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
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined  
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
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

    >>> rtb.gas_ponz2p(poverz=2500, sg=0.75, degf=165)
    2082.5648307800293   
    
    >>> rtb.gas_ponz2p(poverz=[2500,5000], sg=0.75, degf=165)
    array([2082.56483078, 4890.62070847])
    
pyrestoolbox.gas_grad2sg
========================

.. code-block:: python

    gas_grad2sg( grad, p, degf, zmethod='DAK', cmethod='PMC', n2 = 0, co2 = 0, h2s = 0, tc = 0, pc = 0, rtol = 1E-7) -> float

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
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined  
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
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

    >>> rtb.gas_grad2sg(grad=0.0657, p=2500, degf=175)
    0.7500786632299423   
    

pyrestoolbox.gas_dmp
=====================

.. code-block:: python

    gas_dmp(p1, p2, degf, sg, zmethod='DAK', cmethod = 'PMC', n2 = 0, co2 = 0, h2s = 0, tc = 0, pc = 0) -> float

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
   * - n2
     - float
     - Molar fraction of Nitrogen. Defaults to zero if undefined  
   * - co2
     - float
     - Molar fraction of CO2. Defaults to zero if undefined 
   * - h2s
     - float
     - Molar fraction of H2S. Defaults to zero if undefined
   * - tc
     - float
     - Critical gas temperature (deg R). Uses cmethod correlation if not specified  
   * - pc
     - float
     - Critical gas pressure (psia). Uses cmethod correlation if not specified  

Examples:

.. code-block:: python

    >>> rtb.gas_dmp(p1=1000, p2=2000, degf=185, sg=0.78, zmethod='HY', cmethod = 'SUT', n2 = 0.05, co2 = 0.1, h2s = 0.02)
    3690873383.43637  
    
    >>> rtb.gas_dmp(p1=2000, p2=1000, degf=185, sg=0.78, tc = 371, pc = 682)
    -3691052075.812854
        
pyrestoolbox.gas_fws_sg
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

    >>> rtb.gas_fws_sg(sg_g=0.855, cgr=30, api_st=53)
    0.9371015922844881
