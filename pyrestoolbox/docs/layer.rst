===================================
Permeability Layering
===================================

A number of functions used to help characterize heterogeneity in permeability layering and/or create permeability array consistant with observed levels of heterogeneity and observed average permeability.

Heteogeneity is characterized via the Lorenz coefficient, and implemented via one of two methods, both of which can be described by a single Beta value;
  - Exponential ('EXP')
  - Langmuir    ('LANG')

Background on the approach taken as well as the derivation for the Exponential implementation can be found in `this LinkedIn article, <https://www.linkedin.com/pulse/loving-lorenz-new-life-old-parameter-mark-burgoyne/>`_

For the Langmuir formulation:
  - At a specific cumulative Phi.h -> Kh = Phi.h * VL / (Phi.h + PL)
  - Integrating this, subtracting 0.5 and doubling -> Lorenz = (VL - PL * VL * Ln(VL) + PL * VL * Ln(PL) - 0.5) * 2
  - Where PL = 1 / Beta and VL = 1 / Beta + 1
                


pyrestoolbox.layer.lorenz2b
======================

.. code-block:: python

    lorenz2b(lorenz, lrnz_method = 'EXP') -> float

Returns the Beta value consistent with the Lorenz coefficient given, and implementation method selected    

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - lorenz
     - float
     - Lorenz coefficient (0 < lorenz < 1)
   * - lrnz_method
     - float
     - Implementation method. Can be either 'EXP' or 'LANG'. Default is 'EXP'

Examples:

.. code-block:: python

    >>> from pyrestoolbox import layer
    >>> layer.lorenz2b(0.75, lrnz_method = 'LANG')
    16.139518537603912
    
    >>> layer.lorenz2b(0.75)
    7.978108090962671
    
pyrestoolbox.layer.lorenzfromb
======================

.. code-block:: python

    lorenzfromb(B: float, lrnz_method: str = 'EXP') -> float

Returns the Lorenz coefficient consistent with the Beta value given, and implementation method selected       

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - B
     - float
     - Beta value (B > 0)
   * - lrnz_method
     - float
     - Implementation method. Can be either 'EXP' or 'LANG'. Default is 'EXP'

Examples:

.. code-block:: python

    >>> layer.lorenzfromb(16.139518537603912, lrnz_method = 'LANG')
    0.750000182307895
    
    >>> layer.lorenzfromb(7.978108090962671)
    0.7500000108799212
    
pyrestoolbox.layer.lorenz_from_flow_fraction
======================

.. code-block:: python

    lorenz_from_flow_fraction(kh_frac, phih_frac, lrnz_method= 'EXP') -> float

Returns the Lorenz coefficient consistent with observed best flow fraction from a phi_h fraction       

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - kh_frac
     - float
     - The cumulative flow fraction of the best contributing flow unit ( 0 < kh_frac < 1 )
   * - phih_frac
     - float
     - The cumulative porosity thickness fraction of the best contributing flow unit ( phih_frac < kh_frac )
   * - lrnz_method
     - float
     - Implementation method. Can be either 'EXP' or 'LANG'. Default is 'EXP'

Examples:

60% of the observed flow comes from 15% of the net thickness
.. code-block:: python

    >>> lorenz_factor = layer.lorenz_from_flow_fraction(kh_frac=0.6, phih_frac=0.15)
    >>> lorenz_factor
    0.6759312029093838


pyrestoolbox.layer.lorenz_2_flow_frac
======================

.. code-block:: python

    lorenz_2_flow_frac(lorenz, phih_frac, lrnz_method = 'EXP', B = -1) -> float

Returns expected flow fraction from the best phi_h fraction, with a specified Lorenz coefficient.

If B is left default, then it will be calculated. If B is explictly specified > 0, then it will be used instead of the provided lorenz coefficient so as to eliminate repetitive solving for B.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - lorenz
     - float
     - Lorenz coefficient (0 < lorenz < 1). If B is provided, will ignore this parameter to be more efficient. If not, will calculate B from this parameter.
   * - phih_frac
     - float
     - The cumulative porosity thickness fraction of the best contributing flow unit ( 0 < phih_frac < 1 )
   * - lrnz_method
     - float
     - Implementation method. Can be either 'EXP' or 'LANG'. Default is 'EXP'
   * - B
     - float
     - Beta value (B > 0). Will calculate if only lorenz variable defined
     

Examples:

.. code-block:: python

    >>> layer.lorenz_2_flow_frac(lorenz=0.6759312029093838, phih_frac=0.15)
    0.6000001346893536
    
   
       
pyrestoolbox.lorenz_2_layers
======================

.. code-block:: python

    lorenz_2_layers(lorenz, k_avg, nlayers = 1, shuffle = False, lrnz_method = 'EXP', B = -1, phi_h_fracs = []) -> np.ndarray

Returns np.array of permeability values honoring a specified average permeability (assuming equal thickness layers unless list of phi_h_fracs is provided), with degree of heterogeneity consistant with specified Lorenz coefficient and method
        
If B is left default, then it will be calculated. If B is explictly specified > 0, then it will be used instead of the provided lorenz coefficient so as to eliminate repetitive solving for B.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - lorenz
     - float
     - Lorenz coefficient (0 < lorenz < 1). If B is provided, will ignore this parameter to be more efficient. If not, will calculate B from this parameter.
   * - k_avg
     - float
     - The thickness weighted average permeability of all the layers - Sum(kh) / h
   * - nlayers
     - int
     - The number of permeability layers desired (>1 needed unless a list of phi_h_fracs is supplied)
   * - shuffle
     - bool
     - Boolean flag to determine whether to return the permeability array in decreasing order (False), or random order (True). Default False
   * - lrnz_method
     - float
     - Implementation method. Can be either 'EXP' or 'LANG'. Default is 'EXP'
   * - B
     - float
     - Beta value (B > 0). Will calculate if only lorenz variable defined
   * - phi_h_fracs
     - list
     - Optional ability to specify a sorted list of phi_h fractions to calculate permeabilities for. If this list does not add to unity, then one additional layer permeability will be returned. The list needs to be in sorted order of best flow capacity to worst. If list adds to more than 1, it will be normalized
     

Examples:

.. code-block:: python

    >>> layer.lorenz_2_layers(lorenz = 0.67, nlayers = 5, k_avg = 10, shuffle = True)
    array([10.58944038,  0.29499066, 34.9323596 ,  3.21009656,  0.9731128 ])
    
    >>> layer.lorenz_2_layers(lorenz = 0.67, k_avg = 10, phi_h_fracs=[0.05, 0.5])
    array([51.72990694, 14.12556056,  0.77938749]) 
    
   

