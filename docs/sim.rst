===================================
Simulation Helpers
===================================

Function List
=============

.. list-table:: Simulation Helpers
   :widths: 15 40
   :header-rows: 1

   * - Task
     - Function
   * - Extracting cells with convergence issues from IX PRT file
     - `pyrestoolbox.simtools.ix_extract_problem_cells`_  
   * - Solves Van Everdingin & Hurst Constant Terminal Rate solution via inverse Laplace transform and optionally writes out ECLIPSE styled AQUTAB include file
     - `pyrestoolbox.simtools.influence_tables`_
   * - Performs a recursive ECLIPSE deck zip/check for all INCLUDE files
     - `pyrestoolbox.simtools.zip_check_ecl`_
   * - Solves the Rachford Rice equation
     - `pyrestoolbox.simtools.rr_solver`_
     

pyrestoolbox.simtools.ix_extract_problem_cells
======================

.. code-block:: python

    ix_extract_problem_cells(filename = '', silent = False) -> list

Processes Intersect PRT file to extract convergence issue information. Prints a summary of worst offenders to terminal (if silent == False), and returns a list of sorted dataframes summarising all entities in final convergence report rows in the PRT file.

List returned is [well_pressure_df, grid_pressure_df, sat_change_df, comp_change_df]  

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - filename
     - str
     - If empty, will search local directory for PRT file and present list to select from if more than one exists. If a filename is furnished, or only one file exists, then no selection will be presented
   * - silent
     - bool
     - False will return only the list of dataframes, with nothing echoed to the terminal. True will print summary of worst entities to the terminal as well as returning the list.

Examples:

.. code-block:: python

    >>> results = rtb.simtools.ix_extract_problem_cells()
    >>> wells, grid_pres, grid_sat, grid_comp = results
    >>> grid_sat
    Processing TEST.PRT
    
    Issue Type                 Total Instances  Most Frequent Actor      Instances
    -----------------------  -----------------  ---------------------  -----------
    Well Pressure Change                     0  None                             0
    Grid Pressure Change                    64  32,125,4                         2
    Grid Saturation Change                 310  32,121,24                       13
    Grid Composition Change               1627  35,212,25                      544 
    
.. image:: https://github.com/mwburgoyne/pyResToolbox/blob/main/docs/img/grid_sat_df.png
    :alt: DSorted ataFrame of grid blocks with saturation related convergence issues    
    
    
    
pyrestoolbox.simtools.influence_tables
======================

.. code-block:: python

    influence_tables(ReDs, min_td = 0.01, max_td = 200, n_incr = 20, M = 8, export = False)-> tuple

Solves Van Everdingin & Hurst Constant Terminal Rate solution via inverse Laplace transform and optionally writes out ECLIPSE styled AQUTAB include file. 

Returns a tuple of;
    1. Dimensionless time list
    2. list of lists of dimensionless pressures at each dimensionless time for each dimensionless radius

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - ReDs
     - list
     - A list of dimensionless radii > 1.0. These are the ratios of the exteral radius of the aquifer to the external radius of the reservoir (or internal radius of the aquifer)
   * - min_td
     - float
     - Minimum dimensionless time. Default = 0.01
   * - max_td
     - float
     - Maximum dimensionless time. Default = 200
   * - n_incr
     - int
     - Number of log transformed increments to split dimensionless time into. Default = 20
   * - M
     - int
     - Laplace invesrion accuracy. Higher = more accurate, but more time. Generally 6-12 is good range. Default = 8
   * - export
     - bool
     - Boolean value that controls whether an include file with 'INFLUENCE.INC' name is created. Default: False 

Examples:

.. code-block:: python

    >>> from pyrestoolbox import pyrestoolbox as rtb
    >>> import matplotlib.pyplot as plt
    >>> ReDs = [1.5, 2, 3, 5, 10, 25, 1000]
    >>> tds, pds = rtb.simtools.influence_tables(ReDs=ReDs, export=True)
    >>> 
    >>> for p, pd in enumerate(pds):
    >>>     plt.plot(tds, pd, label = str(ReDs[p]))
    >>>     
    >>> plt.xscale('log')
    >>> plt.yscale('log')
    >>> plt.legend(loc='upper left')
    >>> plt.grid(which='both')
    >>> plt.xlabel('Dimensionless Time (tD)')
    >>> plt.ylabel('Dimensionless Pressure Drop (PD)')
    >>> plt.title('Constant Terminal Rate Solution')
    >>> plt.show()

.. image:: https://github.com/mwburgoyne/pyResToolbox/blob/main/docs/img/influence.png
    :alt: Constant Terminal Rate influence tables


pyrestoolbox.simtools.zip_check_ecl
======================

.. code-block:: python

    zip_check_ecl(deck, to_zip = True)

Performs a recursive ECLIPSE deck zip/check. 
Crawls through all INCLUDE files referenced in a deck, including an unlimited number of subdirectories and nested INCLUDE references, 
and (a) checks that all include files exist, then optionally (b) creates a zip file of all required files.
It does NOT zip any files that are in a higher directory than the .DATA file, but it does flag any such files so users can manually include them

Prints to console list of any missing files, and if desired creates a zip archive of all files required

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - deck
     - str
     - ECLIPSE deck File name including extension of the base ECLIPSE data file
   * - to_zip
     - bool
     - True will create a zip archive of DATA file and all associated INCLUDE files from the same directory or below. False will only return summary of whether INCLUDE files are complete.

Examples:

.. code-block:: python

    >>> zip_check_ecl = rtb.simtools.zip_check_ecl('FIELD_A.DATA')

pyrestoolbox.simtools.rr_solver
======================

.. code-block:: python

    rr_solver(zi, ki)

Solves for the root of the Rachford-Rice equation using a method that gracefully handles catastrophic roundoff errors.
The method is outlined in 2022 'Fluid Phase Equilibria' paper by M. Nielsen & H. Lia

Returns a tuple of results

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - zi
     - np.array
     - Molar composition (Percent, Fraction or Amounts - will be normalized)
   * - ki
     - np.array
     - K-values of the respective molar species

.. list-table:: Output Tuple
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - N_it
     - Integer
     - Number of iterations required to solve
   * - yi
     - np.array
     - Vapor mole fraction compositions
   * - xi
     - np.array
     - Liquid mole fraction compositions
   * - V
     - Float
     - Vapor molar fraction   
   * - L
     - Float
     - Liquid molar fraction
     
Examples:

.. code-block:: python

    >>> rtb.simtools.rr_solver(zi =np.array([0.7, 0.15, 0.1, 0.05]), ki = np.array([50, 5, 0.5, 0.01]))
    (6,
    array([0.7406252 , 0.1570315 , 0.09469948, 0.00764382]),
    array([0.0148125 , 0.0314063 , 0.18939896, 0.76438224]),
    0.9440279802330239,
    0.05597201976697608)
    
  