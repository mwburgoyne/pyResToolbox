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
   * - Performs a recursive deck zip/check for all INCLUDE files in ECL and IX files
     - `pyrestoolbox.simtools.zip_check_sim_deck`_
   * - Solves the Rachford Rice equation
     - `pyrestoolbox.simtools.rr_solver`_
   * - Generates ECLIPSE style relative permeability tables
     - `pyrestoolbox.simtools.rel_perm_table`_  
     

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


pyrestoolbox.simtools.zip_check_sim_deck
======================

.. code-block:: python

    zip_check_sim_deck(files2scrape = [], tozip = True, console_summary = True)

Performs a recursive zip/check on one or more ECL/IX decks. 
Crawls through all INCLUDE files referenced in a deck, including an unlimited number of subdirectories and nested INCLUDE references, 
and (a) checks that all include files exist, then optionally (b) creates a zip file of all required files.
It does NOT zip any files that are in a higher directory than the .DATA file, but it does flag any such files so users can manually include them

If files2scrape list is specified: 
 - No user prompt for any options will be given

If files2scrape list is NOT specified: 
 - User will be prompted to select one or more decks by their index number
 - Optionally creates a zip archive of available associated files

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - files2scrape
     - list
     - A list of file names to scrape. These must be .DATA and/or .AFI files. If no (or empty) list is passed, then user will be prompted to select one or more from the files existing in the current directory
   * - to_zip
     - bool
     - True will create a zip archive of DATA file and all associated INCLUDE files from the same directory or below. False will only return summary of whether INCLUDE files are complete.
   * - console_summary
     - bool
     - True will return verbose summary to terminal. False will return only list of missing files with no terminal echo.

Examples:

.. code-block:: python

    >>> rtb.simtools.zip_check_sim_deck(['FIELD_A.DATA', 'FIELD_B.afi'], console_summary=False)
    ['INCLUDE/GridOpts.inc', 'INCLUDE/ZCORN_COORD.GRDECL', 'EPS.ixf']
    
    >>> rtb.simtools.zip_check_sim_deck()
      Index  File Name
    -------  ------------
          0  FIELD_A.DATA
          1  FIELD_B.afi
     
    Please choose index(s) of file to parse separated by commas (0 - 1) :0,1
    Zip or Check files? (Z/c): c
    Scanning through: FIELD_A.DATA
    Scanning through: FIELD_B.afi
    
    ****** MISSING FILES ******
    
    INCLUDE/GridOpts.inc from FIELD_A.DATA
    INCLUDE/ZCORN_COORD.GRDECL from FIELD_A.DATA
    EPS.ixf from FIELD_B.DATA
    

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
    

pyrestoolbox.simtools.rel_perm_table
======================

.. code-block:: python

    rel_perm_table(rows, krtable='SWOF', krfamily='COR', kromax=1, krgmax=1, krwmax=1, swc=0, swcr=0, sorg=0, sorw=0, sgcr=0, no=1, nw=1, ng=1, Lw=1, Ew=1, Tw=1, Lo=1, Eo=1, To=1, Lg=1, Eg=1, Tg=1, export=False)-> pd.DataFrame:
  

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - rows
     - int
     - Integer value specifying the number of table rows desired
   * - krtable
     - str or kr_table class object
     - A string or kr_table Enum class that specifies one of three table type choices - SWOF, SGOF or SGWFN. Default is 'SWOF'.
   * - krfamily
     - str or kr_table class object
     - A string or kr_family Enum class that specifies one of two curve function choices - COR or LET. Default is 'COR'
   * - kromax
     - float
     - Max Kr relative to oil. Default value = 1
   * - krgmax
     - float
     - Max Kr relative to gas. Default value = 1
   * - krwmax
     - float
     - Max Kr relative to water with a second phase present. Default value = 1
   * - swc
     - float
     - Minimum water saturation. Default value = 0
   * - swcr
     - float
     - Maximum water saturation for imobile water. Default value = 0
   * - sorg
     - float
     - Maximum oil saturation relative to gas for imobile oil. Default value = 0
   * - sorw
     - float
     - Maximum oil saturation relative to water for imobile oil. Default value = 0
   * - sgcr
     - float
     - Maximum gas saturation relative to water for imobile gas. Default value = 0
   * - no, nw, ng
     - float
     - Corey exponents to oil, water and gas respectively. Default values = 1
   * - Lw, Ew, Tw
     - float
     - LET exponents to water. Default values = 1
   * - Lo, Eo, To
     - float
     - LET exponents to oil. Default values = 1
   * - Lg, Eg, Tg
     - float
     - LET exponents to gas. Default values = 1   
   * - export
     - bool
     - Boolean flag that controls whether an include file with same name as krtable is created. Default: False

.. list-table:: Method Variables & Class Objects
   :widths: 10 15 40
   :header-rows: 1

   * - Class Variable
     - Class Object 
     - Class Description & Options
   * - krfamily
     - kr_family
     - A string or kr_family Enum class that specifies one of two curve function choices. Defaults to 'COR'. 
       Options are:
        + 'COR': Corey Curve function
        + 'LET': LET Relative permeability function
   * - krtable
     - kr_table
     - A string or kr_table Enum class that specifies one of three table type choices. Default is 'SWOF'.
       Options are:
        + SWOF: Water / Oil table
        + SGOF: Gas / Oil table
        + SGFN: Gas / Water table
          
          
Examples:
    >>> from pyrestoolbox import pyrestoolbox as rtb
    >>> import matplotlib.pyplot as plt
    >>> df = rtb.simtools.rel_perm_table(rows=25, krtable='SGOF', krfamily='LET', kromax =1, krgmax =1, swc =0.2, sorg =0.15, Lo=2.5, Eo = 1.25, To = 1.75, Lg = 1.2, Eg = 1.5, Tg = 2.0)
    >>> plt.plot(df['Sg'], df['Krgo'], c = 'r', label='Gas')
    >>> plt.plot(df['Sg'], df['Krog'], c = 'g', label='Oil')
    >>> plt.title('SGOF Gas Oil LET Relative Permeability Curves')
    >>> plt.xlabel('Sg')
    >>> plt.ylabel('Kr')
    >>> plt.legend()
    >>> plt.grid('both')
    >>> plt.plot()

.. image:: https://github.com/mwburgoyne/pyResToolbox/blob/main/docs/img/sgof.png
    :alt: SGOF Relative Permeability Curves

.. code-block:: python

    >>> df = rtb.simtools.rel_perm_table(rows=25, krtable='SWOF', kromax =1, krwmax =0.25, swc =0.15, swcr = 0.2, sorw =0.15, no=2.5, nw=1.5)
    >>> plt.plot(df['Sw'], df['Krow'], c = 'g', label='Oil')
    >>> plt.plot(df['Sw'], df['Krwo'], c = 'b', label='Water')
    >>> plt.title('SWOF Water Oil Corey Relative Permeability Curves')
    >>> plt.xlabel('Sw')
    >>> plt.ylabel('Kr')
    >>> plt.legend()
    >>> plt.grid('both')
    >>> plt.plot()
    
.. image:: https://github.com/mwburgoyne/pyResToolbox/blob/main/docs/img/swof.png
    :alt: SWOF Relative Permeability Curves
    
   
    
