===================================
Simulation Helpers
===================================


pyrestoolbox.ix_extract_problem_cells
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

    >>> results = rtb.ix_extract_problem_cells()
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
    
    
    
pyrestoolbox.influence_tables
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

    >>> import matplotlib.pyplot as plt
    >>> ReDs = [1.5, 2, 3, 5, 10, 25, 1000]
    >>> tds, pds = rtb.influence_tables(ReDs=ReDs, export=True)
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
    
