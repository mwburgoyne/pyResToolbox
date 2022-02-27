===================================
Relative Permeability Tables Generation
===================================

Returns ECLIPSE style relative permeability tables for
  - Corey or LET family
  - SWOF, SGOF or SGFN table type

While there are many parameters, users need only define those relevant to their table / family selection

Function returns a DataFrame of results, and optionally exports to ECLIPSE styled include file.


pyrestoolbox.rel_perm_table
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


Examples:
    >>> from pyrestoolbox import pyrestoolbox as rtb
    >>> import matplotlib.pyplot as plt
    >>> df = rtb.rel_perm_table(rows=25, krtable='SGOF', krfamily='LET', kromax =1, krgmax =1, swc =0.2, sorg =0.15, Lo=2.5, Eo = 1.25, To = 1.75, Lg = 1.2, Eg = 1.5, Tg = 2.0)
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

    >>> df = rtb.rel_perm_table(rows=25, krtable='SWOF', kromax =1, krwmax =0.25, swc =0.15, swcr = 0.2, sorw =0.15, no=2.5, nw=1.5)
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
          