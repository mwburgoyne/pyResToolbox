===================================
Sensitivity Analysis
===================================

Framework for parameter sweeps and tornado-chart sensitivity analysis. Works with any callable function — not limited to pyResToolbox functions.


pyrestoolbox.sensitivity.sweep
======================

.. code-block:: python

    sweep(func, base_kwargs, vary_param, vary_values, result_key=None) -> SweepResult

Vary one parameter across a range, collecting results. Useful for generating sensitivity curves (e.g. Z-factor vs pressure).

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - func
     - callable
     - Function to evaluate
   * - base_kwargs
     - dict
     - Base keyword arguments for the function
   * - vary_param
     - str
     - Name of the parameter to vary (must be a key in base_kwargs)
   * - vary_values
     - list
     - Values to assign to vary_param
   * - result_key
     - str, optional
     - Key or attribute to extract from each result. If None, the raw result is stored

Examples:

.. code-block:: python

    >>> from pyrestoolbox import sensitivity
    >>> from pyrestoolbox import gas
    >>> s = sensitivity.sweep(
    ...     func=gas.gas_z,
    ...     base_kwargs=dict(p=2000, sg=0.7, degf=200),
    ...     vary_param='p',
    ...     vary_values=[1000, 2000, 3000, 4000]
    ... )
    >>> s.param
    'p'
    >>> s.results[0]
    0.9260188251531628
    >>> s.results[3]
    0.9411603827082226


pyrestoolbox.sensitivity.tornado
======================

.. code-block:: python

    tornado(func, base_kwargs, ranges, result_key=None) -> TornadoResult

Compute tornado-chart sensitivities for multiple parameters. Each parameter is varied to its low and high value while others remain at base. Results are sorted by decreasing sensitivity.

Sensitivity is defined as \|high_result - low_result\| / \|base_result\|.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - func
     - callable
     - Function to evaluate
   * - base_kwargs
     - dict
     - Base keyword arguments for the function
   * - ranges
     - dict
     - Mapping of parameter names to (low, high) tuples
   * - result_key
     - str, optional
     - Key or attribute to extract a scalar from each result

Examples:

.. code-block:: python

    >>> t = sensitivity.tornado(
    ...     func=gas.gas_z,
    ...     base_kwargs=dict(p=2000, sg=0.7, degf=200),
    ...     ranges={'p': (1000, 4000), 'sg': (0.6, 0.8), 'degf': (150, 250)},
    ... )
    >>> t.base_result
    0.8886670011404194
    >>> t.entries[0].param
    'sg'
    >>> t.entries[0].sensitivity
    0.0886413428394408
    >>> t.entries[1].param
    'degf'
    >>> t.entries[2].param
    'p'


Class Objects
======================

.. list-table::
   :widths: 15 40
   :header-rows: 1

   * - Class
     - Description
   * - SweepResult
     - Result from ``sweep()``. Attributes: param (str), values (list), results (list)
   * - TornadoEntry
     - Single parameter sensitivity. Attributes: param, low_value, high_value, low_result, high_result, sensitivity
   * - TornadoResult
     - Result from ``tornado()``. Attributes: base_result (float), entries (list of TornadoEntry, sorted by decreasing sensitivity)
