===================================
Method Recommendations
===================================

Decision-tree engine for selecting appropriate PVT, Z-factor, and VLP correlations based on fluid composition and well geometry.


pyrestoolbox.recommend.recommend_gas_methods
======================

.. code-block:: python

    recommend_gas_methods(sg=0.65, co2=0, h2s=0, n2=0, h2=0) -> dict

Recommend Z-factor and critical property methods for a gas composition. Returns a dict with keys ``'zmethod'`` and ``'cmethod'``, each mapping to a ``MethodRecommendation``.

Decision logic:

- H2 present → BNS mandatory (only method with H2 support)
- Inerts > 55% → BNS recommended (5-component PR-EOS handles extreme compositions)
- CO2 > 10% or H2S > 10% → DAK/PMC (moderate impurity corrections)
- Otherwise → DAK/PMC (standard defaults)

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - sg
     - float
     - Gas specific gravity (default 0.65)
   * - co2
     - float
     - CO2 mole fraction (default 0)
   * - h2s
     - float
     - H2S mole fraction (default 0)
   * - n2
     - float
     - N2 mole fraction (default 0)
   * - h2
     - float
     - H2 mole fraction (default 0)

.. list-table:: Returns (dict)
   :widths: 10 15 40
   :header-rows: 1

   * - Key
     - Type
     - Description
   * - 'zmethod'
     - MethodRecommendation
     - Recommended Z-factor method
   * - 'cmethod'
     - MethodRecommendation
     - Recommended critical property method

Examples:

.. code-block:: python

    >>> from pyrestoolbox import recommend
    >>> r = recommend.recommend_gas_methods()
    >>> r['zmethod'].recommended
    'DAK'
    >>> r['zmethod'].alternatives
    ['HY', 'WYW', 'BNS']

    >>> r = recommend.recommend_gas_methods(h2=0.1)
    >>> r['zmethod'].recommended
    'BNS'
    >>> r['zmethod'].mandatory
    True

    >>> r = recommend.recommend_gas_methods(co2=0.3, n2=0.3)
    >>> r['zmethod'].recommended
    'BNS'


pyrestoolbox.recommend.recommend_oil_methods
======================

.. code-block:: python

    recommend_oil_methods(api=35.0) -> dict

Recommend oil PVT correlation methods based on API gravity. Returns a dict with keys ``'pbmethod'``, ``'rsmethod'``, ``'bomethod'``.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - api
     - float
     - Oil API gravity (default 35.0)

.. list-table:: Returns (dict)
   :widths: 10 15 40
   :header-rows: 1

   * - Key
     - Type
     - Description
   * - 'pbmethod'
     - MethodRecommendation
     - Recommended bubble point method
   * - 'rsmethod'
     - MethodRecommendation
     - Recommended solution GOR method
   * - 'bomethod'
     - MethodRecommendation
     - Recommended oil FVF method

Examples:

.. code-block:: python

    >>> r = recommend.recommend_oil_methods(api=35)
    >>> r['pbmethod'].recommended
    'VELAR'
    >>> r['rsmethod'].recommended
    'VELAR'
    >>> r['bomethod'].recommended
    'MCAIN'
    >>> r['pbmethod'].alternatives
    ['VALMC', 'STAN']


pyrestoolbox.recommend.recommend_vlp_method
======================

.. code-block:: python

    recommend_vlp_method(deviation=0, well_type='gas') -> dict

Recommend VLP multiphase flow correlation based on well deviation. Returns a dict with key ``'vlp_method'``.

HB and Gray were developed for vertical flow data and become unreliable beyond ~30 deg deviation. BB and WG are suitable for all inclinations.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - deviation
     - float
     - Maximum wellbore deviation from vertical in degrees (default 0)
   * - well_type
     - str
     - 'gas' or 'oil' (default 'gas')

.. list-table:: Returns (dict)
   :widths: 10 15 40
   :header-rows: 1

   * - Key
     - Type
     - Description
   * - 'vlp_method'
     - MethodRecommendation
     - Recommended VLP correlation

Examples:

.. code-block:: python

    >>> r = recommend.recommend_vlp_method(deviation=0)
    >>> r['vlp_method'].recommended
    'BB'
    >>> r['vlp_method'].alternatives
    ['HB', 'WG', 'GRAY']

    >>> r = recommend.recommend_vlp_method(deviation=60)
    >>> r['vlp_method'].recommended
    'BB'
    >>> r['vlp_method'].alternatives
    ['WG']


pyrestoolbox.recommend.recommend_methods
======================

.. code-block:: python

    recommend_methods(sg=0.65, co2=0, h2s=0, n2=0, h2=0,
                      api=None, deviation=0, well_type='gas') -> dict

Master recommendation function combining gas, oil, and VLP recommendations. Oil recommendations are included only when ``api`` is provided.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - sg
     - float
     - Gas specific gravity (default 0.65)
   * - co2
     - float
     - CO2 mole fraction (default 0)
   * - h2s
     - float
     - H2S mole fraction (default 0)
   * - n2
     - float
     - N2 mole fraction (default 0)
   * - h2
     - float
     - H2 mole fraction (default 0)
   * - api
     - float, optional
     - Oil API gravity. If provided, oil method recommendations are included
   * - deviation
     - float
     - Maximum wellbore deviation from vertical in degrees (default 0)
   * - well_type
     - str
     - 'gas' or 'oil' (default 'gas')

.. list-table:: Returns (dict)
   :widths: 10 15 40
   :header-rows: 1

   * - Key
     - Type
     - Description
   * - 'zmethod'
     - MethodRecommendation
     - Recommended Z-factor method
   * - 'cmethod'
     - MethodRecommendation
     - Recommended critical property method
   * - 'vlp_method'
     - MethodRecommendation
     - Recommended VLP correlation
   * - 'pbmethod'
     - MethodRecommendation
     - Recommended bubble point method (only when api is provided)
   * - 'rsmethod'
     - MethodRecommendation
     - Recommended solution GOR method (only when api is provided)
   * - 'bomethod'
     - MethodRecommendation
     - Recommended oil FVF method (only when api is provided)

Examples:

.. code-block:: python

    >>> r = recommend.recommend_methods(sg=0.7, co2=0.05, api=30, deviation=45)
    >>> sorted(r.keys())
    ['bomethod', 'cmethod', 'pbmethod', 'rsmethod', 'vlp_method', 'zmethod']
    >>> r['zmethod'].recommended
    'DAK'
    >>> r['vlp_method'].recommended
    'BB'
    >>> r['vlp_method'].alternatives
    ['WG']


Class Objects
======================

.. list-table::
   :widths: 15 40
   :header-rows: 1

   * - Class
     - Description
   * - MethodRecommendation
     - Single method recommendation. Attributes: category, recommended (str), rationale (str), alternatives (list of str), mandatory (bool)
