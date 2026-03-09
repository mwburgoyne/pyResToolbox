===================================
Decline Curve Analysis
===================================

Functions for fitting and forecasting production decline using Arps (exponential, hyperbolic, harmonic) and Duong models. Includes rate-vs-cumulative fitting, secondary phase ratio forecasting, uptime inference, and EUR estimation.

All functions are unit-agnostic — they operate on the numerical values you provide. If you pass rates in stb/d and time in months, the fitted ``qi`` comes back in stb/d and ``di`` in 1/month. No unit conversions are performed internally.


pyrestoolbox.dca.arps_rate
======================

.. code-block:: python

    arps_rate(qi, di, b, t) -> float or np.ndarray

Arps decline rate at time t.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - qi
     - float
     - Initial rate (volume/time)
   * - di
     - float
     - Initial decline rate (1/time). Must be > 0
   * - b
     - float
     - Arps b-factor. b=0: exponential, 0<b<1: hyperbolic, b=1: harmonic
   * - t
     - float or array
     - Time

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float or np.ndarray
     - Rate at time t. Returns same type as input t

Examples:

.. code-block:: python

    >>> from pyrestoolbox import dca
    >>> dca.arps_rate(qi=1000, di=0.1, b=0, t=10)
    367.87944117144235

    >>> dca.arps_rate(qi=1000, di=0.1, b=0.5, t=10)
    444.44444444444446

    >>> dca.arps_rate(qi=1000, di=0.1, b=1.0, t=10)
    500.0

    >>> dca.arps_rate(qi=1000, di=0.1, b=0, t=[0, 5, 10])
    array([1000.        ,  606.53065971,  367.87944117])


pyrestoolbox.dca.arps_cum
======================

.. code-block:: python

    arps_cum(qi, di, b, t) -> float or np.ndarray

Arps cumulative production at time t.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - qi
     - float
     - Initial rate (volume/time)
   * - di
     - float
     - Initial decline rate (1/time). Must be > 0
   * - b
     - float
     - Arps b-factor. b=0: exponential, 0<b<1: hyperbolic, b=1: harmonic
   * - t
     - float or array
     - Time

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float or np.ndarray
     - Cumulative production at time t. Returns same type as input t

Examples:

.. code-block:: python

    >>> dca.arps_cum(qi=1000, di=0.1, b=0, t=10)
    6321.205588285577

    >>> dca.arps_cum(qi=1000, di=0.1, b=1.0, t=10)
    6931.471805599453


pyrestoolbox.dca.duong_rate
======================

.. code-block:: python

    duong_rate(qi, a, m, t) -> float or np.ndarray

Duong decline rate for unconventional reservoirs. q(t) = qi * t^(-m) * exp(a/(1-m) * (t^(1-m) - 1)).

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - qi
     - float
     - Rate at t=1 (volume/time)
   * - a
     - float
     - Duong 'a' parameter. Must be > 0
   * - m
     - float
     - Duong 'm' parameter. Must be > 1
   * - t
     - float or array
     - Time (must be > 0)

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float or np.ndarray
     - Rate at time t. Returns same type as input t

Examples:

.. code-block:: python

    >>> dca.duong_rate(qi=500, a=1.5, m=1.2, t=1.0)
    500.0

    >>> dca.duong_rate(qi=500, a=1.5, m=1.2, t=10.0)
    502.3644755843416


pyrestoolbox.dca.eur
======================

.. code-block:: python

    eur(qi, di, b, q_min) -> float

Estimated ultimate recovery for Arps decline (cumulative production when rate reaches q_min).

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - qi
     - float
     - Initial rate
   * - di
     - float
     - Initial decline rate (1/time)
   * - b
     - float
     - Arps b-factor (0 to 1)
   * - q_min
     - float
     - Economic limit rate

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float
     - EUR (cumulative production when rate reaches q_min)

Examples:

.. code-block:: python

    >>> dca.eur(qi=1000, di=0.1, b=0, q_min=10)
    9900.0

    >>> dca.eur(qi=1000, di=0.1, b=0.5, q_min=10)
    18000.0


pyrestoolbox.dca.fit_decline
======================

.. code-block:: python

    fit_decline(t, q, method='best', t_start=None, t_end=None) -> DeclineResult

Fit a decline model to time-domain production data. Optional ``t_start`` and ``t_end`` parameters restrict fitting to a time window. Data outside the window is excluded, and the window is shifted to start at t=0, so the returned ``qi`` represents the rate at the window start.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - t
     - array-like
     - Time array
   * - q
     - array-like
     - Rate array (must be > 0)
   * - method
     - str
     - 'exponential', 'harmonic', 'hyperbolic', 'duong', or 'best' (default). 'best' tries all four and returns the highest R-squared
   * - t_start
     - float, optional
     - Start of fitting window (inclusive). Data before t_start is excluded
   * - t_end
     - float, optional
     - End of fitting window (inclusive). Data after t_end is excluded

.. list-table:: Returns (DeclineResult)
   :widths: 10 15 40
   :header-rows: 1

   * - Attribute
     - Type
     - Description
   * - method
     - str
     - Decline model name ('exponential', 'harmonic', 'hyperbolic', or 'duong')
   * - qi
     - float
     - Initial rate (volume/time)
   * - di
     - float
     - Initial decline rate (1/time). 0 for Duong
   * - b
     - float
     - Arps b-factor. 0 for exponential, 1 for harmonic. 0 for Duong
   * - a
     - float
     - Duong 'a' parameter. 0 for Arps models
   * - m
     - float
     - Duong 'm' parameter. 0 for Arps models
   * - r_squared
     - float
     - Coefficient of determination of the fit
   * - residuals
     - np.ndarray
     - Residual array (observed - predicted)

Examples:

.. code-block:: python

    >>> import numpy as np
    >>> t = np.arange(1, 51, dtype=float)
    >>> q = 1000 * np.exp(-0.05 * t)
    >>> result = dca.fit_decline(t, q, method='exponential')
    >>> result.method
    'exponential'
    >>> result.qi
    1000.0000000000007
    >>> result.di
    0.05000000000000006
    >>> result.r_squared
    1.0

Windowed fitting example:

.. code-block:: python

    >>> t = np.arange(1, 101, dtype=float)
    >>> q = 1000 * np.exp(-0.05 * t)
    >>> result = dca.fit_decline(t, q, method='exponential', t_start=20, t_end=60)
    >>> result.qi
    367.87944117144156
    >>> result.di
    0.049999999999999975
    >>> result.r_squared
    1.0


pyrestoolbox.dca.fit_decline_cum
======================

.. code-block:: python

    fit_decline_cum(Np, q, method='best', t_calendar=None, Np_start=None, Np_end=None) -> DeclineResult

Fit a decline model to rate-vs-cumulative data, eliminating time from the Arps equations. The returned ``qi`` and ``di`` are identical to time-domain parameters, so the result works directly with ``arps_rate()`` and ``forecast()``.

Supports exponential, harmonic, and hyperbolic models. Duong is excluded (no analytical q-vs-Np form) and raises ``ValueError``.

When ``t_calendar`` is provided, per-interval uptime fractions are inferred by comparing calendar-average rates to fitted capacity rates. Results are stored in ``uptime_mean`` and ``uptime_history`` on the returned ``DeclineResult``.

Optional ``Np_start`` and ``Np_end`` parameters restrict fitting to a cumulative production window. The window is shifted to start at Np=0, so the returned ``qi`` represents the rate at the window start.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - Np
     - array-like
     - Cumulative production array
   * - q
     - array-like
     - Rate array (must be > 0)
   * - method
     - str
     - 'exponential', 'harmonic', 'hyperbolic', or 'best' (default). Passing 'duong' raises ValueError
   * - t_calendar
     - array-like, optional
     - Calendar time stamps corresponding to each data point. When provided, uptime fractions are computed
   * - Np_start
     - float, optional
     - Start of fitting window on cumulative axis (inclusive)
   * - Np_end
     - float, optional
     - End of fitting window on cumulative axis (inclusive)

.. list-table:: Returns (DeclineResult)
   :widths: 10 15 40
   :header-rows: 1

   * - Attribute
     - Type
     - Description
   * - method
     - str
     - Decline model name ('exponential', 'harmonic', or 'hyperbolic')
   * - qi
     - float
     - Initial rate (volume/time)
   * - di
     - float
     - Initial decline rate (1/time)
   * - b
     - float
     - Arps b-factor. 0 for exponential, 1 for harmonic
   * - r_squared
     - float
     - Coefficient of determination of the fit
   * - residuals
     - np.ndarray
     - Residual array (observed - predicted)
   * - uptime_mean
     - float or None
     - Mean uptime fraction when t_calendar is provided
   * - uptime_history
     - np.ndarray or None
     - Per-interval uptime fractions when t_calendar is provided

Examples:

.. code-block:: python

    >>> qi, di = 1000.0, 0.05
    >>> t = np.arange(1, 51, dtype=float)
    >>> q = qi * np.exp(-di * t)
    >>> Np = np.array([float(dca.arps_cum(qi, di, 0, ti)) for ti in t])
    >>> result = dca.fit_decline_cum(Np, q, method='exponential')
    >>> result.method
    'exponential'
    >>> result.qi
    1000.0000000000002
    >>> result.di
    0.05000000000000004
    >>> result.r_squared
    1.0

Windowed cumulative fitting example:

.. code-block:: python

    >>> t = np.arange(1, 101, dtype=float)
    >>> q = 1000 * np.exp(-0.05 * t)
    >>> Np = np.array([float(dca.arps_cum(1000, 0.05, 0, ti)) for ti in t])
    >>> result = dca.fit_decline_cum(Np, q, method='exponential', Np_start=5000, Np_end=15000)
    >>> result.qi
    740.8182206817182
    >>> result.di
    0.05000000000000003
    >>> result.r_squared
    1.0


pyrestoolbox.dca.fit_ratio
======================

.. code-block:: python

    fit_ratio(x, ratio, method='best', domain='cum') -> RatioResult

Fit a ratio model (e.g. GOR, WOR) to data. Four models available: linear (R = a + b*x), exponential (R = a*exp(b*x)), power (R = a*x^b), and logistic (R = Rmax/(1 + c*exp(-b*x))).

The ``domain`` parameter ('cum' or 'time') is stored in the result and controls how ``forecast()`` evaluates the ratio — against cumulative production or time.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - x
     - array-like
     - Independent variable (cumulative production or time)
   * - ratio
     - array-like
     - Ratio values (e.g. GOR, WOR). Must be > 0 for exponential/power methods
   * - method
     - str
     - 'linear', 'exponential', 'power', 'logistic', or 'best' (default)
   * - domain
     - str
     - 'cum' or 'time' — stored in result for use by forecast()

.. list-table:: Returns (RatioResult)
   :widths: 10 15 40
   :header-rows: 1

   * - Attribute
     - Type
     - Description
   * - method
     - str
     - Ratio model name ('linear', 'exponential', 'power', or 'logistic')
   * - a
     - float
     - Primary parameter (intercept / coefficient / Rmax for logistic)
   * - b
     - float
     - Slope / exponent / growth rate
   * - c
     - float
     - Logistic offset parameter (only used for logistic model, 0 otherwise)
   * - domain
     - str
     - 'time' or 'cum' — tells forecast() which x-axis to use
   * - r_squared
     - float
     - Coefficient of determination of the fit
   * - residuals
     - np.ndarray
     - Residual array (observed - predicted)

Examples:

.. code-block:: python

    >>> x = np.arange(1, 51, dtype=float)
    >>> ratio = 0.5 + 0.02 * x
    >>> rr = dca.fit_ratio(x, ratio, method='linear', domain='cum')
    >>> rr.method
    'linear'
    >>> rr.a
    0.5000000000000001
    >>> rr.b
    0.019999999999999993
    >>> rr.r_squared
    1.0


pyrestoolbox.dca.ratio_forecast
======================

.. code-block:: python

    ratio_forecast(result, x) -> float or np.ndarray

Evaluate a fitted ratio model at given x values.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - result
     - RatioResult
     - Fitted ratio result from fit_ratio()
   * - x
     - float or array-like
     - Independent variable values to evaluate at

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float or np.ndarray
     - Ratio value at x. Returns same type as input x

Examples:

.. code-block:: python

    >>> rr = dca.RatioResult(method='linear', a=1.0, b=0.5, domain='cum')
    >>> dca.ratio_forecast(rr, 10.0)
    6.0


pyrestoolbox.dca.forecast
======================

.. code-block:: python

    forecast(result, t_end, dt=1.0, q_min=0.0, uptime=1.0, ratios=None) -> ForecastResult

Generate a rate and cumulative forecast from a fitted decline model. Supports uptime scaling and secondary phase ratio forecasting.

The ``uptime`` parameter scales capacity rate to calendar-effective rate (q_calendar = q_capacity * uptime). Default 1.0 preserves backward compatibility.

The ``ratios`` parameter accepts a dict mapping names to ``RatioResult`` objects. For each entry, secondary phase rate and cumulative are computed based on the ratio's domain ('cum' or 'time').

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - result
     - DeclineResult
     - Fitted decline result from fit_decline() or fit_decline_cum()
   * - t_end
     - float
     - End time for the forecast
   * - dt
     - float
     - Time step (default 1.0)
   * - q_min
     - float
     - Economic limit rate (default 0, no cutoff)
   * - uptime
     - float
     - Uptime fraction, 0 to 1 (default 1.0)
   * - ratios
     - dict, optional
     - Maps names to RatioResult objects for secondary phase forecasting

.. list-table:: Returns (ForecastResult)
   :widths: 10 15 40
   :header-rows: 1

   * - Attribute
     - Type
     - Description
   * - t
     - np.ndarray
     - Time array
   * - q
     - np.ndarray
     - Rate array (calendar-effective, i.e. capacity * uptime)
   * - Qcum
     - np.ndarray
     - Cumulative production array
   * - eur
     - float
     - Estimated ultimate recovery (final cumulative value)
   * - secondary
     - dict or None
     - Per-name dict with 'ratio', 'rate', 'cum' arrays when ratios are provided

Examples:

.. code-block:: python

    >>> dr = dca.DeclineResult(method='exponential', qi=1000, di=0.05, b=0)
    >>> fc = dca.forecast(dr, t_end=50, dt=1.0)
    >>> fc.eur
    17903.167013322283

    >>> fc2 = dca.forecast(dr, t_end=50, dt=1.0, uptime=0.8)
    >>> fc2.eur
    14322.533610657827

    >>> rr = dca.RatioResult(method='linear', a=0.5, b=0.001, domain='cum')
    >>> fc3 = dca.forecast(dr, t_end=50, dt=1.0, ratios={'GOR': rr})
    >>> fc3.secondary['GOR']['ratio'][0]
    1.451229424500714


Class Objects
======================

.. list-table::
   :widths: 15 40
   :header-rows: 1

   * - Class
     - Description
   * - DeclineResult
     - Result from ``fit_decline()`` or ``fit_decline_cum()``. See `Returns table above <#pyrestoolbox-dca-fit-decline>`__ for full attribute details
   * - ForecastResult
     - Result from ``forecast()``. See `Returns table above <#pyrestoolbox-dca-forecast>`__ for full attribute details
   * - RatioResult
     - Result from ``fit_ratio()``. See `Returns table above <#pyrestoolbox-dca-fit-ratio>`__ for full attribute details
