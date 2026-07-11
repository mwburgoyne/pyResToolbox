===================================
Decline Curve Analysis
===================================

Functions for fitting and forecasting production decline using Arps (exponential, hyperbolic, harmonic), modified hyperbolic (hyperbolic-to-exponential), two-segment hyperbolic (transient-to-boundary-dominated) and Duong models. Includes rate-vs-cumulative fitting, EUR-constrained type-curve generation, secondary phase ratio forecasting, uptime inference, and EUR estimation.

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
     - Arps b-factor. b=0: exponential, 0<b<1: hyperbolic, b=1: harmonic, b>1: super-hyperbolic (transient flow, common in unconventionals)
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

    >>> dca.arps_rate(qi=1000, di=0.005, b=2.0, t=100)
    707.1067811865474

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
     - Arps b-factor. b=0: exponential, 0<b<1: hyperbolic, b=1: harmonic, b>1: super-hyperbolic (cumulative is unbounded as t increases unless paired with a terminal decline - see mh_cum)
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


pyrestoolbox.dca.mh_rate
======================

.. code-block:: python

    mh_rate(qi, di, t, b=2.0, dterm=0.0) -> float or np.ndarray

Modified hyperbolic (two-segment Arps) decline rate, after Robertson (1988), SPE-18731-MS. Hyperbolic decline with b-factor ``b`` until the nominal decline D(t) = di / (1 + b*di*t) falls to ``dterm``, exponential decline at ``dterm`` thereafter. Rate and nominal decline are both continuous at the switch, which occurs at t_sw = (1/dterm - 1/di) / b.

This is the standard form for unconventional wells, where transient linear flow gives an early-time b-factor above 1 (default 2.0) that must transition to boundary-dominated behavior for cumulative volumes to stay finite. With ``dterm=0`` (or ``b=0``) the curve is a single Arps segment identical to ``arps_rate``.

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
   * - t
     - float or array
     - Time
   * - b
     - float
     - First-segment Arps b-factor (default 2.0, transient flow). Any b >= 0 is accepted
   * - dterm
     - float
     - Terminal nominal decline rate (1/time). Must be < di. Default 0 disables the terminal segment

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float or np.ndarray
     - Rate at time t. Returns same type as input t

Examples (qi = 1000 stb/d, di = 0.005/day nominal, default b = 2, terminal decline 0.0005/day; the switch falls at t = 900 days):

.. code-block:: python

    >>> dca.mh_rate(qi=1000, di=0.005, t=100, dterm=0.0005)
    707.1067811865474

    >>> dca.mh_rate(qi=1000, di=0.005, t=2000, dterm=0.0005)
    182.44754964045953

    >>> dca.mh_rate(qi=1000, di=0.005, t=[100, 900, 2000], dterm=0.0005)
    array([707.10678119, 316.22776602, 182.44754964])


pyrestoolbox.dca.mh_cum
======================

.. code-block:: python

    mh_cum(qi, di, t, b=2.0, dterm=0.0) -> float or np.ndarray

Modified hyperbolic (two-segment Arps) cumulative production: the analytic piecewise integral of ``mh_rate`` (Arps cumulative to the switch time, exponential-segment cumulative beyond it).

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
   * - t
     - float or array
     - Time
   * - b
     - float
     - First-segment Arps b-factor (default 2.0, transient flow)
   * - dterm
     - float
     - Terminal nominal decline rate (1/time). Must be < di. Default 0 disables the terminal segment

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

    >>> dca.mh_cum(qi=1000, di=0.005, t=2000, dterm=0.0005)
    700015.9647864327


pyrestoolbox.dca.mh_eur
======================

.. code-block:: python

    mh_eur(qi, di, q_min, b=2.0, dterm=0.0) -> float

Estimated ultimate recovery for modified hyperbolic decline (cumulative production when rate reaches q_min). If abandonment falls within the hyperbolic segment the result equals ``eur()``; beyond the switch the exponential-segment integral (q_sw - q_min)/dterm is added to the switch-point cumulative.

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
   * - q_min
     - float
     - Economic limit rate. Must be > 0 when b >= 1 and dterm = 0 (single-segment cumulative is unbounded)
   * - b
     - float
     - First-segment Arps b-factor (default 2.0, transient flow)
   * - dterm
     - float
     - Terminal nominal decline rate (1/time). Must be < di. Default 0 disables the terminal segment

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

    >>> dca.mh_eur(qi=1000, di=0.005, q_min=5.0, dterm=0.0005)
    1054911.0640673516


pyrestoolbox.dca.hyp2_rate
======================

.. code-block:: python

    hyp2_rate(qi, di, t, telf, b1=2.0, b2=0.5, dterm=0.0) -> float or np.ndarray

Two-segment hyperbolic decline rate with the transition at a specified time telf. Hyperbolic decline with b-factor ``b1`` (default 2.0, transient linear flow) to telf, then a second hyperbolic segment with b-factor ``b2`` anchored so that rate and nominal decline (log-slope) are identical to the first segment at telf: q_sw = qi/(1 + b1*di*telf)^(1/b1) and D_sw = di/(1 + b1*di*telf). This is the sharp-transition form of the transient hyperbolic model (Fulford and Blasingame 2013, SPE-167242-MS) and the two-stage type-curve construction of Burgoyne (2016).

An optional terminal decline ``dterm`` converts the second segment to exponential once its nominal decline falls to dterm, exactly as in ``mh_rate``.

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
     - Initial nominal decline rate (1/time). Must be > 0
   * - t
     - float or array
     - Time
   * - telf
     - float
     - Transition time (end of linear flow), same units as t
   * - b1
     - float
     - First-segment Arps b-factor (default 2.0, transient flow)
   * - b2
     - float
     - Second-segment Arps b-factor (default 0.5, boundary-dominated). Must not exceed b1
   * - dterm
     - float
     - Terminal nominal decline rate (1/time) applied to the second segment. Default 0 disables it

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - float or np.ndarray
     - Rate at time t. Returns same type as input t

Examples (qi = 1000 stb/d, di = 0.005/day nominal, transition at 730 days, default b1 = 2 and b2 = 0.5):

.. code-block:: python

    >>> dca.hyp2_rate(qi=1000, di=0.005, t=100, telf=730)
    707.1067811865474

    >>> dca.hyp2_rate(qi=1000, di=0.005, t=2000, telf=730)
    181.59828808766977

    >>> dca.hyp2_rate(qi=1000, di=0.005, t=[100, 730, 2000], telf=730)
    array([707.10678119, 347.10506725, 181.59828809])


pyrestoolbox.dca.hyp2_cum
======================

.. code-block:: python

    hyp2_cum(qi, di, t, telf, b1=2.0, b2=0.5, dterm=0.0) -> float or np.ndarray

Two-segment hyperbolic cumulative production: the analytic piecewise integral of ``hyp2_rate``. Inputs are identical to ``hyp2_rate``.

Examples:

.. code-block:: python

    >>> dca.hyp2_cum(qi=1000, di=0.005, t=2000, telf=730)
    695047.0925841478


pyrestoolbox.dca.hyp2_eur
======================

.. code-block:: python

    hyp2_eur(qi, di, q_min, telf, b1=2.0, b2=0.5, dterm=0.0) -> float

Estimated ultimate recovery for two-segment hyperbolic decline (cumulative production when rate reaches q_min). If abandonment falls within the first segment the result equals ``eur()``. ``q_min`` must be > 0 when b2 >= 1 and dterm = 0 (unbounded cumulative).

Examples:

.. code-block:: python

    >>> dca.hyp2_eur(qi=1000, di=0.005, q_min=10.0, telf=730)
    1332983.3654469885


pyrestoolbox.dca.hyp2_from_eur
======================

.. code-block:: python

    hyp2_from_eur(eur_target, qi, q_min, b1=2.0, b2=0.5,
                  telf=None, d_sw=None, eur_frac=None, rate_frac=None) -> DeclineResult

Solve the first-segment nominal decline ``di`` so that a two-segment hyperbolic delivers ``eur_target`` when the rate reaches ``q_min``. The transition is specified by exactly one of:

- ``telf`` - transition at a fixed time
- ``d_sw`` - transition when the nominal decline falls to d_sw (1/time)
- ``eur_frac`` - transition after this fraction of eur_target is produced
- ``rate_frac`` - transition when the rate falls to rate_frac * qi

Port of the two-stage type-curve generator logic (Burgoyne 2016-2018), with the original bisection replaced by a bracketed brentq root find that keeps the transition time self-consistent at every iteration.

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - eur_target
     - float
     - Target EUR (volume) at abandonment rate q_min
   * - qi
     - float
     - Initial rate (volume/time)
   * - q_min
     - float
     - Abandonment rate. Must satisfy 0 < q_min < qi
   * - b1
     - float
     - First-segment Arps b-factor (default 2.0). Must be > 0
   * - b2
     - float
     - Second-segment Arps b-factor (default 0.5). Must not exceed b1
   * - telf / d_sw / eur_frac / rate_frac
     - float, optional
     - Transition specification - supply exactly one

.. list-table:: Returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * -
     - DeclineResult
     - method='hyp2' with the solved di, b=b1, b2 and resolved telf; works directly with forecast(), hyp2_rate() and hyp2_eur(). For the d_sw specification, if the solved decline never reaches d_sw before abandonment the result is single-segment (method='hyperbolic')

Examples (target EUR of 1.5e6 mscf at 100 mscf/d abandonment from 1500 mscf/d initial rate, transition at 730 days):

.. code-block:: python

    >>> r = dca.hyp2_from_eur(1.5e6, 1500.0, 100.0, telf=730)
    >>> round(r.di, 8)
    0.00520362
    >>> dca.hyp2_eur(r.qi, r.di, 100.0, r.telf, r.b, r.b2)
    1500000.0000000002


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
     - Arps b-factor. b > 1 permitted (EUR remains finite for any q_min > 0); see mh_eur for the two-segment form usually used with transient-flow b-factors
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
     - 'exponential', 'harmonic', 'hyperbolic', 'duong', 'mh' (modified hyperbolic), 'hyp2' (two-segment hyperbolic), or 'best' (default). 'best' tries exponential, harmonic, hyperbolic and Duong and returns the highest R-squared; 'mh' and 'hyp2' are explicit opt-in only, since their extra parameters would dominate the R-squared comparison
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
     - Decline model name ('exponential', 'harmonic', 'hyperbolic', 'mh', 'hyp2', or 'duong')
   * - qi
     - float
     - Initial rate (volume/time)
   * - di
     - float
     - Initial decline rate (1/time). 0 for Duong
   * - b
     - float
     - Arps b-factor (first-segment b1 for 'hyp2'). 0 for exponential, 1 for harmonic. 0 for Duong
   * - a
     - float
     - Duong 'a' parameter. 0 for Arps models
   * - m
     - float
     - Duong 'm' parameter. 0 for Arps models
   * - dterm
     - float
     - Terminal decline rate (1/time) for 'mh' and 'hyp2'. 0 for other models
   * - b2
     - float
     - Second-segment b-factor for 'hyp2'. 0 for other models
   * - telf
     - float
     - Transition time for 'hyp2'. 0 for other models
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

Modified hyperbolic fitting example (monthly rates over 8 years from an unconventional well with transient-flow b = 1.6 and a 0.0006/day terminal decline):

.. code-block:: python

    >>> t = np.arange(15, 2900, 30.0)
    >>> q = dca.mh_rate(1200, 0.008, t, b=1.6, dterm=0.0006)
    >>> result = dca.fit_decline(t, q, method='mh')
    >>> result.method
    'mh'
    >>> round(result.qi, 2), round(result.di, 5), round(result.b, 3), round(result.dterm, 5)
    (1200.0, 0.008, 1.6, 0.0006)
    >>> result.r_squared
    1.0

Two-segment hyperbolic fitting example (transient b1 = 1.8 transitioning to boundary-dominated b2 = 0.5 at 730 days):

.. code-block:: python

    >>> q = dca.hyp2_rate(1500.0, 0.0052, t, 730.0, b1=1.8, b2=0.5)
    >>> result = dca.fit_decline(t, q, method='hyp2')
    >>> result.method
    'hyp2'
    >>> round(result.qi, 1), round(result.di, 5), round(result.b, 3), round(result.b2, 3), round(result.telf, 1)
    (1500.0, 0.0052, 1.8, 0.5, 730.0)
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
     - 'exponential', 'harmonic', 'hyperbolic', or 'best' (default). Passing 'duong' raises ValueError; 'mh' is not supported in the cumulative domain
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
    18358.300027522022

    >>> fc2 = dca.forecast(dr, t_end=50, dt=1.0, uptime=0.8)
    >>> fc2.eur
    14686.640022017618

    >>> rr = dca.RatioResult(method='linear', a=0.5, b=0.001, domain='cum')
    >>> fc3 = dca.forecast(dr, t_end=50, dt=1.0, ratios={'GOR': rr})
    >>> fc3.secondary['GOR']['ratio'][0]
    1.4754115099857197


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
