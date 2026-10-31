=================================
Capacitance-Resistance Models
=================================

Capacitance-resistance models (CRM) for waterfloods, after Sayarpour, Zuluaga, Kabir & Lake (2009), The use of capacitance-resistance models for rapid estimation of waterflood performance and optimization, JPSE 69, 227-238 (doi 10.1016/j.petrol.2009.09.006). Injection rates are the input signal and total-liquid production the response; the fit gives each injector-producer connectivity, each producer's time constant and, when producer BHP is supplied, its productivity index. Field applications are in Sayarpour, Kabir & Lake (2009), SPE REE 12(6), 853-864 (doi 10.2118/114983-PA).

All functions are unit-agnostic. Injection and production must be in the same reservoir-volume rate unit (rb/d for both, say), since f is a volume fraction; tau comes back in the time unit of ``t``, and J in rate per pressure unit of ``bhp``.


pyrestoolbox.crm.crm_fit
========================

.. code-block:: python

    crm_fit(t, inj, prod, bhp=None, model='CRMP') -> CRMResult

Producer j obeys dq_j/dt + q_j / tau_j = (1 / tau_j) sum_i f_ij i_i(t) - J_j dpwf_j/dt (Eq. 1). With injection constant and BHP linear over each interval, the exact solution (Eq. 5) is the one-step recursion

.. code-block:: text

    q_j(t_k) = q_j(t_k-1) a + (1 - a) [sum_i f_ij I_i(k) - J_j tau_j dpwf_j(k) / dt_k],   a = exp(-dt_k / tau_j)

``'CRMP'`` fits f_ij, tau_j, q_j(t0) and (with ``bhp``) J_j for every producer. ``'CRMT'`` lumps the field into one tank with f = 1 (Eq. 8) and fits tau and q(t0) only; BHP is not used.

The fit is least squares with f_ij in [0, 1], q0 >= 0, J >= 0 and, for each injector, sum over producers of f_ij <= 1 (Eq. 20). For fixed tau the model is linear in the other parameters, so each producer is solved exactly by bounded least squares inside a one-dimensional tau search (variable projection). When Eq. (20) binds, one Lagrange multiplier per injector couples the producers. Sayarpour et al. minimise absolute error with a generalised reduced gradient solver and several starts; least squares from a global tau grid replaces both. They recommend at least four production points per unknown: CRMP has n_inj + 2 unknowns per producer (n_inj + 3 with BHP).

.. list-table:: Inputs
   :widths: 10 15 40
   :header-rows: 1

   * - Parameter
     - Type
     - Description
   * - t
     - array-like
     - Times, strictly increasing (any unit)
   * - inj
     - array-like
     - Injection rates, shape (n_times,) or (n_times, n_injectors). inj[k] is the average over (t[k-1], t[k]]; inj[0] is unused
   * - prod
     - array-like
     - Total-liquid production rates, shape (n_times,) or (n_times, n_producers), same rate unit as inj
   * - bhp
     - array-like
     - Optional producer bottomhole pressures, same shape as prod (any pressure unit). Adds J per producer
   * - model
     - str
     - 'CRMP' (default) or 'CRMT'

.. list-table:: CRMResult attributes and methods
   :widths: 15 45
   :header-rows: 1

   * - Name
     - Description
   * - f
     - Connectivities, shape (n_injectors, n_producers); f[i, j] is the steady-state fraction of injector i's rate produced by j. CRMT: [[1.0]]
   * - tau, q0
     - Time constant and fitted initial rate per producer
   * - J
     - Productivity index per producer, None without bhp
   * - r_squared, q_fit
     - Per-producer R-squared of the liquid-rate match, and the fitted rates
   * - predict(t, inj, bhp=None)
     - Rates over a schedule starting at the history's t[0] (history plus any future period); bhp is required exactly when the fit used it
   * - pore_volume(ct)
     - Drainage pore volume per producer from tau = ct Vp / J (Eq. 2), in rate x time units; needs J

Example: two injectors and two producers with known connectivities, recovered exactly from noise-free data, then a forecast with injector 2 shut in and injector 1 raised:

.. code-block:: python

    >>> import numpy as np
    >>> from pyrestoolbox import crm
    >>> t = np.arange(0, 720, 30.0)
    >>> inj = np.column_stack([1000 + 400 * np.sin(t / 90), 800 + 300 * np.cos(t / 60)])
    >>> f, tau = np.array([[0.6, 0.3], [0.2, 0.7]]), np.array([45.0, 120.0])
    >>> q = np.zeros((len(t), 2)); q[0] = [600, 500]
    >>> for k in range(1, len(t)):
    ...     a = np.exp(-30 / tau); q[k] = q[k - 1] * a + (1 - a) * (inj[k] @ f)
    >>> r = crm.crm_fit(t, inj, q)
    >>> np.round(r.f, 3).tolist(), np.round(r.tau, 1).tolist()
    ([[0.6, 0.3], [0.2, 0.7]], [45.0, 120.0])
    >>> future = np.vstack([inj, [1500, 0], [1500, 0]])
    >>> np.round(r.predict(np.r_[t, 720, 750], future)[-2:], 1).tolist()
    [[927.2, 774.7], [914.0, 702.9]]


pyrestoolbox.crm.crm_oil_cut_fit
================================

.. code-block:: python

    crm_oil_cut_fit(cwi, fo) -> dict

Fits the power-law oil-cut model f_o = 1 / (1 + alpha CWI^beta) (Eq. 15) by least squares on ln(WOR) = ln(alpha) + beta ln(CWI), using the points with cwi > 0 and 0 < fo < 1. Sayarpour et al. caution that it holds only on the straight-line part of a log-log WOR vs cumulative-injection plot (mature, high-watercut floods), so pass only that section. For a CRMP producer, CWI is its connected cumulative injection, the time integral of sum_i f_ij i_i (Eq. 18); oil rate is then f_o times the CRM liquid rate.

.. list-table:: Inputs and returns
   :widths: 10 15 40
   :header-rows: 1

   * - Name
     - Type
     - Description
   * - cwi
     - array-like
     - Cumulative water injected (any volume unit; alpha depends on it)
   * - fo
     - array-like
     - Oil cut, fraction of total liquid, at each cwi
   * - returns
     - dict
     - 'alpha', 'beta', 'r_squared' (in ln WOR) and 'n_points'

.. code-block:: python

    >>> r = crm.crm_oil_cut_fit([1e5, 3e5, 1e6, 3e6], [0.5, 0.33, 0.2, 0.12])
    >>> round(r['alpha'], 6), round(r['beta'], 4), round(r['r_squared'], 4)
    (0.001249, 0.5833, 0.9987)
