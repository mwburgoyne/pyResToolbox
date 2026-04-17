Changelist in 3.3.0:

- **BNS Z-factor / critical-property coupling**: When either ``zmethod`` or ``cmethod`` is ``BNS``, both are now forced to ``BNS`` for thermodynamic consistency with a ``UserWarning`` naming the overruled counterpart. ``h2 > 0`` continues to auto-select BNS silently. Non-BNS methods (e.g. ``DAK`` + ``SUT``) remain freely mixable. Implemented via a single ``_resolve_methods`` helper applied at all gas public entry points plus ``GasPVT.__init__``.

- **User-supplied ``tc`` / ``pc`` honored across every method** (including BNS). Semantics:

  - ``SUT`` / ``PMC``: supplied values replace the *mixture* pseudo-critical Tc/Pc.
  - ``BNS``: supplied values replace only the *inert-free hydrocarbon* pseudo-component Tc/Pc. Inert Tc/Pc (CO2, H2S, N2, H2) remain at BNS internal constants, and the BNS 5-component PR-EOS mixes them with the user-supplied HC Tc/Pc.

  This is reflected in ``gas_tc_pc``, ``gas_z``, ``gas_ug``, ``gas_cg``, ``gas_bg``, ``gas_den``, ``gas_dmp``, ``gas_ponz2p``, ``gas_grad2sg``, and ``GasPVT`` (new ``tc=0, pc=0`` kwargs on the class constructor).

- **Rust parity for user Tc/Pc**: The Rust batch paths (``dak_zfactor_batch``, ``hy_zfactor_batch``, ``bns_zfactor_batch``, ``gas_ug_lbc``/``gas_ug_lbc_batch``, ``gas_dmp_rust``, ``gas_ponz2p_rust``) now accept optional ``tc_user`` / ``pc_user`` parameters and apply the same override semantics as Python. Removes the previous Python fallback for BNS+user-Tc/Pc; Rust path is now always exercised.

- 716 validation tests (up from 701 in 3.2.0). New coverage for BNS coupling warnings, ``GasPVT`` user Tc/Pc, and Rust-vs-Python parity for user Tc/Pc.

Changelist in 3.2.0:

- **DCA bug fixes**:

  - ``dca.duong_cum`` — fixed linspace-inversion bug where cumulative volume integrated over a descending axis for ``t < 0.001`` (small-t inputs) and returned a negative value. Lower bound of the integration is now ``min(0.001, t * 0.001)`` so the grid is always ascending.
  - ``dca.forecast`` — now validates ``dt > 0``, ``t_end > 0`` and ``uptime`` in ``(0, 1]`` at entry with clear ``ValueError`` messages. Previously ``dt = 0`` raised an opaque ``ZeroDivisionError`` and out-of-range ``uptime`` was silently accepted.

- **Rachford-Rice solver consolidation**: ``simtools.rr_solver`` is now a thin wrapper that delegates to the canonical ``pyrestoolbox.brine._lib_vle_engine.rr_solver`` (Nielsen & Lia 2022). Removes ~100 lines of duplicated iterative code plus the dead ``ensure_numpy_array`` helper. Inputs validated for length/sum before delegation; behaviour is preserved (``EPS_T=1e-15``, ``max_iter=100``).

- **Rust Sechenov fallback guard**: ``src/vle/mod.rs::flash_tp_rust`` no longer silently substitutes its S&W Eq 8 ``ks`` fallback when a caller passes all-ones ``gamma`` with ``salinity > 0``. The caller-supplied ``gamma`` is now always trusted, eliminating any risk that Python ``framework='proposed'`` calls bypass the specialised Dubessy/Akinfiev/Li/Mao-Duan/Duan-Sun ``ks`` models. ``calc_equilibrium_rust`` retains the S&W Eq 8 path but now carries a prominent doc warning. Python path unchanged (Python always passes the correct ``gamma``).

- **Convergence flag on ``CO2_Brine_Mixture``**: New ``.converged`` attribute on the class. ``True`` after a successful Spycher-Pruess fugacity iteration, ``False`` when the 100-iteration limit is hit (matching the existing ``RuntimeWarning``). Lets downstream callers detect non-convergence programmatically.

- **``sensitivity.tornado`` robustness**: Raises ``ValueError`` if ``base_result`` is not finite (NaN/Inf), and if any ``ranges[param]`` has ``lo > hi``. Previously returned ``nan`` or ``inf`` sensitivities that silently corrupted tornado plots.

- **``layer`` module dedup**: Five copies of the EXP/LANG dispatch (B-clamp, flow-fraction evaluation) consolidated into three private helpers (``_clamp_b``, ``_b_max``, ``_flow_fraction_at_x``). No behaviour change.

- **``recommend`` module docstrings**: ``sg`` on ``recommend_gas_methods`` and ``well_type`` on ``recommend_vlp_method`` are currently unused by the decision logic. Docstrings now flag them as reserved for future rules. Signatures preserved for backward compatibility.

- 701 validation tests (up from 696 in 3.1.5).

Changelist in 3.1.5:

- **Agent-friendly UX**:

  - ``validate_methods`` invalid-method errors now list valid options (e.g. ``Invalid zmethod: 'NOSUCH'. Valid options: ['DAK', 'HY', 'WYW', 'BNS', 'BUR']``). New ``validate_choice`` helper used by all ``nodal`` public entry points (``fbhp``, ``outflow_curve``, ``ipr_curve``, ``operating_point``) to validate ``well_type``.
  - ``simtools.zip_check_sim_deck`` and ``simtools.ix_extract_problem_cells`` accept a ``non_interactive=True`` kwarg that raises a ``ValueError`` instead of prompting on ``input()``. Safe for scripts and agents without stdin.
  - ``simtools.make_vfpinj``/``make_vfpprod`` BHP-failure warnings now use ``warnings.warn`` instead of ``print``.
  - ``DeclineResult``, ``ForecastResult``, ``RatioResult`` ``__repr__`` now summarise array fields as ``ndarray(shape=..., dtype=...)`` so printing a result doesn't flood agent transcripts.
  - ``nodal`` unit-validation errors (``Reservoir.__init__``, ``WellSegment.__init__``) echo the user's original value and unit (e.g. ``got -1 m``) rather than the post-conversion internal number.

- **Release-blocking numerical fixes**:

  - **Garcia CO2-brine density**: Algebraic reformulation of Eq 18 (``brine.garciaDensity`` and ``SoreideWhitson._calc_properties``) removes the ``xCO2 → 1`` singularity. Finite rho at ``xCO2 = 1`` equals ``MwG / vPhi``. Mathematically identical to the old formula for ``xCO2 < 1`` (regression: 1e-12).
  - **``oil.Rs_velarde`` at atmospheric ``pb``**: Now returns ``0.0`` when ``pb <= psc`` instead of emitting NaN from a ``0/0`` division.
  - **``oil.sg_evolved_gas`` silent NaN**: Now calls ``validate_pe_inputs`` at entry; zero pressure or zero ``sg_sp`` raises ``ValueError`` instead of returning NaN.

- 696 validation tests (up from 691 in 3.1.4).

Changelist in 3.1.4:

- **Tier 4 brine improvements**: Adaptive VLE damping in ``flash_tp`` and ``calc_water_content_with_kij`` (replaces fixed 0.7/0.9 factors), ``V2_inf`` cached via ``functools.lru_cache(256)`` in Plyasunov model, ``build_kij_matrix`` cached per ``(T_K, mode)`` on VLE engine instance. Rust VLE flash updated with matching adaptive damping. All three brine models (``CH4_Brine``, ``CO2_Brine_Mixture``, ``SoreideWhitson``) now accept ``p``/``degf``/``wt`` and ``pres``/``temp``/``ppm`` parameter aliases.
- **Oil module split**: Monolithic ``oil.py`` (2100 lines) split into 10 sub-files: ``_constants``, ``_utils``, ``_correlations``, ``_density``, ``_compressibility``, ``_rate``, ``_separator``, ``_harmonize``, ``_tables``, ``_pvt_class``. All public API preserved via ``oil/__init__.py`` re-exports.
- **Named constants**: ~100 correlation coefficient groups extracted to module-level named constants with paper citations across gas, oil, nodal, and brine modules. Replaces ~290 inline magic numbers.
- **VLP scaffold deduplication**: Extracted shared ``_segment_march()`` from 4 gas VLP methods, reducing ~450 lines of duplicate segment-loop boilerplate. Each method now provides a gradient callback only.
- **``@rust_accelerated`` decorator**: New decorator in ``_accelerator.py`` for simple Python/Rust dispatch. Applied to 8 nodal VLP functions where signatures match directly.
- **OilPVT array support**: ``rs()``, ``bo()``, ``density()``, and ``viscosity()`` methods on ``OilPVT`` now accept scalar, list, or numpy array pressure inputs and return matching output shapes.
- **Correlation validity warnings**: Optional warnings for out-of-range inputs on DAK/HY/WYW Z-factor (Tr/Ppr bounds), Standing/VALMC bubble point (API/T ranges), and Beggs-Robinson viscosity (API/T/Rs ranges).
- **``validate_pe_inputs``**: Extended to all oil public functions, gas rate/hydrate, brine init, and nodal ``fbhp()``.
- **Enum aliasing fix**: ``z_method.BNS`` is now canonical (value 3), ``z_method.BUR`` is alias. ``@unique`` decorator applied. Dispatch dicts updated.
- **Error-case tests**: 61 new ``pytest.raises`` tests across all modules covering invalid inputs, out-of-range parameters, and edge cases.
- **Rust Hagedorn-Brown fix**: Fixed temperature discretization and low-rate threshold in Rust HB VLP implementation.
- 691 validation tests (up from 630 in 3.1.3).

Changelist in 3.1.3:

- **Z-factor Z*=Z-B reformulation**: BNS Peng-Robinson cubic solver now uses Aaron Zick's Z* = Z - B variable substitution, which bounds all physical roots to (0, 1] and guarantees Halley solver convergence. Eliminates rare failures at extreme pressures. Falls back to Cardano analytically if needed.
- **Rust batch vectorization**: Rust-accelerated ``gas_z()`` and ``gas_ug()`` now use batch dispatch (single FFI call for all pressures) instead of per-pressure scalar loops. Precomputes critical properties, BIPs, and LBC mixture parameters once per call. New Rust functions: ``bns_zfactor_batch``, ``dak_zfactor_batch``, ``hy_zfactor_batch``, ``gas_ug_lge_batch``, ``gas_ug_lbc_batch``. LBC viscosity batch is 9.3x faster than scalar; full BNS pipeline (Z + viscosity) is 2x faster than pure Python.
- **Rust GWR inverse Laplace transform**: Rust/MPFR-accelerated Gaver-Wynn-Rho algorithm and Bessel functions for influence table generation. Switched dependency from ``gwr_inversion`` to ``ilt-inversion``.
- **oil_co() Pb discontinuity fix**: Fixed oil compressibility discontinuity at bubble point in black oil tables. Added ``undersaturated_only`` mode.
- 630 validation tests (up from 603 in 3.0.5).

Changelist in 3.0.5:

- **Rust acceleration (optional)**: Computationally intensive algorithms now have optional Rust-compiled acceleration. When the compiled extension is present it is used automatically with no API changes; when unavailable, all functions fall back silently to pure Python. Accelerated functions include all 8 VLP segment loops (4 methods x gas/oil), gas Z-factor (DAK, HY, BNS), gas viscosity (LGE, LBC), gas pseudopressure, oil density (SWMH), oil FVF (McCain), DCA hyperbolic grid search (``fit_decline``, ``fit_decline_cum``), and the material balance regression objective. Set ``PYRESTOOLBOX_NO_RUST=1`` to force pure Python; set ``PYRESTOOLBOX_RETRY_RUST=1`` to retry after a blocked probe. Use ``from pyrestoolbox._accelerator import get_status`` for programmatic status checks.
- **brine_props()**: Compressibility return changed from a scalar (saturated only) to a ``[cw_usat, cw_sat]`` list. The undersaturated value (Spivey Eq 4.32) is the isothermal compressibility at constant dissolved gas content. The saturated value (Spivey Eq 4.35) is a pseudo-compressibility of the brine and differentially evolved gas system. Previously only the saturated value was returned.
- **oil_co()**: Changed from saturated pseudo-compressibility (``Co = -1/Bo * (dBo/dp - Bg * dRs/dp)``) to undersaturated compressibility (``Co = -1/Bo * dBo/dp`` at constant Rs). Rs is held at the equilibrium value for the specified pressure, yielding the isothermal liquid-phase compressibility without mixing in differentially evolved gas volume. Values below Pb are now smaller and physically consistent with above-Pb values.
- **oil_co() co_sat parameter**: New ``co_sat=False`` parameter. When ``True``, returns ``[co_usat, co_sat]`` list. Saturated compressibility uses Perrine's definition: ``co_sat = -(1/Bo)*dBo/dp + (Bg/Bo)*dRs/dp``, a pseudo-compressibility including gas evolution effects. Above Pb, both values are equal. Backward compatible — default returns a float.
- **oil_bt()**: New function returning total two-phase oil FVF: ``Bt = Bo + (Rsi - Rs) * Bg``. Above Pb returns Bo. Useful for material balance and reservoir voidage calculations.
- **oil_matbal() metric cw/cf fix**: Fixed bug where ``cw`` and ``cf`` in 1/bar were not converted to 1/psi when ``metric=True``, causing the Efw term to be off by ~14.5x. Regression bounds for ``cw``/``cf`` are also now correctly unit-converted.
- **CO2_Brine_Mixture.Cf_sat** and **SoreideWhitson.Cf_sat**: Documentation clarified that these are pseudo-compressibilities representing the average compressibility of the brine and differentially evolved gas system.
- **make_pvtw_table()**: Table now includes both ``Cw_usat`` and ``Cw_sat`` columns. ``cw_ref`` is now a ``[usat, sat]`` list. PVTW keyword export uses undersaturated compressibility.
- **make_bot_og()**: Water compressibility (``cw`` key) now uses undersaturated value from ``brine_props()``.
- 603 validation tests (up from 588 in 3.0.4).

Changelist in 3.0.4:

- **VLP performance**: Eliminated duplicate Z-factor calculations in all 8 VLP method functions. ``_gas_viscosity()`` now accepts a pre-computed Z-factor, avoiding a redundant Hall-Yarborough solve on every segment iteration. Combined with pre-computing Sutton critical properties (Tc/Pc) once per VLP function call instead of recalculating on every segment step. Delivers ~11% speedup on ``operating_point()`` and ``outflow_curve()`` calls.
- **gas_cg() performance**: Batched the two separate ``gas_z()`` calls (at p and p+1) into a single vectorized call with concatenated pressure arrays, delivering ~43% speedup on gas compressibility calculations.
- Removed commented-out dead code (unused Vasquez-Beggs Rs correlation, debug print statement) from oil module.
- **dca module**: New decline curve analysis module with Arps (exponential, hyperbolic, harmonic) and Duong rate/cumulative functions, EUR calculation, model fitting (``fit_decline()`` with 'best' auto-selection), and forecasting (``forecast()``). Returns ``DeclineResult`` and ``ForecastResult`` dataclasses.
- **matbal module**: New material balance module. ``gas_matbal()`` performs P/Z linear regression for OGIP estimation, returning ``GasMatbalResult`` with fitted slope, intercept, R-squared, and Z-factor at initial pressure. ``oil_matbal()`` implements Havlena-Odeh oil material balance for OOIP estimation with drive index decomposition (DDI/SDI/CDI), gas cap support, and full PVT output.
- **recommend module**: New method recommendation engine. ``recommend_gas_methods()`` uses a decision tree based on composition (H2 presence, inert content, CO2/H2S levels) to recommend Z-factor and critical property methods. ``recommend_oil_methods()`` recommends Pb/Rs/Bo correlations. ``recommend_vlp_method()`` filters VLP methods by well deviation. ``recommend_methods()`` combines all three. Returns ``MethodRecommendation`` dataclasses with rationale and alternatives.
- **sensitivity module**: New sensitivity analysis framework. ``sweep()`` varies one parameter across a range collecting results. ``tornado()`` computes tornado-chart sensitivities for multiple parameters, sorted by impact. Returns ``SweepResult`` and ``TornadoResult`` dataclasses.
- **NodalResult**: ``outflow_curve()``, ``ipr_curve()``, and ``operating_point()`` now return ``NodalResult`` (a ``dict`` subclass) supporting both dict-style (``result['rate']``) and attribute-style (``result.rate``) access. Fully backward compatible.
- **fit_decline_cum()**: New function for fitting decline models to rate-vs-cumulative data, eliminating time from Arps equations. Supports exponential, harmonic, and hyperbolic models (Duong excluded — no analytical q-vs-Np form). Returned ``DeclineResult`` parameters are identical to time-domain fits, so results work directly with ``arps_rate()`` and ``forecast()``. Optional ``t_calendar`` parameter infers per-interval uptime fractions by comparing calendar-average rates to fitted capacity rates, populating ``uptime_mean`` and ``uptime_history`` on the result.
- **fit_ratio() / ratio_forecast()**: New functions for fitting secondary phase ratio models (e.g. GOR, WOR) to production data. Four models: linear, exponential, power, and logistic. Returns ``RatioResult`` dataclass with ``domain`` field ('cum' or 'time') controlling how ``forecast()`` evaluates the ratio. ``ratio_forecast()`` evaluates fitted models at arbitrary x values.
- **forecast() extensions**: New ``uptime`` parameter (default 1.0) scales capacity rate to calendar-effective rate. New ``ratios`` parameter accepts a dict of ``RatioResult`` objects for secondary phase forecasting — each ratio is evaluated against cumulative or time per its domain, producing per-phase rate and cumulative arrays in ``ForecastResult.secondary``. Fully backward compatible.
- **RatioResult**: New dataclass for ratio fitting results with method, parameters (a, b, c), domain, R-squared, and residuals.
- **gas_matbal() aquifer support**: New ``Wp``, ``Bw``, and ``We`` parameters for water production and aquifer influx. Cole plot diagnostics (``F/Et`` vs ``Gp``) computed for all cases. When ``We`` is provided, Havlena-Odeh forced-through-origin regression determines OGIP. ``GasMatbalResult`` extended with ``bg``, ``F``, ``Et``, ``cole_F_over_Et``, and ``method`` fields. Fully backward compatible.
- **fit_decline()** / **fit_decline_cum()**: New ``t_start``/``t_end`` and ``Np_start``/``Np_end`` parameters for windowed fitting. Data outside the window is excluded and the window is shifted to start at zero, so the returned ``qi`` represents the rate at the window start.
- **oil_matbal() regression**: New ``regress`` parameter accepts a dict of parameter names and bounds (e.g. ``{'m': (0, 2), 'cf': (1e-6, 10e-6)}``). Optimizes to minimize the coefficient of variation of OOIP estimates across time steps using ``scipy.optimize.minimize`` with L-BFGS-B. Allowed keys: 'm', 'cf', 'cw', 'sw_i'. Results stored in ``OilMatbalResult.regressed``.
- **oil_matbal() tabulated PVT**: New ``pvt_table`` parameter accepts ``{'p': [...], 'Rs': [...], 'Bo': [...], 'Bg': [...]}`` for user-supplied PVT. When provided, ``api`` and ``sg_sp`` are not required. Bubble point and Rs at Pb are inferred from the table when not explicitly specified.
- **gas_matbal() tabulated PVT**: New ``pvt_table`` parameter accepts either ``{'p': [...], 'Z': [...]}`` or ``{'p': [...], 'Bg': [...]}``. Z and Bg are inter-converted internally. Providing both 'Z' and 'Bg' raises ValueError.
- **RANSAC linear regression**: All linear regression in DCA fitting (``fit_decline()``, ``fit_decline_cum()``, ``fit_ratio()``) and gas material balance (``gas_matbal()`` P/Z regression) now uses RANSAC (Random Sample Consensus) with MAD-based outlier detection. Outlier-contaminated data produces robust fits; clean data reproduces ordinary least squares exactly. New ``ransac_linreg()`` utility in ``shared_fns``.
- **Linearized hyperbolic fitting**: ``_fit_hyperbolic()`` and ``_fit_hyperbolic_cum()`` replaced scipy ``curve_fit`` with a grid-search-over-b linearization approach. For each trial b, the Arps equation becomes linear, and RANSAC regression recovers qi and di algebraically. Eliminates local minima from nonlinear optimization.
- **Documentation**: Standardized Returns tables across all 11 RST documentation files. Every function now has a structured Returns table (between Inputs and Examples) documenting return type, attributes/keys, and descriptions. Covers simple scalars, tuples, dicts, dataclass/class results, and arrays. Class Objects sections updated with cross-references to inline Returns tables.
- **Documentation restructuring**: Reorganized RST files for logical flow. ``nodal.rst`` restructured with expanded intro, Getting Started example, Function List moved before API reference, and method/suitability/unit reference sections moved to end. ``simtools.rst`` intro added with grouped Function List (5 domain categories). Module overview paragraphs added to ``gas.rst``, ``oil.rst``, ``matbal.rst``, and ``library.rst``.
- 588 validation tests (up from 580).


Changelist in 3.0.3:

- **gas_hydrate()**: New function for gas hydrate formation prediction and thermodynamic inhibitor calculations. Returns a ``HydrateResult`` dataclass with hydrate formation temperature (HFT), hydrate formation pressure (HFP), subcooling, hydrate window assessment, inhibitor temperature depression, and required inhibitor concentration. Two HFT correlations: Motiee (1991) and Towler & Mokhatab (2005). Inhibitor depression via Østergaard et al. (2005). Supports five inhibitor types (MeOH, MEG, DEG, TEG, EtOH). Full Eclipse METRIC unit support.
- **gas_hydrate() water content**: Computes equilibrium water content of the gas. Uses the SoreideWhitson VLE model when gas composition (``co2``, ``h2s``, ``n2``, ``h2``) is provided, otherwise uses the Danesh correlation. Separate ``p_res``/``degf_res`` parameters allow specifying reservoir conditions for water content (where gas equilibrated with water) independently of the hydrate assessment point.
- **gas_hydrate() inhibitor capping**: ``required_inhibitor_wt_pct`` is now capped at the physical maximum for each inhibitor type (MEOH: 25%, MEG: 70%, DEG: 70%, TEG: 50%, ETOH: 30%). New ``max_inhibitor_wt_pct`` and ``inhibitor_underdosed`` fields indicate whether the required concentration exceeds the achievable maximum.
- **gas_hydrate() water balance**: Full water balance between reservoir and operating conditions. Reports vaporized water at both reservoir P,T (``water_vaporized_res``) and operating P,T (``water_vaporized_op``), condensed water (``water_condensed``), free water (``free_water``), and total liquid water (``total_liquid_water``). Only liquid water (condensed + free) needs inhibitor treatment.
- **gas_hydrate() injection rate**: New ``inhibitor_mass_rate`` (lb/MMscf | kg/sm3) and ``inhibitor_vol_rate`` (gal/MMscf | L/sm3) fields compute the required inhibitor injection rate from total liquid water (condensed + ``additional_water``) and required concentration.
- **gas_hydrate() default method**: Changed default ``hydmethod`` from ``'MOTIEE'`` to ``'TOWLER'`` (Towler & Mokhatab 2005).
- **hyd_method enum**: New enum for hydrate formation correlation selection (MOTIEE, TOWLER).
- **inhibitor enum**: New enum for thermodynamic hydrate inhibitor selection (MEOH, MEG, DEG, TEG, ETOH).
- New unit conversion constants in ``constants``: ``LB_TO_KG``, ``GAL_TO_LITER``, ``LB_PER_MMSCF_TO_KG_PER_SM3``, ``GAL_PER_MMSCF_TO_L_PER_SM3`` and their inverses.
- **oil_co() bugfix**: Fixed three bugs in oil compressibility numerical derivative. (1) Derivative stencil could cross the bubble point, mixing saturated and undersaturated physics and producing negative Co. (2) Near-Pb one-sided derivative used half the step span of the symmetric derivative without normalization, causing a 2x discontinuity. (3) The derivative formula omitted division by the pressure step size, making results proportional to dp rather than true 1/psi derivatives. All three bugs are now fixed; Co is positive at all pressures with smooth transitions across Pb.
- **_build_bot_tables() bugfix**: Fixed incorrect ``rsb`` parameter passed to ``oil_co()`` during black oil table generation. Was passing current pressure Rs (``rss[-1]``) instead of the actual bubble-point Rs, causing oil_co to compute derivatives around the wrong bubble point at sub-Pb pressures.
- Fixed shadowed duplicate ``test_doc_make_pvtw_table`` in test suite (renamed to ``test_doc_make_pvtw_table_keys``).
- Fixed ``MANIFEST.in`` to reference ``README.rst`` instead of non-existent ``README.md``.
- 416 validation tests (up from 375 in 3.0.1).


Changelist in 3.0.1:

- **OilPVT auto-harmonization**: ``rsb`` is now optional (default 0). When ``degf`` is provided (> 0), the constructor calls ``oil_harmonize()`` internally to resolve consistent Pb, Rsb, rsb_frac, and vis_frac from a single call. Accepts ``uo_target`` and ``p_uo`` for viscosity tuning. ``OilPVT.from_harmonize()`` remains as a deprecated thin wrapper.
- **oil_harmonize()**: New function (replaces ``oil_harmonize_pb_rsb``) that resolves consistent Pb, Rsb, rsb_frac, and ``vis_frac`` viscosity scaling factor from user inputs. Accepts ``uo_target`` and ``p_uo`` parameters to compute vis_frac = uo_target / uo_corr. The deprecated ``oil_harmonize_pb_rsb()`` wrapper remains for backward compatibility, returning the original 3-tuple.
- **OilPVT.vis_frac and OilPVT.rsb_frac**: New ``vis_frac`` and ``rsb_frac`` parameters on ``OilPVT`` constructor (both default 1.0). All ``viscosity()`` outputs are multiplied by vis_frac, and ``rs()`` applies rsb_frac scaling. Both factors flow through to VLP segment calculations (``fbhp()``, ``outflow_curve()``, ``operating_point()``) and ``make_bot_og()`` BOT generation.
- **oil_rate_radial / oil_rate_linear**: Accept ``oil_pvt`` and ``degf`` parameters. When an ``OilPVT`` object is provided, uo, bo, and pb are extracted automatically, and Vogel correction is enabled.
- **gas_rate_radial / gas_rate_linear**: Accept ``gas_pvt`` parameter. When a ``GasPVT`` object is provided, sg, composition, method choices, and pre-computed Tc/Pc are extracted automatically.
- **make_bot_og()**: New ``vis_frac`` parameter (default 1.0) that scales all oil viscosity values in generated black oil tables (PVDO, PVTO). Results dict now includes ``vis_frac`` key.
- **Completion.geometry_at_md()**: New method returning wellbore geometry (TVD, ID, deviation, roughness) at any measured depth along the completion. Supports both oilfield and metric unit systems.
- **Completion.profile()**: New method returning a pandas DataFrame of the wellbore profile at all segment boundaries, including crossover rows where geometry changes. Columns: MD, TVD, Deviation, ID, Roughness.
- Added explicit ``__all__`` exports and module docstrings to all modules for improved discoverability via ``dir()`` and ``help()``.
- ``OilPVT`` docs moved to ``oil.rst``; ``GasPVT`` docs moved to ``gas.rst``. Cross-references remain in ``nodal.rst``.
- 375 validation tests (up from 318 in 3.0.0).


Changelist in 3.0.0:

- **Nodal Analysis module**: Added VLP (Vertical Lift Performance), IPR (Inflow Performance Relationship), and operating point analysis with four multiphase flow correlations (Hagedorn-Brown, Woldesemayat-Ghajar, Gray, Beggs & Brill). Supports gas and oil wells with multi-segment deviated/horizontal wellbores via the new ``WellSegment`` and ``Completion`` classes. PVT convenience wrappers ``GasPVT`` and ``OilPVT`` for use with nodal and VFP functions.
- **VFP Table Generation**: New ``make_vfpprod()`` and ``make_vfpinj()`` functions in simtools for generating Eclipse VFPPROD and VFPINJ lift curve tables directly from wellbore geometry and fluid properties using the nodal VLP correlations.
- **Jerauld Relative Permeability Model**: Added the Jerauld (Arco) two-parameter kr model (``krfamily='JER'``) as a third option alongside Corey and LET in ``rel_perm_table()``.
- **Relative Permeability Curve Fitting**: New ``fit_rel_perm()`` function for fitting Corey, LET, or Jerauld models to measured kr data using least-squares optimization. ``fit_rel_perm_best()`` tries all three models and returns the best fit.
- **LET Physicality Check**: New ``is_let_physical()`` function to verify monotonicity and concavity of LET curves.
- **Simulation Workflow Consolidation**: ``make_bot_og()`` (black oil tables) and ``make_pvtw_table()`` (water PVT tables) now accessible via the ``simtools`` module, consolidating simulation-oriented functions in one place. Original locations (``oil.make_bot_og``, ``brine.make_pvtw_table``) remain as backward-compatible wrappers.
- **Multicomponent Gas-Saturated Brine**: New ``SoreideWhitson`` class for multi-gas brine properties using the Soreide-Whitson VLE framework with Garcia/Plyasunov density corrections and calibrated viscosity corrections. Supports CH4, C2-C4, CO2, H2S, N2, and H2.
- **IAPWS-IF97 Freshwater Density**: All brine models now use IAPWS-IF97 Region 1 for freshwater density base, improving accuracy across temperature and pressure ranges.
- **BNS Z-Factor Improvements**: Fugacity-based root selection for sub-critical conditions, vectorized Halley cubic solver, and tuned LBC viscosity model for hydrogen-containing mixtures.
- **Pseudopressure Performance**: Batch Gauss-Legendre quadrature replacing scipy integration for gas pseudopressure calculations, eliminating scipy dependency from the gas module.
- **Eclipse METRIC unit support**: All public PVT, flow rate, and simulation table functions now accept ``metric=False`` parameter. When ``metric=True``, inputs and outputs use Eclipse METRIC units (barsa, deg C, m, sm3/d, sm3/sm3, kg/m3, 1/bar). Applies to all gas, oil, brine, nodal, and simtools functions including ``GasPVT``, ``OilPVT``, ``WellSegment``, ``Completion``, ``Reservoir``, VFP table generation, and black oil table generation. Standard volumes always reference oilfield standard conditions (60 deg F, 14.696 psia).
- Numerous bugfixes, test suite expansion (318 tests), and code hardening across all modules.


Changelist in 2.2:

- Bugfixes.


Changelist in 2.1.3:

- Updated viscosity parameters for BUR method.


Changelist in 2.1.2:

- Fixed bug in implementation of Velarde, Blasingame & McCain Oil Rs calculation.


Changelist in 2.1.0:

- Fixed variable Typing issue that caused problems with Python 3.9 and older.
- Added reference to the Burgoyne ('BUR') methods for gas Z-Factor and critical property correlation


Changelist in 2.0.0:

- Modified the new Z-Factor method, 'BUR', now a tuned five component Peng Robinson method that is fast and stable and able to handle up to 100% of CO2, H2S, N2 or H2 as well as natural gas. Viscosities are calculated with a tuned LBC model.
- Refactored all code to split into modules for ease of future maintenance

Changelist in 1.4.4:

- Added in new Z-Factor method, 'BUR', which is a tuned five component Peng Robinson method that is fast and stable 

Changelist in 1.4.2:

- Corrected CO2 solubility calculations when two roots in CO2 liquid phase

Changelist in 1.4.1:

- Added calculation of Ezrokhi coefficients for brine density and viscosity with dissolved CO2

Changelist in 1.4.0:

- Introduced CO2 saturated brine calculations using Spycher & Pruess modified SRK EOS method
- Rectified an error introduced in Gas Z-Factor calculations due to errant indentation

Changelist in 1.3.9:

- Tweaks to speed DAK and Hall & Yarborough Z-Factor calculations

Changelist in 1.3.8:

- Fix bug in Hall & Yarborough Z-Factor algorithm

Changelist in 1.3.5:

- Fix bug in ECL deck zip/check recursion


Changelist in 1.3.4:

- Extend ECL deck zip/check function to handle IX formatted decks, and support zipping multiple decks at once.


Changelist in 1.3.2:

- Added robust Rachford Rice solver in Simulation Helpers
- Moved relative permeability functions and simulation helpers to seperate .simtools module


Changelist in v1.2.0:

- Added Component Critical Property Library


Changelist in v1.1.4:

- Attempting to fix reported issue with required dependencies not installing correctly


Changelist in v1.1:

- Fix API to SG calculation (141.4 vs 141.5)
- Added lower limit to first rho_po estimate for Oil Density with McCain method to avoid negative values with high Rs
- Added oil_sg and oil_api functions
- Modified HY Z-Factor solve algorithm to improve robustness
- Modified DAK Z-Factor solve algorithm to improve robustness
- Added Gas Z-Factor correlation from Wang, Ye & Wu (2021)
- Removed 'LIN' Z-Factor method due to significant errors above 12,000 psi. Use WYW method instead if speed needed.