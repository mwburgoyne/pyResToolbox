Changelist in 3.7.3:

- **Gas-free brine viscosity is rebuilt on measured data and is now MULTI-SALT** (``brine.viscosity_route.brine_viscosity``, with ``brine.salt_ratio`` for the ratio alone). The chain is ``mu_water(IAPWS-2008) x salt ratio x pressure factor``. The water leg is IAPWS-2008 (Huber et al. 2009), verified to under 0.034 ppm against that paper's own verification table. The salt ratio is the ion-additive modified Jones-Dole model of Appelo and Parkhurst as implemented in PHREEQC, covering any combination of the ``pitzer.dat`` ions rather than NaCl alone, and reproducing PHREEQC - the reference implementation of the same model - to 0.0015% mean over 145 cases. Kestin's (1981) measured pressure dependence is then grafted on, normalised to 0.1 MPa and evaluated at the brine's ionic strength as an NaCl-equivalent molality; neither salt model carries any pressure term of its own and the measured one is not small, moving the NaCl ratio +3.2% from 0.1 to 35 MPa at 20 degC and 4 molal and reversing sign above about 80 degC. Against Kestin's measured viscosities the chain scores **0.300% mean / 0.78% max on NaCl** and **0.764% / 4.62% on KCl**, over 25-150 degC, 0.1-35 MPa and 0.5-5 molal. **A salinity given with no species named means NaCl**, on every input form (``m=``, ``S=``, ``salts=``, ``composition=``), matching the density side. The previous Mao-Duan chain is retained as ``route='mao_duan'`` and reproduces pre-3.7.3 results exactly; it is NaCl only and raises on any other salt rather than silently treating it as NaCl, which would cost 15.9% mean and 51.2% max on KCl. **Expect existing NaCl viscosities to move by up to 1%** below 5 molal and 150 degC: that is the route change, not a regression. Three limits ship with it and should be quoted. Ion **mixing is not validated against measurement** - every measured brine scored is a single salt, because no measured mixed-brine viscosity was found, so multi-salt rests on reproducing PHREEQC, which is verification of the implementation and not validation of the model. Above 35 MPa the pressure factor is **held, not extrapolated**, so the chain is pressure-blind from 35 to 100 MPa. And the K+ worst case of 4.6% sits at 25-50 degC and 4-5 molal, the same cold concentrated corner where Na+ is worst at 0.78%. The ``pitzer.dat`` per-ion parameters are embedded in ``brine/_lib_pitzer_params.py`` rather than read from disk, since a pip install has no PHREEQC; that table is generated and verified bit-identical to a live parse. Developed in the companion brine_props study, where the ported chain is checked bit-identical to the reference implementation over 2332 cases.

- **Salt effect on the dissolved-gas apparent molar volume is now applied** (``brine.vphi_route.salt_shift``). ``dV = -0.5914*m/(1 + 0.0416*m)`` cm3/mol, giving -0.57 at 1 mol/kg and -2.45 at 5 mol/kg. It is fixed entirely from Tiepel and Gubbins (1972) dilatometry, with no parameter fitted to any brine-density data, and Enns, O'Sullivan and Smith, and Heusler and Gaiser agree on sign and magnitude. An earlier assessment withheld it because the two available brine-density inversions appeared to need opposite signs; that was an artefact of comparing absolute ``V_phi`` error, which mixes each data set's absolute offset with its salinity slope. On the slope alone all four sources agree in sign (-0.6 per molal from direct measurement, -1.65 from Calabrese, -0.33 and -1.71 from the two Yan intervals), and the adopted value is the conservative end of that range. Freshwater results are unchanged exactly; gas-saturated brine densities shift by up to a few hundredths of a percent, and the ``brine.rst`` examples move accordingly. Pass ``salt_effect=False`` to ``brine.vphi_route.V_phi`` to recover the previous behaviour.

- **C3H8 gains its own volume shift** (``brine.pr_vphi``), so it no longer falls back to the Plyasunov correlation. ``VSHIFT['C3H8'] = -0.112963`` gives ``V_phi`` = 72.85 cm3/mol at 298 K and 20 MPa. This entry is weaker than the other six and is labelled as such: no densimetry for propane in water exists, so it is not fitted to a data set but set from the mean of the only two direct 298 K determinations, Moore (1982) 70.7 and Zhou and Battino (2001) 75.0 cm3/mol, whose implied shifts are -0.074763 and -0.151164. Those two disagree by 6.1%, which is the real uncertainty on this gas. It is adopted because the fallback it replaces gave 66.99, below both of them. C3H8 joins ``pr_vphi.UNCALIBRATED_IN_T`` for a stronger reason than N2, H2 and C2H6: it has no temperature-resolved data at all, so its temperature behaviour is the equation of state unaided. **NC4H10 also gains a shift**, ``+0.110924``, from Moore (1982) Table I, 76.6 cm3/mol, the only direct measurement of that quantity. Two things about it are unusual and are flagged in the source: the shift is positive, alone among the eight, and it breaks Moore's own homologous series, whose CH2 increments run +18.4, +17.8 and then +5.9. An earlier note in this changelist said butane was absent from the component set; that was wrong, and was a naming mismatch between ``NC4H10`` here and ``nC4H10`` in the S&W tables.

Changelist in 3.7.2:

- **Dissolved-gas apparent molar volume now comes from the Soreide-Whitson modified PR route plus a volume shift, by default** (``brine.pr_vphi``, exposed as ``brine.V_phi_pr``; selected via ``SoreideWhitson(vphi_route=...)``). ``V2_inf`` is obtained from the exact infinite-dilution relation ``Vbar_2 = -(dP/dn2)_{T,V,n1} / (dP/dV)_{T,n}`` evaluated on the water-rich liquid root of the S&W binary, so no flash is needed and the dissolved gas never has a phase root of its own - it enters only through its co-volume and the cross term ``a12 = sqrt(a1 a2)(1 - kij_aq)``. A single dimensionless Peneloux volume shift per gas (``s = c/b``, the standard VSHIFT convention) closes the remaining offset: CH4 -0.109632, CO2 -0.037913, H2S -0.078975, N2 -0.155288, H2 -0.177625, C2H6 -0.073142. Because Peneloux translation gives ``V2_inf(shifted) = V2_inf(EOS) - s*b`` exactly, the shift is a pure offset with no temperature or pressure shape, and a temperature-dependent shift was tested and rejected as overfitting (on 21 held-out Murphy & Gaines H2S points a constant shift scored 0.6% mean absolute error against 1.3% for ``s(T)``). The shifts are calibrated to **direct volumetric measurement only** - Hnedkovsky (1996) vibrating-tube densimetry, Moore (1982) densimetry, Tiepel and Gubbins (1972) dilatometry, Bignell (1984), Murphy and Gaines (1974) float densities, Zhou and Battino (2001) - excluding determinations that inferred a volume rather than measuring one (Enns, O'Sullivan, Heusler) and excluding V_phi inverted from brine density (Calabrese, Yan), which scatters five to seven times wider. This route is now the **default**, with Plyasunov retained as the reference and the fallback: it needs one number per gas against Plyasunov's 35 eight-significant-digit coefficients per gas plus IAPWS-IF97, it matches or beats Plyasunov on five of six gases against the calibration densimetry (H2S 0.5% against 1.4% mean absolute error, CH4 1.4 against 1.5, N2 3.3 against 4.9, H2 5.1 against 6.2, C2H6 3.0 against 3.4; CO2 0.9 against 0.7), and it reproduces the H2S temperature trend that the Plyasunov correlation misses (Murphy and Gaines measured +0.69 cm3/mol between their cold and warm groups; this route gives +0.53, Plyasunov +0.02). Validity is guarded to 273.15-623.15 K (the IAPWS-IF97 Region 1 ceiling), 100 MPa, above the water saturation pressure, and the six calibrated gases. Three ceilings should not be conflated: ``T_MAX`` 623.15 K is where the arithmetic stops, ``T_CALIBRATED_MAX`` 473.15 K is where the fitted shift stops being fitted, and ``T_VOUCHED_MAX`` 450 K (about 350 degF) is the only one that is an accuracy claim. Note that N2, H2 and C2H6 have no temperature-resolved calibration data, so their temperature behaviour is the EOS unaided (flagged in ``pr_vphi.UNCALIBRATED_IN_T``). Developed and validated in the companion brine_props study.

- **``SoreideWhitson`` gains a ``vphi_route`` argument** (``'auto'`` default, ``'pr'``, ``'plyasunov'``) selecting the source of the dissolved-gas molar volume, with the dispatch policy in ``brine.vphi_route``. ``'auto'`` uses the PR route wherever it is calibrated and falls back to Plyasunov only for C3H8 and NC4H10, which have no fitted shift, and where the quantity is undefined (outside 273.15-623.15 K, above 100 MPa, below the water saturation pressure). There is deliberately **no temperature handover**: an intermediate design switched routes at 473 K, which introduced a discontinuity in the middle of the delivered answer for no gain, and inside the range where accuracy is claimed the two routes agree on gas-saturated brine density to 0.12% at worst. ``vphi_route='plyasunov'`` reproduces pre-3.7.2 densities exactly. Gas-saturated brine densities shift by roughly 0.01-0.02% under the new default (the ``brine.rst`` pure-CO2 field example moves from 0.973227 to 0.973112 sg), and ``Cf_sat`` shifts about 1.2% because it is a difference of two densities.

- **H2S dissolved-gas viscosity coefficient corrected from 1.64 to 1.79.** Murphy and Gaines (1974) report viscosity ratios but no dissolved mole fractions, so the Islam-Carlson coefficient has to be calibrated against an assumed solubility. It was previously calibrated using this library's own Soreide-Whitson flash, which gives ``x_H2S`` 4.5 to 9.2% higher than Burgess and Germann (1969) at the five measurement conditions - and Burgess and Germann is the source Murphy and Gaines state they used. Since the coefficient scales as the reciprocal of ``x``, the mismatch biased it low by the same factor. Re-basing on Burgess and Germann Table 5 gives implied coefficients of 2.31, 1.54, 1.98, 1.31 and 0.27 across their five points; the adopted value remains the mean of the four below 35 degC, which is now 1.79. The H2S viscosity factor at ``x = 0.03`` moves from 1.0469 to 1.0512. One caveat: Burgess and Germann's table starts at 30 degC, so the three 28.1 degC points sit just below its lowest row, Murphy and Gaines having deliberately worked close to the hydrate boundary.

- **Garcia mixing-rule solvent basis corrected in both brine paths.** The apparent molar volume of a dissolved gas is defined against the whole gas-free brine, ``V_phi = (V_solution - V_brine)/n_gas``, so ``n_gas*V_phi`` must be added to the volume of water *and* salt. ``SoreideWhitson`` paired ``MWWAT`` with a brine density and a salt-free mole fraction, which omitted the salt from both the solvent mass and the solvent volume and inflated the whole dissolved-gas density effect by a factor approaching ``1/(1-S)``: +5.6% at 1 mol/kg NaCl, +11% at 2 mol/kg, +29% at 5 mol/kg. ``CO2_Brine_Mixture.garciaDensity`` had the mirror-image error, combining the paired-NaCl brine mole weight with a fully-ionised mole fraction, which overstated the solvent mass by 8% at 5 molal and understated the CO2 densification to match. Both now use the unambiguous per-kilogram-of-water form ``rho = (Wb + m_gas*M_gas)/(Wb/rho_brine + m_gas*V_phi)`` with ``Wb = 1000 + m_NaCl*M_NaCl``, so the two paths finally agree on a basis (they differed by 7.5% of the gas effect at 1 molal and 39% at 5 molal). Freshwater results are unchanged to machine precision; gas-saturated brine densities shift with salinity (the ``brine.rst`` 3% NaCl examples move from 0.973377 to 0.973227 sg for CO2 and 0.964114 to 0.964234 for CH4). Undersaturated compressibility uses the same basis at P+1.

- **CH4 dissolved-gas viscosity rebuilt as a single saturating expression.** The correction was a linear-in-x ramp spliced to a temperature-only plateau by a ``min()``, four parameters in total, built from the three summary plateau ratios quoted in the text of Ostermann (SPE 14211). Reading that paper's Tables 1a-1c (23 measurements at 100/150/250 degF over 500-7000 psi, each with Rsw tabulated) overturned both premises of that construction: the plateau is not flat, drifting up 1.0-1.6 percentage points between 2000 and 7000 psi, and Rsw is tabulated so the dissolved amount need not come from an EOS. All 23 points are now fitted by one expression, ``mu_sat/mu_free = 1 + A*exp(B/T_K)*x/(K + x)`` with A = 1.71739196e-3, B = 1239.77535 K and K = 1.52860547e-3: three parameters, no cap, Arrhenius in temperature (activated flow, hydration structure destroyed thermally; B corresponds to 10.3 kJ/mol, about half a water hydrogen bond) and Langmuir in concentration (finite structurable capacity). It is linear in x in the dilute limit as required and saturates naturally in both variables. Fit quality on the 23 measurements: 0.72 percentage points RMS with three parameters, against 0.71 for the superseded four-parameter spliced form and 1.15 for an unsaturated linear form, so the saturation is real and nothing is lost by dropping the switch. The quadratic in degF printed in the paper is not used anywhere: three coefficients on three summary values is an exact interpolation with no degrees of freedom, and it turns upward past its 273.5 degF vertex. Gas-saturated brine viscosities with dissolved CH4 change modestly, most at low dissolved amounts where the previous ramp shape was assumed rather than fitted.

- **IAPWS-IF97 saturation-line guard.** ``plyasunov.iapws_if97`` gains ``p_sat_if97(T)`` (Region 4 basic equation, reproducing the IF97 Table 35 check values to 2e-9 relative), and ``V2_inf`` now raises when the pressure is below the water saturation pressure, where the solvent is vapour and the infinite-dilution partial molar volume is undefined. The Region 1 density and compressibility functions deliberately do **not** enforce the saturation line, because they are used to obtain an atmospheric-reference brine density at reservoir temperature, which is a legitimate metastable extrapolation.

- **CO2 dissolved-gas viscosity correction moved to Calabrese et al. (2019) Eq. 25** in both ``CO2_Brine_Mixture`` and ``SoreideWhitson``, replacing the temperature-independent Islam-Carlson (2012) factor ``1 + 4.65*x^1.0134``. The Calabrese increment ``ln(mu/mu_brine) = 65.560*exp(-2.468*(T_K/142 - 1))*x_CO2`` is fitted to 415 NaCl-brine points over 274-449 K to 100 MPa (AARD 0.9%) and is independent of salt type and molality (valid to m = 6). Islam-Carlson matches it only near the 25-35 degC window it was calibrated in, and overstates the viscosity increment by 1.7x at 122 degF, 3.5x at 200 degF and 9.4x at 302 degF. Gas-saturated brine viscosities with dissolved CO2 change accordingly: unchanged near 30 degC, materially lower at reservoir temperatures (e.g. at 200 degF and x_CO2 = 0.01 the correction falls from +4.4% to +1.3%).

- **CH4 dissolved-gas viscosity correction gains a sub-plateau ramp** in ``SoreideWhitson``. The Ostermann (SPE 14211/15081) plateau ratios were refitted to the same functional form as the CO2 correction, ``ln(mu/mu_brine) = 93.6918*exp(-0.957408*(T_K/142 - 1))*x_CH4``, using S&W-computed x_CH4 at the 2000 psi plateau onset (slope residuals within 3.4%); the factor is capped at the measured Ostermann plateau. At and above plateau saturation results are unchanged; below ~2000 psi CH4 partial pressure the correction now scales with dissolved CH4 instead of applying the full plateau step (e.g. pure CH4 at 500 psia / 150 degF: +1.3% instead of +4.4%).

Changelist in 3.7.1:

- **Multicomponent flash feed fix in the brine module.** ``calc_gas_brine_equilibrium(method='flash')`` (used by ``SoreideWhitson`` for mixed-gas VLE) flashed a hardcoded 95% water / 5% gas feed. Preferential dissolution of the more-soluble components depleted the equilibrium vapor off the requested dry-gas composition - for a 73/27 CO2/CH4 gas at 50 MPa the dry-basis vapor CO2 slid to about 0.65, biasing dissolved CO2 roughly 10% low and dissolved CH4 high. The ``y_*`` inputs are now interpreted as the dry-basis equilibrium vapor composition: the feed starts gas-excess (50/50) and retries water-rich (0.9/0.95/0.98) near the water boiling curve, where a two-phase solution requires the feed water fraction to lie between the liquid and vapor water mole fractions; both flashes are checked for converged two-phase results. Single-gas (binary) results are unchanged (< 0.04%; the binary tie-line is feed-independent). Validated against ternary CH4+CO2+H2O VLE data (Qin 2008; Al Ghafri 2014; 37 points): dissolved-CO2 MARE improves from 26% to 16%. The ``SoreideWhitson.y`` docstring now states that it equals the input dry-gas composition, which the flash honours by construction. Mixed-gas frozen baselines and doc examples were re-pinned to the corrected values (the ``brine.rst`` mixed-gas example ``Rs_total`` moves from 6.31 to 8.51 sm3/sm3 - the old value had H2S and CO2 depleted by the feed artifact).

Changelist in 3.7.0:

- **Two-segment decline models for unconventional wells** in the ``dca`` module. Both are continuous in rate and nominal decline (log-slope) at the transition, ship with full RST documentation and doc-example tests, and are unit-agnostic like the rest of the module:

  - **Modified hyperbolic** (``mh_rate``, ``mh_cum``, ``mh_eur``): hyperbolic decline with b-factor ``b`` (default 2.0, transient flow) until the nominal decline falls to a terminal decline ``dterm``, exponential thereafter (Robertson 1988, SPE-18731-MS). Analytic piecewise cumulative and EUR.
  - **Two-segment hyperbolic** (``hyp2_rate``, ``hyp2_cum``, ``hyp2_eur``): first-segment ``b1`` (default 2.0) to a specified transition time ``telf``, then a second hyperbolic segment ``b2`` (default 0.5) anchored so rate and slope match the first segment at ``telf`` - the sharp-transition form of the transient hyperbolic model (Fulford and Blasingame 2013, SPE-167242-MS). An optional ``dterm`` gives the second segment an exponential tail.
  - **EUR-constrained type-curve generator** (``hyp2_from_eur``): solves the initial nominal decline so a two-segment hyperbolic delivers a target EUR at a given abandonment rate, with the transition specified by exactly one of fixed time, switch decline rate, EUR fraction or rate fraction. Ports the Burgoyne (2016-2018) two-stage generator logic, with the original bisection replaced by a bracketed brentq root find. Three modes reproduce the reference workbook within its own convergence tolerance; the switch-decline mode differs by about 0.2% because the workbook's VBA froze the transition time at its initial guess (delivering 8.811 BCF against an 8.8 BCF target in the reference case), whereas the port keeps the transition self-consistent and returns the target EUR exactly.

- **Arps b-factors above 1 are now accepted** by ``arps_rate``, ``arps_cum`` and ``eur`` (super-hyperbolic/transient flow; previously raised ``ValueError``). EUR remains finite for any positive economic limit; cumulative is unbounded in time unless paired with a terminal decline.

- **Fitting and forecasting support**: ``fit_decline`` gains explicit ``method='mh'`` and ``method='hyp2'`` fits via bounded curve_fit (excluded from ``'best'`` so existing model selections are unchanged - their extra parameters would otherwise dominate the R-squared comparison). ``DeclineResult`` gains ``dterm``, ``b2`` and ``telf`` fields, and ``forecast()`` dispatches both new models with analytic cumulatives.

- **CI**: GitHub Actions runtimes moved to Node 24 (checkout v5, setup-python v6, upload-artifact v6, download-artifact v7, action-gh-release v3).

Changelist in 3.6.0:

- **Correctness release.** A full literature-verified review of every calculation module found and fixed a set of long-standing correlation bugs. Several outputs change materially - notably Motiee hydrate temperatures, Gray and Beggs-Brill VLP pressures, linear gas rates, stock-tank GOR and undersaturated oil viscosity. Each fix below was confirmed against the original paper (or two independent secondary sources) before the change was made. All frozen regression baselines and documentation examples were re-pinned to the corrected values.

- **Gas Z-factor fixes**:

  - **Hall-Yarborough now actually converges.** The Newton loop exited after exactly one step (a loop-variable bug present since the original implementation), so HY returned "initial guess + one step", in error by up to 14% near the critical region (Tpr 1.05, Ppr 2). The solver now iterates to a 5e-4 relative step with a warning on non-convergence. HY Z values shift accordingly (mid-range shifts are small, e.g. 0.86774 to 0.86775 at 2000 psia / sg 0.75 / 200 degF).
  - **DAK Newton derivative corrected** (the R5 term derivative used ``rhor**3`` where the algebra gives ``rhor**2 - (A11*rhor^2)^2`` terms). Converged Z values shift by < 1e-6; convergence is restored to quadratic (about 5 iterations instead of 18 near-critical). A non-convergence warning was added.
  - **WYW Z-factor method removed from the public API.** The Wang, Ye & Wu (2021) explicit fit behaves incorrectly at low Ppr (Z tending to 0.944 as Ppr tends to 0) and existed mainly as the HY initial guess, which it remains internally. ``zmethod='WYW'`` now raises ``ValueError``.
  - **Single-sided ``tc``/``pc`` overrides honoured on Rust paths.** Effective tc/pc are now resolved once via ``gas_tc_pc`` and forwarded to the Rust kernels, so supplying only one override gives identical answers with and without the accelerator (previously the Rust path silently ignored a lone override). The DAK and HY array paths also now route through the Rust batch kernels for all cmethods, roughly 2x faster than the previous per-element dispatch and bit-identical.

- **``gas_rate_linear`` over-predicted by exactly 2 pi.** The linear-flow branch reused the radial constant 1422, which embeds the 2 pi of radial geometry. A separate linear-flow constant (2 pi x 1422) is now used; linear gas rates drop by a factor of ~6.28 (doc example: 8.20 to 1.31 Mscf/d). Radial rates are unchanged.

- **Motiee hydrate temperature corrected to the published psia/degF basis.** The correlation was evaluated in kPa/degC with a transcribed intercept, reading roughly 30-36 degF hot (1000 psia, sg 0.65: 96.2 degF, vs Towler 61.1 and the Katz chart 62-66). The verified Motiee (1991) form now gives 60.2 degF at the same point. Hydrate formation pressures and inhibitor dose estimates from ``hydmethod='MOTIEE'`` shift substantially; the default TOWLER method was always correct.

- **Oil correlation fixes**:

  - **``oil_rs_st`` returned ln(Rst).** Valko-McCain (2003) Eq 3-2 defines ln(Rst) as the polynomial; the missing exponential is now applied (doc example: 4.18 to 65.13 scf/stb, consistent with the 0.1618 x Rsp heuristic). Confirmed against the original paper, including the psia separator-pressure basis (nomenclature, p168). Metric output is now also converted correctly.
  - **Petrosky-Farshad undersaturated viscosity used natural log where McCain Eq 3.24b specifies log10**, producing unphysical viscosity growth above Pb (6x over 6000 psi in the worst verified case). PVTO undersaturated viscosities and viscosity-dependent doc examples shift.
  - **Vogel IPR edge cases.** The Pb clamp only applied to array ``pwf`` in ``oil_rate_radial`` and never in ``oil_rate_linear``, so scalar undersaturated calls under-predicted by ~22%; saturated reservoirs (pr <= pb) could return negative rates at pwf = pr. Both functions now share an elementwise helper: composite Darcy+Vogel above Pb, pure Vogel driven by pr when saturated (rate is zero at pwf = pr).
  - **``oil_rs_bub`` VALMC inversion is now analytic** (quadratic then cubic root selection) with machine-precision round-trip, replacing a secant loop that diverged silently above the correlation's representable Pb maximum (returning e.g. rsb of 476,000). Requests beyond the maximum now raise ``ValueError``. The docstring's claimed STAN fallback (which never existed) is removed.
  - **``oil_rs`` with ``rsmethod='STAN'`` now honours user ``rsb``/``pb``** (scaled like VALMC), removing an Rs discontinuity at Pb. ``oil_bo`` with ``bomethod='STAN'`` now applies the McCain cofb undersaturated correction above Pb instead of staying constant.
  - **Oil density above Pb** now evaluates the bubble-point anchor with ``rsb`` (matching the Rust path and McCain's formulation), keeping density continuous at Pb; ``oil_co`` uses a one-sided difference near Pb so compressibility stays on a single smooth branch. ``oil_deno`` clamps the temperature term so inputs below 60 degF return real values. **Minor breaking change**: the unused ``pi`` parameter was removed from ``oil_co``.
  - ``import pyrestoolbox.oil`` is now lazy (pandas/tabulate/gas/brine load only when needed), about 4x faster.

- **Brine fixes**:

  - **Salinity-method routing in ``SoreideWhitson``.** The documented ``'sechenov'`` and ``'auto'`` options silently applied no salting-out at all (freshwater answers at any salinity); they now normalise to ``'gamma_phi'``. ``framework='proposed'`` + ``'embedded'`` claimed a gamma_phi fallback that never happened - it now happens. ``framework='dropin'`` + ``'embedded'`` now wires the delta-kij tables into the flash. ``framework='sw_original'`` with the default ``'gamma_phi'`` double-counted salting-out (CH4 solubility ~3x understated at 20 wt%); the Sechenov gamma is now skipped for gases whose kij already embeds salinity. The engine raises ``ValueError`` on unrecognised salinity methods. Regression tests pin x_CH4 at 20 wt% NaCl for every framework/method combination.
  - **Spycher-Pruess 99-109 degC blend corrected** (verified against the 2004 and 2010 papers): the blending weights were inverted (returning the high-temperature model at 99 degC and vice versa), and the Python fugacity blend was a no-op through list aliasing, which was also a real Python/Rust divergence inside the band. Both fixed in both languages; xCO2 is continuous across the band boundaries.
  - **Saturated brine compressibility derivative** (Spivey Eq 4.33 chain) corrected: the Sechenov pressure-derivative was lumped into the denominator instead of subtracted, a 3.5-6% error in saturated cw for saline brines.
  - **Rs mole-basis corrections**: ``CO2_Brine_Mixture`` converted Spycher's ionised-basis xCO2 with paired-NaCl brine moles (Rs ~ -6.7% at 20 wt%); ``SoreideWhitson`` converted the salt-free flash x with brine moles (Rs ~ +7.7% at 20 wt%). Both reduce to the previous behaviour at zero salinity.
  - **Constant consistency**: 273 to 273.15 K in the ``brine_props`` chain, the 0.1015 MPa typos to 0.1013, standard-conditions CH4 molar volume at 288.706 K / 0.101325 MPa, and unified MW and Tb constants that had drifted between modules. Freshwater viscosity shifts ~0.14%; the two brine models now agree on identical inputs.
  - **No more NaN cascade below the water vapour pressure**: ``brine_props`` at low pressure / high temperature now returns finite values with zero dissolved methane. Pressures above the IAPWS-IF97 Region 1 limit (100 MPa = 14,503 psia) raise a clear error naming the limit.
  - Removed ~850 lines of dead code from the VLE engine and salting library (unrouted kij variants, superseded methods, demo blocks).

- **Nodal VLP fixes** (each applied identically in Python and the Rust accelerator):

  - **Gray (1974) holdup rewritten to the published API 14B form.** The previous implementation built its velocity group without surface tension and applied the R/(R+1) transform twice, over-predicting holdup ~30x in mist flow and putting Gray BHPs far out of family (condensate-well example: 1940 to 1067 psia, now beside HB/WG/BB at 1077-1209). The effective (wet film) roughness is also now the published piecewise form.
  - **Beggs-Brill liquid velocity number** used sigma in lbf/ft where the 1.938 constant requires dyne/cm (~11x inflation corrupting the inclination correction; BHP previously fell as water rate rose). **The flow-pattern map** now matches the published revised map (two regions were misclassified as intermittent instead of distributed), and **the Payne et al (1979) 0.924 correction** now applies to all uphill patterns including distributed. BB BHPs shift ~1% in typical cases (doc example: 1213 to 1224 psia).
  - **Hagedorn-Brown two-phase Reynolds constant** corrected from 96778 to 1900.8 (published 2.2e-2 with mass rate in lbm/day). Mostly masked in fully rough turbulence; HB BHPs shift on smooth pipe (doc example: 954.8 to 962.1 psia).
  - **Divergence guard**: an infeasible injection case previously returned a large negative BHP silently (e.g. -239,000 psia). The march now raises ``RuntimeError`` when pressure falls below atmospheric, identically from the Python and Rust paths.
  - The low-rate static-column fallback now warns when liquid rates are being ignored; ``fbhp``/``outflow_curve``/``operating_point`` take ``gsg`` from a supplied ``gas_pvt``; ``fthp``'s ``thp_min`` metric sentinel is fixed (default is now ``None``); negative rates raise ``ValueError``.
  - **Internal**: the Rust VLP port now mirrors the Python ``_segment_march`` scaffold (one shared march, single ``calc_segments``/``condensate_dropout``, ~75 named constants with citations), removing ~1,200 duplicated lines - the structure that bred the 3.4.0/3.5.0 constant-drift bugs. All eight entry points reproduce pure-Python to better than 1e-9.

- **``lorenz_2_layers``** built its saturation grid with ``np.arange``, which appended a spurious extra layer for 28 of nlayers in 2-200 (e.g. 7 layers averaging k=86 when asked for 6 averaging 100). Now ``np.linspace``; a regression sweep pins nlayers 2-50.

- **``oil_matbal`` metric ``Gi``** was understated by exactly 5.6146x (sm3 multiplying an internal rb/scf Bg), inflating OOIP ~34% on gas-injection cases in metric units. Field-unit results were always correct; a field-vs-metric equivalence test now pins this.

- **simtools fixes**: ``rel_perm_table`` now honours ``sgcr`` for SGOF and returns exactly the documented row count in all configurations (two off-by-one cases fixed); ``make_bot_og`` in metric mode now actually converts the Bg column it relabels to rm3/sm3 (~178x), no longer overwrites the harmonised ``sg_sp``, and reports the harmonised ``pb``/``rsb`` in its results dict; ``zip_check_sim_deck`` no-files path raises ``FileNotFoundError`` instead of ``TypeError``; ``zip_check_sim_deck`` and ``ix_extract_problem_cells`` default to ``non_interactive=True`` (**behaviour change** - pass ``False`` for the old prompts). **Internal**: ``simtools.py`` is split into domain sub-files (``_decks``, ``_relperm``, ``_aquifer``, ``_vfp``, ``_pvt_tables``) with the public API unchanged.

- **DCA**: ``forecast()`` cumulative/EUR now uses the analytic cumulative functions (uptime-scaled) instead of a right-rectangle rate sum that biased EUR ~ -1.7% low at dt=1; ``fit_decline(method='best')`` compares R2 on the common t > 0 subset so the comparison is apples-to-apples when t contains 0.

- **Package infrastructure**: ``pyrestoolbox.__version__`` now exists (single-sourced from pyproject.toml) and the Rust extension reports its own version via ``get_status()``; exceptions raised inside Rust-accelerated functions now propagate instead of silently falling back to Python; ``validate_pe_inputs`` rejects NaN and gains a scalar fast path ~40x faster; passing an Enum of the wrong class raises a descriptive ``ValueError`` instead of ``KeyError``; ``bisect_solve`` raises when the root is not bracketed; a new GitHub Actions workflow runs the full suite (pure-Python and Rust-accelerated, including parity tests) on every push and pull request; ``run_all_tests.py`` now delegates to pytest and runs everything.

- 870 validation tests (up from 838 in 3.5.0), including 113 Python/Rust parity tests.

Changelist in 3.5.0:

- **Gas non-Darcy & partial-penetration pseudoskins** - three new public functions in ``pyrestoolbox.gas`` plus two matching ``GasPVT`` convenience methods. Oilfield + metric unit support throughout.

  - ``gas.gas_hvf_beta(k, method='FK' | 'JONES' | 'TCK', phi=0, metric=False)`` - Forchheimer high-velocity-flow coefficient β. Three correlations: ``'FK'`` (default, log-log fit of Firoozabadi & Katz 1979 JPT Feb, consolidated-rock β(k) chart), ``'JONES'`` (Jones 1987 SPE-16949), ``'TCK'`` (Tek, Coats, Katz 1962 JPT Jul, requires ``phi``). Returns β in 1/ft (or 1/m if ``metric=True``).
  - ``gas.gas_non_darcy_skin(qg, k, h_perf, rw, mug, sg, krg=1.0, beta_method='FK', phi=0, metric=False)`` - returns ``{'beta', 'D', 'S_hvf'}`` where D uses the standard ``2.222e-15 · β · γg · k / (μg · h · rw)`` form (Jones 1987; Odeh, Moreland & Schueler 1975 JPT Dec). Output ``D`` is in day/MSCF (or day/sm3 if metric) and can be passed directly to ``gas_rate_radial(D=...)``. ``krg < 1`` evaluates β at the damaged-zone permeability ``k' = k·krg`` for a pessimistic S\ :sub:`hvf`.
  - ``gas.gas_partial_penetration_skin(htot, htop, hbot, rw, kh_kv=10)`` - partial-penetration pseudoskin by direct evaluation of the Streltsova-Adams (1979) SPE-7486 analytical series. Bessel K\ :sub:`0` from ``scipy.special.k0``; vectorised summation to 20,000 terms (sub-millisecond). Warns when the series tail suggests incomplete convergence for very thin wellbores. Unit-agnostic - only ratios enter the formula, so inputs share any consistent length unit.
  - ``GasPVT.non_darcy_skin(qg, p, degf, k, h_perf, rw, krg=1, beta_method='FK', phi=0)`` and ``GasPVT.partial_penetration_skin(htot, htop, hbot, rw, kh_kv=10)`` - the former auto-computes μ\ :sub:`g` from stored PVT state; both honour the ``GasPVT.metric`` flag.

  18 new unit tests in ``test_gas.py`` (all three β correlations, metric round-trip, damaged-zone, series convergence, geometry validation, ``GasPVT`` delegation). 12 new doc-example tests in ``test_doc_examples.py`` pin every RST example value. Reference S\ :sub:`p` values verified against a 50,000-term direct evaluation of the Streltsova-Adams series.

- **PVTO export format fix** (``simtools.make_bot_og(pvto=True, export=True)``). The exported ``PVTO.INC`` placed Eclipse's ``/`` stem-terminator on every saturated row, prematurely closing each stem before its undersaturated continuation rows. A strict parser would then read each continuation row as an unset stem and either error or drop the data. The saturated row now omits ``/`` when undersaturated rows follow; the terminator moves to the LAST undersaturated row of the stem, matching the reference example in the Eclipse manual (``PVTO`` keyword). Stems with no undersaturated extension still carry ``/`` on the saturated row, and the final null-record ``/`` on its own line still terminates the table. Underlying ``usat`` data returned in the results dict was already correct - only the export layout changed.

- **Removed dead module** ``pyrestoolbox/oil/oil.py``. Pre-April-2026 monolith superseded by the ``_correlations.py`` / ``_density.py`` / ``_tables.py`` / etc. sub-modules. Nothing in the package or tests imported from it. No public-API change - ``from pyrestoolbox import oil`` and ``oil.<func>`` keep working via ``oil/__init__.py`` re-exports.

- **Brine framework correctness (Rust path)**. ``SoreideWhitson(framework='dropin')`` and ``framework='sw_original'`` were silently computed as ``'proposed'`` whenever the Rust accelerator was installed: the Rust flash hardcodes the proposed-framework MC-3 water alpha and kij\ :sub:`AQ`, with no framework parameter. The Rust dispatch in ``brine._lib_vle_engine`` (``flash_tp`` and ``calc_equilibrium``) is now gated on ``framework == 'proposed'``; the other two frameworks take the framework-aware Python path. Pure-Python results were always correct. Adds a regression test exercising both non-proposed frameworks.

- **Brine convergence flag (Rust path)**. ``flash_tp_rust`` and ``co2_brine_solubility_rust`` now return the convergence bool they compute, so ``SoreideWhitson``/``CO2_Brine_Mixture`` ``.converged`` (and the ``converged_aq``/``converged_na`` flash flags) report the real outcome instead of being hardcoded ``True`` when Rust is active.

- **Oil ``oil_harmonize`` default ``pbmethod`` unified to VALMC** (was ``VELAR``), matching ``oil_pbub`` / ``oil_rs`` / ``oil_co`` / ``OilPVT``. Previously ``oil_harmonize`` and the other oil functions returned mutually inconsistent Pb/Rsb on identical inputs. The deprecated ``OilPVT.from_harmonize`` default was aligned too. **Behaviour change**: default ``oil_harmonize`` output shifts for callers that relied on the old default (e.g. the ``rsb_frac`` returned when both ``pb`` and ``rsb`` are supplied). Pass ``pbmethod='VELAR'`` explicitly to retain the old behaviour.

- **Nodal HB VLP Rust/Python parity**. Hampson-Brill superficial-gas-velocity and gas-density constants in the Rust accelerator (``35.3741`` molar-volume factor, ``10.73`` vs ``10.732``, ``29.0`` vs ``28.97`` air MW) are aligned with the Python path, making all four VLP methods bit-exact between Rust and Python. The HB-oil/HB-gas parity tests are tightened from ``rtol=1e-3`` to ``1e-4``. Affected HB doc-example BHPs shifted by ~0.05-0.25% to the now-identical Rust/Python value (e.g. ``fbhp`` oil 1771.47 → 1770.58 psi); RST examples updated to match.

- **``simtools.rr_solver`` single-phase guard**. Single-phase feeds (all-liquid or all-vapor) previously fell through to the two-phase Nielsen-Lia solver and returned ``inf``/``nan`` with only a divide-by-zero warning. They are now detected up front and return the trivial split (V=0 or V=1); non-positive K-values raise ``ValueError``.

- **``nodal.fthp`` input validation**. ``fthp`` now validates ``well_type`` / ``vlpmethod`` / ``bhp`` at entry like the other public nodal functions, instead of only erroring transitively through the inner ``fbhp`` solve.

- **``validate_pe_inputs`` array safety**. The ``sg`` and inert-fraction checks are now array-safe (previously ``sg <= 0`` raised an ambiguous-truth-value error on array inputs).

- **Brine undersaturated-compressibility density** now uses the same algebraically reformulated, non-singular Garcia Eq. 18 as the saturated-density step (removes the ``xCO2 → 1`` singularity from the P+1 density used in ``Cf_usat``). Mathematically identical for ``xCO2 < 1``.

- **Removed dead Rust exports** ``oil_bo_mccain_rust`` and ``calc_equilibrium_rust`` (and their now-orphaned internals), which had no Python callers. New Rust-vs-Python parity tests added for ``gas_ponz2p`` and ``simtools.influence_tables`` - the two live accelerated paths that previously had no parity coverage.

- **``oil_matbal`` aquifer water influx**. New optional ``We`` parameter (cumulative aquifer influx in reservoir volume, rb | rm3) brings oil material balance to parity with ``gas_matbal``. The influx is subtracted from underground withdrawal before estimating OOIP (Havlena-Odeh ``F - We = N·[Eo + m·Eg + (1+m)·Efw]``) - feeding both the Python and Rust regression objectives - and a Water Drive Index (``'WDI'``) is added to ``drive_indices`` so the indices sum to 1 under water drive. Fully backward-compatible: with ``We`` omitted, ``WDI`` is all-zero and OOIP/DDI/SDI/CDI are unchanged.

- **Internal: hydrate code split out of ``gas.py``**. ``gas_hydrate``, ``HydrateResult``, the HFT/HFP/Østergaard helpers and their constants now live in ``pyrestoolbox/gas/_hydrate.py`` (≈480 lines carved off the ``gas.py`` monolith). Public API is unchanged - ``gas.gas_hydrate`` and ``gas.HydrateResult`` are re-exported. ``_hydrate`` imports ``gas_water_content`` lazily, so there is no circular import.

- 838 validation tests (up from 812 in 3.4.0).

Changelist in 3.4.0:

- **BREAKING — Brine metric default standardization**: ``CO2_Brine_Mixture`` and ``SoreideWhitson`` constructors now default to ``metric=False`` (oilfield units) to match ``brine_props`` and every other pyrestoolbox API. Previously they defaulted to ``metric=True``. Existing callers that relied on the old default must either pass ``metric=True`` explicitly or switch their input units to psia / degF. RST and notebook examples that relied on the default have been updated to pass ``metric=True`` explicitly.

- **nodal.outflow_curve parameter unification**: ``n_rates`` has been renamed to ``n_points`` to match ``ipr_curve`` / ``operating_point``. ``n_rates`` remains accepted as a deprecated alias (takes precedence when both are passed) — existing callers keep working without changes.

- **Type hints on public APIs**: ``nodal`` (fbhp / fthp / outflow_curve / ipr_curve / operating_point), ``oil`` (oil_co/bt return types, OilPVT constructor + methods, oil_rate_radial/linear), and ``dca`` (arps_*, duong_*, eur, fit_*, ratio_forecast, forecast) public signatures now carry type annotations. IDE autocomplete and static analysis now work out of the box for the primary user surface. ``matbal`` / ``simtools`` / ``plyasunov`` type-hint backfill deferred.

- **Rust-vs-Python parity harness**: 63 new parametrized tests added to ``test_rust_acceleration.py``. Coverage:

  - Nodal gas VLP: 4 methods × 6 (THP, rate) points on a vertical well; WG/BB × 2 points on a deviated well (45-degree lower segment).
  - Nodal oil VLP: 4 methods × 5 (THP, rate) points on a vertical well; WG/BB × 2 points on a deviated well.
  - SoreideWhitson: pure CO2 across 4 (P, T) × 3 salinity points; sour-gas mix (CO2 + H2S + N2) at 2 points; natural-gas-only at 2 points.

  Purpose is prophylactic — catches future one-sided edits between Python and Rust. Tolerance is ``RTOL_MEDIUM = 1e-4``.

- **Rust HB-oil VLP bugfix (correctness)**. Three dimensional bugs in ``src/vlp/segment_oil.rs`` (HB oil only — WG/GRAY/BB oil were coded correctly):

  1. ``min_l`` formula (``rho_g * 62.37 / (rho_g * 62.37 + lsg * 62.37)`` simplified to ``rho_g_lbft3 / (rho_g_lbft3 + lsg_SG)``) mixed lb/cuft with dimensionless SG. Forced ~35% minimum liquid holdup at low rate when Python computed ~2%. Dominant contributor to the divergence.
  2. ``rho_avg`` averaged gas vs stock-tank oil density, not the oil+water mixture density that includes live-oil McCain correction and water cut.
  3. ``ul`` superficial liquid velocity used stock-tank volumetric flow (``ql * 5.615 / 86400 / area``) instead of reservoir-conditions volumetric (``mflow_l / rho_l / area``). Missed the live-oil expansion / water-cut weighting.

  Combined effect was up to ~75% over-prediction of BHP for HB oil wells at low rate / low THP. Python path was always correct; bug was isolated to the Rust accelerator. The ``hb_holdup`` helper signature changed from taking ``lsg`` (implicit ``62.4 * lsg`` internally) to taking ``rho_l`` directly; gas-segment callers updated to pass ``62.4 * lsg_loc`` to preserve prior behaviour. Parity harness drift reduced from 75% to <0.1% on the affected grid.

- **Nodal ``operating_point`` non-monotonic VLP handling**. HB oil VLP curves can be non-monotonic at very low rate (spurious near-shut-in high BHP from holdup correlations). The previous bisection bracketed ``[min_rate, AOF]`` and failed when both endpoints had the same sign despite a crossing existing in between. Replaced with a scan-then-bisect: sample 25 rates, locate all sign changes in the error function, and bisect in the highest-rate bracket — the physical operating point. Sets ``converged=False`` if no sign change is found.

- **HB oil doc-example values updated**. Three expected values changed due to the Rust HB-oil bugfix (Python values unchanged, Rust values now match Python):

  - ``fbhp(thp=200, ..., vlpmethod='HB', well_type='oil', qt_stbpd=2000, ...)``: 2271.72 psi → 1771.47 psi
  - ``fbhp`` with ``oil_pvt=opvt``: 2273.72 psi → 1772.24 psi
  - ``operating_point`` (HB oil, thp=200, reservoir pr=3000): rate 1391.4 → 2019.0 stb/d, bhp 2206.3 → 1778.5 psi

  RST and notebook examples updated to match the corrected values.

- **nodal.fthp**: New reverse-solve function. ``fthp(bhp, completion, ...)`` returns the tubing head pressure that produces the specified BHP under the given VLP correlation. Wraps ``bisect_solve`` over ``fbhp``, accepts the same flow/well parameters as ``fbhp``, and has matching units/metric handling. Users no longer need to hand-roll a bisection for wellhead back-calculation workflows.

- **nodal.fbhp ``return_profile=True``**: New kwarg. When enabled, ``fbhp`` returns a ``NodalResult`` with per-segment-boundary ``md``, ``tvd``, and ``p`` arrays (plus scalar ``bhp``) instead of just the bottom-hole scalar. Exposes the wellbore pressure traverse without forking the internal segment march. Metric-aware outputs.

- **nodal.operating_point extensions**:

  - New ``injection=False`` kwarg, forwarded to internal ``fbhp`` and ``outflow_curve`` calls so injection wells can be solved directly through ``operating_point``. Previously ``injection=True`` was honoured only by direct ``fbhp`` calls.
  - New ``converged`` key in the returned ``NodalResult``. ``True`` when the VLP/IPR bisection succeeded, ``False`` when it fell back to ``rate=0, bhp=pr`` (no intersection). Surfaces a condition that previously looked like a plausible zero-rate answer.

- **nodal.outflow_curve dict key**: Returned ``NodalResult`` now carries both ``'rate'`` and ``'rates'`` entries (same list). ``'rate'`` matches ``ipr_curve`` and ``operating_point``; ``'rates'`` is kept for backward compatibility.

- **SoreideWhitson framework / salinity_method exposed**: New ``framework`` (``'proposed'`` | ``'sw_original'`` | ``'dropin'``) and ``salinity_method`` (``'gamma_phi'`` | ``'embedded'`` | ``'explicit'`` | ``'sechenov'`` | ``'auto'``) constructor kwargs, forwarded to the VLE engine. Defaults unchanged (``'proposed'`` + ``'gamma_phi'``). Combining ``framework='proposed'`` with ``salinity_method='embedded'`` now emits a warning because the engine falls back to ``gamma_phi`` for the ``proposed`` kij set.

- **Brine ``ppm`` / ``wt`` conflict detection**: ``brine_props``, ``CO2_Brine_Mixture``, and ``SoreideWhitson`` now raise ``ValueError("Supply either ... not both.")`` when both aliases are passed. Previously one silently won with no feedback.

- **Oil API cleanup**:

  - Removed private-helper re-exports ``_cofb_mccain``, ``_perrine_co_sat``, ``_resolve_pb_rsb``, ``_build_bot_tables``, ``_format_bot_results`` from ``pyrestoolbox.oil`` public namespace. Still reachable via ``pyrestoolbox.oil._tables`` / ``._density`` / ``._compressibility`` for advanced users; ``simtools.make_bot_og`` updated internally.
  - Dead imports removed across ``_compressibility``, ``_correlations``, ``_density``, ``_harmonize``, ``_tables`` sub-files (11 symbols total).
  - Valko-McCain (2003) coefficient matrices for ``sg_st_gas`` and ``oil_rs_st`` extracted to ``_constants.py`` with paper citations, per CLAUDE.md named-constants rule.

- **Oil documentation sync**: ``docs/oil.rst`` updated — ``oil_deno`` parameter listed as ``sg_o`` (matches code; was incorrectly ``sg_sto``). ``make_bot_og`` result-dict table now documents ``vis_frac`` (was previously omitted despite being returned).

- **DCA internals**: 20 previously inline magic numbers extracted to named constants at the top of ``dca.py`` (Arps b-grid 0.05/0.96/0.01, Duong trap bounds and curve_fit bounds, logistic fit bounds, hyperbolic numerical floor). ``_build_decline_result`` helper deduplicates R-squared/residual/DeclineResult construction across the six ``_fit_*`` helpers (~60 lines removed).

- 726 validation tests (up from 716 in 3.3.0). New coverage for ``fthp`` roundtrip, ``fbhp`` pressure-traverse return, ``operating_point.converged``, ``SoreideWhitson`` framework/salinity_method validation, and the three ``ppm``/``wt`` conflict raises.

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