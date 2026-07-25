# API/UX — oil

Scope: Public-surface ergonomics only. Correctness, refactoring, and duplication are covered by other passes. Focus is on signatures, naming, defaults, docstrings, RST alignment, enum UX, OilPVT API, harmonization clarity, and `__init__.py` re-exports. Review limited to items not already marked DONE in `.claude/docs/review-2026-04-12.md`.

## Critical

- **`oil.rst` `oil_deno` parameter name wrong.** `pyrestoolbox/docs/oil.rst:780`, `785`, `1201` document the parameter as `sg_sto` but the code parameter is `sg_o` (see `oil/_density.py:48`, `oil/_correlations.py:373`). Users following the docs will hit `TypeError: unexpected keyword argument 'sg_sto'`. Also the doc text at line 782 ("Specific gravity of stock tank liquid (rel water). Will calculate from api if not specified") should be reworded – spec says "If both defined api value will prevail" which means sg_sto/sg_o is **overwritten** by api, not "used if api missing". Fix: rename `sg_sto` → `sg_o` in RST, and clarify that api takes precedence when both are supplied.

- **`make_bot_og` result dict: docs and code use different keys.** `oil.rst:1157` documents `'rsb_scale'` and code at `oil/_tables.py:259` returns `"rsb_scale"` — consistent. BUT the harmonize return tuple key is named `rsb_frac` (`_harmonize.py:57`), `OilPVT` attribute is `rsb_frac`, and the `vis_frac` attribute/tuple field is also `..._frac`. So the BOT dict key `rsb_scale` is the single outlier. Either rename to `rsb_frac` (breaking) or add an alias. At minimum, call this out in the RST description so users who chain `oil_harmonize` → `make_bot_og` aren't confused. `vis_frac` **is** returned in the BOT dict (`_tables.py:260`) but NOT listed in the RST result-dict table (`oil.rst:1126-1162`) — add it.

## High

- **`pbmethod` default varies across functions for no obvious reason.** `oil_pbub` default `VALMC`; `oil_rs` default `VALMC`; `oil_co` default `VALMC`; `oil_bt` default `VALMC`; but `oil_harmonize` default `VELAR` (`_harmonize.py:27`); `simtools.make_bot_og` default `VELAR` (`simtools.py:1758`); `OilPVT` default `VALMC` (`_pvt_class.py:39`); `OilPVT.from_harmonize` default `VELAR` (`_pvt_class.py:75`). Pick one. A user who harmonizes with VELAR then computes Bo with VALMC gets silently inconsistent results across a workflow. Recommend: standardize on `VALMC` (already the module-wide default for oil_pbub/oil_rs/oil_co/oil_bt) and document why the few call-sites differ — or make them match.

- **`OilPVT.density()` and `OilPVT.bo()` ignore the class's `denomethod`.** The class constructor accepts `rsmethod`, `pbmethod`, `bomethod` but NOT `denomethod`, yet `OilPVT.density()` calls `oil_deno(...)` without passing any denomethod — which is fine today because SWMH is the only option, but if a second density method is ever added the class silently won't track it. Either add `denomethod` to the constructor now, or document clearly that SWMH is fixed. Same consideration for `zmethod`/`cmethod` (currently unused since OilPVT doesn't expose Bt/Co, so perhaps not yet relevant).

- **`oil_harmonize` requires `degf > 0` but the parameter defaults to 0.** `_harmonize.py:69` calls `validate_pe_inputs(degf=degf)` which will raise if degf ≤ 0, yet `degf=0` is the default (`_harmonize.py:20`). Any call that omits degf errors out. Either make degf a required positional argument, or document explicitly "degf: REQUIRED (no default)". Same pattern in the `OilPVT` constructor which auto-harmonizes "when degf > 0" (line 45) — rename the docstring contract to be explicit about mandatory-vs-optional.

- **`OilPVT` constructor: `pb` is required positional but `rsb=0` is optional with no indication either is needed.** If `degf=0` (the default) and `rsb<=0`, the constructor raises `"Either rsb or degf must be specified"` (`_pvt_class.py:53`). But if pb is a bad value (e.g., 0 accidentally) no error fires until a downstream call. The error message should also be rephrased — "Either rsb or degf must be specified" is misleading; the actual rule is "rsb required unless degf>0 triggers auto-harmonization". Suggest: `"rsb must be > 0, OR degf > 0 to auto-harmonize rsb from pb"`.

- **`_utils.get_real_part` has no leading underscore — looks public.** It is an internal helper (strips imaginary part of complex numbers from Velarde correlations). It's not in `__all__` but it IS imported across sub-files, so renaming to `_get_real_part` would be a surgical change. Current name suggests it's a public utility. Not re-exported from `__init__.py`, so low-impact, but inconsistent with codebase convention.

- **`__init__.py` re-exports three private BOT helpers (`_resolve_pb_rsb`, `_build_bot_tables`, `_format_bot_results`) and two private numeric helpers (`_cofb_mccain`, `_perrine_co_sat`).** Lines 73, 75-76. These are imported because `simtools.make_bot_og` uses them, but re-exporting private names from a package `__init__.py` pollutes the public namespace (`oil._cofb_mccain` is now reachable). Consider keeping them reachable via `oil._tables._build_bot_tables` (`simtools` already accesses them through `_oil_impl`, so no change needed there) and removing from `oil/__init__.py`.

## Medium

- **`oil_rs` signature: `sg_g=0` is absent but other oil PVT functions include it.** Compare `oil_rs(api, degf, sg_sp, p, pb=0, rsb=0, ...)` (`_correlations.py:246`) vs `oil_pbub(api, degf, rsb, sg_g=0, sg_sp=0, ...)` (`_correlations.py:30`) and `oil_bo(p, pb, degf, rs, rsb, sg_o, sg_g=0, sg_sp=0, ...)`. `oil_rs` hardcodes `check_sgs(sg_g=0, sg_sp=sg_sp)` (line 292), so users who have `sg_g` but not `sg_sp` cannot call `oil_rs` without first converting. Add `sg_g: float = 0` for consistency.

- **Parameter-order inconsistency: `oil_bo` puts `sg_o` before `sg_g/sg_sp` as a required positional.** `oil_bo(p, pb, degf, rs, rsb, sg_o, sg_g=0, sg_sp=0, ...)`. Every other function takes sg_g/sg_sp as kwargs with defaults. Users calling `oil_bo(p, pb, degf, rs, rsb, sg_g=0.75, sg_o=0.82)` works via kwargs but unfamiliar users hit `TypeError: missing argument 'sg_o'`. Consider making `sg_o=0` optional (consistent with `oil_deno` which accepts either sg_o or api).

- **`oil_deno` default `pb=1e6` is a magic number.** `_density.py:47, 56`. The intent ("not used for densities below Pb") is documented but `1e6` psia is an arbitrary large value. Use `pb=0` and branch on `pb <= 0`, matching the pattern used in every other oil function (oil_co, oil_bt, make_bot_og).

- **Docstring RST signature strings are outdated.**
  - `oil.rst:551` shows `oil_co(..., undersaturated_only` absent — the actual oil_co signature includes `undersaturated_only: bool = False` (`_compressibility.py:61`), documented later in the prose but missing from the pseudo-signature line. Similarly `pi` kwarg exists (`_compressibility.py:59`) but isn't shown.
  - `oil.rst:354` — `oil_pbub` signature shown uses `pbmethod ='VALMC'` (correct) but missing a visual `sg_o` context.
  - `oil.rst:1041` — `make_bot_og` signature is missing the `pvto`, `vis_frac`, and `metric` kwargs (all present in code at `simtools.py:1760-1762`).

- **Error messages could be more actionable.** Examples:
  - `"Need valid values for rs, api, sg_g and degf for Standing Pb calculation"` (`_correlations.py:64`) — should tell the user which one was zero.
  - `"Either oil_pvt or both uo and bo must be specified"` (`_rate.py:84`) — fine, but no info on which (uo vs bo) was missing.
  - `"Could not solve Pb & Rsb for these combination of inputs"` (`_harmonize.py:114`) — print the inputs in the message so users can reproduce/debug.

- **`oil_rs_bub` hidden fallback to Velarde for extreme inputs isn't documented in the RST.** The RST says "will revert to the Standing method" (`oil.rst:421`), but the code uses Velarde as the initial guess then Newton iteration (`_correlations.py:188`). Fix: align docs with code.

- **`oil_rate_radial`/`oil_rate_linear` implicit behavior when `oil_pvt` is provided.** Passing `oil_pvt` silently switches `vogel=True` (`_rate.py:82, 170`) and overrides user-specified `pb`. This is noted in docs but is a surprising side-effect; consider raising a warning when the user explicitly passed `vogel=False` and `oil_pvt` is not None, so they know Vogel was auto-enabled.

- **`_perrine_co_sat` uses positional args while every other internal helper uses kwargs.** `_compressibility.py:22-23` has 13 positional parameters. If someone adds a parameter later the call sites at lines 166 and 258 both break silently. Use kwargs.

## Low

- **`oil_harmonize_pb_rsb` deprecation.** The wrapper is marked deprecated in the docstring but does not emit `DeprecationWarning`. Add a `warnings.warn(..., DeprecationWarning, stacklevel=2)`.

- **`OilPVT.from_harmonize` deprecation.** Same — marked deprecated in docstring (`_pvt_class.py:77`) but no `DeprecationWarning`.

- **`oil.make_bot_og` is deprecated (per `_tables.py:270`) but only prose mentions it — no `DeprecationWarning`.** The RST notes "v3.0: see also simtools.make_bot_og" (`oil.rst:167, 1034`) but running `oil.make_bot_og()` doesn't warn.

- **`oil_rs_st` parameter `degf_sp` is lowercase-f but accepts only FIELD temperature when `metric=False`** — consistent with other functions, but the name `degf_sp` visually reads as "always degF". This is a codebase-wide convention (not unique to oil) so flagging only because `oil.rst:323` labels it "Separator temperature (deg F, or deg C if metric=True)". OK as-is but worth documenting the naming rule once somewhere.

- **Twu `oil_twu_props` docstring double-documents return type.** `_utils.py:72-74` — the "Returns Tuple of sg, tb ..." sentence is repeated twice (lines 72-74). Trim the duplication.

- **`sg_evolved_gas`, `sg_st_gas`, `sgg_wt_avg` have no `metric` param** (`oil.rst:130` notes this explicitly). Fine design choice, but inconsistent API footprint — `oil_rs_st` (a sibling separator function) DOES accept metric. Pick one rule: either all separator utilities accept metric, or none do.

- **`oil.rst` has a heading-underline length inconsistency** (`oil_ja_sg`, `oil_rs_st`, `oil_twu_props` sections). Underlines shorter than heading text — Sphinx may emit warnings. E.g. `oil.rst:185-186`: heading is `pyrestoolbox.oil.oil_ja_sg` (27 chars) but the `=====` underline below is only 22 chars. Same at lines 226-227, 295-296. Cosmetic but will raise Sphinx build warnings.

- **`oil.rst` "Oil Bubble Point Pressure" section lists `pbmethod` default as `'VALMC'` in the code-block at line 354 but describes defaults as `'VELAR'` for `pbmethod` in the class-options table at line 35.** Two conflicting defaults listed in one doc page. The actual code default in `oil_pbub` is `VALMC`. Fix the table at line 35.

- **`co_method` enum has only one value (`EXPLT`).** Noted in the RST (line 49-51). Having an enum for a single-option method is API overhead that forces users to pass a method string. Consider keeping `comethod` kwarg but omitting the enum until a 2nd method is added; or leave as-is for forward compat (current choice — fine, just verbose).

- **`oil_deno` accepts both `sg_o` and `api` with "both defined, api prevails" semantics.** `_density.py:184`: `else: sg_o = oil_sg(api)` — silently overwrites user's sg_o. Either warn when both are provided and differ, or make it an XOR error. Currently a user who provides inconsistent sg_o and api gets no indication their sg_o was ignored.

- **`OilPVT.rs()`, `.bo()`, `.density()`, `.viscosity()` accept array `p` but require scalar `degf`.** Documented in each docstring but surprising — most array-capable APIs in pyrestoolbox allow either. Not a bug; worth a one-line callout in the RST class overview.
