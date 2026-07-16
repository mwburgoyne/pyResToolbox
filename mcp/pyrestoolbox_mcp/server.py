"""MCP server exposing the pyrestoolbox reservoir engineering library.

Meta-tool design: rather than one MCP tool per library function (125+), a
small set of tools lets an agent discover and call the whole public API:

    list_functions  - module map, or the public functions of one module
    describe        - signature + docstring for module.function
    call            - invoke any public function with keyword arguments
    recommend_methods - structured correlation-method recommendations
    co2_brine_props / sw_brine_props / fit_and_forecast - one-shot wrappers
                      around the stateful classes

The full RST documentation shipped inside the pyrestoolbox wheel is exposed
as MCP resources (docs://gas, docs://oil, ...).
"""

import importlib
import inspect

from mcp.server.fastmcp import FastMCP

from pyrestoolbox_mcp.serialize import to_jsonable

# Modules whose public functions are callable through this server.
ALLOWED_MODULES = {
    "gas": "Gas PVT (Z-factor, viscosity, Bg, Cg, density, pseudopressure), flow rates, pseudoskins, hydrates",
    "oil": "Oil PVT (Pb, Rs, Bo, density, viscosity, compressibility), flow rates, harmonization (scalar inputs only)",
    "brine": "Brine properties: CH4-saturated, CO2-saturated, multicomponent gas-saturated",
    "plyasunov": "Infinite-dilution partial molar volumes of dissolved gases in water (SI units: K, MPa, cm3/mol)",
    "nodal": "VLP/IPR nodal analysis: fbhp, fthp, outflow_curve, ipr_curve, operating_point",
    "dca": "Decline curve analysis: Arps, Duong, modified hyperbolic, two-segment hyperbolic, fitting, EUR",
    "matbal": "Material balance: P/Z gas and Havlena-Odeh oil",
    "simtools": "Simulation helpers: rel-perm tables, aquifer influence, VFP/PVT tables, deck checking",
    "layer": "Lorenz heterogeneity and permeability layering",
    "library": "Component critical property database",
    "recommend": "Correlation method recommendation engine",
    "sensitivity": "Parameter sweeps and tornado analysis",
}

mcp = FastMCP(
    "pyrestoolbox",
    instructions=(
        "Petroleum/reservoir engineering calculations from the pyrestoolbox "
        "library. Start with list_functions() for the module map, "
        "describe('module.function') before calling anything, and "
        "recommend_methods(...) if unsure which correlation method suits a "
        "fluid. All inputs default to oilfield units (psia, deg F, ft, mD, "
        "cP); pass metric=True for Eclipse METRIC units where supported. "
        "Full module documentation is available as docs:// resources."
    ),
)


def _get_module(module: str):
    if module not in ALLOWED_MODULES:
        raise ValueError(f"Unknown module: {module!r}. Available: {sorted(ALLOWED_MODULES)}")
    return importlib.import_module(f"pyrestoolbox.{module}")


def _public_names(mod):
    names = getattr(mod, "__all__", None)
    if names is None:
        names = [n for n in dir(mod) if not n.startswith("_")]
    return [n for n in names if not n.startswith("_")]


def _resolve(name: str):
    """Resolve 'module.function' to (module, attribute)."""
    if "." not in name:
        raise ValueError(
            f"Use 'module.function' form, e.g. 'gas.gas_z'. Got: {name!r}. "
            f"Modules: {sorted(ALLOWED_MODULES)}"
        )
    module, _, attr = name.partition(".")
    mod = _get_module(module)
    if attr.startswith("_") or attr not in _public_names(mod):
        raise ValueError(f"{attr!r} is not a public name in pyrestoolbox.{module}")
    return mod, getattr(mod, attr)


@mcp.tool()
def list_functions(module: str = "") -> dict:
    """List pyrestoolbox modules, or the public functions of one module.

    With no argument, returns {module: description}. With a module name
    (e.g. 'gas'), returns {function: one-line summary} for every public
    function and class in that module.
    """
    if not module:
        return dict(ALLOWED_MODULES)
    mod = _get_module(module)
    out = {}
    for name in _public_names(mod):
        obj = getattr(mod, name, None)
        if callable(obj):
            doc = inspect.getdoc(obj) or ""
            out[name] = doc.split("\n", 1)[0].strip()
    return out


@mcp.tool()
def describe(name: str) -> str:
    """Signature and full docstring for a public function or class.

    name is 'module.function', e.g. 'gas.gas_z' or 'dca.fit_decline'.
    Also see the docs://<module> resources for worked examples.
    """
    _, obj = _resolve(name)
    target = obj.__init__ if inspect.isclass(obj) else obj
    try:
        raw = inspect.signature(target)
        # Strip type annotations - docstrings carry the units/types, and
        # numpy typing unions are unreadable noise in a tool result.
        params = [p.replace(annotation=inspect.Parameter.empty)
                  for p in raw.parameters.values() if p.name != "self"]
        sig = str(raw.replace(parameters=params,
                              return_annotation=inspect.Signature.empty))
    except (TypeError, ValueError):
        sig = "(...)"
    kind = "class" if inspect.isclass(obj) else "function"
    doc = inspect.getdoc(obj) or "(no docstring)"
    return f"{kind} {name}{sig}\n\n{doc}"


@mcp.tool()
def call(name: str, arguments: dict = {}) -> dict:
    """Call any public pyrestoolbox function with keyword arguments.

    name is 'module.function', e.g. 'gas.gas_z'. arguments is a JSON object
    of keyword arguments, e.g. {"p": 2000, "sg": 0.75, "degf": 200}.
    Correlation methods are strings (e.g. "zmethod": "DAK"); an invalid
    method string returns an error listing the valid options. Gas and brine
    functions accept lists for pressure-like inputs; oil functions are
    scalar-only. Classes cannot be called here - use the one-shot wrapper
    tools (co2_brine_props, sw_brine_props, fit_and_forecast) instead.
    """
    _, obj = _resolve(name)
    if inspect.isclass(obj):
        raise ValueError(
            f"{name} is a class. Use describe({name!r}) to inspect it, or a "
            "one-shot wrapper tool (co2_brine_props, sw_brine_props, "
            "fit_and_forecast) if one covers it."
        )
    if not callable(obj):
        raise ValueError(f"{name} is not callable (it is {type(obj).__name__})")
    result = obj(**arguments)
    return {"result": to_jsonable(result)}


@mcp.tool()
def recommend_methods(sg: float = 0.65, co2: float = 0, h2s: float = 0,
                      n2: float = 0, h2: float = 0, api: float | None = None,
                      deviation: float = 0, well_type: str = "gas") -> dict:
    """Recommend correlation methods for a fluid/well, with rationale.

    Composition fractions are molar. api: stock tank oil gravity (include to
    get oil PVT method recommendations). deviation: max wellbore deviation
    from vertical (degrees). well_type: 'gas' or 'oil'. Returns one
    recommendation per category (zmethod, cmethod, vlp_method, ...), each
    with the recommended method, rationale, alternatives and whether the
    choice is mandatory.
    """
    from pyrestoolbox import recommend
    recs = recommend.recommend_methods(sg=sg, co2=co2, h2s=h2s, n2=n2, h2=h2,
                                       api=api, deviation=deviation,
                                       well_type=well_type)
    return to_jsonable(recs)


@mcp.tool()
def co2_brine_props(pres: float, temp: float, ppm: float = 0,
                    cw_sat: bool = False, metric: bool = False) -> dict:
    """CO2-saturated brine properties (Spycher-Pruess mutual solubility).

    pres: pressure (psia, or barsa if metric=True). temp: temperature
    (deg F, or deg C if metric=True). ppm: NaCl-equivalent salinity (weight
    salt per million weight of brine). cw_sat=True also computes the
    saturated pseudo-compressibility. Returns aqueous/vapor compositions,
    gas density (g/cm3), brine density/viscosity/Bw triplets
    [CO2-saturated, pure brine, freshwater], Rs and compressibilities.
    """
    from pyrestoolbox import brine
    mix = brine.CO2_Brine_Mixture(pres=pres, temp=temp, ppm=ppm,
                                  cw_sat=cw_sat, metric=metric)
    keys = ("x", "y", "xSalt", "rhoGas", "bDen", "bVis", "bVisblty", "bw",
            "Rs", "Cf_usat", "Cf_sat")
    return {k: to_jsonable(getattr(mix, k, None)) for k in keys}


@mcp.tool()
def sw_brine_props(pres: float, temp: float, ppm: float = 0,
                   y_CO2: float = 0, y_H2S: float = 0, y_N2: float = 0,
                   y_H2: float = 0, sg: float = 0.65, cw_sat: bool = False,
                   metric: bool = False) -> dict:
    """Multicomponent gas-saturated brine properties (Soreide-Whitson VLE).

    pres: pressure (psia, or barsa if metric=True). temp: temperature
    (deg F, or deg C if metric=True). ppm: NaCl-equivalent salinity.
    y_*: molar fractions of CO2/H2S/N2/H2 in the equilibrium gas; the
    remainder is hydrocarbon gas of specific gravity sg. Returns aqueous
    composition, brine density/viscosity triplets, Rs, water content and
    compressibilities.
    """
    from pyrestoolbox import brine
    sw = brine.SoreideWhitson(pres=pres, temp=temp, ppm=ppm, y_CO2=y_CO2,
                              y_H2S=y_H2S, y_N2=y_N2, y_H2=y_H2, sg=sg,
                              cw_sat=cw_sat, metric=metric)
    keys = ("x", "bDen", "bVis", "bw", "Rs", "Rs_total", "water_content",
            "Cf_usat", "Cf_sat")
    return {k: to_jsonable(getattr(sw, k, None)) for k in keys}


@mcp.tool()
def fit_and_forecast(t: list, q: list, t_end: float, method: str = "best",
                     dt: float = 1.0, q_min: float = 0.0) -> dict:
    """Fit a decline model to rate-time data and forecast to t_end.

    t: time array. q: rate array (same length; any consistent units - the
    dca module is unit-agnostic). method: 'best' (excludes 'mh'/'hyp2'),
    'exponential', 'harmonic', 'hyperbolic', 'duong', 'mh' or 'hyp2'.
    Forecast runs from the fitted model at dt spacing until t_end or until
    rate falls to q_min. Returns fitted parameters (qi, di, b, ...),
    r_squared, EUR and the forecast t/q/cumulative arrays.
    """
    from pyrestoolbox import dca
    fit = dca.fit_decline(t=t, q=q, method=method)
    fc = dca.forecast(fit, t_end=t_end, dt=dt, q_min=q_min)
    return {
        "fit": {k: to_jsonable(getattr(fit, k)) for k in
                ("method", "qi", "di", "b", "a", "m", "dterm", "b2", "telf",
                 "r_squared")},
        "forecast": {"t": to_jsonable(fc.t), "q": to_jsonable(fc.q),
                     "cum": to_jsonable(fc.Qcum), "eur": to_jsonable(fc.eur)},
    }


def _register_doc_resources():
    """Expose every RST file shipped in the pyrestoolbox wheel as docs://<name>."""
    import os
    import pyrestoolbox
    docs = pyrestoolbox.docs_dir()
    for fname in sorted(os.listdir(docs)):
        if not fname.endswith(".rst"):
            continue
        stem = fname[:-4]
        path = os.path.join(docs, fname)

        def _make_reader(p):
            def _read() -> str:
                with open(p, encoding="utf-8") as f:
                    return f.read()
            return _read

        mcp.resource(
            f"docs://{stem}",
            name=f"pyrestoolbox {stem} documentation",
            description=f"RST documentation with worked examples: {stem}",
            mime_type="text/x-rst",
        )(_make_reader(path))


_register_doc_resources()


def build_parser():
    import argparse
    parser = argparse.ArgumentParser(
        prog="pyrestoolbox-mcp",
        description="MCP server exposing the pyrestoolbox library",
    )
    parser.add_argument("--transport", choices=["stdio", "streamable-http"],
                        default="stdio",
                        help="stdio for local hosts (default); streamable-http "
                             "for remote deployment (serves on /mcp)")
    parser.add_argument("--host", default="127.0.0.1",
                        help="bind address for streamable-http (default "
                             "127.0.0.1; use 0.0.0.0 in a container)")
    parser.add_argument("--port", type=int, default=8000,
                        help="port for streamable-http (default 8000)")
    return parser


def main():
    args = build_parser().parse_args()
    if args.transport == "streamable-http":
        mcp.settings.host = args.host
        mcp.settings.port = args.port
        # Independent request handling; no server-side session affinity needed
        # behind a load balancer.
        mcp.settings.stateless_http = True
    mcp.run(transport=args.transport)


if __name__ == "__main__":
    main()
