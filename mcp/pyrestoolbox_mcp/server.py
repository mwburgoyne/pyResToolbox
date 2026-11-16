"""MCP server exposing the pyrestoolbox reservoir engineering library.

Meta-tool design: rather than one MCP tool per library function (125+), a
small set of tools lets an agent discover and call the whole public API:

    search          - engineering task in plain words -> matching functions
    read_docs       - a module's full RST documentation with worked examples
    list_functions  - module map, or the public functions of one module
    describe        - signature + docstring for module.function
    call            - invoke any public function with keyword arguments
    recommend_methods - structured correlation-method recommendations
    co2_brine_props / sw_brine_props / fit_and_forecast - one-shot wrappers
                      around the stateful classes

The full RST documentation shipped inside the pyrestoolbox wheel is exposed
as MCP resources (docs://gas, docs://oil, ...).
"""

import functools
import importlib
import inspect
import json
import os
import re
import threading
import warnings

import anyio
from mcp.types import ToolAnnotations

try:  # mcp >= 2 renamed FastMCP to MCPServer
    from mcp.server.mcpserver import MCPServer as _MCPServer
    from mcp.server.mcpserver.exceptions import ToolError
    _MCP_V2 = True
except ImportError:  # mcp 1.x
    from mcp.server.fastmcp import FastMCP as _MCPServer
    from mcp.server.fastmcp.exceptions import ToolError
    _MCP_V2 = False

from pyrestoolbox_mcp import __version__
from pyrestoolbox_mcp.serialize import to_jsonable

# Server policy, set from the command line in main().
# File I/O: a few simtools functions read or write files on the host running
# the server. Refused unless the operator opts in, because over HTTP the caller
# is not the host's owner.
ALLOW_FILE_IO = False
# Largest tool result returned, in characters of JSON. A 20,000-row table is
# about 3 MB, which no agent context can use; the error tells it to ask for less.
MAX_RESULT_CHARS = 500_000

_FILE_IO_FUNCTIONS = {"simtools.ix_extract_problem_cells", "simtools.zip_check_sim_deck"}
_FILE_IO_ARGS = {"export", "filename", "files2scrape", "tozip"}

# Modules whose public functions are callable through this server.
ALLOWED_MODULES = {
    "gas": "Gas PVT (Z-factor, viscosity, Bg, Cg, density, pseudopressure), flow rates, pseudoskins, hydrates",
    "oil": "Oil PVT (Pb, Rs, Bo, density, viscosity, compressibility), flow rates, harmonization (scalar inputs only)",
    "brine": "Brine properties: CH4-saturated, CO2-saturated, multicomponent gas-saturated",
    "plyasunov": "Infinite-dilution partial molar volumes of dissolved gases in water (SI units: K, MPa, cm3/mol)",
    "nodal": "VLP/IPR nodal analysis: fbhp, fthp, outflow_curve, ipr_curve, operating_point; liquid-loading critical rate, API RP 14E erosional velocity, wellhead chokes (Gilbert family, Sachdeva)",
    "pipeline": "Gas pipeline flow (one segment): General Flow (Colebrook, AGA), Weymouth, Panhandle A/B, IGT; solves p1, p2 or rate",
    "crm": "Waterflood capacitance-resistance models (CRMP/CRMT): connectivity and time-constant fit from injection/production history, prediction, oil cut",
    "dca": "Decline curve analysis: Arps, Duong, modified hyperbolic, two-segment hyperbolic, fitting, EUR",
    "matbal": "Material balance: P/Z gas and Havlena-Odeh oil, Carter-Tracy aquifer influx",
    "simtools": "Simulation helpers: rel-perm tables, aquifer influence, VFP/PVT tables, deck checking",
    "layer": "Lorenz heterogeneity and permeability layering",
    "library": "Component critical property database",
    "recommend": "Correlation method recommendation engine",
    "sensitivity": "Parameter sweeps and tornado analysis",
}

mcp = _MCPServer(
    "pyrestoolbox",
    instructions=(
        "Petroleum/reservoir engineering calculations from the pyrestoolbox "
        "library. Start with list_functions() for the module map, "
        "describe('module.function') before calling anything, and "
        "recommend_methods(...) if unsure which correlation method suits a "
        "fluid. All inputs default to oilfield units (psia, deg F, ft, mD, "
        "cP); pass metric=True for Eclipse METRIC units where supported. "
        "IMPORTANT: parameter NAMES keep their oilfield spelling under "
        "metric=True. `degf` then takes deg C, `p`/`pr`/`pwf` take barsa, and "
        "lengths take m; the docstring unit notation `(deg F | deg C)` gives "
        "the field value first and the metric value second. Passing a "
        "Fahrenheit value to `degf` with metric=True converts it twice and "
        "returns a wrong answer silently. Standard volumes are always at "
        "60 deg F and 14.696 psia, so metric sm3 is on that basis. "
        "Gas rates are Mscf/d everywhere (sm3/d with metric=True). "
        "Functions that take objects (Completion, Reservoir, GasPVT, OilPVT, "
        "decline results) accept them through call() as JSON objects of their "
        "constructor keywords; see the call tool description. "
        "search('liquid loading') finds functions by engineering task; "
        "read_docs('gas') returns a module's full documentation with worked "
        "examples (also served as docs:// resources). Errors carry the "
        "library's message (valid method names, missing arguments, range "
        "limits): read it and retry. Results may carry `notes` (library "
        "validity warnings): pass them on with the number."
    ),
    **({"version": __version__} if _MCP_V2 else {}),
)
if not _MCP_V2:  # 1.x takes no version argument and would report the SDK's own
    mcp._mcp_server.version = __version__

# Every tool is a pure calculation: no side effects, same answer for the same
# input, no network. Hosts can show or auto-approve on these hints.
_READ_ONLY = ToolAnnotations(readOnlyHint=True, destructiveHint=False,
                             idempotentHint=True, openWorldHint=False)

# One library call at a time per process. warnings.catch_warnings swaps
# process-global state, so two concurrent calls would trade or lose each
# other's notes; the library's pure-Python paths hold the GIL anyway. Scale
# out with processes or replicas, not threads.
_LIB_LOCK = threading.Lock()


def _check_size(name, out):
    n = len(out) if isinstance(out, str) else len(json.dumps(out, default=str))
    if n > MAX_RESULT_CHARS:
        raise ValueError(
            f"{name} result is {n:,} characters of JSON, over this server's "
            f"{MAX_RESULT_CHARS:,} limit. Ask for fewer rows or points "
            "(e.g. rows, nrows, n_points), or a narrower range."
        )


def _tool(fn):
    """Register fn as an MCP tool and return fn itself (tests call it directly).

    The registered adapter runs fn in a worker thread under the library lock,
    so a slow call never stalls the event loop (pings, other clients), and
    re-raises any failure as ToolError. The SDK passes a ToolError's message to
    the model; mcp 2.x replaces every other exception with a generic "Error
    executing tool", which would hide the valid-options list an agent needs to
    correct itself.
    """
    @functools.wraps(fn)
    async def adapter(*args, **kwargs):
        def work():
            with _LIB_LOCK:
                out = fn(*args, **kwargs)
            _check_size(fn.__name__, out)
            return out
        try:
            return await anyio.to_thread.run_sync(work)
        except ToolError:
            raise
        except Exception as exc:
            raise ToolError(f"{type(exc).__name__}: {exc}") from exc

    mcp.tool(annotations=_READ_ONLY)(adapter)
    return fn


def _engine() -> dict:
    """Which library produced a number: package version and computation path."""
    import pyrestoolbox
    from pyrestoolbox._accelerator import get_status
    st = get_status()
    return {"pyrestoolbox": pyrestoolbox.__version__,
            "rust": bool(st["rust_available"] and not st["forced_python"])}


def _run(fn, *args, **kwargs):
    """Call fn, returning (result, notes) with every library warning captured.

    Validity and range warnings would otherwise go to the server's stderr and
    never reach the agent, which would quote the number without its caveat.
    """
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = fn(*args, **kwargs)
    return result, list(dict.fromkeys(str(w.message) for w in caught))


def _envelope(payload: dict, notes: list) -> dict:
    """Attach captured warnings (when any) and the engine stamp to a tool result."""
    if notes:
        payload["notes"] = notes
    payload["engine"] = _engine()
    return payload


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


@_tool
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


@_tool
def describe(name: str) -> str:
    """Signature and full docstring for a public function or class.

    name is 'module.function', e.g. 'gas.gas_z' or 'dca.fit_decline'.
    Also see the docs://<module> resources for worked examples.

    Units live in the docstring, not the signature, which is why annotations
    are stripped below. A parameter documented `(deg F | deg C)` takes the
    first under the default and the second when metric=True; the parameter
    name itself never changes.
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


@_tool
def call(name: str, arguments: dict | None = None) -> dict:
    """Call any public pyrestoolbox function with keyword arguments.

    name is 'module.function', e.g. 'gas.gas_z'. arguments is a JSON object
    of keyword arguments, e.g. {"p": 2000, "sg": 0.75, "degf": 200}.
    That example is oilfield units. With {"metric": true} the same keys take
    metric values, so `degf` would be deg C and `p` barsa; the names do not
    change. Call describe(name) first and read its unit notation: `(deg F |
    deg C)` means field first, metric second.
    Correlation methods are strings (e.g. "zmethod": "DAK"); an invalid
    method string returns an error listing the valid options. Gas and brine
    functions accept lists for pressure-like inputs; oil functions are
    scalar-only. Gas rates are Mscf/d (sm3/d with metric=True).

    Object arguments are passed as JSON objects of their constructor keywords
    and built here (describe the class for its keywords):
      completion -> nodal.Completion   (a "segments" list of objects -> WellSegment)
      reservoir  -> nodal.Reservoir
      gas_pvt    -> gas.GasPVT
      oil_pvt    -> oil.OilPVT
      result     -> dca.DeclineResult (forecast) or dca.RatioResult (ratio_forecast);
                    pass the object a fit_decline / fit_ratio call returned
      ratios     -> {name: RatioResult object}
      func       -> a 'module.function' name (sensitivity.sweep / tornado); objects
                    inside base_kwargs are built by the same rules
    e.g. call('nodal.fbhp', {"thp": 500, "well_type": "gas", "qg_mscfd": 5000,
    "completion": {"tid": 2.441, "length": 10000, "tht": 100, "bht": 200}}).
    Classes themselves cannot be called here; one-shot wrappers
    (co2_brine_props, sw_brine_props, fit_and_forecast) cover the common cases.
    Tables come back inline; `export` and other file arguments are refused
    unless the server was started with --allow-file-io.
    """
    arguments = arguments or {}
    _, obj = _resolve(name)
    if inspect.isclass(obj):
        raise ValueError(
            f"{name} is a class. Use describe({name!r}) to inspect it, or a "
            "one-shot wrapper tool (co2_brine_props, sw_brine_props, "
            "fit_and_forecast) if one covers it."
        )
    if not callable(obj):
        raise ValueError(f"{name} is not callable (it is {type(obj).__name__})")
    _check_file_io(name, arguments)
    result, notes = _run(obj, **_build_objects(obj, arguments))
    return _envelope({"result": to_jsonable(result)}, notes)


def _check_file_io(name, arguments):
    if ALLOW_FILE_IO:
        return
    hint = ("This server does not touch files on its host (start it with "
            "--allow-file-io for local use).")
    if name in _FILE_IO_FUNCTIONS:
        raise ValueError(f"{name} reads files on the server host. {hint}")
    used = sorted(k for k in _FILE_IO_ARGS & set(arguments) if arguments[k])
    if used:
        raise ValueError(f"{used} would read or write files on the server host. "
                         f"Omit them: the table is returned in the result. {hint}")


def _library_class(dotted: str):
    module, cls = dotted.split(".")
    return getattr(_get_module(module), cls)


# Keyword -> library class built from a JSON object. The library names these
# parameters consistently, so one table covers every function that takes them.
_OBJECT_PARAMS = {
    "completion": "nodal.Completion",
    "reservoir": "nodal.Reservoir",
    "gas_pvt": "gas.GasPVT",
    "oil_pvt": "oil.OilPVT",
}


def _build_objects(fn, arguments: dict) -> dict:
    """Turn JSON objects in `arguments` into the library objects fn expects.

    Raises ValueError naming the parameter and class when construction fails,
    so an agent learns which keywords to fix instead of an AttributeError from
    deep inside the library.
    """
    try:
        annotations = {k: str(v.annotation) for k, v in inspect.signature(fn).parameters.items()}
    except (TypeError, ValueError):
        annotations = {}
    out = {}
    for key, value in arguments.items():
        out[key] = _build_one(key, value, annotations.get(key, ""))
    return out


def _build_one(key, value, annotation):
    if key == "func" and isinstance(value, str):
        _, fn = _resolve(value)
        return fn
    if key == "base_kwargs" and isinstance(value, dict):
        return {k: _build_one(k, v, "") for k, v in value.items()}
    if not isinstance(value, dict):
        return value
    if key == "ratios":
        return {name: _construct("dca.RatioResult", key, spec) for name, spec in value.items()}
    if key == "result":
        dotted = "dca.RatioResult" if "RatioResult" in annotation else "dca.DeclineResult"
        return _construct(dotted, key, value)
    if key in _OBJECT_PARAMS:
        if key == "completion" and isinstance(value.get("segments"), list):
            value = dict(value, segments=[
                _construct("nodal.WellSegment", "completion.segments", seg)
                if isinstance(seg, dict) else seg for seg in value["segments"]])
        return _construct(_OBJECT_PARAMS[key], key, value)
    return value


def _construct(dotted, key, spec):
    cls = _library_class(dotted)
    try:
        return cls(**spec)
    except TypeError as exc:
        raise ValueError(
            f"Could not build {dotted} for '{key}' from {sorted(spec)}: {exc}. "
            f"Use describe({dotted!r}) for its keywords."
        ) from exc


@_tool
def recommend_methods(sg: float = 0.65, co2: float = 0, h2s: float = 0,
                      n2: float = 0, h2: float = 0, api: float | None = None,
                      deviation: float = 0, well_type: str = "gas", he: float = 0) -> dict:
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
                                       well_type=well_type, he=he)
    return to_jsonable(recs)


@_tool
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
    mix, notes = _run(brine.CO2_Brine_Mixture, pres=pres, temp=temp, ppm=ppm,
                      cw_sat=cw_sat, metric=metric)
    keys = ("x", "y", "xSalt", "rhoGas", "bDen", "bVis", "bVisblty", "bw",
            "Rs", "Cf_usat", "Cf_sat")
    return _envelope({k: to_jsonable(getattr(mix, k, None)) for k in keys}, notes)


@_tool
def sw_brine_props(pres: float, temp: float, ppm: float = 0,
                   y_CO2: float = 0, y_H2S: float = 0, y_N2: float = 0,
                   y_H2: float = 0, sg: float = 0.65, cw_sat: bool = False,
                   metric: bool = False, y_He: float = 0) -> dict:
    """Multicomponent gas-saturated brine properties (Soreide-Whitson VLE).

    pres: pressure (psia, or barsa if metric=True). temp: temperature
    (deg F, or deg C if metric=True). ppm: NaCl-equivalent salinity.
    y_*: molar fractions of CO2/H2S/N2/H2/He in the equilibrium gas; the
    remainder is hydrocarbon gas of specific gravity sg. Returns aqueous
    composition, brine density/viscosity triplets, Rs, water content and
    compressibilities.
    """
    from pyrestoolbox import brine
    sw, notes = _run(brine.SoreideWhitson, pres=pres, temp=temp, ppm=ppm, y_CO2=y_CO2,
                     y_H2S=y_H2S, y_N2=y_N2, y_H2=y_H2, sg=sg,
                     cw_sat=cw_sat, metric=metric, y_He=y_He)
    keys = ("x", "bDen", "bVis", "bw", "Rs", "Rs_total", "water_content",
            "Cf_usat", "Cf_sat")
    return _envelope({k: to_jsonable(getattr(sw, k, None)) for k in keys}, notes)


@_tool
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
    fit, notes = _run(dca.fit_decline, t=t, q=q, method=method)
    fc, notes2 = _run(dca.forecast, fit, t_end=t_end, dt=dt, q_min=q_min)
    return _envelope({
        "fit": {k: to_jsonable(getattr(fit, k)) for k in
                ("method", "qi", "di", "b", "a", "m", "dterm", "b2", "telf",
                 "r_squared")},
        "forecast": {"t": to_jsonable(fc.t), "q": to_jsonable(fc.q),
                     "cum": to_jsonable(fc.Qcum), "eur": to_jsonable(fc.eur)},
    }, list(dict.fromkeys(notes + notes2)))


def _doc_paths() -> dict:
    """{stem: path} for every RST file shipped in the pyrestoolbox wheel."""
    import pyrestoolbox
    docs = pyrestoolbox.docs_dir()
    return {f[:-4]: os.path.join(docs, f) for f in sorted(os.listdir(docs)) if f.endswith(".rst")}


@_tool
def read_docs(topic: str = "index") -> str:
    """Full documentation for one module: function list, signatures, input
    and return tables, and worked examples whose printed outputs the test
    suite verifies.

    topic is a module name ('gas', 'oil', 'nodal', ...), 'index' for the map
    of engineering tasks to functions, or 'changelist' for release history.
    Same content as the docs:// resources, for hosts that do not surface
    resources. Prefer describe() for one function; read_docs for the
    background, method choices and examples.
    """
    paths = _doc_paths()
    if topic not in paths:
        raise ValueError(f"Unknown topic {topic!r}. Topics: {sorted(paths)}")
    with open(paths[topic], encoding="utf-8") as f:
        return f.read()


_STOPWORDS = {"a", "an", "and", "the", "of", "for", "to", "in", "on", "from",
              "with", "by", "or", "at", "is", "how", "what", "do", "i", "my"}


# Engineering words that the docstrings do not use, mapped to ones they do.
_SYNONYMS = {
    "place": "ogip ooip material balance", "ogip": "material balance gas",
    "ooip": "material balance oil", "gip": "material balance gas",
    "reserves": "eur material balance", "bhp": "bottomhole flowing",
    "pwf": "bottomhole flowing ipr", "vlp": "bottomhole flowing tubing",
    "tubing": "vlp fbhp", "lift": "vlp fbhp", "solubility": "saturated rs",
    "dissolved": "saturated rs", "pvt": "properties table",
    "loading": "critical rate", "deliverability": "aof backpressure",
    "kr": "relative permeability", "relperm": "relative permeability",
    "forecast": "decline forecast", "z": "z factor deviation",
}

# The one-shot tools, so a search can land on them as well as on call().
_WRAPPER_ENTRIES = [
    ("tool:co2_brine_props", "CO2-saturated brine: CO2 solubility, density, viscosity, Bw, Rs",
     "CO2 solubility in brine, CO2-water mutual solubility"),
    ("tool:sw_brine_props", "Gas-saturated brine for CH4/CO2/H2S/N2/H2/He mixtures: solubility, density, Rs, water content",
     "Gas solubility in brine, sour gas, hydrogen storage"),
    ("tool:fit_and_forecast", "Fit a decline curve to rate-time data and forecast with EUR",
     "Decline curve fit and production forecast"),
]


def _words(text):
    return set(re.findall(r"[a-z0-9]+", text.lower()))


def _search_entry(name, summary, task):
    strong = _words(re.sub(r"[._:]", " ", name) + " " + summary)
    return (name, summary, task, strong, _words(task) - strong)


@functools.lru_cache(maxsize=1)
def _search_index():
    """[(name, summary, task, strong_words, task_words)] over every public function.

    The word set joins the function name, its docstring summary and every
    task row of docs/index.rst that names it, so a query in engineering words
    ("liquid loading") finds a function whose name does not contain them
    (nodal.critical_gas_rate). Classes are left out: call() cannot run them.
    """
    tasks = {}
    with open(_doc_paths()["index"], encoding="utf-8") as f:
        for line in f:
            cols = re.split(r"\s{2,}", line.strip())
            if len(cols) == 2 and "." in cols[1] and not cols[0].startswith("="):
                for fn_name in re.split(r",\s*", cols[1]):
                    tasks.setdefault(fn_name.strip(), []).append(cols[0])
    entries = []
    for module in ALLOWED_MODULES:
        mod = _get_module(module)
        for n in _public_names(mod):
            obj = getattr(mod, n, None)
            if not callable(obj) or inspect.isclass(obj):
                continue
            name = f"{module}.{n}"
            summary = (inspect.getdoc(obj) or "").split("\n", 1)[0].strip()
            task = "; ".join(tasks.get(name, []))
            entries.append(_search_entry(name, summary, task))
    for name, summary, task in _WRAPPER_ENTRIES:
        entries.append(_search_entry(name, summary, task))
    return entries


@_tool
def search(query: str, limit: int = 12) -> dict:
    """Find functions for an engineering task, best match first.

    query is plain words, e.g. 'liquid loading', 'hydrate inhibitor dose',
    'z factor sour gas', 'black oil table', 'erosional velocity'. Matches the
    documentation's task index plus every public function's name and
    summary. Returns {matches: [{name, summary, task}]}; then describe(name)
    and call it.
    Names starting 'tool:' are the one-shot tools of this server.
    """
    import math
    words = [w for w in re.findall(r"[a-z0-9]+", query.lower()) if w not in _STOPWORDS]
    if not words:
        raise ValueError("Give a few words describing the task, e.g. 'gas viscosity'.")
    words += [x for w in words for x in _SYNONYMS.get(w, "").split()]
    index = _search_index()

    def has(wordset, w):
        # Prefix match ('inhibitor' finds 'inhibitors'), but only for words of
        # three or more letters, so 'z' matches the token 'z' and nothing else.
        return w in wordset or (len(w) > 2 and any(x.startswith(w) for x in wordset))

    df = {w: sum(has(e[3] | e[4], w) for e in index) for w in set(words)}
    scored = []
    for name, summary, task, strong, task_words in index:
        # Rare words carry the meaning ('loading' outranks 'gas'). A match in
        # the name or summary counts double one in a shared task-index row.
        score = 0.0
        for w in words:
            idf = math.log((1 + len(index)) / (1 + df[w]))
            score += idf if has(strong, w) else 0.5 * idf if has(task_words, w) else 0
        if score > 0:
            scored.append((score, name, summary, task))
    scored.sort(key=lambda s: (-s[0], s[1]))
    return {"matches": [{"name": n, "summary": s, "task": t}
                        for _, n, s, t in scored[:max(1, limit)]]}


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
    parser.add_argument("--allow-file-io", action="store_true",
                        help="let simtools functions read or write files on this "
                             "host (export=True, deck scrapers); local use only")
    parser.add_argument("--max-result-chars", type=int, default=MAX_RESULT_CHARS,
                        help="largest tool result returned, in characters of JSON "
                             f"(default {MAX_RESULT_CHARS:,})")
    parser.add_argument("--version", action="version",
                        version=f"pyrestoolbox-mcp {__version__}")
    return parser


def main():
    global ALLOW_FILE_IO, MAX_RESULT_CHARS
    args = build_parser().parse_args()
    ALLOW_FILE_IO = args.allow_file_io
    MAX_RESULT_CHARS = args.max_result_chars
    if args.transport != "streamable-http":
        mcp.run(transport=args.transport)
        return
    # stateless_http: independent request handling; no server-side session
    # affinity needed behind a load balancer.
    http = {"host": args.host, "port": args.port, "stateless_http": True}
    if _MCP_V2:
        mcp.run(transport=args.transport, **http)
    else:  # 1.x reads these from settings, not run()
        for key, value in http.items():
            setattr(mcp.settings, key, value)
        mcp.run(transport=args.transport)


if __name__ == "__main__":
    main()
