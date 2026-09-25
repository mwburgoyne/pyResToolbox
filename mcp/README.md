# pyrestoolbox-mcp

MCP (Model Context Protocol) server exposing the
[pyrestoolbox](https://github.com/mwburgoyne/pyResToolbox) reservoir
engineering library to agents without a Python execution environment
(Claude Desktop, claude.ai connectors, other MCP hosts).

If your agent CAN run Python, prefer `pip install pyrestoolbox` directly -
the library ships its full documentation inside the wheel
(`pyrestoolbox.docs_dir()`), and native code execution is more flexible than
any tool surface.

## Design

A small meta-tool surface instead of one MCP tool per library function:

| Tool | Purpose |
|---|---|
| `list_functions(module?)` | Module map, or public functions of one module |
| `describe(name)` | Signature + docstring for `module.function` |
| `call(name, arguments)` | Invoke any public function with keyword args |
| `recommend_methods(...)` | Correlation-method recommendations with rationale |
| `co2_brine_props(...)` | One-shot CO2-saturated brine properties |
| `sw_brine_props(...)` | One-shot multicomponent gas-saturated brine properties |
| `fit_and_forecast(...)` | Fit a decline model and forecast in one call |

Every module's RST documentation (shipped in the pyrestoolbox wheel) is
exposed as an MCP resource: `docs://gas`, `docs://oil`, `docs://index`, ...

Correlation methods are passed as strings (`"zmethod": "DAK"`); an invalid
method string returns an error listing the valid options, so agents
self-correct in one round trip. Results are JSON: numpy arrays become lists,
DataFrames become `{columns, records}`, non-finite floats become null.

## Install and configure

```
pip install pyrestoolbox-mcp        # once published; or: pip install ./mcp
```

Claude Code:

```
claude mcp add pyrestoolbox -- pyrestoolbox-mcp
```

Claude Desktop (`claude_desktop_config.json`):

```json
{
  "mcpServers": {
    "pyrestoolbox": {
      "command": "pyrestoolbox-mcp"
    }
  }
}
```

## Remote deployment (enterprise LLM platforms)

The server also runs as a remote MCP server over streamable HTTP:

```
pyrestoolbox-mcp --transport streamable-http --host 0.0.0.0 --port 8000
```

It is stateless (no session affinity needed behind a load balancer) and
serves the MCP endpoint at `/mcp`. A `Dockerfile` is included; put an API
gateway in front for authentication.

## Tests

```
PYTHONPATH=<repo>:<repo>/mcp python3 -m pytest mcp/tests/ -q
```

## Notes

- Units default to oilfield (psia, deg F, ft, mD, cP); pass `metric: true`
  where the function supports it. The plyasunov module is SI (K, MPa).
- The oil module is scalar-only; gas and brine functions accept lists.
- Gas rates are Mscf/d everywhere (sm3/d with `metric: true`).
- Functions that take objects are callable through `call` with those objects
  written as JSON objects of their constructor keywords: `completion`
  (`nodal.Completion`, with an optional `segments` list of `WellSegment`
  objects), `reservoir`, `gas_pvt`, `oil_pvt`, `result` (the object a
  `fit_decline` or `fit_ratio` call returned), `ratios`, and `func` as a
  `'module.function'` name for `sensitivity.sweep`/`tornado`. For example
  `call('nodal.fbhp', {"thp": 500, "well_type": "gas", "qg_mscfd": 5000,
  "completion": {"tid": 2.441, "length": 10000, "tht": 100, "bht": 200}})`.
  Classes themselves are not called directly; the one-shot wrapper tools
  cover the common stateful workflows.
- File-writing simtools functions keep their `export=False` defaults, so
  calls return table text rather than touching disk.
