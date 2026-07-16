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
gateway in front for authentication. See `DEPLOYMENT.md` for the platform
team one-pager.

## Tests

```
PYTHONPATH=<repo>:<repo>/mcp python3 -m pytest mcp/tests/ -q
```

## Notes

- Units default to oilfield (psia, deg F, ft, mD, cP); pass `metric: true`
  where the function supports it. The plyasunov module is SI (K, MPa).
- The oil module is scalar-only; gas and brine functions accept lists.
- Stateful classes (`GasPVT`, `CO2_Brine_Mixture`, `SoreideWhitson`,
  `DeclineResult`) are not exposed through `call` - the one-shot wrapper
  tools cover the common workflows.
- File-writing simtools functions keep their `export=False` defaults, so
  calls return table text rather than touching disk.
