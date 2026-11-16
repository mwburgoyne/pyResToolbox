# pyrestoolbox-mcp

A [Model Context Protocol](https://modelcontextprotocol.io) server that gives AI agents (Claude Code, Claude Desktop, Cursor, VS Code Copilot, any MCP host) the whole [pyrestoolbox](https://github.com/mwburgoyne/pyResToolbox) reservoir engineering library: gas, oil and brine PVT, VLP/IPR and nodal analysis, chokes and liquid loading, deliverability tests, gas pipelines, decline curves, material balance, capacitance-resistance models, simulator tables and more.

Requires Python 3.10 or later; installs pyrestoolbox 3.10.1 or later and the official `mcp` SDK (1.22 or later, 1.x and 2.x both supported).

## Set up

With [uv](https://docs.astral.sh/uv/) installed there is nothing to install first; `uvx` fetches and caches the server.

Claude Code:

```
claude mcp add pyrestoolbox -- uvx pyrestoolbox-mcp
```

Claude Desktop (`claude_desktop_config.json`) or Cursor (`.cursor/mcp.json`):

```json
{
  "mcpServers": {
    "pyrestoolbox": {"command": "uvx", "args": ["pyrestoolbox-mcp"]}
  }
}
```

VS Code (`.vscode/mcp.json`):

```json
{
  "servers": {
    "pyrestoolbox": {"command": "uvx", "args": ["pyrestoolbox-mcp"]}
  }
}
```

Without uv: `pip install pyrestoolbox-mcp`, then use `pyrestoolbox-mcp` as the command with no arguments. Desktop apps launched from a dock or start menu often do not see your shell's PATH; if the host reports that it cannot start the server, give the full path to `uvx` or `pyrestoolbox-mcp` (`which uvx`, or `where uvx` on Windows).

To pin the version for repeatable numbers: `uvx pyrestoolbox-mcp@0.1.0` and, if the library must stay fixed too, `uvx --with pyrestoolbox==3.10.1 pyrestoolbox-mcp`.

## What the agent gets

A small set of tools that reach all 120-odd public functions, instead of one tool per function:

| Tool | Purpose |
|---|---|
| `search(query)` | Engineering task in plain words ("liquid loading", "hydrate inhibitor dose", "gas in place") to the matching functions, best first |
| `describe(name)` | Signature and full docstring, with units, for `module.function` |
| `call(name, arguments)` | Run any public function with keyword arguments |
| `read_docs(topic)` | A module's full documentation with worked examples; `index` maps tasks to functions |
| `list_functions(module?)` | Module map, or the functions of one module |
| `recommend_methods(...)` | Which correlations suit a fluid and well, with rationale |
| `co2_brine_props(...)` | CO2-saturated brine in one call (Spycher-Pruess) |
| `sw_brine_props(...)` | Gas-saturated brine for CH4/CO2/H2S/N2/H2/He mixtures (Soreide-Whitson) |
| `fit_and_forecast(...)` | Fit a decline curve and forecast with EUR in one call |

The same documentation is also served as MCP resources (`docs://gas`, `docs://index`, ...).

A typical exchange: the agent calls `search("flowing bottomhole pressure")`, gets `nodal.fbhp`, calls `describe("nodal.fbhp")`, then

```
call("nodal.fbhp", {"thp": 500, "well_type": "gas", "qg_mscfd": 5000,
                    "completion": {"tid": 2.441, "length": 10000, "tht": 100, "bht": 200}})
```

and gets `{"result": 935.8, "engine": {"pyrestoolbox": "3.10.1", "rust": true}}`. Every result names the library version and whether the compiled Rust path or pure Python produced it.

Behaviour an agent can rely on:

- **Errors are readable.** A bad input returns the library's own message, e.g. `Invalid zmethod: 'FOO'. Valid options: ['DAK', 'HY', 'BNS']` or `gas_z() missing 2 required positional arguments: 'sg' and 'degf'`, so the agent corrects itself in one round trip.
- **Caveats travel with numbers.** Library warnings (correlation range, extrapolation) come back in a `notes` list on the result instead of disappearing into a server log.
- **Units.** Oilfield by default (psia, deg F, ft, mD, cP); `"metric": true` switches to Eclipse METRIC (barsa, deg C, m). Parameter *names* keep their oilfield spelling, so `degf` takes deg C under metric. Docstrings write units as `(deg F | deg C)`, field first. Gas rates are Mscf/d (sm3/d metric); standard conditions are 60 deg F and 14.696 psia.
- **Objects as JSON.** Functions that take a `Completion`, `Reservoir`, `GasPVT`, `OilPVT`, a decline result or a `module.function` name accept them as JSON objects of their constructor keywords, as in the `completion` above.
- **Results are JSON.** Arrays become lists, tables become `{columns, records}`, non-finite numbers become null.

## Remote deployment

```
pyrestoolbox-mcp --transport streamable-http --host 0.0.0.0 --port 8000
```

Serves the MCP endpoint at `/mcp`, stateless, so no session affinity is needed behind a load balancer. A `Dockerfile` is included (runs as a non-root user, pins the library version). The server has no authentication of its own: put a gateway in front of it.

One calculation runs at a time per process (the library's warning capture is process-global); scale with more processes or replicas.

## Safety

Every tool is a pure calculation, annotated read-only, idempotent and closed-world, so hosts can auto-approve them.

- **No file access by default.** A few `simtools` functions write Eclipse include files (`export=True`) or read simulator decks and print files. The server refuses those arguments and functions unless started with `--allow-file-io`, which is meant for local use only. Tables always come back in the result.
- **Bounded results.** A result over 500,000 characters of JSON is refused with a message asking for fewer rows or points; change the limit with `--max-result-chars`.
- **Public API only.** `call` reaches public functions of the library's calculation modules; private names, classes and anything outside the library are rejected.

## Tests

From the repository root, with an `mcp` SDK installed:

```
PYTHONPATH=$PWD:$PWD/mcp python -m pytest mcp/tests/ -q
```

## Licence

GPL-3.0-or-later, as for pyrestoolbox.
