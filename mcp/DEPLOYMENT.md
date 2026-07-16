# pyrestoolbox-mcp: deployment note for the Mawson platform team

## What this is

An MCP (Model Context Protocol) server exposing pyResToolbox, a petroleum
engineering calculation library (gas/oil/brine PVT, nodal analysis, decline
curve analysis, material balance). Source: https://github.com/mwburgoyne/pyResToolbox
(server code in `mcp/`). Contact: Mark Burgoyne.

It gives an LLM seven tools: three discovery/dispatch tools that reach the
library's full public API (about 125 functions), a correlation-method
recommender, and three one-shot calculation wrappers. Module documentation
with worked examples is served as MCP resources. Wrong inputs return errors
that list the valid options, so the model self-corrects without help.

## Security properties

- Pure computation. The server makes no outbound network calls, reads no
  external data sources, and writes nothing to disk. Requests contain only
  numeric/string calculation inputs; responses are JSON numbers and text.
- No credentials, no state, no user data. Nothing is retained between
  requests.
- Licence: GPL-3.0. Internal use is not distribution, so no obligations are
  triggered, and the source is public in any case.

## How to deploy

1. Build the container from the repo's `mcp/` directory (Dockerfile
   included; pyresToolbox version is pinned via a build arg). No base-image
   requirements beyond python:3.11-slim.
2. Run it anywhere inside the Santos network. It listens on port 8000 and
   serves MCP streamable HTTP at `/mcp`. It is stateless: replicas can sit
   behind a load balancer without session affinity.
3. Put your standard API gateway in front for authentication. The server
   itself does not authenticate; it expects to live behind the gateway.
4. Register the endpoint URL with Mawson as a remote MCP server.

If Mawson does not support remote MCP servers and needs an OpenAPI tool
spec instead, tell us; the tool layer is plain typed Python and a REST
wrapper is straightforward.

## Sizing

Single-digit milliseconds to low seconds per calculation, CPU-bound, modest
memory (numpy/scipy stack, roughly 500 MB resident). One small container is
enough to start; scale horizontally if usage grows.

## Updates

Version bumps are deliberate: rebuild the image with a new
`PYRESTOOLBOX_VERSION` build arg after the library publishes a release.
The calculation library is regression-tested (900+ tests, frozen numerical
baselines, documentation examples verified in CI).
