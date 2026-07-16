# pyrestoolbox-mcp: remote MCP server (streamable HTTP)
# Build from the mcp/ directory:  docker build -t pyrestoolbox-mcp .
# Run:                            docker run -p 8000:8000 pyrestoolbox-mcp
# MCP endpoint:                   http://<host>:8000/mcp

FROM python:3.11-slim

# Pin the library version for reproducible corporate deployments; bump
# deliberately and rebuild.
ARG PYRESTOOLBOX_VERSION=3.7.0

WORKDIR /app
COPY . /app/mcp

RUN pip install --no-cache-dir "pyrestoolbox==${PYRESTOOLBOX_VERSION}" /app/mcp

# Non-root runtime user
RUN useradd --create-home --shell /usr/sbin/nologin mcp
USER mcp

EXPOSE 8000

CMD ["pyrestoolbox-mcp", "--transport", "streamable-http", "--host", "0.0.0.0", "--port", "8000"]
