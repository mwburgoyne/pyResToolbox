"""Smoke-test an installed pyrestoolbox-mcp over stdio with raw JSON-RPC.

Usage: python mcp_stdio_smoke.py <command> [args...]
Independent of the mcp SDK client API, so it checks the server the way any
MCP host would: initialize, tools/list, a search, a library call, and a bad
call whose error message must reach the client intact.
"""
import json
import subprocess
import sys

EXPECTED_TOOLS = {"search", "read_docs", "list_functions", "describe", "call",
                  "recommend_methods", "co2_brine_props", "sw_brine_props", "fit_and_forecast"}


def main():
    proc = subprocess.Popen(sys.argv[1:], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            text=True, encoding="utf-8")

    def send(msg):
        proc.stdin.write(json.dumps(msg) + "\n")
        proc.stdin.flush()

    def reply(msg_id):
        for line in proc.stdout:
            msg = json.loads(line)
            if msg.get("id") == msg_id:
                if "error" in msg:
                    raise SystemExit(f"error reply: {msg['error']}")
                return msg["result"]
        raise SystemExit(f"server exited before replying to id {msg_id}")

    try:
        send({"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {
            "protocolVersion": "2025-06-18", "capabilities": {},
            "clientInfo": {"name": "smoke", "version": "0"}}})
        init = reply(1)
        send({"jsonrpc": "2.0", "method": "notifications/initialized"})
        send({"jsonrpc": "2.0", "id": 2, "method": "tools/list"})
        names = {t["name"] for t in reply(2)["tools"]}
        missing = EXPECTED_TOOLS - names
        if missing:
            raise SystemExit(f"tools missing: {sorted(missing)}")
        send({"jsonrpc": "2.0", "id": 3, "method": "tools/call", "params": {
            "name": "call", "arguments": {"name": "gas.gas_z",
                                          "arguments": {"p": 2000, "sg": 0.75, "degf": 180}}}})
        res = reply(3)
        if res.get("isError"):
            raise SystemExit(f"tools/call failed: {res}")
        text = res["content"][0]["text"]
        payload = json.loads(text)
        send({"jsonrpc": "2.0", "id": 4, "method": "tools/call", "params": {
            "name": "search", "arguments": {"query": "liquid loading"}}})
        hits = json.loads(reply(4)["content"][0]["text"])["matches"]
        if hits[0]["name"] != "nodal.critical_gas_rate":
            raise SystemExit(f"search returned {hits[:3]}")
        send({"jsonrpc": "2.0", "id": 5, "method": "tools/call", "params": {
            "name": "call", "arguments": {"name": "gas.gas_z", "arguments": {
                "p": 2000, "sg": 0.75, "degf": 180, "zmethod": "FOO"}}}})
        bad = reply(5)
        if not bad.get("isError") or "Valid options" not in bad["content"][0]["text"]:
            raise SystemExit(f"error message did not reach the client: {bad}")
        print(f"server {init['serverInfo']}; {len(names)} tools; gas_z -> {payload}; "
              f"search ok; error passthrough ok")
    finally:
        proc.stdin.close()
        proc.terminate()
        proc.wait(timeout=10)


if __name__ == "__main__":
    main()
