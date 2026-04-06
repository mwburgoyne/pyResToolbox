#!/usr/bin/env python3
"""
Benchmark: Rust-accelerated vs pure-Python pyResToolbox functions.

Runs 8 test cases with BNS method, comparing results, accuracy and speedup.
Pure-Python runs in a subprocess with PYRESTOOLBOX_NO_RUST=1 (set before
import since the accelerator probes at import time).
"""

import json
import os
import subprocess
import sys
import time

import numpy as np
from tabulate import tabulate


# ── Helper: run a test case in a subprocess ──────────────────────────────
def _run_in_subprocess(test_code: str, env_override: dict | None = None) -> dict:
    """Execute test_code in a fresh Python process and return JSON result."""
    env = os.environ.copy()
    if env_override:
        env.update(env_override)

    wrapper = f"""
import json, time, sys, os
sys.path.insert(0, {os.getcwd()!r})
{test_code}
"""
    result = subprocess.run(
        [sys.executable, "-c", wrapper],
        capture_output=True, text=True, env=env, timeout=300
    )
    if result.returncode != 0:
        return {"error": result.stderr.strip(), "elapsed": 0, "value": None}
    try:
        return json.loads(result.stdout.strip().split("\n")[-1])
    except (json.JSONDecodeError, IndexError):
        return {"error": f"Bad output: {result.stdout[:200]}", "elapsed": 0, "value": None}


# ── Test case definitions ────────────────────────────────────────────────

PREAMBLE = """
from pyrestoolbox import gas, oil
from pyrestoolbox.nodal import Completion, fbhp
from pyrestoolbox import simtools
from pyrestoolbox.classes import z_method, c_method
import numpy as np
import time
import json
"""

# (label, function_call_string, test_code)
TESTS = [
    (
        "(a) Single Z-factor",
        "gas.gas_z(p=3000, sg=0.7, degf=200, zmethod=z_method.BNS, cmethod=c_method.BNS)",
        """
z = gas.gas_z(p=3000, sg=0.7, degf=200, zmethod=z_method.BNS, cmethod=c_method.BNS)
print(json.dumps({"elapsed": _elapsed, "value": float(z)}))
""",
    ),
    (
        "(b) 100x Z-factors",
        "[gas.gas_z(p=p, sg=0.7, degf=200, zmethod=z_method.BNS, cmethod=c_method.BNS)\n  for p in np.linspace(500, 5000, 100)]",
        """
pressures = np.linspace(500, 5000, 100).tolist()
results = [gas.gas_z(p=p, sg=0.7, degf=200, zmethod=z_method.BNS, cmethod=c_method.BNS) for p in pressures]
print(json.dumps({"elapsed": _elapsed, "value": float(np.mean(results))}))
""",
    ),
    (
        "(c) Single viscosity",
        "gas.gas_ug(p=3000, sg=0.7, degf=200, zmethod=z_method.BNS, cmethod=c_method.BNS)",
        """
ug = gas.gas_ug(p=3000, sg=0.7, degf=200, zmethod=z_method.BNS, cmethod=c_method.BNS)
print(json.dumps({"elapsed": _elapsed, "value": float(ug)}))
""",
    ),
    (
        "(d) 100x viscosities",
        "[gas.gas_ug(p=p, sg=0.7, degf=200, zmethod=z_method.BNS, cmethod=c_method.BNS)\n  for p in np.linspace(500, 5000, 100)]",
        """
pressures = np.linspace(500, 5000, 100).tolist()
results = [gas.gas_ug(p=p, sg=0.7, degf=200, zmethod=z_method.BNS, cmethod=c_method.BNS) for p in pressures]
print(json.dumps({"elapsed": _elapsed, "value": float(np.mean(results))}))
""",
    ),
    (
        "(e) Pseudopressure",
        "gas.gas_dmp(p1=500, p2=4000, degf=200, sg=0.7, zmethod=z_method.BNS, cmethod=c_method.BNS)",
        """
pp = gas.gas_dmp(p1=500, p2=4000, degf=200, sg=0.7, zmethod=z_method.BNS, cmethod=c_method.BNS)
print(json.dumps({"elapsed": _elapsed, "value": float(pp)}))
""",
    ),
    (
        "(f) Gas outflow",
        "fbhp(thp=500, completion=Completion(2.441, 10000, 100, 200),\n  vlpmethod='WG', well_type='gas', qg_mmscfd=5, cgr=10,\n  qw_bwpd=10, api=45, gsg=0.65)",
        """
c = Completion(tid=2.441, length=10000, tht=100, bht=200)
bhp = fbhp(thp=500, completion=c, vlpmethod='WG', well_type='gas',
           qg_mmscfd=5.0, cgr=10, qw_bwpd=10, api=45, gsg=0.65, wsg=1.07)
print(json.dumps({"elapsed": _elapsed, "value": float(bhp)}))
""",
    ),
    (
        "(g) Oil outflow",
        "fbhp(thp=200, completion=Completion(2.441, 8000, 100, 180),\n  vlpmethod='HB', well_type='oil', qt_stbpd=2000, gor=800,\n  wc=0.3, pb=2500, rsb=500, api=35, sgsp=0.65, gsg=0.65)",
        """
c = Completion(tid=2.441, length=8000, tht=100, bht=180)
bhp = fbhp(thp=200, completion=c, vlpmethod='HB', well_type='oil',
           qt_stbpd=2000, gor=800, wc=0.3, pb=2500, rsb=500,
           api=35, sgsp=0.65, gsg=0.65)
print(json.dumps({"elapsed": _elapsed, "value": float(bhp)}))
""",
    ),
    (
        "(h) CT Influence table",
        "simtools.influence_tables(ReDs=[2, 5, 10],\n  min_td=0.01, max_td=10, n_incr=10, M=7)",
        """
tds, pds = simtools.influence_tables(ReDs=[2, 5, 10], min_td=0.01, max_td=10, n_incr=10, M=7, export=False)
val = sum(sum(pd) for pd in pds)
print(json.dumps({"elapsed": _elapsed, "value": float(val)}))
""",
    ),
]


def _wrap_timed(test_code: str) -> str:
    """Wrap test code with timing. Replaces _elapsed placeholder."""
    lines = test_code.strip().split("\n")
    print_line = lines[-1]
    setup_lines = "\n".join(lines[:-1])
    return f"""
{PREAMBLE}
t0 = time.perf_counter()
{setup_lines}
_elapsed = time.perf_counter() - t0
{print_line}
"""


def main():
    print("=" * 78)
    print("  pyResToolbox Benchmark: Rust-accelerated vs Pure Python (BNS)")
    print("=" * 78)
    print()

    # Confirm Rust status
    rust_check = _run_in_subprocess(
        PREAMBLE + """
from pyrestoolbox._accelerator import get_status
s = get_status()
print(json.dumps({"elapsed": 0, "value": 1 if s['rust_available'] else 0}))
"""
    )
    rust_avail = rust_check.get("value", 0) == 1
    print(f"  Rust acceleration available: {rust_avail}")
    if not rust_avail:
        print("  WARNING: Rust not available - benchmark will show Python vs Python")
    print()

    rows = []
    for label, call_str, test_code in TESTS:
        full_code = _wrap_timed(test_code)

        # Run with Rust (default)
        print(f"  Running {label} (Rust)...", end="", flush=True)
        r_rust = _run_in_subprocess(full_code)
        if r_rust.get("error"):
            print(f" ERROR: {r_rust['error'][:80]}")
            rows.append([label, call_str, "ERROR", "", "", "", "", ""])
            continue
        print(f" {r_rust['elapsed']:.4f}s", flush=True)

        # Run pure Python
        print(f"  Running {label} (Python)...", end="", flush=True)
        r_py = _run_in_subprocess(full_code, env_override={"PYRESTOOLBOX_NO_RUST": "1"})
        if r_py.get("error"):
            print(f" ERROR: {r_py['error'][:80]}")
            rows.append([label, call_str, f"{r_rust['elapsed']:.4f}", "ERROR",
                         f"{r_rust['value']:.6g}", "", "", ""])
            continue
        print(f" {r_py['elapsed']:.4f}s", flush=True)

        speedup = r_py["elapsed"] / r_rust["elapsed"] if r_rust["elapsed"] > 1e-9 else float("inf")
        v_rust, v_py = r_rust["value"], r_py["value"]
        if v_rust is not None and v_py is not None and abs(v_py) > 1e-30:
            rel_err = abs(v_rust - v_py) / abs(v_py)
            accuracy = "Exact" if rel_err < 1e-12 else f"{rel_err:.2e}"
        else:
            accuracy = "N/A"

        rows.append([
            label,
            call_str,
            f"{r_rust['elapsed']:.4f}",
            f"{r_py['elapsed']:.4f}",
            f"{v_rust:.6g}",
            f"{v_py:.6g}",
            f"{speedup:.1f}x",
            accuracy,
        ])

    print()
    print(tabulate(
        rows,
        headers=["Test", "Function Call", "Rust (s)", "Python (s)",
                 "Rust Result", "Python Result", "Speedup", "Rel Error"],
        tablefmt="pretty",
        colalign=("left", "left", "right", "right", "right", "right", "right", "center"),
    ))


if __name__ == "__main__":
    main()
