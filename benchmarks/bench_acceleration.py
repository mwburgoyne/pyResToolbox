"""
Benchmark Rust vs Python Z-factor, viscosity, pseudopressure, and critical property implementations.

Run with:
    PYTHONPATH=/home/mark/projects python3 benchmarks/bench_acceleration.py
"""
import subprocess
import sys
import os
import json

PYTHONPATH = "/home/mark/projects"
WORKER = os.path.join(os.path.dirname(__file__), "_bench_worker.py")


def run_bench(mode, env_overrides=None):
    """Run benchmark in a clean subprocess."""
    env = os.environ.copy()
    env["PYTHONPATH"] = PYTHONPATH
    env["BENCH_MODE"] = mode
    if env_overrides:
        env.update(env_overrides)
    result = subprocess.run(
        [sys.executable, WORKER],
        capture_output=True, text=True, env=env,
        timeout=600,
    )
    print(result.stdout)
    if result.returncode != 0:
        print("STDERR:", result.stderr[:1000])
        return {}

    for line in result.stdout.splitlines():
        if line.startswith("RESULTS:"):
            return json.loads(line[8:])
    return {}


print("=" * 70)
print("PHASE 1: Pure Python benchmarks (PYRESTOOLBOX_NO_RUST=1)")
print("=" * 70)
py_results = run_bench("Python", {"PYRESTOOLBOX_NO_RUST": "1"})

print("\n" + "=" * 70)
print("PHASE 2: Rust-accelerated benchmarks")
print("=" * 70)
rs_results = run_bench("Rust")

# Summary
if py_results and rs_results:
    print("\n" + "=" * 70)
    print("SPEEDUP SUMMARY (via public API)")
    print("=" * 70)
    print(f"{'Method':<35} {'Python':>10} {'Rust':>10} {'Speedup':>10}")
    print("-" * 67)

    labels = {
        "dak_scalar": "DAK Z scalar (10k calls)",
        "hy_scalar": "HY Z scalar (10k calls)",
        "bns_scalar": "BNS Z scalar (10k calls)",
        "dak_array": "DAK Z 50-pt array (200x)",
        "bns_array": "BNS Z 50-pt array (200x)",
        "sut_tc_pc": "Sutton+WA Tc/Pc (10k)",
        "bns_tc_pc": "BNS Tc/Pc (10k)",
        "ug_lge_scalar": "LGE viscosity scalar (10k)",
        "ug_lbc_scalar": "LBC viscosity scalar (10k)",
        "ug_lge_array": "LGE visc 50-pt array (200x)",
        "ug_lbc_array": "LBC visc 50-pt array (200x)",
        "dmp_dak": "Pseudopressure DAK/SUT (1k)",
        "dmp_hy": "Pseudopressure HY/SUT (1k)",
        "dmp_bns": "Pseudopressure BNS (1k)",
        "dmp_bns_h2": "Pseudopressure BNS+H2 (1k)",
        "vlp_hb_gas": "VLP HB Gas (50 calls)",
        "vlp_wg_gas": "VLP WG Gas (50 calls)",
        "vlp_gray_gas": "VLP Gray Gas (50 calls)",
        "vlp_bb_gas": "VLP BB Gas (50 calls)",
        "dca_hyp_time": "DCA Hyp fit time (200 calls)",
        "dca_hyp_cum": "DCA Hyp fit cum (200 calls)",
        "oil_deno": "Oil density SWMH (5k calls)",
        "oil_bo": "Oil Bo McCain (5k calls)",
        "sp_co2_brine": "Spycher-Pruess CO2 (500 calls)",
        "vle_equilibrium": "VLE equilibrium (500 calls)",
    }

    keys_order = [
        "dak_scalar", "hy_scalar", "bns_scalar", "dak_array", "bns_array",
        "sut_tc_pc", "bns_tc_pc",
        "ug_lge_scalar", "ug_lbc_scalar", "ug_lge_array", "ug_lbc_array",
        "dmp_dak", "dmp_hy", "dmp_bns", "dmp_bns_h2",
        "vlp_hb_gas", "vlp_wg_gas", "vlp_gray_gas", "vlp_bb_gas",
        "dca_hyp_time", "dca_hyp_cum",
        "oil_deno", "oil_bo",
        "sp_co2_brine", "vle_equilibrium",
    ]

    for key in keys_order:
        if key in py_results and key in rs_results:
            py_t = py_results[key]
            rs_t = rs_results[key]
            speedup = py_t / rs_t
            print(f"{labels[key]:<35} {py_t:>9.4f}s {rs_t:>9.4f}s {speedup:>9.1f}x")
else:
    print("\nCould not compute summary - benchmark results missing.")
