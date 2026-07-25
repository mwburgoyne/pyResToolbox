#!/usr/bin/env python3
"""
BNS Viscosity Benchmark: Rust-accelerated vs Pure Python

Replicates the benchmarking approach from Peter Kirkham's analysis at
https://pkirkham.github.io/analysis/viscosity-bns/

Runs randomised BNS viscosity calculations in 5-second windows, comparing
Rust-accelerated pyResToolbox against the pure-Python fallback. Reports
calculations per second for single and batch (1000-pressure) modes, and
produces a logarithmic bar chart.
"""

import json
import os
import subprocess
import sys
import time

# ── Benchmark code executed inside subprocesses ──────────────────────────

BENCHMARK_CODE = r'''
import time
import random
import numpy as np
from pyrestoolbox import gas
from pyrestoolbox.classes import z_method, c_method

P_SC = 14.7    # psia
T_SC = 60.0    # degF
DURATION = 5.0 # seconds per benchmark window

def random_between(a, b):
    return a + (b - a) * random.random()

def run_single(duration=DURATION):
    """Single-pressure calculations, one at a time."""
    count = 0
    start = time.time()
    while time.time() - start < duration:
        p = random_between(P_SC, 15000.0)
        t = random_between(T_SC, 300.0)
        sg = random_between(0.65, 0.9)
        if random.random() < 0.3:
            inert_frac = random.random()
            co2 = random_between(0.0, inert_frac)
            h2s = random_between(0.0, 1.0 - co2)
            n2  = random_between(0.0, 1.0 - co2 - h2s)
            h2  = random_between(0.0, 1.0 - co2 - h2s - n2)
        else:
            co2 = h2s = n2 = h2 = 0.0

        mu = gas.gas_ug(p=p, sg=sg, degf=t,
                        zmethod=z_method.BNS, cmethod=c_method.BNS,
                        co2=co2, h2s=h2s, n2=n2, h2=h2)
        assert np.all(np.isfinite(mu)) and np.all(mu > 0)
        count += 1

    elapsed = time.time() - start
    return count / elapsed

def run_batch(n_samples=1000, duration=DURATION):
    """Batch-pressure calculations, n_samples pressures per iteration."""
    rng = np.random.default_rng()
    count = 0
    start = time.time()
    while time.time() - start < duration:
        p = rng.uniform(P_SC, 15000.0, size=n_samples)
        t = random_between(T_SC, 300.0)
        sg = random_between(0.65, 0.9)
        if random.random() < 0.3:
            inert_frac = random.random()
            co2 = random_between(0.0, inert_frac)
            h2s = random_between(0.0, 1.0 - co2)
            n2  = random_between(0.0, 1.0 - co2 - h2s)
            h2  = random_between(0.0, 1.0 - co2 - h2s - n2)
        else:
            co2 = h2s = n2 = h2 = 0.0

        mu = gas.gas_ug(p=p, sg=sg, degf=t,
                        zmethod=z_method.BNS, cmethod=c_method.BNS,
                        co2=co2, h2s=h2s, n2=n2, h2=h2)
        assert np.all(np.isfinite(mu)) and np.all(mu > 0)
        count += n_samples

    elapsed = time.time() - start
    return count / elapsed

import json
single_rate = run_single()
batch_rate = run_batch()
print(json.dumps({"single": single_rate, "batch": batch_rate}))
'''


def run_subprocess(env_override=None):
    """Run the benchmark in a subprocess, return dict with single/batch rates."""
    env = os.environ.copy()
    if env_override:
        env.update(env_override)

    code = f"import sys, os\nsys.path.insert(0, {os.getcwd()!r})\n{BENCHMARK_CODE}"
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True, text=True, env=env, timeout=120
    )
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr.strip()[:200]}")
        return None
    try:
        return json.loads(result.stdout.strip().split("\n")[-1])
    except (json.JSONDecodeError, IndexError):
        print(f"  Bad output: {result.stdout[:200]}")
        return None


def main():
    print("=" * 70)
    print("  BNS Viscosity Benchmark: Rust vs Pure Python")
    print("  (after P. Kirkham, https://pkirkham.github.io/analysis/viscosity-bns/)")
    print("=" * 70)
    print()
    print("  Parameters:")
    print("    Pressure:    14.7 – 15,000 psia (uniform random)")
    print("    Temperature: 60 – 300 °F (uniform random)")
    print("    SG:          0.65 – 0.9 (uniform random)")
    print("    Inerts:      30% chance of CO2/H2S/N2/H2")
    print("    Duration:    5 s per test window")
    print("    Batch size:  1,000 pressures per call")
    print()

    # ── Rust-accelerated run ─────────────────────────────────────────────
    print("  Running Rust-accelerated benchmark...", flush=True)
    rust = run_subprocess()
    if rust is None:
        print("  Rust benchmark failed!")
        return

    print(f"    Single: {rust['single']:>10,.0f} calc/s")
    print(f"    Batch:  {rust['batch']:>10,.0f} calc/s")
    print()

    # ── Pure Python run ──────────────────────────────────────────────────
    print("  Running pure-Python benchmark...", flush=True)
    python = run_subprocess(env_override={"PYRESTOOLBOX_NO_RUST": "1"})
    if python is None:
        print("  Python benchmark failed!")
        return

    print(f"    Single: {python['single']:>10,.0f} calc/s")
    print(f"    Batch:  {python['batch']:>10,.0f} calc/s")
    print()

    # ── Summary table ────────────────────────────────────────────────────
    speedup_single = rust["single"] / python["single"] if python["single"] > 0 else 0
    speedup_batch = rust["batch"] / python["batch"] if python["batch"] > 0 else 0

    print("  ┌─────────────────┬──────────────┬──────────────┬──────────┐")
    print("  │ Mode            │ Rust (calc/s)│Python (calc/s)│ Speedup  │")
    print("  ├─────────────────┼──────────────┼──────────────┼──────────┤")
    print(f"  │ Single pressure │ {rust['single']:>12,.0f}│{python['single']:>14,.0f}│ {speedup_single:>7.1f}x │")
    print(f"  │ 1000x batch     │ {rust['batch']:>12,.0f}│{python['batch']:>14,.0f}│ {speedup_batch:>7.1f}x │")
    print("  └─────────────────┴──────────────┴──────────────┴──────────┘")
    print()

    # ── Bar chart ────────────────────────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        labels = ["Python\nSingle", "Rust\nSingle", "Python\nBatch (1000)", "Rust\nBatch (1000)"]
        rates = [python["single"], rust["single"], python["batch"], rust["batch"]]
        colors = ["#4878CF", "#D65F5F", "#4878CF", "#D65F5F"]
        hatches = ["", "", "//", "//"]

        fig, ax = plt.subplots(figsize=(8, 5))
        bars = ax.bar(labels, rates, color=colors, edgecolor="black", linewidth=0.8)
        for bar, h in zip(bars, hatches):
            bar.set_hatch(h)

        ax.set_yscale("log")
        ax.set_ylabel("BNS Viscosity Calculations per Second")
        ax.set_title("BNS Viscosity Benchmark: Rust-accelerated vs Pure Python\n"
                      "(P: 14.7–15,000 psia, T: 60–300 °F, SG: 0.65–0.9, 30% inerts)")
        ax.grid(axis="y", which="both", alpha=0.3)

        # Annotate bars with values
        for bar, rate in zip(bars, rates):
            ax.text(bar.get_x() + bar.get_width() / 2, rate * 1.3,
                    f"{rate:,.0f}/s", ha="center", va="bottom", fontsize=10, fontweight="bold")

        # Add speedup annotations
        ax.annotate(f"{speedup_single:.1f}x", xy=(1, rust["single"]),
                    xytext=(0.5, rust["single"] * 3),
                    arrowprops=dict(arrowstyle="->", color="green", lw=1.5),
                    fontsize=11, color="green", fontweight="bold", ha="center")
        ax.annotate(f"{speedup_batch:.1f}x", xy=(3, rust["batch"]),
                    xytext=(2.5, rust["batch"] * 3),
                    arrowprops=dict(arrowstyle="->", color="green", lw=1.5),
                    fontsize=11, color="green", fontweight="bold", ha="center")

        plt.tight_layout()
        outpath = os.path.join(os.path.dirname(__file__) or ".", "benchmark_viscosity.png")
        plt.savefig(outpath, dpi=150)
        print(f"  Chart saved to: {outpath}")
    except ImportError:
        print("  (matplotlib not available — skipping chart)")


if __name__ == "__main__":
    main()
