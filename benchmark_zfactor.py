"""Benchmark: Halley (Z*=Z-B) vs Cardano analytical solver for PR-EOS cubic.

Tests both pure Python and Rust-accelerated paths.
"""

import time
import numpy as np
import sys
sys.path.insert(0, '/home/mark/projects')

from pyrestoolbox.gas.gas import _cardano_cubic, _halley_cubic_vec

# ---------- Setup: generate realistic cubic coefficients ----------
# Use BNS-style PR-EOS coefficients for a range of pressures
# Simulate a typical gas: sg=0.7, T=200F, pressures 14.7 to 15000 psia
np.random.seed(42)

N_POINTS = 10_000  # Number of pressure points to solve

# Realistic ranges for A and B from PR-EOS mixing rules
# A ~ 0.01 to 5.0, B ~ 0.001 to 0.5
A = np.random.uniform(0.01, 5.0, N_POINTS)
B = np.random.uniform(0.005, 0.5, N_POINTS)

# PR cubic coefficients: Z^3 - (1-B)*Z^2 + (A-3B^2-2B)*Z - (AB-B^2-B^3) = 0
c2 = -(1.0 - B)
c1 = A - 3.0 * B**2 - 2.0 * B
c0 = -(A * B - B**2 - B**3)

# ---------- Pure Python benchmarks ----------
N_REPEATS = 20

print(f"Benchmarking Z-factor solvers: {N_POINTS} cubic equations, {N_REPEATS} repeats")
print("=" * 70)

# --- Halley (vectorized, Z* reformulation) ---
times_halley = []
for _ in range(N_REPEATS):
    t0 = time.perf_counter()
    Z_halley = _halley_cubic_vec(c2.copy(), c1.copy(), c0.copy(), A=A, B=B)
    t1 = time.perf_counter()
    times_halley.append(t1 - t0)

halley_median = np.median(times_halley) * 1000
halley_min = np.min(times_halley) * 1000

# --- Cardano (element-by-element, analytical) ---
times_cardano = []
for _ in range(N_REPEATS):
    t0 = time.perf_counter()
    Z_cardano = np.array([_cardano_cubic(c2[i], c1[i], c0[i], flag=1) for i in range(N_POINTS)])
    t1 = time.perf_counter()
    times_cardano.append(t1 - t0)

cardano_median = np.median(times_cardano) * 1000
cardano_min = np.min(times_cardano) * 1000

# --- Halley without Z* (legacy path, for comparison) ---
times_halley_legacy = []
for _ in range(N_REPEATS):
    t0 = time.perf_counter()
    Z_halley_leg = _halley_cubic_vec(c2.copy(), c1.copy(), c0.copy())
    t1 = time.perf_counter()
    times_halley_legacy.append(t1 - t0)

halley_leg_median = np.median(times_halley_legacy) * 1000
halley_leg_min = np.min(times_halley_legacy) * 1000

print("\n--- Pure Python Results ---")
print(f"{'Method':<30} {'Median (ms)':>12} {'Min (ms)':>12} {'Speedup':>10}")
print("-" * 70)
print(f"{'Halley Z* (vectorized)':<30} {halley_median:>12.2f} {halley_min:>12.2f} {'baseline':>10}")
print(f"{'Halley legacy (vectorized)':<30} {halley_leg_median:>12.2f} {halley_leg_min:>12.2f} {halley_min/halley_leg_min:>10.2f}x")
print(f"{'Cardano (scalar loop)':<30} {cardano_median:>12.2f} {cardano_min:>12.2f} {halley_min/cardano_min:>10.2f}x")

# Verify agreement
max_diff = np.max(np.abs(Z_halley - Z_cardano))
print(f"\nMax |Z_halley - Z_cardano|: {max_diff:.2e}")

# Count Halley fallbacks (where Halley needed Cardano backup)
d2 = 4.0 * B - 1.0
d1_star = A + 2.0 * B * (B - 2.0)
d0_star = -2.0 * B * B
Zs_check = Z_halley - B
f_check = Zs_check**3 + d2 * Zs_check**2 + d1_star * Zs_check + d0_star
n_fallback = np.sum((np.abs(f_check) > 1e-6) | (Zs_check <= 0))
print(f"Halley Z* fallbacks to Cardano: {n_fallback}/{N_POINTS}")

# ---------- Rust benchmarks ----------
print("\n" + "=" * 70)
try:
    from pyrestoolbox.gas.gas import _rust_module, RUST_AVAILABLE
    if not RUST_AVAILABLE:
        raise ImportError("Rust module not available")
    print("Rust module loaded successfully")

    # For Rust, we need to go through the full BNS path which includes
    # mixing rules. Let's benchmark the full gas_z call instead.
    from pyrestoolbox import gas as gmod

    # Test parameters
    sg = 0.7
    degf = 200.0
    co2, h2s, n2, h2 = 0.01, 0.0, 0.02, 0.0
    pressures = np.linspace(14.7, 15000, N_POINTS).tolist()

    # --- Rust BNS path ---
    # Force Rust
    times_rust = []
    for _ in range(N_REPEATS):
        t0 = time.perf_counter()
        z_rust = np.array([_rust_module.bns_zfactor_full(
            float(pi), float(degf), float(sg),
            float(co2), float(h2s), float(n2), float(h2)
        ) for pi in pressures])
        t1 = time.perf_counter()
        times_rust.append(t1 - t0)

    rust_median = np.median(times_rust) * 1000
    rust_min = np.min(times_rust) * 1000

    # --- Python BNS path (force by temporarily hiding Rust) ---
    import pyrestoolbox.gas.gas as gas_mod
    orig_rust = gas_mod.RUST_AVAILABLE
    gas_mod.RUST_AVAILABLE = False

    times_python_full = []
    for _ in range(N_REPEATS):
        t0 = time.perf_counter()
        z_python = gmod.gas_z(pressures, sg, degf, zmethod='BNS',
                              co2=co2, h2s=h2s, n2=n2, h2=h2)
        t1 = time.perf_counter()
        times_python_full.append(t1 - t0)

    gas_mod.RUST_AVAILABLE = orig_rust

    python_full_median = np.median(times_python_full) * 1000
    python_full_min = np.min(times_python_full) * 1000

    # --- Rust with Cardano-only path ---
    # Check if Rust module has a cardano-only entry point
    has_cardano_only = hasattr(_rust_module, 'bns_zfactor_full_cardano')

    print("\n--- Full BNS Z-Factor Pipeline (mixing rules + solve + fugacity) ---")
    print(f"{'Method':<30} {'Median (ms)':>12} {'Min (ms)':>12} {'Speedup':>10}")
    print("-" * 70)
    print(f"{'Python Halley Z*':<30} {python_full_median:>12.2f} {python_full_min:>12.2f} {'baseline':>10}")
    print(f"{'Rust Halley Z* (scalar loop)':<30} {rust_median:>12.2f} {rust_min:>12.2f} {python_full_min/rust_min:>10.1f}x")

    max_diff_full = np.max(np.abs(np.array(z_rust) - np.array(z_python)))
    print(f"\nMax |Z_rust - Z_python|: {max_diff_full:.2e}")

    if has_cardano_only:
        times_rust_cardano = []
        for _ in range(N_REPEATS):
            t0 = time.perf_counter()
            z_rc = np.array([_rust_module.bns_zfactor_full_cardano(
                float(pi), float(degf), float(sg),
                float(co2), float(h2s), float(n2), float(h2)
            ) for pi in pressures])
            t1 = time.perf_counter()
            times_rust_cardano.append(t1 - t0)
        rc_median = np.median(times_rust_cardano) * 1000
        rc_min = np.min(times_rust_cardano) * 1000
        print(f"{'Rust Cardano (scalar loop)':<30} {rc_median:>12.2f} {rc_min:>12.2f} {python_full_min/rc_min:>10.1f}x")

except ImportError as e:
    print(f"Rust module not available: {e}")
    print("Skipping Rust benchmarks.")

# ---------- Iteration count analysis ----------
print("\n" + "=" * 70)
print("--- Halley Z* Iteration Count Analysis ---")

# Count how many iterations Halley needs
iter_counts = np.zeros(N_POINTS, dtype=int)
d2 = 4.0 * B - 1.0
d1_star = A + 2.0 * B * (B - 2.0)
d0_star = -2.0 * B * B

Zs = np.ones(N_POINTS)
for it in range(50):
    f = Zs**3 + d2 * Zs**2 + d1_star * Zs + d0_star
    fp = 3.0 * Zs**2 + 2.0 * d2 * Zs + d1_star
    fpp = 6.0 * Zs + 2.0 * d2
    safe_fp = np.where(np.abs(fp) < 1e-30, 1e-30, fp)
    dZ = f / safe_fp
    denom = safe_fp - 0.5 * dZ * fpp
    denom = np.where(np.abs(denom) < 1e-30, 1e-30, denom)
    dZ = f / denom
    Zs -= dZ

    # Mark converged
    not_converged = np.abs(dZ) >= 1e-12
    iter_counts[not_converged] = it + 1

# The ones that converged at iteration X have iter_counts = X
# The ones still at 0 converged in 1 iteration
iter_counts[iter_counts == 0] = 1

print(f"Mean iterations: {np.mean(iter_counts):.1f}")
print(f"Max iterations:  {np.max(iter_counts)}")
print(f"Distribution:")
for i in range(1, min(np.max(iter_counts) + 1, 15)):
    count = np.sum(iter_counts == i)
    if count > 0:
        print(f"  {i} iterations: {count:>6} ({100*count/N_POINTS:.1f}%)")
