"""Benchmark: Batch vs scalar Rust dispatch for gas_z and gas_ug."""

import time
import numpy as np
import sys
sys.path.insert(0, '/home/mark/projects')

from pyrestoolbox._native import (
    bns_zfactor_full, bns_zfactor_batch,
    dak_zfactor_full, dak_zfactor_batch,
    hy_zfactor_batch, hall_yarborough_zfactor_full,
    gas_ug_lge, gas_ug_lge_batch,
    gas_ug_lbc, gas_ug_lbc_batch,
)
from pyrestoolbox import gas as gmod
import pyrestoolbox.gas.gas as gas_mod

N_POINTS = 10_000
N_REPEATS = 20
pressures = np.linspace(14.7, 15000, N_POINTS)
p_list = pressures.tolist()

# Common parameters
sg, degf = 0.7, 200.0
co2, h2s, n2, h2 = 0.01, 0.0, 0.02, 0.0

print(f"Batch vs Scalar Rust dispatch: {N_POINTS} pressure points, {N_REPEATS} repeats")
print("=" * 80)

def bench(name, fn, repeats=N_REPEATS):
    times = []
    result = None
    for _ in range(repeats):
        t0 = time.perf_counter()
        result = fn()
        t1 = time.perf_counter()
        times.append(t1 - t0)
    return np.median(times) * 1000, np.min(times) * 1000, result

# ==================== Z-Factor ====================
print("\n--- Z-Factor: BNS Method ---")
_, t_scalar, z_scalar = bench("BNS scalar", lambda: np.array([
    bns_zfactor_full(float(pi), degf, sg, co2, h2s, n2, h2) for pi in pressures
]))
_, t_batch, z_batch = bench("BNS batch", lambda: np.array(
    bns_zfactor_batch(p_list, degf, sg, co2, h2s, n2, h2)
))

# Also benchmark pure Python (Rust disabled)
orig_rust = gas_mod.RUST_AVAILABLE
gas_mod.RUST_AVAILABLE = False
_, t_python, z_python = bench("BNS Python", lambda: gmod.gas_z(
    pressures, sg, degf, zmethod='BNS', co2=co2, h2s=h2s, n2=n2, h2=h2
))
gas_mod.RUST_AVAILABLE = orig_rust

print(f"  {'Method':<30} {'Min (ms)':>10} {'Speedup':>10}")
print(f"  {'-'*55}")
print(f"  {'Python (NumPy vectorized)':<30} {t_python:>10.2f} {'baseline':>10}")
print(f"  {'Rust scalar (N FFI calls)':<30} {t_scalar:>10.2f} {t_python/t_scalar:>10.1f}x")
print(f"  {'Rust batch  (1 FFI call)':<30} {t_batch:>10.2f} {t_python/t_batch:>10.1f}x")
print(f"  Batch vs scalar speedup: {t_scalar/t_batch:.1f}x")
print(f"  Max |batch - scalar|: {np.max(np.abs(z_batch - z_scalar)):.2e}")

# DAK
print("\n--- Z-Factor: DAK Method (SUT critical props) ---")
_, t_dak_scalar, z_dak_s = bench("DAK scalar", lambda: np.array([
    dak_zfactor_full(float(pi), degf, sg, co2, h2s, n2) for pi in pressures
]))
_, t_dak_batch, z_dak_b = bench("DAK batch", lambda: np.array(
    dak_zfactor_batch(p_list, degf, sg, co2, h2s, n2)
))
gas_mod.RUST_AVAILABLE = False
_, t_dak_py, _ = bench("DAK Python", lambda: gmod.gas_z(
    pressures, sg, degf, zmethod='DAK', cmethod='SUT', co2=co2, h2s=h2s, n2=n2
))
gas_mod.RUST_AVAILABLE = orig_rust

print(f"  {'Method':<30} {'Min (ms)':>10} {'Speedup':>10}")
print(f"  {'-'*55}")
print(f"  {'Python (NumPy vectorized)':<30} {t_dak_py:>10.2f} {'baseline':>10}")
print(f"  {'Rust scalar (N FFI calls)':<30} {t_dak_scalar:>10.2f} {t_dak_py/t_dak_scalar:>10.1f}x")
print(f"  {'Rust batch  (1 FFI call)':<30} {t_dak_batch:>10.2f} {t_dak_py/t_dak_batch:>10.1f}x")
print(f"  Batch vs scalar speedup: {t_dak_scalar/t_dak_batch:.1f}x")
print(f"  Max |batch - scalar|: {np.max(np.abs(z_dak_b - z_dak_s)):.2e}")

# HY
print("\n--- Z-Factor: HY Method (SUT critical props) ---")
_, t_hy_scalar, z_hy_s = bench("HY scalar", lambda: np.array([
    hall_yarborough_zfactor_full(float(pi), degf, sg, co2, h2s, n2) for pi in pressures
]))
_, t_hy_batch, z_hy_b = bench("HY batch", lambda: np.array(
    hy_zfactor_batch(p_list, degf, sg, co2, h2s, n2)
))
gas_mod.RUST_AVAILABLE = False
_, t_hy_py, _ = bench("HY Python", lambda: gmod.gas_z(
    pressures, sg, degf, zmethod='HY', cmethod='SUT', co2=co2, h2s=h2s, n2=n2
))
gas_mod.RUST_AVAILABLE = orig_rust

print(f"  {'Method':<30} {'Min (ms)':>10} {'Speedup':>10}")
print(f"  {'-'*55}")
print(f"  {'Python (NumPy vectorized)':<30} {t_hy_py:>10.2f} {'baseline':>10}")
print(f"  {'Rust scalar (N FFI calls)':<30} {t_hy_scalar:>10.2f} {t_hy_py/t_hy_scalar:>10.1f}x")
print(f"  {'Rust batch  (1 FFI call)':<30} {t_hy_batch:>10.2f} {t_hy_py/t_hy_batch:>10.1f}x")
print(f"  Batch vs scalar speedup: {t_hy_scalar/t_hy_batch:.1f}x")
print(f"  Max |batch - scalar|: {np.max(np.abs(z_hy_b - z_hy_s)):.2e}")

# ==================== Viscosity ====================
print("\n--- Viscosity: LGE Method ---")
z_dak = np.array(dak_zfactor_batch(p_list, degf, sg, co2, h2s, n2))
_, t_lge_scalar, ug_lge_s = bench("LGE scalar", lambda: np.array([
    gas_ug_lge(float(pi), sg, degf, float(zi)) for pi, zi in zip(pressures, z_dak)
]))
_, t_lge_batch, ug_lge_b = bench("LGE batch", lambda: np.array(
    gas_ug_lge_batch(p_list, z_dak.tolist(), sg, degf)
))
print(f"  {'Method':<30} {'Min (ms)':>10} {'Speedup':>10}")
print(f"  {'-'*55}")
print(f"  {'Rust scalar (N FFI calls)':<30} {t_lge_scalar:>10.2f} {'baseline':>10}")
print(f"  {'Rust batch  (1 FFI call)':<30} {t_lge_batch:>10.2f} {t_lge_scalar/t_lge_batch:>10.1f}x")
print(f"  Max |batch - scalar|: {np.max(np.abs(ug_lge_b - ug_lge_s)):.2e}")

print("\n--- Viscosity: LBC Method (BNS) ---")
z_bns = np.array(bns_zfactor_batch(p_list, degf, sg, co2, h2s, n2, h2))
_, t_lbc_scalar, ug_lbc_s = bench("LBC scalar", lambda: np.array([
    gas_ug_lbc(float(pi), sg, degf, co2, h2s, n2, h2, float(zi)) for pi, zi in zip(pressures, z_bns)
]))
_, t_lbc_batch, ug_lbc_b = bench("LBC batch", lambda: np.array(
    gas_ug_lbc_batch(p_list, z_bns.tolist(), sg, degf, co2, h2s, n2, h2)
))
print(f"  {'Method':<30} {'Min (ms)':>10} {'Speedup':>10}")
print(f"  {'-'*55}")
print(f"  {'Rust scalar (N FFI calls)':<30} {t_lbc_scalar:>10.2f} {'baseline':>10}")
print(f"  {'Rust batch  (1 FFI call)':<30} {t_lbc_batch:>10.2f} {t_lbc_scalar/t_lbc_batch:>10.1f}x")
print(f"  Max |batch - scalar|: {np.max(np.abs(ug_lbc_b - ug_lbc_s)):.2e}")

# ==================== Full Pipeline (gas_z + gas_ug via public API) ====================
print("\n" + "=" * 80)
print("--- Full Pipeline: gas_z() + gas_ug() via public API ---")

# With batch Rust (current)
gas_mod.RUST_AVAILABLE = True
_, t_full_batch, _ = bench("Batch Rust", lambda: (
    gmod.gas_z(pressures, sg, degf, zmethod='BNS', co2=co2, h2s=h2s, n2=n2, h2=h2),
    gmod.gas_ug(pressures, sg, degf, zmethod='BNS', co2=co2, h2s=h2s, n2=n2, h2=h2),
))

# Pure Python
gas_mod.RUST_AVAILABLE = False
_, t_full_python, _ = bench("Python", lambda: (
    gmod.gas_z(pressures, sg, degf, zmethod='BNS', co2=co2, h2s=h2s, n2=n2, h2=h2),
    gmod.gas_ug(pressures, sg, degf, zmethod='BNS', co2=co2, h2s=h2s, n2=n2, h2=h2),
))
gas_mod.RUST_AVAILABLE = orig_rust

print(f"  {'Method':<30} {'Min (ms)':>10} {'Speedup':>10}")
print(f"  {'-'*55}")
print(f"  {'Python (NumPy vectorized)':<30} {t_full_python:>10.2f} {'baseline':>10}")
print(f"  {'Rust batch (2 FFI calls)':<30} {t_full_batch:>10.2f} {t_full_python/t_full_batch:>10.1f}x")
