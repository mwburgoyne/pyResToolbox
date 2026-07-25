"""Worker script for benchmarks - run in subprocess."""
import sys
import os

# Ensure we use the development version, not the installed one
_dev_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _dev_root)

import time
import json
import numpy as np


def bench(name, func, args, n=100_000):
    for _ in range(100):
        func(*args)
    start = time.perf_counter()
    for _ in range(n):
        func(*args)
    elapsed = time.perf_counter() - start
    rate = n / elapsed
    print(f"  {name}: {elapsed:.4f}s ({rate:,.0f} calls/sec)")
    return elapsed


def bench_array(name, func, args, n=1000):
    for _ in range(10):
        func(*args)
    start = time.perf_counter()
    for _ in range(n):
        func(*args)
    elapsed = time.perf_counter() - start
    rate = n / elapsed
    print(f"  {name}: {elapsed:.4f}s ({rate:,.0f} calls/sec)")
    return elapsed


from pyrestoolbox import gas
from pyrestoolbox._accelerator import RUST_AVAILABLE

mode = os.environ.get("BENCH_MODE", "unknown")
print(f"Mode: {mode}, Rust available: {RUST_AVAILABLE}")

N_SCALAR = 10_000
N_ARRAY = 200
p_array = np.linspace(500, 5000, 50)

results = {}

# ===== Z-Factor benchmarks =====
print(f"\n--- DAK Z-Factor (scalar, {N_SCALAR:,} calls) ---")
results["dak_scalar"] = bench(f"{mode} DAK", gas.gas_z,
    (2000.0, 0.75, 200.0, "DAK", "SUT"), n=N_SCALAR)

print(f"\n--- HY Z-Factor (scalar, {N_SCALAR:,} calls) ---")
results["hy_scalar"] = bench(f"{mode} HY", gas.gas_z,
    (2000.0, 0.75, 200.0, "HY", "SUT"), n=N_SCALAR)

print(f"\n--- BNS Z-Factor (scalar, {N_SCALAR:,} calls) ---")
results["bns_scalar"] = bench(f"{mode} BNS", gas.gas_z,
    (2000.0, 0.75, 200.0, "BNS", "BNS", 0.05, 0.02, 0.01, 0.0), n=N_SCALAR)

print(f"\n--- DAK Z-Factor (50-element array, {N_ARRAY:,} calls) ---")
results["dak_array"] = bench_array(f"{mode} DAK array", gas.gas_z,
    (p_array, 0.75, 200.0, "DAK", "SUT"), n=N_ARRAY)

print(f"\n--- BNS Z-Factor (50-element array, {N_ARRAY:,} calls) ---")
results["bns_array"] = bench_array(f"{mode} BNS array", gas.gas_z,
    (p_array, 0.75, 200.0, "BNS", "BNS", 0.05, 0.02, 0.01, 0.0), n=N_ARRAY)

# ===== Critical Properties benchmarks =====
print(f"\n--- Sutton+WA Critical Properties ({N_SCALAR:,} calls) ---")
results["sut_tc_pc"] = bench(f"{mode} Sutton+WA", gas.gas_tc_pc,
    (0.75, 0.05, 0.02, 0.01, 0.0, "SUT"), n=N_SCALAR)

print(f"\n--- BNS Critical Properties ({N_SCALAR:,} calls) ---")
results["bns_tc_pc"] = bench(f"{mode} BNS Tc/Pc", gas.gas_tc_pc,
    (0.75, 0.05, 0.02, 0.01, 0.0, "BNS"), n=N_SCALAR)

# ===== Gas Viscosity benchmarks =====
print(f"\n--- LGE Viscosity (scalar, {N_SCALAR:,} calls) ---")
results["ug_lge_scalar"] = bench(f"{mode} LGE ug", gas.gas_ug,
    (2000.0, 0.75, 200.0, "DAK", "SUT"), n=N_SCALAR)

print(f"\n--- LBC Viscosity (scalar, {N_SCALAR:,} calls) ---")
results["ug_lbc_scalar"] = bench(f"{mode} LBC ug", gas.gas_ug,
    (2000.0, 0.75, 200.0, "BNS", "BNS", 0.05, 0.02, 0.01, 0.0), n=N_SCALAR)

print(f"\n--- LGE Viscosity (50-element array, {N_ARRAY:,} calls) ---")
results["ug_lge_array"] = bench_array(f"{mode} LGE ug array", gas.gas_ug,
    (p_array, 0.75, 200.0, "DAK", "SUT"), n=N_ARRAY)

print(f"\n--- LBC Viscosity (50-element array, {N_ARRAY:,} calls) ---")
results["ug_lbc_array"] = bench_array(f"{mode} LBC ug array", gas.gas_ug,
    (p_array, 0.75, 200.0, "BNS", "BNS", 0.05, 0.02, 0.01, 0.0), n=N_ARRAY)

# ===== Pseudopressure benchmarks =====
N_DMP = 1_000
print(f"\n--- DAK/SUT Pseudopressure ({N_DMP:,} calls) ---")
results["dmp_dak"] = bench(f"{mode} dmp DAK/SUT", gas.gas_dmp,
    (1000.0, 3000.0, 200.0, 0.75, "DAK", "SUT"), n=N_DMP)

print(f"\n--- HY/SUT Pseudopressure ({N_DMP:,} calls) ---")
results["dmp_hy"] = bench(f"{mode} dmp HY/SUT", gas.gas_dmp,
    (1000.0, 3000.0, 200.0, 0.75, "HY", "SUT"), n=N_DMP)

print(f"\n--- BNS Pseudopressure ({N_DMP:,} calls) ---")
results["dmp_bns"] = bench(f"{mode} dmp BNS", gas.gas_dmp,
    (1000.0, 3000.0, 200.0, 0.75, "BNS", "BNS", 0.05, 0.02, 0.01, 0.0), n=N_DMP)

print(f"\n--- BNS+H2 Pseudopressure ({N_DMP:,} calls) ---")
results["dmp_bns_h2"] = bench(f"{mode} dmp BNS+H2", gas.gas_dmp,
    (1000.0, 4000.0, 180.0, 0.70, "BNS", "BNS", 0.03, 0.01, 0.01, 0.10), n=N_DMP)

# ===== Phase 2: VLP benchmarks =====
from pyrestoolbox.nodal.nodal import _hb_fbhp_gas, _wg_fbhp_gas, _gray_fbhp_gas, _bb_core_gas

N_VLP = 50
print(f"\n--- HB FBHP Gas ({N_VLP:,} calls) ---")
results["vlp_hb_gas"] = bench(f"{mode} HB gas", _hb_fbhp_gas,
    (300.0, 35.0, 0.65, 2.441, 0.0006, 9000.0, 120.0, 200.0, 0.0,
     5.0, 0.0, 0.0, 0.0, False, 0.0, 0.0), n=N_VLP)

print(f"\n--- WG FBHP Gas ({N_VLP:,} calls) ---")
results["vlp_wg_gas"] = bench(f"{mode} WG gas", _wg_fbhp_gas,
    (300.0, 35.0, 0.65, 2.441, 0.0006, 9000.0, 120.0, 200.0, 0.0,
     5.0, 0.0, 0.0, 0.0, False, 0.0, 0.0), n=N_VLP)

print(f"\n--- Gray FBHP Gas ({N_VLP:,} calls) ---")
results["vlp_gray_gas"] = bench(f"{mode} Gray gas", _gray_fbhp_gas,
    (300.0, 35.0, 0.65, 2.441, 0.0006, 9000.0, 120.0, 200.0, 0.0,
     5.0, 0.0, 0.0, 0.0, False, 0.0, 0.0), n=N_VLP)

print(f"\n--- BB FBHP Gas ({N_VLP:,} calls) ---")
results["vlp_bb_gas"] = bench(f"{mode} BB gas", _bb_core_gas,
    (300.0, 35.0, 0.65, 2.441, 0.0006, 9000.0, 120.0, 200.0, 0.0,
     5.0, 0.0, 0.0, 0.0, False, 0.0, 0.0), n=N_VLP)

# ===== Phase 2: DCA Hyperbolic benchmarks =====
from pyrestoolbox.dca.dca import _fit_hyperbolic, _fit_hyperbolic_cum

N_DCA = 200
qi_true, di_true, b_true = 1000.0, 0.05, 0.5
t_dca = np.arange(1, 61, dtype=float)
q_dca = qi_true / (1.0 + b_true * di_true * t_dca) ** (1.0 / b_true)

print(f"\n--- Hyperbolic fit time-domain ({N_DCA:,} calls) ---")
results["dca_hyp_time"] = bench(f"{mode} hyp fit", _fit_hyperbolic,
    (t_dca, q_dca), n=N_DCA)

# Cumulative production for cum-domain fit
np_cum = np.cumsum(q_dca)
print(f"\n--- Hyperbolic fit cum-domain ({N_DCA:,} calls) ---")
results["dca_hyp_cum"] = bench(f"{mode} hyp cum", _fit_hyperbolic_cum,
    (np_cum, q_dca), n=N_DCA)

# ===== Phase 2: Oil Density benchmarks =====
from pyrestoolbox import oil

N_OIL = 5_000
print(f"\n--- Oil density SWMH ({N_OIL:,} calls) ---")
results["oil_deno"] = bench(f"{mode} oil_deno", oil.oil_deno,
    (2000.0, 200.0, 500.0, 600.0, 0.8, 0.8, 3000.0, 0.0, 35.0), n=N_OIL)

print(f"\n--- Oil Bo McCain ({N_OIL:,} calls) ---")
results["oil_bo"] = bench(f"{mode} oil_bo", oil.oil_bo,
    (2000.0, 3000.0, 200.0, 500.0, 600.0, oil.oil_sg(35.0), 0.8, 0.8), n=N_OIL)

# ===== Phase 2: Spycher-Pruess CO2-Brine benchmarks =====
from pyrestoolbox import brine

N_SP = 500
print(f"\n--- Spycher-Pruess CO2-Brine ({N_SP:,} calls) ---")
results["sp_co2_brine"] = bench(f"{mode} SP CO2", lambda: brine.CO2_Brine_Mixture(
    pres=200.0, temp=135.0, ppm=30000.0, metric=True),
    (), n=N_SP)

# ===== Phase 2: VLE Flash benchmarks =====
from pyrestoolbox.brine._lib_vle_engine import SWMultiComponentFlash

N_VLE = 500
flash_obj = SWMultiComponentFlash(['H2O', 'CH4', 'CO2'],
    salinity_molal=0.529, framework='proposed')
z_vle = np.array([0.95, 0.04, 0.01])
T_vle = 373.15  # 100 degC
P_vle = 200e5   # 200 bar

print(f"\n--- VLE calc_equilibrium ({N_VLE:,} calls) ---")
results["vle_equilibrium"] = bench(f"{mode} VLE equil", flash_obj.calc_equilibrium,
    (T_vle, P_vle, z_vle), n=N_VLE)

print(f"RESULTS:{json.dumps(results)}")
