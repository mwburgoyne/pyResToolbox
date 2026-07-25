# pyResToolbox Rust Acceleration - Claude Code Instructions

## Objective

Add optional Rust-compiled acceleration modules to pyResToolbox for computationally expensive functions. The Rust layer must be completely invisible to users - if compiled extensions are available and loadable, they are used automatically. If not (missing wheel, unsupported platform, or corporate runtime protection like Airlock blocking non-whitelisted executables), the existing pure Python implementations are used with zero user intervention required.

## Architecture Principles

1. **Pure Python is the contract.** Every function must have a complete, correct, tested pure Python implementation. Rust is an optimisation, never a requirement.
2. **Single import-time probe.** Detect Rust availability once at import, not on every function call.
3. **Fail safe, fail quiet.** If the Rust module cannot load for any reason (ImportError, OSError, permission denied, segfault on load, Airlock intercept), log a single debug-level message and continue with Python. Never raise, never warn to stderr, never print.
4. **Numerical equivalence.** Rust and Python implementations must produce identical results within floating-point tolerance (relative error < 1e-10 for Z-factor and critical properties, < 1e-8 for integrated quantities like pseudo-pressure).
5. **No new user-facing dependencies.** A `pip install pyrestoolbox` from a wheel must work exactly as today. The Rust extension is bundled inside the wheel as a compiled `.so`/`.pyd` - users never see Cargo, rustc, or any Rust toolchain.

## Technology Stack

- **PyO3** (latest stable) for Rust-Python bindings
- **maturin** as the PEP 517 build backend (replaces setuptools/setup.py for the compiled extension only)
- **cibuildwheel** in GitHub Actions for cross-platform wheel builds
- **numpy** interop via PyO3's numpy bindings where array I/O is needed (VLP pressure traverses return arrays)

## Project Structure

```
pyrestoolbox/
├── Cargo.toml                          # Rust crate configuration
├── pyproject.toml                      # maturin build backend config
├── src/                                # Rust source
│   ├── lib.rs                          # PyO3 module definition, exports
│   ├── zfactor/
│   │   ├── mod.rs
│   │   ├── dak.rs                      # Dranchuk-Abou-Kassem (1975)
│   │   ├── hall_yarborough.rs          # Hall-Yarborough (1973)
│   │   └── bns.rs                      # Borges-Experiment, Experiment-Table,
│   │                                   #   standing-based five-component mixing
│   ├── critical_properties/
│   │   ├── mod.rs
│   │   ├── sutton.rs                   # Sutton (1985) pseudocritical correlations
│   │   ├── wichert_aziz.rs            # Wichert-Aziz (1972) acid gas correction
│   │   │                               #   applied to Sutton Ppc/Tpc for use
│   │   │                               #   with DAK and Hall-Yarborough
│   │   └── bns_mixing.rs              # BNS five-component (CO2/H2S/N2/H2/HC)
│   │                                   #   pseudocritical property mixing rules
│   │                                   #   (does NOT use Wichert-Aziz)
│   ├── gas_viscosity/
│   │   ├── mod.rs
│   │   ├── lee_gonzalez_eakin.rs       # Lee-Gonzalez-Eakin (1966)
│   │   └── carr.rs                     # Carr-Kobayashi-Burrows (1954)
│   ├── pseudopressure.rs               # Numerical integration of 2p/(mu*Z)
│   └── vlp/
│       ├── mod.rs
│       ├── beggs_brill.rs              # Beggs & Brill (1973) revised
│       ├── hagedorn_brown.rs           # Hagedorn & Brown (1965)
│       ├── gray.rs                     # Gray (1978) - gas wells
│       └── pressure_traverse.rs        # Marching algorithm, calls correlations
├── pyrestoolbox/                       # Python package (existing)
│   ├── __init__.py                     # Existing, modified to import _accel
│   ├── _accelerator.py                 # NEW - Rust probe and dispatch layer
│   ├── gas.py                          # Existing gas property functions
│   ├── zfactor.py                      # Existing Z-factor implementations
│   └── ...
└── tests/
    ├── test_rust_parity.py             # Numerical equivalence tests
    └── test_fallback.py                # Fallback behaviour tests
```

## The Accelerator Module: `_accelerator.py`

This is the critical piece. It must handle every failure mode gracefully, including corporate endpoint protection software (Airlock, CrowdStrike, Carbon Black, etc.) that intercepts and blocks execution of non-whitelisted shared libraries.

### Why a sentinel file?

The try/except probe runs once per Python process (module-level code, cached in `sys.modules`). The per-call dispatch is just `if RUST_AVAILABLE:` - a boolean lookup, essentially free. So the cost of re-probing is not performance.

The real problem is **security noise**. If Airlock blocks the `.so`/`.pyd` load, it typically logs a security event. Every new Python process that imports pyResToolbox triggers another alert. Restart a Jupyter kernel ten times while debugging, that is ten Airlock alerts. Corporate IT starts asking questions, or worse, escalates.

The solution: on first runtime failure (not ImportError - that just means no wheel, which is benign), write a small sentinel file. On subsequent imports, check the sentinel before attempting the dlopen. The sentinel is keyed to the specific `.so`/`.pyd` file's identity (path + mtime + size), so installing a new wheel automatically invalidates the cache and retries.

### Sentinel file location and lifecycle

The sentinel lives in the platform-appropriate user cache directory:
- Linux/macOS: `~/.cache/pyrestoolbox/rust_blocked.json`
- Windows: `%LOCALAPPDATA%\pyrestoolbox\rust_blocked.json`

The sentinel is a small JSON file containing:
```json
{
    "extension_path": "/path/to/_native.cpython-312-x86_64-linux-gnu.so",
    "extension_mtime": 1711900000.0,
    "extension_size": 524288,
    "failure_reason": "OSError (possible runtime protection): ...",
    "blocked_at": "2026-03-31T14:22:00Z"
}
```

**Invalidation triggers** (sentinel is ignored and probe retried):
1. Extension file path, mtime, or size has changed (new wheel installed)
2. `PYRESTOOLBOX_RETRY_RUST=1` environment variable is set
3. Sentinel file is manually deleted
4. Sentinel file is corrupt or unreadable (fail open - retry the probe)

**What does NOT write a sentinel:**
- `ImportError` (no compiled extension present at all - this is expected for sdist installs and doesn't generate security events, so no need to cache)
- `PYRESTOOLBOX_NO_RUST=1` (explicit env var override doesn't need persistence)

**What DOES write a sentinel:**
- `OSError`, `PermissionError`, or any unexpected exception during load or smoke test - these are the runtime-protection signatures

```python
"""
Rust acceleration layer for pyResToolbox.

Attempts to load compiled Rust extensions at import time.
Falls back to pure Python silently on any failure:
  - Platform without a pre-built wheel (ImportError)
  - Rust extension not installed / sdist install without Rust toolchain
  - Corporate runtime protection (Airlock, CrowdStrike, etc.) blocking
    execution of non-whitelisted .so/.pyd files (OSError, PermissionError)
  - Any other unexpected error during shared library load

On runtime-protection failures (OSError, PermissionError, etc.), a
sentinel file is written to the user cache directory. Subsequent imports
skip the probe entirely, avoiding repeated security alerts. The sentinel
auto-invalidates when a new wheel is installed (file identity changes).

The probe runs at most once per Python process. After that,
RUST_AVAILABLE is a module-level boolean and function references
point to either Rust or Python implementations.
"""

import json
import logging
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

logger = logging.getLogger(__name__)

RUST_AVAILABLE: bool = False
_rust_module = None
_failure_reason: str = ""


# ---------------------------------------------------------------------------
# Cache directory helpers
# ---------------------------------------------------------------------------

def _get_cache_dir() -> Path:
    """Return platform-appropriate user cache directory for pyResToolbox."""
    if sys.platform == "win32":
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Caches"
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
    return base / "pyrestoolbox"


def _get_extension_identity() -> dict | None:
    """Get identity (path, mtime, size) of the compiled extension file.

    Returns None if the extension file does not exist on disk (sdist
    install, no wheel). This is distinct from the extension existing
    but being blocked at runtime.
    """
    try:
        import importlib.util
        spec = importlib.util.find_spec("pyrestoolbox._native")
        if spec is None or spec.origin is None:
            return None
        ext_path = Path(spec.origin)
        if not ext_path.exists():
            return None
        stat = ext_path.stat()
        return {
            "extension_path": str(ext_path),
            "extension_mtime": stat.st_mtime,
            "extension_size": stat.st_size,
        }
    except Exception:
        return None


def _sentinel_path() -> Path:
    return _get_cache_dir() / "rust_blocked.json"


def _read_sentinel() -> dict | None:
    """Read and validate the sentinel file.

    Returns the sentinel dict if valid and current, None otherwise.
    A None return means "go ahead and probe".
    """
    try:
        sentinel_file = _sentinel_path()
        if not sentinel_file.exists():
            return None
        data = json.loads(sentinel_file.read_text(encoding="utf-8"))
        # Validate required fields
        if not all(k in data for k in ("extension_path", "extension_mtime", "extension_size")):
            return None
        return data
    except Exception:
        # Corrupt sentinel - ignore it, retry the probe
        return None


def _sentinel_matches_current(sentinel: dict, current_identity: dict) -> bool:
    """Check whether sentinel was written for the currently installed extension."""
    return (
        sentinel["extension_path"] == current_identity["extension_path"]
        and sentinel["extension_mtime"] == current_identity["extension_mtime"]
        and sentinel["extension_size"] == current_identity["extension_size"]
    )


def _write_sentinel(identity: dict, reason: str) -> None:
    """Write a sentinel file recording that Rust loading was blocked.

    Best-effort only. If we can't write the file (read-only filesystem,
    permissions, etc.), we just move on. The worst case is that the probe
    retries next time, which is the pre-sentinel behaviour anyway.
    """
    try:
        cache_dir = _get_cache_dir()
        cache_dir.mkdir(parents=True, exist_ok=True)
        sentinel_data = {
            **identity,
            "failure_reason": reason,
            "blocked_at": datetime.now(timezone.utc).isoformat(),
        }
        _sentinel_path().write_text(
            json.dumps(sentinel_data, indent=2),
            encoding="utf-8",
        )
        logger.debug("Wrote Rust-blocked sentinel to %s", _sentinel_path())
    except Exception as e:
        logger.debug("Could not write sentinel file: %s", e)


def _clear_sentinel() -> None:
    """Remove the sentinel file. Called on successful Rust load."""
    try:
        sentinel_file = _sentinel_path()
        if sentinel_file.exists():
            sentinel_file.unlink()
            logger.debug("Cleared stale Rust-blocked sentinel")
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Main probe logic
# ---------------------------------------------------------------------------

# Environment variable overrides
_force_python = os.environ.get("PYRESTOOLBOX_NO_RUST", "").strip() in ("1", "true", "yes")
_force_retry = os.environ.get("PYRESTOOLBOX_RETRY_RUST", "").strip() in ("1", "true", "yes")

if _force_python:
    _failure_reason = "disabled via PYRESTOOLBOX_NO_RUST environment variable"
    logger.debug("pyResToolbox Rust acceleration disabled by environment variable")
else:
    # Step 1: Check if the compiled extension file even exists
    _ext_identity = _get_extension_identity()

    if _ext_identity is None:
        # No compiled extension on disk. Pure sdist install or unsupported
        # platform. No sentinel needed - ImportError is benign.
        _failure_reason = "compiled extension not found (sdist install or unsupported platform)"
        logger.debug("No Rust extension found on disk. Using pure Python.")

    else:
        # Extension file exists. Check sentinel before attempting load.
        _skip_probe = False

        if not _force_retry:
            _sentinel = _read_sentinel()
            if _sentinel is not None and _sentinel_matches_current(_sentinel, _ext_identity):
                # Sentinel is valid and matches current extension - skip probe
                _skip_probe = True
                _failure_reason = (
                    f"skipped probe (previously blocked: {_sentinel.get('failure_reason', 'unknown')}). "
                    f"Set PYRESTOOLBOX_RETRY_RUST=1 to retry, or delete {_sentinel_path()}"
                )
                logger.debug(
                    "Rust probe skipped - sentinel indicates previous block. "
                    "Set PYRESTOOLBOX_RETRY_RUST=1 or delete %s to retry.",
                    _sentinel_path(),
                )

        if _skip_probe:
            pass  # RUST_AVAILABLE stays False, _failure_reason already set

        else:
            # Actually attempt to load and run the extension
            try:
                from pyrestoolbox import _native

                # Smoke test: call a trivial function to verify the shared
                # library actually executes. This catches runtime protection
                # that allows dlopen but intercepts code execution.
                _native._smoke_test()

                _rust_module = _native
                RUST_AVAILABLE = True
                logger.debug("pyResToolbox Rust acceleration loaded successfully")

                # Clear any stale sentinel from a previous blocked install
                _clear_sentinel()

            except ImportError as e:
                # Extension file exists on disk but can't be imported.
                # Unusual but possible (e.g. Python version mismatch in
                # the .so filename). Don't write sentinel - this isn't
                # a security event.
                _failure_reason = f"ImportError: {e}"
                logger.debug(
                    "Rust extension import failed (ImportError: %s). Using pure Python.", e
                )

            except (OSError, PermissionError) as e:
                # Runtime protection signature. Write sentinel to avoid
                # repeated security events.
                _failure_reason = f"{type(e).__name__} (possible runtime protection): {e}"
                logger.debug(
                    "Rust extension blocked at runtime (%s: %s). "
                    "Writing sentinel and using pure Python.", type(e).__name__, e
                )
                _write_sentinel(_ext_identity, _failure_reason)

            except Exception as e:
                # Unexpected failure. Write sentinel as a precaution -
                # if it's transient, the user can retry with env var.
                _failure_reason = f"Unexpected: {type(e).__name__}: {e}"
                logger.debug(
                    "Rust extension failed unexpectedly (%s: %s). "
                    "Writing sentinel and using pure Python.",
                    type(e).__name__, e
                )
                _write_sentinel(_ext_identity, _failure_reason)


def get_status() -> dict:
    """Return diagnostic information about Rust acceleration status.

    Useful for debugging and support. Users can call:
        import pyrestoolbox as rtb
        print(rtb.acceleration_status())

    Returns:
        dict with keys:
            rust_available (bool): True if Rust backend is active
            failure_reason (str): Explanation if not available
            forced_python (bool): True if PYRESTOOLBOX_NO_RUST is set
            sentinel_path (str): Location of the sentinel file
            sentinel_exists (bool): Whether a sentinel file is present
    """
    return {
        "rust_available": RUST_AVAILABLE,
        "failure_reason": _failure_reason if not RUST_AVAILABLE else "",
        "forced_python": _force_python,
        "sentinel_path": str(_sentinel_path()),
        "sentinel_exists": _sentinel_path().exists(),
    }


def clear_block() -> dict:
    """Clear the Rust-blocked sentinel and re-probe.

    Call this if you believe the blocking condition has been resolved
    (e.g. Airlock policy updated, new wheel installed) and you want
    to retry without restarting Python with PYRESTOOLBOX_RETRY_RUST=1.

    Note: This clears the sentinel but cannot reload the extension
    within a running process (Python caches module imports). You will
    need to restart the Python process for the retry to take effect.

    Returns:
        dict: Updated status after clearing the sentinel.
    """
    _clear_sentinel()
    return {
        "sentinel_cleared": True,
        "note": "Restart Python process to retry Rust extension loading",
        **get_status(),
    }
```
```

## How Functions Dispatch (Pattern for Every Accelerated Function)

Each Python module that has Rust-accelerated functions should follow this pattern. The dispatch happens at function level, not module level, so that partial Rust availability is possible if needed in future.

Example for `zfactor.py`:

```python
from pyrestoolbox._accelerator import RUST_AVAILABLE, _rust_module

def _dak_zfactor_python(Pr: float, Tr: float) -> float:
    """Pure Python DAK Z-factor. The reference implementation."""
    # ... existing implementation unchanged ...

def dak_zfactor(Pr: float, Tr: float) -> float:
    """Dranchuk-Abou-Kassem (1975) Z-factor correlation.

    Uses Rust acceleration if available, otherwise pure Python.
    """
    if RUST_AVAILABLE:
        return _rust_module.dak_zfactor(Pr, Tr)
    return _dak_zfactor_python(Pr, Tr)
```

Do NOT use a decorator or metaclass for this. Keep it explicit and readable. The `if RUST_AVAILABLE` branch prediction will be essentially free after the first call.

## Rust Implementation Requirements

### `lib.rs` - Module Entry Point

```rust
use pyo3::prelude::*;

mod zfactor;
mod critical_properties;
mod gas_viscosity;
mod pseudopressure;
mod vlp;

/// Smoke test function called during import-time probe.
/// Must execute actual compiled code (not just return).
/// Validates the shared library is genuinely executable,
/// catching runtime protection that blocks code execution
/// after allowing dlopen.
#[pyfunction]
fn _smoke_test() -> PyResult<()> {
    // Do a trivial but non-optimisable-away computation
    let _x: f64 = (2.0_f64).sqrt();
    Ok(())
}

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_smoke_test, m)?)?;

    // Z-factor functions (correlation only - caller supplies Pr, Tr)
    m.add_function(wrap_pyfunction!(zfactor::dak_zfactor, m)?)?;
    m.add_function(wrap_pyfunction!(zfactor::hall_yarborough_zfactor, m)?)?;
    m.add_function(wrap_pyfunction!(zfactor::bns_zfactor, m)?)?;

    // Critical properties and corrections
    m.add_function(wrap_pyfunction!(critical_properties::sutton_pseudocritical, m)?)?;
    m.add_function(wrap_pyfunction!(critical_properties::wichert_aziz_correction, m)?)?;
    m.add_function(wrap_pyfunction!(critical_properties::bns_pseudocritical, m)?)?;

    // Convenience: full pipeline functions that chain critical properties
    // through to Z-factor in a single Rust call (no Python round-trips)
    m.add_function(wrap_pyfunction!(zfactor::dak_zfactor_full, m)?)?;
    m.add_function(wrap_pyfunction!(zfactor::hall_yarborough_zfactor_full, m)?)?;
    m.add_function(wrap_pyfunction!(zfactor::bns_zfactor_full, m)?)?;

    // Gas viscosity
    m.add_function(wrap_pyfunction!(gas_viscosity::lee_gonzalez_eakin, m)?)?;

    // Pseudo-pressure
    m.add_function(wrap_pyfunction!(pseudopressure::pseudo_pressure, m)?)?;
    m.add_function(wrap_pyfunction!(pseudopressure::pseudo_pressure_array, m)?)?;

    // VLP
    m.add_function(wrap_pyfunction!(vlp::pressure_traverse_beggs_brill, m)?)?;
    m.add_function(wrap_pyfunction!(vlp::pressure_traverse_hagedorn_brown, m)?)?;
    m.add_function(wrap_pyfunction!(vlp::pressure_traverse_gray, m)?)?;

    Ok(())
}
```

### Function Signatures (Rust Side)

All Rust functions exposed via PyO3 must accept and return standard Python types (f64, Vec<f64>, tuples). No custom Rust structs in the public API.

There are two levels of functions for Z-factor:

1. **Correlation-only** functions (accept Pr, Tr directly). Useful when the caller has already computed reduced properties.
2. **Full pipeline** functions (accept p, T, sg, composition). These chain critical properties, corrections, and the Z-factor correlation in a single Rust call with no Python round-trips. This is where the real speedup lives for pseudo-pressure and VLP, because the entire chain stays in Rust.

```rust
// ===================================================================
// Critical properties and corrections
// ===================================================================

// Sutton (1985) pseudocritical properties from gas specific gravity.
// Used as the first step for DAK and Hall-Yarborough methods.
#[pyfunction]
fn sutton_pseudocritical(
    sg: f64,               // Gas specific gravity (air = 1.0)
) -> PyResult<(f64, f64)> {
    // Returns (Ppc_psia, Tpc_degR)
    // ...
}

// Wichert-Aziz (1972) acid gas correction to Sutton pseudocriticals.
// Applied AFTER Sutton, BEFORE computing Pr and Tr for DAK or HY.
// NOT used with BNS (BNS has its own mixing rules).
#[pyfunction]
fn wichert_aziz_correction(
    ppc: f64,              // Uncorrected Ppc from Sutton, psia
    tpc: f64,              // Uncorrected Tpc from Sutton, deg R
    co2_frac: f64,         // CO2 mole fraction
    h2s_frac: f64,         // H2S mole fraction
) -> PyResult<(f64, f64)> {
    // Returns (Ppc_corrected, Tpc_corrected)
    // ...
}

// BNS five-component pseudocritical properties.
// Uses its own mixing rules for CO2, H2S, N2, H2, and HC.
// Does NOT apply Wichert-Aziz - the acid gas effects are handled
// internally by the BNS mixing rule coefficients.
#[pyfunction]
fn bns_pseudocritical(
    sg: f64,               // Hydrocarbon specific gravity
    co2_frac: f64,         // CO2 mole fraction
    h2s_frac: f64,         // H2S mole fraction
    n2_frac: f64,          // N2 mole fraction
    h2_frac: f64,          // H2 mole fraction
) -> PyResult<(f64, f64)> {
    // Returns (Ppc_psia, Tpc_degR) via BNS mixing rules
    // ...
}

// ===================================================================
// Z-factor correlations (reduced-property input)
// ===================================================================

// DAK and HY accept reduced properties directly.
// Caller is responsible for computing Pr = p/Ppc, Tr = T/Tpc
// using whatever critical property method they prefer.
#[pyfunction]
fn dak_zfactor(Pr: f64, Tr: f64) -> PyResult<f64> { /* ... */ }

#[pyfunction]
fn hall_yarborough_zfactor(Pr: f64, Tr: f64) -> PyResult<f64> { /* ... */ }

// BNS Z-factor correlation. Also accepts reduced properties,
// but these MUST be computed via bns_pseudocritical, not Sutton.
#[pyfunction]
fn bns_zfactor(Pr: f64, Tr: f64) -> PyResult<f64> { /* ... */ }

// ===================================================================
// Full pipeline functions (pressure, temperature, composition in;
// Z-factor out). These keep the entire calculation chain in Rust.
// ===================================================================

// DAK full pipeline: Sutton -> Wichert-Aziz -> DAK
// This is the high-value function for pseudo-pressure and VLP,
// because calling it 500 times in Rust avoids 500 Python round-trips.
#[pyfunction]
fn dak_zfactor_full(
    p_psia: f64,           // Pressure, psia
    t_degf: f64,           // Temperature, deg F
    sg: f64,               // Gas specific gravity
    co2_frac: f64,         // CO2 mole fraction
    h2s_frac: f64,         // H2S mole fraction
    n2_frac: f64,          // N2 mole fraction (applied via Kay's rule
                           //   adjustment to Sutton, or ignored if zero)
) -> PyResult<f64> {
    // Internally:
    //   1. Sutton pseudocriticals from sg
    //   2. Wichert-Aziz correction for CO2 + H2S
    //   3. Compute Pr = p / Ppc_corrected, Tr = T_rankine / Tpc_corrected
    //   4. DAK correlation
    // Returns Z
    // ...
}

// Hall-Yarborough full pipeline: Sutton -> Wichert-Aziz -> HY
#[pyfunction]
fn hall_yarborough_zfactor_full(
    p_psia: f64,
    t_degf: f64,
    sg: f64,
    co2_frac: f64,
    h2s_frac: f64,
    n2_frac: f64,
) -> PyResult<f64> {
    // Same chain as dak_zfactor_full but uses Hall-Yarborough at step 4
    // ...
}

// Note: BNS does NOT need a separate "full" wrapper because its
// bns_pseudocritical + bns_zfactor chain is already self-contained.
// But for API consistency, expose one anyway:
#[pyfunction]
fn bns_zfactor_full(
    p_psia: f64,
    t_degf: f64,
    sg: f64,
    co2_frac: f64,
    h2s_frac: f64,
    n2_frac: f64,
    h2_frac: f64,          // H2 only relevant for BNS
) -> PyResult<f64> {
    // Internally:
    //   1. BNS pseudocritical mixing (NOT Sutton, NOT Wichert-Aziz)
    //   2. Compute Pr, Tr from BNS pseudocriticals
    //   3. BNS Z-factor correlation
    // Returns Z
    // ...
}

// ===================================================================
// Pseudo-pressure
// ===================================================================

// Pseudo-pressure returning array. The zmethod parameter selects
// the full pipeline function internally, so the entire integration
// loop (hundreds of pressure steps, each calling Z-factor and
// viscosity) stays in Rust.
#[pyfunction]
fn pseudo_pressure_array(
    p_array: Vec<f64>,     // Pressure array, psia
    t_degf: f64,           // Temperature, deg F
    sg: f64,               // Gas specific gravity
    co2_frac: f64,
    h2s_frac: f64,
    n2_frac: f64,
    h2_frac: f64,          // Ignored unless zmethod = "BNS"
    p_base: f64,           // Base pressure for integration, psia
    zmethod: &str,         // "DAK", "HY", or "BNS"
) -> PyResult<Vec<f64>> {
    // Returns pseudo-pressure at each pressure in p_array
    // Internally calls dak_zfactor_full, hy_zfactor_full, or
    // bns_zfactor_full depending on zmethod - no Python involved.
    // ...
}

// ===================================================================
// VLP pressure traverse
// ===================================================================

#[pyfunction]
fn pressure_traverse_beggs_brill(
    p_start: f64,          // Starting pressure, psia
    t_start: f64,          // Starting temperature, deg F
    depth_array: Vec<f64>, // Measured depth array, ft
    inclination: Vec<f64>, // Well inclination at each depth, degrees
    d_tubing: f64,         // Tubing ID, inches
    roughness: f64,        // Pipe roughness, inches
    q_oil: f64,            // Oil rate, STB/d
    q_gas: f64,            // Gas rate, Mscf/d
    q_water: f64,          // Water rate, STB/d
    api: f64,              // Oil API gravity
    sg_gas: f64,           // Gas specific gravity
    sg_water: f64,         // Water specific gravity
    direction: &str,       // "up" or "down"
) -> PyResult<Vec<f64>> {
    // Returns pressure at each depth point
    // ...
}
```

### Critical Property Chain Summary

The relationship between critical property methods and Z-factor correlations:

```
DAK / Hall-Yarborough path:
  sg ──> Sutton (1985) ──> Ppc, Tpc
                              │
  co2_frac, h2s_frac ──> Wichert-Aziz (1972) ──> Ppc', Tpc'
                                                      │
  p, T ──────────────────────────────────────> Pr, Tr ──> DAK or HY ──> Z

BNS path:
  sg, co2, h2s, n2, h2 ──> BNS mixing rules ──> Ppc, Tpc
                                                      │
  p, T ──────────────────────────────────────> Pr, Tr ──> BNS ──> Z

  (No Wichert-Aziz. BNS mixing rules handle acid gas
   effects through their own fitted coefficients.)
```
```

### Numerical Precision Rules

- All internal arithmetic in f64 (never f32)
- Newton-Raphson convergence tolerance: 1e-12 (tighter than Python, costs nothing in Rust)
- Maximum iterations for Z-factor solvers: 100 (with convergence check)
- Pseudo-pressure integration: use adaptive Simpson's rule or Gauss-Legendre quadrature, not fixed-step trapezoidal
- Sutton + Wichert-Aziz chain (DAK/HY path): maintain full precision through the correction; do not round intermediate Ppc/Tpc values
- BNS mixing rules: maintain full precision through the five-component mixing; the fitted coefficients are sensitive to truncation

### Error Handling (Rust Side)

Return `PyResult<T>` from all exported functions. Map Rust errors to Python exceptions:

```rust
use pyo3::exceptions::PyValueError;

if pr <= 0.0 {
    return Err(PyValueError::new_err(
        format!("Reduced pressure must be positive, got {}", pr)
    ));
}

if iterations >= max_iter {
    return Err(PyValueError::new_err(
        format!("Z-factor failed to converge after {} iterations (Pr={}, Tr={})",
                max_iter, pr, tr)
    ));
}
```

These must match the exceptions the Python implementations raise, so the user experience is identical regardless of backend.

## BNS Z-Factor Specific Requirements

The BNS correlation handles four non-hydrocarbon components: CO2, H2S, N2, and H2. The implementation must include:

1. **BNS pseudocritical property mixing rules** for the five-component gas (four non-HC + hydrocarbon fraction). These mixing rules have their own fitted coefficients that account for acid gas effects internally.
2. **The BNS Z-factor correlation itself** using the BNS-mixed pseudocritical properties
3. All functions must accept individual mole fractions (not a dictionary or struct) for clean PyO3 binding

**BNS does NOT use Wichert-Aziz.** The Wichert-Aziz correction is only applied in the Sutton pseudocritical path (used by DAK and Hall-Yarborough). The BNS mixing rules handle CO2 and H2S effects through their own regression coefficients. Do not apply Wichert-Aziz anywhere in the BNS chain.

Port the existing Python BNS implementation exactly. Do not simplify, approximate, or "improve" the correlation. Numerical equivalence with the Python version is mandatory.

## DAK / Hall-Yarborough Critical Property Chain

Unlike BNS, the DAK and Hall-Yarborough correlations are pure Z(Pr, Tr) functions that require external critical property computation. The full pipeline for these methods is:

1. **Sutton (1985)** pseudocritical properties from gas specific gravity
2. **Wichert-Aziz (1972)** acid gas correction applied to Sutton Ppc and Tpc, using CO2 and H2S mole fractions
3. Compute reduced properties: Pr = p / Ppc_corrected, Tr = T_rankine / Tpc_corrected
4. Apply the DAK or Hall-Yarborough Z-factor correlation

The `_full` pipeline functions (`dak_zfactor_full`, `hall_yarborough_zfactor_full`) perform all four steps in a single Rust call. This is critical for performance in pseudo-pressure integration where the function is called hundreds of times per evaluation.

## Build Configuration

### `Cargo.toml`

```toml
[package]
name = "pyrestoolbox-native"
version = "0.1.0"
edition = "2021"

[lib]
name = "_native"
crate-type = ["cdylib"]

[dependencies]
pyo3 = { version = "0.23", features = ["extension-module"] }

[profile.release]
opt-level = 3
lto = "fat"
codegen-units = 1
strip = true
```

### `pyproject.toml` Changes

The existing pyproject.toml needs modification to use maturin as the build backend while preserving the pure Python fallback for sdist installs.

```toml
[build-system]
requires = ["maturin>=1.7,<2.0"]
build-backend = "maturin"

[tool.maturin]
# Build the Rust extension as a submodule inside the Python package
module-name = "pyrestoolbox._native"
# Include the Python source in wheels and sdists
python-source = "."
# Features for the Rust build
features = ["pyo3/extension-module"]

# When building from sdist without Rust toolchain, maturin will fail.
# This is expected - the CI builds wheels for all platforms.
# Users installing from sdist get pure Python via the fallback.
```

**Important:** The sdist must still install successfully without Rust. Configure the build so that if `maturin` or `rustc` is unavailable, the package installs with just the Python files. This requires a fallback `setup.py` or alternative build config. The recommended approach:

- Primary path: maturin builds wheel with compiled extension
- Fallback path: provide a `setup.py` that installs only the Python package without the `_native` module. Use a try/except in `setup.py` or a conditional build backend.

One clean way to handle this:

```python
# setup.py (fallback for sdist installs without Rust)
from setuptools import setup, find_packages

setup(
    name="pyrestoolbox",
    packages=find_packages(),
    # ... existing setup kwargs ...
)
```

And in `pyproject.toml`, use `maturin` for wheel builds but include `setup.py` in the sdist for fallback installs. The CI workflow always uses maturin. End users installing from sdist (rare) get setuptools and pure Python.

## GitHub Actions CI Workflow

### `.github/workflows/build-wheels.yml`

```yaml
name: Build and Publish Wheels

on:
  push:
    tags:
      - 'v*'
  workflow_dispatch:

permissions:
  contents: read

jobs:
  build-wheels:
    name: Build wheels - ${{ matrix.os }}
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
        os: [ubuntu-latest, macos-13, macos-14, windows-latest]
        # macos-13 = x86_64, macos-14 = arm64 (Apple Silicon)

    steps:
      - uses: actions/checkout@v4

      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: '3.12'

      - name: Install Rust
        uses: dtolnay/rust-toolchain@stable

      - name: Build wheels
        uses: PyO3/maturin-action@v1
        with:
          command: build
          args: --release --out dist
          # manylinux for broad Linux compatibility
          manylinux: auto
          # Build for multiple Python versions
          target: ${{ matrix.os == 'macos-14' && 'aarch64-apple-darwin' || '' }}

      - name: Upload wheels
        uses: actions/upload-artifact@v4
        with:
          name: wheels-${{ matrix.os }}
          path: dist/*.whl

  build-wheels-linux-multi:
    # Additional Linux targets if needed (aarch64, etc.)
    name: Build wheels - Linux ${{ matrix.target }}
    runs-on: ubuntu-latest
    strategy:
      matrix:
        target: [x86_64, aarch64]
    steps:
      - uses: actions/checkout@v4
      - name: Build wheels
        uses: PyO3/maturin-action@v1
        with:
          target: ${{ matrix.target }}
          command: build
          args: --release --out dist
          manylinux: auto
      - uses: actions/upload-artifact@v4
        with:
          name: wheels-linux-${{ matrix.target }}
          path: dist/*.whl

  build-sdist:
    name: Build sdist
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - name: Build sdist
        uses: PyO3/maturin-action@v1
        with:
          command: sdist
          args: --out dist
      - uses: actions/upload-artifact@v4
        with:
          name: sdist
          path: dist/*.tar.gz

  # Target Python versions: 3.9, 3.10, 3.11, 3.12, 3.13
  # maturin-action handles the matrix automatically via
  # the Cargo.toml and pyproject.toml configuration

  publish:
    name: Publish to PyPI
    runs-on: ubuntu-latest
    needs: [build-wheels, build-wheels-linux-multi, build-sdist]
    # Only publish on tagged releases
    if: startsWith(github.ref, 'refs/tags/v')
    environment: pypi
    permissions:
      id-token: write  # For trusted publishing
    steps:
      - uses: actions/download-artifact@v4
        with:
          path: dist
          merge-multiple: true
      - name: Publish to PyPI
        uses: pypa/gh-action-pypi-publish@release/v1
        with:
          packages-dir: dist/
```

This produces approximately 15-20 wheels covering:
- Linux x86_64 (manylinux)
- Linux aarch64 (manylinux)
- macOS x86_64
- macOS arm64 (Apple Silicon)
- Windows x86_64
- Each across Python 3.9 through 3.13

## Testing Requirements

### `tests/test_rust_parity.py`

Test numerical equivalence between Rust and Python for every accelerated function. These tests must pass in CI on all platforms.

```python
"""
Numerical parity tests: Rust vs Python implementations.

These tests import both implementations explicitly and compare
results across a range of inputs, including edge cases.
"""
import pytest
import math
from pyrestoolbox._accelerator import RUST_AVAILABLE, _rust_module
from pyrestoolbox.zfactor import (
    _dak_zfactor_python,
    _hall_yarborough_zfactor_python,
    _bns_zfactor_python,
)

# Skip entire module if Rust is not available in this environment
pytestmark = pytest.mark.skipif(
    not RUST_AVAILABLE,
    reason="Rust extension not available"
)

# Test vectors: (Pr, Tr, expected_Z_approximate)
# Cover low, medium, high reduced pressure and temperature
DAK_TEST_VECTORS = [
    (0.5, 1.05),
    (1.0, 1.2),
    (2.0, 1.5),
    (5.0, 2.0),
    (8.0, 1.4),
    (15.0, 2.5),
    # Near-critical (challenging convergence)
    (1.0, 1.01),
    (0.1, 1.05),
]

@pytest.mark.parametrize("Pr, Tr", DAK_TEST_VECTORS)
def test_dak_parity(Pr, Tr):
    py_result = _dak_zfactor_python(Pr, Tr)
    rust_result = _rust_module.dak_zfactor(Pr, Tr)
    assert math.isclose(py_result, rust_result, rel_tol=1e-10), (
        f"DAK parity failure at Pr={Pr}, Tr={Tr}: "
        f"Python={py_result}, Rust={rust_result}"
    )

# Similar parametrised tests for:
# - Hall-Yarborough
# - BNS (with various compositions including high CO2, high H2S, H2-bearing)
# - Pseudo-pressure (compare arrays element-wise)
# - VLP pressure traverse (compare full pressure profile)
# - Critical properties (Sutton, BNS mixing)
# - Gas viscosity

# Edge cases to include:
# - Pure methane (all non-HC fractions = 0)
# - High CO2 (30%+ CO2)
# - High H2S (30%+ H2S)
# - H2-bearing gas (5-20% H2, relevant for energy transition)
# - Near-critical conditions (Pr ~ 1.0, Tr ~ 1.0-1.05)
# - Very high pressure (Pr > 15)
# - Low pressure (Pr < 0.2)
```

### `tests/test_fallback.py`

Test that the fallback and sentinel mechanisms work correctly. These tests use subprocesses to get clean module state (Python caches modules in `sys.modules`, so `importlib.reload` is unreliable for testing import-time logic).

```python
"""
Fallback and sentinel behaviour tests.

Verify that pyResToolbox:
  - Works correctly when Rust acceleration is unavailable
  - Environment variable overrides work
  - Sentinel file is written on runtime-protection failures
  - Sentinel file prevents repeated probe attempts
  - Sentinel auto-invalidates when a new wheel is installed
  - PYRESTOOLBOX_RETRY_RUST=1 overrides the sentinel
  - clear_block() removes the sentinel
"""
import json
import os
import subprocess
import sys
import tempfile

import pytest


def _run_in_subprocess(code: str, env_overrides: dict = None) -> subprocess.CompletedProcess:
    """Run Python code in a clean subprocess with optional env var overrides."""
    env = os.environ.copy()
    if env_overrides:
        env.update(env_overrides)
    return subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True, text=True, env=env,
    )


def test_env_var_disables_rust():
    """PYRESTOOLBOX_NO_RUST=1 forces pure Python."""
    result = _run_in_subprocess("""
import pyrestoolbox._accelerator as acc
assert not acc.RUST_AVAILABLE
status = acc.get_status()
assert status['forced_python'] is True
assert 'environment variable' in status['failure_reason']
print('PASS')
""", env_overrides={"PYRESTOOLBOX_NO_RUST": "1"})
    assert "PASS" in result.stdout, result.stderr


def test_functions_work_without_rust():
    """All public functions produce valid results without Rust."""
    result = _run_in_subprocess("""
import pyrestoolbox as rtb

# Test Z-factor
z = rtb.zfactor.dak_zfactor(2.0, 1.5)
assert 0.0 < z < 2.0, f"Z-factor out of range: {z}"

# Test BNS Z-factor
z_bns = rtb.zfactor.bns_zfactor(
    p_psia=2000, t_degf=200, sg=0.7,
    co2_frac=0.05, h2s_frac=0.02,
    n2_frac=0.01, h2_frac=0.0
)
assert 0.0 < z_bns < 2.0, f"BNS Z-factor out of range: {z_bns}"

print('PASS')
""", env_overrides={"PYRESTOOLBOX_NO_RUST": "1"})
    assert "PASS" in result.stdout, result.stderr


class TestSentinel:
    """Tests for the sentinel file mechanism.

    These tests use a temporary directory as the cache location
    to avoid polluting the real user cache.
    """

    def test_sentinel_written_on_oserror(self, tmp_path):
        """Sentinel file is created when the extension raises OSError."""
        sentinel_dir = tmp_path / "pyrestoolbox"
        result = _run_in_subprocess(f"""
import sys
import json
from pathlib import Path
from unittest.mock import patch, MagicMock

# Patch the cache dir to our temp location
import pyrestoolbox._accelerator as acc

# Simulate: extension file exists but raises OSError on import
fake_identity = {{
    "extension_path": "/fake/path/_native.so",
    "extension_mtime": 1234567890.0,
    "extension_size": 999,
}}

with patch.object(acc, '_get_cache_dir', return_value=Path("{sentinel_dir}")), \\
     patch.object(acc, '_get_extension_identity', return_value=fake_identity):
    acc._write_sentinel(fake_identity, "OSError: blocked by Airlock")

sentinel_file = Path("{sentinel_dir}") / "rust_blocked.json"
assert sentinel_file.exists(), "Sentinel file not created"

data = json.loads(sentinel_file.read_text())
assert data["extension_path"] == "/fake/path/_native.so"
assert "Airlock" in data["failure_reason"]
assert "blocked_at" in data

print('PASS')
""")
        assert "PASS" in result.stdout, result.stderr

    def test_sentinel_skips_probe(self, tmp_path):
        """When a valid sentinel exists, the probe is skipped entirely."""
        # Pre-create a sentinel file
        sentinel_dir = tmp_path / "pyrestoolbox"
        sentinel_dir.mkdir(parents=True)
        sentinel_file = sentinel_dir / "rust_blocked.json"
        sentinel_data = {
            "extension_path": "/fake/_native.so",
            "extension_mtime": 1234567890.0,
            "extension_size": 999,
            "failure_reason": "OSError: test block",
            "blocked_at": "2026-03-31T00:00:00+00:00",
        }
        sentinel_file.write_text(json.dumps(sentinel_data))

        result = _run_in_subprocess(f"""
from pathlib import Path
from unittest.mock import patch

import pyrestoolbox._accelerator as acc

fake_identity = {{
    "extension_path": "/fake/_native.so",
    "extension_mtime": 1234567890.0,
    "extension_size": 999,
}}

# Reload with patched paths and matching identity
with patch.object(acc, '_get_cache_dir', return_value=Path("{sentinel_dir}")), \\
     patch.object(acc, '_get_extension_identity', return_value=fake_identity):
    sentinel = acc._read_sentinel()
    assert sentinel is not None, "Sentinel should be readable"
    assert acc._sentinel_matches_current(sentinel, fake_identity), "Sentinel should match"

print('PASS')
""")
        assert "PASS" in result.stdout, result.stderr

    def test_sentinel_invalidated_by_new_wheel(self, tmp_path):
        """Sentinel is ignored when extension file identity changes (new install)."""
        sentinel_dir = tmp_path / "pyrestoolbox"
        sentinel_dir.mkdir(parents=True)
        sentinel_file = sentinel_dir / "rust_blocked.json"
        sentinel_data = {
            "extension_path": "/fake/_native.so",
            "extension_mtime": 1234567890.0,  # Old mtime
            "extension_size": 999,
            "failure_reason": "OSError: old block",
            "blocked_at": "2026-03-01T00:00:00+00:00",
        }
        sentinel_file.write_text(json.dumps(sentinel_data))

        result = _run_in_subprocess(f"""
from pathlib import Path
from unittest.mock import patch

import pyrestoolbox._accelerator as acc

# New wheel has different mtime and size
new_identity = {{
    "extension_path": "/fake/_native.so",
    "extension_mtime": 9999999999.0,
    "extension_size": 1234,
}}

with patch.object(acc, '_get_cache_dir', return_value=Path("{sentinel_dir}")):
    sentinel = acc._read_sentinel()
    assert sentinel is not None, "Sentinel file should exist"
    assert not acc._sentinel_matches_current(sentinel, new_identity), \\
        "Sentinel should NOT match new extension identity"

print('PASS')
""")
        assert "PASS" in result.stdout, result.stderr

    def test_retry_env_var_overrides_sentinel(self, tmp_path):
        """PYRESTOOLBOX_RETRY_RUST=1 forces a probe even with valid sentinel."""
        # This is tested via the _force_retry flag in the main logic.
        # A full integration test would require a real compiled extension.
        result = _run_in_subprocess("""
import os
os.environ['PYRESTOOLBOX_RETRY_RUST'] = '1'
# Just verify the flag is read correctly
assert os.environ.get('PYRESTOOLBOX_RETRY_RUST', '').strip() in ('1', 'true', 'yes')
print('PASS')
""")
        assert "PASS" in result.stdout, result.stderr

    def test_clear_block_removes_sentinel(self, tmp_path):
        """clear_block() removes the sentinel file."""
        sentinel_dir = tmp_path / "pyrestoolbox"
        sentinel_dir.mkdir(parents=True)
        sentinel_file = sentinel_dir / "rust_blocked.json"
        sentinel_file.write_text('{"test": true}')

        result = _run_in_subprocess(f"""
from pathlib import Path
from unittest.mock import patch

import pyrestoolbox._accelerator as acc

sentinel_file = Path("{sentinel_file}")
assert sentinel_file.exists(), "Sentinel should exist before clear"

with patch.object(acc, '_sentinel_path', return_value=sentinel_file):
    result = acc.clear_block()
    assert result['sentinel_cleared'] is True

assert not sentinel_file.exists(), "Sentinel should be deleted after clear_block()"
print('PASS')
""")
        assert "PASS" in result.stdout, result.stderr

    def test_corrupt_sentinel_triggers_retry(self, tmp_path):
        """Corrupt or unparseable sentinel file causes a fresh probe."""
        sentinel_dir = tmp_path / "pyrestoolbox"
        sentinel_dir.mkdir(parents=True)
        sentinel_file = sentinel_dir / "rust_blocked.json"
        sentinel_file.write_text("this is not json {{{")

        result = _run_in_subprocess(f"""
from pathlib import Path
from unittest.mock import patch

import pyrestoolbox._accelerator as acc

with patch.object(acc, '_get_cache_dir', return_value=Path("{sentinel_dir}")):
    sentinel = acc._read_sentinel()
    assert sentinel is None, "Corrupt sentinel should return None (triggering retry)"

print('PASS')
""")
        assert "PASS" in result.stdout, result.stderr
```

## Implementation Order

Work through these in sequence. Each step should be committed and tested before moving to the next.

### Phase 1: Scaffolding
1. Create `Cargo.toml`, update `pyproject.toml` for maturin
2. Create `src/lib.rs` with just the `_smoke_test` function
3. Create `_accelerator.py` with the full probe and sentinel logic
4. Verify: `maturin develop` builds and imports successfully
5. Verify: `PYRESTOOLBOX_NO_RUST=1` disables it
6. Verify: deleting the `.so`/`.pyd` triggers graceful fallback
7. Verify sentinel lifecycle:
   - Simulate an OSError and confirm sentinel file is written
   - Confirm subsequent import skips the probe (check debug log output)
   - Confirm `PYRESTOOLBOX_RETRY_RUST=1` overrides the sentinel
   - Confirm sentinel auto-invalidates after `maturin develop` rebuilds (new mtime/size)
   - Confirm `clear_block()` removes the sentinel file
   - Confirm corrupt sentinel file triggers a fresh probe

### Phase 2: Critical Properties and Z-Factor Correlations
1. Port Sutton pseudocritical correlations to Rust
2. Port Wichert-Aziz acid gas correction to Rust
3. Port DAK Z-factor (Pr, Tr version) to Rust, expose via PyO3
4. Build `dak_zfactor_full` pipeline (Sutton -> Wichert-Aziz -> DAK), expose via PyO3
5. Wire into `zfactor.py` dispatch
6. Parity tests (both reduced-property and full-pipeline versions)
7. Port Hall-Yarborough (Pr, Tr version and `_full` pipeline), same process
8. Port BNS pseudocritical mixing rules to Rust (no Wichert-Aziz)
9. Port BNS Z-factor correlation and `bns_zfactor_full` pipeline
10. Parity tests for BNS (including high-CO2, high-H2S, and H2-bearing compositions)
11. Benchmark: time 100,000 Z-factor evaluations per method, Python vs Rust

### Phase 3: Gas Viscosity and Pseudo-Pressure
1. Port Lee-Gonzalez-Eakin viscosity to Rust
2. Implement pseudo-pressure integration in Rust (calling Rust Z-factor and viscosity internally, no Python round-trips)
3. Expose both single-point and array versions
4. Parity tests with tight tolerances on the integrated quantity
5. Benchmark: pseudo-pressure table generation (500 pressure points)

### Phase 4: VLP Correlations
1. Port Beggs & Brill revised correlation to Rust
2. Port pressure traverse marching algorithm
3. Port Hagedorn & Brown
4. Port Gray (gas wells)
5. Parity tests across representative well configurations
6. Benchmark: single well VLP curve (200 depth increments, 20 rate steps)

### Phase 5: CI and Packaging
1. Set up GitHub Actions workflow (use template above)
2. Build wheels for all target platforms
3. Test wheel installation on clean VMs/containers for each platform
4. Configure PyPI trusted publishing
5. Update README with acceleration status badge / documentation
6. Tag and release

## Public API Additions

Add to the top-level `pyrestoolbox` namespace:

```python
def acceleration_status() -> dict:
    """Check whether Rust acceleration is active.

    Returns a dict with:
        rust_available (bool): True if Rust backend is loaded
        failure_reason (str): Empty if available, otherwise explains why not
        forced_python (bool): True if PYRESTOOLBOX_NO_RUST is set
        sentinel_path (str): Location of the sentinel cache file
        sentinel_exists (bool): Whether a block sentinel is present

    Example:
        >>> import pyrestoolbox as rtb
        >>> rtb.acceleration_status()
        {'rust_available': True, 'failure_reason': '', 'forced_python': False,
         'sentinel_path': '/home/user/.cache/pyrestoolbox/rust_blocked.json',
         'sentinel_exists': False}
    """
    from pyrestoolbox._accelerator import get_status
    return get_status()


def clear_rust_block() -> dict:
    """Clear the Rust-blocked sentinel file.

    Call this if you believe the blocking condition has been resolved
    (e.g. security policy updated, new wheel installed) and want to
    retry Rust acceleration.

    Note: You must restart the Python process after calling this for
    the retry to take effect (Python caches module imports).

    Example:
        >>> import pyrestoolbox as rtb
        >>> rtb.clear_rust_block()
        {'sentinel_cleared': True, 'note': 'Restart Python process to retry...', ...}
    """
    from pyrestoolbox._accelerator import clear_block
    return clear_block()
```

This gives users (and support/debugging) three ways to manage the Rust backend:
1. `rtb.acceleration_status()` - check what is happening and why
2. `rtb.clear_rust_block()` - clear a stale sentinel from within Python
3. `PYRESTOOLBOX_RETRY_RUST=1` - retry on next import (e.g. from shell)
4. `PYRESTOOLBOX_NO_RUST=1` - force pure Python permanently (e.g. corporate env var policy)

## Things NOT To Do

- Do not add `rust` or `cargo` to the pip install requirements
- Do not require users to have a Rust toolchain installed
- Do not print warnings to stdout or stderr about Rust availability
- Do not add a `[rust]` extras_require - the extension is always bundled in wheels, always absent in sdist-without-Rust
- Do not use `ctypes` or `cffi` as an alternative to PyO3 - the ergonomics are worse and the maintenance burden is higher
- Do not attempt WASM compilation as an alternative - the performance overhead defeats the purpose
- Do not change any existing function signatures in the public API
- Do not add numpy as a hard dependency just for the Rust layer (use it only where VLP functions already return arrays)
- Do not use `warnings.warn()` for Rust unavailability - use `logging.debug()` only
- Do not attempt to detect specific security products by name (Airlock, etc.) - just catch the exception classes they raise
- Do not write the sentinel file on `ImportError` - that just means no wheel is installed, which is benign and does not generate security events
- Do not store the sentinel inside the package directory (it may be read-only, and `pip install --upgrade` would delete it anyway). Use the user cache directory.
- Do not add `platformdirs` as a dependency just for the cache path - use the manual platform detection shown in `_get_cache_dir()`
- Do not make the sentinel write a hard requirement - if the cache directory is not writable, skip it silently. The worst case is the probe runs again next time, which is the pre-sentinel behaviour

## Benchmarking Script

Include `benchmarks/bench_acceleration.py` for validating performance gains:

```python
"""
Benchmark Rust vs Python implementations.
Run with: python benchmarks/bench_acceleration.py
"""
import time
import os

def bench(name, func, args, n=100_000):
    start = time.perf_counter()
    for _ in range(n):
        func(*args)
    elapsed = time.perf_counter() - start
    rate = n / elapsed
    print(f"  {name}: {elapsed:.3f}s ({rate:,.0f} calls/sec)")
    return elapsed

print("=== Z-Factor (DAK) ===")
from pyrestoolbox.zfactor import _dak_zfactor_python, dak_zfactor
from pyrestoolbox._accelerator import RUST_AVAILABLE

py_time = bench("Python", _dak_zfactor_python, (2.0, 1.5))
if RUST_AVAILABLE:
    rs_time = bench("Rust  ", dak_zfactor, (2.0, 1.5))
    print(f"  Speedup: {py_time/rs_time:.1f}x")
else:
    print("  Rust not available - skipping")

# Repeat for pseudo-pressure, VLP, etc.
```

## Summary of Platform Coverage

| Platform | Architecture | How Covered |
|----------|-------------|-------------|
| Linux | x86_64 | manylinux wheel via CI |
| Linux | aarch64 | manylinux wheel via CI (cross-compile) |
| macOS | x86_64 | Native wheel via CI (macos-13 runner) |
| macOS | arm64 | Native wheel via CI (macos-14 runner) |
| Windows | x86_64 | Native wheel via CI |
| Other | Any | sdist install, pure Python fallback |
| Corporate locked-down | Any | Wheel installs but .so/.pyd blocked at runtime, pure Python fallback |

Users on any platform get `pip install pyrestoolbox` working first time. The compiled extension is a silent bonus on supported platforms with permissive runtime environments.
