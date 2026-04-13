"""
Rust acceleration layer for pyResToolbox.

Attempts to load compiled Rust extensions at import time.
Falls back to pure Python silently on any failure.
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


def _get_cache_dir() -> Path:
    if sys.platform == "win32":
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Caches"
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
    return base / "pyrestoolbox"


def _get_extension_identity():
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


def _read_sentinel():
    try:
        sentinel_file = _sentinel_path()
        if not sentinel_file.exists():
            return None
        data = json.loads(sentinel_file.read_text(encoding="utf-8"))
        if not all(k in data for k in ("extension_path", "extension_mtime", "extension_size")):
            return None
        return data
    except Exception:
        return None


def _sentinel_matches_current(sentinel, current_identity):
    return (
        sentinel["extension_path"] == current_identity["extension_path"]
        and sentinel["extension_mtime"] == current_identity["extension_mtime"]
        and sentinel["extension_size"] == current_identity["extension_size"]
    )


def _write_sentinel(identity, reason):
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


def _clear_sentinel():
    try:
        sentinel_file = _sentinel_path()
        if sentinel_file.exists():
            sentinel_file.unlink()
            logger.debug("Cleared stale Rust-blocked sentinel")
    except Exception:
        pass


# --- Main probe logic ---
_force_python = os.environ.get("PYRESTOOLBOX_NO_RUST", "").strip() in ("1", "true", "yes")
_force_retry = os.environ.get("PYRESTOOLBOX_RETRY_RUST", "").strip() in ("1", "true", "yes")

if _force_python:
    _failure_reason = "disabled via PYRESTOOLBOX_NO_RUST environment variable"
    logger.debug("pyResToolbox Rust acceleration disabled by environment variable")
else:
    _ext_identity = _get_extension_identity()

    if _ext_identity is None:
        _failure_reason = "compiled extension not found"
        logger.debug("No Rust extension found on disk. Using pure Python.")
    else:
        _skip_probe = False

        if not _force_retry:
            _sentinel = _read_sentinel()
            if _sentinel is not None and _sentinel_matches_current(_sentinel, _ext_identity):
                _skip_probe = True
                _failure_reason = (
                    f"skipped probe (previously blocked: {_sentinel.get('failure_reason', 'unknown')}). "
                    f"Set PYRESTOOLBOX_RETRY_RUST=1 to retry."
                )
                logger.debug("Rust probe skipped - sentinel indicates previous block.")

        if _skip_probe:
            pass
        else:
            try:
                from pyrestoolbox import _native

                _native._smoke_test()

                _rust_module = _native
                RUST_AVAILABLE = True
                logger.debug("pyResToolbox Rust acceleration loaded successfully")
                _clear_sentinel()

            except ImportError as e:
                _failure_reason = f"ImportError: {e}"
                logger.debug("Rust extension import failed: %s", e)

            except (OSError, PermissionError) as e:
                _failure_reason = f"{type(e).__name__}: {e}"
                logger.debug("Rust extension blocked: %s", e)
                _write_sentinel(_ext_identity, _failure_reason)

            except Exception as e:
                _failure_reason = f"Unexpected: {type(e).__name__}: {e}"
                logger.debug("Rust extension failed: %s", e)
                _write_sentinel(_ext_identity, _failure_reason)


def get_status():
    status = {
        "rust_available": RUST_AVAILABLE,
        "failure_reason": _failure_reason if not RUST_AVAILABLE else "",
        "forced_python": _force_python,
        "sentinel_path": str(_sentinel_path()),
        "sentinel_exists": _sentinel_path().exists(),
    }
    if RUST_AVAILABLE:
        status["rust_version"] = get_rust_version()
    return status


def get_rust_version():
    """Return Rust extension version string, or None if unavailable."""
    if not RUST_AVAILABLE or _rust_module is None:
        return None
    return getattr(_rust_module, '__version__', getattr(_rust_module, 'version', None))


def clear_block():
    _clear_sentinel()
    return {
        "sentinel_cleared": True,
        "note": "Restart Python process to retry Rust extension loading",
        **get_status(),
    }


def rust_accelerated(rust_fn_name):
    """Decorator that dispatches to a Rust implementation when available.

    The decorated function is the pure-Python fallback. When Rust is available,
    the decorator calls ``_rust_module.<rust_fn_name>`` with the same positional
    and keyword arguments. On ImportError or AttributeError (missing function),
    it falls back to the Python implementation transparently.

    Usage::

        @rust_accelerated('hb_fbhp_gas_rust')
        def _hb_fbhp_gas(thp, api, gsg, ...):
            return _segment_march_gas(...)  # pure-Python path
    """
    import functools

    def decorator(fn):
        @functools.wraps(fn)
        def wrapper(*args, **kwargs):
            if RUST_AVAILABLE:
                try:
                    rust_fn = getattr(_rust_module, rust_fn_name)
                    return rust_fn(*args, **kwargs)
                except (ImportError, AttributeError):
                    pass
            return fn(*args, **kwargs)
        wrapper._rust_fn_name = rust_fn_name
        wrapper._python_fn = fn
        return wrapper
    return decorator
