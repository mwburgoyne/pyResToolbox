#!/usr/bin/env python3
"""Build a pure-Python wheel (no Rust extension) as a platform fallback.

Usage:
    python build_pure_python.py

Produces a universal wheel in dist/ that works on any platform.
Users get the pure-Python implementation; the Rust accelerator is not included.
Pre-built Rust-accelerated wheels (from CI) take priority when available.
"""

import subprocess
import sys
import tempfile
import shutil
from pathlib import Path

ROOT = Path(__file__).parent

def main():
    # Write a temporary pyproject.toml for setuptools (pure Python)
    # Strategy: keep the full [project] section from the real pyproject.toml
    # (so modern setuptools gets name, version, dependencies, etc.) but swap
    # [build-system] from maturin to setuptools.
    import re
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)

        real_content = (ROOT / "pyproject.toml").read_text()

        # Replace [build-system] section (everything up to next [section] or EOF)
        pp_content = re.sub(
            r'^\[build-system\].*?(?=^\[(?!build-system)|\Z)',
            '[build-system]\n'
            'requires = ["setuptools", "wheel"]\n'
            'build-backend = "setuptools.build_meta"\n\n',
            real_content,
            flags=re.MULTILINE | re.DOTALL,
        )
        # Remove [tool.maturin] section (not relevant for setuptools)
        pp_content = re.sub(
            r'^\[tool\.maturin\].*?(?=^\[(?!tool\.maturin)|\Z)',
            '',
            pp_content,
            flags=re.MULTILINE | re.DOTALL,
        )

        pp_toml = tmp / "pyproject.toml"
        pp_toml.write_text(pp_content)

        # Build using setup.py + the temporary pyproject.toml
        # We swap pyproject.toml temporarily
        real_toml = ROOT / "pyproject.toml"
        backup = ROOT / "pyproject.toml.maturin"

        shutil.copy2(real_toml, backup)
        shutil.copy2(pp_toml, real_toml)

        try:
            dist_dir = ROOT / "dist"
            dist_dir.mkdir(exist_ok=True)
            subprocess.check_call(
                [sys.executable, "-m", "build", "--wheel", "--outdir", str(dist_dir)],
                cwd=str(ROOT),
            )
        finally:
            # Restore maturin pyproject.toml
            shutil.copy2(backup, real_toml)
            backup.unlink()

    print("\nPure-Python wheel built in dist/")


if __name__ == "__main__":
    main()
