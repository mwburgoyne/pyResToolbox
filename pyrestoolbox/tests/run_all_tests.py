#!/usr/bin/env python3
"""
Master test runner for pyResToolbox validation suite.

Thin delegation to pytest so that class-based and parametrised tests run too.
Run from the repo root:  python3 pyrestoolbox/tests/run_all_tests.py
Or with pytest directly: python3 -m pytest pyrestoolbox/tests/ -v

Exit code is non-zero if any test fails (pytest exit code passed through).
"""

import os
import sys


def main():
    tests_dir = os.path.dirname(os.path.abspath(__file__))
    # Ensure the repo root is importable so `import pyrestoolbox` resolves
    repo_root = os.path.dirname(os.path.dirname(tests_dir))
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    import pytest
    return pytest.main(["-q", tests_dir])


if __name__ == '__main__':
    sys.exit(main())
