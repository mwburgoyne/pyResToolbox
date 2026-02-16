#!/usr/bin/env python3
"""
Master test runner for pyResToolbox validation suite.
Run from project root: PYTHONPATH=/home/mark/projects python3 tests/run_all_tests.py
Or with pytest:        PYTHONPATH=/home/mark/projects python3 -m pytest tests/ -v
"""

import sys
import os
import importlib.util
import traceback

# Ensure project parent is on path for pyrestoolbox imports
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(project_root))  # parent of pyResToolbox

def load_test_module(filepath, module_name):
    """Load a test module from filepath"""
    spec = importlib.util.spec_from_file_location(module_name, filepath)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def run_module_tests(filepath, module_name):
    """Run all test_ functions in a test module and return (passed, failed, errors)"""
    try:
        mod = load_test_module(filepath, module_name)
    except Exception as e:
        print(f"  ERROR: Could not import {module_name}: {e}")
        traceback.print_exc()
        return 0, 1, [(module_name, str(e))]

    tests = [(k, v) for k, v in vars(mod).items() if k.startswith('test_') and callable(v)]
    passed, failed, errors = 0, 0, []

    for name, func in sorted(tests):
        try:
            func()
            passed += 1
            print(f"    PASS: {name}")
        except Exception as e:
            failed += 1
            errors.append((name, str(e)))
            print(f"    FAIL: {name}")
            print(f"          {e}")

    return passed, failed, errors


def main():
    print("=" * 70)
    print("pyResToolbox VALIDATION TEST SUITE")
    print("=" * 70)

    tests_dir = os.path.dirname(os.path.abspath(__file__))

    test_modules = [
        ('test_gas.py', 'Gas Module'),
        ('test_oil.py', 'Oil Module'),
        ('test_brine.py', 'Brine Module'),
        ('test_layer.py', 'Layer Module'),
        ('test_simtools.py', 'SimTools Module'),
        ('test_doc_examples.py', 'Documentation Examples'),
    ]

    total_passed = 0
    total_failed = 0
    all_errors = []

    for filename, display_name in test_modules:
        filepath = os.path.join(tests_dir, filename)
        module_name = filename.replace('.py', '')
        print(f"\n--- {display_name} ---")
        p, f, e = run_module_tests(filepath, module_name)
        total_passed += p
        total_failed += f
        all_errors.extend(e)
        print(f"  Subtotal: {p} passed, {f} failed")

    print(f"\n{'=' * 70}")
    print(f"TOTAL: {total_passed} passed, {total_failed} failed out of {total_passed + total_failed}")

    if all_errors:
        print(f"\nFailed tests:")
        for name, msg in all_errors:
            print(f"  - {name}: {msg}")

    print("=" * 70)
    return 1 if total_failed > 0 else 0


if __name__ == '__main__':
    sys.exit(main())
