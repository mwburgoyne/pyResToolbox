#!/usr/bin/env python3
"""Tests for the sensitivity module."""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from pyrestoolbox import sensitivity


def _simple_func(x=1.0, y=2.0):
    """Test function: returns x * y."""
    return x * y


def _dict_func(x=1.0, y=2.0):
    """Test function returning a dict."""
    return {'product': x * y, 'sum': x + y}


# ======================== sweep tests ========================

def test_sweep_basic():
    result = sensitivity.sweep(_simple_func, {'x': 2.0, 'y': 3.0}, 'x', [1.0, 2.0, 3.0])
    assert isinstance(result, sensitivity.SweepResult)
    assert result.param == 'x'
    assert len(result.values) == 3
    assert result.results == [3.0, 6.0, 9.0]


def test_sweep_with_result_key():
    result = sensitivity.sweep(_dict_func, {'x': 2.0, 'y': 3.0}, 'x', [1.0, 2.0],
                               result_key='product')
    assert result.results == [3.0, 6.0]


def test_sweep_single_value():
    result = sensitivity.sweep(_simple_func, {'x': 2.0, 'y': 3.0}, 'x', [5.0])
    assert result.results == [15.0]


def test_sweep_invalid_param():
    try:
        sensitivity.sweep(_simple_func, {'x': 2.0, 'y': 3.0}, 'z', [1.0])
        assert False, "Should have raised ValueError"
    except ValueError:
        pass


# ======================== tornado tests ========================

def test_tornado_basic():
    result = sensitivity.tornado(_simple_func, {'x': 2.0, 'y': 3.0},
                                 {'x': (1.0, 3.0), 'y': (2.0, 4.0)})
    assert isinstance(result, sensitivity.TornadoResult)
    assert result.base_result == 6.0
    assert len(result.entries) == 2
    # Both entries should be TornadoEntry
    for entry in result.entries:
        assert isinstance(entry, sensitivity.TornadoEntry)


def test_tornado_sensitivity_ordering():
    result = sensitivity.tornado(_simple_func, {'x': 2.0, 'y': 3.0},
                                 {'x': (1.0, 3.0), 'y': (1.0, 5.0)})
    # y range is larger (1-5 vs 1-3), so y should have higher sensitivity
    # x: low=1*3=3, high=3*3=9, sensitivity=|9-3|/6=1.0
    # y: low=2*1=2, high=2*5=10, sensitivity=|10-2|/6=1.333
    assert result.entries[0].param == 'y'
    assert result.entries[1].param == 'x'


def test_tornado_with_result_key():
    result = sensitivity.tornado(_dict_func, {'x': 2.0, 'y': 3.0},
                                 {'x': (1.0, 3.0)}, result_key='sum')
    assert result.base_result == 5.0
    assert len(result.entries) == 1
    assert result.entries[0].low_result == 4.0  # 1+3
    assert result.entries[0].high_result == 6.0  # 3+3


def test_tornado_invalid_param():
    try:
        sensitivity.tornado(_simple_func, {'x': 2.0}, {'z': (1.0, 3.0)})
        assert False, "Should have raised ValueError"
    except ValueError:
        pass


def test_tornado_zero_base():
    """When base_result is 0, sensitivity uses absolute difference."""
    def zero_func(x=0):
        return x
    result = sensitivity.tornado(zero_func, {'x': 0}, {'x': (-1.0, 1.0)})
    assert result.base_result == 0.0
    assert result.entries[0].sensitivity == 2.0


def test_tornado_single_param():
    result = sensitivity.tornado(_simple_func, {'x': 2.0, 'y': 3.0},
                                 {'x': (1.0, 4.0)})
    assert len(result.entries) == 1
    e = result.entries[0]
    assert e.param == 'x'
    assert e.low_value == 1.0
    assert e.high_value == 4.0
    assert e.low_result == 3.0   # 1*3
    assert e.high_result == 12.0  # 4*3


def test_sweep_result_attributes():
    result = sensitivity.sweep(_simple_func, {'x': 2.0, 'y': 3.0}, 'y', [1.0, 2.0, 3.0])
    assert result.param == 'y'
    assert result.values == [1.0, 2.0, 3.0]
    assert result.results == [2.0, 4.0, 6.0]


def test_tornado_entry_sensitivity_value():
    result = sensitivity.tornado(_simple_func, {'x': 5.0, 'y': 2.0},
                                 {'x': (4.0, 6.0)})
    e = result.entries[0]
    # base = 10, low = 8, high = 12, sensitivity = |12-8|/10 = 0.4
    assert abs(e.sensitivity - 0.4) < 1e-10
