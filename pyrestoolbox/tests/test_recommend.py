#!/usr/bin/env python3
"""Tests for the recommend module."""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from pyrestoolbox import recommend


# ======================== recommend_gas_methods tests ========================

def test_gas_clean():
    recs = recommend.recommend_gas_methods(sg=0.65)
    assert recs['zmethod'].recommended == 'DAK'
    assert recs['cmethod'].recommended == 'PMC'
    assert not recs['zmethod'].mandatory


def test_gas_h2():
    recs = recommend.recommend_gas_methods(sg=0.65, h2=0.05)
    assert recs['zmethod'].recommended == 'BNS'
    assert recs['zmethod'].mandatory
    assert recs['cmethod'].recommended == 'BNS'
    assert recs['cmethod'].mandatory


def test_gas_high_inerts():
    recs = recommend.recommend_gas_methods(sg=0.65, co2=0.30, n2=0.30)
    assert recs['zmethod'].recommended == 'BNS'
    assert not recs['zmethod'].mandatory  # recommended but not mandatory


def test_gas_moderate_co2():
    recs = recommend.recommend_gas_methods(sg=0.65, co2=0.15)
    assert recs['zmethod'].recommended == 'DAK'
    assert recs['cmethod'].recommended == 'PMC'
    assert 'BNS' in recs['zmethod'].alternatives


def test_gas_moderate_h2s():
    recs = recommend.recommend_gas_methods(sg=0.7, h2s=0.12)
    assert recs['zmethod'].recommended == 'DAK'
    assert 'BNS' in recs['zmethod'].alternatives


# ======================== recommend_oil_methods tests ========================

def test_oil_default():
    recs = recommend.recommend_oil_methods(api=35.0)
    assert recs['pbmethod'].recommended == 'VELAR'
    assert recs['rsmethod'].recommended == 'VELAR'
    assert recs['bomethod'].recommended == 'MCAIN'


def test_oil_heavy():
    recs = recommend.recommend_oil_methods(api=8.0)
    assert recs['pbmethod'].recommended == 'VELAR'
    assert 'heavy' in recs['pbmethod'].rationale.lower() or 'Heavy' in recs['pbmethod'].rationale


def test_oil_light():
    recs = recommend.recommend_oil_methods(api=55.0)
    assert recs['pbmethod'].recommended == 'VELAR'
    assert 'condensate' in recs['pbmethod'].rationale.lower() or 'Light' in recs['pbmethod'].rationale


# ======================== recommend_vlp_method tests ========================

def test_vlp_vertical():
    recs = recommend.recommend_vlp_method(deviation=0)
    assert recs['vlp_method'].recommended == 'BB'
    assert len(recs['vlp_method'].alternatives) == 3  # HB, WG, GRAY all OK


def test_vlp_deviated():
    recs = recommend.recommend_vlp_method(deviation=60)
    assert recs['vlp_method'].recommended == 'BB'
    assert 'WG' in recs['vlp_method'].alternatives
    assert 'HB' not in recs['vlp_method'].alternatives  # HB not suitable
    assert 'GRAY' not in recs['vlp_method'].alternatives  # Gray not suitable


def test_vlp_30deg_boundary():
    recs = recommend.recommend_vlp_method(deviation=30)
    assert len(recs['vlp_method'].alternatives) == 3  # All methods OK at exactly 30


# ======================== recommend_methods tests ========================

def test_master_gas_only():
    recs = recommend.recommend_methods(sg=0.65, h2=0.1)
    assert 'zmethod' in recs
    assert 'cmethod' in recs
    assert 'vlp_method' in recs
    assert 'pbmethod' not in recs  # No API specified


def test_master_with_oil():
    recs = recommend.recommend_methods(sg=0.65, api=30.0, deviation=45)
    assert 'zmethod' in recs
    assert 'pbmethod' in recs
    assert 'vlp_method' in recs
    # Deviated well: only BB and WG suitable
    assert len(recs['vlp_method'].alternatives) == 1  # Only WG
