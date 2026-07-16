"""Tests for the pyrestoolbox MCP server.

Run from the repo root:
    PYTHONPATH=/home/mark/projects/pyResToolbox:/home/mark/projects/pyResToolbox/mcp \
        python3 -m pytest mcp/tests/ -q

Most tests exercise the tool functions directly (FastMCP's @tool decorator
returns the undecorated function). test_client_session_* run the full MCP
protocol over an in-process transport.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from pyrestoolbox_mcp import server as srv
from pyrestoolbox_mcp.serialize import to_jsonable

RTOL = 1e-6


# ----------------------------------------------------------------- serialize

def test_serialize_numpy_and_nonfinite():
    import numpy as np
    assert to_jsonable(np.float64(1.5)) == 1.5
    assert type(to_jsonable(np.float64(1.5))) is float
    assert to_jsonable(np.array([1.0, 2.0])) == [1.0, 2.0]
    assert to_jsonable(float('nan')) is None
    assert to_jsonable(float('inf')) is None
    assert to_jsonable({'a': (1, 2)}) == {'a': [1, 2]}


def test_serialize_dataframe():
    import pandas as pd
    df = pd.DataFrame({'a': [1, 2], 'b': [3.0, 4.0]})
    out = to_jsonable(df)
    assert out['columns'] == ['a', 'b']
    assert out['records'][0] == {'a': 1, 'b': 3.0}


# ------------------------------------------------------------ discovery tools

def test_list_functions_modules():
    mods = srv.list_functions()
    assert 'gas' in mods and 'dca' in mods and 'plyasunov' in mods


def test_list_functions_gas():
    fns = srv.list_functions('gas')
    assert 'gas_z' in fns and 'gas_hydrate' in fns
    assert fns['gas_z']  # has a one-line summary


def test_list_functions_unknown_module():
    with pytest.raises(ValueError, match='Unknown module'):
        srv.list_functions('nope')


def test_describe_function_and_class():
    d = srv.describe('gas.gas_z')
    assert d.startswith('function gas.gas_z(')
    assert 'zmethod' in d
    assert 'numpy._typing' not in d  # annotations stripped
    c = srv.describe('gas.GasPVT')
    assert c.startswith('class gas.GasPVT(')
    assert 'self' not in c.split('\n')[0]


def test_describe_rejects_bad_names():
    with pytest.raises(ValueError, match='module.function'):
        srv.describe('gas_z')
    with pytest.raises(ValueError, match='not a public name'):
        srv.describe('gas._resolve_methods')


# ------------------------------------------------------------------ call tool

def test_call_gas_z_scalar():
    r = srv.call('gas.gas_z', {'p': 2000, 'sg': 0.75, 'degf': 200})
    assert abs(r['result'] - 0.8682653947215557) / 0.8682653947215557 < RTOL


def test_call_gas_z_list_input():
    r = srv.call('gas.gas_z', {'p': [1000, 2000], 'sg': 0.75, 'degf': 200})
    assert isinstance(r['result'], list) and len(r['result']) == 2


def test_call_invalid_method_lists_options():
    with pytest.raises(ValueError, match='Valid options'):
        srv.call('gas.gas_z', {'p': 2000, 'sg': 0.75, 'degf': 200,
                               'zmethod': 'NOPE'})


def test_call_rejects_class():
    with pytest.raises(ValueError, match='is a class'):
        srv.call('gas.GasPVT', {'sg': 0.75})


def test_call_rejects_private():
    with pytest.raises(ValueError, match='not a public name'):
        srv.call('gas._resolve_methods', {})


# ----------------------------------------------------------- composite tools

def test_recommend_methods_h2_mandatory_bns():
    rec = srv.recommend_methods(h2=0.1)
    assert rec['zmethod']['recommended'] == 'BNS'
    assert rec['zmethod']['mandatory'] is True
    assert 'rationale' in rec['zmethod']


def test_co2_brine_props_matches_doc_example():
    out = srv.co2_brine_props(pres=5000, temp=275, ppm=30000)
    # brine.rst worked example values
    assert abs(out['bw'][0] - 1.1091672843736888) < 1e-6
    assert abs(out['x'][0] - 0.02431225) < 1e-6
    assert type(out['bw'][0]) is float  # json-clean


def test_fit_and_forecast_exponential():
    import numpy as np
    t = list(np.arange(0, 36.0))
    q = [1000 * np.exp(-0.05 * x) for x in t]
    out = srv.fit_and_forecast(t=t, q=q, t_end=120, method='exponential')
    assert out['fit']['method'] == 'exponential'
    assert abs(out['fit']['di'] - 0.05) < 1e-6
    assert out['fit']['r_squared'] > 0.9999
    assert len(out['forecast']['t']) == len(out['forecast']['q'])
    assert out['forecast']['eur'] > 0


# ------------------------------------------------- full MCP protocol session

@pytest.mark.anyio
async def test_client_session_tools_and_resources():
    try:
        from mcp.shared.memory import create_connected_server_and_client_session
    except ImportError:
        pytest.skip('mcp.shared.memory helper not available in this SDK version')

    async with create_connected_server_and_client_session(
            srv.mcp._mcp_server) as client:
        tools = await client.list_tools()
        names = {t.name for t in tools.tools}
        assert {'list_functions', 'describe', 'call', 'recommend_methods',
                'co2_brine_props', 'sw_brine_props',
                'fit_and_forecast'} <= names

        result = await client.call_tool(
            'call', {'name': 'gas.gas_z',
                     'arguments': {'p': 2000, 'sg': 0.75, 'degf': 200}})
        assert not result.isError
        assert '0.868265' in result.content[0].text

        resources = await client.list_resources()
        uris = {str(r.uri) for r in resources.resources}
        assert 'docs://gas' in uris and 'docs://index' in uris

        doc = await client.read_resource('docs://gas')
        assert 'gas_z' in doc.contents[0].text


@pytest.fixture
def anyio_backend():
    return 'asyncio'


# ------------------------------------------------------------------ CLI args

def test_cli_parser_defaults_and_http():
    p = srv.build_parser()
    args = p.parse_args([])
    assert args.transport == 'stdio'
    args = p.parse_args(['--transport', 'streamable-http',
                         '--host', '0.0.0.0', '--port', '9000'])
    assert (args.transport, args.host, args.port) == ('streamable-http',
                                                      '0.0.0.0', 9000)
