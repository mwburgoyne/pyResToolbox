#!/usr/bin/env python3
"""
Tests for every documented code example in the RST documentation files.
Ensures that documented API calls actually work and return expected results.
Run with: PYTHONPATH=/home/mark/projects python3 -m pytest tests/test_doc_examples.py -v
Or standalone: PYTHONPATH=/home/mark/projects python3 tests/test_doc_examples.py
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
import pyrestoolbox.gas as gas
import pyrestoolbox.oil as oil
import pyrestoolbox.brine as brine
import pyrestoolbox.layer as layer
import pyrestoolbox.nodal as nodal
import pyrestoolbox.simtools as simtools
import pyrestoolbox.library as library

RTOL = 1e-4  # Relative tolerance for floating point comparisons

# =============================================================================
# Gas Module Documentation Examples (docs/gas.rst)
# =============================================================================

def test_doc_gas_z_dak_pmc():
    """gas.rst: gas_z with DAK/PMC for sg=0.68"""
    result = gas.gas_z(p=2350, sg=0.68, degf=180, zmethod='DAK', cmethod='PMC')
    assert isinstance(result, float)
    assert abs(result - 0.8785399927100872) / 0.8785399927100872 < RTOL

def test_doc_gas_z_bur_co2():
    """gas.rst: gas_z for pure CO2 with BUR method"""
    result = gas.gas_z(p=2350, sg=0.68, degf=180, co2=1.0, zmethod='BUR', cmethod='BUR')
    assert isinstance(result, float)
    assert abs(result - 0.5258309021348752) / 0.5258309021348752 < RTOL

def test_doc_gas_sg():
    """gas.rst: gas_sg for mixture with H2"""
    result = gas.gas_sg(hc_mw=19.0, co2=0.05, h2s=0.10, n2=0, h2=0.20)
    assert isinstance(result, float)
    assert abs(result - 0.6338246461857093) / 0.6338246461857093 < RTOL

def test_doc_gas_z_bur_mixture():
    """gas.rst: gas_z BUR for complex mixture"""
    gsg = gas.gas_sg(hc_mw=19.0, co2=0.05, h2s=0.10, n2=0, h2=0.20)
    result = gas.gas_z(p=2350, sg=gsg, degf=180, co2=0.05, h2s=0.10, n2=0, h2=0.20, zmethod='BUR', cmethod='BUR')
    assert isinstance(result, float)
    assert abs(result - 0.9048153036714465) / 0.9048153036714465 < RTOL

def test_doc_gas_tc_pc_pmc():
    """gas.rst: gas_tc_pc with PMC defaults"""
    tc, pc = gas.gas_tc_pc(sg=0.7, co2=0.15)
    assert abs(tc - 363.9387708314338) / 363.9387708314338 < RTOL
    assert abs(pc - 738.3190067714969) / 738.3190067714969 < RTOL

def test_doc_gas_tc_pc_sut_fixed_tc():
    """gas.rst: gas_tc_pc with SUT and fixed tc"""
    tc, pc = gas.gas_tc_pc(sg=0.7, co2=0.15, tc=365, cmethod='SUT')
    assert tc == 365  # Should be returned unchanged
    assert abs(pc - 709.2356299485114) / 709.2356299485114 < RTOL

def test_doc_gas_z_n2_co2():
    """gas.rst: gas_z with N2 and CO2"""
    result = gas.gas_z(p=1000, sg=0.75, degf=160, n2=0.02, co2=0.17)
    assert isinstance(result, float)
    assert abs(result - 0.9138558878125714) / 0.9138558878125714 < RTOL

def test_doc_gas_z_hy():
    """gas.rst: gas_z with HY method"""
    result = gas.gas_z(p=1000, sg=0.75, degf=160, n2=0.02, co2=0.17, zmethod='HY')
    assert isinstance(result, float)
    assert abs(result - 0.9142136711443208) / 0.9142136711443208 < RTOL

def test_doc_gas_z_sut_array():
    """gas.rst: gas_z array with SUT"""
    result = gas.gas_z(p=[1000, 2000], sg=0.75, degf=160, cmethod='SUT', n2=0.02, co2=0.17)
    assert isinstance(result, np.ndarray)
    assert len(result) == 2
    expected = np.array([0.91900003, 0.87160514])
    np.testing.assert_allclose(result, expected, rtol=RTOL)

def test_doc_gas_ug_hy_sut():
    """gas.rst: gas_ug with HY/SUT"""
    result = gas.gas_ug(p=1000, sg=0.75, degf=180, zmethod='HY', cmethod='SUT')
    assert isinstance(result, float)
    assert abs(result - 0.014118890100250796) / 0.014118890100250796 < RTOL

def test_doc_gas_ug_default():
    """gas.rst: gas_ug with defaults"""
    result = gas.gas_ug(p=1000, sg=0.75, degf=180)
    assert isinstance(result, float)
    assert abs(result - 0.014110092961853301) / 0.014110092961853301 < RTOL

def test_doc_gas_cg_scalar():
    """gas.rst: gas_cg scalar"""
    result = gas.gas_cg(p=2000, sg=0.68, degf=120, co2=0.05)
    assert isinstance(result, float)
    assert abs(result - 0.0005374854430839333) / 0.0005374854430839333 < RTOL

def test_doc_gas_cg_array():
    """gas.rst: gas_cg array"""
    result = gas.gas_cg(p=np.array([1000, 2000]), sg=0.68, degf=120, co2=0.05)
    assert isinstance(result, np.ndarray)
    expected = np.array([0.00110369, 0.00053749])
    np.testing.assert_allclose(result, expected, rtol=RTOL)

def test_doc_gas_bg_scalar():
    """gas.rst: gas_bg scalar"""
    result = gas.gas_bg(p=3000, sg=0.78, degf=240)
    assert isinstance(result, float)
    assert abs(result - 0.005927563975073749) / 0.005927563975073749 < RTOL

def test_doc_gas_bg_array_inverse():
    """gas.rst: 1/gas_bg array"""
    result = 1 / gas.gas_bg(p=[3000, 5000], sg=0.78, degf=240)
    assert isinstance(result, np.ndarray)
    expected = np.array([168.70336688, 249.54573283])
    np.testing.assert_allclose(result, expected, rtol=RTOL)

def test_doc_gas_den():
    """gas.rst: gas_den with impurities"""
    result = gas.gas_den(p=2000, sg=0.75, degf=150, zmethod='HY', cmethod='SUT', n2=0.02, co2=0.15, h2s=0.02)
    assert isinstance(result, float)
    assert abs(result - 7.736656004563576) / 7.736656004563576 < RTOL

def test_doc_gas_water_content():
    """gas.rst: gas_water_content"""
    result = gas.gas_water_content(p=1500, degf=165)
    assert isinstance(result, float)
    assert abs(result - 0.6474226409378979) / 0.6474226409378979 < RTOL

def test_doc_gas_water_content_salinity():
    """gas.rst: gas_water_content with salinity"""
    result = gas.gas_water_content(p=1500, degf=165, salinity=5)
    assert isinstance(result, float)
    assert result < gas.gas_water_content(p=1500, degf=165)  # salinity reduces water content
    assert abs(result - 0.628635730743162) / 0.628635730743162 < RTOL

def test_doc_gas_sg_pure_methane():
    """gas.rst: gas_sg for pure methane"""
    result = gas.gas_sg(hc_mw=16.043, co2=0, h2s=0, n2=0, h2=0)
    assert isinstance(result, float)
    assert abs(result - 0.5537797721781152) / 0.5537797721781152 < RTOL

def test_bns_recommended_for_high_inerts():
    """gas.rst: BNS method is strongly recommended for high CO2/H2S/N2/H2 gases.

    Standard correlations (DAK/HY/WYW with PMC/SUT) were developed for sweet natural
    gas and become unreliable at high impurity concentrations. The BNS method (tuned
    5-component Peng-Robinson EOS) handles the full range from pure hydrocarbon to
    100% inerts. This test demonstrates the BNS method working correctly for:
    - Pure CO2 (z should be << 1 at moderate pressure, sub-critical)
    - Pure H2S (should produce valid z-factor)
    - Pure N2 (should be near 1.0 at moderate conditions)
    - High H2 blend (auto-selects BNS when h2 > 0)
    - Mixed high-inert gas (50%+ inerts)

    Standard methods will either fail or produce unreliable results for these cases.
    """
    # Pure CO2 at 2000 psia, 150 degF - BNS handles supercritical CO2 correctly
    z_co2 = gas.gas_z(p=2000, sg=0.65, degf=150, co2=1.0, zmethod='BNS', cmethod='BNS')
    assert 0.2 < z_co2 < 0.8, f"Pure CO2 Z-factor should be well below 1.0, got {z_co2}"

    # Pure H2S at 2000 psia, 200 degF
    z_h2s = gas.gas_z(p=2000, sg=0.65, degf=200, h2s=1.0, zmethod='BNS', cmethod='BNS')
    assert 0.2 < z_h2s < 0.9, f"Pure H2S Z-factor should be well below 1.0, got {z_h2s}"

    # Pure N2 at 2000 psia, 150 degF - N2 is supercritical, z near 1
    z_n2 = gas.gas_z(p=2000, sg=0.65, degf=150, n2=1.0, zmethod='BNS', cmethod='BNS')
    assert 0.8 < z_n2 < 1.2, f"Pure N2 Z-factor should be near 1.0, got {z_n2}"

    # High H2 blend (20% H2) - auto-selects BNS when h2 > 0
    gsg = gas.gas_sg(hc_mw=19.0, co2=0.05, h2s=0.10, n2=0, h2=0.20)
    z_h2_blend = gas.gas_z(p=2000, sg=gsg, degf=180, co2=0.05, h2s=0.10, n2=0, h2=0.20)
    assert 0.5 < z_h2_blend < 1.1, f"H2 blend Z-factor should be reasonable, got {z_h2_blend}"

    # Mixed high-inert gas: 30% CO2, 10% H2S, 15% N2 = 55% inerts
    gsg_mix = gas.gas_sg(hc_mw=18.0, co2=0.30, h2s=0.10, n2=0.15, h2=0)
    z_mix = gas.gas_z(p=3000, sg=gsg_mix, degf=200, co2=0.30, h2s=0.10, n2=0.15, zmethod='BNS', cmethod='BNS')
    assert 0.3 < z_mix < 1.0, f"High-inert mix Z-factor should be valid, got {z_mix}"

    # Verify BNS viscosity also works for these cases (LBC correlation)
    ug_co2 = gas.gas_ug(p=2000, sg=0.65, degf=150, co2=1.0, zmethod='BNS', cmethod='BNS')
    assert ug_co2 > 0, "BNS viscosity for pure CO2 must be positive"

    ug_h2_blend = gas.gas_ug(p=2000, sg=gsg, degf=180, co2=0.05, h2s=0.10, n2=0, h2=0.20, zmethod='BNS', cmethod='BNS')
    assert ug_h2_blend > 0, "BNS viscosity for H2 blend must be positive"

def test_doc_gas_ponz2p_scalar():
    """gas.rst: gas_ponz2p scalar"""
    result = gas.gas_ponz2p(poverz=2500, sg=0.75, degf=165)
    assert isinstance(result, float)
    assert abs(result - 2081.5489292144775) / 2081.5489292144775 < RTOL

def test_doc_gas_ponz2p_array():
    """gas.rst: gas_ponz2p array"""
    result = gas.gas_ponz2p(poverz=[2500, 5000], sg=0.75, degf=165)
    assert isinstance(result, np.ndarray)
    expected = np.array([2081.54892921, 4856.97983205])
    np.testing.assert_allclose(result, expected, rtol=RTOL)

def test_doc_gas_grad2sg():
    """gas.rst: gas_grad2sg"""
    result = gas.gas_grad2sg(grad=0.0657, p=2500, degf=175)
    assert isinstance(result, float)
    assert abs(result - 0.7495803994806547) / 0.7495803994806547 < RTOL

def test_doc_gas_dmp_positive():
    """gas.rst: gas_dmp positive (p1 < p2)"""
    result = gas.gas_dmp(p1=1000, p2=2000, degf=185, sg=0.78, zmethod='HY', cmethod='SUT', n2=0.05, co2=0.1, h2s=0.02)
    assert isinstance(result, float)
    assert result > 0
    assert abs(result - 213690308.9907268) / 213690308.9907268 < RTOL

def test_doc_gas_dmp_negative():
    """gas.rst: gas_dmp negative (p1 > p2) with fixed tc/pc"""
    result = gas.gas_dmp(p1=2000, p2=1000, degf=185, sg=0.78, tc=371, pc=682)
    assert isinstance(result, float)
    assert result < 0
    assert abs(result - (-213713909.36339885)) / 213713909.36339885 < RTOL

def test_doc_gas_fws_sg():
    """gas.rst: gas_fws_sg"""
    result = gas.gas_fws_sg(sg_g=0.855, cgr=30, api_st=53)
    assert isinstance(result, float)
    assert abs(result - 0.937116010334538) / 0.937116010334538 < RTOL

def test_doc_gas_rate_radial_scalar():
    """gas.rst: gas_rate_radial scalar"""
    result = gas.gas_rate_radial(k=5, h=50, pr=2000, pwf=750, r_w=0.3, r_ext=1500, degf=180, sg=0.75, D=0.01, S=5)
    assert isinstance(result, float)
    assert abs(result - 2078.9101970773477) / 2078.9101970773477 < RTOL

def test_doc_gas_rate_radial_array():
    """gas.rst: gas_rate_radial array"""
    result = gas.gas_rate_radial(k=1, h=50, pr=[2000, 1000], pwf=750, r_w=0.3, r_ext=1500, degf=180, sg=0.75, D=0.01, S=5)
    assert isinstance(result, np.ndarray)
    expected = np.array([704.29202227, 135.05317439])
    np.testing.assert_allclose(result, expected, rtol=RTOL)

def test_doc_gas_rate_linear_scalar():
    """gas.rst: gas_rate_linear scalar"""
    result = gas.gas_rate_linear(k=0.1, area=50, length=200, pr=2000, pwf=250, degf=180, sg=0.8)
    assert isinstance(result, float)
    assert abs(result - 8.202200317597859) / 8.202200317597859 < RTOL

def test_doc_gas_rate_linear_array():
    """gas.rst: gas_rate_linear array"""
    result = gas.gas_rate_linear(k=0.1, area=50, length=200, pr=[2000, 1000, 500], pwf=250, degf=180, sg=0.8)
    assert isinstance(result, np.ndarray)
    expected = np.array([8.20220032, 2.10691337, 0.42685002])
    np.testing.assert_allclose(result, expected, rtol=RTOL)

# =============================================================================
# Oil Module Documentation Examples (docs/oil.rst)
# =============================================================================

def test_doc_oil_ja_sg():
    """oil.rst: oil_ja_sg"""
    result = oil.oil_ja_sg(mw=150, ja=0.5)
    assert isinstance(result, float)
    assert abs(result - 0.8583666666666667) / 0.8583666666666667 < RTOL

def test_doc_oil_twu_props():
    """oil.rst: oil_twu_props"""
    result = oil.oil_twu_props(mw=225, ja=0.5)
    assert isinstance(result, tuple)
    assert len(result) == 5
    expected = (0.8954444444444445, 1068.3961103813851, 1422.4620493584146, 264.23402773211745, 13.498328588856445)
    for r, e in zip(result, expected):
        assert abs(float(r) - e) / abs(e) < RTOL

def test_doc_oil_rs_st():
    """oil.rst: oil_rs_st"""
    result = oil.oil_rs_st(psp=114.7, degf_sp=80, api=38)
    assert isinstance(result, float)
    assert abs(result - 4.176458005559282) / 4.176458005559282 < RTOL

def test_doc_oil_pbub_valmc():
    """oil.rst: oil_pbub with VALMC default"""
    result = oil.oil_pbub(api=43, degf=185, rsb=2350, sg_g=0.72)
    assert isinstance(result, float)
    assert abs(result - 5199.2406069808885) / 5199.2406069808885 < RTOL

def test_doc_oil_pbub_stan():
    """oil.rst: oil_pbub with Standing via sg_sp"""
    result = oil.oil_pbub(api=43, degf=185, rsb=2350, sg_sp=0.72, pbmethod='STAN')
    assert isinstance(result, float)
    assert abs(result - 6390.281894698239) / 6390.281894698239 < RTOL

def test_doc_oil_pbub_class_object():
    """oil.rst: oil_pbub using class object"""
    result = oil.oil_pbub(api=43, degf=185, rsb=2350, sg_g=0.72, pbmethod=oil.pb_method.STAN)
    assert isinstance(result, float)
    assert result > 0

def test_doc_oil_pbub_metric():
    """oil.rst: oil_pbub metric example"""
    result = oil.oil_pbub(api=43, degf=85, rsb=418.8, sg_g=0.72, metric=True)
    assert abs(result - 358.5685338835858) / 358.5685338835858 < RTOL

def test_doc_oil_rs_bub():
    """oil.rst: oil_rs_bub"""
    result = oil.oil_rs_bub(api=43, degf=185, pb=5179.5, sg_sp=0.72)
    assert isinstance(result, float)
    assert abs(result - 1872.666133282599) / 1872.666133282599 < RTOL

def test_doc_oil_rs_with_pb_rsb():
    """oil.rst: oil_rs with both pb and rsb"""
    result = oil.oil_rs(api=43, degf=185, sg_sp=0.72, p=3000, pb=5179.5, rsb=2370)
    assert isinstance(result, float)
    assert abs(result - 1017.9424383646037) / 1017.9424383646037 < RTOL

def test_doc_oil_rs_with_rsb_only():
    """oil.rst: oil_rs with rsb only"""
    result = oil.oil_rs(api=43, degf=185, sg_sp=0.72, p=3000, rsb=2370)
    assert isinstance(result, float)
    assert abs(result - 1010.0669567201218) / 1010.0669567201218 < RTOL

def test_doc_oil_rs_with_pb_only():
    """oil.rst: oil_rs with pb only"""
    result = oil.oil_rs(api=43, degf=185, sg_sp=0.72, p=3000, pb=5180)
    assert isinstance(result, float)
    assert abs(result - 804.2857187814161) / 804.2857187814161 < RTOL

def test_doc_oil_rs_stan():
    """oil.rst: oil_rs with Standing method"""
    result = oil.oil_rs(api=43, degf=185, sg_sp=0.72, p=3000, pb=5180, rsmethod='STAN')
    assert isinstance(result, float)
    assert abs(result - 947.1133546937306) / 947.1133546937306 < RTOL

def test_doc_oil_co_above_pb():
    """oil.rst: oil_co above bubble point"""
    result = oil.oil_co(p=4500, api=47, degf=180, sg_sp=0.72, rsb=2750)
    assert isinstance(result, float)
    assert abs(result - 0.0007587726853322233) / 0.0007587726853322233 < RTOL

def test_doc_oil_co_below_pb():
    """oil.rst: oil_co below bubble point"""
    result = oil.oil_co(p=2000, api=47, degf=180, sg_sp=0.72, rsb=2750, pb=4945)
    assert isinstance(result, float)
    assert abs(result - 0.0009245540028053584) / 0.0009245540028053584 < RTOL

def test_doc_oil_deno():
    """oil.rst: oil_deno"""
    result = oil.oil_deno(p=2000, degf=165, rs=1000, rsb=2000, sg_g=0.72, api=38)
    assert isinstance(result, float)
    assert abs(result - 40.98349866963842) / 40.98349866963842 < RTOL

def test_doc_oil_bo_mcain():
    """oil.rst: oil_bo with McCain default"""
    result = oil.oil_bo(p=2000, pb=3000, degf=165, rs=1000, rsb=2000, sg_o=0.8, sg_g=0.68)
    assert isinstance(result, float)
    assert abs(result - 1.5075107735318138) / 1.5075107735318138 < RTOL

def test_doc_oil_bo_stan():
    """oil.rst: oil_bo with Standing"""
    result = oil.oil_bo(p=2000, pb=3000, degf=165, rs=1000, rsb=2000, sg_o=0.8, sg_g=0.68, bomethod='STAN')
    assert isinstance(result, float)
    assert abs(result - 1.5393786735904431) / 1.5393786735904431 < RTOL

def test_doc_oil_viso():
    """oil.rst: oil_viso"""
    result = oil.oil_viso(p=2000, api=38, degf=165, pb=3500, rs=1000)
    assert isinstance(result, float)
    assert abs(result - 0.416858469042502) / 0.416858469042502 < RTOL

def test_doc_sg_evolved_gas():
    """oil.rst: sg_evolved_gas"""
    result = oil.sg_evolved_gas(p=2000, degf=185, rsb=2370, api=43, sg_sp=0.72)
    assert isinstance(result, float)
    assert abs(result - 0.7872810977386344) / 0.7872810977386344 < RTOL

def test_doc_sg_st_gas():
    """oil.rst: sg_st_gas"""
    result = oil.sg_st_gas(114.7, rsp=1500, api=42, sg_sp=0.72, degf_sp=80)
    assert isinstance(result, float)
    assert abs(result - 1.1923932340625523) / 1.1923932340625523 < RTOL

def test_doc_sgg_wt_avg():
    """oil.rst: sgg_wt_avg"""
    result = oil.sgg_wt_avg(sg_sp=0.72, rsp=1000, sg_st=1.1, rst=5)
    assert isinstance(result, float)
    assert abs(result - 0.7218905472636816) / 0.7218905472636816 < RTOL

def test_doc_oil_api():
    """oil.rst: oil_api"""
    result = oil.oil_api(sg_value=0.82)
    assert isinstance(result, float)
    assert abs(result - 41.0609756097561) / 41.0609756097561 < RTOL

def test_doc_oil_sg():
    """oil.rst: oil_sg"""
    result = oil.oil_sg(api_value=45)
    assert isinstance(result, float)
    assert abs(result - 0.8016997167138811) / 0.8016997167138811 < RTOL

def test_doc_oil_rate_radial_scalar():
    """oil.rst: oil_rate_radial scalar with Vogel"""
    result = oil.oil_rate_radial(k=20, h=20, pr=1500, pwf=250, r_w=0.3, r_ext=1500, uo=0.8, bo=1.4, vogel=True, pb=1800)
    assert isinstance(result, float)
    assert abs(result - 213.8147848023242) / 213.8147848023242 < RTOL

def test_doc_oil_rate_radial_array():
    """oil.rst: oil_rate_radial array with Vogel"""
    result = oil.oil_rate_radial(k=20, h=20, pr=[1500, 2000], pwf=250, r_w=0.3, r_ext=1500, uo=0.8, bo=1.4, vogel=True, pb=1800)
    assert isinstance(result, np.ndarray)
    expected = np.array([213.8147848, 376.58731835])
    np.testing.assert_allclose(result, expected, rtol=RTOL)

def test_doc_oil_rate_linear_scalar():
    """oil.rst: oil_rate_linear scalar"""
    result = oil.oil_rate_linear(k=0.1, area=15000, pr=3000, pwf=500, length=500, uo=0.4, bo=1.5)
    assert isinstance(result, float)
    assert abs(result - 14.08521246363274) / 14.08521246363274 < RTOL

def test_doc_oil_rate_linear_array():
    """oil.rst: oil_rate_linear array"""
    result = oil.oil_rate_linear(k=[0.1, 1, 5, 10], area=15000, pr=3000, pwf=500, length=500, uo=0.4, bo=1.5)
    assert isinstance(result, np.ndarray)
    expected = np.array([14.08521246, 140.85212464, 704.26062318, 1408.52124636])
    np.testing.assert_allclose(result, expected, rtol=RTOL)

def test_doc_make_bot_og():
    """oil.rst: make_bot_og returns correct structure"""
    results = oil.make_bot_og(pvto=False, pi=4000, api=38, degf=175, sg_g=0.68, pmax=5500, pb=4500, nrows=10, export=False)
    assert isinstance(results, dict)
    for key in ['bot', 'deno', 'deng', 'denw', 'cw', 'uw', 'pb', 'rsb', 'rsb_scale', 'usat']:
        assert key in results, f"Missing key: {key}"
    assert results['pb'] == 4500
    assert results['bot'].shape[0] == 10

# =============================================================================
# Brine Module Documentation Examples (docs/brine.rst)
# =============================================================================

def test_doc_make_pvtw_table():
    """brine.rst: make_pvtw_table"""
    result = brine.make_pvtw_table(pi=3000, degf=200, wt=0, ch4_sat=0)
    assert isinstance(result, dict)
    assert abs(result['bw_ref'] - 1.027589195773527) / 1.027589195773527 < RTOL
    assert abs(result['cw_ref'] - 3.0887176266534516e-06) / 3.0887176266534516e-06 < RTOL
    assert abs(result['visw_ref'] - 0.3083544960904146) / 0.3083544960904146 < RTOL
    assert 'table' in result

def test_doc_brine_props():
    """brine.rst: brine_props"""
    bw, lsg, visw, cw, rsw = brine.brine_props(p=160, degf=135, wt=1.5, ch4_sat=1.0)
    assert abs(bw - 1.0152007040056148) / 1.0152007040056148 < RTOL
    assert abs(lsg - 0.9950108179684295) / 0.9950108179684295 < RTOL
    assert abs(visw - 0.4994004662758671) / 0.4994004662758671 < RTOL
    assert abs(cw - 0.0001539690974662865) / 0.0001539690974662865 < RTOL
    assert abs(rsw - 1.2540982731813703) / 1.2540982731813703 < RTOL

def test_doc_co2_brine_field():
    """brine.rst: CO2_Brine_Mixture field units"""
    mix = brine.CO2_Brine_Mixture(pres=5000, temp=275, ppm=30000, metric=False)
    assert isinstance(mix.bw, list)
    assert len(mix.bw) == 3
    assert abs(float(mix.bw[0]) - 1.108578337107381) / 1.108578337107381 < RTOL
    assert isinstance(mix.x, np.ndarray)
    assert abs(mix.x[0] - 0.02431225) / 0.02431225 < RTOL

def test_doc_co2_brine_metric():
    """brine.rst: CO2_Brine_Mixture metric units"""
    mix = brine.CO2_Brine_Mixture(pres=175, temp=85)
    assert abs(mix.Rs - 24.743651168969475) / 24.743651168969475 < RTOL

def test_doc_make_pvtw_table():
    """brine.rst: make_pvtw_table"""
    result = brine.make_pvtw_table(pi=3000, degf=200, wt=0, ch4_sat=0)
    assert isinstance(result, dict)
    for key in ['table', 'pref', 'bw_ref', 'cw_ref', 'visw_ref', 'rsw_ref', 'den_ref']:
        assert key in result, f"Missing key: {key}"
    assert abs(result['bw_ref'] - 1.027589195773527) / 1.027589195773527 < RTOL

def test_doc_sw_pure_co2_field():
    """brine.rst: SoreideWhitson pure CO2 field units"""
    mix = brine.SoreideWhitson(pres=5000, temp=275, ppm=30000, y_CO2=1.0, metric=False)
    assert isinstance(mix.bDen, list) and len(mix.bDen) == 3
    assert abs(mix.bDen[0] - 0.9733769457162755) / 0.9733769457162755 < RTOL
    assert abs(mix.Rs['CO2'] - 140.90858294709142) / 140.90858294709142 < RTOL
    assert abs(mix.bw[0] - 1.0968991160573776) / 1.0968991160573776 < RTOL

def test_doc_sw_pure_ch4_field():
    """brine.rst: SoreideWhitson pure CH4 field units"""
    mix = brine.SoreideWhitson(pres=5000, temp=275, ppm=30000, y_CO2=0, sg=0.554, metric=False)
    assert abs(mix.Rs['CH4'] - 21.414423540331008) / 21.414423540331008 < RTOL
    assert abs(mix.bDen[0] - 0.9641137202631425) / 0.9641137202631425 < RTOL

def test_doc_sw_mixed_gas_metric():
    """brine.rst: SoreideWhitson mixed gas metric"""
    mix = brine.SoreideWhitson(pres=200, temp=80, ppm=10000, y_CO2=0.1, y_H2S=0.05, sg=0.7, metric=True)
    assert abs(mix.Rs_total - 6.32875877743837) / 6.32875877743837 < RTOL
    assert abs(mix.bDen[0] - 0.9854845215724698) / 0.9854845215724698 < RTOL
    assert 'CO2' in mix.gas_comp and 'H2S' in mix.gas_comp and 'CH4' in mix.gas_comp

def test_doc_sw_co2_freshwater_cw_sat():
    """brine.rst: SoreideWhitson pure CO2 freshwater with Cf_sat"""
    mix = brine.SoreideWhitson(pres=175, temp=85, ppm=0, y_CO2=1.0, metric=True, cw_sat=True)
    assert abs(mix.Rs_total - 24.188037633302223) / 24.188037633302223 < RTOL
    assert abs(mix.Cf_sat - 0.00016012590421810821) / 0.00016012590421810821 < RTOL
    assert isinstance(mix.water_content, dict)
    assert abs(mix.water_content['stb_mmscf'] - 1.923030543083137) / 1.923030543083137 < RTOL

# =============================================================================
# Layer Module Documentation Examples (docs/layer.rst)
# =============================================================================

def test_doc_lorenz2b_lang():
    """layer.rst: lorenz2b with Langmuir"""
    result = layer.lorenz2b(0.75, lrnz_method='LANG')
    assert abs(result - 16.139518537603912) / 16.139518537603912 < RTOL

def test_doc_lorenz2b_exp():
    """layer.rst: lorenz2b with Exponential default"""
    result = layer.lorenz2b(0.75)
    assert abs(result - 7.978108090962671) / 7.978108090962671 < RTOL

def test_doc_lorenzfromb_lang():
    """layer.rst: lorenzfromb with Langmuir"""
    result = layer.lorenzfromb(16.139518537603912, lrnz_method='LANG')
    assert abs(result - 0.750000182307895) / 0.750000182307895 < RTOL

def test_doc_lorenzfromb_exp():
    """layer.rst: lorenzfromb with Exponential default"""
    result = layer.lorenzfromb(7.978108090962671)
    assert abs(result - 0.7500000108799212) / 0.7500000108799212 < RTOL

def test_doc_lorenz_from_flow_fraction():
    """layer.rst: lorenz_from_flow_fraction"""
    result = layer.lorenz_from_flow_fraction(kh_frac=0.6, phih_frac=0.15)
    assert abs(result - 0.6759312029093838) / 0.6759312029093838 < RTOL

def test_doc_lorenz_2_flow_frac():
    """layer.rst: lorenz_2_flow_frac"""
    result = layer.lorenz_2_flow_frac(lorenz=0.6759312029093838, phih_frac=0.15)
    assert abs(result - 0.6000001346893536) / 0.6000001346893536 < RTOL

def test_doc_lorenz_2_layers_nlayers():
    """layer.rst: lorenz_2_layers with nlayers"""
    result = layer.lorenz_2_layers(lorenz=0.67, nlayers=5, k_avg=10, shuffle=False)
    assert isinstance(result, np.ndarray)
    assert len(result) == 5
    # With shuffle=False, should be sorted descending
    assert all(result[i] >= result[i+1] for i in range(len(result)-1))
    # Average should be close to 10
    assert abs(np.mean(result) - 10) / 10 < 0.01
    expected = np.array([34.9323596, 10.58944038, 3.21009656, 0.9731128, 0.29499066])
    np.testing.assert_allclose(result, expected, rtol=RTOL)

def test_doc_lorenz_2_layers_phi_h_fracs():
    """layer.rst: lorenz_2_layers with phi_h_fracs"""
    result = layer.lorenz_2_layers(lorenz=0.67, k_avg=10, phi_h_fracs=[0.05, 0.5])
    assert isinstance(result, np.ndarray)
    assert len(result) == 3
    expected = np.array([51.72990694, 14.12556056, 0.77938749])
    np.testing.assert_allclose(result, expected, rtol=RTOL)

# =============================================================================
# SimTools Module Documentation Examples (docs/simtools.rst)
# =============================================================================

def test_doc_rel_perm_table_sgof():
    """simtools.rst: rel_perm_table SGOF with LET"""
    df = simtools.rel_perm_table(rows=25, krtable='SGOF', krfamily='LET', kromax=1, krgmax=1, swc=0.2, sorg=0.15, Lo=2.5, Eo=1.25, To=1.75, Lg=1.2, Eg=1.5, Tg=2.0)
    assert 'Sg' in df.columns
    assert 'Krgo' in df.columns
    assert 'Krog' in df.columns
    assert df.shape[0] >= 25

def test_doc_rel_perm_table_swof():
    """simtools.rst: rel_perm_table SWOF with Corey"""
    df = simtools.rel_perm_table(rows=25, krtable='SWOF', kromax=1, krwmax=0.25, swc=0.15, swcr=0.2, sorw=0.15, no=2.5, nw=1.5)
    assert 'Sw' in df.columns
    assert 'Krwo' in df.columns
    assert 'Krow' in df.columns
    assert df.shape[0] == 25

def test_doc_rr_solver():
    """simtools.rst: rr_solver"""
    n_it, yi, xi, V, L = simtools.rr_solver(zi=np.array([0.7, 0.15, 0.1, 0.05]), ki=np.array([50, 5, 0.5, 0.01]))
    assert n_it == 6
    assert abs(V - 0.9440279802330239) / 0.9440279802330239 < RTOL
    assert abs(L - 0.05597201976697608) / 0.05597201976697608 < RTOL
    expected_yi = np.array([0.7406252, 0.1570315, 0.09469948, 0.00764382])
    expected_xi = np.array([0.0148125, 0.0314063, 0.18939896, 0.76438224])
    np.testing.assert_allclose(yi, expected_yi, rtol=RTOL)
    np.testing.assert_allclose(xi, expected_xi, rtol=RTOL)

# =============================================================================
# Library Module Documentation Examples (docs/library.rst)
# =============================================================================

def test_doc_library_prop_pc():
    """library.rst: library.prop for CH4 Pc_psia"""
    result = library.prop(comp='CH4', prop='Pc_psia')
    assert result == 667.029

def test_doc_library_prop_vtran_pr79():
    """library.rst: library.prop for C3 VTran PR79"""
    result = library.prop(comp='C3', prop='VTran')
    assert result == -0.06381

def test_doc_library_prop_vtran_srk():
    """library.rst: library.prop for C3 VTran SRK"""
    result = library.prop(comp='C3', prop='VTran', model='SRK')
    assert result == 0.09075

def test_doc_library_components():
    """library.rst: library.components is a list with CH4"""
    assert isinstance(library.components, list)
    assert 'CH4' in library.components

def test_doc_library_property_list():
    """library.rst: library.property_list"""
    assert isinstance(library.property_list, list)
    assert 'MW' in library.property_list
    assert 'Tc_R' in library.property_list
    assert 'Pc_psia' in library.property_list

def test_doc_library_models():
    """library.rst: library.models"""
    assert library.models == ['PR79', 'PR77', 'SRK', 'RK']


# =============================================================================
# Nodal Module Documentation Examples (docs/nodal.rst)
# =============================================================================

def test_doc_nodal_gas_pvt_z():
    """nodal.rst: GasPVT z-factor"""
    gpvt = gas.GasPVT(sg=0.65, co2=0.1)
    result = gpvt.z(2000, 180)
    assert abs(result - 0.9026719828498643) / 0.9026719828498643 < RTOL

def test_doc_nodal_gas_pvt_viscosity():
    """nodal.rst: GasPVT viscosity"""
    gpvt = gas.GasPVT(sg=0.65, co2=0.1)
    result = gpvt.viscosity(2000, 180)
    assert abs(result - 0.016666761192334678) / 0.016666761192334678 < RTOL

def test_doc_nodal_gas_pvt_density():
    """nodal.rst: GasPVT density"""
    gpvt = gas.GasPVT(sg=0.65, co2=0.1)
    result = gpvt.density(2000, 180)
    assert abs(result - 6.077743278791424) / 6.077743278791424 < RTOL

def test_doc_nodal_gas_pvt_bg():
    """nodal.rst: GasPVT bg"""
    gpvt = gas.GasPVT(sg=0.65, co2=0.1)
    result = gpvt.bg(2000, 180)
    assert abs(result - 0.008164459661048012) / 0.008164459661048012 < RTOL

def test_doc_nodal_oil_pvt_rs():
    """nodal.rst: OilPVT rs"""
    opvt = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    result = opvt.rs(2000, 180)
    assert abs(result - 403.58333168879415) / 403.58333168879415 < RTOL

def test_doc_nodal_oil_pvt_bo():
    """nodal.rst: OilPVT bo"""
    opvt = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    result = opvt.bo(2000, 180)
    assert abs(result - 1.22370082673546) / 1.22370082673546 < RTOL

def test_doc_nodal_oil_pvt_density():
    """nodal.rst: OilPVT density"""
    opvt = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    result = opvt.density(2000, 180)
    assert abs(result - 46.23700811760461) / 46.23700811760461 < RTOL

def test_doc_nodal_oil_pvt_viscosity():
    """nodal.rst: OilPVT viscosity"""
    opvt = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    result = opvt.viscosity(2000, 180)
    assert abs(result - 0.7187504436478858) / 0.7187504436478858 < RTOL

def test_doc_nodal_gas_pvt_metric_z():
    """nodal.rst: GasPVT metric z-factor"""
    gpvt_m = gas.GasPVT(sg=0.65, co2=0.1, metric=True)
    result = gpvt_m.z(137.9, 82.2)
    assert abs(result - 0.9026433033953588) / 0.9026433033953588 < RTOL

def test_doc_nodal_gas_pvt_metric_density():
    """nodal.rst: GasPVT metric density"""
    gpvt_m = gas.GasPVT(sg=0.65, co2=0.1, metric=True)
    result = gpvt_m.density(137.9, 82.2)
    assert abs(result - 97.36871728783241) / 97.36871728783241 < RTOL

def test_doc_nodal_oil_pvt_metric_rs():
    """nodal.rst: OilPVT metric rs"""
    opvt_m = oil.OilPVT(api=35, sg_sp=0.65, pb=172.4, rsb=89, metric=True)
    result = opvt_m.rs(137.9, 82.2)
    assert abs(result - 71.82727018664512) / 71.82727018664512 < RTOL

def test_doc_nodal_oil_pvt_metric_density():
    """nodal.rst: OilPVT metric density"""
    opvt_m = oil.OilPVT(api=35, sg_sp=0.65, pb=172.4, rsb=89, metric=True)
    result = opvt_m.density(137.9, 82.2)
    assert abs(result - 740.7086089268661) / 740.7086089268661 < RTOL

def test_doc_nodal_completion_basic():
    """nodal.rst: Completion with no casing"""
    c = nodal.Completion(tid=2.441, length=10000, tht=100, bht=200)
    assert c.has_casing_section == False

def test_doc_nodal_completion_casing():
    """nodal.rst: Completion with casing section"""
    c2 = nodal.Completion(tid=2.441, length=9000, tht=100, bht=200, cid=6.184, mpd=10000)
    assert c2.has_casing_section == True
    assert c2.casing_length == 1000

def test_doc_nodal_reservoir():
    """nodal.rst: Reservoir constructor"""
    r = nodal.Reservoir(pr=3000, degf=200, k=10, h=50, re=1500, rw=0.35, S=2, D=0.001)
    assert r.pr == 3000

def test_doc_nodal_fbhp_gas():
    """nodal.rst: fbhp gas well HB"""
    c = nodal.Completion(tid=2.441, length=10000, tht=100, bht=200)
    result = nodal.fbhp(thp=500, completion=c, vlpmethod='HB', well_type='gas',
                        qg_mmscfd=5.0, gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    assert abs(result - 952.6868477414688) / 952.6868477414688 < RTOL

def test_doc_nodal_fbhp_gas_wg():
    """nodal.rst: fbhp gas well WG"""
    c = nodal.Completion(tid=2.441, length=10000, tht=100, bht=200)
    result = nodal.fbhp(thp=500, completion=c, vlpmethod='WG', well_type='gas',
                        qg_mmscfd=5.0, gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    assert abs(result - 1172.8626065704736) / 1172.8626065704736 < RTOL

def test_doc_nodal_fbhp_gas_gray():
    """nodal.rst: fbhp gas well GRAY"""
    c = nodal.Completion(tid=2.441, length=10000, tht=100, bht=200)
    result = nodal.fbhp(thp=500, completion=c, vlpmethod='GRAY', well_type='gas',
                        qg_mmscfd=5.0, gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    assert abs(result - 1940.034905804135) / 1940.034905804135 < RTOL

def test_doc_nodal_fbhp_gas_bb():
    """nodal.rst: fbhp gas well BB"""
    c = nodal.Completion(tid=2.441, length=10000, tht=100, bht=200)
    result = nodal.fbhp(thp=500, completion=c, vlpmethod='BB', well_type='gas',
                        qg_mmscfd=5.0, gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    assert abs(result - 1213.2399909265969) / 1213.2399909265969 < RTOL

def test_doc_nodal_fbhp_oil():
    """nodal.rst: fbhp oil well HB"""
    c = nodal.Completion(tid=2.441, length=8000, tht=100, bht=180)
    result = nodal.fbhp(thp=200, completion=c, vlpmethod='HB', well_type='oil',
                        qt_stbpd=2000, gor=800, wc=0.3, gsg=0.65,
                        pb=2500, rsb=500, sgsp=0.65, api=35)
    assert abs(result - 2256.2340921828286) / 2256.2340921828286 < RTOL

def test_doc_nodal_fbhp_oil_pvt():
    """nodal.rst: fbhp oil well with OilPVT"""
    c = nodal.Completion(tid=2.441, length=8000, tht=100, bht=180)
    opvt = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    result = nodal.fbhp(thp=200, completion=c, vlpmethod='HB', well_type='oil',
                        oil_pvt=opvt, qt_stbpd=2000, gor=800, wc=0.3, gsg=0.65)
    assert abs(result - 2258.2220198140935) / 2258.2220198140935 < RTOL

def test_doc_nodal_outflow_curve():
    """nodal.rst: outflow_curve gas"""
    c = nodal.Completion(tid=2.441, length=10000, tht=100, bht=200)
    result = nodal.outflow_curve(thp=500, completion=c, vlpmethod='HB',
                                  well_type='gas', rates=[2.0, 5.0, 10.0, 15.0, 20.0], gsg=0.65)
    assert result['rates'] == [2.0, 5.0, 10.0, 15.0, 20.0]
    expected_bhp = [676.8, 925.7, 1498.7, 2121.5, 2757.8]
    for actual, expected in zip(result['bhp'], expected_bhp):
        assert abs(round(actual, 1) - expected) < 0.2, f"BHP mismatch: {round(actual, 1)} vs {expected}"

def test_doc_nodal_ipr_curve():
    """nodal.rst: ipr_curve gas"""
    r = nodal.Reservoir(pr=3000, degf=200, k=10, h=50, re=1500, rw=0.35, S=2, D=0.001)
    ipr = nodal.ipr_curve(r, well_type='gas', gsg=0.65, n_points=5)
    expected_pwf = [14.7, 761.0, 1507.4, 2253.7, 3000.0]
    for actual, expected in zip(ipr['pwf'], expected_pwf):
        assert abs(round(actual, 1) - expected) < 0.2
    expected_rate = [13456.5, 12812.2, 10861.0, 7290.7, 0.0]
    for actual, expected in zip(ipr['rate'], expected_rate):
        assert abs(round(actual, 1) - expected) < 0.2

def test_doc_nodal_operating_point_gas():
    """nodal.rst: operating_point gas"""
    c = nodal.Completion(tid=2.441, length=10000, tht=100, bht=200)
    r = nodal.Reservoir(pr=3000, degf=200, k=10, h=50, re=1500, rw=0.35, S=2, D=0.001)
    result = nodal.operating_point(thp=500, completion=c, reservoir=r,
                                    vlpmethod='HB', well_type='gas', gsg=0.65)
    assert abs(round(result['rate'], 2) - 10.61) < 0.02
    assert abs(round(result['bhp'], 1) - 1573.7) < 0.5

def test_doc_nodal_operating_point_oil():
    """nodal.rst: operating_point oil"""
    c = nodal.Completion(tid=2.441, length=8000, tht=100, bht=180)
    r = nodal.Reservoir(pr=3000, degf=180, k=50, h=30, re=1000, rw=0.35)
    opvt = oil.OilPVT(api=35, sg_sp=0.65, pb=2500, rsb=500)
    result = nodal.operating_point(thp=200, completion=c, reservoir=r,
                                    vlpmethod='HB', well_type='oil',
                                    oil_pvt=opvt, gor=800, wc=0.3, gsg=0.65)
    assert abs(round(result['rate'], 1) - 1412.6) < 1.0
    assert abs(round(result['bhp'], 1) - 2193.0) < 1.0


def test_doc_nodal_wellsegment_vertical():
    """nodal.rst: WellSegment vertical tvd"""
    seg = nodal.WellSegment(md=10000, id=2.441, deviation=0)
    assert round(seg.tvd, 1) == 10000.0

def test_doc_nodal_wellsegment_deviated():
    """nodal.rst: WellSegment 45-degree tvd"""
    seg = nodal.WellSegment(md=10000, id=2.441, deviation=45)
    assert round(seg.tvd, 1) == 7071.1

def test_doc_nodal_completion_segments():
    """nodal.rst: Completion with segments and total_md/total_tvd"""
    segs = [nodal.WellSegment(md=5000, id=2.441, deviation=0),
            nodal.WellSegment(md=3000, id=2.441, deviation=45),
            nodal.WellSegment(md=2000, id=2.441, deviation=60)]
    c3 = nodal.Completion(segments=segs, tht=100, bht=250)
    assert c3.total_md == 10000
    assert round(c3.total_tvd, 1) == 8121.3

def test_doc_nodal_geometry_at_md_deviated():
    """nodal.rst: geometry_at_md on multi-segment deviated well"""
    segs = [nodal.WellSegment(md=5000, id=2.441, deviation=0),
            nodal.WellSegment(md=3000, id=2.441, deviation=45),
            nodal.WellSegment(md=2000, id=4.0, deviation=60)]
    c = nodal.Completion(segments=segs, tht=100, bht=250)
    g = c.geometry_at_md(6500)
    assert round(g['tvd'], 1) == 6060.7
    assert g['id'] == 2.441
    assert g['deviation'] == 45

def test_doc_nodal_geometry_at_md_casing():
    """nodal.rst: geometry_at_md on legacy completion with casing"""
    c2 = nodal.Completion(tid=2.441, length=9000, tht=100, bht=200, cid=6.184, mpd=10000)
    g2 = c2.geometry_at_md(9500)
    assert g2['id'] == 6.184
    assert g2['tvd'] == 9500.0

def test_doc_nodal_profile_deviated():
    """nodal.rst: profile on multi-segment deviated well"""
    segs = [nodal.WellSegment(md=5000, id=2.441, deviation=0),
            nodal.WellSegment(md=3000, id=2.441, deviation=45),
            nodal.WellSegment(md=2000, id=4.0, deviation=60)]
    c = nodal.Completion(segments=segs, tht=100, bht=250)
    df = c.profile()
    assert len(df) == 6
    assert list(df.columns) == ['MD', 'TVD', 'Deviation', 'ID', 'Roughness']

def test_doc_nodal_profile_casing():
    """nodal.rst: profile on legacy completion with casing"""
    c2 = nodal.Completion(tid=2.441, length=9000, tht=100, bht=200, cid=6.184, mpd=10000)
    df2 = c2.profile()
    assert len(df2) == 4
    assert df2['ID'].iloc[1] == 2.441
    assert df2['ID'].iloc[2] == 6.184

def test_doc_nodal_fbhp_deviated():
    """nodal.rst: fbhp deviated well using WellSegment"""
    segs = [nodal.WellSegment(md=5000, id=2.441, deviation=0),
            nodal.WellSegment(md=5000, id=2.441, deviation=45)]
    c_dev = nodal.Completion(segments=segs, tht=100, bht=200)
    result = nodal.fbhp(thp=500, completion=c_dev, vlpmethod='HB', well_type='gas',
                        qg_mmscfd=5.0, gsg=0.65, cgr=10, qw_bwpd=10, api=45, oil_vis=1.0)
    assert abs(result - 923.092017723091) < 0.01


# =============================================================================
# Main runner
# =============================================================================

if __name__ == '__main__':
    import traceback
    tests = [(k, v) for k, v in sorted(globals().items()) if k.startswith('test_') and callable(v)]
    passed = failed = 0
    errors = []
    for name, func in tests:
        try:
            func()
            passed += 1
            print(f"  PASS: {name}")
        except Exception as e:
            failed += 1
            errors.append((name, str(e)))
            print(f"  FAIL: {name}")
            print(f"        {e}")
            traceback.print_exc()

    print(f"\n{'='*60}")
    print(f"TOTAL: {passed} passed, {failed} failed out of {passed + failed}")
    if errors:
        print(f"\nFailed tests:")
        for name, msg in errors:
            print(f"  - {name}: {msg}")
    print("="*60)
