#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""
    pyResToolbox - A collection of Reservoir Engineering Utilities
              Copyright (C) 2022, Mark Burgoyne

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    The GNU General Public License can be found in the LICENSE directory,
    and at  <https://www.gnu.org/licenses/>.

          Contact author at mark.w.burgoyne@gmail.com

Method recommendation engine for selecting appropriate correlations.

Functions
---------
recommend_gas_methods   Recommend Z-factor and critical property methods
recommend_oil_methods   Recommend oil PVT correlation methods
recommend_vlp_method    Recommend VLP multiphase flow correlation
recommend_methods       Master function combining all recommendations

Classes
-------
MethodRecommendation    Single method recommendation with rationale
"""

__all__ = [
    'recommend_gas_methods', 'recommend_oil_methods',
    'recommend_vlp_method', 'recommend_methods',
    'MethodRecommendation',
]

from dataclasses import dataclass, field
from typing import List, Optional, Dict


@dataclass
class MethodRecommendation:
    """Single method recommendation.

    Attributes
    ----------
    category : str
        Category of the recommendation (e.g. 'zmethod', 'vlp_method').
    recommended : str
        Recommended method name.
    rationale : str
        Explanation for the recommendation.
    alternatives : list of str
        Other viable methods.
    mandatory : bool
        True if the recommended method is the only valid choice.
    """
    category: str
    recommended: str
    rationale: str
    alternatives: List[str] = field(default_factory=list)
    mandatory: bool = False


def recommend_gas_methods(sg: float = 0.65, co2: float = 0, h2s: float = 0,
                          n2: float = 0, h2: float = 0) -> Dict[str, MethodRecommendation]:
    """Recommend Z-factor and critical property methods for a gas composition.

    Parameters
    ----------
    sg : float
        Gas specific gravity (default 0.65). Currently unused by the decision
        logic — accepted for API consistency and reserved for future rules
        (e.g. heavy gas / condensate-laden streams). Breaking removal deferred.
    co2 : float
        CO2 mole fraction (default 0).
    h2s : float
        H2S mole fraction (default 0).
    n2 : float
        N2 mole fraction (default 0).
    h2 : float
        H2 mole fraction (default 0).

    Returns
    -------
    dict of str to MethodRecommendation
        Keys: 'zmethod', 'cmethod'.
    """
    _ = sg  # reserved, see docstring
    inerts = co2 + h2s + n2 + h2
    recs = {}

    if h2 > 0:
        recs['zmethod'] = MethodRecommendation(
            category='zmethod',
            recommended='BNS',
            rationale=f'H2 present ({h2:.1%}). BNS is the only method with H2 support.',
            alternatives=[],
            mandatory=True,
        )
        recs['cmethod'] = MethodRecommendation(
            category='cmethod',
            recommended='BNS',
            rationale='BNS critical properties required when using BNS Z-factor method.',
            alternatives=[],
            mandatory=True,
        )
    elif inerts > 0.55:
        recs['zmethod'] = MethodRecommendation(
            category='zmethod',
            recommended='BNS',
            rationale=f'High inert content ({inerts:.1%} > 55%). BNS 5-component PR-EOS handles extreme compositions.',
            alternatives=['DAK'],
        )
        recs['cmethod'] = MethodRecommendation(
            category='cmethod',
            recommended='BNS',
            rationale='BNS critical properties recommended for consistency with BNS Z-factor.',
            alternatives=['PMC'],
        )
    elif co2 > 0.10 or h2s > 0.10:
        recs['zmethod'] = MethodRecommendation(
            category='zmethod',
            recommended='DAK',
            rationale=f'Elevated CO2 ({co2:.1%}) or H2S ({h2s:.1%}). DAK with PMC handles moderate impurities well.',
            alternatives=['BNS', 'HY'],
        )
        recs['cmethod'] = MethodRecommendation(
            category='cmethod',
            recommended='PMC',
            rationale='Piper-McCain-Corredor includes CO2/H2S corrections.',
            alternatives=['SUT', 'BNS'],
        )
    else:
        recs['zmethod'] = MethodRecommendation(
            category='zmethod',
            recommended='DAK',
            rationale='Clean gas. DAK (Dranchuk & Abou-Kassem) is the standard default.',
            alternatives=['HY', 'WYW', 'BNS'],
        )
        recs['cmethod'] = MethodRecommendation(
            category='cmethod',
            recommended='PMC',
            rationale='Piper-McCain-Corredor is the standard default for critical properties.',
            alternatives=['SUT', 'BNS'],
        )

    return recs


def recommend_oil_methods(api: float = 35.0) -> Dict[str, MethodRecommendation]:
    """Recommend oil PVT correlation methods.

    Parameters
    ----------
    api : float
        Oil API gravity (default 35.0).

    Returns
    -------
    dict of str to MethodRecommendation
        Keys: 'pbmethod', 'rsmethod', 'bomethod'.
    """
    recs = {}

    if api < 10:
        rationale_pb = f'Heavy oil (API={api:.1f}). Velarde-Blasingame-McCain is most robust for heavy oils.'
    elif api > 50:
        rationale_pb = f'Light/condensate (API={api:.1f}). Velarde-Blasingame-McCain provides good coverage.'
    else:
        rationale_pb = f'Medium oil (API={api:.1f}). Velarde-Blasingame-McCain is the recommended default.'

    recs['pbmethod'] = MethodRecommendation(
        category='pbmethod',
        recommended='VELAR',
        rationale=rationale_pb,
        alternatives=['VALMC', 'STAN'],
    )
    recs['rsmethod'] = MethodRecommendation(
        category='rsmethod',
        recommended='VELAR',
        rationale='Velarde Rs is the recommended default, consistent with VELAR bubble point.',
        alternatives=['STAN', 'VALMC'],
    )
    recs['bomethod'] = MethodRecommendation(
        category='bomethod',
        recommended='MCAIN',
        rationale='McCain density-based FVF is the recommended default.',
        alternatives=['STAN'],
    )

    return recs


def recommend_vlp_method(deviation: float = 0,
                         well_type: str = 'gas') -> Dict[str, MethodRecommendation]:
    """Recommend VLP multiphase flow correlation.

    Parameters
    ----------
    deviation : float
        Maximum wellbore deviation from vertical in degrees (default 0).
    well_type : str
        'gas' or 'oil' (default 'gas'). Currently unused by the decision logic
        — accepted for API consistency and reserved for future fluid-specific
        recommendations. Breaking removal deferred.

    Returns
    -------
    dict of str to MethodRecommendation
        Key: 'vlp_method'.
    """
    _ = well_type  # reserved, see docstring
    recs = {}

    if deviation <= 30:
        recs['vlp_method'] = MethodRecommendation(
            category='vlp_method',
            recommended='BB',
            rationale=f'Vertical/near-vertical well (dev={deviation:.0f} deg). All four methods are suitable. '
                       'Beggs & Brill is the most widely used general-purpose correlation.',
            alternatives=['HB', 'WG', 'GRAY'],
        )
    else:
        recs['vlp_method'] = MethodRecommendation(
            category='vlp_method',
            recommended='BB',
            rationale=f'Deviated well (dev={deviation:.0f} deg > 30 deg). '
                       'HB and Gray were developed for vertical flow data and become unreliable beyond ~30 deg. '
                       'BB and WG are suitable for all inclinations.',
            alternatives=['WG'],
        )

    return recs


def recommend_methods(sg: float = 0.65, co2: float = 0, h2s: float = 0,
                      n2: float = 0, h2: float = 0,
                      api: Optional[float] = None,
                      deviation: float = 0,
                      well_type: str = 'gas') -> Dict[str, MethodRecommendation]:
    """Master recommendation function combining gas, oil, and VLP recommendations.

    Parameters
    ----------
    sg : float
        Gas specific gravity (default 0.65).
    co2 : float
        CO2 mole fraction (default 0).
    h2s : float
        H2S mole fraction (default 0).
    n2 : float
        N2 mole fraction (default 0).
    h2 : float
        H2 mole fraction (default 0).
    api : float, optional
        Oil API gravity. If provided, oil method recommendations are included.
    deviation : float
        Maximum wellbore deviation from vertical in degrees (default 0).
    well_type : str
        'gas' or 'oil' (default 'gas').

    Returns
    -------
    dict of str to MethodRecommendation
    """
    recs = {}
    recs.update(recommend_gas_methods(sg=sg, co2=co2, h2s=h2s, n2=n2, h2=h2))

    if api is not None:
        recs.update(recommend_oil_methods(api=api))

    recs.update(recommend_vlp_method(deviation=deviation, well_type=well_type))

    return recs
