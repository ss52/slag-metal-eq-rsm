"""Equilibrium constants and the slag-Fe-O oxygen-potential closure (7.7)."""

from __future__ import annotations

import math

from .constants import R
from .data import (
    DELTA_G_CO,
    DELTA_G_CR_DISSOLUTION_J,
    DELTA_G_CR_SOLID_CAL,
    DELTA_G_FEO,
)


def _delta_g(coeffs: tuple[float, float], T: float) -> float:
    return coeffs[0] + coeffs[1] * T


def k_FeO(T: float) -> float:
    """Fe(l) + 1/2 O2 = FeO(l):  K_FeO = a_FeO / (a_Fe * sqrt(P_O2))."""
    return math.exp(-_delta_g(DELTA_G_FEO, T) / (R * T))


def k_CO(T: float) -> float:
    """[C](1wt%) + 1/2 O2 = CO(g):  K_CO = P_CO / (a_C * sqrt(P_O2)).

    Coefficients already on the Henrian 1 wt% standard state (PLAN.md 6.4b).
    """
    return math.exp(-_delta_g(DELTA_G_CO, T) / (R * T))


def k_Cr(T: float) -> float:
    """2 CrO1.5 + [Cr]1wt% = 3 CrO on the dissolved-metal standard state."""
    dg_oxide_j = 4.184 * _delta_g(DELTA_G_CR_SOLID_CAL, T)
    dg_dissolution_j = _delta_g(DELTA_G_CR_DISSOLUTION_J, T)
    return math.exp(-(dg_oxide_j - dg_dissolution_j) / (R * T))


def p_o2_from_slag(a_FeO: float, T: float, a_Fe: float) -> float:
    """Oxygen potential from the Fe-O equilibrium buffer in the slag (7.7).

    P_O2 = (a_FeO / (K_FeO * a_Fe))^2
    """
    return (a_FeO / (k_FeO(T) * a_Fe)) ** 2
