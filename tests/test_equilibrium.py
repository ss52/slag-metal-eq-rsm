"""Equilibrium constants and the two-DeltaG decomposition (PLAN.md test 4)."""

from __future__ import annotations

import math

import pytest

from slag_model.constants import R
from slag_model.data import DELTA_G_C_DISSOLUTION, DELTA_G_CO_GRAPHITE
from slag_model.equilibrium import k_CO, k_Cr, k_FeO

T = 1823.15


def test_k_FeO_range():
    k = k_FeO(T)
    assert 1.3e4 <= k <= 1.6e4, f"K_FeO = {k:.4e}, expected [1.3e4, 1.6e4]"


def test_k_CO_range():
    k = k_CO(T)
    assert 1.5e6 <= k <= 1.8e6, f"K_CO = {k:.4e}, expected [1.5e6, 1.8e6]"


@pytest.mark.parametrize("temperature", [1773.15, 1823.15, 1873.15])
def test_k_Cr_source_reaction(temperature):
    dg = 4.184 * (25690.0 - 13.36 * temperature) - (19246.0 - 46.86 * temperature)
    assert k_Cr(temperature) == pytest.approx(math.exp(-dg / (R * temperature)), rel=1e-12)


def test_k_CO_approx_value():
    assert k_CO(T) == pytest.approx(1.66e6, rel=0.02)


def test_k_FeO_approx_value():
    assert k_FeO(T) == pytest.approx(1.45e4, rel=0.02)


def test_k_CO_henrian_decomposition():
    """Verify K_CO matches the two-DeltaG decomposition (PLAN.md 6.4b).

    K_CO(Henrian) = exp(-(DG_graphite - DG_dissolution)/RT)
    """
    dg_graphite = DELTA_G_CO_GRAPHITE[0] + DELTA_G_CO_GRAPHITE[1] * T
    dg_diss = DELTA_G_C_DISSOLUTION[0] + DELTA_G_C_DISSOLUTION[1] * T
    k_decomp = math.exp(-(dg_graphite - dg_diss) / (R * T))
    # Should equal the direct K_CO = exp(-(-134300 - 45.40*T)/RT)
    assert k_decomp == pytest.approx(k_CO(T), rel=1e-10)
