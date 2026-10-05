"""WIPF regression and carbon fixed-point (PLAN.md test 3)."""

from __future__ import annotations

import math

import pytest

from slag_model.equilibrium import k_CO
from slag_model.metal import (
    balanced_carbon,
    carbon_fixed_point,
    dissolved_oxygen,
    wipf_factors,
)

T = 1823.15
METAL = {"Cr": 0.038, "Mn": 0.007, "P": 0.069}
C_VAL = 0.018


def test_wipf_log_fC():
    f = wipf_factors(METAL, C_VAL)
    assert math.log10(f["C"]) == pytest.approx(0.0051, abs=5e-4)


def test_wipf_log_fCr():
    f = wipf_factors(METAL, C_VAL)
    assert math.log10(f["Cr"]) == pytest.approx(-0.0058, abs=5e-4)


def test_wipf_log_fMn():
    f = wipf_factors(METAL, C_VAL)
    assert math.log10(f["Mn"]) == pytest.approx(-0.0015015, rel=0.0, abs=1e-12)


def test_wipf_log_fP():
    f = wipf_factors(METAL, C_VAL)
    assert math.log10(f["P"]) == pytest.approx(0.0012, abs=5e-4)


def test_wipf_manganese_carbon_factor_uses_compiled_parameter():
    f = wipf_factors({}, C_wtpc=0.5)
    assert f["Mn"] == pytest.approx(10.0 ** (-0.07 * 0.5), rel=1e-13)
    assert f["C"] == pytest.approx(10.0 ** (0.14 * 0.5), rel=1e-13)


def test_dissolved_oxygen_accounts_for_iron_activity():
    expected = (0.3 / 0.8) * 10.0 ** (-6320.0 / 1823.15 + 2.734)
    assert dissolved_oxygen(0.3, 1823.15, a_Fe=0.8) == pytest.approx(expected, rel=1e-14)
    assert dissolved_oxygen(0.3, 1823.15) == pytest.approx(
        0.3 * 10.0 ** (-6320.0 / 1823.15 + 2.734), rel=1e-14
    )


def test_balanced_carbon_fixed_point():
    """C-O balance converges and is self-consistent."""
    P_O2 = 3e-10
    K_CO = k_CO(T)
    C_converged, f_C = carbon_fixed_point(1.0, K_CO, P_O2, METAL, passes=3)
    assert 0.02 < C_converged < 0.06
    # self-consistency: recompute C from the final f_C
    C_check = balanced_carbon(1.0, K_CO, P_O2, f_C)
    assert C_check == pytest.approx(C_converged, rel=1e-6)
