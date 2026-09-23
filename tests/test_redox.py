"""Fixed-P_O2 redox regression and cubic solver (PLAN.md test 4, 5)."""

from __future__ import annotations

import math
import pathlib

import pytest

from slag_model.constants import R
from slag_model.io import load_input
from slag_model.redox import cr_ratio, fe_ratio_ban_ya
from slag_model.solver import solve_redox_fixed_po2

ROOT = pathlib.Path(__file__).resolve().parents[1]
SLAG_PATH = ROOT / "examples" / "eaf_slag_01.json"


def test_fe_ratio_ban_ya():
    """Eq. 24 at fixed conditions."""
    T = 1823.15
    P_O2 = 1e-9
    g_FeO = 1.754
    g_FeO1_5 = 0.5512
    r = fe_ratio_ban_ya(T, P_O2, g_FeO, g_FeO1_5)
    assert r == pytest.approx(0.131, abs=0.002)


def test_cr_ratio_cubic():
    """Exact cubic r^3/(1+r) = rhs at a known value."""
    r_target = 1.25
    rhs = r_target**3 / (1.0 + r_target)
    assert cr_ratio(rhs) == pytest.approx(r_target, abs=1e-8)


def test_cr_ratio_zero_rhs():
    assert cr_ratio(0.0) == 0.0


def test_cr_ratio_monotonic():
    assert cr_ratio(0.1) < cr_ratio(1.0) < cr_ratio(10.0)


def test_fe_split_regression():
    """Fixed-P_O2 Fe ratio satisfies Eq. 24 using both RS gammas."""
    cfg = load_input(SLAG_PATH)
    cfg.options.cross_terms = "major5"
    cfg.options.ti_handling = "as_excel"
    cfg.options.sio2_conversion = "workbook"
    state = solve_redox_fixed_po2(cfg, P_O2=1e-9, C=0.018)
    rt = R * cfg.temperature_K
    expected = fe_ratio_ban_ya(
        cfg.temperature_K,
        1e-9,
        math.exp(state.rtln_gamma_rs["Fe2+"] / rt),
        math.exp(state.rtln_gamma_rs["Fe3+"] / rt),
    )
    assert state.r_Fe == pytest.approx(expected, rel=1e-8)


def test_cr_split_regression():
    """Fixed-P_O2 chromium ratio satisfies the Henrian equilibrium identity."""
    cfg = load_input(SLAG_PATH)
    cfg.options.cross_terms = "major5"
    cfg.options.ti_handling = "as_excel"
    cfg.options.sio2_conversion = "workbook"
    state = solve_redox_fixed_po2(cfg, P_O2=1e-9, C=0.018)
    gamma = state.activities.gamma_conventional_by_species
    rhs = (
        state.k_Cr
        * state.a_Cr
        * gamma["CrO1.5"] ** 2
        * state.N
        / (gamma["CrO"] ** 3 * state.split.n_totals["Cr"])
    )
    assert state.r_Cr**3 / (1.0 + state.r_Cr) == pytest.approx(rhs, rel=1e-8)
