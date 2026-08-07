"""Balanced solver integration tests (PLAN.md test 5, 6)."""

from __future__ import annotations

import pathlib

import pytest

from slag_model.io import load_input

ROOT = pathlib.Path(__file__).resolve().parents[1]
SLAG_PATH = ROOT / "examples" / "eaf_slag_01.json"


@pytest.fixture(scope="module")
def default_result():
    cfg = load_input(SLAG_PATH)
    from slag_model.solver import solve

    return solve(cfg)


def test_converges(default_result):
    assert default_result.iterations < 200


def test_sum_X(default_result):
    assert sum(default_result.X.values()) == pytest.approx(1.0, abs=1e-12)


def test_mass_balance_Fe(default_result):
    n = default_result.split.n_cations
    tot = default_result.split.n_totals["Fe"]
    assert n["Fe2+"] + n["Fe3+"] == pytest.approx(tot, abs=1e-12)


def test_mass_balance_Cr(default_result):
    n = default_result.split.n_cations
    tot = default_result.split.n_totals["Cr"]
    assert n["Cr2+"] + n["Cr3+"] == pytest.approx(tot, abs=1e-12)


def test_P_O2_range(default_result):
    assert 2e-10 <= default_result.P_O2_atm <= 6e-10, (
        f"P_O2 = {default_result.P_O2_atm:.4e}, expected [2e-10, 6e-10]"
    )


def test_r_Fe_range(default_result):
    assert 0.05 <= default_result.r_Fe <= 0.25, (
        f"r_Fe = {default_result.r_Fe:.6f}, expected [0.05, 0.25]"
    )


def test_balanced_C_range(default_result):
    assert 0.025 <= default_result.C_wtpc <= 0.045, (
        f"[%C] = {default_result.C_wtpc:.4f}, expected [0.025, 0.045]"
    )


def test_dissolved_O_range(default_result):
    assert 0.04 <= default_result.O_wtpc <= 0.06, (
        f"[%O] = {default_result.O_wtpc:.4f}, expected [0.04, 0.06]"
    )


def test_QK_FeO_unity(default_result):
    assert default_result.q_over_k["FeO"] == pytest.approx(1.0, abs=1e-6)


def test_QK_CO_unity(default_result):
    assert default_result.q_over_k["CO"] == pytest.approx(1.0, abs=1e-6)
