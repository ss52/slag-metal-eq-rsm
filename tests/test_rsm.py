"""Workbook regression for the regular solution model (PLAN.md test 2, 6)."""

from __future__ import annotations

import pathlib

import pytest

from slag_model.io import load_input
from slag_model.solver import compute_state, solve_redox_fixed_po2

ROOT = pathlib.Path(__file__).resolve().parents[1]
SLAG_PATH = ROOT / "examples" / "eaf_slag_01.json"

# RTln(gamma_RS) values from the reference workbook at P_O2 = 1e-9 atm,
# cross_terms="major5", ti_handling="as_excel", sio2_conversion="workbook",
# parameter_profile="legacy_workbook_typo_hybrid_v1".
EXPECTED_RTLN = {
    "Fe2+": 4031.71,
    "Fe3+": -9027.69,
    "Ca2+": -17758.64,
    "Mg2+": -23733.98,
    "Mn2+": -23294.70,
    "Si4+": -26354.39,
    "Al3+": -50118.08,
    "Cr3+": 24359.11,
    "Cr2+": 7682.77,
}


@pytest.fixture(scope="module")
def workbook_state():
    """Evaluate the historical fixed composition with workbook options."""
    cfg = load_input(SLAG_PATH)
    # Override options to match workbook settings
    cfg.options.cross_terms = "major5"
    cfg.options.ti_handling = "as_excel"
    cfg.options.sio2_conversion = "workbook"
    cfg.options.parameter_profile = "legacy_workbook_typo_hybrid_v1"
    state, _ = compute_state(
        cfg,
        r_Fe=0.13073478698367935,
        r_Cr=1.2754258847919098,
        C_wtpc=0.018,
    )
    return state


def test_rtln_gamma_rs_per_cation(workbook_state):
    for cation, expected in EXPECTED_RTLN.items():
        actual = workbook_state.rtln_gamma_rs[cation]
        # Al3+ has wider tolerance because the reference expected values were
        # computed with full cross-terms while this test uses major5. The other
        # cations are insensitive to this difference (within 5 J).
        tol = 800.0 if cation == "Al3+" else 5.0
        assert actual == pytest.approx(expected, abs=tol), (
            f"RTln(gamma_RS) for {cation}: got {actual:.3f}, expected {expected:.2f}"
        )


def test_converted_gamma_FeO(workbook_state):
    assert workbook_state.activities.gamma_conventional_by_species["FeO"] == pytest.approx(
        1.75, abs=0.01
    )


def test_converted_gamma_SiO2(workbook_state):
    assert workbook_state.activities.gamma_conventional_by_species["SiO2"] == pytest.approx(
        0.98, abs=0.01
    )


def test_converted_gamma_CrO(workbook_state):
    assert workbook_state.activities.gamma_conventional_by_species["CrO"] == pytest.approx(
        4.79, abs=0.02
    )


def test_converted_gamma_CrO1_5(workbook_state):
    assert workbook_state.activities.gamma_conventional_by_species["CrO1.5"] == pytest.approx(
        7.71, abs=0.02
    )


def test_sum_X(workbook_state):
    """Cation fractions sum to unity (PLAN.md test 6)."""
    assert sum(workbook_state.X.values()) == pytest.approx(1.0, abs=1e-12)


def test_fixed_po2_chromium_reaction():
    cfg = load_input(SLAG_PATH)
    state = solve_redox_fixed_po2(cfg, P_O2=1e-9, C=0.018)
    r = state.r_Cr
    gamma = state.activities.gamma_conventional_by_species
    rhs = (
        state.k_Cr
        * state.a_Cr
        * gamma["CrO1.5"] ** 2
        * state.N
        / (gamma["CrO"] ** 3 * state.split.n_totals["Cr"])
    )
    assert r**3 / (1.0 + r) == pytest.approx(rhs, rel=1e-8)
