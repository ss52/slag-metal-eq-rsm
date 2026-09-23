"""Paper-derived validation tests for the thermodynamic model.

Known defects are strict xfails: they document the required scientific behavior
without making the repository's ordinary test run red before the fixes land.
Remove each xfail marker together with the corresponding production-code fix.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from slag_model.constants import R
from slag_model.data import ALPHA, CATION_ORDER, CONVERSION_DEFAULTS, EIJ, SIO2_CONVERSIONS
from slag_model.equilibrium import k_Cr
from slag_model.io import load_input
from slag_model.redox import fe_ratio_ban_ya
from slag_model.rsm import convert_gammas, rsm_gamma_rtln
from slag_model.solver import compute_state, solve_redox_fixed_po2

REVIEW_DIR = Path(__file__).resolve().parent
REPO_ROOT = REVIEW_DIR.parents[1]
CASES = json.loads((REVIEW_DIR / "paper_cases.json").read_text(encoding="utf-8"))
REFERENCE_INPUT = REPO_ROOT / "examples" / "eaf_slag_01.json"
IDX = {name: index for index, name in enumerate(CATION_ORDER)}


@pytest.mark.parametrize(
    ("oxide", "expected"),
    CASES["banya_1993"]["conversion_A_plus_B_T_J_per_mol"].items(),
)
def test_banya_table_3_conversion_coefficients(oxide, expected):
    assert CONVERSION_DEFAULTS[oxide] == pytest.approx(expected, abs=0.0)


@pytest.mark.parametrize(
    ("oxide", "expected"),
    CASES["xiao_holappa_reuter_2002"]["conversion_A_plus_B_T_J_per_mol"].items(),
)
def test_xiao_2002_table_iv_conversion_coefficients(oxide, expected):
    actual = SIO2_CONVERSIONS["workbook"] if oxide == "SiO2" else CONVERSION_DEFAULTS[oxide]
    assert actual == pytest.approx(expected, abs=0.0)


@pytest.mark.parametrize(
    ("pair", "expected"),
    [
        item
        for item in CASES["xiao_holappa_reuter_2002"]["interaction_energy_J_per_mol"].items()
        if item[0] != "Ca2+|Si4+"
    ],
)
def test_xiao_2002_table_iii_chromium_interactions(pair, expected):
    left, right = pair.split("|")
    assert ALPHA[IDX[left], IDX[right]] == pytest.approx(expected, abs=0.0)
    assert ALPHA[IDX[right], IDX[left]] == pytest.approx(expected, abs=0.0)


@pytest.mark.xfail(
    strict=True,
    reason="Known source-data defect: Fe3+-Ca2+ is transcribed as -96810 instead of -95810 J/mol.",
)
def test_banya_table_2_fe3_ca_interaction_energy():
    expected = CASES["banya_1993"]["interaction_energy_J_per_mol"]["Fe3+|Ca2+"]
    assert ALPHA[IDX["Fe3+"], IDX["Ca2+"]] == pytest.approx(expected, abs=0.0)


@pytest.mark.xfail(
    strict=True,
    reason="Known defect: no selectable Xiao-2002 parameter set.",
)
def test_xiao_2002_table_iii_ca_si_interaction_energy():
    expected = CASES["xiao_holappa_reuter_2002"]["interaction_energy_J_per_mol"]["Ca2+|Si4+"]
    assert ALPHA[IDX["Ca2+"], IDX["Si4+"]] == pytest.approx(expected, abs=0.0)


def test_banya_eq_24_is_coupled_with_regular_solution_gammas():
    config = load_input(REFERENCE_INPUT)
    state, targets = compute_state(config, r_Fe=0.10, r_Cr=1.0, C_wtpc=0.02)
    gamma_fe2_rs = math.exp(state.rtln_gamma_rs["Fe2+"] / (R * config.temperature_K))
    gamma_fe3_rs = math.exp(state.rtln_gamma_rs["Fe3+"] / (R * config.temperature_K))
    expected = fe_ratio_ban_ya(
        config.temperature_K,
        state.P_O2_atm,
        gamma_fe2_rs,
        gamma_fe3_rs,
    )
    assert targets["r_Fe"] == pytest.approx(expected, rel=1e-12)


def test_fixed_po2_fe_redox_uses_two_rs_gammas():
    cfg = load_input(REFERENCE_INPUT)
    state = solve_redox_fixed_po2(cfg, P_O2=1e-9, C=0.018)
    rt = R * cfg.temperature_K
    expected = fe_ratio_ban_ya(
        cfg.temperature_K,
        1e-9,
        math.exp(state.rtln_gamma_rs["Fe2+"] / rt),
        math.exp(state.rtln_gamma_rs["Fe3+"] / rt),
    )
    assert state.r_Fe == pytest.approx(expected, rel=1e-8)


def test_xiao_chromium_equilibrium_is_on_library_metal_standard_state():
    temperature = 1823.15
    source = CASES["xiao_holappa_1993"]
    a_cal, b_cal = source["chromium_redox"]["delta_g_A_plus_B_T_cal_per_mol"]
    a_dis, b_dis = source["chromium_dissolution_for_library_standard"][
        "delta_g_A_plus_B_T_J_per_mol"
    ]
    delta_g_oxide_j = 4.184 * (a_cal + b_cal * temperature)
    delta_g_dissolution_j = a_dis + b_dis * temperature
    expected = math.exp(-(delta_g_oxide_j - delta_g_dissolution_j) / (R * temperature))
    assert k_Cr(temperature) == pytest.approx(expected, rel=1e-12)


@pytest.mark.xfail(
    strict=True,
    reason="Known defect: P2O5 is treated as one PO2.5 unit.",
)
def test_banya_p2o5_formula_unit_activity_relation():
    temperature = 1873.0
    x_po2_5 = 0.02
    rtln_gamma_po2_5_rs = 4000.0
    gamma, delta_g = convert_gammas({"P5+": rtln_gamma_po2_5_rs}, temperature)
    library_activity = gamma["P2O5"] * x_po2_5
    pseudo_activity = math.exp(rtln_gamma_po2_5_rs / (R * temperature)) * x_po2_5
    expected = pseudo_activity**2 * math.exp(delta_g["P2O5"] / (R * temperature))
    assert library_activity == pytest.approx(expected, rel=1e-12, abs=0.0)


@pytest.mark.xfail(
    strict=True,
    reason="Known defect: e_Mn^C is transcribed as -0.012.",
)
def test_sigworth_elliott_manganese_carbon_interaction_parameter():
    assert EIJ[IDX_METAL("Mn"), IDX_METAL("C")] == pytest.approx(-0.07, abs=0.0)


def IDX_METAL(species: str) -> int:
    """Return the stable index used by the library's four-species metal matrix."""
    from slag_model.data import METAL_SPECIES

    return METAL_SPECIES.index(species)


@pytest.mark.parametrize(
    "case",
    CASES["xiao_holappa_reuter_2002"]["experimental_cases_1873_K"],
    ids=lambda case: case["id"],
)
def test_xiao_2002_table_i_experimental_activity_envelope(case):
    """Check the paper's stated approximate agreement, not exact curve-fitting.

    The factor-of-2.1 envelope is intentionally broad and is declared in the
    fixture. Exact tests above protect transcription and equation identities.
    """
    mol_percent = case["mol_percent"]
    total = sum(mol_percent.values())
    x = {
        "Ca2+": mol_percent["CaO"] / total,
        "Si4+": mol_percent["SiO2"] / total,
        "Cr2+": mol_percent["CrO"] / total,
        "Cr3+": mol_percent["CrO1.5"] / total,
    }
    rtln = rsm_gamma_rtln(x, "full")
    gamma, _ = convert_gammas(rtln, 1873.0, sio2_conversion="workbook")
    calculated = {
        "CrO": gamma["CrO"] * x["Cr2+"],
        "CrO1.5": gamma["CrO1.5"] * x["Cr3+"],
    }
    factor = 2.1
    for oxide, measured in case["measured_activity"].items():
        ratio = calculated[oxide] / measured
        assert 1.0 / factor <= ratio <= factor, (
            f"{case['id']} {oxide}: calculated={calculated[oxide]:.6g}, "
            f"measured={measured:.6g}, ratio={ratio:.3f}"
        )
