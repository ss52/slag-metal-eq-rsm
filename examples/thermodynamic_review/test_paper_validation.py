"""Paper-derived source checks and isolated activity comparisons.

Known defects are strict xfails: they document the required scientific behavior
without making the repository's ordinary test run red before the fixes land.
Remove each xfail marker together with the corresponding production-code fix.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from slag_model.constants import MOLAR_MASS_OXIDE, R
from slag_model.data import ALPHA, CATION_ORDER, CONVERSION_DEFAULTS, EIJ, SIO2_CONVERSIONS
from slag_model.equilibrium import k_Cr
from slag_model.io import load_input
from slag_model.redox import fe_ratio_ban_ya
from slag_model.rsm import build_slag_activities, rsm_gamma_rtln
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


def test_chromium_thermodynamic_coefficients_have_separate_verified_sources():
    sources = {source["id"]: source for source in CASES["sources"]}
    assert {"xiao_holappa_1995", "xiao_kou_fang_2018"} <= sources.keys()

    oxide_source = sources["xiao_holappa_1995"]
    dissolution_source = sources["xiao_kou_fang_2018"]
    sigworth_source = sources["sigworth_elliott_1974"]
    assert "1995" in oxide_source["citation"]
    assert "2018" in dissolution_source["citation"]
    assert any("p. 322" in use and "Eqs. 9" in use for use in oxide_source["used_for"])
    assert any("p. 419" in use and "Table 1" in use for use in dissolution_source["used_for"])
    assert any("Wagner" in use for use in sigworth_source["used_for"])
    assert not any("dissolv" in use.lower() for use in sigworth_source["used_for"])

    oxide = CASES["xiao_holappa_1995"]["chromium_redox"]
    dissolution = CASES["xiao_kou_fang_2018"]["chromium_dissolution_for_library_standard"]
    assert oxide["delta_g_A_plus_B_T_cal_per_mol"] == [25690.0, -13.36]
    assert dissolution["delta_g_A_plus_B_T_J_per_mol"] == [19246.0, -46.86]


def test_xiao_chromium_equilibrium_is_on_library_metal_standard_state():
    temperature = 1823.15
    oxide_source = CASES["xiao_holappa_1995"]
    dissolution_source = CASES["xiao_kou_fang_2018"]
    a_cal, b_cal = oxide_source["chromium_redox"]["delta_g_A_plus_B_T_cal_per_mol"]
    a_dis, b_dis = dissolution_source["chromium_dissolution_for_library_standard"][
        "delta_g_A_plus_B_T_J_per_mol"
    ]
    delta_g_oxide_j = 4.184 * (a_cal + b_cal * temperature)
    delta_g_dissolution_j = a_dis + b_dis * temperature
    expected = math.exp(-(delta_g_oxide_j - delta_g_dissolution_j) / (R * temperature))
    assert k_Cr(temperature) == pytest.approx(expected, rel=1e-12)


def test_banya_p2o5_formula_unit_activity_relation():
    temperature = 1873.0
    x_po2_5 = 0.02
    rtln_gamma_po2_5_rs = 4000.0
    activities = build_slag_activities(
        {"P5+": rtln_gamma_po2_5_rs},
        {"P5+": x_po2_5},
        temperature,
    )
    library_activity = activities.a_conventional_by_species["P2O5"]
    pseudo_activity = math.exp(rtln_gamma_po2_5_rs / (R * temperature)) * x_po2_5
    a, b = CASES["banya_1993"]["conversion_A_plus_B_T_J_per_mol"]["P2O5"]
    expected = pseudo_activity**2 * math.exp((a + b * temperature) / (R * temperature))
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


def _banya_1985_prediction(case):
    wt = case["reported_wt_percent"]
    n = {
        "Fe2+": wt["FeO"] / MOLAR_MASS_OXIDE["FeO"],
        "Fe3+": 2.0 * wt["Fe2O3"] / MOLAR_MASS_OXIDE["Fe2O3"],
        "Si4+": wt["SiO2"] / MOLAR_MASS_OXIDE["SiO2"],
        "Mn2+": wt["MnO"] / MOLAR_MASS_OXIDE["MnO"],
    }
    total = sum(n.values())
    x = {cation: amount / total for cation, amount in n.items()}
    full_x = {cation: x.get(cation, 0.0) for cation in CATION_ORDER}
    rtln = rsm_gamma_rtln(x, "full")
    activities = build_slag_activities(rtln, full_x, 1723.15, sio2_conversion="banya")
    return activities.a_conventional_by_species["FeO"], x["Fe2+"] + x["Fe3+"], n["Fe3+"] / n["Fe2+"]


def _xiao_2002_prediction(case):
    mol_percent = case["mol_percent"]
    total = sum(mol_percent.values())
    x = {
        "Ca2+": mol_percent["CaO"] / total,
        "Si4+": mol_percent["SiO2"] / total,
        "Cr2+": mol_percent["CrO"] / total,
        "Cr3+": mol_percent["CrO1.5"] / total,
    }
    rtln = rsm_gamma_rtln(x, "full")
    full_x = {cation: x.get(cation, 0.0) for cation in CATION_ORDER}
    activities = build_slag_activities(rtln, full_x, 1873.0, sio2_conversion="workbook")
    return {oxide: activities.a_conventional_by_species[oxide] for oxide in ("CrO", "CrO1.5")}, x


@pytest.mark.parametrize(
    "case",
    CASES["banya_1985"]["experimental_cases_1723_15_K"],
    ids=lambda case: case["id"],
)
def test_banya_1985_table_i_oxide_split_matches_printed_fe_ratio(case):
    _, _, calculated_ratio = _banya_1985_prediction(case)
    assert calculated_ratio == pytest.approx(case["reported_Fe3_to_Fe2"], abs=0.001)


def test_banya_1985_table_i_cases_are_available():
    """Keep a fixed source-based Ban-ya 1985 sample for activity comparison."""
    cases = CASES["banya_1985"]["experimental_cases_1723_15_K"]
    assert [case["id"] for case in cases] == ["101", "301", "501", "701", "901"]


@pytest.mark.parametrize(
    "case",
    [
        case
        for case in CASES["banya_1985"]["experimental_cases_1723_15_K"]
        if case["validation_role"] == "project_screen"
    ],
    ids=lambda case: case["id"],
)
def test_banya_1985_screened_rows_match_reported_fe_to_activity(case):
    calculated, _, _ = _banya_1985_prediction(case)
    measured = case["measured_activity_Fe_tO"]
    tolerance = CASES["validation_policy"]["banya_activity_relative_error_screen"]
    relative_error = calculated / measured - 1.0
    assert abs(relative_error) <= tolerance, (
        f"{case['id']}: calculated={calculated:.6f}, measured={measured:.3f}, "
        f"relative error={relative_error:+.2%}"
    )


@pytest.mark.parametrize(
    "case",
    CASES["banya_1985"]["experimental_cases_1723_15_K"],
    ids=lambda case: case["id"],
)
def test_banya_1985_calculated_activity_matches_diagnostic_snapshot(case):
    calculated, _, _ = _banya_1985_prediction(case)
    expected = CASES["model_diagnostic_snapshots"]["banya_1985_1723_15_K"]["activity_Fe_tO"][
        case["id"]
    ]
    assert calculated == pytest.approx(expected, abs=5e-7)


def test_banya_1985_high_iron_row_is_stress_diagnostic_only():
    case = CASES["banya_1985"]["experimental_cases_1723_15_K"][-1]
    _, x_fe_total, _ = _banya_1985_prediction(case)
    assert case["id"] == "901"
    assert case["validation_role"] == "stress_diagnostic"
    assert x_fe_total == pytest.approx(0.84236, abs=5e-6)
    assert x_fe_total > 0.7


def _xiao_2002_pure_solid_cr_k():
    temperature = 1873.0
    a_cal, b_cal = CASES["xiao_holappa_1995"]["chromium_redox"]["delta_g_A_plus_B_T_cal_per_mol"]
    delta_g_j_per_mol = 4.184 * (a_cal + b_cal * temperature)
    return math.exp(-delta_g_j_per_mol / (R * temperature))


def _xiao_2002_measured_redox_q_over_k(case):
    activity = case["measured_activity"]
    q = activity["CrO"] ** 3 / activity["CrO1.5"] ** 2
    return q / _xiao_2002_pure_solid_cr_k()


@pytest.mark.parametrize(
    "case",
    CASES["xiao_holappa_reuter_2002"]["experimental_cases_1873_K"],
    ids=lambda case: case["id"],
)
def test_xiao_2002_measured_activity_identity_matches_reported_gamma(case):
    _, x = _xiao_2002_prediction(case)
    for oxide, cation in (("CrO", "Cr2+"), ("CrO1.5", "Cr3+")):
        reconstructed = case["measured_gamma"][oxide] * x[cation]
        measured = case["measured_activity"][oxide]
        assert reconstructed == pytest.approx(
            measured,
            abs=CASES["validation_policy"]["xiao_source_activity_identity_abs_tolerance"],
        )


@pytest.mark.parametrize(
    "case",
    CASES["xiao_holappa_reuter_2002"]["experimental_cases_1873_K"],
    ids=lambda case: case["id"],
)
def test_xiao_2002_calculated_activities_match_diagnostic_snapshot(case):
    calculated, _ = _xiao_2002_prediction(case)
    expected = CASES["model_diagnostic_snapshots"]["xiao_2002_1873_K"]["activities"][case["id"]]
    for oxide in ("CrO", "CrO1.5"):
        assert calculated[oxide] == pytest.approx(expected[oxide], abs=5e-7)


@pytest.mark.parametrize(
    "case",
    CASES["xiao_holappa_reuter_2002"]["experimental_cases_1873_K"],
    ids=lambda case: case["id"],
)
def test_xiao_2002_calculated_redox_q_over_k_matches_snapshot(case):
    calculated, _ = _xiao_2002_prediction(case)
    q = calculated["CrO"] ** 3 / calculated["CrO1.5"] ** 2
    q_over_k = q / _xiao_2002_pure_solid_cr_k()
    expected = CASES["model_diagnostic_snapshots"]["xiao_2002_1873_K"]["q_over_k_pure_solid_cr"][
        case["id"]
    ]
    assert q_over_k == pytest.approx(expected, abs=0.001)


@pytest.mark.parametrize(
    "case",
    CASES["xiao_holappa_reuter_2002"]["experimental_cases_1873_K"],
    ids=lambda case: case["id"],
)
def test_xiao_2002_published_activity_redox_diagnostic(case):
    q_over_k = _xiao_2002_measured_redox_q_over_k(case)
    expected = CASES["xiao_2002_measured_activity_q_over_k_snapshot"][case["id"]]
    assert q_over_k == pytest.approx(expected, abs=0.001)
    if case["id"] == "CSC7":
        assert q_over_k > 1.7
    else:
        assert 0.97 <= q_over_k <= 1.03
