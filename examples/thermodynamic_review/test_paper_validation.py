"""Source-exact paper data and explicitly in-sample activity comparisons."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from slag_model import data
from slag_model.constants import R
from slag_model.data import CATION_ORDER, CONVERSION_DEFAULTS, EIJ, SIO2_CONVERSIONS
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


def test_paper_pack_uses_only_permitted_oxide_sources_and_labels_external_sources():
    sources = {source["id"]: source for source in CASES["sources"]}
    assert set(sources) == {
        "banya_1993",
        "xiao_holappa_1995",
        "xiao_kou_fang_2018",
        "sigworth_elliott_1974",
    }
    assert sources["banya_1993"]["scope"] == "oxide_model"
    assert sources["xiao_holappa_1995"]["scope"] == "oxide_model_and_in_sample_comparison"
    assert sources["xiao_kou_fang_2018"]["scope"] == "metal_side_coefficient_only"
    assert sources["sigworth_elliott_1974"]["scope"] == "metal_side_coefficients_only"


@pytest.mark.parametrize(
    ("oxide", "expected"),
    CASES["banya_1993"]["conversion_A_plus_B_T_J_per_mol"].items(),
)
def test_banya_1993_table_3_conversion_coefficients(oxide, expected):
    assert CONVERSION_DEFAULTS[oxide] == tuple(expected)


def test_banya_1993_silica_conversion_is_labeled_cristobalite():
    source_values = CASES["banya_1993"]["silica_conversions_A_plus_B_T_J_per_mol"]
    assert source_values["SiO2(beta-cristobalite)=SiO2(RS)"] == [27030.0, -1.983]
    assert data.SIO2_CONVERSION_STANDARD_STATES["banya"] == (
        "SiO2(beta-cristobalite), Ban-ya (1993) Table 3"
    )


def test_banya_1993_table_1_interactions_are_transcribed_exactly():
    expected = {"Ca2+|Si4+": -133890.0, "Fe3+|Ca2+": -95810.0}
    assert CASES["banya_1993"]["interaction_energy_J_per_mol"] == expected


def test_banya_1993_table_1_fe3_ca_interaction_is_exact():
    assert data.PARAMETER_PROFILES["banya93_xiao95_v1"][IDX["Fe3+"], IDX["Ca2+"]] == -95810.0


def test_xiao_1995_table_2_chromium_interactions_match_model_matrix():
    interactions = CASES["xiao_holappa_1995"]["interaction_energy_J_per_mol"]
    assert set(interactions) == {
        "Ca2+|Si4+",
        "Ca2+|Cr2+",
        "Ca2+|Cr3+",
        "Si4+|Cr2+",
        "Si4+|Cr3+",
        "Cr2+|Cr3+",
        "Cr2+|Mg2+",
        "Cr3+|Mg2+",
        "Cr2+|Al3+",
        "Cr3+|Al3+",
    }
    for pair, expected in interactions.items():
        left, right = pair.split("|")
        if pair == "Ca2+|Si4+":
            continue
        assert data.PARAMETER_PROFILES["banya93_xiao95_v1"][IDX[left], IDX[right]] == expected
        assert data.PARAMETER_PROFILES["banya93_xiao95_v1"][IDX[right], IDX[left]] == expected


def test_xiao_1995_table_2_ca_si_interaction_is_the_selected_default():
    assert data.PARAMETER_PROFILES["banya93_xiao95_v1"][IDX["Ca2+"], IDX["Si4+"]] == -139100.0


def test_banya_ca_si_comparison_and_workbook_legacy_profiles_are_distinct():
    banya_comparison = data.PARAMETER_PROFILES["banya93_casi_xiao95cr_hybrid_v1"]
    workbook_legacy = data.PARAMETER_PROFILES["legacy_workbook_typo_hybrid_v1"]
    assert banya_comparison[IDX["Fe3+"], IDX["Ca2+"]] == -95810.0
    assert banya_comparison[IDX["Ca2+"], IDX["Si4+"]] == -133890.0
    assert workbook_legacy[IDX["Fe3+"], IDX["Ca2+"]] == -96810.0
    assert workbook_legacy[IDX["Ca2+"], IDX["Si4+"]] == -133890.0


def test_xiao_1995_table_3_all_chromium_conversion_reactions_are_source_exact():
    expected = {
        "CrO(liq)=CrO(RS)": [77150.0, -33.5],
        "CrO1.5(s)=CrO1.5(RS)": [74967.0, -37.5],
        "CrO1.5(liq)=CrO1.5(RS)": [10152.0, -12.6],
        "CrO1.5(s)=CrO1.5(liq)": [64815.0, -24.9],
    }
    assert CASES["xiao_holappa_1995"]["conversion_A_plus_B_T_J_per_mol"] == expected
    assert getattr(data, "CHROMIUM_OXIDE_CONVERSION_REACTIONS", {}) == {
        reaction: tuple(coefficients) for reaction, coefficients in expected.items()
    }
    assert CONVERSION_DEFAULTS["CrO"] == (77150.0, -33.5)
    assert CONVERSION_DEFAULTS["CrO1.5"] == (74967.0, -37.5)
    assert data.CONVERSION_STANDARD_STATES["CrO"] == ("CrO(liquid), Xiao & Holappa (1995) Table 3")
    assert data.CONVERSION_STANDARD_STATES["CrO1.5"] == (
        "CrO1.5(solid), Xiao & Holappa (1995) Table 3"
    )


def test_xiao_1995_silica_conversion_uses_the_correct_paper_label():
    assert CASES["xiao_holappa_1995"]["silica_conversion_A_plus_B_T_J_per_mol"] == [
        51346.0,
        -13.88,
    ]
    assert SIO2_CONVERSIONS["workbook"] == (51346.0, -13.88)
    assert data.SIO2_CONVERSION_STANDARD_STATES["workbook"] == (
        "SiO2(solid), polymorph unspecified; Xiao & Holappa (1995), silica-in-lime-silica section"
    )


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


def test_compute_state_uses_the_selected_parameter_profile_without_global_mutation():
    cfg = load_input(REFERENCE_INPUT)
    cfg.options.parameter_profile = "banya93_casi_xiao95cr_hybrid_v1"
    state, _ = compute_state(cfg, r_Fe=0.10, r_Cr=1.0, C_wtpc=0.02)
    expected = rsm_gamma_rtln(
        state.X, cfg.options.cross_terms, parameter_profile=cfg.options.parameter_profile
    )
    default = rsm_gamma_rtln(state.X, cfg.options.cross_terms)
    assert state.rtln_gamma_rs == pytest.approx(expected, rel=1e-12, abs=1e-9)
    assert expected["Ca2+"] != default["Ca2+"]


def test_chromium_coefficients_keep_oxide_and_metal_sources_separate():
    sources = {source["id"]: source for source in CASES["sources"]}
    oxide_source = sources["xiao_holappa_1995"]
    dissolution_source = sources["xiao_kou_fang_2018"]
    sigworth_source = sources["sigworth_elliott_1974"]
    assert any("p. 322" in use and "Eqs. 9" in use for use in oxide_source["used_for"])
    assert any("p. 419" in use and "Table 1" in use for use in dissolution_source["used_for"])
    assert any("Wagner" in use for use in sigworth_source["used_for"])
    assert not any("dissolv" in use.lower() for use in sigworth_source["used_for"])

    oxide = CASES["xiao_holappa_1995"]["chromium_redox"]
    dissolution = CASES["xiao_kou_fang_2018"]["chromium_dissolution_for_library_standard"]
    assert oxide["delta_g_A_plus_B_T_cal_per_mol"] == [25690.0, -13.36]
    assert dissolution["delta_g_A_plus_B_T_J_per_mol"] == [19246.0, -46.86]


def test_xiao_chromium_equilibrium_uses_the_library_metal_standard_state():
    temperature = 1823.15
    oxide = CASES["xiao_holappa_1995"]["chromium_redox"]
    dissolution = CASES["xiao_kou_fang_2018"]["chromium_dissolution_for_library_standard"]
    a_cal, b_cal = oxide["delta_g_A_plus_B_T_cal_per_mol"]
    a_dis, b_dis = dissolution["delta_g_A_plus_B_T_J_per_mol"]
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


def test_compiled_manganese_carbon_parameter_is_separate():
    row, column = "Mn", "C"
    from slag_model.data import METAL_SPECIES

    i = METAL_SPECIES.index(row)
    j = METAL_SPECIES.index(column)
    assert EIJ[i, j] == -0.07


_XIAO1995_TABLE1 = [
    (
        {"CrO": 14.16, "CrO1.5": 67.54, "SiO2": 18.30},
        {"CrO": 0.91, "CrO1.5": 0.95},
        0.12,
        0.17,
        -12.06,
    ),
    (
        {"CrO": 8.55, "CrO1.5": 75.93, "SiO2": 15.52},
        {"CrO": 0.94, "CrO1.5": 1.00},
        0.06,
        0.10,
        -12.03,
    ),
    (
        {"CrO": 19.48, "CrO1.5": 59.32, "SiO2": 21.20},
        {"CrO": 0.80, "CrO1.5": 0.79},
        0.19,
        0.25,
        -12.17,
    ),
    (
        {"CrO": 11.99, "CrO1.5": 61.31, "SiO2": 26.70},
        {"CrO": 0.72, "CrO1.5": 0.67},
        0.37,
        0.16,
        -12.26,
    ),
    (
        {"CrO": 18.88, "CrO1.5": 54.59, "SiO2": 26.53},
        {"CrO": 0.66, "CrO1.5": 0.59},
        0.36,
        0.26,
        -12.33,
    ),
    (
        {"CrO": 19.46, "CrO1.5": 45.01, "SiO2": 35.53},
        {"CrO": 0.54, "CrO1.5": 0.43},
        0.66,
        0.30,
        -12.51,
    ),
    (
        {"CrO": 10.48, "CrO1.5": 50.82, "SiO2": 38.70},
        {"CrO": 0.55, "CrO1.5": 0.44},
        0.75,
        0.17,
        -12.50,
    ),
    (
        {"CrO": 13.16, "CrO1.5": 63.37, "SiO2": 23.47},
        {"CrO": 0.77, "CrO1.5": 0.74},
        0.26,
        0.17,
        -12.20,
    ),
    (
        {"CrO": 13.35, "CrO1.5": 57.48, "SiO2": 29.17},
        {"CrO": 0.62, "CrO1.5": 0.53},
        0.46,
        0.19,
        -12.39,
    ),
    (
        {"CrO": 7.23, "CrO1.5": 55.74, "SiO2": 37.03},
        {"CrO": 0.60, "CrO1.5": 0.51},
        0.70,
        0.11,
        -12.42,
    ),
    (
        {"CrO": 13.10, "CrO1.5": 47.40, "SiO2": 39.50},
        {"CrO": 0.55, "CrO1.5": 0.44},
        0.77,
        0.22,
        -12.50,
    ),
]


def test_xiao_1995_table_1_all_eleven_rows_are_transcribed_exactly():
    source = CASES["xiao_holappa_1995"]
    rows = source["experimental_cases_1873_K"]
    assert len(rows) == 11
    for row, (mol_percent, measured_activity, siO2_activity, cr2_fraction, log_po2) in zip(
        rows, _XIAO1995_TABLE1, strict=True
    ):
        assert row["mol_percent"] == mol_percent
        assert row["measured_activity"] == measured_activity
        assert row["gibbs_duhem_calculated_activity_SiO2"] == siO2_activity
        assert row["reported_Cr2plus_fraction"] == cr2_fraction
        assert row["reported_log10_P_O2_atm"] == log_po2


def test_xiao_1995_pure_solid_chromium_redox_identity_within_printed_rounding():
    source = CASES["xiao_holappa_1995"]
    assert source["measured_activity_decimal_places"] == 2
    a_cal, b_cal = source["chromium_redox"]["delta_g_A_plus_B_T_cal_per_mol"]
    temperature = source["temperature_K"]
    k_solid_cr = math.exp(-4.184 * (a_cal + b_cal * temperature) / (R * temperature))
    half_printed_step = 0.5 * 10 ** (-source["measured_activity_decimal_places"])

    for case in source["experimental_cases_1873_K"]:
        activity = case["measured_activity"]
        a_cro = activity["CrO"]
        a_cro15 = activity["CrO1.5"]
        q_min = (a_cro - half_printed_step) ** 3 / (a_cro15 + half_printed_step) ** 2
        q_max = (a_cro + half_printed_step) ** 3 / (a_cro15 - half_printed_step) ** 2
        assert q_min <= k_solid_cr <= q_max, (
            f"{case['id']}: K={k_solid_cr:.6f}, source-rounding interval=[{q_min:.6f}, {q_max:.6f}]"
        )


def _xiao1995_model_activities(case):
    mol_percent = case["mol_percent"]
    total = sum(mol_percent.values())
    x_by_cation = {
        "Si4+": mol_percent["SiO2"] / total,
        "Cr2+": mol_percent["CrO"] / total,
        "Cr3+": mol_percent["CrO1.5"] / total,
    }
    full_x = {cation: x_by_cation.get(cation, 0.0) for cation in CATION_ORDER}
    rtln = rsm_gamma_rtln(x_by_cation, "full")
    activities = build_slag_activities(rtln, full_x, 1873.0, sio2_conversion="workbook")
    return {
        species: activities.a_conventional_by_species[species]
        for species in ("CrO", "CrO1.5", "SiO2")
    }


def test_xiao_1995_table_1_comparison_is_labeled_in_sample_and_not_thresholded():
    source = CASES["xiao_holappa_1995"]
    assert source["experimental_data_role"] == (
        "in_sample_parameter_assessment_data; not independent validation"
    )
    assert source["model_activity_comparison_basis"] == (
        "regular-solution calculation using Table 2 parameters and Table 3 conversions"
    )
    assert source["activity_provenance"] == {
        "CrO": "EMF-derived activity measurement",
        "CrO1.5": "EMF-derived activity measurement",
        "SiO2": "calculated from Gibbs-Duhem, not directly measured",
    }
    assert "activity_acceptance_tolerance" not in CASES["validation_policy"]
    expected = CASES["model_diagnostic_snapshots"]["xiao_holappa_1995_1873_K"]["activities"]
    for case in source["experimental_cases_1873_K"]:
        calculated = _xiao1995_model_activities(case)
        assert calculated == pytest.approx(expected[case["id"]], abs=5e-7)


def test_banya_ca_si_analytical_case_is_not_measured_or_independent_validation():
    case = CASES["analytical_cases"]["banya_1993_ca_si_equal_cations"]
    assert case["role"] == "equation_derived_analytical_illustration_not_measured_or_validation"
    assert "measured_activity" not in case
    assert case["cation_fractions"] == {"Ca2+": 0.5, "Si4+": 0.5}
    X = case["cation_fractions"]
    assert "parameter_profile" in __import__("inspect").signature(rsm_gamma_rtln).parameters
    rtln_gamma = rsm_gamma_rtln(X, "full", parameter_profile=case["parameter_profile"])
    full_x = {cation: X.get(cation, 0.0) for cation in CATION_ORDER}
    activities = build_slag_activities(
        rtln_gamma,
        full_x,
        case["temperature_K"],
        sio2_conversion="banya",
    )
    assert activities.a_conventional_by_species["SiO2"] == pytest.approx(
        case["expected_activity"]["SiO2"], rel=1e-10
    )
    assert activities.a_conventional_by_species["CaO"] == pytest.approx(
        case["expected_activity"]["CaO"], rel=1e-10
    )
