"""Input validation tests (PLAN.md section 5)."""

from __future__ import annotations

import json
import math

import pytest

from slag_model.constants import R
from slag_model.data import CONVERSION_STANDARD_STATES
from slag_model.io import InputError, build_report, parse_input
from slag_model.solver import solve

VALID_JSON = {
    "temperature_K": 1823.15,
    "P_CO_atm": 1.0,
    "slag_wtpc": {
        "SiO2": 27.495,
        "CaO": 35.597,
        "MgO": 1.999,
        "Al2O3": 7.973,
        "MnO": 4.762,
        "TiO2": 0.725,
        "P2O5": 0.0,
        "FeO_total": 18.851,
        "Cr2O3_total": 2.58,
    },
    "metal_wtpc": {"Cr": 0.038, "Mn": 0.007, "P": 0.069},
    "options": {
        "cross_terms": "full",
        "ti_handling": "exclude_renormalize",
        "sio2_conversion": "workbook",
    },
}


def test_valid_input():
    cfg = parse_input(json.dumps(VALID_JSON))
    assert cfg.temperature_K == 1823.15
    assert cfg.P_CO_atm == 1.0
    assert cfg.slag_wtpc["SiO2"] == 27.495
    assert cfg.metal_wtpc["Cr"] == 0.038


def test_defaults():
    minimal = {
        "slag_wtpc": {
            "SiO2": 40.0,
            "CaO": 40.0,
            "MgO": 5.0,
            "Al2O3": 5.0,
            "FeO_total": 10.0,
            "Cr2O3_total": 0.01,
        },
        "metal_wtpc": {},
    }
    cfg = parse_input(json.dumps(minimal))
    assert cfg.temperature_K == 1823.15
    assert cfg.P_CO_atm == 1.0
    assert cfg.metal_wtpc == {"Cr": 0.0, "Mn": 0.0, "P": 0.0}


def test_unknown_top_key():
    d = {**VALID_JSON, "bogus": 42}
    with pytest.raises(InputError, match="unknown top-level"):
        parse_input(json.dumps(d))


def test_exactly_one_iron():
    d = {**VALID_JSON}
    d["slag_wtpc"] = {k: v for k, v in d["slag_wtpc"].items() if k != "FeO_total"}
    with pytest.raises(InputError, match="exactly one of FeO_total"):
        parse_input(json.dumps(d))


def test_exactly_one_chromium():
    d = {**VALID_JSON}
    d["slag_wtpc"] = {k: v for k, v in d["slag_wtpc"].items() if k != "Cr2O3_total"}
    with pytest.raises(InputError, match="exactly one of Cr2O3_total"):
        parse_input(json.dumps(d))


def test_negative_slag_value():
    d = {**VALID_JSON, "slag_wtpc": {**VALID_JSON["slag_wtpc"], "SiO2": -1.0}}
    with pytest.raises(InputError, match="non-negative"):
        parse_input(json.dumps(d))


def test_unknown_slag_key():
    d = {**VALID_JSON, "slag_wtpc": {**VALID_JSON["slag_wtpc"], "BOGUS": 1.0}}
    with pytest.raises(InputError, match="unknown slag keys"):
        parse_input(json.dumps(d))


def test_unknown_metal_key():
    d = {**VALID_JSON, "metal_wtpc": {"Sn": 0.1}}
    with pytest.raises(InputError, match="unknown metal keys"):
        parse_input(json.dumps(d))


def test_invalid_cross_terms():
    d = {**VALID_JSON, "options": {"cross_terms": "bogus"}}
    with pytest.raises(InputError, match="cross_terms"):
        parse_input(json.dumps(d))


def test_slag_sum_warning():
    d = {**VALID_JSON, "slag_wtpc": {**VALID_JSON["slag_wtpc"], "SiO2": 10.0}}
    cfg = parse_input(json.dumps(d))
    assert any("slag wt% sum" in w for w in cfg.warnings)


def test_al2o3_conversion_warning():
    cfg = parse_input(json.dumps(VALID_JSON))
    assert any("Al2O3" in w for w in cfg.warnings)


def test_fe_total_input_form():
    d = {**VALID_JSON, "slag_wtpc": {**VALID_JSON["slag_wtpc"]}}
    del d["slag_wtpc"]["FeO_total"]
    d["slag_wtpc"]["Fe_total"] = 14.68
    d["options"] = {**d["options"], "iron_input": "Fe_total"}
    cfg = parse_input(json.dumps(d))
    assert cfg.options.iron_input == "Fe_total"


def test_report_v2_distinguishes_rs_and_formula_unit_activities():
    slag = {**VALID_JSON["slag_wtpc"], "CaO": 34.597, "P2O5": 1.0}
    cfg = parse_input(json.dumps({**VALID_JSON, "slag_wtpc": slag}))
    report = build_report(cfg, solve(cfg))

    assert report["schema_version"] == 2
    rows = {row["cation"]: row for row in report["solution"]["components"]}
    phosphorus = rows["P5+"]
    assert phosphorus["rs_species"] == "PO2.5"
    assert phosphorus["conventional"]["species"] == "P2O5"
    assert phosphorus["conventional"]["standard_state"] == CONVERSION_STANDARD_STATES["P2O5"]
    assert phosphorus["conventional"]["rs_units_per_species"] == 2
    assert "gamma" not in phosphorus["conventional"]
    assert "gamma_fraction_basis" not in phosphorus["conventional"]
    assert "gamma" not in phosphorus
    assert "a" not in phosphorus
    assert phosphorus["conventional"]["a"] == pytest.approx(
        phosphorus["a_RS"] ** 2
        * math.exp((52720.0 - 230.706 * cfg.temperature_K) / (R * cfg.temperature_K)),
        rel=1e-12,
    )
    assert {
        "rs_species",
        "X_cation",
        "RTln_gamma_RS_J_per_mol_cation",
        "gamma_RS",
        "a_RS",
    }.issubset(phosphorus)

    aluminum = rows["Al3+"]
    assert aluminum["rs_species"] == "AlO1.5"
    assert aluminum["conventional"] is None
    assert aluminum["conventional_unavailable_reason"] == "no documented conversion"
    assert "gamma" not in aluminum
    assert "a" not in aluminum

    iron = rows["Fe2+"]["conventional"]
    assert iron["gamma_fraction_basis"] == "X_cation"
    assert "gamma" in iron


def test_as_excel_titanium_report_is_denominator_only_without_activity():
    cfg = parse_input(json.dumps(VALID_JSON))
    cfg.options.ti_handling = "as_excel"
    report = build_report(cfg, solve(cfg))

    titanium = next(row for row in report["solution"]["components"] if row.get("cation") == "Ti4+")
    assert titanium["model_status"] == "denominator_only_legacy"
    assert not {
        "RTln_gamma_RS_J_per_mol_cation",
        "gamma_RS",
        "a_RS",
        "conventional",
        "conventional_unavailable_reason",
        "gamma",
        "a",
    }.intersection(titanium)


@pytest.mark.parametrize("option_name", ["fe2o3_conversion", "al2o3_conversion"])
def test_custom_conversion_requires_standard_state_migration(option_name):
    options = {
        **VALID_JSON["options"],
        option_name: {"A": 100.0, "B": -0.1},
    }
    payload = {**VALID_JSON, "options": options}

    with pytest.raises(InputError, match="migration"):
        parse_input(json.dumps(payload))


@pytest.mark.parametrize("option_name", ["fe2o3_conversion", "al2o3_conversion"])
def test_labeled_custom_conversion_round_trips_through_report(option_name):
    conversion = {
        "A": 100.0,
        "B": -0.1,
        "standard_state": "custom reference state, specified by caller",
    }
    options = {**VALID_JSON["options"], option_name: conversion}
    cfg = parse_input(json.dumps({**VALID_JSON, "options": options}))
    report = build_report(cfg, solve(cfg))

    assert report["input"]["options"][option_name] == conversion
    reparsed = parse_input(json.dumps(report["input"]))
    assert getattr(reparsed.options, option_name) == getattr(cfg.options, option_name)
