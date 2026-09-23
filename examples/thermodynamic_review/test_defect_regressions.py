"""Executable reproductions for non-paper defects found during the audit."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from slag_model.equilibrium import k_CO, k_FeO
from slag_model.io import InputError, load_input, parse_input
from slag_model.redox import cr_ratio
from slag_model.rsm import split_slag
from slag_model.solver import solve_redox_fixed_po2

REPO_ROOT = Path(__file__).resolve().parents[2]
REFERENCE_INPUT = REPO_ROOT / "examples" / "eaf_slag_01.json"


def _payload() -> dict:
    return json.loads(REFERENCE_INPUT.read_text(encoding="utf-8"))


@pytest.mark.xfail(strict=True, reason="Known input defect: non-finite floats are accepted.")
def test_nonfinite_input_is_rejected():
    payload = _payload()
    payload["P_CO_atm"] = "NaN"
    with pytest.raises(InputError, match="finite"):
        parse_input(json.dumps(payload))


@pytest.mark.xfail(
    strict=True, reason="Known root-bracketing defect: roots below 1e-6 clamp to the lower bound."
)
def test_chromium_ratio_keeps_trace_roots():
    rhs = 1.0e-24
    ratio = cr_ratio(rhs)
    reconstructed = ratio**3 / (1.0 + ratio)
    assert reconstructed == pytest.approx(rhs, rel=1e-8, abs=0.0)


@pytest.mark.xfail(
    strict=True,
    reason="Known defect: fixed-P_O2 diagnostics use a different pressure.",
)
def test_fixed_po2_diagnostics_use_reported_pressure():
    config = load_input(REFERENCE_INPUT)
    state = solve_redox_fixed_po2(config, P_O2=1.0e-9, C=0.018)
    q_fe = (state.a_FeO / (state.a_Fe * math.sqrt(state.P_O2_atm))) / k_FeO(config.temperature_K)
    q_co = (config.P_CO_atm / (state.a_C * math.sqrt(state.P_O2_atm))) / k_CO(config.temperature_K)
    assert state.q_over_k["FeO"] == pytest.approx(q_fe, rel=1e-12)
    assert state.q_over_k["CO"] == pytest.approx(q_co, rel=1e-12)


def test_custom_conversion_report_round_trip():
    payload = _payload()
    payload["options"]["fe2o3_conversion"] = {
        "A": 100.0,
        "B": -0.1,
        "standard_state": "custom FeO1.5 reference state",
    }
    payload["options"]["al2o3_conversion"] = {
        "A": 200.0,
        "B": -0.2,
        "standard_state": "custom Al2O3 reference state",
    }
    config = parse_input(json.dumps(payload))
    payload["options"] = config.options.as_dict()
    reparsed = parse_input(json.dumps(payload))
    assert reparsed.options.fe2o3_conversion == config.options.fe2o3_conversion
    assert reparsed.options.al2o3_conversion == config.options.al2o3_conversion


@pytest.mark.xfail(
    strict=True,
    reason="Known validation defect: an empty options array is silently treated as an object.",
)
def test_options_array_is_rejected():
    payload = _payload()
    payload["options"] = []
    with pytest.raises(InputError, match="JSON object"):
        parse_input(json.dumps(payload))


@pytest.mark.xfail(
    strict=True,
    reason="Known defect: elemental totals produce a false oxide-sum warning.",
)
def test_element_total_inputs_do_not_trigger_false_slag_sum_warning():
    payload = _payload()
    slag = payload["slag_wtpc"]
    slag["Fe_total"] = slag.pop("FeO_total") * 55.845 / 71.844
    slag["Cr_total"] = slag.pop("Cr2O3_total") * (2.0 * 51.9961) / 151.9904
    payload["options"]["iron_input"] = "Fe_total"
    payload["options"]["chromium_input"] = "Cr_total"
    config = parse_input(json.dumps(payload))
    assert not any("slag wt% sum" in warning for warning in config.warnings)


@pytest.mark.xfail(
    strict=True,
    reason="Known cancellation defect: trace oxidized/reduced cations are lost by subtraction.",
)
def test_trace_redox_species_are_not_rounded_to_zero():
    config = load_input(REFERENCE_INPUT)
    split = split_slag(
        config.slag_wtpc,
        config.options.iron_input,
        config.options.chromium_input,
        r_Fe=1.0e-20,
        r_Cr=1.0e-20,
    )
    assert split.n_cations["Fe3+"] > 0.0
    assert split.n_cations["Cr2+"] > 0.0


@pytest.mark.xfail(
    strict=True,
    reason="Known defect: malformed JSON escapes the InputError contract.",
)
def test_malformed_json_raises_input_error():
    with pytest.raises(InputError):
        parse_input("{")
