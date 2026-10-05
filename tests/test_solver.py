"""Balanced solver integration tests (PLAN.md test 5, 6)."""

from __future__ import annotations

import math
import pathlib

import pytest

from slag_model.io import load_input
from slag_model.redox import cr_ratio, fe_ratio_ban_ya
from slag_model.rsm import SlagActivities
from slag_model.solver import compute_state, solve

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


def test_xfe_oxygen_estimate_uses_result_iron_activity():
    cfg = load_input(SLAG_PATH)
    cfg.options.a_Fe = "xfe"
    result = solve(cfg)
    expected = (result.a_FeO / result.a_Fe) * 10.0 ** (-6320.0 / cfg.temperature_K + 2.734)
    assert result.O_wtpc == pytest.approx(expected, rel=1e-12)
    assert result.O_FeO_equivalent_wtpc == pytest.approx(expected, rel=1e-12)
    assert result.calculation_mode == "coupled"


def test_QK_FeO_unity(default_result):
    assert default_result.q_over_k["FeO"] == pytest.approx(1.0, abs=1e-6)


def test_QK_CO_unity(default_result):
    assert default_result.q_over_k["CO"] == pytest.approx(1.0, abs=1e-6)


def test_solve_result_uses_explicit_activity_bases(default_result):
    assert isinstance(default_result.activities, SlagActivities)
    assert default_result.a_FeO == pytest.approx(
        default_result.activities.a_conventional_by_species["FeO"]
    )


def test_compute_state_uses_rs_fe_and_conventional_chromium_activities():
    cfg = load_input(SLAG_PATH)
    result, targets = compute_state(cfg, r_Fe=0.10, r_Cr=1.0, C_wtpc=0.02)
    activities = result.activities

    expected_fe = fe_ratio_ban_ya(
        cfg.temperature_K,
        result.P_O2_atm,
        activities.gamma_rs_by_cation["Fe2+"],
        activities.gamma_rs_by_cation["Fe3+"],
    )
    assert targets["r_Fe"] == pytest.approx(expected_fe, rel=1e-12)

    cr_rhs = (
        result.k_Cr
        * result.a_Cr
        * activities.gamma_conventional_by_species["CrO1.5"] ** 2
        * result.N
        / (activities.gamma_conventional_by_species["CrO"] ** 3 * result.split.n_totals["Cr"])
    )
    assert targets["r_Cr"] == pytest.approx(cr_ratio(cr_rhs), rel=1e-12)


def test_solve_rejects_nonfinite_trailing_target(monkeypatch):
    import slag_model.solver as solver_module

    cfg = load_input(SLAG_PATH)
    valid_state, _ = compute_state(cfg, r_Fe=0.10, r_Cr=1.0, C_wtpc=0.02)
    monkeypatch.setattr(
        solver_module,
        "compute_state",
        lambda *_args: (valid_state, {"r_Fe": 0.10, "r_Cr": 1.0, "C": float("nan")}),
    )

    with pytest.raises(RuntimeError) as exc_info:
        solver_module.solve(cfg)
    assert type(exc_info.value).__name__ == "NumericalStateError"


def test_fixed_po2_nonconvergence_state_keeps_fixed_pressure_diagnostics(monkeypatch):
    import slag_model.solver as solver_module
    from slag_model.equilibrium import k_CO, k_FeO
    from slag_model.solver import NonConvergenceError, solve_redox_fixed_po2

    cfg = load_input(SLAG_PATH)
    imposed_pressure = 1.0e-9
    monkeypatch.setattr(solver_module, "CONVERGENCE_TOL", -1.0)
    with pytest.raises(NonConvergenceError) as exc_info:
        solve_redox_fixed_po2(cfg, P_O2=imposed_pressure, C=0.018)

    state = exc_info.value.state
    expected_fe = (state.a_FeO / (state.a_Fe * math.sqrt(imposed_pressure))) / k_FeO(
        cfg.temperature_K
    )
    expected_co = (cfg.P_CO_atm / (state.a_C * math.sqrt(imposed_pressure))) / k_CO(
        cfg.temperature_K
    )
    assert state.P_O2_atm == imposed_pressure
    assert state.calculation_mode == "fixed_P_O2"
    assert state.O_wtpc is None
    assert state.q_over_k["FeO"] == pytest.approx(expected_fe, rel=1e-12)
    assert state.q_over_k["CO"] == pytest.approx(expected_co, rel=1e-12)


def test_reference_case_after_p0_corrections():
    cfg = load_input(SLAG_PATH)
    result = solve(cfg)

    # Regression sentinel for the selected default parameter profile, not empirical validation.
    assert cfg.options.parameter_profile == "banya93_xiao95_v1"
    assert result.k_Cr == pytest.approx(0.00878747397894681, rel=1e-10)
    assert result.r_Fe == pytest.approx(0.13100535169281785, rel=2e-10)
    assert result.r_Cr == pytest.approx(0.2491105269636371, rel=2e-10)
    assert result.C_wtpc == pytest.approx(0.03079811331754755, rel=2e-10)
    assert result.P_O2_atm == pytest.approx(3.718198690650037e-10, rel=2e-10)

    cations = result.split.n_cations
    assert cations["Fe2+"] + cations["Fe3+"] == pytest.approx(
        result.split.n_totals["Fe"], abs=1e-12
    )
    assert cations["Cr2+"] + cations["Cr3+"] == pytest.approx(
        result.split.n_totals["Cr"], abs=1e-12
    )
    assert cfg.options.ti_handling == "exclude_renormalize"
    assert sum(result.X.values()) == pytest.approx(1.0, abs=1e-12)

    activities = result.activities
    q_over_k_fe = activities.a_conventional_by_species["FeO"] / (
        result.a_Fe * math.sqrt(result.P_O2_atm) * result.k_FeO
    )
    q_over_k_co = cfg.P_CO_atm / (result.a_C * math.sqrt(result.P_O2_atm) * result.k_CO)
    q_over_k_cr = activities.a_conventional_by_species["CrO"] ** 3 / (
        activities.a_conventional_by_species["CrO1.5"] ** 2 * result.a_Cr * result.k_Cr
    )
    assert q_over_k_fe == pytest.approx(1.0, abs=1e-6)
    assert q_over_k_co == pytest.approx(1.0, abs=1e-6)
    assert q_over_k_cr == pytest.approx(1.0, abs=1e-6)
