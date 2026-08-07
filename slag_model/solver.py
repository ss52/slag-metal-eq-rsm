"""Damped fixed-point solver (PLAN.md section 8) for the coupled slag-metal model."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .data import CATION_TO_OXIDE
from .equilibrium import k_CO, k_Cr, k_FeO, p_o2_from_slag
from .metal import (
    activity_C,
    activity_Cr,
    balanced_carbon,
    dissolved_oxygen,
    wipf_factors,
    x_fe_metal,
)
from .redox import cr_ratio, fe_ratio_ban_ya
from .rsm import (
    SlagSplit,
    cation_fractions,
    convert_gammas,
    rsm_gamma_rtln,
    split_slag,
)

MAX_ITERATIONS: int = 200
DAMPING: float = 0.4
CONVERGENCE_TOL: float = 1.0e-8


class NonConvergenceError(RuntimeError):
    """Raised when the fixed-point iteration fails to converge."""

    def __init__(self, iterations: int, state: SolveResult, message: str):
        super().__init__(message)
        self.iterations = iterations
        self.state = state


@dataclass
class SolveResult:
    """Full converged state of a slag-metal solve (full or fixed-P_O2)."""

    r_Fe: float
    r_Cr: float
    C_wtpc: float
    P_O2_atm: float
    a_Fe: float
    a_FeO: float
    iterations: int
    split: SlagSplit
    X: dict[str, float]
    N: float
    rtln_gamma_rs: dict[str, float]
    delta_g_conv: dict[str, float]
    gamma: dict[str, float]
    a_slag: dict[str, float]
    f_metal: dict[str, float]
    a_C: float
    a_Cr: float
    O_wtpc: float
    k_FeO: float
    k_CO: float
    k_Cr: float
    q_over_k: dict[str, float | None] = field(default_factory=dict)


def compute_state(
    cfg, r_Fe: float, r_Cr: float, C_wtpc: float
) -> tuple[SolveResult, dict[str, float]]:
    """Evaluate all slag/metal quantities at (r_Fe, r_Cr, C) and return them plus
    the next-iteration targets (r_Fe_n, r_Cr_n, C_n).
    """
    o = cfg.options
    T = cfg.temperature_K

    split = split_slag(cfg.slag_wtpc, o.iron_input, o.chromium_input, r_Fe, r_Cr)
    X, N = cation_fractions(split.n_cations, o.ti_handling)
    rtln_rs = rsm_gamma_rtln(X, o.cross_terms)
    gamma, delta_g = convert_gammas(
        rtln_rs,
        T,
        sio2_conversion=o.sio2_conversion,
        fe2o3_conversion=o.fe2o3_conversion,
        al2o3_conversion=o.al2o3_conversion,
    )
    if "Ti4+" in X:
        gamma["TiO2"] = 1.0
        delta_g["TiO2"] = 0.0

    a_FeO = gamma["FeO"] * X["Fe2+"]
    a_Fe = 1.0 if o.a_Fe == "unity" else x_fe_metal(cfg.metal_wtpc, C_wtpc)

    k_FeO_val = k_FeO(T)
    P_O2 = p_o2_from_slag(a_FeO, T, a_Fe)

    r_Fe_n = fe_ratio_ban_ya(T, P_O2, gamma["FeO"], gamma["FeO1.5"])

    f = wipf_factors(cfg.metal_wtpc, C_wtpc)
    a_C = activity_C(f, C_wtpc)
    a_Cr = activity_Cr(f, cfg.metal_wtpc)

    n_Cr_tot = split.n_totals["Cr"]
    if n_Cr_tot > 0.0:
        rhs = k_Cr(T) * a_Cr * gamma["CrO1.5"] ** 2 * N / (gamma["CrO"] ** 3 * n_Cr_tot)
        r_Cr_n = cr_ratio(rhs)
    else:
        r_Cr_n = 1.0

    k_CO_val = k_CO(T)
    C_n = balanced_carbon(cfg.P_CO_atm, k_CO_val, P_O2, f["C"])

    a_slag: dict[str, float] = {}
    for cation, oxide in CATION_TO_OXIDE.items():
        if cation in X:
            a_slag[oxide] = gamma[oxide] * X[cation]
    if "Ti4+" in X:
        a_slag["TiO2"] = X["Ti4+"]

    q_over_k: dict[str, float | None] = {
        "FeO": (a_FeO / (a_Fe * math.sqrt(P_O2))) / k_FeO_val,
        "CO": (cfg.P_CO_atm / (a_C * math.sqrt(P_O2))) / k_CO_val,
    }
    if a_Cr > 0.0 and a_slag.get("CrO", 0.0) > 0.0:
        q_over_k["Cr"] = (a_slag["CrO"] ** 3 / (a_slag["CrO1.5"] ** 2 * a_Cr)) / k_Cr(T)
    else:
        q_over_k["Cr"] = None

    result = SolveResult(
        r_Fe=r_Fe,
        r_Cr=r_Cr,
        C_wtpc=C_wtpc,
        P_O2_atm=P_O2,
        a_Fe=a_Fe,
        a_FeO=a_FeO,
        iterations=0,
        split=split,
        X=X,
        N=N,
        rtln_gamma_rs=rtln_rs,
        delta_g_conv=delta_g,
        gamma=gamma,
        a_slag=a_slag,
        f_metal=f,
        a_C=a_C,
        a_Cr=a_Cr,
        O_wtpc=dissolved_oxygen(a_FeO, T),
        k_FeO=k_FeO_val,
        k_CO=k_CO_val,
        k_Cr=k_Cr(T),
        q_over_k=q_over_k,
    )
    targets = {"r_Fe": r_Fe_n, "r_Cr": r_Cr_n, "C": C_n}
    return result, targets


def _relative_changes(targets: dict[str, float], current: dict[str, float]) -> dict[str, float]:
    return {k: abs(targets[k] - current[k]) / max(abs(current[k]), 1.0e-12) for k in targets}


def solve(cfg) -> SolveResult:
    """Damped fixed-point driver with the slag-derived P_O2 closure (sec. 8)."""
    r_Fe, r_Cr, C = 0.10, 1.0, 0.02
    state: SolveResult | None = None
    targets: dict[str, float] = {"r_Fe": r_Fe, "r_Cr": r_Cr, "C": C}
    changes: dict[str, float] = {"r_Fe": 1.0, "r_Cr": 1.0, "C": 1.0}

    for it in range(1, MAX_ITERATIONS + 1):
        state, targets = compute_state(cfg, r_Fe, r_Cr, C)
        current = {"r_Fe": r_Fe, "r_Cr": r_Cr, "C": C}
        changes = _relative_changes(targets, current)
        if max(changes.values()) < CONVERGENCE_TOL:
            final, _ = compute_state(cfg, targets["r_Fe"], targets["r_Cr"], targets["C"])
            final.iterations = it
            return final
        r_Fe += DAMPING * (targets["r_Fe"] - r_Fe)
        r_Cr += DAMPING * (targets["r_Cr"] - r_Cr)
        C += DAMPING * (targets["C"] - C)

    assert state is not None
    raise NonConvergenceError(
        MAX_ITERATIONS,
        state,
        f"did not converge in {MAX_ITERATIONS} iterations; "
        f"last targets={targets}, last changes={changes}",
    )


def solve_redox_fixed_po2(cfg, P_O2: float, C: float = 0.018) -> SolveResult:
    """Damped fixed-point solve of the redox system at a *fixed* oxygen potential.

    Used by the workbook regression tests (PLAN.md sec. 10 item 2): iron and
    chromium ratios are brought to consistency at the given P_O2 while the rest
    of the slag composition is updated accordingly. P_CO is irrelevant here.
    """
    r_Fe, r_Cr = 0.10, 1.0
    it = 0
    current_state: SolveResult | None = None

    for it in range(1, MAX_ITERATIONS + 1):
        current_state, _ = compute_state(cfg, r_Fe, r_Cr, C)
        # r_Fe target from Ban-ya Eq. 24 at the EXTERNAL P_O2 (overwrite the
        # slag-derived one used inside compute_state).
        r_Fe_n = fe_ratio_ban_ya(
            cfg.temperature_K,
            P_O2,
            current_state.gamma["FeO"],
            current_state.gamma["FeO1.5"],
        )
        n_Cr_tot = current_state.split.n_totals["Cr"]
        if n_Cr_tot > 0.0:
            rhs = (
                k_Cr(cfg.temperature_K)
                * current_state.a_Cr
                * current_state.gamma["CrO1.5"] ** 2
                * current_state.N
                / (current_state.gamma["CrO"] ** 3 * n_Cr_tot)
            )
            r_Cr_n = cr_ratio(rhs)
        else:
            r_Cr_n = 1.0

        rel_fe = abs(r_Fe_n - r_Fe) / max(abs(r_Fe), 1.0e-12)
        rel_cr = abs(r_Cr_n - r_Cr) / max(abs(r_Cr), 1.0e-12)
        if rel_fe < CONVERGENCE_TOL and rel_cr < CONVERGENCE_TOL:
            final, _ = compute_state(cfg, r_Fe_n, r_Cr_n, C)
            final.iterations = it
            final.P_O2_atm = P_O2  # diagnostic honesty
            return final
        r_Fe += DAMPING * (r_Fe_n - r_Fe)
        r_Cr += DAMPING * (r_Cr_n - r_Cr)

    assert current_state is not None
    raise NonConvergenceError(
        it,
        current_state,
        f"did not converge in {MAX_ITERATIONS} iterations at fixed P_O2={P_O2}",
    )
