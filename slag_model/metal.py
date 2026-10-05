"""Metal phase: WIPF interaction parameters, activities, carbon balance (7.7-7.10)."""

from __future__ import annotations

import math

from .constants import MOLAR_MASS_METAL
from .data import (
    EIJ,
    METAL_INPUT_KEYS,
    METAL_SPECIES,
    O_SAT_CONSTANT,
    O_SAT_T_COEF,
)


def wipf_factors(metal_wtpc: dict[str, float], C_wtpc: float) -> dict[str, float]:
    """First-order Wagner interaction parameter formalism (PLAN.md 6.5, 7.10).

    log10 f_i = sum_j e_i^j * [%j]
    """
    comp = {"C": C_wtpc, **{k: metal_wtpc.get(k, 0.0) for k in METAL_INPUT_KEYS}}
    f: dict[str, float] = {}
    for i, si in enumerate(METAL_SPECIES):
        logf = sum(EIJ[i, j] * comp[sj] for j, sj in enumerate(METAL_SPECIES))
        f[si] = 10.0**logf
    return f


def activity_C(f: dict[str, float], C_wtpc: float) -> float:
    """Henrian carbon activity a_C = f_C * [%C]."""
    return f["C"] * C_wtpc


def activity_Cr(f: dict[str, float], metal_wtpc: dict[str, float]) -> float:
    """Henrian chromium activity a_Cr = f_Cr * [%Cr]."""
    return f["Cr"] * metal_wtpc.get("Cr", 0.0)


def balanced_carbon(P_CO_atm: float, k_CO_val: float, P_O2: float, f_C: float) -> float:
    """Steel carbon in equilibrium with P_O2 via the C-O balance (PLAN.md 7.8).

    [%C] = P_CO / (K_CO * f_C * sqrt(P_O2))
    """
    return P_CO_atm / (k_CO_val * f_C * math.sqrt(P_O2))


def carbon_fixed_point(
    P_CO_atm: float,
    k_CO_val: float,
    P_O2: float,
    metal_wtpc: dict[str, float],
    passes: int = 3,
) -> tuple[float, float]:
    """Solve the C-O balance self-consistently for f_C (which depends on [%C]).

    f_C depends weakly on [%C] via e_C^C, so a few fixed-point passes starting
    from the WIPF factors evaluated at the current iterate converge quickly.
    """
    C = 0.02
    for _ in range(passes):
        C = balanced_carbon(P_CO_atm, k_CO_val, P_O2, wipf_factors(metal_wtpc, C)["C"])
    return C, wipf_factors(metal_wtpc, C)["C"]


def dissolved_oxygen(a_FeO: float, T: float, a_Fe: float = 1.0) -> float:
    """Approximate FeO-equilibrium oxygen diagnostic with f_O assumed unity."""
    return (a_FeO / a_Fe) * 10.0 ** (O_SAT_T_COEF / T + O_SAT_CONSTANT)


def x_fe_metal(metal_wtpc: dict[str, float], C_wtpc: float) -> float:
    """Mole fraction of Fe in the steel (for the a_Fe="xfe" option)."""
    comp = {"C": C_wtpc, **{k: metal_wtpc.get(k, 0.0) for k in METAL_INPUT_KEYS}}
    fe_wt = 100.0 - sum(comp.values())
    if fe_wt <= 0.0:
        raise ValueError("metal composition must leave a positive Fe balance")
    n_Fe = fe_wt / MOLAR_MASS_METAL["Fe"]
    n = {s: comp[s] / MOLAR_MASS_METAL[s] for s in ("C", "Cr", "Mn", "P")}
    return n_Fe / (n_Fe + sum(n.values()))
