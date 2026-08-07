"""Regular solution model: split, cation fractions, RSM gammas, conversion."""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .constants import CATIONS_PER_FORMULA, MOLAR_MASS_ELEMENT, MOLAR_MASS_OXIDE, R
from .data import (
    ALPHA,
    CATION_ORDER,
    CATION_TO_OXIDE,
    CONVERSION_DEFAULTS,
    INPUT_OXIDE_KEYS,
    OXIDE_TO_CATION,
    SIO2_CONVERSIONS,
)


@dataclass
class SlagSplit:
    """Cation molar amounts (per 100 g of as-received slag) and oxide wt% split."""

    n_cations: dict[str, float]
    n_totals: dict[str, float]
    wt_split: dict[str, float]


def split_slag(
    slag_wtpc: dict[str, float],
    iron_input: str,
    chromium_input: str,
    r_Fe: float,
    r_Cr: float,
) -> SlagSplit:
    """Split analytical slag input into cation molar amounts (PLAN.md 7.1).

    Amounts are moles per 100 g of as-received slag. Iron is split into Fe2+/Fe3+
    via r_Fe = n_Fe3+/n_Fe2+; chromium into Cr3+/Cr2+ via r_Cr = n_Cr2+/n_Cr3+.
    """
    n: dict[str, float] = {}
    for oxide in INPUT_OXIDE_KEYS:
        w = slag_wtpc.get(oxide, 0.0)
        cation = OXIDE_TO_CATION.get(oxide)
        if cation is None:
            # TiO2 maps to Ti4+ which is not a regular-solution cation
            cation = "Ti4+"
        n[cation] = CATIONS_PER_FORMULA[oxide] * w / MOLAR_MASS_OXIDE[oxide]

    w_fe = slag_wtpc[iron_input]
    if iron_input == "Fe_total":
        n_Fe_tot = w_fe / MOLAR_MASS_ELEMENT["Fe"]
    else:  # FeO_total
        n_Fe_tot = w_fe / MOLAR_MASS_OXIDE["FeO"]
    n["Fe2+"] = n_Fe_tot / (1.0 + r_Fe)
    n["Fe3+"] = n_Fe_tot - n["Fe2+"]

    w_cr = slag_wtpc[chromium_input]
    if chromium_input == "Cr_total":
        n_Cr_tot = w_cr / MOLAR_MASS_ELEMENT["Cr"]
    else:  # Cr2O3_total -> 2 cations per formula
        n_Cr_tot = w_cr / (MOLAR_MASS_OXIDE["Cr2O3"] / CATIONS_PER_FORMULA["Cr2O3"])
    n["Cr3+"] = n_Cr_tot / (1.0 + r_Cr)
    n["Cr2+"] = n_Cr_tot - n["Cr3+"]

    # Oxide wt% table (per 100 g basis). Start from the non-redox input values,
    # then add the redox-split oxides.
    wt: dict[str, float] = {}
    for oxide in INPUT_OXIDE_KEYS:
        wt[oxide] = slag_wtpc.get(oxide, 0.0)
    wt["FeO"] = n["Fe2+"] * MOLAR_MASS_OXIDE["FeO"]
    wt["Fe2O3"] = n["Fe3+"] * MOLAR_MASS_OXIDE["Fe2O3"] / CATIONS_PER_FORMULA["Fe2O3"]
    wt["CrO"] = n["Cr2+"] * MOLAR_MASS_OXIDE["CrO"]
    wt["Cr2O3"] = n["Cr3+"] * MOLAR_MASS_OXIDE["Cr2O3"] / CATIONS_PER_FORMULA["Cr2O3"]

    return SlagSplit(
        n_cations=n,
        n_totals={"Fe": n_Fe_tot, "Cr": n_Cr_tot},
        wt_split=wt,
    )


def cation_fractions(
    n_cations: dict[str, float], ti_handling: str
) -> tuple[dict[str, float], float]:
    """Cation mole fractions over the active cation set (PLAN.md 7.2, 11).

    With ti_handling="as_excel" the Ti4+ cation is included in N with zero alpha
    (reproduces the reference workbook); with "exclude_renormalize" (default) Ti
    is dropped from the cation sum and the remaining fractions renormalize.
    """
    if ti_handling == "as_excel":
        included = list(CATION_ORDER) + ["Ti4+"]
    elif ti_handling == "exclude_renormalize":
        included = list(CATION_ORDER)
    else:
        raise ValueError(f"unknown ti_handling: {ti_handling!r}")

    N = sum(n_cations.get(c, 0.0) for c in included)
    if N <= 0.0:
        raise ValueError("cation total is non-positive")
    X = {c: n_cations.get(c, 0.0) / N for c in included}
    return X, N


# Major5 cross-term pair set (exists only to reproduce the workbook regression).
MAJOR5_PAIRS: list[tuple[str, str]] = [
    ("Ca2+", "Si4+"),
    ("Ca2+", "Al3+"),
    ("Si4+", "Al3+"),
    ("Mg2+", "Si4+"),
    ("Ca2+", "Mg2+"),
]


def _pair_indices(cross_terms: str) -> list[tuple[int, int]]:
    """Return list of (j, k) cross-term pair indices into CATION_ORDER."""
    idx = {c: i for i, c in enumerate(CATION_ORDER)}
    n = len(CATION_ORDER)
    if cross_terms == "full":
        return [(j, k) for j in range(n) for k in range(j + 1, n)]
    if cross_terms == "major5":
        return [(idx[a], idx[b]) for (a, b) in MAJOR5_PAIRS]
    raise ValueError(f"unknown cross_terms: {cross_terms!r}")


def rsm_gamma_rtln(X: dict[str, float], cross_terms: str) -> dict[str, float]:
    """RSM activity coefficients on the regular-solution scale (PLAN.md 7.3).

    Returns R*T*ln(gamma_i^RS) in J/mol for each cation in CATION_ORDER.
    """
    Xv = np.array([X.get(c, 0.0) for c in CATION_ORDER])
    pairs = _pair_indices(cross_terms)

    out: dict[str, float] = {}
    for i, cation in enumerate(CATION_ORDER):
        # sum_j alpha_ij * X_j^2 (alpha_ii = 0 so the i term drops out)
        s = float(np.sum(ALPHA[i, :] * Xv**2))
        t = 0.0
        for j, k in pairs:
            if i in (j, k):
                continue
            t += (ALPHA[i, j] + ALPHA[i, k] - ALPHA[j, k]) * Xv[j] * Xv[k]
        out[cation] = s + t
    return out


def conversion_parameters(
    oxide: str,
    sio2_conversion: str,
    fe2o3_conversion: tuple[float, float] | None,
    al2o3_conversion: tuple[float, float] | None,
) -> tuple[float, float]:
    """Return (A, B) for the standard-state conversion of `oxide`."""
    if oxide == "SiO2":
        return SIO2_CONVERSIONS[sio2_conversion]
    if oxide == "FeO1.5":
        return fe2o3_conversion if fe2o3_conversion is not None else (0.0, 0.0)
    if oxide == "Al2O3":
        return al2o3_conversion if al2o3_conversion is not None else (0.0, 0.0)
    return CONVERSION_DEFAULTS.get(oxide, (0.0, 0.0))


def convert_gammas(
    rtln_gamma_rs: dict[str, float],
    T: float,
    sio2_conversion: str = "workbook",
    fe2o3_conversion: tuple[float, float] | None = None,
    al2o3_conversion: tuple[float, float] | None = None,
) -> tuple[dict[str, float], dict[str, float]]:
    """Convert R.S. gammas to the pure-stable-oxide standard state (PLAN.md 7.4).

    Returns (gamma_final keyed by oxide, DeltaG_conv keyed by oxide).
    """
    gamma: dict[str, float] = {}
    delta_g: dict[str, float] = {}
    for cation, rtln in rtln_gamma_rs.items():
        oxide = CATION_TO_OXIDE[cation]
        A, B = conversion_parameters(oxide, sio2_conversion, fe2o3_conversion, al2o3_conversion)
        dg = A + B * T
        delta_g[oxide] = dg
        gamma[oxide] = math.exp((rtln + dg) / (R * T))
    return gamma, delta_g
