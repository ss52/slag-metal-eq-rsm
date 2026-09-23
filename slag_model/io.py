"""JSON input parsing, validation, and report building (PLAN.md sections 5, 9)."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

from .data import (
    CHROMIUM_INPUT_KEYS,
    IRON_INPUT_KEYS,
    METAL_INPUT_KEYS,
)
from .rsm import ConversionSpec

KNOWN_TOP_KEYS = {
    "temperature_K",
    "P_CO_atm",
    "slag_wtpc",
    "metal_wtpc",
    "options",
}

INPUT_OXIDE_KEYS = {"SiO2", "CaO", "MgO", "Al2O3", "MnO", "TiO2", "P2O5"}
ALL_INPUT_OXIDE_KEYS = INPUT_OXIDE_KEYS | set(IRON_INPUT_KEYS) | set(CHROMIUM_INPUT_KEYS)

OPTION_KEYS = {
    "cross_terms",
    "ti_handling",
    "sio2_conversion",
    "fe2o3_conversion",
    "al2o3_conversion",
    "a_Fe",
    "iron_input",
    "chromium_input",
}


class InputError(ValueError):
    """Raised when the input JSON violates the specification."""


@dataclass
class Options:
    cross_terms: str = "full"
    ti_handling: str = "exclude_renormalize"
    sio2_conversion: str = "workbook"
    fe2o3_conversion: ConversionSpec | None = None
    al2o3_conversion: ConversionSpec | None = None
    a_Fe: str = "unity"
    iron_input: str = "FeO_total"
    chromium_input: str = "Cr2O3_total"

    def as_dict(self) -> dict:
        return {
            "cross_terms": self.cross_terms,
            "ti_handling": self.ti_handling,
            "sio2_conversion": self.sio2_conversion,
            "fe2o3_conversion": _conversion_as_dict(self.fe2o3_conversion),
            "al2o3_conversion": _conversion_as_dict(self.al2o3_conversion),
            "a_Fe": self.a_Fe,
            "iron_input": self.iron_input,
            "chromium_input": self.chromium_input,
        }


@dataclass
class SlagModelConfig:
    temperature_K: float
    P_CO_atm: float
    slag_wtpc: dict[str, float]
    metal_wtpc: dict[str, float]
    options: Options
    warnings: list[str] = field(default_factory=list)


def _conversion_as_dict(conversion: ConversionSpec | None) -> dict[str, float | str] | str:
    if conversion is None:
        return "none"
    return {
        "A": conversion.A,
        "B": conversion.B,
        "standard_state": conversion.standard_state,
    }


def _parse_custom_conversion(v, name: str) -> ConversionSpec | None:
    if v is None or v == "none":
        return None
    if isinstance(v, dict):
        required = {"A", "B", "standard_state"}
        missing = required - set(v)
        if missing:
            raise InputError(
                f"{name}: migration required; custom conversions must include numeric "
                "'A' and 'B' plus a non-empty 'standard_state' label"
            )
        extra = set(v) - required
        if extra:
            raise InputError(f"{name} has unknown fields: {sorted(extra)}")
        A, B = v["A"], v["B"]
        if any(isinstance(value, bool) or not isinstance(value, (int, float)) for value in (A, B)):
            raise InputError(f"{name}: 'A' and 'B' must be numbers")
        standard_state = v["standard_state"]
        if not isinstance(standard_state, str) or not standard_state.strip():
            raise InputError(f"{name}: 'standard_state' must be a non-empty string")
        return ConversionSpec(float(A), float(B), standard_state)
    raise InputError(
        f"{name} must be 'none' or an object with numeric 'A', numeric 'B', "
        "and a non-empty 'standard_state'"
    )


def _parse_options(raw: dict | None, warnings: list[str]) -> tuple[Options, set[str]]:
    raw = raw or {}
    if not isinstance(raw, dict):
        raise InputError("options must be a JSON object")
    unknown = set(raw) - OPTION_KEYS
    if unknown:
        raise InputError(f"unknown option keys: {sorted(unknown)}")

    cross_terms = raw.get("cross_terms", "full")
    if cross_terms not in ("full", "major5"):
        raise InputError(f"cross_terms must be 'full' or 'major5', got {cross_terms!r}")

    ti_handling = raw.get("ti_handling", "exclude_renormalize")
    if ti_handling not in ("exclude_renormalize", "as_excel"):
        raise InputError(
            f"ti_handling must be 'exclude_renormalize' or 'as_excel', got {ti_handling!r}"
        )

    sio2_conversion = raw.get("sio2_conversion", "workbook")
    if sio2_conversion not in ("workbook", "banya"):
        raise InputError(f"sio2_conversion must be 'workbook' or 'banya', got {sio2_conversion!r}")

    a_Fe = raw.get("a_Fe", "unity")
    if a_Fe not in ("unity", "xfe"):
        raise InputError(f"a_Fe must be 'unity' or 'xfe', got {a_Fe!r}")

    fe2o3 = _parse_custom_conversion(raw.get("fe2o3_conversion"), "fe2o3_conversion")
    al2o3 = _parse_custom_conversion(raw.get("al2o3_conversion"), "al2o3_conversion")

    if al2o3 is None:
        warnings.append(
            "Al2O3 conventional activity is unavailable because no documented "
            "conversion is configured."
        )

    iron_input = raw.get("iron_input", "FeO_total")
    if iron_input not in ("FeO_total", "Fe_total"):
        raise InputError(f"iron_input must be 'FeO_total' or 'Fe_total', got {iron_input!r}")

    chromium_input = raw.get("chromium_input", "Cr2O3_total")
    if chromium_input not in ("Cr2O3_total", "Cr_total"):
        raise InputError(
            f"chromium_input must be 'Cr2O3_total' or 'Cr_total', got {chromium_input!r}"
        )

    return Options(
        cross_terms=cross_terms,
        ti_handling=ti_handling,
        sio2_conversion=sio2_conversion,
        fe2o3_conversion=fe2o3,
        al2o3_conversion=al2o3,
        a_Fe=a_Fe,
        iron_input=iron_input,
        chromium_input=chromium_input,
    ), set(raw.keys())


def _validate_number(v, name: str, allow_negative: bool = False):
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        raise InputError(f"{name} must be a number, got {type(v).__name__}")
    if not allow_negative and v < 0:
        raise InputError(f"{name} must be non-negative, got {v}")


def parse_input(text: str) -> SlagModelConfig:
    """Parse and validate a JSON input string (PLAN.md section 5)."""
    data = json.loads(text)
    if not isinstance(data, dict):
        raise InputError("input JSON must be a top-level object")

    unknown_top = set(data) - KNOWN_TOP_KEYS
    if unknown_top:
        raise InputError(f"unknown top-level keys: {sorted(unknown_top)}")

    # temperature / P_CO
    T = float(data.get("temperature_K", 1823.15))
    P_CO = float(data.get("P_CO_atm", 1.0))
    _validate_number(T, "temperature_K")
    _validate_number(P_CO, "P_CO_atm")
    if T <= 0:
        raise InputError("temperature_K must be positive")
    if P_CO <= 0:
        raise InputError("P_CO_atm must be positive")

    # slag composition
    slag = data.get("slag_wtpc")
    if not isinstance(slag, dict) or len(slag) == 0:
        raise InputError("slag_wtpc is required and must be a non-empty JSON object")
    unknown_ox = set(slag) - ALL_INPUT_OXIDE_KEYS
    if unknown_ox:
        raise InputError(f"unknown slag keys: {sorted(unknown_ox)}")
    for k, v in slag.items():
        _validate_number(v, f"slag_wtpc.{k}")
        if v < 0:
            raise InputError(f"slag component {k} must be non-negative, got {v}")

    n_iron = [k for k in IRON_INPUT_KEYS if k in slag]
    if len(n_iron) != 1:
        raise InputError("exactly one of FeO_total / Fe_total must be present")
    n_cr = [k for k in CHROMIUM_INPUT_KEYS if k in slag]
    if len(n_cr) != 1:
        raise InputError("exactly one of Cr2O3_total / Cr_total must be present")
    iron_key = n_iron[0]
    cr_key = n_cr[0]
    if slag[iron_key] <= 0:
        raise InputError("slag iron content must be positive (slag must contain Fe)")

    total_wt = sum(slag.values())
    warnings: list[str] = []
    if not (99.0 <= total_wt <= 101.0):
        warnings.append(f"slag wt% sum is {total_wt:.2f} (expected 100 +/- 1.0)")

    # metal composition
    metal_raw = data.get("metal_wtpc", {})
    if not isinstance(metal_raw, dict):
        raise InputError("metal_wtpc must be a JSON object")
    unknown_m = set(metal_raw) - set(METAL_INPUT_KEYS)
    if unknown_m:
        raise InputError(f"unknown metal keys: {sorted(unknown_m)}")
    metal: dict[str, float] = {}
    for k in METAL_INPUT_KEYS:
        v = metal_raw.get(k, 0.0)
        _validate_number(v, f"metal_wtpc.{k}")
        metal[k] = float(v)

    # options
    options, explicit_opts = _parse_options(data.get("options"), warnings)

    # derive / validate iron_input and chromium_input vs the keys actually present
    if "iron_input" in explicit_opts:
        if options.iron_input != iron_key:
            raise InputError(
                f"options.iron_input {options.iron_input!r} inconsistent with "
                f"present input key {iron_key!r}"
            )
    else:
        options.iron_input = iron_key
    if "chromium_input" in explicit_opts:
        if options.chromium_input != cr_key:
            raise InputError(
                f"options.chromium_input {options.chromium_input!r} inconsistent with "
                f"present input key {cr_key!r}"
            )
    else:
        options.chromium_input = cr_key

    return SlagModelConfig(
        temperature_K=T,
        P_CO_atm=P_CO,
        slag_wtpc=dict(slag),
        metal_wtpc=metal,
        options=options,
        warnings=warnings,
    )


def load_input(path: str | Path) -> SlagModelConfig:
    """Load and validate a JSON input file."""
    p = Path(path)
    if not p.exists():
        raise InputError(f"input file not found: {p}")
    text = p.read_text(encoding="utf-8")
    return parse_input(text)


def build_report(cfg: SlagModelConfig, result) -> dict:
    """Build the v2 report by translating the converged result without recalculation."""
    from .data import CATION_TO_OXIDE, CATION_TO_RS_SPECIES

    comps = []
    for cation, oxide in CATION_TO_OXIDE.items():
        if cation not in result.X:
            continue
        activities = result.activities
        if oxide in activities.a_conventional_by_species:
            conventional = {
                "species": oxide,
                "standard_state": activities.standard_state_by_species[oxide],
                "rs_units_per_species": 2 if oxide in {"P2O5", "Al2O3"} else 1,
                "a": activities.a_conventional_by_species[oxide],
                "DeltaG_conversion_J_per_mol_species": (
                    activities.delta_g_conversion_J_per_mol_species[oxide]
                ),
            }
            if oxide in activities.gamma_conventional_by_species:
                conventional["gamma"] = activities.gamma_conventional_by_species[oxide]
                conventional["gamma_fraction_basis"] = "X_cation"
            unavailable_reason = None
        else:
            conventional = None
            unavailable_reason = "no documented conversion"
        comps.append(
            {
                "cation": cation,
                "rs_species": CATION_TO_RS_SPECIES[cation],
                "X_cation": result.X[cation],
                "RTln_gamma_RS_J_per_mol_cation": result.rtln_gamma_rs[cation],
                "gamma_RS": activities.gamma_rs_by_cation[cation],
                "a_RS": activities.a_rs_by_cation[cation],
                "conventional": conventional,
                "conventional_unavailable_reason": unavailable_reason,
            }
        )
    if "Ti4+" in result.X:
        comps.append(
            {
                "cation": "Ti4+",
                "rs_species": None,
                "X_cation": result.X["Ti4+"],
                "model_status": "denominator_only_legacy",
            }
        )

    eq = {
        "K_FeO": result.k_FeO,
        "K_CO": result.k_CO,
        "K_Cr": result.k_Cr,
        "Q_over_K_FeO": result.q_over_k.get("FeO"),
        "Q_over_K_CO": result.q_over_k.get("CO"),
        "Q_over_K_Cr": result.q_over_k.get("Cr"),
    }

    return {
        "schema_version": 2,
        "input": {
            "temperature_K": cfg.temperature_K,
            "P_CO_atm": cfg.P_CO_atm,
            "slag_wtpc": cfg.slag_wtpc,
            "metal_wtpc": cfg.metal_wtpc,
            "options": cfg.options.as_dict(),
        },
        "solution": {
            "P_O2_atm": result.P_O2_atm,
            "log10_P_O2": math.log10(result.P_O2_atm),
            "r_Fe": result.r_Fe,
            "r_Cr": result.r_Cr,
            "iterations": result.iterations,
            "split_wtpc": result.split.wt_split,
            "components": comps,
        },
        "metal": {
            "C_wtpc": result.C_wtpc,
            "O_wtpc": result.O_wtpc,
            "f_C": result.f_metal["C"],
            "f_Cr": result.f_metal["Cr"],
            "f_Mn": result.f_metal["Mn"],
            "f_P": result.f_metal["P"],
            "a_C": result.a_C,
            "a_Cr": result.a_Cr,
        },
        "equilibrium": eq,
        "warnings": cfg.warnings,
    }
