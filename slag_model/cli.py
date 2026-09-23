"""CLI entry point (PLAN.md section 3, 9)."""

from __future__ import annotations

import json
import sys

from .data import CATION_TO_OXIDE
from .io import InputError, build_report, load_input
from .solver import NonConvergenceError, solve


def _print_report(report: dict) -> None:
    sol = report["solution"]
    print("=" * 78)
    print("Slag Activity Coefficient Model  (Ban-ya RSM + oxygen-potential balance)")
    print("=" * 78)

    inp = report["input"]
    print(f"  Temperature : {inp['temperature_K']:.2f} K")
    print(f"  P_CO        : {inp['P_CO_atm']} atm")
    print()
    print("  Slag input (wt%):")
    for k, v in inp["slag_wtpc"].items():
        if v != 0.0 or k not in ("P2O5",):
            print(f"    {k:12s} {v:10.4f}")
    print("  Metal input (wt%):")
    for k, v in inp["metal_wtpc"].items():
        print(f"    {k:12s} {v:10.4f}")
    print()
    print(f"  P_O2 = {sol['P_O2_atm']:.6e} atm  (log10 = {sol['log10_P_O2']:.4f})")
    print(
        f"  r_Fe = {sol['r_Fe']:.6f}   r_Cr = {sol['r_Cr']:.6f}   iterations = {sol['iterations']}"
    )
    print()

    print("  Split slag composition (wt%):")
    for k, v in sol["split_wtpc"].items():
        print(f"    {k:12s} {v:10.4f}")
    print()

    print("  Component activities:")
    for c in sol["components"]:
        if c.get("model_status") == "denominator_only_legacy":
            print(
                f"    {c['cation']:5s} X_cation={c['X_cation']:.6f}  "
                "denominator-only legacy; activity unavailable"
            )
            continue
        print(
            f"    {c['cation']:5s} ({c['rs_species']})  "
            f"X_cation={c['X_cation']:.6f}  "
            f"RTln(gamma_RS)={c['RTln_gamma_RS_J_per_mol_cation']:.3f} J/mol cation  "
            f"gamma_RS={c['gamma_RS']:.4f}  a_RS={c['a_RS']:.6g}"
        )
        conventional = c["conventional"]
        if conventional is None:
            print(
                f"      conventional {CATION_TO_OXIDE[c['cation']]} unavailable: "
                f"{c['conventional_unavailable_reason']}"
            )
            continue
        gamma_text = (
            f", gamma={conventional['gamma']:.6g} (X_cation basis)"
            if "gamma" in conventional
            else ""
        )
        print(
            f"      conventional {conventional['species']} "
            f"[{conventional['standard_state']}; "
            f"{conventional['rs_units_per_species']} R.S. unit(s)]: "
            f"a={conventional['a']:.6g}{gamma_text}"
        )
    print()

    m = report["metal"]
    print("  Metal:")
    print(f"    [%C] = {m['C_wtpc']:.4f}    [%O] = {m['O_wtpc']:.4f}")
    print(
        f"    f_C = {m['f_C']:.5f}    f_Cr = {m['f_Cr']:.5f}"
        f"    f_Mn = {m['f_Mn']:.5f}    f_P = {m['f_P']:.5f}"
    )
    print(f"    a_C = {m['a_C']:.6f}    a_Cr = {m['a_Cr']:.6f}")
    print()

    eq = report["equilibrium"]
    print("  Equilibrium:")
    print(f"    K_FeO = {eq['K_FeO']:.4e}    K_CO = {eq['K_CO']:.4e}    K_Cr = {eq['K_Cr']:.4f}")
    print(
        f"    Q/K :  Fe-O = {eq['Q_over_K_FeO']:.9f}"
        f"    C-O = {eq['Q_over_K_CO']:.9f}"
        f"    Cr  = {eq['Q_over_K_Cr']}"
    )

    for w in report["warnings"]:
        print(f"  WARNING: {w}")
    print()


def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        prog="slag_model",
        description="EAF slag activity coefficient model (Ban-ya RSM)",
    )
    parser.add_argument("input", help="path to JSON input file")
    parser.add_argument(
        "-o", "--output", default=None, help="path to write JSON report (default: none)"
    )
    args = parser.parse_args(argv)

    try:
        cfg = load_input(args.input)
        for w in cfg.warnings:
            print(f"warning: {w}", file=sys.stderr)
        result = solve(cfg)
    except (InputError, NonConvergenceError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    report = build_report(cfg, result)
    _print_report(report)

    if args.output:
        with open(args.output, "w", encoding="utf-8") as fh:
            json.dump(report, fh, indent=2)
        print(f"report written to {args.output}")

    return 0
