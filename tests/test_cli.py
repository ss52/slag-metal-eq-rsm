"""CLI output safety tests."""

from __future__ import annotations

from pathlib import Path

from slag_model import cli
from slag_model.io import build_report, load_input
from slag_model.solver import solve_redox_fixed_po2

ROOT = Path(__file__).resolve().parents[1]
REFERENCE_INPUT = ROOT / "examples" / "eaf_slag_01.json"


def test_nonfinite_report_does_not_create_output_file(tmp_path, monkeypatch, capsys):
    output = tmp_path / "report.json"
    build_report = cli.build_report

    def report_with_nan(cfg, result):
        report = build_report(cfg, result)
        report["solution"]["r_Fe"] = float("nan")
        return report

    monkeypatch.setattr(cli, "build_report", report_with_nan)

    status = cli.main([str(REFERENCE_INPUT), "--output", str(output)])

    assert status != 0
    assert not output.exists()
    assert "error:" in capsys.readouterr().err


def test_fixed_pressure_cli_labels_unavailable_oxygen_and_diagnostics(capsys):
    cfg = load_input(REFERENCE_INPUT)
    result = solve_redox_fixed_po2(cfg, P_O2=1.0e-9, C=0.018)

    cli._print_report(build_report(cfg, result))

    output = capsys.readouterr().out
    assert "calculation mode: fixed_P_O2" in output
    assert "[%O] = unavailable" in output
    assert "FeO-equivalent [%O]" in output
    assert "not enforced" in output
    assert "Q/K values are diagnostics" in output
