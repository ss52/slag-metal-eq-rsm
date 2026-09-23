"""CLI output safety tests."""

from __future__ import annotations

from pathlib import Path

from slag_model import cli

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
