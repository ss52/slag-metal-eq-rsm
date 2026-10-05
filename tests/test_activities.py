"""Formula-unit activity relationships for slag components."""

from __future__ import annotations

import math

import pytest

import slag_model.rsm as rsm
from slag_model.constants import R
from slag_model.data import CATION_TO_RS_SPECIES


def build_slag_activities(*args, **kwargs):
    """Call the new builder through the current RSM module."""
    return rsm.build_slag_activities(*args, **kwargs)


def test_banya_p2o5_formula_unit_activity_relation():
    T = 1873.0
    x = {"P5+": 0.02}
    rtln = {"P5+": 4000.0}
    result = build_slag_activities(rtln, x, T)
    a_rs = 0.02 * math.exp(4000.0 / (R * T))
    dg = 52720.0 - 230.706 * T
    assert result.a_rs_by_cation["P5+"] == pytest.approx(a_rs, rel=1e-12)
    assert result.a_conventional_by_species["P2O5"] == pytest.approx(
        a_rs**2 * math.exp(dg / (R * T)), rel=1e-12
    )
    assert "P2O5" not in result.gamma_conventional_by_species


def test_zero_phosphorus_and_unconverted_alumina():
    result = build_slag_activities({"P5+": 0.0, "Al3+": 0.0}, {"P5+": 0.0, "Al3+": 0.1}, 1873.0)
    assert result.a_conventional_by_species["P2O5"] == 0.0
    assert "Al2O3" not in result.a_conventional_by_species
    assert result.a_rs_by_cation["Al3+"] == pytest.approx(0.1)


def test_one_cation_iron_conversion_matches_existing_coefficient():
    T = 1873.0
    x = {"Fe2+": 0.2}
    rtln = {"Fe2+": -1200.0}
    result = build_slag_activities(rtln, x, T)
    a_rs = 0.2 * math.exp(-1200.0 / (R * T))
    dg = -8540.0 + 7.142 * T
    assert result.a_conventional_by_species["FeO"] == pytest.approx(
        a_rs * math.exp(dg / (R * T)), rel=1e-12
    )
    assert result.gamma_conventional_by_species["FeO"] == pytest.approx(
        math.exp(-1200.0 / (R * T)) * math.exp(dg / (R * T)), rel=1e-12
    )
    assert result.standard_state_by_species["FeO"] == (
        "Ban-ya (1993) Table 3: Fe_tO(l) equilibrated with Fe"
    )


def test_custom_alumina_conversion_uses_two_rs_units():
    T = 1873.0
    spec = rsm.ConversionSpec(A=100.0, B=-0.1, standard_state="Al2O3(custom-solid)")
    result = build_slag_activities({"Al3+": 0.0}, {"Al3+": 0.1}, T, al2o3_conversion=spec)
    assert result.a_conventional_by_species["Al2O3"] == pytest.approx(
        0.1**2 * math.exp((100.0 - 0.1 * T) / (R * T))
    )
    assert "Al2O3" not in result.gamma_conventional_by_species


def test_custom_iron_oxide_conversion_uses_one_rs_unit():
    T = 1873.0
    spec = rsm.ConversionSpec(A=250.0, B=-0.2, standard_state="FeO1.5(custom-solid)")
    result = build_slag_activities({"Fe3+": 0.0}, {"Fe3+": 0.15}, T, fe2o3_conversion=spec)
    assert result.a_conventional_by_species["FeO1.5"] == pytest.approx(
        0.15 * math.exp((250.0 - 0.2 * T) / (R * T))
    )
    assert result.gamma_conventional_by_species["FeO1.5"] == pytest.approx(
        math.exp((250.0 - 0.2 * T) / (R * T))
    )
    assert result.standard_state_by_species["FeO1.5"] == "FeO1.5(custom-solid)"


def test_source_qualified_states_preserve_paper_phases_and_rs_species():
    cations = {
        "Fe2+": 0.1,
        "Fe3+": 0.1,
        "Ca2+": 0.1,
        "Mg2+": 0.1,
        "Mn2+": 0.1,
        "Si4+": 0.1,
        "Al3+": 0.1,
        "P5+": 0.1,
        "Cr3+": 0.1,
        "Cr2+": 0.1,
    }
    result = build_slag_activities(
        dict.fromkeys(cations, 0.0), cations, 1873.0, sio2_conversion="banya"
    )
    assert CATION_TO_RS_SPECIES["Al3+"] == "AlO1.5"
    assert CATION_TO_RS_SPECIES["P5+"] == "PO2.5"
    assert result.standard_state_by_species["SiO2"] == (
        "SiO2(beta-cristobalite), Ban-ya (1993) Table 3"
    )
    assert result.standard_state_by_species["P2O5"] == ("P2O5(l), Ban-ya (1993) Table 3")
    assert result.standard_state_by_species["FeO"] == (
        "Ban-ya (1993) Table 3: Fe_tO(l) equilibrated with Fe"
    )
    assert result.standard_state_by_species["CaO"] == ("CaO(s), Ban-ya (1993) Table 3")
    assert result.standard_state_by_species["MgO"] == ("MgO(s), Ban-ya (1993) Table 3")
    assert result.standard_state_by_species["MnO"] == ("MnO(s), Ban-ya (1993) Table 3")
    assert result.standard_state_by_species["CrO"] == ("CrO(liquid), Xiao & Holappa (1995) Table 3")
    assert result.standard_state_by_species["CrO1.5"] == (
        "CrO1.5(solid), Xiao & Holappa (1995) Table 3"
    )
    workbook_silica = build_slag_activities(
        {"Si4+": 0.0}, {"Si4+": 0.1}, 1873.0, sio2_conversion="workbook"
    )
    assert workbook_silica.standard_state_by_species["SiO2"] == (
        "SiO2(solid), polymorph unspecified; Xiao & Holappa (1995), silica-in-lime-silica section"
    )
