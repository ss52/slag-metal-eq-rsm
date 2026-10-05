"""Tests for data constants: alpha matrix, molar masses (PLAN.md test 1)."""

from __future__ import annotations

import numpy as np
import pytest

from slag_model import data
from slag_model.constants import (
    CATIONS_PER_FORMULA,
    MOLAR_MASS_ELEMENT,
    MOLAR_MASS_OXIDE,
)
from slag_model.data import ALPHA, CATION_ORDER, CATION_TO_OXIDE, EIJ, METAL_SPECIES


def test_alpha_symmetric():
    assert ALPHA.shape == (10, 10)
    assert np.allclose(ALPHA, ALPHA.T)


def test_alpha_zero_diagonal():
    assert np.all(np.diag(ALPHA) == 0.0)


def test_alpha_fe_al_typofix():
    """Verify alpha(Fe2+,Al3+) = alpha(Al3+,Fe2+) = -41000 (not -410)."""
    i = CATION_ORDER.index("Fe2+")
    j = CATION_ORDER.index("Al3+")
    assert ALPHA[i, j] == -41000.0
    assert ALPHA[j, i] == -41000.0


def test_molar_masses_positive():
    assert all(v > 0.0 for v in MOLAR_MASS_OXIDE.values())
    assert all(v > 0.0 for v in MOLAR_MASS_ELEMENT.values())


def test_cations_per_formula_positive():
    assert all(v > 0 for v in CATIONS_PER_FORMULA.values())


def test_cation_order_length():
    assert len(CATION_ORDER) == 10


def test_eij_shape():
    assert EIJ.shape == (4, 4)
    assert len(METAL_SPECIES) == 4


def test_cation_to_oxide_covers_order():
    for cation in CATION_ORDER:
        assert cation in CATION_TO_OXIDE


def test_named_parameter_profiles_have_source_exact_values_and_immutable_matrices():
    assert hasattr(data, "DEFAULT_PARAMETER_PROFILE")
    assert data.DEFAULT_PARAMETER_PROFILE == "banya93_xiao95_v1"
    assert set(data.PARAMETER_PROFILES) == {
        "banya93_xiao95_v1",
        "banya93_casi_xiao95cr_hybrid_v1",
        "legacy_workbook_typo_hybrid_v1",
    }
    i_fe3, i_ca = CATION_ORDER.index("Fe3+"), CATION_ORDER.index("Ca2+")
    i_si = CATION_ORDER.index("Si4+")
    assert data.PARAMETER_PROFILES[data.DEFAULT_PARAMETER_PROFILE][i_fe3, i_ca] == -95810.0
    assert data.PARAMETER_PROFILES[data.DEFAULT_PARAMETER_PROFILE][i_ca, i_si] == -139100.0
    assert data.PARAMETER_PROFILES["banya93_casi_xiao95cr_hybrid_v1"][i_fe3, i_ca] == -95810.0
    assert data.PARAMETER_PROFILES["banya93_casi_xiao95cr_hybrid_v1"][i_ca, i_si] == -133890.0
    assert data.PARAMETER_PROFILES["legacy_workbook_typo_hybrid_v1"][i_fe3, i_ca] == -96810.0
    assert data.PARAMETER_PROFILES["legacy_workbook_typo_hybrid_v1"][i_ca, i_si] == -133890.0
    with pytest.raises(ValueError):
        data.PARAMETER_PROFILES[data.DEFAULT_PARAMETER_PROFILE][i_fe3, i_ca] = 0.0
