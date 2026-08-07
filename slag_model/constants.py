"""Physical constants and molar masses (PLAN.md section 4, 6.1)."""

from __future__ import annotations

R: float = 8.314  # J/mol/K

# Molar masses of oxides (g/mol)
MOLAR_MASS_OXIDE: dict[str, float] = {
    "SiO2": 60.084,
    "CaO": 56.077,
    "MgO": 40.304,
    "Al2O3": 101.961,
    "MnO": 70.937,
    "TiO2": 79.866,
    "P2O5": 141.943,
    "FeO": 71.846,
    "Fe2O3": 159.688,
    "CrO": 67.996,
    "Cr2O3": 151.990,
}

# Molar masses of pure elements (g/mol), for total-element input forms
MOLAR_MASS_ELEMENT: dict[str, float] = {
    "Fe": 55.847,
    "Cr": 51.996,
}

# Molar masses of steel solutes (g/mol)
MOLAR_MASS_METAL: dict[str, float] = {
    "Fe": 55.847,
    "Cr": 51.996,
    "Mn": 54.938,
    "P": 30.974,
    "C": 12.011,
}

# Cations per formula unit of each oxide
CATIONS_PER_FORMULA: dict[str, int] = {
    "SiO2": 1,
    "CaO": 1,
    "MgO": 1,
    "Al2O3": 2,
    "MnO": 1,
    "TiO2": 1,
    "P2O5": 2,
    "FeO": 1,
    "Fe2O3": 2,
    "CrO": 1,
    "Cr2O3": 2,
}
