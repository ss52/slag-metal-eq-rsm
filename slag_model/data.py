"""Thermodynamic data (PLAN.md section 6): alpha matrix, conversions, e_ij, DeltaG."""

from __future__ import annotations

from types import MappingProxyType

import numpy as np

# Cation ordering used throughout the regular solution model
CATION_ORDER: list[str] = [
    "Fe2+",
    "Fe3+",
    "Ca2+",
    "Mg2+",
    "Mn2+",
    "Si4+",
    "Al3+",
    "P5+",
    "Cr3+",
    "Cr2+",
]

# Cation interaction energy alpha_ij (J/mol), symmetric, zero diagonal.
# The default uses Ban-ya (1993) Table 1 base pairs (including the corrected
# Fe3+-Ca2+ value) and Xiao & Holappa (1995) Table 2 values where selected,
# including Ca2+-Si4+ and the chromium-bearing terms.
# Order: [Fe2+, Fe3+, Ca2+, Mg2+, Mn2+, Si4+, Al3+, P5+, Cr3+, Cr2+]
ALPHA: np.ndarray = np.array(
    [
        [
            0.0,
            -18660.0,
            -31380.0,
            33470.0,
            7110.0,
            -41840.0,
            -41000.0,
            -31380.0,
            0.0,
            0.0,
        ],
        [
            -18660.0,
            0.0,
            -95810.0,
            -2930.0,
            -56480.0,
            32640.0,
            -161080.0,
            14640.0,
            0.0,
            0.0,
        ],
        [
            -31380.0,
            -95810.0,
            0.0,
            -100420.0,
            -92050.0,
            -139100.0,
            -154810.0,
            -251040.0,
            44235.0,
            -6140.0,
        ],
        [
            33470.0,
            -2930.0,
            -100420.0,
            0.0,
            61920.0,
            -66940.0,
            -71130.0,
            -37660.0,
            28085.0,
            5520.0,
        ],
        [
            7110.0,
            -56480.0,
            -92050.0,
            61920.0,
            0.0,
            -75310.0,
            -83680.0,
            -84940.0,
            0.0,
            0.0,
        ],
        [
            -41840.0,
            32640.0,
            -139100.0,
            -66940.0,
            -75310.0,
            0.0,
            -127610.0,
            83680.0,
            -48975.0,
            -60540.0,
        ],
        [
            -41000.0,
            -161080.0,
            -154810.0,
            -71130.0,
            -83680.0,
            -127610.0,
            0.0,
            -261500.0,
            -46210.0,
            -30285.0,
        ],
        [
            -31380.0,
            14640.0,
            -251040.0,
            -37660.0,
            -84940.0,
            83680.0,
            -261500.0,
            0.0,
            0.0,
            0.0,
        ],
        [0.0, 0.0, 44235.0, 28085.0, 0.0, -48975.0, -46210.0, 0.0, 0.0, 32710.0],
        [0.0, 0.0, -6140.0, 5520.0, 0.0, -60540.0, -30285.0, 0.0, 32710.0, 0.0],
    ],
    dtype=float,
)

DEFAULT_PARAMETER_PROFILE = "banya93_xiao95_v1"
BANYA_CASI_PARAMETER_PROFILE = "banya93_casi_xiao95cr_hybrid_v1"
LEGACY_WORKBOOK_PARAMETER_PROFILE = "legacy_workbook_typo_hybrid_v1"

_banya_casi_matrix = ALPHA.copy()
_i_ca = CATION_ORDER.index("Ca2+")
_i_si = CATION_ORDER.index("Si4+")
_banya_casi_matrix[_i_ca, _i_si] = -133890.0
_banya_casi_matrix[_i_si, _i_ca] = -133890.0

_legacy_workbook_matrix = _banya_casi_matrix.copy()
_i_fe3 = CATION_ORDER.index("Fe3+")
_legacy_workbook_matrix[_i_fe3, _i_ca] = -96810.0
_legacy_workbook_matrix[_i_ca, _i_fe3] = -96810.0

for _matrix in (ALPHA, _banya_casi_matrix, _legacy_workbook_matrix):
    _matrix.setflags(write=False)

PARAMETER_PROFILES: dict[str, np.ndarray] = MappingProxyType(
    {
        DEFAULT_PARAMETER_PROFILE: ALPHA,
        BANYA_CASI_PARAMETER_PROFILE: _banya_casi_matrix,
        LEGACY_WORKBOOK_PARAMETER_PROFILE: _legacy_workbook_matrix,
    }
)
PARAMETER_PROFILE_METADATA: dict[str, dict[str, str | int]] = {
    DEFAULT_PARAMETER_PROFILE: {
        "id": DEFAULT_PARAMETER_PROFILE,
        "version": 1,
        "name": "Ban-ya 1993 base with Xiao-Holappa 1995 selected parameters",
        "source_basis": (
            "Ban-ya 1993 Table 1 major-pair data with corrected Fe3+-Ca2+; "
            "Xiao & Holappa 1995 Table 2 Ca2+-Si4+ and chromium-bearing terms."
        ),
    },
    BANYA_CASI_PARAMETER_PROFILE: {
        "id": BANYA_CASI_PARAMETER_PROFILE,
        "version": 1,
        "name": "Ban-ya Ca-Si comparison with Xiao-Holappa chromium terms (hybrid)",
        "source_basis": (
            "Ban-ya 1993 Table 1 Ca2+-Si4+ and corrected Fe3+-Ca2+; "
            "Xiao & Holappa 1995 Table 2 chromium-bearing terms. This is not a pure Ban-ya set."
        ),
    },
    LEGACY_WORKBOOK_PARAMETER_PROFILE: {
        "id": LEGACY_WORKBOOK_PARAMETER_PROFILE,
        "version": 1,
        "name": "Historical workbook interaction snapshot (hybrid, with Fe3+-Ca2+ typo)",
        "source_basis": (
            "Historical regression only: Ban-ya Ca2+-Si4+, Xiao & Holappa chromium terms, "
            "and the old -96810 J/mol Fe3+-Ca2+ transcription. Not a paper source set."
        ),
    },
}

# Backward-compatible module-level name points only to the scientific default.
# Parameter selection at runtime uses the immutable profile mapping above.
ALPHA = PARAMETER_PROFILES[DEFAULT_PARAMETER_PROFILE]

# Mapping between a cation and the oxide label used for its activity / standard state.
CATION_TO_OXIDE: dict[str, str] = {
    "Fe2+": "FeO",
    "Fe3+": "FeO1.5",
    "Ca2+": "CaO",
    "Mg2+": "MgO",
    "Mn2+": "MnO",
    "Si4+": "SiO2",
    "Al3+": "Al2O3",
    "P5+": "P2O5",
    "Cr3+": "CrO1.5",
    "Cr2+": "CrO",
}
OXIDE_TO_CATION: dict[str, str] = {v: k for k, v in CATION_TO_OXIDE.items()}

# R.S. pseudo-species are defined per cation and are distinct from the
# conventional oxide formula units used by CATION_TO_OXIDE.
CATION_TO_RS_SPECIES: dict[str, str] = {
    "Fe2+": "FeO",
    "Fe3+": "FeO1.5",
    "Ca2+": "CaO",
    "Mg2+": "MgO",
    "Mn2+": "MnO",
    "Si4+": "SiO2",
    "Al3+": "AlO1.5",
    "P5+": "PO2.5",
    "Cr3+": "CrO1.5",
    "Cr2+": "CrO",
}

# Input composition keys (analytical slag form)
INPUT_OXIDE_KEYS: list[str] = ["SiO2", "CaO", "MgO", "Al2O3", "MnO", "TiO2", "P2O5"]
IRON_INPUT_KEYS: list[str] = ["FeO_total", "Fe_total"]
CHROMIUM_INPUT_KEYS: list[str] = ["Cr2O3_total", "Cr_total"]
METAL_INPUT_KEYS: list[str] = ["Cr", "Mn", "P"]

# Standard-state conversion factors: Delta G_conv = A + B*T (J/mol of the
# conventional species). P2O5 is a two-cation formula-unit conversion; it must
# not be interpreted as a conventional gamma on the cation-fraction basis.
CHROMIUM_OXIDE_CONVERSION_REACTIONS: dict[str, tuple[float, float]] = {
    # Xiao & Holappa (1995) Table 3; all energies are J/mol.
    "CrO(liq)=CrO(RS)": (77150.0, -33.5),
    "CrO1.5(s)=CrO1.5(RS)": (74967.0, -37.5),
    "CrO1.5(liq)=CrO1.5(RS)": (10152.0, -12.6),
    "CrO1.5(s)=CrO1.5(liq)": (64815.0, -24.9),
}
CONVERSION_DEFAULTS: dict[str, tuple[float, float]] = {
    "FeO": (-8540.0, 7.142),
    "CaO": (18160.0, -23.309),
    "MgO": (34350.0, -16.736),
    "MnO": (-32470.0, 26.143),
    "P2O5": (52720.0, -230.706),
    "CrO": CHROMIUM_OXIDE_CONVERSION_REACTIONS["CrO(liq)=CrO(RS)"],
    "CrO1.5": CHROMIUM_OXIDE_CONVERSION_REACTIONS["CrO1.5(s)=CrO1.5(RS)"],
}
SIO2_CONVERSIONS: dict[str, tuple[float, float]] = {
    "workbook": (51346.0, -13.88),
    "banya": (27030.0, -1.983),
}

# Conventional reference states from the cited conversion tables. Keep the
# source's wording when it does not identify a pure oxide phase; do not infer
# a phase from a species formula.
CONVERSION_STANDARD_STATES: dict[str, str] = {
    "FeO": "Ban-ya (1993) Table 3: Fe_tO(l) equilibrated with Fe",
    "CaO": "CaO(s), Ban-ya (1993) Table 3",
    "MgO": "MgO(s), Ban-ya (1993) Table 3",
    "MnO": "MnO(s), Ban-ya (1993) Table 3",
    "P2O5": "P2O5(l), Ban-ya (1993) Table 3",
    "CrO": "CrO(liquid), Xiao & Holappa (1995) Table 3",
    "CrO1.5": "CrO1.5(solid), Xiao & Holappa (1995) Table 3",
}
SIO2_CONVERSION_STANDARD_STATES: dict[str, str] = {
    "banya": "SiO2(beta-cristobalite), Ban-ya (1993) Table 3",
    "workbook": (
        "SiO2(solid), polymorph unspecified; Xiao & Holappa (1995), silica-in-lime-silica section"
    ),
}
# The explicit activity builder leaves FeO1.5 and Al2O3 unavailable unless a
# labeled custom conversion is supplied.

# Equilibrium reaction coefficients (Delta G = A + B*T, J/mol)
DELTA_G_FEO: tuple[float, float] = (-232600.0, 47.9)  # Fe(l) + 1/2 O2 = FeO(l)
DELTA_G_CO: tuple[float, float] = (-134300.0, -45.40)  # [C](1wt%) + 1/2 O2 = CO(g)
DELTA_G_CO_GRAPHITE: tuple[float, float] = (-111710.0, -87.66)  # C(gr) + 1/2 O2 = CO(g)
DELTA_G_C_DISSOLUTION: tuple[float, float] = (22590.0, -42.26)  # C(gr) = [C]_1wt%
# Chromium equilibrium source reactions.
# Xiao & Holappa, INFACON VII (1995), p. 322, Eqs. 9-10:
# 2 CrO1.5 + Cr(s) = 3 CrO, with Delta G in cal/mol reaction.
DELTA_G_CR_SOLID_CAL: tuple[float, float] = (25690.0, -13.36)
# Cr(s) = [Cr] on the Henrian 1 mass-% standard, per mol Cr. Xiao, Kou & Fang
# (2018), p. 419, Table 1, lists (19246 - 46.86*T) J/mol Cr, the rounded
# cal-to-J transcription of 4600 - 11.20*T cal/g-atom in Stefanescu & Katz,
# ASM Handbook vol. 15 (2008), p. 47, Table 3. Do not attribute it to
# Sigworth & Elliott (1974).
DELTA_G_CR_DISSOLUTION_J: tuple[float, float] = (19246.0, -46.86)

# Wagner first-order interaction parameters e_i^j (1600 degC values) on the
# Henrian 1 wt% scale. Rows/cols in METAL_SPECIES order.
# The compiled e_Mn^C value is -0.07 at 1873 K: PAN, Archives of Foundry
# Engineering 24(2) (2024), Table 1 p. 112, citing Chen (2010), Common Charts
# and Databook for Steelmaking. This is a compilation source, not direct
# verification of the original Sigworth & Elliott table.
METAL_SPECIES: list[str] = ["C", "Cr", "Mn", "P"]
EIJ: np.ndarray = np.array(
    [
        [0.14, -0.024, -0.012, 0.051],
        [-0.12, 0.0, 0.0, -0.053],
        [-0.07, 0.0, 0.0, -0.0035],
        [0.13, -0.03, 0.0, 0.0],
    ],
    dtype=float,
)

# Fe redox parameters (Ban-ya Eq. 24):
#   log10(Fe3+/Fe2+) = 6625/T - 2.77 + 0.25*log10(P_O2)
#                       + log10(gamma_FeO) - log10(gamma_FeO1.5)
BANYA_T_COEF: float = 6625.0
BANYA_CONSTANT: float = -2.77
BANYA_PO2_COEF: float = 0.25

# Approximate dissolved-oxygen diagnostic for liquid Fe, assuming f_O = 1:
#   [%O] = (a_FeO/a_Fe) * 10**(-6320/T + 2.734)
O_SAT_T_COEF: float = -6320.0
O_SAT_CONSTANT: float = 2.734
