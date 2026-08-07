# Slag Activity Model

Thermodynamic model for **electric-arc-furnace (EAF) slag–metal equilibrium**.

The package computes activity coefficients and activities of all slag oxide
components at a given temperature using **Ban-ya's quadratic formalism**
(a cation regular-solution model, RSM), coupled to the steel bath through a
**single shared oxygen potential**. Iron and chromium are entered analytically
as one total each and are split automatically into their valence states
(Fe²⁺/Fe³⁺ and Cr²⁺/Cr³⁺) by the redox equilibria.

[![Tests](https://img.shields.io/badge/tests-56%20passed-brightgreen)](#tests)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Python](https://img.shields.io/badge/python-%3E%3D3.13-blue)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/license-MIT-yellow)](LICENSE)

---

## Physical basis

The model closes the slag–metal system with **one** oxygen potential $P_{O_2}$:

1. **Slag side (RSM).** The slag is treated as a regular solution of cations
   (Ban-ya 1993). Oxide-component activity coefficients follow from cation
   fractions $X_i$ and pair interaction energies $\alpha_{ij}$:

$$RT\ln\gamma_i^{RS} = \sum_j \alpha_{ij} X_j^2 + \sum_{j} \sum_{k} (\alpha_{ij} + \alpha_{ik} - \alpha_{jk}) X_j X_k$$

   A linear conversion $\Delta G_{conv} = A + BT$ shifts the raw
   regular-solution coefficients onto the **pure-stable-oxide (Raoultian)
   standard state** used by thermodynamic software such as HSC Chemistry.

2. **Oxygen potential.** $P_{O_2}$ is *not* an input — the FeO-rich slag
   buffers it through the Fe–O equilibrium:

$$P_{O_2} = \left( \frac{a_{FeO}}{K_{FeO}\, a_{Fe}} \right)^2$$

3. **Iron redox.** The Fe³⁺/Fe²⁺ ratio is slaved to $P_{O_2}$
   (Ban-ya Eq. 24).

4. **Chromium redox.** The Cr²⁺/Cr³⁺ ratio is pinned by the chromium activity
   in the steel via $2\,\underline{CrO_{1.5}} + [Cr] = 3\,\underline{CrO}$
   (Xiao & Holappa) — oxygen cancels, so no $P_{O_2}$ is needed. The model
   solves the **exact cubic** relation rather than the approximate square-root
   fixed point found in spreadsheet implementations.

5. **Metal side (WIPF).** Solute activity coefficients come from **Wagner's
   first-order interaction-parameter formalism** (coefficients from
   Sigworth & Elliott). The carbon content that is in equilibrium with the slag
   follows from the C–O reaction at the same $P_{O_2}$ — slag and steel are
   in balance *by construction*.

Because the oxygen potential is read off the slag and the balanced carbon is then computed from it, the coupled system converges in typically **5–30 damped fixed-point iterations**.

---

## Installation

Requires **Python ≥ 3.13** and [uv](https://docs.astral.sh/uv/).

```bash
git clone <repository-url>
cd "2026 TD calc project"
uv sync
```

`uv sync` creates a virtual environment and installs the runtime dependency
(`numpy`) plus the dev tools (`pytest`, `ruff`).

## Usage

Run the model on the reference EAF slag and write a JSON report:

```bash
uv run python -m slag_model examples/eaf_slag_01.json -o report.json
```

or via the installed console script:

```bash
uv run slag-model examples/eaf_slag_01.json -o report.json
```

The program prints a formatted report to the terminal and, with `-o`, writes
the full machine-readable report (input echo, redox splits, per-component
$\gamma$ and $a$, metal balance, equilibrium diagnostics $Q/K$) to JSON.

## Input JSON reference

Top-level object. All keys are optional unless noted.

| Key            | Type   | Default   | Description |
|----------------|--------|-----------|-------------|
| `temperature_K`| number | `1823.15` | Bath temperature in kelvin. Must be > 0. |
| `P_CO_atm`     | number | `1.0`     | CO partial pressure in atm for the C–O balance. Must be > 0. |
| `slag_wtpc`    | object | —         | **Required.** Slag composition in wt% (analytical form). See below. |
| `metal_wtpc`   | object | `{}`      | Steel solute contents in wt%. Missing keys default to 0.0. |
| `options`      | object | `{}`      | Model configuration flags (see the Options table). |

### `slag_wtpc` fields

| Key           | Description |
|---------------|-------------|
| `SiO2`        | wt% silica. |
| `CaO`         | wt% lime. |
| `MgO`         | wt% magnesia. |
| `Al2O3`       | wt% alumina. |
| `MnO`         | wt% manganese oxide. |
| `TiO2`        | wt% titania. (Ti⁴⁺ is *not* an RSM cation — see `ti_handling`.) |
| `P2O5`        | wt% phosphate. |
| `FeO_total`   | **Iron total as FeO** (wt%). Mutually exclusive with `Fe_total`. Exactly one must be present and > 0. |
| `Fe_total`    | **Iron total as elemental Fe** (wt%). Mutually exclusive with `FeO_total`. |
| `Cr2O3_total` | **Chromium total as Cr₂O₃** (wt%). Mutually exclusive with `Cr_total`. Exactly one must be present. |
| `Cr_total`    | **Chromium total as elemental Cr** (wt%). Mutually exclusive with `Cr2O3_total`. |

Iron and chromium are entered as *single analytical totals*; the model splits
them into Fe²⁺/Fe³⁺ and Cr²⁺/Cr³⁺ during the solve.

### `metal_wtpc` fields

| Key  | Description |
|------|-------------|
| `Cr` | wt% chromium dissolved in the steel. |
| `Mn` | wt% manganese. |
| `P`  | wt% phosphorus. |

(Steel carbon is **not** an input — it is the *output* of the C–O balance.)

### Options

| Option             | Values                                   | Default                | Meaning |
|--------------------|------------------------------------------|------------------------|---------|
| `cross_terms`      | `full` / `major5`                        | `full`                 | RSM cross-term set. `major5` exists **only** to reproduce the reference workbook in tests. |
| `ti_handling`      | `exclude_renormalize` / `as_excel`       | `exclude_renormalize`  | `as_excel` counts Ti⁴⁺ in the cation sum with zero α (reproduces the workbook); default drops TiO₂ and renormalises (~0.5 % difference). |
| `sio2_conversion`  | `workbook` / `banya`                     | `workbook`             | SiO₂ standard-state conversion. The 2.6 kJ difference moves γ(SiO₂) by ~1.19×. |
| `fe2o3_conversion` | `none` / `{ "A": ..., "B": ... }`        | `none`                 | Optional custom ΔG conversion for FeO₁.₅ (`ΔG = A + B·T`, J/mol). |
| `al2o3_conversion` | `none` / `{ "A": ..., "B": ... }`        | `none`                 | **Data gap** — see below. |
| `a_Fe`             | `unity` / `xfe`                          | `unity`                | Iron activity: 1.0 (pure liquid Fe) or the metal mole fraction. |
| `iron_input`       | `FeO_total` / `Fe_total`                 | `FeO_total`            | Iron input form (must match the key used in `slag_wtpc`). |
| `chromium_input`   | `Cr2O3_total` / `Cr_total`               | `Cr2O3_total`          | Chromium input form (must match the key used in `slag_wtpc`). |

### Validation rules

* Exactly one of `FeO_total` / `Fe_total` and one of `Cr2O3_total` / `Cr_total`
  must be present; unknown keys, negative values, or an inconsistent
  `iron_input`/`chromium_input` flag raise an `InputError`.
* Slag wt% sum should be 100 ± 1, otherwise a warning is emitted. (After the
  redox split the oxide mass can drift slightly from 100 g because Fe₂O₃/Cr₂O₃
  carry more oxygen than FeO/CrO — this is expected, not an error.)
* Missing `metal_wtpc` components default to 0.0.

## Output JSON reference

The report written with `-o` has five top-level blocks.

### `input`

An exact echo of the resolved inputs: `temperature_K`, `P_CO_atm`,
`slag_wtpc`, `metal_wtpc`, and `options` (with all defaults filled in and the
custom conversions shown as `none` when unset).

### `solution`

The converged slag state.

| Key            | Description |
|----------------|-------------|
| `P_O2_atm`     | Equilibrium oxygen potential set by the slag (atm). |
| `log10_P_O2`   | `log10(P_O2_atm)`. |
| `r_Fe`         | Converged ratio $r_{Fe} = n_{Fe^{3+}} / n_{Fe^{2+}}$. |
| `r_Cr`         | Converged ratio $r_{Cr} = n_{Cr^{2+}} / n_{Cr^{3+}}$. |
| `iterations`   | Number of damped fixed-point iterations to convergence. |
| `split_wtpc`   | Slag composition (wt%) **after** redox splitting: `FeO`, `Fe2O3`, `CrO`, `Cr2O3` plus the un-split oxides. |
| `components`   | Array of per-oxide rows (see below). |

Each entry in `components`:

| Field            | Description |
|------------------|-------------|
| `oxide`          | Oxide label (e.g. `FeO`, `SiO2`, `CrO1.5`). |
| `X`              | Cation mole fraction $X_i$ of the corresponding cation. |
| `RTln_gamma_RS`  | $RT\ln\gamma_i^{RS}$ on the regular-solution scale (J/mol). |
| `DeltaG_conv`    | Standard-state conversion $\Delta G_{conv} = A + BT$ (J/mol). |
| `gamma`          | Final activity coefficient $\gamma_i$ (pure-stable-oxide scale). |
| `a`              | Activity $a_i = \gamma_i X_i$. |

### `metal`

The steel bath in equilibrium with the slag.

| Field    | Description |
|----------|-------------|
| `C_wtpc` | Balanced dissolved carbon $[\%C]$ from the C–O equilibrium (wt%). |
| `O_wtpc` | Dissolved oxygen $[\%O]$ diagnostic $= a_{FeO}\cdot[\%O]_{sat}$ (wt%). |
| `f_C`    | Henrian activity coefficient of carbon (WIPF). |
| `f_Cr`   | Henrian activity coefficient of chromium. |
| `f_Mn`   | Henrian activity coefficient of manganese. |
| `f_P`    | Henrian activity coefficient of phosphorus. |
| `a_C`    | Carbon activity $a_C = f_C[\%C]$ (Henrian 1 wt% scale). |
| `a_Cr`   | Chromium activity $a_{Cr} = f_{Cr}[\%Cr]$. |

### `equilibrium`

| Field         | Description |
|---------------|-------------|
| `K_FeO`       | Equilibrium constant of $Fe(l) + \tfrac12 O_2 = FeO(l)$ at $T$. |
| `K_CO`        | Equilibrium constant of $[C] + \tfrac12 O_2 = CO(g)$ (Henrian C). |
| `K_Cr`        | Equilibrium constant of $2\,CrO_{1.5} + [Cr] = 3\,CrO$ at $T$. |
| `Q_over_K_FeO`| Reaction-quotient / $K$ check for Fe–O (equals 1.0 by construction). |
| `Q_over_K_CO` | Reaction-quotient / $K$ check for C–O (equals 1.0 by construction). |
| `Q_over_K_Cr` | Reaction-quotient / $K$ check for Cr redox (`null` if no Cr present). |

### `warnings`

Array of human-readable warning strings emitted during parsing (e.g. slag wt%
sum out of range, the `Al2O3` data-gap notice). Empty when there are no
warnings.

---

## Project layout

```
slag_model/
  constants.py     # R, molar masses, cations per formula unit
  data.py          # alpha matrix, conversion factors, WIPF e_ij, Delta G coeffs
  rsm.py           # split, cation fractions, RSM gammas, standard-state conversion
  redox.py         # Fe split (Ban-ya Eq. 24), Cr split (exact cubic)
  metal.py         # WIPF, activities, balanced carbon, dissolved oxygen
  equilibrium.py   # K_FeO, K_CO, K_Cr, P_O2 closure
  solver.py        # damped fixed-point driver
  io.py            # JSON input parsing / validation, report builder
  cli.py           # entry point: python -m slag_model input.json [-o report.json]
examples/
  eaf_slag_01.json # reference EAF slag
tests/             # pytest suite (56 tests)
PLAN.md            # full technical specification (single source of truth)
```

The authoritative, self-contained specification — including every equation,
the full $\alpha_{ij}$ matrix, all conversion factors, and the acceptance
criteria — is in **[PLAN.md](PLAN.md)**.

## Tests

```bash
uv run pytest -v
```

The suite covers data integrity, a **workbook regression** (RSM `RT·ln γ` per cation, converted
γ, redox ratios), the WIPF regression, equilibrium-constant ranges including
the Henrian two-ΔG decomposition of $K_{CO}$ (the standard-state trap), and
the balanced solve (mass balance, $Q/K = 1$, physical ranges).

## Code quality

Linted and formatted with [Ruff](https://github.com/astral-sh/ruff):

```bash
uv run ruff check .
uv run ruff format --check .
```

## Known data gaps

These are **explicit** — the model does not invent values:

1. **Al₂O₃ conversion factor** is TBD (default `A=0, B=0`, so γ(Al₂O₃) is on
   the regular-solution scale and a warning is emitted).
2. **FeO₁.₅ conversion** defaults to none (pure Fe₂O₃ is nearly solid at
   1823 K, $T_m \approx 1838$ K).
3. **Fe–Cr, Mn–Cr, P–Cr pair interactions are 0** — a limit of the Xiao
   chromium parameter set, which only covers Ca/Mg/Si/Al/Cr partners.
4. **WIPF $e_i^j$** are 1600 °C values used at 1550 °C; second-order terms are
   neglected (acceptable at the dilute C ≈ 0.04 %, Cr ≈ 0.04 % levels).

## References

1. **S. Ban-ya**, *Mathematical Expression of Slag-Metal Reactions in
   Steelmaking Process by Quadratic Formalism Based on the Regular Solution
   Model*, ISIJ International **33** (1993), No. 1, 2–11.
   <https://www.jstage.jst.go.jp/article/isijinternational1989/33/1/33_1_2/_article>
   — the quadratic (regular-solution) formalism, the α-matrix, the Fe-redox
   Eq. 24, and the oxide conversion factors.

2. **S. Ban-ya, F. Ishii, Y. Iguchi, T. Nagasaka**, *Regular Solution Model
   for the Equilibrium of Distribution of Oxygen between Liquid Iron and
   Slag*, Tetsu-to-Hagané **71** (1985), No. 7, 853–860.
   <https://www.jstage.jst.go.jp/article/tetsutohagane1955/71/7/71_7_853/_article>
   — the regular-solution treatment of the Fe–O equilibrium / oxygen potential
   of iron-bearing slags.

3. **Y. Xiao, L. Holappa**, *Determination of Activities in Slags Containing
   Chromium Oxides*, ISIJ International **33** (1993), No. 1, 66–74.
   <https://www.jstage.jst.go.jp/article/isijinternational1989/33/1/33_1_66/_article>
   — the chromium redox equilibrium, the CrO / CrO₁.₅ conversion factors, and
   the Cr-cation interaction parameters.

4. **Y. Xiao, L. Holappa, M.A. Reuter**, *Oxidation State and Activities of
   Chromium Oxides in CaO-SiO₂-CrOₓ Slag System*, Metallurgical and Materials
   Transactions B **33** (2002), 595–603.
   <https://link.springer.com/article/10.1007/s11663-002-0039-9>
   — experimental validation of the chromium oxidation state and activities.

5. **G.K. Sigworth, J.F. Elliott**, *The Thermodynamics of Liquid Dilute Iron
   Alloys*, Metal Science **8** (1974), No. 1, 298–310.
   <https://doi.org/10.1179/msc.1974.8.1.298>
   — source of the first-order Wagner interaction parameters $e_i^j$ used for
   the steel bath.

## License

[MIT](LICENSE)
