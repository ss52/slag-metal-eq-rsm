# Slag Activity Model

Thermodynamic model for **electric-arc-furnace (EAF) slag–metal equilibrium**.

The package computes activity coefficients and activities of all slag oxide
components at a given temperature using **Ban-ya's quadratic formalism**
(a cation regular-solution model, RSM), coupled to the steel bath through a
**single shared oxygen potential**. Iron and chromium are entered analytically
as one total each and are split automatically into their valence states
(Fe²⁺/Fe³⁺ and Cr²⁺/Cr³⁺) by the redox equilibria.

[![Tests](https://img.shields.io/badge/tests-86%20passed-brightgreen)](#tests)
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

   The model reports both regular-solution cation activities and conventional
   oxide activities when a documented conversion is available. Ban-ya Eq. 24
   uses the two Fe regular-solution coefficients; the Fe–O oxygen-potential
   closure uses conventional $a_{FeO}$.

2. **Oxygen potential.** $P_{O_2}$ is *not* an input — the FeO-rich slag
   buffers it through the Fe–O equilibrium:

$$P_{O_2} = \left( \frac{a_{FeO}}{K_{FeO}\, a_{Fe}} \right)^2$$

3. **Iron redox.** The Fe³⁺/Fe²⁺ ratio is slaved to $P_{O_2}$
   (Ban-ya Eq. 24).

4. **Chromium redox.** The Cr²⁺/Cr³⁺ ratio is pinned by the chromium activity
   in the steel via $2\,\underline{CrO_{1.5}} + [Cr]_{1\,wt\%} = 3\,\underline{CrO}$
   (oxide reaction from Xiao & Holappa 1995; dissolution coefficient from Xiao,
   Kou & Fang 2018) — oxygen cancels, so no $P_{O_2}$ is needed. The model solves
   the **exact cubic** relation rather than the approximate square-root fixed
   point found in spreadsheet implementations.

5. **Metal side (WIPF).** Solute activity coefficients come from **Wagner's
   first-order interaction-parameter formalism** (coefficients from
   Sigworth & Elliott). The carbon content that is in equilibrium with the slag
   follows from the C–O reaction at the same $P_{O_2}$ — slag and steel are
   in balance *by construction*.

Because the oxygen potential is read off the slag and the balanced carbon is then computed from it, the coupled system converges in typically **5–30 damped fixed-point iterations**.

### Corrected equilibrium equations and numerical validity

Xiao and Holappa (1995, p. 322, Eqs. 9–10) report the oxide reaction with pure
solid chromium and its free energy in cal/mol. Xiao, Kou, and Fang (2018, p. 419,
Table 1) give the dissolution coefficient for chromium on the Henrian 1 wt%
standard state. The solver uses that dissolved-chromium standard state, so the
reaction energy is converted as

$$\Delta G^\circ_{1\,wt\%Cr}
=4.184(25690-13.36T)-(19246-46.86T)\ \mathrm{J\ mol^{-1}},\qquad
K_{Cr}=\exp\!\left(-\frac{\Delta G^\circ_{1\,wt\%Cr}}{RT}\right).$$

Ban-ya Eq. 24 uses both Fe coefficients on the regular-solution basis:

$$\log_{10}\!\left(\frac{n_{Fe^{3+}}}{n_{Fe^{2+}}}\right)
=\frac{6625}{T}-2.77+0.25\log_{10}P_{O_2}
+\log_{10}\gamma_{FeO}^{RS}-\log_{10}\gamma_{FeO_{1.5}}^{RS}.$$

The oxygen-potential closure still uses conventional $a_{FeO}$. For the
formula-unit conversion $P_2O_5(l)=2PO_{2.5}(RS)$, the reported activity is
$a_{P_2O_5}= (a_{PO_{2.5}}^{RS})^2\exp[(52720-230.706T)/(RT)]$;
conventional $a_{Al_2O_3}$ is unavailable unless a labeled custom conversion is
provided. The corresponding regular-solution cation activities remain
available.

Numeric inputs and solver states must be finite. Non-finite or undefined
iteration values raise an error, and JSON reports are serialized with strict
JSON number handling before an output file is opened, so a failed report cannot
leave a partial file.

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
| `fe2o3_conversion` | `none` / `{ "A": number, "B": number, "standard_state": string }` | `none` | Optional custom ΔG conversion for FeO₁.₅ (`ΔG = A + B·T`, J/mol of FeO₁.₅); the caller must identify its reference state. |
| `al2o3_conversion` | `none` / `{ "A": number, "B": number, "standard_state": string }` | `none` | Optional custom conversion for `Al2O3(reference) = 2 AlO1.5(RS)`; A and B are J/mol of Al₂O₃ formula unit and require a reference-state label. |
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
* Custom conversion objects must include numeric `A`, numeric `B`, and a
  non-empty `standard_state` label. Legacy `A`/`B`-only objects require this
  label to be added before they can be read.

## Output JSON reference

The report written with `-o` has `schema_version: 2`, followed by `input`,
`solution`, `metal`, `equilibrium`, and `warnings`.

### `input`

An exact echo of the resolved inputs: `temperature_K`, `P_CO_atm`,
`slag_wtpc`, `metal_wtpc`, and `options` (with all defaults filled in and the
custom conversions shown as `none` when unset or as labeled `{ "A", "B",
"standard_state" }` objects). The options object can be passed back to the
input parser without changing either custom conversion.

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
| `components`   | Array of per-cation rows distinguishing regular-solution and conventional activities (see below). |

Each modeled cation entry in `components`:

| Field | Description |
|-------|-------------|
| `cation` | Cation identity, such as `Fe2+` or `P5+`. |
| `rs_species` | Regular-solution species per cation; P⁵⁺ uses `PO2.5` and Al³⁺ uses `AlO1.5`. |
| `X_cation` | Cation mole fraction on the RSM basis. |
| `RTln_gamma_RS_J_per_mol_cation` | $RT\ln\gamma_i^{RS}$, in J/mol of cation. |
| `gamma_RS` | Regular-solution coefficient for that cation. |
| `a_RS` | Regular-solution activity, $X_{cation}\gamma_{RS}$. |
| `conventional` | `null` when no conversion is documented; otherwise the conventional species, its `standard_state`, number of RSM units per formula unit, conversion energy, and activity. |
| `conventional_unavailable_reason` | Explanation when `conventional` is `null`. |

The nested `conventional` object includes `gamma` and
`gamma_fraction_basis: "X_cation"` only for one-cation species. `P2O5` and
custom `Al2O3` conversions each use two RSM units and do not report a
conventional gamma on the cation-fraction basis. By default, conventional
activities for `FeO1.5` and `Al2O3` are unavailable. The optional legacy Ti
row appears only with `ti_handling: "as_excel"`; it has
`model_status: "denominator_only_legacy"` and carries no activity.

The `standard_state` field keeps the source wording. For example, the workbook
silica option is reported as `SiO2(s), polymorph unspecified; Xiao, Holappa &
Reuter (2002) Table IV`; the model does not guess an unstated crystal phase.

### Migrating from report v1

Read `schema_version` before consuming a report. In v2, component fields
`oxide`, `X`, `RTln_gamma_RS`, `DeltaG_conv`, `gamma`, and `a` are removed.
Use `cation`, `rs_species`, `X_cation`, `RTln_gamma_RS_J_per_mol_cation`,
`gamma_RS`, and `a_RS` for the RSM result. Read the nested `conventional`
object for a conventional oxide activity; it can be `null`, so consumers must
handle unavailability instead of treating it as zero. Custom conversion
inputs and echoed options now use labeled objects containing `A`, `B`, and
`standard_state`.

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

The suite covers source-equation checks for Cr equilibrium and Fe redox,
formula-unit activity conversions, numeric input/output safety, legacy workbook
regressions of the interaction matrix, and the reference coupled solve. The
paper-derived cases are an independent experimental envelope, not an oracle for
the parameter matrix.

## Code quality

Linted and formatted with [Ruff](https://github.com/astral-sh/ruff):

```bash
uv run ruff check .
uv run ruff format --check .
```

## Known data gaps

These are **explicit** — the model does not invent values:

1. **Al₂O₃ conventional activity** has no default conversion and is reported as
   unavailable. The regular-solution `AlO1.5` activity remains available; a
   conventional value requires a labeled custom conversion.
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
   — the CrO / CrO₁.₅ conversion factors and the Cr-cation interaction
   parameters.

4. **Y. Xiao, L. Holappa**, *Thermodynamics of Slags Containing Chromium Oxides*,
   INFACON 7, Trondheim (1995), 319–328.
   <https://www.pyrometallurgy.co.za/InfaconVII/319-Xiao.pdf>
   — the oxide reaction coefficients in Eqs. 9–10 (p. 322).

5. **S. Xiao, Q. Kou, X. Fang**, *Mutual Calculation between Standard Dissolved
   Gibbs Free Energy and Differential Dissolution Enthalpy Based on a Model of
   Dilute Solution*, Chinese Journal of Process Engineering **18** (2018), No. 2,
   417–421, doi:10.12034/j.issn.1009-606X.217281.
   <https://jproeng.ipe.ac.cn/CN/abstract/abstract3039.shtml>
   — the chromium dissolution coefficient in Table 1 (p. 419).

6. **Y. Xiao, L. Holappa, M.A. Reuter**, *Oxidation State and Activities of
   Chromium Oxides in CaO-SiO₂-CrOₓ Slag System*, Metallurgical and Materials
   Transactions B **33** (2002), 595–603.
   <https://link.springer.com/article/10.1007/s11663-002-0039-9>
   — experimental validation of the chromium oxidation state and activities.

7. **G.K. Sigworth, J.F. Elliott**, *The Thermodynamics of Liquid Dilute Iron
   Alloys*, Metal Science **8** (1974), No. 1, 298–310.
   <https://doi.org/10.1179/msc.1974.8.1.298>
   — source of the first-order Wagner interaction parameters $e_i^j$ used for
   the steel bath.

## License

[MIT](LICENSE)
