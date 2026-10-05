# Slag Activity Model

Thermodynamic calculations for **prescribed EAF-type slag and metal compositions**.

The package computes activity coefficients and activities of modeled slag oxide
components at a given temperature using **Ban-ya's quadratic formalism**
(a cation regular-solution model, RSM). It uses a **Fe-buffer oxygen potential**
for Fe redox and carbon, with a separate chromium comproportionation constraint.
Iron and chromium are entered analytically as one total each and are split
automatically into their valence states (Fe²⁺/Fe³⁺ and Cr²⁺/Cr³⁺) by these
selected relations.

[Tests](#tests)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Python](https://img.shields.io/badge/python-%3E%3D3.13-blue)](https://www.python.org)
[![License: MIT](https://img.shields.io/badge/license-MIT-yellow)](LICENSE)

---

## Physical basis

In coupled mode, the model links Fe–O, Fe redox and carbon through the Fe-buffer
oxygen potential $P_{O_2}$, and separately solves chromium comproportionation
while retaining the prescribed compositions:

1. **Slag side (RSM).** The slag is treated as a regular solution of cations
   (Ban-ya 1993). Oxide-component activity coefficients follow from cation
   fractions $X_i$ and pair interaction energies $\alpha_{ij}$:

$$RT\ln\gamma_i^{RS} = \sum_j \alpha_{ij} X_j^2 + \sum_{\substack{j<k \\ j,k\ne i}} (\alpha_{ij} + \alpha_{ik} - \alpha_{jk}) X_j X_k$$

   The model reports both regular-solution cation activities and conventional
   oxide activities when a documented conversion is available. Ban-ya Eq. 24
   uses the two Fe regular-solution coefficients; the Fe–O oxygen-potential
   closure uses conventional $a_{FeO}$.

2. **Oxygen potential.** In coupled mode, $P_{O_2}$ is calculated from the
   slag's FeO activity through the Fe–O equilibrium:

$$P_{O_2} = \left( \frac{a_{FeO}}{K_{FeO}\, a_{Fe}} \right)^2$$

3. **Iron redox.** The Fe³⁺/Fe²⁺ ratio is slaved to $P_{O_2}$
   (Ban-ya Eq. 24).

4. **Chromium comproportionation.** The Cr²⁺/Cr³⁺ ratio is pinned by the
   chromium activity in the steel via $2\,\underline{CrO_{1.5}} + [Cr]_{1\,wt\%} = 3\,\underline{CrO}$
   (oxide reaction from Xiao & Holappa 1995; dissolution coefficient from Xiao,
   Kou & Fang 2018) — oxygen cancels, so no $P_{O_2}$ is needed. The model solves
   the **exact cubic** relation rather than the approximate square-root fixed
   point found in spreadsheet implementations.

   The independent equilibrium $CrO_{1.5}(s)=CrO(l)+\tfrac14 O_2(g)$
   ([Xiao and Holappa 1995, p. 324, Eqs. 14–15](https://www.pyrometallurgy.co.za/InfaconVII/319-Xiao.pdf))
   is not enforced. The chromium valence constraint therefore need not be
   consistent with the Fe-buffer oxygen potential; the calculation is not a
   complete slag–metal redox equilibrium.

   For the default reference EAF case, the Fe-buffer $P_{O_2}$ is about
   $3.72\times10^{-10}$ atm, whereas Eq. 15 with the same oxide reference
   states and calculated activities implies $9.06\times10^{-10}$ atm, about
   2.44 times higher. The Cr–O $Q/K$ at the reported pressure is about 0.80.
   This illustrates an unenforced source-equation discrepancy, not solver
   roundoff or measured-activity fit scatter; these values are neither an
   empirical measurement nor an accuracy tolerance.

5. **Metal side (WIPF).** Solute activity coefficients come from **Wagner's
   first-order interaction-parameter formalism**, using auxiliary metal-side
   data discussed below. Carbon is calculated from the C–O reaction at the same
   $P_{O_2}$. The selected FeO, CO, chromium comproportionation and Fe-redox
   relations are therefore satisfied at the prescribed compositions, within
   solver tolerance.

Slag elemental Fe/Cr totals, metal Cr/Mn/P contents, temperature and $P_{CO}$
remain prescribed. The model updates Fe/Cr valence ratios and carbon; it has
no phase masses, slag/metal mass ratio, gas inventory or coupled element
transfer. It does not predict closed-system equilibration or Mn/P partition
coefficients. The coupled calculation uses damped fixed-point iteration.

The separate `solve_redox_fixed_po2` Python function imposes $P_{O_2}$ and carbon.
It enforces Fe redox and chromium comproportionation; FeO and CO equilibrium
are not enforced. The output identifies this mode and evaluates every reported
$Q/K$ at the imposed pressure.

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

## Scientific applicability and tested evidence

The checks below establish specific computational and source comparisons.
They do not establish a validated temperature/composition box or experimental
validation of the complete coupled calculation.

| Tested calculation | Temperature and composition | What the evidence establishes |
|---|---|---|
| Reference EAF solve, [`eaf_slag_01.json`](examples/eaf_slag_01.json), default `banya93_xiao95_v1` profile | 1823.15 K, $P_{CO}=1$ atm; prescribed metal Cr 0.038, Mn 0.007, P 0.069 wt%; solved C about 0.0308 wt% | Computational regression and selected equilibrium identities at this input; no experimental full-solver validation. |
| Chromium-free equimolar CaO–SiO₂ illustration, `banya93_casi_xiao95cr_hybrid_v1` with Ban-ya silica conversion | 1823 K, $X_{Ca}=X_{Si}=0.5$, $\alpha_{Ca,Si}=-133890$ J/mol; β-cristobalite reference | Equation-derived $a_{SiO_2}=0.2575019792$, $a_{CaO}=0.01103118178$; no measured activity row. |
| All 11 Xiao and Holappa (1995) Table 1 SiO₂–CrO–CrO₁.₅ rows, default profile | 1873 K, source-reported valence splits, equilibrium with solid Cr | Isolated slag-activity comparisons using parameter-assessment data; in-sample, with measured Cr oxide and Gibbs–Duhem-derived silica activities distinguished below. |

The reference slag assay is SiO₂ 27.495, CaO 35.597, MgO 1.999, Al₂O₃ 7.973,
MnO 4.762, TiO₂ 0.725, total Fe as FeO 18.851, total Cr as Cr₂O₃ 2.580, and
P₂O₅ 0 wt%. This single assay is a tested point, not a general EAF applicability
envelope. Exact paper compositions and the analytical illustration are in
[`paper_cases.json`](examples/thermodynamic_review/paper_cases.json).

Ban-ya's regular-solution treatment is an approximate liquid-slag model whose
applicability depends on composition. This implementation does not predict
solid phases, spinel formation or phase separation. Xiao's chromium parameters
cover the studied Ca–Si–Mg–Al–Cr systems. Fe–Cr, Mn–Cr and P–Cr slag pair
interactions are assumed zero because parameters are unavailable; Fe-rich
chromium slags therefore involve extrapolation, not measured zero interactions.

Ti is excluded and the active cations renormalized by default, an implementation
approximation. Legacy `as_excel` is denominator-only dilution without a physical
Ti activity; Ban-ya lists nonzero Ti interactions that this implementation does
not use. No universal safe Ti wt% limit has been established.

The metal model uses first-order coefficients tabulated at 1873 K approximately
at the 1823.15 K reference temperature. Second-order terms and $f_O$ corrections
are absent. Carbon-rich hot metal (for example, 3–5 wt% C) and the complete
coupled model have not been empirically validated here. The FeO/CO equilibrium,
carbon-dissolution and oxygen-saturation fits are inherited auxiliary data with
incomplete directly verified provenance; they are not derived from the two
oxide papers.

Extreme chromium roots outside the current $10^{-6}$–$10^6$ bracket and trace
species near floating-point roundoff remain deferred numerical findings. Results
in those regions require further tests; passing the tested cases does not
establish their reliability.

## Comparisons with published slag activities

These comparisons exercise the isolated slag-activity calculation. They do not
test the complete coupled EAF slag–metal solve. Relative differences are
calculated as `calculated / reported - 1`.

Ban-ya (1993) provides equations, interaction energies and conversion factors,
but no source-exact measured numerical row for this implementation. Its checks
therefore cover Eq. 24, selected Table 1 interaction energies, and Table 3
conversion factors; they do not claim empirical activity validation.

### Xiao and Holappa 1995, Table 1

The 11 SiO₂–CrO–CrO₁.₅ compositions and chromium oxide activities are from
[Table 1](https://www.pyrometallurgy.co.za/InfaconVII/319-Xiao.pdf), at 1873 K
in equilibrium with solid chromium. The paper reports the chromium oxide
activities from EMF/oxygen-potential measurements. Its listed silica activities
were calculated by Gibbs–Duhem, so they are shown separately and are not counted
as direct measurements. The reported Cr²⁺ fraction and log $P_{O_2}$ are also
kept in the source fixture; they are not predictions from this isolated activity
comparison.

| Case | Reported → calculated $a_{CrO}$ (relative difference) | Reported → calculated $a_{CrO_{1.5}}$ (relative difference) | Gibbs–Duhem $a_{SiO_2}$ → calculated |
|---|---:|---:|---:|
| T1-01 | 0.91 → 0.966670 (+6.23%) | 0.95 → 0.924719 (−2.66%) | 0.12 → 0.085732 (−28.56%) |
| T1-02 | 0.94 → 0.773519 (−17.71%) | 1.00 → 1.005703 (+0.57%) | 0.06 → 0.069266 (+15.44%) |
| T1-03 | 0.80 → 1.024726 (+28.09%) | 0.79 → 0.849820 (+7.57%) | 0.19 → 0.107210 (−43.57%) |
| T1-04 | 0.72 → 0.630419 (−12.44%) | 0.67 → 0.749387 (+11.85%) | 0.37 → 0.201461 (−45.55%) |
| T1-05 | 0.66 → 0.824443 (+24.92%) | 0.59 → 0.736625 (+24.85%) | 0.36 → 0.179744 (−50.07%) |
| T1-06 | 0.54 → 0.571285 (+5.79%) | 0.43 → 0.540390 (+25.67%) | 0.66 → 0.371085 (−43.77%) |
| T1-07 | 0.55 → 0.331747 (−39.68%) | 0.44 → 0.493709 (+12.21%) | 0.75 → 0.515357 (−31.29%) |
| T1-08 | 0.77 → 0.762132 (−1.02%) | 0.74 → 0.817495 (+10.47%) | 0.26 → 0.147549 (−43.25%) |
| T1-09 | 0.62 → 0.607896 (−1.95%) | 0.53 → 0.691102 (+30.40%) | 0.46 → 0.243303 (−47.11%) |
| T1-10 | 0.60 → 0.271959 (−54.67%) | 0.51 → 0.535309 (+4.96%) | 0.70 → 0.481352 (−31.24%) |
| T1-11 | 0.55 → 0.372342 (−32.30%) | 0.44 → 0.472193 (+7.32%) | 0.77 → 0.526504 (−31.62%) |

The Table 2 parameters were assessed using these Table 1 measurements together
with earlier chromium-slag data (pp. 325–326), so these are in-sample
comparisons, not independent validation. The model's mean absolute relative
differences are 20.4% for CrO, 12.6% for CrO₁.₅ and 37.4% against the
Gibbs–Duhem silica values. The discrepancies are not uniform: for example,
calculated CrO activity in T1-10 is 54.67% below the reported value. The paper's
Fig. 4 also shows scatter between calculated and measured chromium activities
across systems; these numerical differences are reported as-is, without an
activity-error acceptance threshold.

At 1873 K, the printed CrO and CrO₁.₅ values give $Q/K$ values from 0.981 to
1.028 when $Q=a_{CrO}^3/a_{CrO_{1.5}}^2$ and $K$ uses the pure-solid-Cr reaction
in Xiao and Holappa Eqs. 9–10. All eleven rows' intervals obtained from their
printed two-decimal activity precision include the calculated $K$; this is a
rounding-interval check on the source data, not a model-to-measurement pass
criterion. The formation relation for liquid CrO cited on p. 321 was available
only from 1665–1750 °C and was extrapolated by the authors to their measurements
at 1600 °C, one limitation on treating the tabulated values as precision ground
truth. Exact compositions and source-reported redox values are in
[`paper_cases.json`](examples/thermodynamic_review/paper_cases.json).

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
$\gamma$ and $a$, prescribed metal and solved carbon, equilibrium diagnostics
$Q/K$) to JSON.

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
| `parameter_profile` | `banya93_xiao95_v1` / `banya93_casi_xiao95cr_hybrid_v1` / `legacy_workbook_typo_hybrid_v1` | `banya93_xiao95_v1` | Named, versioned immutable cation-interaction matrix; echoed in each v2 report. |
| `ti_handling`      | `exclude_renormalize` / `as_excel`       | `exclude_renormalize`  | Default excludes Ti and renormalizes active cations; `as_excel` emulates denominator-only workbook dilution without Ti interaction terms or physical Ti activity. |
| `sio2_conversion`  | `workbook` / `banya`                     | `workbook`             | SiO₂ standard-state conversion. The 2.6 kJ difference moves γ(SiO₂) by ~1.19×. |
| `fe2o3_conversion` | `none` / `{ "A": number, "B": number, "standard_state": string }` | `none` | Optional custom ΔG conversion for FeO₁.₅ (`ΔG = A + B·T`, J/mol of FeO₁.₅); the caller must identify its reference state. |
| `al2o3_conversion` | `none` / `{ "A": number, "B": number, "standard_state": string }` | `none` | Optional custom conversion for `Al2O3(reference) = 2 AlO1.5(RS)`; A and B are J/mol of Al₂O₃ formula unit and require a reference-state label. |
| `a_Fe`             | `unity` / `xfe`                          | `unity`                | Iron activity: 1.0 (pure liquid Fe) or the metal mole fraction. |
| `iron_input`       | `FeO_total` / `Fe_total`                 | `FeO_total`            | Iron input form (must match the key used in `slag_wtpc`). |
| `chromium_input`   | `Cr2O3_total` / `Cr_total`               | `Cr2O3_total`          | Chromium input form (must match the key used in `slag_wtpc`). |

The default profile uses Ban-ya (1993) Table 1 base interactions, including
the corrected Fe³⁺–Ca²⁺ value −95,810 J/mol, and Xiao and Holappa (1995) Table 2
values for Ca²⁺–Si⁴⁺ (−139,100 J/mol) and the chromium-bearing terms. The
`banya93_casi_xiao95cr_hybrid_v1` comparison changes only Ca²⁺–Si⁴⁺ to Ban-ya's
−133,890 J/mol; its chromium terms remain from Xiao and Holappa, so it is a
hybrid rather than a pure Ban-ya set. The `legacy_workbook_typo_hybrid_v1`
profile preserves the historical workbook snapshot, including its incorrect
Fe³⁺–Ca²⁺ value −96,810 J/mol, and is not a source parameter set. Reports record
the profile ID and version under both `input.options.parameter_profile` and
`parameter_profile`.

### Validation rules

* Exactly one of `FeO_total` / `Fe_total` and one of `Cr2O3_total` / `Cr_total`
  must be present; unknown keys, negative values, or an inconsistent
  `iron_input`/`chromium_input` flag raise an `InputError`.
* The current 100 ± 1 wt% warning sums raw input values. It can give a false
  warning for an equivalent assay expressed with elemental Fe/Cr totals; that
  input-basis defect remains open. After redox splitting, oxide mass can change
  because the valence states carry different oxygen amounts; the calculation
  does not close an oxygen inventory.
* Missing `metal_wtpc` components default to 0.0.
* Custom conversion objects must include numeric `A`, numeric `B`, and a
  non-empty `standard_state` label. Legacy `A`/`B`-only objects require this
  label to be added before they can be read.

## Output JSON reference

The report written with `-o` has `schema_version: 2`, a root-level
`parameter_profile` identity, then `input`, `solution`, `metal`, `equilibrium`,
and `warnings`. The selected profile ID is also preserved in `input.options` so
the report input can be parsed again without losing that choice.

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
| `calculation_mode` | `coupled` or `fixed_P_O2`; identifies which relations are enforced. |
| `P_O2_atm`     | Oxygen potential calculated from FeO in coupled mode, or imposed in fixed-$P_{O_2}$ mode (atm). |
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
silica option is reported as `SiO2(solid), polymorph unspecified; Xiao & Holappa
(1995), silica-in-lime-silica section`; the model does not guess an unstated
crystal phase.

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

Prescribed metal solutes, calculated activities and mode-dependent carbon and
oxygen diagnostics.

| Field    | Description |
|----------|-------------|
| `C_wtpc` | Carbon from C–O equilibrium in coupled mode, or prescribed carbon in fixed-$P_{O_2}$ mode (wt%). |
| `O_wtpc` | Approximate FeO-equilibrium oxygen estimate in coupled mode; `null` in fixed-$P_{O_2}$ mode. |
| `O_FeO_equivalent_wtpc` | Separately named approximate FeO-equivalent estimate in either mode; it need not describe oxygen at the imposed pressure. |
| `oxygen_diagnostic` | Status and basis: `approximate_FeO_equilibrium` for coupled mode or `FeO_equivalent_only` for fixed mode. |
| `f_C`    | Henrian activity coefficient of carbon (WIPF). |
| `f_Cr`   | Henrian activity coefficient of chromium. |
| `f_Mn`   | Henrian activity coefficient of manganese. |
| `f_P`    | Henrian activity coefficient of phosphorus. |
| `a_C`    | Carbon activity $a_C = f_C[\%C]$ (Henrian 1 wt% scale). |
| `a_Cr`   | Chromium activity $a_{Cr} = f_{Cr}[\%Cr]$. |

The oxygen estimate is explicitly approximate:

$$[\%O]_{FeO\ equivalent}
=\frac{a_{FeO}}{a_{Fe}}\,10^{-6320/T+2.734}.$$

It uses the selected iron activity, assumes $f_O=1$, and retains the inherited
empirical saturation fit. Coupled mode reports the same estimate under both
oxygen fields. Fixed mode cannot infer dissolved oxygen from this FeO relation
at its imposed pressure because FeO equilibrium is not enforced.

### `equilibrium`

| Field         | Description |
|---------------|-------------|
| `K_FeO`       | Equilibrium constant of $Fe(l) + \tfrac12 O_2 = FeO(l)$ at $T$. |
| `K_CO`        | Equilibrium constant of $[C] + \tfrac12 O_2 = CO(g)$ (Henrian C). |
| `K_Cr`        | Equilibrium constant of chromium comproportionation $2\,CrO_{1.5} + [Cr] = 3\,CrO$ at $T$. |
| `enforced_relations` | `['FeO', 'CO', 'Cr', 'Fe_redox']` in coupled mode; `['Cr', 'Fe_redox']` in fixed-$P_{O_2}$ mode. |
| `Q_over_K_FeO`| Fe–O quotient at the reported pressure; approximately 1 within tolerance in coupled mode, diagnostic only in fixed mode. |
| `Q_over_K_CO` | C–O quotient at the reported pressure; approximately 1 within tolerance in coupled mode, diagnostic only in fixed mode. |
| `Q_over_K_Cr` | Reaction-quotient / $K$ check for the enforced Cr comproportionation relation only; does not certify Cr–O equilibrium (`null` if no Cr present). |

An enforced-relation identity near unity checks the calculation's consistency;
it does not establish experimental accuracy. A non-unity FeO or CO quotient in
fixed mode describes the imposed state and is not a failure to enforce a
relation that this mode does not solve.

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
tests/             # source and computational regression checks
PLAN.md            # original implementation specification
```

[PLAN.md](PLAN.md) records the original implementation specification. Use this
README, the named parameter-profile metadata and the current source fixtures
for the corrected equations, active datasets and evidence limitations.

## Tests

```bash
uv run pytest -v
```

The suite covers source-equation checks for Cr equilibrium and Fe redox,
formula-unit activity conversions, numeric input/output safety, legacy workbook
regressions of the interaction matrix, and the reference coupled solve. The
paper-review pack checks Ban-ya 1993 equations and constants, preserves all
eleven Xiao and Holappa 1995 Table 1 rows, and reports their in-sample model
differences without an activity-error acceptance threshold.

The separate paper-review pack can be run with:

```bash
uv run pytest examples/thermodynamic_review -p no:cacheprovider -ra
```

Ban-ya 1993 has no source-exact measured numerical row in this pack. The Xiao
comparison distinguishes measured Cr oxide activities from the paper's
Gibbs–Duhem-derived SiO₂ values and flags that Table 1 was used in parameter
assessment. Neither paper pack result is presented as independent or full-solver
validation. A separate chromium-free CaO–SiO₂ illustration at 1823 K uses
Ban-ya's −133,890 J/mol Ca–Si interaction and β-cristobalite conversion to
reproduce hand-calculated activities; it is equation-derived, not experimental
evidence. The reference solve is retained as a regression sentinel for the
selected parameter profile, not as empirical validation.

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
2. **FeO₁.₅ conventional activity** has no default documented conversion and is
   unavailable unless a labeled custom conversion is provided.
3. **Fe–Cr, Mn–Cr, P–Cr pair interactions are assumed zero** because parameters
   are unavailable. The Xiao and Holappa (1995) chromium set covers the studied
   Ca/Mg/Si/Al/Cr partners; the zero assumptions are not measurements.
4. **WIPF $e_i^j$** are 1873 K values approximately applied at 1823.15 K.
   Neglecting higher-order terms has not been validated for carbon-rich hot
   metal. The corrected $e_{Mn}^{C}=-0.07$ is from Lei, Fu and Xiong (2024),
   Table 1, p. 112, which cites Chen's 2010 steelmaking databook; the original
   Sigworth and Elliott table was not directly verified for this coefficient.
   The distinct $e_C^{Mn}=-0.012$ is retained; Wagner coefficients are asymmetric.

## Related literature (context only)

The oxide parameters and selected paper-activity rows above come from the two
primary sources, Ban-ya (1993) and Xiao and Holappa (1995). The studies below
help interpret chromium activity and redox behavior; none supplies oxide-model
coefficients or a selected benchmark row in this project. Comparing numerical
activities across them requires matching temperature, slag and metal composition,
oxygen potential, oxide phase, and activity standard state.

- **Xiao and Holappa (1993)** measured CrO and CrO₁.₅ activities and chromium
  valence in CaO–SiO₂-based slags, including MgO and Al₂O₃ additions, by an EMF
  method. This is a companion study by the authors of the 1995 source, not an
  independent model validation. [ISIJ International **33**, 66–74](https://doi.org/10.2355/isijinternational.33.66).
- **Park, Min, and Rhee (1998)** measured chromium activity in liquid Fe–Cr–C
  at 1823 K, finding Henrian behavior up to 3 mass% Cr. It informs the metal
  side only within that dilute-Cr range.
  [ISIJ International **38**, 1287–1291](https://doi.org/10.2355/isijinternational.38.1287).
- **Itoh, Nagasaka, and Hino (2000)** determined the dissolved Cr–O equilibrium
  and a Wagner Cr–O interaction parameter in Fe–Cr melts saturated with pure
  solid Cr₂O₃ at 1823–1923 K.
  [ISIJ International **40**, 1051–1058](https://doi.org/10.2355/isijinternational.40.1051).
- **Kimoto, Itoh, Nagasaka, and Hino (2002)** studied Fe–Cr–O equilibrium with
  FeO·Cr₂O₃ solid solution. They report that below approximately 7 mass% Cr in
  the metal, this oxide phase replaces pure Cr₂O₃ at equilibrium.
  [ISIJ International **42**, 23–32](https://doi.org/10.2355/isijinternational.42.23).
- **Dong, Wang, and Seetharaman (2009)** measured CrO activity in
  CaO–SiO₂–MgO–Al₂O₃ slags at 1803–1923 K by gas–slag equilibrium. Their oxygen
  partial pressures of 10⁻³–10⁻⁵ Pa differ from strongly reducing cases.
  [steel research international **80**, 202–208](https://doi.org/10.2374/SRI08SP124).
- **Shi, Li, Liu, and Kobayashi (2024)** measured chromium-oxide activities in
  CaO–SiO₂–MgO–Al₂O₃–MnO–CaF₂ slags at 1823 K and low initial CrO₁.₅ content
  (at most 0.5 mass%) to avoid spinel formation. The extra MnO/CaF₂ and dilute
  chromium range limit direct comparison with the selected 1995 rows.
  [ISIJ International **64**, 893–900](https://doi.org/10.2355/isijinternational.ISIJINT-2023-214).

## References

1. **S. Ban-ya**, *Mathematical Expression of Slag-Metal Reactions in
   Steelmaking Process by Quadratic Formalism Based on the Regular Solution
   Model*, ISIJ International **33** (1993), No. 1, 2–11.
   <https://www.jstage.jst.go.jp/article/isijinternational1989/33/1/33_1_2/_article>
   — the quadratic (regular-solution) formalism, the α-matrix, the Fe-redox
   Eq. 24, and the oxide conversion factors.

2. **Y. Xiao, L. Holappa**, *Thermodynamics of Slags Containing Chromium Oxides*,
   INFACON 7, Trondheim (1995), 319–328.
   <https://www.pyrometallurgy.co.za/InfaconVII/319-Xiao.pdf>
   — Table 1 chromium oxide activities and compositions, Table 2 assessed
   chromium interactions, Table 3 oxide conversions, and the oxide reaction
   coefficients in Eqs. 9–10 (p. 322).

3. **S. Xiao, Q. Kou, X. Fang**, *Mutual Calculation between Standard Dissolved
   Gibbs Free Energy and Differential Dissolution Enthalpy Based on a Model of
   Dilute Solution*, Chinese Journal of Process Engineering **18** (2018), No. 2,
   417–421, doi:10.12034/j.issn.1009-606X.217281.
   <https://jproeng.ipe.ac.cn/CN/abstract/abstract3039.shtml>
   — an external metal-side chromium dissolution coefficient in Table 1 (p. 419).

4. **G.K. Sigworth, J.F. Elliott**, *The Thermodynamics of Liquid Dilute Iron
   Alloys*, Metal Science **8** (1974), No. 1, 298–310.
   <https://doi.org/10.1179/msc.1974.8.1.298>
   — a historical reference for first-order Wagner metal-side interaction
   parameters; the complete original table has not been directly reverified
   in this audit.

5. **J. Lei, Y. Fu, L. Xiong**, *Study on Mn Volatilization Behavior During
   Vacuum Melting of High-manganese Steel*, Archives of Foundry Engineering
   **24** (2024), No. 2, 110–116, doi:10.24425/afe.2024.149277.
   [Table 1, p. 112](https://journals.pan.pl/Content/131630/AFE%202_2024_13-Final.pdf)
   — the checked compilation for $e_{Mn}^{C}=-0.07$ at 1873 K; its reference
   [22] is Chen (2010), *Common Charts and Databook for Steelmaking*.

## License

[MIT](LICENSE)
