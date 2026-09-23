# P0 Thermodynamic Corrections Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Correct the four critical thermodynamic and numeric-safety defects so the coupled result and v2 activity report have consistent, declared standard states.

**Architecture:** Keep the existing solver and parameter matrix, but make the reaction basis explicit in equilibrium data and separate RS from conventional slag activities in one value object. Each solver equation consumes only the activity basis stated by its source; v2 output reports the same distinctions, and finite-number guards prevent false success.

**Tech Stack:** Python 3.13, uv, NumPy, pytest, Ruff, local `.venv`.

**Spec:** `examples/thermodynamic_review/P0_DESIGN.md`

## Global Constraints

- Implement C1-C4 from `examples/thermodynamic_review/REVIEW.md` only; keep normal- and low-priority corrections separate.
- Use the cited primary papers and `examples/thermodynamic_review/paper_cases.json` for equations and literal coefficients, not the workbook as an oracle.
- Thermodynamic energies are J/mol of the written reaction at K; dissolved Cr uses the Henrian 1 mass-% standard state.
- The v2 report has `schema_version: 2`, no ambiguous legacy `gamma`/`a` activity fields, and no non-standard JSON numbers.
- Preserve unrelated worktree changes; do not edit the workbook or change the interaction/Wagner parameter matrix in this plan.
- Do not use a numerical match to the workbook to choose or fit any conversion coefficient.

## Review Focus

- Bare JSON `NaN`/`Infinity`, quoted non-finite strings, and huge integers: Task 5 tests require `InputError` with the field name.
- A trailing `NaN` target in solver convergence order: Task 5 test requires `NumericalStateError`, not a returned result.
- Zero P2O5 input and nonzero P2O5 input: Task 3 tests require zero and the two-power formula respectively.
- Default Al2O3 and legacy custom `A,B`-only Al2O3 input: Tasks 3-4 tests require unavailable conventional activity and an explicit migration error.
- Ti `as_excel` output and serializer failure: Tasks 4-5 tests require an explicit denominator-only status and no partial report file.

## File map and interfaces

- `slag_model/data.py`: separately named, unit-bearing Cr reaction/dissolution coefficients and mappings from cations to RS species and documented conventional references.
- `slag_model/equilibrium.py`: `k_Cr(T)` on the same 1 mass-% Cr basis as the solver's metal activity.
- `slag_model/redox.py`: Ban-ya Fe-ratio signature whose coefficient names say `RS`.
- `slag_model/rsm.py`: `ConversionSpec` and `SlagActivities`; `build_slag_activities(...)` performs all standard-state conversions and formula-unit stoichiometry.
- `slag_model/solver.py`: consumes explicit activity maps; owns finite iteration/state checks and `NumericalStateError`.
- `slag_model/io.py`: validates finite input and custom conversion schema; translates `SlagActivities` to v2 report without recomputing thermodynamics.
- `slag_model/cli.py`: displays v2 labels and pre-serializes strict JSON before opening the output file.
- `tests/` and `examples/thermodynamic_review/`: paper-equation, integration, input, and CLI tests. Existing workbook assertions affected by C1/C2 are replaced with source-equation assertions rather than widened.

---

### Task 1: Put chromium equilibrium on the dissolved-metal basis

**Files:**
- Modify: `slag_model/data.py` (Cr free-energy constants)
- Modify: `slag_model/equilibrium.py` (`k_Cr`)
- Modify: `tests/test_equilibrium.py` (remove incorrect 0.916 expectations)
- Modify: `examples/thermodynamic_review/test_paper_validation.py` (remove C1 `xfail`)

**Interfaces:**
- Consumes: `R`, temperature in K, the published solid-Cr oxide reaction and Cr dissolution reactions.
- Produces: `k_Cr(T: float) -> float`, dimensionless for `2 CrO1.5 + [Cr]1wt = 3 CrO`.

- [ ] **Step 1: Verify the source ledger.** Recheck the 1995 Xiao-Holappa reaction, its cal/mol unit and solid-Cr basis, plus the 1 mass-% Cr dissolution relation and its J/mol basis. Record page/equation locators beside the constants. If the exact `19246-46.86T` relation cannot be established for this reference state, stop this task and ask for a source decision; do not substitute a workbook-fitted value.

- [ ] **Step 2: Make the existing paper test fail normally.** Remove only its strict `xfail` marker and run:

```powershell
uv run pytest examples/thermodynamic_review/test_paper_validation.py::test_xiao_chromium_equilibrium_is_on_library_metal_standard_state -q -p no:cacheprovider
```

Expected: FAIL because current `k_Cr(1823.15)` is about 0.915831 instead of 0.00878747.

- [ ] **Step 3: Replace the ambiguous constant and calculate the compound reaction.** Use separate literals and the exact source-to-library conversion:

```python
# data.py; reaction 2 CrO1.5 + Cr(s) = 3 CrO, cal/mol reaction.
DELTA_G_CR_SOLID_CAL = (25690.0, -13.36)
# Cr(s) = [Cr] on the Henrian 1 mass-% scale, J/mol Cr.
DELTA_G_CR_DISSOLUTION_J = (19246.0, -46.86)

# equilibrium.py
def k_Cr(T: float) -> float:
    dg_oxide_j = 4.184 * _delta_g(DELTA_G_CR_SOLID_CAL, T)
    dg_dissolution_j = _delta_g(DELTA_G_CR_DISSOLUTION_J, T)
    return math.exp(-(dg_oxide_j - dg_dissolution_j) / (R * T))
```

Remove the old `DELTA_G_CR` symbol or make any retained alias impossible to mistake for a J/mol Henrian reaction.

- [ ] **Step 4: Replace workbook-derived `k_Cr` checks.** In `tests/test_equilibrium.py`, assert the independent compound-reaction expression at multiple temperatures, rather than `0.90 <= K <= 0.93` or `K≈0.916`:

```python
@pytest.mark.parametrize("temperature", [1773.15, 1823.15, 1873.15])
def test_k_Cr_source_reaction(temperature):
    dg = 4.184 * (25690.0 - 13.36 * temperature) - (19246.0 - 46.86 * temperature)
    assert k_Cr(temperature) == pytest.approx(math.exp(-dg / (R * temperature)), rel=1e-12)
```

- [ ] **Step 5: Run the focused tests.** Task 2 updates the workbook-coupled redox assertions before a combined C1+C2 commit; do not treat their interim failures as scientific evidence against the corrected `k_Cr`.

```powershell
uv run pytest tests/test_equilibrium.py examples/thermodynamic_review/test_paper_validation.py::test_xiao_chromium_equilibrium_is_on_library_metal_standard_state -q -p no:cacheprovider
```

Expected: focused tests PASS. Keep these changes together with Task 2 for the redox commit so no commit knowingly leaves the workbook-coupled test suite red.

### Task 2: Keep both Fe redox coefficients on the RS basis

**Files:**
- Modify: `slag_model/redox.py` (parameter names/docstring)
- Modify: `slag_model/solver.py` (both Ban-ya call sites)
- Modify: `tests/test_rsm.py` (replace affected workbook redox expectations)
- Modify: `examples/thermodynamic_review/test_paper_validation.py` (remove C2 `xfail`, add fixed-pressure test)

**Interfaces:**
- Consumes: `RTln_gamma_rs` from `rsm_gamma_rtln`, `R`, `T`, and `P_O2`.
- Produces: `fe_ratio_ban_ya(T, P_O2, gamma_FeO_RS, gamma_FeO1_5_RS)`; the FeO oxygen buffer still consumes converted `a_FeO`.

- [ ] **Step 1: Unmark and run the existing coupled Eq. 24 test.**

```powershell
uv run pytest examples/thermodynamic_review/test_paper_validation.py::test_banya_eq_24_is_coupled_with_regular_solution_gammas -q -p no:cacheprovider
```

Expected: FAIL at the FeO conversion factor.

- [ ] **Step 2: Add a failing fixed-pressure Eq. 24 test.** Place it beside the coupled paper test:

```python
def test_fixed_po2_fe_redox_uses_two_rs_gammas():
    cfg = load_input(REFERENCE_INPUT)
    state = solve_redox_fixed_po2(cfg, P_O2=1e-9, C=0.018)
    rt = R * cfg.temperature_K
    expected = fe_ratio_ban_ya(
        cfg.temperature_K,
        1e-9,
        math.exp(state.rtln_gamma_rs["Fe2+"] / rt),
        math.exp(state.rtln_gamma_rs["Fe3+"] / rt),
    )
    assert state.r_Fe == pytest.approx(expected, rel=1e-8)
```

Import `solve_redox_fixed_po2` in that file. Run this test and confirm FAIL before changing the solver.

- [ ] **Step 3: Fix both callers and clarify the signature.** Until Task 4 switches the solver to `SlagActivities`, calculate both coefficients directly from the same RS map:

```python
rt = R * T
gamma_feo_rs = math.exp(rtln_rs["Fe2+"] / rt)
gamma_feo1_5_rs = math.exp(rtln_rs["Fe3+"] / rt)
r_Fe_n = fe_ratio_ban_ya(T, P_O2, gamma_feo_rs, gamma_feo1_5_rs)
```

In `solve_redox_fixed_po2`, use `current_state.rtln_gamma_rs` and `cfg.temperature_K` identically. Rename `fe_ratio_ban_ya` arguments to `gamma_FeO_RS` and `gamma_FeO1_5_RS`; retain converted FeO only for `p_o2_from_slag`.

- [ ] **Step 4: Decouple the workbook coefficient regression from the wrong redox solution.** The current workbook fixture's historical fixed-pressure result is `r_Fe=0.13073478698367935`, `r_Cr=1.2754258847919098`. Preserve its RSM coefficient checks by evaluating this *fixed composition*, not resolving redox after the equations change:

```python
@pytest.fixture(scope="module")
def workbook_state():
    cfg = load_input(SLAG_PATH)
    cfg.options.cross_terms = "major5"
    cfg.options.ti_handling = "as_excel"
    cfg.options.sio2_conversion = "workbook"
    state, _ = compute_state(
        cfg, r_Fe=0.13073478698367935,
        r_Cr=1.2754258847919098, C_wtpc=0.018,
    )
    return state
```

Import `compute_state` for this fixture. Remove `test_r_Fe` and `test_r_Cr` as workbook-oracle assertions. Add a separate fixed-pressure chromium equilibrium identity test on an actual `solve_redox_fixed_po2` result:

```python
def test_fixed_po2_chromium_reaction():
    cfg = load_input(SLAG_PATH)
    state = solve_redox_fixed_po2(cfg, P_O2=1e-9, C=0.018)
    r = state.r_Cr
    rhs = (
        state.k_Cr * state.a_Cr * state.gamma["CrO1.5"] ** 2 * state.N
        / (state.gamma["CrO"] ** 3 * state.split.n_totals["Cr"])
    )
    assert r**3 / (1.0 + r) == pytest.approx(rhs, rel=1e-8)
```

Import `solve_redox_fixed_po2` for this separate test and use the new paper test from Step 2 for Fe. Task 4 migrates the Cr identity to `state.activities.gamma_conventional_by_species`. Do not relabel the historical fixed-composition fixture as a corrected equilibrium case.

- [ ] **Step 5: Run focused and integrated tests, then commit Task 2.**

```powershell
uv run pytest -q -p no:cacheprovider
git add slag_model/data.py slag_model/equilibrium.py slag_model/redox.py slag_model/solver.py tests/test_equilibrium.py tests/test_rsm.py examples/thermodynamic_review/test_paper_validation.py
git commit -m "fix: align chromium and iron redox standard states"
```

Expected: C1/C2 source-equation tests PASS and the ordinary suite is green. Other strict `xfail` cases remain expected failures, not silently broadened tolerances.

### Task 3: Model RS pseudo-components and formula-unit conversions explicitly

**Files:**
- Modify: `slag_model/data.py` (distinct cation-to-RS-species mapping and reference-state metadata)
- Modify: `slag_model/rsm.py` (new activity value, conversion functions)
- Create: `tests/test_activities.py` (new builder equations)

**Interfaces:**
- Consumes: `rtln_gamma_rs: dict[str,float]`, `X: dict[str,float]`, `T`, and optional `ConversionSpec(A, B, standard_state)`.
- Produces: `build_slag_activities(rtln_gamma_rs, X, T, *, sio2_conversion, fe2o3_conversion, al2o3_conversion) -> SlagActivities`. The existing solver/report remains operational until Task 4 switches it atomically; `convert_gammas` is not removed in this preparatory task.

- [ ] **Step 1: Verify the conventional-reference labels from the source tables.** Transcribe the phase/reference for every current conversion coefficient. Use source-qualified labels where the paper does not define a phase; never infer liquid versus solid from the oxide name. Preserve the selected `sio2_conversion` variant in that label.

- [ ] **Step 2: Write and run failing formula-unit tests.** Create `tests/test_activities.py` with:

```python
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
    result = build_slag_activities(
        {"P5+": 0.0, "Al3+": 0.0}, {"P5+": 0.0, "Al3+": 0.1}, 1873.0
    )
    assert result.a_conventional_by_species["P2O5"] == 0.0
    assert "Al2O3" not in result.a_conventional_by_species
    assert result.a_rs_by_cation["Al3+"] == pytest.approx(0.1)
```

Run the two tests; expected FAIL because the API does not exist yet.

- [ ] **Step 3: Add the types and activity builder.** Define the precise fields, then calculate one-cation and two-cation conversions through separate branches:

```python
@dataclass(frozen=True)
class ConversionSpec:
    A: float
    B: float
    standard_state: str

@dataclass(frozen=True)
class SlagActivities:
    gamma_rs_by_cation: dict[str, float]
    a_rs_by_cation: dict[str, float]
    gamma_conventional_by_species: dict[str, float]
    a_conventional_by_species: dict[str, float]
    delta_g_conversion_J_per_mol_species: dict[str, float]
    standard_state_by_species: dict[str, str]

# Inside build_slag_activities, with rt = R*T:
gamma_rs = {c: math.exp(v / rt) for c, v in rtln_gamma_rs.items()}
a_rs = {c: gamma_rs[c] * X[c] for c in rtln_gamma_rs}
dg_p = 52720.0 - 230.706 * T
a_conventional["P2O5"] = a_rs["P5+"] ** 2 * math.exp(dg_p / rt)
```

Keep the P2O5 coefficients in `data.py` and reference them rather than duplicate literals in production. For optional Al2O3, use `a_rs["Al3+"]**2 * exp((A+B*T)/RT)` and its supplied `standard_state`; for optional FeO1.5, use the one-power formula. Do not put either species in the new conventional maps by default. Leave legacy `convert_gammas` untouched only until Task 4 migrates the solver and report in the same change.

- [ ] **Step 4: Test the builder independently.** Check the one-cation identity for FeO against the existing conversion coefficient and check custom Al2O3 with two powers and an explicit state:

```python
def test_custom_alumina_conversion_uses_two_rs_units():
    T = 1873.0
    spec = ConversionSpec(A=100.0, B=-0.1, standard_state="Al2O3(custom-solid)")
    result = build_slag_activities(
        {"Al3+": 0.0}, {"Al3+": 0.1}, T, al2o3_conversion=spec
    )
    assert result.a_conventional_by_species["Al2O3"] == pytest.approx(
        0.1**2 * math.exp((100.0 - 0.1 * T) / (R * T))
    )
    assert "Al2O3" not in result.gamma_conventional_by_species
```

- [ ] **Step 5: Run builder and existing solver tests, then commit Task 3.** The existing solver still uses its old representation in this intermediate commit, so its existing tests must remain green; the legacy P2O5 paper test remains strict-xfail until Task 4 migrates the caller.

```powershell
uv run pytest tests/test_activities.py tests/test_rsm.py tests/test_solver.py -q -p no:cacheprovider
git add slag_model/data.py slag_model/rsm.py tests/test_activities.py
git commit -m "feat: calculate explicit RS and conventional slag activities"
```

Expected: new builder equations PASS and the existing solver/report tests remain green. C3 is not declared fixed until Task 4.

### Task 4: Publish a versioned, unambiguous v2 report

**Files:**
- Modify: `slag_model/solver.py` (`SolveResult.activities`, explicit activity consumers)
- Modify: `slag_model/io.py` (`Options`, custom conversion parsing, `build_report`)
- Modify: `slag_model/cli.py` (component table)
- Modify: `tests/test_io.py`, `tests/test_solver.py`, and `tests/test_rsm.py` (schema and migrated fields)
- Modify: `examples/thermodynamic_review/test_paper_validation.py` (migrate P2O5 and Cr paper tests; remove C3 `xfail`)
- Modify: `README.md` (v2 JSON contract and migration note)

**Interfaces:**
- Consumes: `build_slag_activities(...)` and `ConversionSpec` from Task 3.
- Produces: `SolveResult.activities: SlagActivities`; `build_report(cfg, result) -> dict` with `schema_version == 2`; custom conversion input objects use `{"A": number, "B": number, "standard_state": nonempty string}`.

- [ ] **Step 1: Add failing schema tests.** Use the reference input and call `solve`/`build_report`:

```python
def test_report_v2_distinguishes_rs_and_formula_unit_activities():
    cfg = parse_input(json.dumps(VALID_JSON))
    report = build_report(cfg, solve(cfg))
    assert report["schema_version"] == 2
    rows = {row["cation"]: row for row in report["solution"]["components"]}
    assert rows["P5+"]["rs_species"] == "PO2.5"
    assert rows["P5+"]["conventional"]["species"] == "P2O5"
    assert rows["P5+"]["conventional"]["rs_units_per_species"] == 2
    assert "gamma" not in rows["P5+"]["conventional"]
    assert rows["Al3+"]["rs_species"] == "AlO1.5"
    assert rows["Al3+"]["conventional"] is None
    assert "gamma" not in rows["P5+"]
    assert "a" not in rows["P5+"]
```

Also test `as_excel` Ti row has `model_status == "denominator_only_legacy"` and no activity. Run tests; expected FAIL against v1.

- [ ] **Step 2: Migrate solver consumers and the paper test.** Use `activities.a_conventional_by_species["FeO"]` for oxygen potential; `activities.gamma_rs_by_cation["Fe2+"]` and `["Fe3+"]` for Fe redox; `activities.gamma_conventional_by_species["CrO"]` and `["CrO1.5"]` for Cr equilibrium. Replace `SolveResult.gamma`, `a_slag`, and `delta_g_conv` with `SolveResult.activities`. The old `convert_gammas` must now be removed or reject multi-cation `P5+`/`Al3+` calls with an explicit migration error. Replace the C3 strict-xfail paper test with the new builder formula assertion; migrate Xiao's experimental activity test to `a_conventional_by_species["CrO"]` and `["CrO1.5"]` without changing published inputs or the declared envelope.

- [ ] **Step 3: Make custom conversion input explicit.** Parse `ConversionSpec` only when `A`, `B`, and a nonempty `standard_state` are present; serialize the same object shape through `Options.as_dict()`. Tests must reject `{"A": 100, "B": -0.1}` for both optional conversions with a migration message, and accept/round-trip the three-field form. The pre-existing strict-xfail round-trip test becomes a passing test if this change resolves it.

- [ ] **Step 4: Build v2 rows from activities without new calculations.** For each modeled cation, emit the RS fields and, when a conventional activity exists, a nested object:

```python
species = CATION_TO_OXIDE[cation]
activities = result.activities
if species in activities.a_conventional_by_species:
    reason = None
    conventional = {
        "species": species,
        "standard_state": activities.standard_state_by_species[species],
        "rs_units_per_species": 2 if species in {"P2O5", "Al2O3"} else 1,
        "DeltaG_conversion_J_per_mol_species": (
            activities.delta_g_conversion_J_per_mol_species[species]
        ),
        "a": activities.a_conventional_by_species[species],
    }
    if species in activities.gamma_conventional_by_species:
        conventional["gamma"] = activities.gamma_conventional_by_species[species]
        conventional["gamma_fraction_basis"] = "X_cation"
else:
    conventional = None
    reason = "no documented conversion"

row = {
    "cation": cation,
    "rs_species": CATION_TO_RS_SPECIES[cation],
    "X_cation": result.X[cation],
    "RTln_gamma_RS_J_per_mol_cation": result.rtln_gamma_rs[cation],
    "gamma_RS": result.activities.gamma_rs_by_cation[cation],
    "a_RS": result.activities.a_rs_by_cation[cation],
    "conventional": conventional,
    "conventional_unavailable_reason": reason,
}
```

The nested object contains `species`, `standard_state`, `rs_units_per_species`, `DeltaG_conversion_J_per_mol_species`, and `a`; add `gamma` and `gamma_fraction_basis: "X_cation"` only for one-cation species. `P2O5` and optional `Al2O3` have two RS units and no `gamma` field. Add `schema_version: 2` at the report root. Print the same distinction in the CLI table; do not print zero for unavailable values.

- [ ] **Step 5: Update documentation and run focused tests.** Document the v1-to-v2 field migration and explicitly label all report standard states, including source-qualified references where phase was not stated. Run:

```powershell
uv run pytest tests/test_activities.py tests/test_io.py tests/test_solver.py tests/test_rsm.py examples/thermodynamic_review/test_paper_validation.py -q -p no:cacheprovider
git add slag_model/solver.py slag_model/io.py slag_model/cli.py tests/test_io.py tests/test_solver.py tests/test_rsm.py examples/thermodynamic_review/test_paper_validation.py README.md
git commit -m "feat: publish explicit activity report schema v2"
```

Expected: report and input round-trip tests PASS, with no legacy `gamma`/`a` rows.

### Task 5: Reject non-finite values at every boundary

**Files:**
- Modify: `slag_model/io.py` (`parse_input`, `load_input`, `build_report`/serialization helper)
- Modify: `slag_model/solver.py` (finite checks and named failure)
- Modify: `slag_model/cli.py` (error handling, strict pre-serialization)
- Modify: `tests/test_io.py` (input cases)
- Create: `tests/test_cli.py` (no invalid/partial report)
- Modify: `examples/thermodynamic_review/test_defect_regressions.py` (remove C4 `xfail`)

**Interfaces:**
- Consumes: scalar input fields, per-iteration state/targets/changes, and report dict.
- Produces: `InputError` for malformed/non-finite input; `NumericalStateError` for invalid calculation; strict JSON text or a controlled CLI error.

- [ ] **Step 1: Add failing input tests.** Parameterize bare constants, quoted strings, huge integers, composition values, metal values, and custom conversion values:

```python
@pytest.mark.parametrize("bad", ["NaN", "Infinity", "-Infinity"])
def test_nonfinite_pressure_string_rejected(bad):
    payload = {**VALID_JSON, "P_CO_atm": bad}
    with pytest.raises(InputError, match="P_CO_atm.*finite"):
        parse_input(json.dumps(payload))

def test_bare_json_nan_rejected():
    with pytest.raises(InputError, match="non-finite"):
        parse_input(json.dumps(VALID_JSON).replace('"P_CO_atm": 1.0', '"P_CO_atm": NaN'))
```

Add corresponding tests for `slag_wtpc.SiO2`, `metal_wtpc.Cr`, custom `A/B`, and an integer large enough to overflow `float`. Remove the C4 strict `xfail` marker and verify these fail before code changes.

- [ ] **Step 2: Normalize and validate numeric input once.** Use a helper that rejects bools, catches `TypeError`, `ValueError`, and `OverflowError` from `float`, then checks `math.isfinite` before sign checks. Preserve existing finite numeric-string support only for temperature and CO pressure. Pass `parse_constant` to `json.loads` so bare `NaN`/`Infinity` becomes `InputError`. Check finite composition sums after addition. Catch malformed JSON as `InputError` while preserving location context.

- [ ] **Step 3: Add a failing false-convergence test.** Monkeypatch `solver.compute_state` to return an otherwise valid state plus targets `{"r_Fe": 0.10, "r_Cr": 1.0, "C": float("nan")}`. `solve(cfg)` must raise `NumericalStateError`, never return a result merely because the NaN is last in `max(changes.values())`.

- [ ] **Step 4: Guard iteration and result values.** Add `NumericalStateError(RuntimeError)` and a single finite checker over input iterates, the `SolveResult` numeric scalars/maps, targets, relative changes, and damped updates. Use it before convergence comparison and before returning. Trap arithmetic overflow/domain errors from `compute_state` and rethrow with iteration and field context; do not turn them into convergence. Apply the same guard to `solve_redox_fixed_po2` without changing its normal-priority diagnostic semantics.

- [ ] **Step 5: Serialize before opening an output file.** In `cli.py`, call `json.dumps(report, indent=2, allow_nan=False)` before displaying success or opening `args.output`. Catch serialization `ValueError` as a controlled numerical/report error; write only already-valid text. Add a CLI test that injects `NaN` into a report and checks nonzero exit and no newly created output file.

- [ ] **Step 6: Run safety and CLI tests, then commit Task 5.**

```powershell
uv run pytest tests/test_io.py tests/test_cli.py examples/thermodynamic_review/test_defect_regressions.py -q -p no:cacheprovider
git add slag_model/io.py slag_model/solver.py slag_model/cli.py tests/test_io.py tests/test_cli.py examples/thermodynamic_review/test_defect_regressions.py
git commit -m "fix: reject non-finite thermodynamic inputs and states"
```

Expected: invalid values fail with named errors; normal-priority strict `xfail` cases remain expected failures unless incidentally fixed and unmarked.

### Task 6: Independent integration acceptance and handoff

**Files:**
- Modify: `tests/test_solver.py` (source-consistent integration sentinel)
- Modify: `README.md` (corrected equations and safety statement)
- Modify: `examples/thermodynamic_review/REVIEW.md` (mark C1-C4 resolved with verification evidence; retain N/low findings)

**Interfaces:**
- Consumes: corrected Tasks 1-5 and the unchanged current parameter matrix.
- Produces: reproducible full-suite evidence and a clear record of what remains outside P0 scope.

- [ ] **Step 1: Add a reference integration test.** Use `examples/eaf_slag_01.json` and compare the coupled solve against the independently derived C1+C2 sentinel under the unchanged matrix:

```python
def test_reference_case_after_p0_corrections():
    cfg = load_input(SLAG_PATH)
    result = solve(cfg)
    assert result.k_Cr == pytest.approx(0.00878747397894681, rel=1e-10)
    assert result.r_Fe == pytest.approx(0.13189054746, rel=2e-5)
    assert result.r_Cr == pytest.approx(0.25229986894, rel=2e-5)
    assert result.C_wtpc == pytest.approx(0.0319280093, rel=2e-5)
    assert result.P_O2_atm == pytest.approx(3.45717e-10, rel=2e-5)
```

The sentinel is specific to the existing hybrid matrix; a later named parameter-set change must version it, not silently adjust tolerance.

- [ ] **Step 2: Check algebraic invariants independently.** Recalculate Fe and Cr cation totals from `split.n_cations`, check `sum(X)` for the default mode, and calculate the Fe, CO, and Cr $Q/K$ values directly from the returned activities and `k` values. These identity checks supplement, but do not replace, the paper tests.

- [ ] **Step 3: Run complete verification from the repository root.**

```powershell
uv run pytest -p no:cacheprovider -v
uv run ruff check --no-cache .
uv run ruff format --check --no-cache .
git status --short
```

Expected: all ordinary and four formerly P0 `xfail` tests PASS; remaining strict `xfail` tests correspond only to documented normal findings. Do not claim a pass without the command output.

- [ ] **Step 4: Update the audit and commit only the P0 handoff.** Record exact test counts and any changed baseline values. Keep N1-N12 open, describe the workbook as legacy, and identify the exact report-schema migration. Commit only files in this plan if commits are authorized in the execution environment.

## Execution handoff

This is a plan, not implementation. Review the design and this plan before execution. The tasks share evolving `SlagActivities` and report interfaces, so **native task-by-task execution** is the least wasteful default; use a fresh whole-branch review at the end because a mistaken standard state would be costly. If independent task reviewers are preferred despite the interface dependency, use the subagent-driven method specified in the header. No production files were changed while writing this plan.
