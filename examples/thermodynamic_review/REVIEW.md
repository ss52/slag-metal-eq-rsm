# Thermodynamic and code audit

Date: 2026-09-23  
Repository: `slag-activity-model` 0.1.0  
Scope: equations, standard states, source-data transcription, numerical solver,
input/output behavior, reference workbook, existing tests, and paper-derived
validation cases.

## Executive conclusion

The implementation is structurally clear and most of the basic algebra is sound,
but the current coupled result is **not yet scientifically reliable**. Three
thermodynamic standard-state/stoichiometry defects materially change results, and
one input-validation defect permits a false successful solve containing `NaN`.

The ordinary project suite passes (`56 passed`), but it does not detect these
problems. Several expectations are derived from the supplied workbook, which
contains source-data and formula errors of its own. The library should not be used
for process decisions until the critical findings below are fixed and the strict
`xfail` tests in this directory are converted to passing assertions.

Correcting only the two coupled redox defects, while leaving the rest of the code
unchanged, changes the reference case as follows:

| Quantity | Current code | Source-consistent Fe/Cr calculation |
|---|---:|---:|
| $K_{Cr}$ at 1823.15 K | 0.915831 | 0.00878747 |
| $n_{Fe^{3+}}/n_{Fe^{2+}}$ | 0.173117 | 0.131891 |
| $n_{Cr^{2+}}/n_{Cr^{3+}}$ | 1.53820 | 0.252300 |
| dissolved C, wt% | 0.0326668 | 0.0319280 |
| $P_{O_2}$, atm | $3.30099\times10^{-10}$ | $3.45717\times10^{-10}$ |

This comparison is an independent calculation made with the source equations; it
is not a proposed new empirical fit.

## What is correct

The following parts agree with the cited equations or satisfy independent balance
checks:

- Oxide-to-cation mole conversion and the Fe/Cr total-cation mass balances.
- Normalization of the active cation fractions when Ti is excluded.
- The quadratic regular-solution expression in `slag_model/rsm.py:126-146`.
- The algebra converting the chromium equilibrium into
  $r_{Cr}^3/(1+r_{Cr})=\mathrm{RHS}$, conditional on a correct $K_{Cr}$ and
  consistent standard states.
- The FeO-based oxygen-potential equation and the C-CO equation, considered in
  isolation.
- Ban-ya Table 3 coefficients for FeO, CaO, MgO, MnO, and P2O5 in the Python
  data table.
- Xiao, Holappa, and Reuter Table IV coefficients for CrO, CrO1.5, and SiO2.
- Five chromium-related Table III interaction energies used by the Python matrix.
- The eight selected Xiao 2002 experimental compositions remain within a factor
  2.1 of the reported CrO and CrO1.5 activities. This deliberately broad envelope
  reflects the paper's stated scatter and approximate model agreement; it is not
  a high-precision acceptance criterion.

Primary sources:

- [Ban-ya 1993](https://www.jstage.jst.go.jp/article/isijinternational1989/33/1/33_1_2/_article)
- [Xiao and Holappa 1993](https://www.jstage.jst.go.jp/article/isijinternational1989/33/1/33_1_66/_article)
- [Xiao, Holappa, and Reuter 2002](https://doi.org/10.1007/s11663-002-0039-9)
- [Sigworth and Elliott 1974](https://doi.org/10.1179/msc.1974.8.1.298)

## Critical findings

### C1. Chromium equilibrium uses the wrong energy units and the wrong metal standard state

Evidence:

- `slag_model/data.py:170` stores `(25690, -13.36)` under a comment that calls
  the reaction `2 CrO1.5 + [Cr] = 3 CrO`.
- `slag_model/equilibrium.py:34` treats those coefficients as J/mol.
- `slag_model/solver.py:102-107` supplies $a_{Cr}=f_{Cr}[\%Cr]$.

The published reaction is

$$2\,CrO_{1.5}+Cr(s)=3\,CrO,$$

with

$$\Delta G^\circ=(25690-13.36T)\ \mathrm{cal\ mol^{-1}},$$

and pure solid chromium as the metal standard state. The code therefore makes two
independent changes without conversion: cal is interpreted as joules and $Cr(s)$
is replaced by dissolved chromium on the 1 mass-% standard state.

For the library's metal activity convention, the required reaction energy is

$$\Delta G^\circ_{1\%Cr}
=4.184(25690-13.36T)-(19246-46.86T)\ \mathrm{J\ mol^{-1}},$$

where the second term is $Cr(s)=[Cr]_{1\%}$. At 1823.15 K this gives
$K_{Cr}=0.00878747$, versus 0.915831 in the library. The reference-case chromium
ratio consequently changes from 1.538 to 0.252 when this and C2 are corrected.

Required action: represent the reaction and its metal standard state explicitly;
convert calories to joules at the data boundary; add a source-equation test for
the compound reaction.

### C2. Ban-ya Fe redox equation is called with mixed standard states

Evidence:

- `slag_model/solver.py:92` correctly uses converted, conventional FeO activity
  for the Fe-FeO oxygen buffer.
- `slag_model/solver.py:98` then passes that converted FeO coefficient to Ban-ya
  Eq. 24, while FeO1.5 remains on the regular-solution scale by default.
- The same problem occurs in `slag_model/solver.py:205-210`.

Ban-ya Eq. 24 uses both $\gamma_{FeO(RS)}$ and
$\gamma_{FeO_{1.5}(RS)}$. The FeO conversion contributes about a factor 1.344 at
1823.15 K, so mixing the standards is not a cosmetic labeling issue. With the
same remaining assumptions, correcting the caller changes the reference-case
$r_{Fe}$ from 0.173117 to 0.131891.

Required action: preserve RS coefficients as first-class values and pass them to
the Fe redox equation; use converted FeO only where the FeO(l) standard is
actually required.

### C3. P2O5 and Al2O3 activities are calculated as one-cation pseudo-components

Evidence:

- `slag_model/data.py:135-136` maps $Al^{3+}$ directly to `Al2O3` and $P^{5+}$
  directly to `P2O5`.
- `slag_model/rsm.py:179-184` applies one conversion exponent.
- `slag_model/solver.py:114-117` reports every oxide activity as $\gamma_iX_i$.

Ban-ya's cation model components are $AlO_{1.5}$ and $PO_{2.5}$. For example,
Table 3 gives

$$P_2O_5(l)=2\,PO_{2.5}(RS),$$

so the formula-unit activity is

$$RT\ln a_{P_2O_5}
=2RT\ln\!\left(\gamma^{RS}_{PO_{2.5}}X_P\right)
+\Delta G^\circ_{conv}.$$

The one-power expression now returned as `a_P2O5` is neither this conventional
formula-unit activity nor clearly labeled pseudo-component activity. The same
stoichiometric problem applies to Al2O3. In addition, the Al2O3 conversion is not
known in the project, so a conventional Al2O3 activity cannot currently be
reported at all.

Required action: expose pseudo-component activities as `a_PO2.5_RS` and
`a_AlO1.5_RS`; calculate formula-unit activities only when a documented
conversion is available; never label an RS pseudo-component as conventional
P2O5 or Al2O3.

### C4. Non-finite input can falsely converge and serialize invalid JSON

Evidence:

- `slag_model/io.py:151-155` checks type and sign but not `math.isfinite`.
- `slag_model/io.py:169-172` accepts strings convertible to `NaN` or infinity.
- `slag_model/solver.py:171-174` can declare convergence because Python's `max`
  can ignore a trailing `NaN` depending on ordering.
- `slag_model/cli.py:104` uses the default JSON encoder, which emits non-standard
  `NaN` tokens.

Reproduction: setting `P_CO_atm` to `"NaN"` returns a nominally converged result
after 22 iterations with `C_wtpc=NaN` and `Q/K_CO=NaN`.

Required action: reject every non-finite scalar at input, check all state and
residual values for finiteness on every iteration, and serialize reports with
`allow_nan=False`.

## Normal-priority findings

### N1. One Ban-ya interaction energy is mistyped

`slag_model/data.py:40-51` uses
$\alpha_{Fe^{3+},Ca^{2+}}=-96810$ J/mol. Ban-ya Table 2 gives
$-95810$ J/mol. The same typo is present in the reference workbook.

### N2. The global interaction matrix is an undocumented hybrid dataset

The interaction matrix uses the Ban-ya parameter set, supplemented only with the
chromium terms from Xiao 2002 Table III. The Ca-Si value is therefore intentionally
the Ban-ya value $-133890$ J/mol rather than Xiao's optimized $-139100$ J/mol. This
extended Ban-ya parameter set should be named, versioned, and identified in the
calculation report so that results are reproducible.

### N3. The Wagner matrix appears to contain a reciprocal-coefficient transcription

`slag_model/data.py:179` uses $e_{Mn}^{C}=-0.012$. Source compilations traceable to
Sigworth and Elliott give approximately $e_{Mn}^{C}=-0.07$ at 1873 K; $-0.012$
is the commonly tabulated $e_C^{Mn}$. Wagner coefficients are not generally
symmetric. This has a small effect for the current low-C/low-Mn reference case,
but it is a source-data error and becomes material outside that composition.

### N4. Chromium root bracketing and redox splitting lose trace species

`slag_model/redox.py:38-46` hard-clamps the positive root to
$[10^{-6},10^6]$. For RHS $=10^{-24}$ it returns $10^{-6}$ although the root is
approximately $10^{-8}$; the reconstructed RHS is wrong by about six orders of
magnitude. Separately, `slag_model/rsm.py:57-58` and `:65-66` compute a trace
species by subtracting two nearly equal floats, which can round it to zero.

Required action: use an adaptive bracket or a monotone analytic/numerical method,
and compute both split species directly from stable ratio formulas.

### N5. Fixed-$P_{O_2}$ diagnostics are inconsistent with the reported state

`solve_redox_fixed_po2` calculates all diagnostics using the internally derived
oxygen pressure, then overwrites only `P_O2_atm` at `slag_model/solver.py:229`.
For the workbook case the report says Fe-O $Q/K=1$ and CO $Q/K=1.902$, while
recalculation at the reported $10^{-9}$ atm gives about 0.551 and 1.048.

Required action: recompute the complete final state and diagnostics under the
fixed external pressure, or clearly omit equations that are not enforced.

### N6. `ti_handling="as_excel"` is denominator-only dilution, not zero-interaction Ti

`slag_model/rsm.py:94-105` adds Ti to the cation denominator, but the activity
coefficient calculation has no Ti row/column; `slag_model/solver.py:88-90` forces
$\gamma_{TiO_2}=1$. A true zero-$\alpha$ Ti component would still receive cross
terms in the multicomponent regular-solution expression. Ban-ya also publishes
several nonzero Ti interactions, so the current behavior is only workbook
emulation and should not be presented as a physical model.

Required action: rename this option to make the emulation explicit, or add Ti as a
fully modeled cation with documented parameters.

### N7. Input and report contracts are not closed under round-trip

- `Options.as_dict()` serializes custom conversion pairs as arrays
  (`slag_model/io.py:54-64`), but `_parse_custom_conversion` accepts only objects
  with `A` and `B` (`slag_model/io.py:81-91`).
- `options: []` is silently treated as default options because
  `raw = raw or {}` precedes the type check (`slag_model/io.py:95`).
- Malformed JSON and null numeric fields raise raw `JSONDecodeError`/`TypeError`,
  while the CLI catches only `InputError` and `NonConvergenceError`.

Required action: define one canonical JSON schema, validate against it, and test
`parse(report["input"])` as an invariant.

### N8. Element-total inputs trigger a false 100 wt% warning

`slag_model/io.py:201-204` sums the raw input values even when Fe and Cr are given
as elemental totals. A composition exactly equivalent to the reference oxide
input is therefore reported as roughly 94.97 wt%. Convert analytical forms to a
common mass basis before the check, or define separate validation rules for the
two input conventions.

### N9. The workbook is a regression artifact, not an authoritative oracle

The supplied workbook contains, among others:

- Fe redox based on converted rather than RS FeO gamma (`Balance!C11`).
- the chromium expression `25690-13.36T` labeled J/mol (`Balance!C33`);
- a P2O5 coefficient `-23.706` instead of `-230.706` (`Slag!B65:C65`);
- MnO conversion cells containing the MgO coefficients (`Slag!B64:C64`);
- an Fe2+-Al3+ displayed matrix entry `-410` while its cross-term table uses
  `-41000` (`Slag!M5:W15`).

The Python code corrected some of these workbook errors but retained others. Tests
must therefore trace expectations to papers, not merely to the spreadsheet.

### N10. Existing acceptance tests are too permissive and partly self-referential

- The full-model reference $r_{Fe}$ is 0.173, while the plan's stated range was
  0.08-0.11; `tests/test_solver.py` accepts 0.05-0.25.
- `tests/test_rsm.py` allows an 800 J/mol tolerance for Al3+ and its comment says
  the expected value uses full cross terms, although the actual full-cross-term
  value differs much more.
- $Q/K=1$ tests for FeO and CO mostly confirm that the code inverted its own
  equations; they are not independent validation.

Required action: separate identity tests, source-transcription tests, independent
paper cases, and empirical acceptance envelopes. Never widen a tolerance to make
an unexplained discrepancy pass.

### N11. Applicability gaps are not surfaced to the caller

Fe-Cr, Mn-Cr, and P-Cr slag interaction energies are zero because data are absent,
not because zero interaction has been established. The approximation matters in a
Fe-rich chromium-bearing slag, but no runtime warning or report metadata identifies
it. Reports should include parameter-set identity, missing-pair warnings, source
references, standard states, and model-validity ranges.

### N12. The bundled HSC PDFs do not validate this model

`Docs/EQ_HSC.pdf` and `Docs/EQ_HSC_examples.pdf` are general HSC Equilibrium
module manuals. They explain Gibbs minimization and generic solution-model use,
but they do not supply an independent EAF slag/hot-metal case or the constants
used here. They are useful background documentation, not a validation source.

## Low-priority findings

Per the requested scope, these were recorded but not investigated further:

- The plan says a default CLI run writes `report.json`; the implementation writes
  only when `-o` is supplied.
- A few names and comments conflate cation pseudo-components with oxide formula
  units even where the numerical result is otherwise unambiguous.
- Documentation still calls the Xiao 2002 SiO2 conversion a workbook-specific or
  uncertain value even though Table IV provides its source.

## Paper-derived test pack

Files:

- `paper_cases.json`: source citations, exact constants, equations, and eight
  experimental cases from Xiao 2002 Table I.
- `test_paper_validation.py`: exact transcription/equation tests plus the declared
  experimental activity envelope.
- `test_defect_regressions.py`: executable reproductions for non-paper numerical
  and input/output defects.

Run from the repository root:

```powershell
uv run pytest examples\thermodynamic_review -p no:cacheprovider -v
```

Current result:

```text
21 passed, 14 xfailed
```

The strict `xfail` markers are intentional. Each corresponds to a known defect;
after fixing that defect, remove its marker so an unexpected pass becomes a normal
required pass. The experimental factor-of-2.1 envelope is not a substitute for the
strict equation and data tests.

## Recommended repair order

1. Fix C1 and explicitly model every reaction's units and standard states.
2. Fix C2 by keeping RS and conventional gammas distinct in the type/API design.
3. Fix C3 by separating cation pseudo-component and formula-unit activities.
4. Add finite-number invariants at input, iteration, result, and JSON boundaries.
5. Introduce named/versioned parameter sets and correct the two confirmed data
   transcription errors.
6. Fix the root bracket, trace-species arithmetic, fixed-pressure diagnostics, and
   JSON round-trip contract.
7. Rebuild the normal test suite around the paper fixtures; retain the workbook
   only as a clearly labeled legacy-regression profile.

## Verification performed

- `uv run pytest -p no:cacheprovider -v` -> 56 passed.
- `uv run ruff check --no-cache .` -> passed before audit artifacts were added.
- `uv run ruff format --check --no-cache .` -> 21 files already formatted before
  audit artifacts were added.
- New audit pack -> 21 passed, 14 expected failures.
- The reference case was recalculated independently with source-consistent Fe and
  Cr equations.
- The supplied workbook was inspected as formulas and cached values, not only as
  rendered output.
- Relevant pages of both bundled HSC PDFs were rendered and visually checked.

No production code was changed in this audit.
