# Thermodynamic and code audit

## Current scientific status (2026-10-05)

This section supersedes the open-finding and applicability claims in the dated
audit below. It records the source/profile corrections and the Task 7 output
contract; final integrated test counts are pending and are not asserted here.
For active equations, report fields and exact tested inputs, use the current
[`README.md`](../../README.md) and [`paper_cases.json`](paper_cases.json).

### Critical

The earlier C1–C4 corrections address chromium reaction units and the dissolved
Cr standard state, consistent RS coefficients in Ban-ya Fe redox, formula-unit
activities for P₂O₅/Al₂O₃, and non-finite input/output handling. These are distinct
from empirical model validation. Conventional Al₂O₃ and FeO₁.₅ activities remain
unavailable without labeled custom conversions; their RS activities are
available. No new critical scientific finding is identified in the present
source assessment at the tested points.

### Normal

- **Source profiles and evidence.** The default `banya93_xiao95_v1` uses the
  corrected Ban-ya Fe³⁺–Ca²⁺ coefficient −95,810 J/mol and Xiao and Holappa
  (1995) Ca²⁺–Si⁴⁺ coefficient −139,100 J/mol. The separate
  `banya93_casi_xiao95cr_hybrid_v1` restores Ban-ya's −133,890 J/mol Ca–Si
  value while retaining Xiao chromium terms. `legacy_workbook_typo_hybrid_v1`
  preserves the historical workbook matrix, including −96,810 J/mol, and is
  explicitly legacy. Reports carry profile identity/version; the old hybrid
  and transcription findings are superseded by these named choices.
- **Mn–C source correction.** Task 7 uses $e_{Mn}^{C}=-0.07$ while retaining
  $e_C^{Mn}=-0.012$. The directly checked source is
  [Lei, Fu and Xiong (2024), Table 1, p. 112](https://journals.pan.pl/Content/131630/AFE%202_2024_13-Final.pdf),
  whose reference [22] is Chen's 2010 steelmaking databook. The original
  Sigworth and Elliott (1974) table was not directly verified for this value;
  Wagner coefficients are asymmetric.
- **Mode-consistent diagnostics.** Task 7 identifies `solution.calculation_mode`
  as `coupled` or `fixed_P_O2`. Every $Q/K$ is evaluated at the reported
  pressure. Coupled mode enforces FeO, CO, Cr comproportionation and Fe-redox
  relations; fixed mode imposes pressure and carbon, enforcing only Cr
  comproportionation and Fe redox. Reports list these under
  `equilibrium.enforced_relations`; fixed-mode FeO/CO quotients are diagnostics,
  not enforced equilibrium claims.
- **Chromium–oxygen scope qualification — documented.** The implemented Cr
  relation is Xiao and Holappa (1995) Eq. 9 comproportionation, not the independent
  $CrO_{1.5}(s)=CrO(l)+\tfrac14 O_2(g)$ equilibrium in
  [p. 324, Eqs. 14–15](https://www.pyrometallurgy.co.za/InfaconVII/319-Xiao.pdf).
  Its valence constraint need not agree with the Fe-buffer oxygen potential.
  In the default reference case, the reported Fe-buffer pressure is about
  $3.72\times10^{-10}$ atm while Eq. 15 with the same oxide reference states
  and calculated activities implies $9.06\times10^{-10}$ atm (2.44 times
  higher), giving Cr–O $Q/K\approx0.80$ at the reported pressure. This is an
  unenforced source-equation discrepancy, not roundoff or empirical fit scatter;
  it is neither a measured validation result nor an accuracy tolerance.
  `Q_over_K_Cr` checks only the implemented comproportionation identity.
  The overbroad common-oxygen redox scope claim is corrected in the
  documentation; no additional model relation, calibration or program behavior
  has been introduced. The implemented numerical relations remain unchanged,
  and the calculation does not establish complete slag–metal redox equilibrium.
- **Approximate oxygen estimate.** The selected iron activity is used in
  $[\%O]_{FeO\ equivalent}=(a_{FeO}/a_{Fe})10^{-6320/T+2.734}$, assuming
  $f_O=1$ and retaining an inherited empirical saturation fit. Coupled mode
  reports this approximate value as both `O_wtpc` and
  `O_FeO_equivalent_wtpc`. Fixed mode reports `O_wtpc=null` and labels the
  separate estimate `FeO_equivalent_only`; it is not dissolved oxygen inferred
  at the imposed pressure. The auxiliary FeO/CO, carbon-dissolution and
  oxygen-saturation fits have incomplete directly verified provenance and do
  not come from the two oxide papers.
- **Prescribed composition and validation scope.** The solver fixes slag
  elemental Fe/Cr totals, metal Cr/Mn/P, temperature and $P_{CO}$, then solves
  valence ratios and carbon. It has no phase masses, slag/metal mass ratio,
  gas inventory or coupled element transfer. Selected equilibrium identities
  therefore do not establish closed-system equilibration or Mn/P partition
  coefficients. The reference EAF assay at 1823.15 K and $P_{CO}=1$ atm,
  with Cr 0.038, Mn 0.007 and P 0.069 wt% in the metal and solved C about
  0.0308 wt%, is computational regression evidence. The exact assay is in
  `examples/eaf_slag_01.json`; it is not an experimental full-solver case.
- **Paper comparison scope.** The active oxide sources are Ban-ya (1993) and
  Xiao and Holappa (1995). Ban-ya's chromium-free equimolar CaO–SiO₂ check at
  1823 K, with −133,890 J/mol Ca–Si and β-cristobalite conversion, gives
  equation-derived $a_{SiO_2}=0.2575019792$ and $a_{CaO}=0.01103118178$.
  No complete source-exact measured Ban-ya 1993 row is available in this pack.
  All eleven Xiao Table 1 SiO₂–CrO–CrO₁.₅ activity comparisons at 1873 K use
  observed valence splits and parameter-assessment data, so they are in-sample.
  Mean absolute relative differences remain 20.4% CrO, 12.6% CrO₁.₅ and
  37.4% against Gibbs–Duhem-derived SiO₂; T1-10 CrO is −54.67%. No empirical
  acceptance envelope or complete coupled-model validation is established.
- **Applicability limitations.** Ban-ya's approximate liquid-slag model has
  composition-dependent applicability; this implementation does not predict
  solid phases, spinel or phase separation. Xiao chromium parameters cover
  studied Ca–Si–Mg–Al–Cr systems; missing Fe–Cr, Mn–Cr and P–Cr interactions
  are assumed zero, so Fe-rich chromium slags are extrapolations. Ti exclusion
  and renormalization is an approximation; `as_excel` is denominator-only
  dilution without Ti activity. The published nonzero Ti interactions are not
  used and no universal safe Ti wt% limit is supported. First-order metal
  coefficients at 1873 K are approximately applied at 1823.15 K; carbon-rich
  hot metal, including 3–5 wt% C, remains unvalidated. These tested points and
  source systems do not define a validated temperature/composition box. This
  disclosure addresses the scientific scope of archived N11; adding separate
  runtime applicability warnings is optional deferred work.
- **Deferred numerical findings.** The chromium ratio root uses a
  $10^{-6}$–$10^6$ bracket and subtraction can lose trace species near
  roundoff. Extreme roots and trace fractions outside the tested examples
  require further tests. These are deferred normal findings under the
  scientific release scope, not release blockers for the stated tested points.
- **Open ordinary input warning.** Element-total assays can still trigger a
  false raw 100 wt% sum warning (archived N8). This affects a warning, not the
  calculation, and is deferred under the focus on meaningful calculated
  results. It is not claimed fixed or treated as a tested-point release
  blocker, unlike the calculation-output errors addressed by Task 7.

### Low

The original audit and its dated verification results below are retained as
history. Their line numbers, test counts, earlier matrices, Ban-ya 1985/Xiao
2002 active-benchmark claims and repair order do not describe the current pack.
The remaining archived silica-polymorph label and other maintenance notes are
low priority and are not evidence of empirical accuracy. No human-prose tests
or new numerical acceptance threshold are introduced by this documentation
update; final integration verification is recorded separately.

## Source-benchmark correction (2026-09-25)

This note supersedes the benchmark-source and validation claims in the dated
audit below. The active pack uses Ban-ya 1993 for equations, interaction energies
and conversions, and all eleven Xiao and Holappa 1995 Table 1 cases. Ban-ya 1993
has no source-exact measured numerical row here. The Xiao Table 1 rows were used
in assessing the paper's Table 2 parameters, so calculated comparisons are
in-sample, not independent validation. The paper's CrO and CrO1.5 activities
come from EMF/oxygen-potential measurements; its SiO2 activities were calculated
by Gibbs-Duhem. Current mean absolute relative differences are 20.4% for CrO,
12.6% for CrO1.5, and 37.4% versus those Gibbs-Duhem SiO2 values; T1-10 CrO is
54.67% below the reported activity. Xiao 1995 p. 321 also notes that the liquid
CrO formation relation was extrapolated from 1665-1750 C to the 1600 C
measurements.

Everything below this notice is an archived audit snapshot kept for history.
Its older source-data, benchmark, interaction-matrix, and open-finding claims
are superseded; use the current README and `paper_cases.json` for active
parameter profiles and source provenance.

Date: 2026-09-23  
Repository: `slag-activity-model` 0.1.0  
Scope: equations, standard states, source-data transcription, numerical solver,
input/output behavior, reference workbook, existing tests, and paper-derived
validation cases.

## Executive conclusion

At the 2026-09-23 audit baseline, three thermodynamic standard-state or
stoichiometry defects materially changed results, and one input-validation defect
permitted a false successful solve containing `NaN`. P0 Tasks 1-6 have since
resolved C1-C4; each resolution and its source-based test evidence is recorded
below. The ordinary suite and paper-derived pack were rerun for the Task 6
handoff. Normal- and low-priority findings remain outside P0 scope.

Correcting only the two coupled redox defects, while leaving the rest of the code
unchanged, changes the reference case as follows:

| Quantity | Audit baseline before P0 | Source-consistent Fe/Cr calculation |
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
- Four preselected Ban-ya 1985 Table 1 rows (101, 301, 501, and 701) are within
  the declared 10% project activity screen. This threshold is not source
  measurement uncertainty. Row 901 is reported only as a stress diagnostic
  because its total iron cation fraction, 0.84236, is outside the stated range.
- Xiao 2002 Table I measured activities reconstruct from the printed activity
  coefficients and compositions within 0.006 absolute. Seven printed rows give
  pure-solid-Cr $Q/K=0.973$–$1.026$; CSC7 gives 1.769 and remains flagged as a
  source-row anomaly. The current model's mean absolute relative differences are
  29.6% for CrO and 39.6% for CrO1.5. The paper gives no quantitative
  activity-error acceptance bound; Fig. 11 shows substantial model-to-measurement
  scatter, so these values are contextual comparisons, not proof of a code defect
  or an acceptance test.

Primary sources:

- [Ban-ya et al. 1985, Table 1, p. 854](https://www.jstage.jst.go.jp/article/tetsutohagane1955/71/7/71_7_853/_pdf)
- [Ban-ya 1993](https://www.jstage.jst.go.jp/article/isijinternational1989/33/1/33_1_2/_article)
- [Xiao and Holappa 1993](https://www.jstage.jst.go.jp/article/isijinternational1989/33/1/33_1_66/_article) — chromium conversion factors and interaction data.
- [Xiao and Holappa 1995, p. 322, Eqs. 9–10](https://www.pyrometallurgy.co.za/InfaconVII/319-Xiao.pdf) — oxide reaction coefficients.
- [Xiao, Kou, and Fang 2018, p. 419, Table 1](https://jproeng.ipe.ac.cn/CN/abstract/abstract3039.shtml) — chromium dissolution coefficient.
- [Xiao, Holappa, and Reuter 2002](https://doi.org/10.1007/s11663-002-0039-9)
- [Sigworth and Elliott 1974](https://doi.org/10.1179/msc.1974.8.1.298) — Wagner first-order interaction parameters.

## Critical findings

### C1. Chromium equilibrium uses the wrong energy units and the wrong metal standard state — RESOLVED

Evidence:

- `slag_model/data.py:170` stores `(25690, -13.36)` under a comment that calls
  the reaction `2 CrO1.5 + [Cr] = 3 CrO`.
- `slag_model/equilibrium.py:34` treats those coefficients as J/mol.
- `slag_model/solver.py:102-107` supplies $a_{Cr}=f_{Cr}[\%Cr]$.

The published reaction is

$$2\,CrO_{1.5}+Cr(s)=3\,CrO,$$

with

$$\Delta G^\circ=(25690-13.36T)\ \mathrm{cal\ mol^{-1}},$$

and pure solid chromium as the metal standard state (Xiao and Holappa 1995,
p. 322, Eqs. 9–10). The dissolution coefficient for $Cr(s)=[Cr]$ on the Henrian
1 mass-% standard state comes from Xiao, Kou, and Fang 2018 (p. 419, Table 1).
At the audit baseline the
code made two independent changes without conversion: cal was interpreted as
joules and $Cr(s)$ was replaced by dissolved chromium on the 1 mass-% standard
state.

For the library's metal activity convention, the required reaction energy is

$$\Delta G^\circ_{1\%Cr}
=4.184(25690-13.36T)-(19246-46.86T)\ \mathrm{J\ mol^{-1}},$$

where the second term is $Cr(s)=[Cr]_{1\%}$. At 1823.15 K this gives
$K_{Cr}=0.00878747$, versus 0.915831 at the audit baseline. The reference-case chromium
ratio consequently changes from 1.538 to 0.252 when this and C2 are corrected.

**Resolution and evidence:** `k_Cr` now combines the cal/mol oxide reaction with
the J/mol dissolution reaction on the Henrian 1 mass-% Cr basis. Verified by
`tests/test_equilibrium.py::test_k_Cr_source_reaction`,
`examples/thermodynamic_review/test_paper_validation.py::test_xiao_chromium_equilibrium_is_on_library_metal_standard_state`,
and the full-reference `tests/test_solver.py::test_reference_case_after_p0_corrections`.

### C2. Ban-ya Fe redox equation is called with mixed standard states — RESOLVED

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

**Resolution and evidence:** Ban-ya Eq. 24 now receives both Fe coefficients
from the RS activity map in the coupled and fixed-$P_{O_2}$ solvers; conventional
$a_{FeO}$ remains in the oxygen-potential closure. Verified by
`examples/thermodynamic_review/test_paper_validation.py::test_banya_eq_24_is_coupled_with_regular_solution_gammas`,
`examples/thermodynamic_review/test_paper_validation.py::test_fixed_po2_fe_redox_uses_two_rs_gammas`,
and the independent Fe quotient identity in
`tests/test_solver.py::test_reference_case_after_p0_corrections`.

### C3. P2O5 and Al2O3 activities are calculated as one-cation pseudo-components — RESOLVED

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

At the audit baseline, the one-power expression returned as `a_P2O5` was neither
this conventional formula-unit activity nor clearly labeled pseudo-component
activity. The same stoichiometric problem applied to Al2O3. In addition, the
Al2O3 conversion was not known in the project, so a conventional Al2O3 activity
could not be reported.

**Resolution and evidence:** RS cation activities and conventional species
activities are now separate. P2O5 uses the two-unit formula conversion; Al2O3
has no default conventional activity and only accepts a labeled custom
conversion. Verified by `tests/test_activities.py::test_banya_p2o5_formula_unit_activity_relation`,
`tests/test_activities.py::test_zero_phosphorus_and_unconverted_alumina`,
`tests/test_activities.py::test_custom_alumina_conversion_uses_two_rs_units`,
and the v2 report tests.

### C4. Non-finite input can falsely converge and serialize invalid JSON — RESOLVED

Evidence:

- `slag_model/io.py:151-155` checks type and sign but not `math.isfinite`.
- `slag_model/io.py:169-172` accepts strings convertible to `NaN` or infinity.
- `slag_model/solver.py:171-174` can declare convergence because Python's `max`
  can ignore a trailing `NaN` depending on ordering.
- `slag_model/cli.py:104` uses the default JSON encoder, which emits non-standard
  `NaN` tokens.

At the audit baseline, setting `P_CO_atm` to `"NaN"` returned a nominally
converged result after 22 iterations with `C_wtpc=NaN` and `Q/K_CO=NaN`.

**Resolution and evidence:** parsing rejects non-finite scalars with field
context; solver state, targets, residuals, and updates are checked before
convergence; the CLI serializes with `allow_nan=False` before opening the output
file. Verified by `tests/test_io.py::test_nonfinite_pressure_string_rejected`,
`tests/test_io.py::test_bare_json_nonfinite_constant_rejected`,
`tests/test_solver.py::test_solve_rejects_nonfinite_trailing_target`, and
`tests/test_cli.py::test_nonfinite_report_does_not_create_output_file`.

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

- `paper_cases.json`: source citations, exact constants and equations, five
  selected Ban-ya 1985 rows, eight Xiao 2002 rows, and model-output snapshots
  explicitly separated from source measurements.
- `test_paper_validation.py`: source-data identity checks, isolated activity
  comparisons, redox diagnostics, and snapshots of current calculated values.
- `test_defect_regressions.py`: executable reproductions for non-paper numerical
  and input/output defects.

Run from the repository root:

```powershell
uv run pytest examples\thermodynamic_review -p no:cacheprovider -v
```

Paper-review pack result (2026-09-24):

```text
70 passed, 7 xfailed
```

The seven strict `xfail` markers correspond to documented normal-priority findings: chromium root bracketing,
fixed-pressure diagnostic pressure, elemental-total warning, trace-species
cancellation, the Fe3+-Ca2+ coefficient, selectable Xiao-2002 parameter sets,
and the Mn-C interaction coefficient. The four P0 critical regressions pass.
The former factor-of-2.1 experimental envelope has been removed because the
paper does not prescribe a quantitative activity-error bound. The recorded
model-to-measurement differences remain contextualized against Fig. 11; they are
not on their own evidence of a code defect or of a model outside the paper's
reported accuracy.

## Recommended repair order

This was the original pre-P0 triage. Items 1-4 are complete under P0 Tasks 1-6;
the normal- and low-priority findings below remain open.

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

## Audit baseline verification (2026-09-23)

- `uv run pytest -p no:cacheprovider -v` -> 56 passed.
- `uv run ruff check --no-cache .` -> passed before audit artifacts were added.
- `uv run ruff format --check --no-cache .` -> 21 files already formatted before
  audit artifacts were added.
- New audit pack -> 21 passed, 14 expected failures at the audit baseline.
- The reference case was recalculated independently with source-consistent Fe and
  Cr equations.
- The supplied workbook was inspected as formulas and cached values, not only as
  rendered output.
- Relevant pages of both bundled HSC PDFs were rendered and visually checked.

No production code was changed during the original audit.

## P0 closeout verification (2026-09-24)

- `uv run pytest -p no:cacheprovider -v` -> 86 passed.
- `uv run pytest examples/thermodynamic_review -ra -p no:cacheprovider` ->
  29 passed, 7 xfailed; all remaining xfails correspond to the normal-priority
  findings listed above.
- `uv run ruff check --no-cache .` -> passed.
- `uv run ruff format --check --no-cache .` -> 27 files already formatted after
  mechanical formatting of the six earlier P0 branch files and `tests/test_solver.py`.
- `git diff --check` -> passed.
- `tests/test_solver.py::test_reference_case_after_p0_corrections` -> passed;
  `K_Cr=0.00878747397894681`, `r_Fe=0.13189054746`,
  `r_Cr=0.25229986894`, `[C]=0.0319280093 wt%`, and
  `P_O2=3.45717e-10 atm` under the unchanged hybrid parameter matrix.

These values are a regression sentinel for the unchanged hybrid matrix, not
empirical validation. A future named parameter-set change must version the
sentinel rather than silently widen or move its tolerances.
