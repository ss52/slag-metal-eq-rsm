# P0 Thermodynamic Corrections Design

Date: 2026-09-23  
Status: draft for user review; no production changes authorized by this document  
Evidence: `examples/thermodynamic_review/REVIEW.md` and its paper-case fixtures

## Intent and scope

Make the coupled slag/hot-metal result scientifically interpretable by fixing the four critical findings C1-C4 in the audit: chromium reaction units and metal standard state, the Fe redox regular-solution basis, formula-unit oxide activity reporting, and non-finite input/result handling. The user chose a versioned v2 JSON report rather than legacy aliases. Success means the source-equation tests pass, the ordinary suite is updated to assert corrected behavior, and no input can produce a successful report containing `NaN` or infinity.

This is a focused repair of the existing Ban-ya/Xiao/Wagner calculation, not a new thermodynamic framework. The normal-priority findings in `REVIEW.md` remain outside this plan. In particular, do not silently alter the interaction matrix, Wagner coefficients, Ti approximation, chromium root solver, or fixed-oxygen-pressure diagnostic semantics. Those need separate reviewed changes. The workbook is a labeled legacy regression aid, not a scientific oracle.

## Selected approach

Three approaches were considered:

1. Patch four expressions in place. Lowest code churn, but the shared `gamma` and `a_slag` dictionaries would continue to mix RS pseudo-components and conventional oxide formula units, making another basis error likely.
2. **Selected:** introduce an explicit activity result with separate RS and conventional maps, make the four scientific corrections, and expose a v2 report that carries species and standard states. This is a moderate, testable interface change confined to the affected modules.
3. Replace all thermodynamic quantities with a general units/reaction/standard-state engine. Stronger long-term abstraction but disproportionate for four P0 repairs and likely to entangle the normal-priority parameter-set work.

## Scientific contract

All energies passed to `exp(-DeltaG/(R*T))` are J/mol of the reaction as written. Temperature is K; pressure is atm relative to the papers' implicit 1-atm gas standard; slag `X` is cation fraction; dissolved metal `a_Cr=f_Cr[%Cr]` is Henrian on the 1 mass-% state.

### Chromium equilibrium

The source oxide reaction is

$$2\,CrO_{1.5}+Cr(s)=3\,CrO,$$

with $\Delta G^\circ=(25690-13.36T)$ cal/mol. The independent dissolution reaction is $Cr(s)=[Cr]_{1\%}$ with $\Delta G^\circ=(19246-46.86T)$ J/mol. The solver uses dissolved chromium, so the reaction energy it must use is

$$\Delta G^\circ_{\mathrm{Henrian}}=4.184(25690-13.36T)-(19246-46.86T)\quad\mathrm{J/mol}.$$

`k_Cr(T)` must use that combined expression. At 1823.15 K its source-equation target is 0.00878747397894681. Keep the oxide and dissolution coefficients separately named and unit-labeled; do not retain a generic `DELTA_G_CR` whose basis is ambiguous. Confirm the cited dissolution value and its reference state in the primary-source ledger before changing production data; if its provenance cannot be confirmed, pause the C1 implementation rather than fitting it to the workbook.

### Iron redox

Ban-ya Eq. 24 uses $\gamma_{FeO}^{RS}$ and $\gamma_{FeO_{1.5}}^{RS}$ together. Both come directly from the RSM `RTln_gamma_rs` values. The separate Fe-FeO oxygen-potential closure continues to use converted $a_{FeO(l)}$. Rename the redox function parameters to make the RS requirement visible at the call site and test both the coupled and fixed-$P_{O_2}$ callers.

### Slag activities

RSM components are per cation. The RS result is $a_i^{RS}=X_i\exp(RT\ln\gamma_i^{RS}/RT)$. For a one-cation oxide with a documented conversion, the conventional activity is $a_j=a_i^{RS}\exp(\Delta G_{conv,j}/RT)$. Keep these in separate maps; no conventional value exists merely because an RS value exists.

For Ban-ya's formula-unit reaction $P_2O_5(l)=2PO_{2.5}(RS)$,

$$a_{P_2O_5(l)}=(a_{PO_{2.5}}^{RS})^2\exp[(52720-230.706T)/(RT)].$$

`gamma_P2O5` is not reported: defining it would require an explicitly chosen mole-fraction basis for formula units, not the cation fraction $X_P$. `AlO1.5(RS)` is always reportable. `a_Al2O3` is unavailable by default because no documented conversion is present. A custom Al2O3 conversion, if retained, has the explicit reaction `Al2O3(reference)=2 AlO1.5(RS)` and its `A,B` are J/mol of Al2O3 formula unit; it requires an input `standard_state` label. Legacy `A,B`-only custom Al2O3 inputs must receive a migration error, not be silently reinterpreted. The optional FeO1.5 custom conversion likewise needs a standard-state label before appearing as a conventional result. No default `A=B=0` can masquerade as a known conventional conversion.

### Numeric validity

Input validation rejects booleans, malformed scalars, and all `NaN`/positive or negative infinity forms, including quoted strings and non-standard bare JSON constants. Every solver iteration checks the state, targets, relative changes, and updated iterates for finite values before the convergence test. A non-finite or mathematically undefined intermediate raises a named numerical failure, never a successful result. The CLI emits strict JSON (`allow_nan=False`) and does not leave a partial output file when serialization fails.

## Internal boundaries and v2 report

`rsm.py` owns a `SlagActivities` value with maps `gamma_rs_by_cation`, `a_rs_by_cation`, `gamma_conventional_by_species`, `a_conventional_by_species`, `delta_g_conversion_J_per_mol_species`, and `standard_state_by_species`. "Species" is deliberate: CrO1.5 is a one-cation model species, whereas P2O5 is a two-cation formula unit. The conventional maps contain only values for which a conversion is defined. `solver.py` consumes this value: RS Fe coefficients for Ban-ya Eq. 24, conventional FeO for oxygen potential, conventional CrO/CrO1.5 for Cr equilibrium. `io.py` only translates that state into a report; it must not perform thermodynamic conversion.

The report has `schema_version: 2`. Each `solution.components` entry identifies `cation`, `rs_species`, `X_cation`, `RTln_gamma_RS_J_per_mol_cation`, `gamma_RS`, and `a_RS`. A nested `conventional` object has `species`, `standard_state`, `rs_units_per_species`, `a`, and `DeltaG_conversion_J_per_mol_species` only when a conversion exists; otherwise it is `null` with a `conventional_unavailable_reason`. A conventional `gamma` is present only for one-cation species whose fraction basis is explicitly `X_cation`; it is absent for P2O5 and Al2O3. The P5+ row therefore contains both `a_PO2.5_RS` through `a_RS` and `a_P2O5(l)` through `conventional.a`; the Al3+ row has `a_AlO1.5_RS` and no default conventional activity. Top-level `equilibrium` and `metal` sections keep their physical meanings, with corrected values. The legacy Ti denominator-only row reports `model_status: denominator_only_legacy` and no invented activity.

The input schema remains composition-compatible. The two optional custom conversion objects acquire a required `standard_state` string when used. The report's `input.options` must be canonical and round-trippable for these options; old tuple-array serialization is not retained in v2. The CLI's human-readable table uses the same labels and prints unavailable conventional results as unavailable, not zero.

## Verification and acceptance

Use the strict paper-derived tests already in `examples/thermodynamic_review/` and remove the four relevant `xfail` markers as their behaviors are implemented. Where a test calls an obsolete public activity API, replace it with the new explicit interface while keeping the same independent paper equation. Add focused tests for both Fe-redox call paths, P2O5 and unavailable/default Al2O3, custom conversion semantics, and each non-finite input/iteration/report boundary. Check the independent 1823.15 K reference calculation after C1 and C2: $r_{Fe}\approx0.131891$, $r_{Cr}\approx0.252300$, $[C]\approx0.0319280$ wt%, and $P_{O_2}\approx3.45717\times10^{-10}$ atm under the current parameter matrix. These are integration sentinels for this parameter set, not empirical validation.

Run from the repository root: `uv run pytest -p no:cacheprovider -v`, `uv run ruff check --no-cache .`, and `uv run ruff format --check --no-cache .`. No P0 `xfail` may remain. Update outdated tests and documentation to the corrected source equations; do not loosen tolerances just to make workbook values pass.

## Sources

- [Ban-ya, ISIJ International 33 (1993), regular-solution model and Eq. 24](https://www.jstage.jst.go.jp/article/isijinternational1989/33/1/33_1_2/_article)
- [Xiao and Holappa, ISIJ International 33 (1993), chromium-bearing slag](https://www.jstage.jst.go.jp/article/isijinternational1989/33/1/33_1_66/_article)
- [Xiao and Holappa, INFACON VII (1995), chromium oxide reaction and cal/mol units](https://www.pyrometallurgy.co.za/InfaconVII/319-Xiao.pdf)
- [Paper-case provenance and executable equations](paper_cases.json)
