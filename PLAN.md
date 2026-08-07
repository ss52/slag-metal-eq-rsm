# PLAN.md - Slag Activity Coefficient Model with Slag-Metal Oxygen-Potential Balance

Spec version 1.0 - 2026-08-04
Context: EAF slag-metal equilibrium model (Ban-ya quadratic formalism / regular
solution model for the slag + Wagner interaction parameter formalism for the
steel). Later stages will couple this engine to HSC Chemistry GIBBS. This file
is the single source of truth for the implementation.

## 0. How to use this document

This is a complete, self-contained specification. Suggested opencode prompt:

  "Implement the Python package specified in PLAN.md exactly. Create the full
   module layout, the example input file, and the pytest suite. Run the tests
   and iterate until they all pass. Do not invent thermodynamic data: where
   the spec marks a value TBD, implement the config flag and default exactly
   as described."

## 1. Purpose

Compute activity coefficients of all slag components of an EAF slag at a given
temperature using Ban-ya's quadratic formalism (regular solution model, RSM),
with automatic redox splitting and slag-metal oxygen-potential balance:

1. The slag is entered in analytical form: iron as ONE number (total iron,
   given as FeO_total or Fe_total) and chromium as ONE number (total chromium,
   given as Cr2O3_total or Cr_total). The model splits iron into FeO/Fe2O3
   (Fe2+/Fe3+) and chromium into CrO/Cr2O3 (Cr2+/Cr3+) according to the
   equilibrium oxygen potential P_O2.
2. P_O2 is NOT an input. It is derived from the slag itself through the Fe-O
   equilibrium (the high-FeO slag buffers the oxygen potential).
3. The steel carbon content that is in equilibrium with this slag is computed
   from the C-O equilibrium at the same P_O2. After the solve, slag and steel
   share exactly one oxygen potential: the two phases are in balance by
   construction. This is the "adjust P_O2 via steel carbon" closure.
4. Metal solute activity coefficients (needed for the Cr redox reaction and
   the carbon balance) come from Wagner's interaction parameter formalism
   (WIPF, first order).

Primary output: activity coefficients gamma_i and activities a_i of all slag
components on the pure-stable-oxide (Raoultian) standard state, the redox
splits, P_O2, and the balanced steel carbon content.

## 2. Physical basis (for the implementer)

- The slag is modeled as a regular solution of CATIONS (Ban-ya 1993, ISIJ
  Int. 33(1):2-11). Oxide component activities are obtained from cation
  fractions and interaction energies alpha_ij.
- Raw RSM activity coefficients refer to a hypothetical regular-solution
  standard state. Conversion factors (Delta G = A + B*T) shift them to the
  pure stable oxide at T (solid SiO2, solid Cr2O3, liquid FeO, ...), i.e. the
  Raoultian standard state used by HSC Chemistry.
- Fe redox: FeO(slag) + 1/4 O2 = FeO1.5(slag). The Fe3+/Fe2+ ratio is slaved
  to P_O2 (Ban-ya Eq. 24). One may fix either the ratio or P_O2, never both.
- Cr redox: 2 CrO1.5(slag) + [Cr](metal) = 3 CrO(slag) (Xiao). Oxygen cancels;
  the ratio is pinned by the chromium activity in the steel. No P_O2 needed.
- ONE oxygen potential: at equilibrium, slag and steel share a single P_O2.
  The slag sets it through a_FeO (Fe-O equilibrium); the steel carbon content
  consistent with it follows from the C-O equilibrium. Numerically the solver
  reads P_O2 off the slag and then computes the balanced C - this is
  mathematically identical to "adjust C until the C-O potential matches the
  slag potential", but converges far better.

## 3. Tech stack and module layout

Python 3.11+, numpy, pytest. Type hints everywhere. No pandas/scipy needed.
Use uv to manage dependencies. Use the latest stable version of Python.

  slag-activity-model/
    PLAN.md                (this file)
    pyproject.toml
    README.md              (short usage: how to run, input format)
    slag_model/
      __init__.py
      constants.py         # R, molar masses, cations per formula unit
      data.py              # alpha matrix, conversion factors, e_ij, DeltaG coeffs
      rsm.py               # cation fractions, RSM gammas, standard-state conversion
      redox.py             # Fe split (Ban-ya Eq. 24), Cr split (exact cubic)
      metal.py             # WIPF, a_Cr, balanced carbon, dissolved oxygen
      equilibrium.py       # K_FeO, K_CO, K_Cr, P_O2 closure
      solver.py            # damped fixed-point driver
      io.py                # JSON input parsing/validation, report builder
      cli.py               # entry: python -m slag_model input.json [-o report.json]
    examples/
      eaf_slag_01.json     # reference slag (see section 5)
    tests/
      test_data.py         # matrix symmetry, constants
      test_rsm.py          # workbook regression (major5 mode), sum X = 1
      test_metal.py        # WIPF regression
      test_redox.py        # fixed-P_O2 redox regression
      test_solver.py       # balanced solve, residuals, mass balance

## 4. Units and constants

- R = 8.314 J/mol/K (exactly 8.314, to reproduce the reference workbook)
- T in K; energies in J/mol; pressures in atm; compositions in wt%
- amounts in mol per 100 g of as-received slag
- ln(10) * R = 19.144 J/mol/K (for log10 conversions; prefer computing
  log10(gamma) directly instead of the workbook's "19.1*T" denominator form)

## 5. Input

JSON file. Example (examples/eaf_slag_01.json) - the reference EAF slag:

  {
    "temperature_K": 1823.15,
    "P_CO_atm": 1.0,
    "slag_wtpc": {
      "SiO2": 27.495, "CaO": 35.597, "MgO": 1.999, "Al2O3": 7.973,
      "MnO": 4.762, "TiO2": 0.725, "P2O5": 0.0,
      "FeO_total": 18.851, "Cr2O3_total": 2.58
    },
    "metal_wtpc": { "Cr": 0.038, "Mn": 0.007, "P": 0.069 },
    "options": {
      "cross_terms": "full",
      "ti_handling": "exclude_renormalize",
      "sio2_conversion": "workbook",
      "iron_input": "FeO_total",
      "chromium_input": "Cr2O3_total",
      "a_Fe": "unity"
    }
  }

Rules:
- Exactly one of FeO_total / Fe_total must be present (same for Cr forms).
- Unknown oxide keys -> error. Negative values -> error. Missing metal keys
  default to 0.0.
- Slag wt% sum should be 100 +/- 1.0; otherwise warn (do not fail). Note:
  after the redox split the oxide mass sum can drift slightly from 100 g per
  100 g basis because Fe2O3/Cr2O3 carry more oxygen than FeO/CrO. This is
  expected, not an error. All cation-fraction math is ratio-based and
  unaffected.
- temperature_K default 1823.15; P_CO_atm default 1.0.

## 6. Thermodynamic data (data.py)

### 6.1 Oxides, molar masses, cations per formula unit

  oxide    M (g/mol)   cations   cation
  SiO2     60.084      1         Si4+
  CaO      56.077      1         Ca2+
  MgO      40.304      1         Mg2+
  Al2O3    101.961     2         Al3+
  MnO      70.937      1         Mn2+
  TiO2     79.866      1         Ti4+   (not an RSM cation; see ti_handling)
  P2O5     141.943     2         P5+
  FeO      71.846      1         Fe2+
  Fe2O3    159.688     2         Fe3+
  CrO      67.996      1         Cr2+
  Cr2O3    151.990     2         Cr3+
  (Fe 55.847, Cr 51.996 for total-element input forms)

### 6.2 Cation interaction energy matrix alpha_ij (J/mol), symmetric

Order: [Fe2+, Fe3+, Ca2+, Mg2+, Mn2+, Si4+, Al3+, P5+, Cr3+, Cr2+]

  Fe2+:      0  -18660  -31380   33470    7110  -41840  -41000  -31380      0      0
  Fe3+:  -18660      0  -96810   -2930  -56480   32640 -161080   14640      0      0
  Ca2+:  -31380  -96810       0 -100420  -92050 -133890 -154810 -251040  44235  -6140
  Mg2+:   33470   -2930 -100420       0   61920  -66940  -71130  -37660  28085   5520
  Mn2+:    7110  -56480  -92050   61920       0  -75310  -83680  -84940      0      0
  Si4+:  -41840   32640 -133890  -66940  -75310       0 -127610   83680 -48975 -60540
  Al3+:  -41000 -161080 -154810  -71130  -83680 -127610       0 -261500 -46210 -30285
  P5+:   -31380   14640 -251040  -37660  -84940   83680 -261500       0      0      0
  Cr3+:       0       0   44235   28085       0  -48975  -46210       0      0  32710
  Cr2+:       0       0   -6140    5520       0  -60540  -30285       0  32710      0

Notes:
- alpha(Fe2+,Al3+) = -41000. The reference workbook contains a typo (-410) in
  one triangle of its matrix; its cross-term formulas use -41000, which is the
  correct value. Implement the matrix symmetric with -41000.
- Zeros for Fe-Cr, Mn-Cr, P-Cr pairs are a DATA GAP (the Cr parameter set only
  covers Ca/Mg/Si/Al/Cr partners). Keep them as 0.0; do not invent values.

### 6.3 Standard-state conversion factors: Delta G_conv = A + B*T (J/mol)

  component    A        B          note
  FeO          -8540    +7.142     FeO(l) = FeO(R.S.) [Ban-ya Table 3]
  SiO2         +51346   -13.88     "workbook" option (default; source: user workbook)
  SiO2         +27030   -1.983     "banya" option: SiO2(beta-cr) = SiO2(R.S.)
  CaO          +18160   -23.309    CaO(s) = CaO(R.S.) [Ban-ya]
  MgO          +34350   -16.736    MgO(s) = MgO(R.S.) [Ban-ya]
  MnO          -32470   +26.143    MnO(s) = MnO(R.S.) [Ban-ya]
  P2O5         +52720   -230.706   P2O5(l) = 2 PO2.5(R.S.) [Ban-ya]
  CrO          +77150   -33.5      CrO(l) = CrO(R.S.) [Xiao]
  CrO1.5       +74967   -37.5      CrO1.5(s) = CrO1.5(R.S.) [Xiao]; solid at 1823 K
  Al2O3        TBD      TBD        DATA GAP - see section 12. Default A=0, B=0 + warn.
  FeO1.5       0        0          default none (workbook behavior); Fe2O3 is nearly
                                   solid at 1823 K (Tm approx 1838 K) - refine later.

Final coefficient:  R*T*ln(gamma_i) = [R*T*ln(gamma_i)]_RS + A_i + B_i*T
Activity:           a_i = gamma_i * X_i   (X_i = cation fraction)

### 6.4 Equilibrium reactions and constants

(a) Fe(l) + 1/2 O2 = FeO(l):   Delta G = -232600 + 47.9*T  J/mol
    K_FeO = a_FeO / (a_Fe * sqrt(P_O2));  K(1823.15 K) approx 1.45e4
    a_Fe = 1.0 (option "unity") or X_Fe of the metal (option "xfe").

(b) [C] + 1/2 O2 = CO(g):      Delta G = -134300 - 45.40*T  J/mol
    K_CO = P_CO / (a_C * sqrt(P_O2));     K(1823.15 K) approx 1.66e6
    STANDARD-STATE WARNING: carbon is on the Henrian 1 wt% scale. The
    coefficients above = DeltaG(C(gr) + 1/2 O2 = CO) = -111710 - 87.66*T
    MINUS DeltaG(C(gr) = [C]_1wt%) = +22590 - 42.26*T. Never combine the
    graphite-based K with a Henrian a_C (error of ~3 orders of magnitude).
    a_C = f_C * [%C].

(c) 2 CrO1.5 + [Cr] = 3 CrO:   Delta G = +25690 - 13.36*T   J/mol  [Xiao]
    K_Cr = a_CrO^3 / (a_CrO1.5^2 * a_Cr); K(1823.15 K) approx 0.916
    a_Cr = f_Cr * [%Cr] (Henrian 1 wt%).

(d) Fe redox (Ban-ya Eq. 24):
    log10(Fe3+/Fe2+) = 6625/T - 2.77 + 0.25*log10(P_O2)
                       + log10(gamma_FeO) - log10(gamma_FeO1.5)
    gamma_FeO = converted coefficient; gamma_FeO1.5 = R.S. value unless an
    fe2o3_conversion is configured.

(e) Oxygen saturation of liquid Fe (for the [%O] diagnostic output):
    log10[%O]_sat = -6320/T + 2.734 ;  [%O] = a_FeO * [%O]_sat

### 6.5 Metal WIPF first-order interaction parameters e_i^j (1600 degC values,
    used at 1550 degC - accepted approximation)

  log10 f_i = sum_j e_i^j * [%j],  i,j in {C, Cr, Mn, P}

  i \ j     C        Cr       Mn       P
  C         0.14    -0.024   -0.012    0.051
  Cr       -0.12     0        0       -0.053
  Mn       -0.012    0        0       -0.0035
  P         0.13    -0.03     0        0

## 7. Equations (rsm.py, redox.py, metal.py, equilibrium.py)

### 7.1 Input splitting (per 100 g as-received slag)
  n_Fe_tot = FeO_total / 71.846            (iron_input = "FeO_total")
           = Fe_total  / 55.847            (iron_input = "Fe_total")
  n_Cr_tot = 2 * Cr2O3_total / 151.990     (chromium_input = "Cr2O3_total")
           = Cr_total / 51.996             (chromium_input = "Cr_total")
  Given r_Fe = n_Fe3+/n_Fe2+:  n_Fe2+ = n_Fe_tot/(1+r_Fe), n_Fe3+ = n_Fe_tot - n_Fe2+
  Given r_Cr = n_Cr2+/n_Cr3+:  n_Cr3+ = n_Cr_tot/(1+r_Cr), n_Cr2+ = n_Cr_tot - n_Cr3+
  Other oxides: n_cat,i = c_i * w_i / M_i

### 7.2 Cation fractions
  N = sum of cation mol over the RSM cation set (per ti_handling, section 11)
  X_i = n_i / N ;  assert sum(X) = 1 within 1e-12

### 7.3 RSM activity coefficients (Ban-ya quadratic formalism)
  R*T*ln(gamma_i^RS) = sum_j alpha_ij * X_j^2
      + sum_{j<k, j!=i, k!=i} (alpha_ij + alpha_ik - alpha_jk) * X_j * X_k
  - cross_terms = "full":   all pairs j<k (36 pairs per component for 10 cations)
  - cross_terms = "major5": ONLY pairs (Ca,Si), (Ca,Al), (Si,Al), (Mg,Si),
    (Ca,Mg) - exists ONLY to reproduce the reference workbook in tests.

### 7.4 Standard-state conversion: see 6.3.

### 7.5 Fe redox: r_Fe_new from Eq. 24 (explicit, see 6.4d).

### 7.6 Cr redox (exact relation - solve, do not use the sqrt fixed-point)
  From K_Cr = (gamma_CrO^3 * n_Cr2+^3) / (gamma_CrO1.5^2 * n_Cr3+^2 * a_Cr * N)
  with n_Cr2+ = r*n_Cr_tot/(1+r), n_Cr3+ = n_Cr_tot/(1+r):
      r^3/(1+r) = K_Cr * a_Cr * gamma_CrO1.5^2 * N / (gamma_CrO^3 * n_Cr_tot)
  LHS is monotonic in r > 0; solve by bisection or Newton (r in [1e-6, 1e6]).

### 7.7 Oxygen potential from the slag (Fe-O closure)
  P_O2 = ( a_FeO / (K_FeO * a_Fe) )^2 ,  a_FeO = gamma_FeO * X_Fe2+

### 7.8 Balanced steel carbon (C-O balance at the same P_O2)
  [%C] = P_CO / (K_CO * f_C * sqrt(P_O2))
  f_C = 10^( e_C^C*[%C] + e_C^Cr*[%Cr] + e_C^Mn*[%Mn] + e_C^P*[%P] )
  f_C depends weakly on [%C]: 2-3 fixed-point passes starting from f_C = 1.

### 7.9 Dissolved oxygen diagnostic: [%O] = a_FeO * 10^(-6320/T + 2.734)

### 7.10 WIPF for f_Cr, f_Mn, f_P, f_C: see 6.5.

## 8. Solver (solver.py)

Damped fixed-point iteration:

  state: r_Fe = 0.10, r_Cr = 1.0, C = 0.02 wt%
  repeat (max 200 iterations, damping w = 0.4):
      n       = split(totals, r_Fe, r_Cr)                    # 7.1
      X, N    = cation_fractions(n, ti_handling)             # 7.2
      g_RS    = rsm_gamma(X, alpha, cross_terms)             # 7.3
      g       = convert(g_RS, T, conversions)                # 7.4
      a_FeO   = g[FeO] * X[Fe2+]
      P_O2    = (a_FeO / (K_FeO(T) * a_Fe))^2                # 7.7
      r_Fe_n  = 10^(6625/T - 2.77 + 0.25*log10(P_O2)
                    + log10(g[FeO]) - log10(g[FeO1.5]))      # 7.5
      f       = wipf(C, Cr, Mn, P)                           # 7.10
      a_Cr    = f[Cr] * Cr_wtpc
      r_Cr_n  = solve_cubic(K_Cr(T)*a_Cr*g[CrO1.5]^2*N
                            / (g[CrO]^3 * n_Cr_tot))         # 7.6
      C_n     = P_CO / (K_CO(T) * f[C] * sqrt(P_O2))         # 7.8
      if max relative change of (r_Fe_n, r_Cr_n, C_n) < 1e-8: converged
      r_Fe += w*(r_Fe_n - r_Fe);  r_Cr += w*(r_Cr_n - r_Cr);  C += w*(C_n - C)
  else: raise NonConvergenceError with full state diagnostics

Typical convergence: 5-30 iterations. Gammas are weak functions of the splits,
so the loop is well behaved. Keep P_O2 in log10 space internally if preferred.

## 9. Output (JSON report + printed table)

  - echo of inputs and options
  - converged P_O2 (atm and log10), r_Fe, r_Cr, iteration count
  - split slag composition (wt%): FeO, Fe2O3, CrO, Cr2O3, and all other oxides
  - per slag component: X_i, RTln(gamma_RS), DeltaG_conv, gamma_final, a_i
  - metal: balanced [%C], [%O], f_C, f_Cr, f_Mn, f_P, a_C, a_Cr
  - equilibrium diagnostics: K_FeO, K_CO, K_Cr at T; reaction-quotient checks
    Q/K for Fe-O and C-O (must equal 1.000 within 1e-6 by construction)

## 10. Tests (tests/)

1. test_data: alpha matrix symmetric, zero diagonal; molar masses positive.
2. test_rsm (WORKBOOK REGRESSION, options cross_terms="major5",
   ti_handling="as_excel", sio2_conversion="workbook", T=1823.15,
   fixed P_O2 = 1e-9 atm, slag = eaf_slag_01.json):
   RT*ln(gamma_RS) totals per cation (J/mol, tol 5 J):
     Fe2+ +4031.71, Fe3+ -9027.69, Ca2+ -17758.64, Mg2+ -23733.98,
     Mn2+ -23294.70, Si4+ -26354.39, Al3+ -50118.08, Cr3+ +24359.11, Cr2+ +7682.77
   converted gammas: FeO 1.75 +/- 0.01, SiO2 0.98 +/- 0.01,
     CrO 4.79 +/- 0.02, CrO1.5 7.71 +/- 0.02
   r_Fe = 0.131 +/- 0.002 ; r_Cr = 1.25 +/- 0.05
   (note: the workbook's own sqrt fixed-point for r_Cr carries a small internal
   inconsistency of ~3%; the exact cubic of section 7.6 is the reference.)
3. test_metal (WIPF regression, C 0.018, Cr 0.038, Mn 0.007, P 0.069 wt%):
   log10 f: C +0.0051, Cr -0.0058, Mn -0.0005, P +0.0012 (tol 5e-4)
4. test_equilibrium: K_FeO in [1.3e4, 1.6e4], K_CO in [1.5e6, 1.8e6],
   K_Cr in [0.90, 0.93] at 1823.15 K. Verify the Henrian K_CO via the
   two-DeltaG decomposition of section 6.4b (standard-state trap test).
5. test_solver (balanced solve, cross_terms="full", defaults):
   converges in < 200 iterations; mass balance n_Fe2+ + n_Fe3+ = n_Fe_tot and
   n_Cr2+ + n_Cr3+ = n_Cr_tot to 1e-12; Q/K = 1 +/- 1e-6 for Fe-O and C-O;
   expected physical ranges: P_O2 in [2e-10, 6e-10] atm, r_Fe in [0.08, 0.11],
   balanced [%C] in [0.025, 0.045], [%O] in [0.04, 0.06].
6. sum(X) = 1 within 1e-12 in every solve.

## 11. Config flags (options block)

  flag              values                                  default
  cross_terms       "full" / "major5"                       "full"
  ti_handling       "exclude_renormalize" / "as_excel"      "exclude_renormalize"
                    ("as_excel": Ti4+ counted in N with alpha_Ti,j = 0,
                     reproduces the workbook; default drops TiO2 from the
                     cation sum and renormalizes - difference approx 0.5%)
  sio2_conversion   "workbook" / "banya"                    "workbook"
  fe2o3_conversion  "none" / custom {A, B}                  "none"
  al2o3_conversion  "none" / custom {A, B}                  "none" (TBD, sec. 12)
  a_Fe              "unity" / "xfe"                         "unity"
  iron_input        "FeO_total" / "Fe_total"                "FeO_total"
  chromium_input    "Cr2O3_total" / "Cr_total"              "Cr2O3_total"

## 12. Known data gaps - DO NOT INVENT VALUES

1. Al2O3 conversion factor (A, B): TBD. Primary source: Ban-ya 1993 Table 3
   (PDF in Zotero, item key JFSC4VDW). Until filled: default "none"
   (A=0, B=0) and emit a warning that gamma_Al2O3 is on the R.S. scale.
2. FeO1.5 conversion: default "none" (workbook behavior). Pure Fe2O3 is
   nearly solid at 1823 K (Tm approx 1838 K); a solid-standard-state
   conversion should be added before HSC coupling.
3. alpha for Fe-Cr, Mn-Cr, P-Cr pairs = 0 (limit of the Xiao parameter set).
4. SiO2 conversion: "workbook" (51346 - 13.88T, 26.0 kJ at 1823 K) vs
   "banya" (27030 - 1.983T, 23.4 kJ). The 2.6 kJ difference moves
   gamma_SiO2 by a factor ~1.19. Source of the workbook fit to be confirmed.
5. e_i^j are 1600 degC values used at 1550 degC; second-order WIPF terms
   neglected (acceptable at these dilute levels: C ~0.04%, Cr ~0.04%).

## 13. Acceptance criteria

- All tests in section 10 pass.
- CLI run: python -m slag_model examples/eaf_slag_01.json prints the full
  report and writes report.json.
- The balanced solve reproduces the workbook regression when the regression
  options are selected, and the physical ranges of test 5 with defaults.

## 14. Out of scope / future extensions

- Fixed-carbon closure (steel C given, slag total FeO floats) - the mirror
  problem; add as closure="metal_sets_PO2" later.
- HSC Chemistry GIBBS coupling (spreadsheet formulas, C# DLL implementing
  IActivityCoefficientModel2, or Python orchestration).
- UIPF / second-order metal formalism; temperature-dependent e_i^j.
- Ti4+ alpha row (then ti_handling="include").
