# Quark splittings: the gg → ZZ test

The quark-splitting result in the [Letter](../manuscript_prl/letter.tex),
with the derivation and checks in the [Supplement](../manuscript_prl/supplement.tex)
(`sec:quarks` and `sec:zz`), predicts the two-loop DPS logarithms for a gluon that
splits into a quark–antiquark pair, with ΔC = 2C_F − C_A = −1/N_c, B = 6C_F and the scale
Q*² = s^{N_c²}/m^{2(N_c²−1)}. This folder tests the prediction in the two-loop gg → ZZ
amplitudes through massless quark loops of von Manteuffel and Tancredi (arXiv:1503.08835). The
results are collected in the Supplement's gg → ZZ comparison and its numerical-fit discussion.

The prediction holds in both colour structures, at β = 1/4, 1/2 and 2/3, for both kinds of
amplitude that carry a DPS term at one loop:

- the amplitudes with a one-loop logarithm, ++++ and ++−− (and their parity images): the
  coefficients of G³ and G²;
- the amplitudes in which J_z changes by two, +−++ and +−−− (and their images), whose one-loop
  forward limit is a constant K: the coefficients of G² and G, and the absence of G³.

The amplitudes with ΔJ_z = 0 and no one-loop logarithm (++00, +−+−) have no logarithm at two
loops either. The tables at the end are written by `make_summary.py` from the fits and the logs.

The one-loop coefficients have closed forms (normalisation of the source, x = (1+β)/2):

    c₁(++++) = −8πi x²,   c₁(++−−) = −8πi (1−x)²,   K(+−++) = K(+−−−) = −4πi x(1−x),

and K(++00) = −8πi x(1−x) for the amplitude with longitudinal bosons, which has no logarithm.
In the normalisation of the paper (a_s per loop) the moduli are twice these. The x², (1−x)²
dependence and the restriction of the logarithm to ++++ and ++−− agree with Gaunt–Stirling
(arXiv:1103.1888: [s − 2M² ± sβ]/s = 2x², 2(1−x)²). An independent full fermion trace now
checks the absolute normalization, including the measure, colour and two orientations:
see the Supplement's "Helicity coefficients and normalisation" subsection.
The repository-only [derivation notes](../../derivations/README.md) also include a restored
standalone trace calculation in `derivations/checks/independent/phase3/zz_dirac_trace.py`;
that four-channel symbolic check has a narrower scope than the original amplitude comparison.
In the retained amplitude comparison, all 36 helicities were checked and the twelve
nonzero comparisons agree with the fits within their window shifts. The section 2 and
section 3 ZZ displays of that reference differ by a factor of two; the trace fixes the
single-box magnitude without assuming that both displays have the same normalization.
The two-loop prediction for the angular constant assumes that K comes from the DPS region. That
the relation holds in the data is evidence for it, not a separate calculation; the closed form
of K alone would not show it, since K(++00) has a closed form too and no logarithms.

The previously reported C_A discrepancy was caused by two interchanged master-family
identifications in our FORM reader. The finite-export label `MN13m` denotes `N13x34x12`,
and `MN13cm` denotes `N13x34`. The mass exchange at fixed t combines the incoming- and
outgoing-leg crossings. Six comparisons with the bare weight-four coefficients independently
establish these identifications. Correcting the mapping restores agreement with the C++
amplitudes and removes the spurious G/p_T² growth. The retained six-coefficient
result is `scripts/fixtures/nonplanar_weight4_residues.json`; the native reader,
its tests, source provenance and commands for regenerating the primary check are
documented in [the mapping-check instructions](../../checks/README_ZZ_CA_MAPPING.md).

The stored-point analysis is complete. The compiled algebraic and FORM comparisons
used by this release are recorded in `logs/` and `points_routeB/`; the corresponding
scripts are retained under `scripts/`. Rechecking these stored records is distinct
from a fresh evaluation of the external amplitude expressions.

## Inputs

From the VVamp project (vvamp.hepforge.org), gluon channel and master integrals; checksums in
`inputs/SHA256SUMS`. The gluon-channel links on the project page returned "not found" on
13 October 2026. Historical notes say that the files used here were "supplied by
the author", but do not unambiguously identify the supplier or delivery route.
The checksums identify the retained versions; they do not establish how those
copies were acquired. Confirm that route before finalising the Data Availability
statement. On 12 November 2026 the essential inputs were recovered publicly with
exactly the retained checksums. The C++ archive is available from
[CERN's source mirror](https://lcgpackages.web.cern.ch/lcgpackages/tarFiles/sources/ggvvamp-1.0.tar.gz).
The three other required files are available at
`https://vvamp.hepforge.org/downloads/FILENAME` (substitute the filenames below).
These recovered downloads establish a reproducible route today; they do not
identify the historical supplier. Third-party files remain downloaded inputs.

| file | used for |
|---|---|
| `ggvvamp-1.0.tar.gz` | route A: the amplitude coefficients (generated C++ code) |
| `gg_projectors-1.0.inc` | identification of the 20 tensor structures |
| `gg_Aj_A_qT-1.0.inc.gz`, `vvpmasters_opt-1.0.m.gz` | route B (cross-check) |
| `vvpmasters_opt-1.0.inc.gz` | optional direct comparison of the FORM and Wolfram master exports; recovered and checksum-pinned on 12 November 2026 |
| `gg_Aj_A-1.0.inc.gz` | independent bare-residue determination of the crossed nonplanar labels |
| `vvpmasters_def-1.0.m`, `vvpmasters_relations-1.0.m`, `*.yaml` | definitions of the underlying integrands and crossings; optional for numerical production |

## Method

**Coefficients (route A).** ggvvamp computes 9 of the 20 coefficients A_j^[A] of the tensor
decomposition, at one and two loops in the qT scheme at μ² = s; the other 11 follow from
crossing relations (ggvvampprec.cpp). Its generated code has three parts: the polylogarithms
(p0002.cpp, through GiNaC), 492 004 lines of arithmetic (p0003–p2463.cpp) and the assignment of
the coefficients (p2464.cpp), which combines the two-loop C_F and C_A parts with C_F = 4/3 and
C_A = 3. `translate_ggvvamp.py` turns them into Python, keeping the C_F and C_A parts separate;
`ggvv.py` runs them with mpmath at any precision, with Li₂, Li₃, Li₄ from mpmath and Li₂,₂ from
its defining series or, near the edge of convergence, the integral representation. Li₂,₂
follows eqs. (5.14), (5.15) of arXiv:1503.04812, Li₂,₂(a,b) = Σ_{i>j} aⁱbʲ/(i²j²) with the
arguments chosen such that |a| < 1, |ab| < 1. At the benchmark point all 401 Li₂,₂ calls of the
code satisfy these conditions in the order given, and only 159 would with the arguments
exchanged (`check_li22_order.py`); at every scan point they hold too (summary.json).

**Tensor structures.** With μ, ν, ρ, σ the indices of gluon 1, gluon 2, Z₃, Z₄:

    T₁ = g^{μν}g^{ρσ},  T₂ = g^{μρ}g^{νσ},  T₃ = g^{μσ}g^{νρ},
    T₄…T₇ = g^{μν} p_a^ρ p_b^σ  with (a,b) = (1,1), (1,2), (2,1), (2,2),
    T₈, T₉ = g^{μρ} p₃^ν p_{1,2}^σ,   T₁₀, T₁₁ = g^{μσ} p₃^ν p_{1,2}^ρ,
    T₁₂, T₁₃ = g^{νρ} p₃^μ p_{1,2}^σ,  T₁₄, T₁₅ = g^{νσ} p₃^μ p_{1,2}^ρ,
    T₁₆ = g^{ρσ} p₃^μ p₃^ν,   T₁₇…T₂₀ = p₃^μ p₃^ν p_a^ρ p_b^σ  with (a,b) as for T₄…T₇.

This basis was not taken from the paper: `check_tensors.py` builds it with explicit vectors in
six dimensions, computes M_ik = T_i†·T_k with the d-dimensional polarisation sums of the paper,
and finds that M⁻¹ equals the projector matrix of `gg_projectors-1.0.inc` at d = 6 to 4 × 10⁻⁴⁰.
It is also consistent with the crossing relations of the code.

**Helicity amplitudes.** On-shell Z bosons of equal mass m = 1, in the centre-of-mass frame with
gluon 1 along +z and Z₃ at polar angle θ in the xz plane; helicity polarisation vectors, which
satisfy the gauge conditions ε₁·p₂ = ε₂·p₁ = 0 of the decomposition; M = Σ_j A_j T_j(ε₁, ε₂, ε₃*,
ε₄*). ΔJ_z = (λ₁ − λ₂) − (λ₃ − λ₄) in the forward direction. Overall couplings and the colour
factor δ^{ab} are stripped, as in the source.

**Forward lines.** As for HH: x̂ = (1−β)/(1+β), s = (1+x̂)²/x̂, −t = x̂ + δ, δ = 10⁻ᵏ,
k = 4, …, 12, G = ½ ln(1 − tu) on the principal branch, p_T² = (tu − 1)x̂/(1+x̂)². The
coefficients carry up to four inverse powers of the Gram determinant, and many digits cancel
towards the forward limit; the points are computed with 30 + 5k digits. One point repeated with
15 more digits (see *Precision*).

**Fits.** At each β, Σ_{j≤J} f_j G^j + p_T² Σ_{l≤P} h_l G^l through J + P + 2 consecutive
points ending at the smallest δ: J = 1, P = 2 at one loop (5 points) and J = 3, P = 3 at two
loops (8 points). The uncertainty quoted is the change when the window moves one point up
(δ = 10⁻⁴ … 10⁻¹¹ instead of 10⁻⁵ … 10⁻¹²). With P = 2 at two loops the deviations of f₂/c₁ of
++++ and ++−− grow to the value given under the first two-loop table (10⁻⁷ in the units of the
source), the effect of the truncated power corrections. `analyse_windows.py` refits the leading
two-loop coefficients of ++++ and +−++ through every window of consecutive points, without p_T²
terms and with P = 0, 1, 2, 3 (`results_windows.json`, Figure 5 of the paper): every model
improves as the window moves to smaller δ, ++++ needs P = 3, and for +−++ P = 2 does somewhat
better than P = 3, which is used for both. `check_fit_models.py` (`results_fit_models.json`) makes
the same choice without the predictions: for all eight amplitudes of the two groups, both colour
factors and both powers, the window shift is smallest for P = 3 in the first group and P = 2 in the
second, the models with the smallest deviations; with P = 3 the deviations are 0.010–0.018 of the
window shift. Least-squares fits of the P = 3 model through one more point lie halfway between the
two windows, about half a window shift from the predictions.

**Approach to the forward limit.** `forward_curve.py` stores (M − K)/K point by point for the
one-loop +−++ amplitude (`results_forward.json`, Sec. 6.2 of the paper). K is taken in closed
form: the fitted K of `results_zz.json` agrees with it to about 10⁻¹⁷ (recorded in the file), but
at δ = 10⁻¹² the difference M − K is only about 10⁻¹² of K, so that a fitted K would limit
(M − K)/(K p_T²) to about five digits. The difference vanishes like p_T² with no logarithm:
Re (M − K)/(K p_T²) is constant for δ ≤ 10⁻⁶ to the spread recorded in the file, as for the
gg → HH amplitude of Sec. 2.5.

**From 90 degrees.** `zz_wide.py` evaluates the amplitudes on the same lines of fixed β away
from the forward limit, at sin²θ = 10^(−J/2) for J = 0 (90°), 1, …, 7, that is from
p_T² = β²/(1−β²) m² down to about 3 × 10⁻⁴ of that, with the working precision of the scan
for the corresponding δ (`points_wide/`, 24 points). `wide_curve.py` joins them to the scan
and stores, for the one-loop ++++ and +−++ amplitudes, M/c₁ − G and M/K from 90° to
δ = 10⁻¹² (`results_wide.json`, Figure 2 of the paper, with the gg → AA amplitudes of P1).

**Prediction in the normalisation of the source.** The source has a factor α_s/(2π) per loop,
The manuscript uses a_s = α_s/(4π), so M^(n) = 2ⁿ × (source) and the kernel prediction with G → G_q,
G_q = G − [(C_A − C_F)/ΔC] L, L = ln[x(1−x)] (L_x in the paper), becomes

    one-loop c₁G:   M^(2) = (c₁/2)[−(4/3)ΔC G_q³ − B G_q²] + O(G)
        C_F part:   f₃/c₁ = −4/3,   f₂/c₁ = −2L − 3
        C_A part:   f₃/c₁ = +2/3,   f₂/c₁ = 2L
    one-loop K:     M^(2) = (K/2)[−4ΔC G_q² + (4ΔC − 2B) G_q] + O(1)
        C_F part:   f₂/K = −4,      f₁/K = −4L − 2
        C_A part:   f₂/K = +2,      f₁/K = 4L − 2

In the paper's normalisation every ratio is twice this (eq. 6.3 and Table 7, written by the
paper's `make_tables.py` from `results_zz.json` and `summary.json`).

The finite remainders are in the qT scheme at μ² = s. For gg → VV a change of the infrared
scheme or of μ adds a constant times the one-loop amplitude (ggvvamp converts to Catani's scheme
by A⁽²⁾ → A⁽²⁾ + (−½π²C_A + iπβ₀)A⁽¹⁾), so none of the coefficients tested here depends on
either. n_f does not appear in the remainders at μ² = s.

## Cross-checks

All numbers below are in `summary.md` and `summary.json`, read from the logs of `make checks`,
`make cxx`, `make li22` and from the route-B points.

- **Compiled C++ code.** `check_cxx.py` compiles the generated files p0003–p2464.cpp unchanged
  (apart from taking C_F and C_A from the harness in p2464.cpp, so that the two parts come out
  separately) in the quad precision of the original code, with a small replacement header
  (`scripts/cxx/`) instead of GiNaC and CLN, and runs them on the polylogarithm values of the
  Python translation. The nine coefficients agree with the translation to 2 × 10⁻³⁰ (one loop)
  and 2 × 10⁻²⁸ (C_F and C_A) at the benchmark point, and to 3 × 10⁻¹⁹, 10⁻¹⁶ and 4 × 10⁻¹⁶ at
  β = 1/2, δ = 10⁻⁴, where quad precision loses digits to the cancellations. This checks the
  translation of the arithmetic and of the assignment, C_A part included, but not of p0002.
- **Li₂,₂.** `check_li22_values.py` compares the Li₂,₂ of the translation, at 30 digits, with an
  independent integral representation (summing over the smaller index first) at 60 digits, for
  every distinct argument pair of the code at the benchmark point and at δ = 10⁻¹² for the three
  velocities (1490 pairs): agreement to 3 × 10⁻²⁵, and to 4 × 10⁻³⁰ at the benchmark point and
  at β = 1/4, 1/2. The largest differences, at β = 2/3, are those of the translation's integral
  representation at 30 digits for pairs with a ≈ 0.96 and b ≈ 10⁻¹³–10⁻¹¹; with 45 digits it
  agrees with the reference to all digits shown (the scan uses 30 + 5k digits).
  `check_li22_scan.py` (`make li22-scan`, 3 November 2026; about 15 CPU hours) compares, at the working
  precision of every point of the scan, of its control and of the blind test at β = 4/5, every real
  argument on the integral branch (7163 in all; none complex), in both channels, with the same integral
  split at points accumulating at t = 1 and evaluated with 15 more digits, as gg → HH now evaluates it
  (`../hh/hh_two_loop/evaluate_hh.py`). The integral branch here is split at fixed points only, which
  for gg → HH, at 3k fewer digits, lost up to 21 digits when ab approaches 1 from below with b > 1. At
  30 + 5k digits the largest difference is 169 units in the last digit (β = 1/4, δ = 10⁻⁵, 55 digits,
  a ≈ 1 − 6 × 10⁻⁶, b ≈ 0.64, where the pole of 1/(1 − at) lies just beyond t = 1), 50 or less at every
  other point and about one unit at most of them: the loss is at most three of the 30 + 5k digits,
  against the k + 30 or so that the cancellations leave, so the stored points stand (App. A of Paper 2).
- **p0002.** `check_p0002.py` transcribes the 1613 statements of p0002.cpp (401 Li₂,₂, 501 Liₙ,
  19 logarithms, 3 ζ values and 689 algebraic statements) independently of
  `translate_ggvvamp.py` and compares them with the generated `p0002.py`: none differ. The
  C++ code takes the real part of every function value (`torealnumeric` in ggvvamp.cpp); the
  translation does the same. This matters because the compiled check above takes p0002 from
  the translation, and route B below checks only the function values that enter the one-loop
  and C_F parts.
- **Route B.** The same coefficients from the FORM file `gg_Aj_A_qT-1.0.inc` with the master
  integrals of `vvpmasters_opt-1.0.m` (`crosscheck_form.py`, `zz_routeB.py`; an independent
  parser, the masters evaluated from their own file, but the same Li₂,₂ routine). At the
  benchmark point of the code (s = 15625, t = −2384.89727896833, p₃² = 6467.215561, p₄² = 625),
  30 digits: all 20 one-loop coefficients agree with route A to 4 × 10⁻²⁷; the two-loop C_F
  parts of A₁, A₃, A₁₀, A₁₇ there, and of A₂ at the point with t and u exchanged, to 3 × 10⁻²⁶.
  At one forward point (β = 1/2, δ = 10⁻⁴, 40 digits) the C_F parts of all 36 helicity
  amplitudes agree to 4 × 10⁻²¹. This checks the polylogarithms (p0002) as well.
- **The C_A comparison after correcting the family mapping.** The C_A terms of A₁, A₃,
  A₁₀ and A₁₇ at the benchmark and A₂ at the crossed benchmark agree with route A to
  2 × 10⁻⁵⁶. All 36 forward helicity amplitudes agree to 3.2 × 10⁻²¹ at β = 1/2,
  δ = 10⁻⁴, using 40 digits, and 8.3 × 10⁻²³ at δ = 10⁻⁵, using 45 digits.
  The corrected values of p_T²|M(++++)| are 0.596, 0.109 and 0.0181 at δ = 10⁻³,
  10⁻⁴ and 10⁻⁵. The earlier inverse-power growth came from our mapping error.
  The earlier records are historical development data outside this release; the
  corrected comparisons are retained in `logs/` and `points_routeB/`.
- One-loop coefficients equal their closed forms to 5 × 10⁻¹⁸ or better at all three
  velocities.

## Rerunning

From this directory (`paper2/zz/` in the P1 package):

    make test                    # stored points and logs -> results_zz.json, summary.md,
                                 # summary.json, results_windows.json, results_forward.json,
                                 # results_wide.json, exact match (seconds)
    make analyse                 # rewrite them

To recompute amplitudes or rerun the checks, first copy the VVamp files into `inputs/` (at least
`ggvvamp-1.0.tar.gz`, `gg_projectors-1.0.inc`, `gg_Aj_A_qT-1.0.inc.gz`,
`vvpmasters_opt-1.0.m.gz`), then

    make inputs                  # checksums, unpack into build/
    make build                   # translate ggvvamp into Python (about 30 s)
    make checks                  # tensor basis, Li22 order, p0002, FORM file at the benchmark point (2 min)
    make cxx                     # compiled C++ code against the translation (about 15 min; g++)
    make li22                    # Li22 against an independent representation (about 45 min)
    make routeb                  # FORM file on the forward line (about 10 min)
    make point BETA=1/2 K=6      # one forward point, 30 + 5K digits (2-10 min)
    make scan                    # all 28 points (about 2.5 CPU hours; existing points are kept)
    make wide                    # the 24 points from 90 degrees (about an hour; existing points are kept)

and `make analyse` to update the summary. The 27 points of the scan and the precision control
are in `points/`, those from 90 degrees in `points_wide/`. Python ≥ 3.10 and mpmath only (and
g++ with libquadmath for `make cxx`); `build/` (about 330 MB, and about 650 MB more after
`make cxx`) is not included and is regenerated by `make inputs build`.

The tabulated results were checked twice by a re-fit of the raw points with separate
code, reproducing the table entries and the quark-kernel predictions. The current
manuscript presentation is the Supplement's gg → ZZ subsection (`sec:zz`).

## Current numerical results

The generated [forward-test results](summary.md) contain the current fit tables,
C++ and FORM comparisons, and corrected crossed-master mapping results.
Regenerate that single report with `make analyse` after updating its inputs.
Keeping the numerical report in one file avoids stale copies of superseded
comparison results in this README.
