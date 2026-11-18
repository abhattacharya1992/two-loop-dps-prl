# Stage B: the one-loop DPS kernel and the boundary terms r₁, r₂

Goal: derive r₁ = r₂ = 0 in the boundary term R̃ ⊃ K̃₀(r₂X² + r₁X) of Paper 2, Sec. 5 (eq:kernel),
X = ln(Q_*²e^{−iπ}b²/b₀²), X_μ = ln(μ²b²/b₀²); r₀ is the by-product that gate C would test.

**Status: r₂ = r₁ = 0 shown at one loop, given the list of regions**, for
any splitting (g → gg, g → qq̄) and any colour-singlet 2 → 1 subprocesses; written into Paper 2 as
appendix E. A separate AI-assisted review of the derivation found no fatal error; its points are
addressed below (BPS decoupling as the argument for step 2, one homogeneous regulator, the Glauber
statements corrected and checked, the claims of the scalar graph narrowed, the b₀ convention). r₀
(gate C) and an end-to-end check of r₁, r₂ in a full model are open; the end-to-end test so far is
the data of Paper 2 (gg → HH, gg → ZZ) and Stage D (planar N = 4).

The boundary conclusion applies to the one-loop kernel correction (the two-loop amplitude).
The ultrasoft Wilson-line cancellation below does not by itself exclude higher-order
collinear or soft rapidity logarithms. The three-loop QCD coefficients in `../r0/theory/`
additionally assume no independent quadratic or higher-power boundary logarithm at that order.

## The result in three steps

1. **Regions.** At leading power in λ ~ p_T/m the extra gluon ℓ of the two-loop DPS term is hard (h),
   collinear to p₁ or p₂ (c₁: (1, λ², λ), c₂), soft (s: (λ, λ, λ)), Glauber (G: (λ², λ², λ)) or
   ultrasoft (us: (λ², λ², λ²)). h is the one-loop correction to H₁ and H₂ and is removed by the
   definition of the kernel. c, s and G resolve b and have virtuality μ_b² = b₀²/b². The us gluon does
   not resolve b; it sees the off-shell daughters through σ₁ = k_⊥²/(x(1−x)p₁⁺) and σ₂ (the light-cone
   energy denominators of the two splittings), so its virtuality is σ₁σ₂ ~ μ_b⁴/Q_*². It is leading
   power: in the explicit two-loop graph λ⁶F → const at ℓ_⊥ ~ λ² (B2.2).
2. **The ultrasoft gluons cancel** (B2.3). The argument is the decoupling of ultrasoft gluons from
   collinear fields (Bauer–Pirjol–Stewart field redefinition): at leading power they couple only
   through n·A_us; after the redefinition they sit in Wilson lines Y_n (Y_n̄) at the points x₁, x₂ where
   the daughters enter the vertices, and the numerators of the collinear subgraphs no longer contain
   them. Ultrasoft momenta do not resolve x₁ − x₂ (transverse b ~ 1/(λ√s), light-cone parts ~ 1/√s),
   so colour conservation, T_k + T_k̄ = T_{p₁}, combines the two Y_n into one in the representation
   of p₁; ⟨Y_n† Y_n̄⟩ at one point has no scale and is 1 to all orders, in any gauge. The one-loop
   illustration: the k⁻ integral of the DPS loop gives, relative to the lowest order and with eikonal
   couplings, attachment to k: 1/(ℓ⁻ − σ₁), to k̄: 1/(ℓ⁻ − σ₁), to the parent p₁ before the splitting:
   1/ℓ⁻ − 1/(ℓ⁻ − σ₁); the colour-weighted sum is T_{p₁}/ℓ⁻. With the full vertices in a covariant gauge
   the numerator also contains terms proportional to the daughters' virtualities (contact terms at
   lowest order); they contribute to single graphs (e.g. a k² term in graph (a) leaves
   −2πi/(x(1−x)(ℓ⁻ − σ₁)), outside the identity) and cancel in the sum, as the field redefinition
   requires — this part rests on the theorem, not on our scripts. The argument uses colour conservation
   at the splitting and the eikonal coupling, neither of which depends on the representation, so it
   holds for gluon and quark splittings alike.
3. **Bookkeeping** (B2.4). Regulators: dimensional regularisation and a rapidity regulator that leaves
   scaleless integrals zero (analytic or η regulator); in this scheme the telescoped us integral and the
   soft exchanges between partons at the same transverse position vanish. (The δ-regulator of B1 is
   used only for the coefficient of the rapidity logarithm, the Fourier transform of 1/ℓ_⊥², which is
   the same in both.) Write each region as (μ²/Λ_R²)^ε (c_R/ε² + d_R/ε + e_R), with Λ_R² = μ_b² for
   c, s, G and Λ² = μ_b⁴/(Q²κ_u) for us, and with the hard scale entering the μ_b-regions only through
   the one rapidity logarithm of a one-loop graph, (g/ε + g₁) ln(Q²/μ_b²) (boost invariance). Requiring b-independent poles
   (IR factorisation of the 2 → 2 amplitude and of H₁, H₂) and expanding in X, X_μ gives

       r₂ = c_u/2,      r₁ = g₁ + d_u + c_u ln κ_u,

   and the coefficients of X_μ² and XX_μ reproduce Δ_C and −2Δ_C. At this one-loop order, only ultrasoft regions can make r₂;
   with them cancelled, r₂ = 0 and r₁ = g₁, the ε⁰ term of the rapidity coefficient at μ = μ_b. That
   coefficient multiplies the Fourier transform of the eikonal 1/ℓ_⊥² in 2 − 2ε dimensions,
   2a_sΔ_C(μ²b²/b₀²)^ε e^{−εγ_E}Γ(−ε) = 2a_sΔ_C(μ²/μ_b²)^ε(−1/ε − ζ₂ε/2 + …): no ε⁰ term in the b₀
   convention, so **r₁ = 0**, as for the one-loop Collins–Soper kernel. The cancellation of the
   regulator dependence between soft and collinear regions forces the collinear coefficient of
   ln(ν/p⁺) to equal the soft one at all orders in ε, so g₁ = 0 holds for the total coefficient.
   r₁ = 0 refers to b₀ = 2e^{−γ_E} in X and X_μ, the convention of Φ_ν (eq. 5.5): replacing b₀ by
   b₀e^{−κ'/2} in both moves 2Δ_C κ' X from the renormalisation-group terms into r₁X without changing
   the amplitude. In the explicit two-loop graph, at fixed ℓ_⊥ ~ λ, the single logarithm of λ has the
   eikonal coefficient to 3.3 × 10⁻³ (B2.2); this tests the region structure of one graph, not g₁
   (an O(ε) property of the transverse integral), not the cancellation (which needs all attachments)
   and not r₁ itself.

Correction to B0 below: its boost-invariance argument for r₂ = 0 overlooked the ultrasoft regions,
which know both p₁⁺ and p₂⁻ (σ₁σ₂ ∝ 1/s) and in a single graph do produce X² (as in the laboratory
triangle, whose "kernel" c₁ + c₂ + us has X² coefficient c_u/2 = −1/2, `kernel_logs.py` part 2). The
cancellation over attachments is what makes r₂ = 0. This also explains why the data of Paper 2 give
r₂ = 0 for gluon (gg → HH) and quark (gg → ZZ) splittings although the colour factors of the
individual ultrasoft pairs, −C_{R₁}, −C_{R₂}, Δ_C/2, Δ_C/2, sum to T_{p₁}·T_{p₂} = −C_a ≠ 0: without
the parent attachments r₂ would be proportional to C_a.

Glauber region (corrected after the review): at leading power the gluon propagator and the p₂
side do not depend on ℓ⁻ there. The parent attachments do not telescope (ℓ_⊥ ~ k_⊥ shifts the
splitting differently), but both of their ℓ⁻ poles lie below the real axis and they fall like 1/ℓ⁻²,
so their ℓ⁻ integral **vanishes** (`us_telescoping.py` (4a): 2 × 10⁻²⁸ relative) — they do not feed
r₀ (an earlier version of this note said they did). The daughter attachments fall like 1/ℓ⁻ at both
ends; with a cutoff symmetric in ℓ_z their ℓ⁻ integral tends to −iπ times the coefficient of the
tail, independent of σ₁ (checked at two values of K_⊥, (4b)) and, by boost invariance, of p₁⁺: no
logarithm of the hard scale. With an asymmetric regulator the tails would give ln p₁⁺: the statement
needs a regulator symmetric in the two light-cone directions, like the η regulator.

## B0. What fixes r₂ and r₁ (structural argument, as first written; see the correction above)

In the paper's basis at fixed μ, X = X_μ + ℓ with ℓ = ln(Q_*²e^{−iπ}/μ²), the one-loop kernel
divided by K̃₀ is α X_μ² + β X_μ ℓ + γ X_μ + r₂ ℓ² + r₁ ℓ + const. So

- r₂ is the coefficient of ln²(Q_*²/μ²) without L_b = X_μ, and r₁ that of ln(Q_*²/μ²) without L_b.
- Each sector depends on the hard scales only as follows: the collinear-1 function on p₁⁺ only
  through the rapidity regulator (boost invariance along n), the collinear-2 function on p₂⁻
  likewise, the soft function on neither. A one-loop graph has a single rapidity divergence, so
  logarithms of p₁⁺, p₂⁻ appear at most linearly. Hence r₂ = 0 at one loop **provided no region
  knows both p₁⁺ and p₂⁻**; the ultrasoft regions do, and r₂ = 0 needs their cancellation (step 2).
- The regulator dependence cancels in the sum; what is left is C(ε) ln(p₁⁺p₂⁻b²/b₀²) plus
  logarithms of the momentum fractions, with C(ε) the (common) coefficient function of the
  rapidity logarithm. r₁ is the O(ε⁰) constant of C(ε) at μ = μ_b: r₁ = 0 iff the one-loop rapidity
  coefficient has no constant term in the b₀ convention (the analogue of the vanishing one-loop
  constant of the Collins–Soper kernel), plus the ultrasoft terms d_u, c_u ln κ_u, which cancel.

## B1. The soft function (done, `soft_function.py`)

Four incoming Wilson lines at (0, 0, b, b); at one loop only the cross pairs (k, r̄), (k̄, r)
contribute, colour Δ_C/2 each. With the δ-regulator (eikonal propagators 1/(ℓ⁻ + iδ₁),
1/(−ℓ⁺ + iδ₂)), the virtual exchange at separation b is done exactly (ℓ⁺ by residues, ℓ⁻ in
closed form, transverse Fourier transform in 2 − 2ε dimensions):

    S⁽¹⁾ = 2a_sΔ_C (μ²b²/b₀²)^ε e^{−εγ_E} Γ(−ε) [ψ(−ε) + γ_E − ln(δ₁δ₂b²/b₀²)]
         = 2a_sΔ_C [ −1/ε² − L_b/ε − L_b²/2 + ln(δ₁δ₂b²/b₀²)(1/ε + L_b) + ζ₂/2 ] + O(ε).

Every step is checked numerically in the script. Remarks: (i) the rapidity logarithm comes with
exactly 1/ε + L_b, no constant (this is g₁ = 0 of step 3); (ii) the constant ζ₂/2 × 2a_sΔ_C is a
first piece of r₀ (in N = 4, κ₁ = ζ₂ in units of a = 2a_sΔ_C: the soft function gives half);
(iii) with real δ's the soft function is real; with the physical sign of the off-shellness of
incoming lines the logarithm ln(δ₁δ₂) acquires the phase, and the Coulomb phase 4iπΔ_C of γ_K is
attached to the hard scale as in X = ln(Q_*²e^{−iπ}b²/b₀²) — a check for gate A3, not needed for r₁, r₂.

## Finding (15–16 October): the ultrasoft region

The daughters k, k̄, r, r̄ are off shell by O(λ²s). For a gluon exchanged between a daughter and
anything else, the leading-power small-ℓ regions are the soft mode (λ, λ, λ), where the daughter
propagators are eikonal and b is resolved (ℓ_⊥b ~ 1), AND the ultrasoft mode (λ², λ², λ²), where
(k + ℓ)² ≈ k² + k⁺ℓ⁻ keeps the off-shellness and b is not resolved. This is the soft region of the
off-shell Sudakov form factor (B2.1). In the DPS kernel it is present graph by graph (B2.2) and
cancels in the sum over attachments to the parent and the two daughters (B2.3). Consequence for the
definitions note: Definition 1 lists the ultrasoft region, with the statement that its sum over
attachments vanishes at one loop (note v4).

## B2.1 Laboratory (done, `lab_triangle.py`)

The one-loop massless triangle with two off-shell legs, p₁² = −a, p₂² = −b, a, b ≪ |s|: the
regions hard, collinear-1, collinear-2 and ultrasoft are each computed in closed form in
d = 4 − 2ε; the (λ, λ, λ) soft and the Glauber regions are scaleless. The poles cancel in the sum,
which equals (1/s)[ln(a/(−s−i0)) ln(b/(−s−i0)) + π²/3]: checked against the exact
Feynman-parameter integral to 1.4 × 10⁻³ (power corrections) for spacelike s, and against the
analytic continuation for timelike s. The Coulomb phase of the timelike case comes entirely from the
ultrasoft region, where the pole ℓ⁺ = b/p₂⁻ lies on the contour.

## B2.2 One explicit two-loop graph (done, `model_numerics.py`)

The cross-pair exchange of the two-loop DPS term (gluon from r̄ to k), scalar lines, numerator
N = k_⊥·(k + ℓ − q)_⊥, x = 0.75, k_⊥ = λ(0.6, 0.3), q = λ(1, 0). The DPS loop is reduced to its
residues (k⁻: D₃, and D₁ for k⁺ < 0; k⁺: D₄, and D₅ for ℓ⁻ < −x; `model_kres.py`, checked against
quadrature), the ℓ⁻ integral is done in closed form (`model_lminus.py`), the ℓ⁺ integral numerically
with finite i0 = 10⁻⁴λ² and windows at the pinches (`model_fast.py`; `model_exact2.py` in mpmath).
The ℓ⁺ partition has breakpoints at the pinches, at B = xℓ⁺ + a₆ = 0 and at the kinks (poles crossing
the end points of the ℓ⁻ ranges), with windows from 10⁻³ to 10⁵ i0 around each (first set: windows
around the pinches only, which missed the real part of F next to t = ℓ_⊥ + k − q = 0). On this
partition the two evaluators (Gauss–Legendre in double precision; tanh-sinh at 30 digits) agree to
10⁻¹¹–10⁻¹² at four points, 1.3 × 10⁻⁶ at a hard point (|ℓ_⊥| = 30λ) and 5.4 × 10⁻⁵ next to t = 0, where
64 Gauss–Legendre nodes are the limit (`model_ref_check.log`).

Bug found and fixed on the way (16 October, 02:00): the pole z₆′ of 1/D₆(u₅) lies **above** the
real ℓ⁻ axis (D₆(u₅) = −B + xa/(x + ℓ⁻) + iη with η > 0 on x + ℓ⁻ < 0), not below as first written. With
the wrong side, z₆′ and the gluon pole formed a fake pinch on the half line ℓ_⊥ ∥ q − k, beyond
q − k, where F grew like 1/distance and the transverse integral did not converge. The earlier test of
the closed form did not see it because the residue at z₆′ is small at its test points. The three
the pre-correction sources are retained in Git history at commit `1259ee1`, rather than in the active source tree. They are not used to generate the published results.

Results (λ ≥ 10⁻³; at λ = 3 × 10⁻⁴ the smallest windows reach the rounding level of double precision):

- ultrasoft transverse momenta ℓ_⊥ = cλ²: λ⁶F tends to a non-zero constant (c = 0.1, 1, 3; changes
  between successive λ shrink like λ), so the us region is leading power in this graph;
- ℓ_⊥ = λℓ̂: λ⁴F = f₀(ℓ̂) + f₁(ℓ̂) ln λ, a single logarithm (the rapidity logarithm, from
  λ² < |ℓ⁺| < 1); the real part has none (changes ≤ 3.4 × 10⁻⁴);
- f₁ is the eikonal prediction −2(2πi)³N(ℓ)/(x²k²a(ℓ)ℓ_⊥²)λ⁴ (the lowest-order integrand at shifted
  momentum times 1/ℓ_⊥², i.e. in b space the Fourier transform of 1/ℓ_⊥² times K̃₀) to 2.5 × 10⁻⁴ …
  3.3 × 10⁻³ at five values of ℓ̂ (power corrections).

Also seen: the Glauber pinch z₆ = z₂ at ℓ⁺ ≈ −a₆/x − a/(1 − x), ℓ⁻ ~ λ²; a second pinch at
k̄⁺ = 1 − x + ℓ⁺ ≈ 1.8λ² with a hard gluon (ℓ ≈ (−(1−x), 1−x, λ)) and soft k̄ and r, which carries
about 1 % of F and no logarithm in this scalar model and is suppressed by the effective vertices in
QCD (the ggH vertex is linear in each gluon momentum, so a soft k̄ with all components ~ λ² costs λ²;
a soft quark brings a spinor suppression); and a Landau singularity where the two
pinches merge, at |t| ≈ 0.2 in absolute units, i.e. at hard ℓ_⊥, power suppressed.

## B2.3 The ultrasoft cancellation (done, `us_telescoping.py`)

Checks: (1) the residue formulas of the three attachments (k, k̄, parent p₁) against direct
quadrature in k⁻ with finite i0, ≤ 1.7 × 10⁻⁶; (2) in the ultrasoft scaling ℓ = λ²ℓ̂ each attachment
tends to its limit and both colour components of the sum (the T_k part w_a(a) + w_c(c), the T_k̄ part
w_b(b) + w_c(c)) tend to tree/ℓ⁻, with corrections falling like λ (6.7 × 10⁻⁵ at λ = 10⁻⁴), while one
attachment alone is O(1) in these units; (3) the cancellation needs the ultrasoft scaling: in the
Glauber scaling the two components deviate by 0.34 and 0.10, independent of λ, and in the soft scaling
the parent attachment is suppressed by λ; (4a) in the Glauber region the ℓ⁻ integral of the parent
attachment vanishes; (4b) that of a daughter attachment, with a symmetric cutoff, tends to −iπ times
its tail coefficient, independent of σ₁; (4c) the parent attachment falls like 1/ℓ⁻², the daughter
attachment like 1/ℓ⁻; (5) the colour bookkeeping: the ultrasoft integral of one pair is
J = −2iπ²/ε² + …; given that an integral with a zero off-shellness is scaleless, the sum over the 3 × 3
attachments to {p₁, k, k̄} × {p₂, r, r̄} vanishes identically, while the daughters alone give
T_{p₁}·T_{p₂} J ≠ 0 (the content is in (1) and (2); (5) only adds up the colours).

## B2.4 The bookkeeping (done, `kernel_logs.py`)

Part 1: the general region ansatz, b-independent poles (g = −c − 2c_u), finite part in X, X_μ:
X_μ²: (c + c_u)/2, XX_μ: −c − c_u, X²: c_u/2, X: g₁ + d_u + c_u κ_u; matching the
renormalisation-group terms gives c + c_u = 2Δ_C and r₂ = c_u/2, r₁ = g₁ + d_u + c_u κ_u (the script's
kappa_u is ln κ_u of step 3). Part 2: on
the laboratory triangle the "kernel" c₁ + c₂ + us has c_u = −1 and X² coefficient −1/2 = c_u/2.
Part 3: the DPS kernel, c_u = d_u = 0 (B2.3), g₁ = 0 (B1): r₂ = r₁ = 0.

## What remains

1. **Write-up**: done in Paper 2 (Sec. 5 and appendix E, 16 October); the definitions note v4
   (Definition 1 with the ultrasoft region and its cancellation; the bookkeeping lemma) is still to do.
1b. **End-to-end check (open)**: a model with all nine attachments and gauge (eikonal) couplings,
   transformed to b space, in which r₁ and r₂ are read off directly; the review regards this as
   what would turn "shown, given the list of regions" into "derived".
2. **Gate C (optional, quantitative)**: r₀ from the finite parts (soft function ζ₂/2, collinear
   functions with the DPS numerators of eq. (2.3) for both helicity components, Glauber parent terms),
   against the G coefficient of M₁⁽²⁾ and the constant of M₂⁽²⁾ in gg → HH and gg → ZZ by colour factor;
   and the N = 4 constant κ₁ = ζ₂. Both done: `../r0/theory/` (Sec. 7.2 and 7.3; κ₁ = ζ₂ on 18 October).
3. Gates A2/A3 at the level of the collinear functions (poles, the phase 4iπΔ_C, the scale Q_*): not
   needed for r₁, r₂, which follow from the structure above; part of gate C.

## Files

| file | what it does |
|---|---|
| `soft_function.py`, `.log` | B1: one-loop soft function of four Wilson lines, δ-regulator, every step checked numerically |
| `lab_triangle.py`, `.log` | B2.1: regions of the off-shell triangle against the exact result |
| `model_kres.py`, `.log` | B2.2: residues of the DPS loop (k⁻ and k⁺) against quadrature |
| `model_lminus.py` | B2.2: closed-form ℓ⁻ integral (z₆′ side corrected) |
| `model_exact2.py` | B2.2: mpmath evaluator of F(ℓ_⊥) (z₆′ side corrected) |
| `model_fast.py` | B2.2: numpy evaluator (pinches, kinks, windows) used for the results |
| `model_ref_check.py`, `.log` | B2.2: numpy against mpmath at six points |
| `model_numerics.py`, `.log` | B2.2: ultrasoft scaling, single logarithm, eikonal coefficient |
| `model_2d.py`, `model_2d_conv.log` | transverse integral of F (written before the z₆′ fix was complete; not used for the results) |
| `us_telescoping.py`, `.log` | B2.3: the ultrasoft cancellation |
| `kernel_logs.py`, `.log` | B2.4: r₂ = c_u/2, r₁ = g₁ + d_u + c_u κ_u; triangle check; r₁ = r₂ = 0 |
| `stage_b_results.py`, `results_stage_b.json` | the numbers quoted in the manuscript, recomputed from the functions above |
| `model_exact.py` | first (scipy) version, superseded |

Run: `python3 soft_function.py`, `lab_triangle.py`, `model_kres.py`, `us_telescoping.py`
(each < 1 min), `kernel_logs.py` (seconds), `model_numerics.py` (about 1 min), `model_ref_check.py`
(about 4 min), `stage_b_results.py` (about 2 min).
