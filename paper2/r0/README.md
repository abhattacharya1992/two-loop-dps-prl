# Finite matching in the two-loop DPS kernel

The finite matching functions $r_0$ are derived analytically from the one-loop
collinear contributions, the soft function and parton self-energy corrections.
The extraction from full two-loop amplitudes and numerical region calculations
provide independent checks. The Letter and Supplemental Material contain the
derivation and comparisons.

The calculations are organised as follows:

- `derivation/`: analytic diagram calculations and checks of the gluon and quark
  splitting contributions, including the helicity and parity dependence.
  Run `make -C paper2/r0/derivation test` from the repository root.
- `data/`: finite matching from the full two-loop amplitudes.
- `theory/`: numerical evaluation of the momentum regions.

The original scripts, numerical results and further checks are retained below.

## `data/`: finite matching from two-loop amplitudes

r₀ enters the coefficient of 𝒢 in 𝓜₁⁽²⁾ (ν = 0) and the constant of 𝓜₂⁽²⁾ (ν = 2). With the hard functions,
the kernel evolution and conversion of the finite remainders to minimal subtraction,
the extracted imaginary parts are consistent with zero and the complex
matching coefficients agree with the analytic expressions to 5.4 × 10⁻¹⁰
(gg → HH) and 1.6 × 10⁻⁹ (gg → ZZ).

    gluon splittings: r₀ = C_A[134/9 − 5ζ₂ − 4 ln x ln(1−x)] − (40/9)T_F n_f + (2/3)(C_A − 2T_F n_f)/F_ν(x)
    quark splittings: r₀ = C_F[12 − 6ζ₂ + 2 ln²(x/(1−x))] + C_A[2 + ζ₂ − 2 ln²x − 2 ln²(1−x)]

Here F₀ = F_S and F₂ = W, as defined in the Letter.

| file | what it does |
|---|---|
| `r0_from_data.py` → `r0_from_data.json` | r₀ from the gg → HH scan points of `../hh/hh_two_loop/points` (Catani's scheme at μ = m_H), N and n_f parts, both components |
| `r0_from_zz.py` → `r0_from_zz.json` | r₀ from the gg → ZZ fits of `../zz/results_zz.json` (q_T scheme at μ² = s), C_F and C_A parts, all eight amplitudes |
| `r0_closed_forms.py` → `r0_closed_forms.json` | the closed forms against both, largest deviation |
| `points_hh/`, `run_hh_new_velocities.py` | new gg → HH points at β = 1/10 and 4/5 (x = 0.55, 0.9), with the scan code of `../hh/hh_two_loop` (needs the files of `make fetch` in `../hh`) |
| `points_zz/`, `zz_jobs_a.txt`, `zz_jobs_b.txt` | new gg → ZZ points at β = 4/5 (x = 0.9), computed with `../zz/scripts/zz_points.py` and these job lists |
| `blind_test_hh.py`, `blind_test_zz.py` → `blind_test_*.json` | the closed forms against fits at the new velocities, which were not used to find them: 2 × 10⁻⁹ (HH), 4 × 10⁻⁹ (ZZ) |
| `points_endpoint/`, `zz_points_acc.py`, `endpoint_test.py` → `endpoint_test.json` | gg → HH and gg → ZZ at β = 19/20 (x = 0.975, where the x-dependent part of the quark-splitting r₀ is 2.4 times its value at x = 0.9): r₀ against the closed forms, every value within a tenth of its fit-window shift, and the precision controls. The ZZ points come from `zz_points_acc.py`, which runs `../zz/scripts/zz_points.py` with Li₂,₂ on its integral branch evaluated by `../hh/hh_two_loop/evaluate_hh.li22`: at this velocity the fixed subdivision of `ggvv.py` loses up to 29 of 90 digits (`points_endpoint/zz/li22_*.json`; `check_li22_scan.py` stops at δ = 10⁻⁴, `points_endpoint/logs/check_li22_scan.log`). `points_endpoint/zz_fixed_split/` holds the same points from the unchanged code; their amplitudes differ by less than 10⁻²². The commands are in the docstring of `endpoint_test.py` |
| `test_data.py` | reruns the six scripts and checks that the stored files are reproduced |

## `theory/`: numerical evaluation of momentum regions

The one-loop kernel at μ = μ_b in impact-parameter space, K̃⁽¹⁾/K̃₀ = J₁ + J₂ + S + Σ (see "One-loop correction to the DPS kernel" in the Supplemental Material): the collinear
functions of the two incoming gluons in Feynman gauge, as one-loop light-front wave functions of the gauge-invariant
fields (χ_n, χ̄_n for g → qq̄, B_⊥ for g → gg) with the η regulator on the Wilson lines; the soft function with the
same regulator; the self-energy corrections to the four parton lines. The light-cone energies are integrated by residues, the
transverse momentum with Feynman parameters, the rest by sector decomposition; the Laurent coefficients in η and ε
come from discrete Fourier sums on circles. At x = 0.1, 0.3, 0.5 the 1/η poles cancel, the real parts of the
poles are those of 2I_K, and the finite parts agree with the closed forms above, colour factor by colour factor
and in both components, to 10⁻⁷ (the corresponding table in the Supplemental Material). The collinear and soft functions and the self-energies are real; the
imaginary part of the kernel is the Glauber region of the cross-pair exchange (`glauber.py`, workbook VII.15,
20 October): the hard scales enter the kernel only through its rapidity logarithm, as ln(Q_*²e^{−iπ}/μ_b²) with
the phase of −s − i0 of the two incoming parents, and the Glauber region is −iπ times the coefficient of that
logarithm, −2πiΔ_C e^{−εγ}Γ(−ε)(μ²b²/b₀²)^ε = 2πiΔ_C[1/ε + X_μ + ε(X_μ²/2 + ζ₂/2)]: the pole of 2I_K, the X_μ
term of the RG bracket, and no finite part at μ = μ_b (the same mechanism as r₁ = 0). It comes from the daughters'
light-cone energy denominators, which put the eikonal pole of the ℓ⁻ integral inside the range (for an incoming
and an outgoing parent, 3 → 3, the region is real); the η-regulated soft and collinear functions carry no phase.
Checked on the explicit two-loop graph of `../kernel` (Re f₀ = −(π/2) Im f₁ at eighteen configurations, to
10⁻⁴), against Im r₀ = 0 of the data (104 fitted entries (including symmetry-related helicity components), within the fit uncertainties), and in planar N = 4 in
2 → 4 kinematics, where the single poles cancel only with the kernel's iπ/ε and Dixon–Esterlis' two-loop single
logarithm, −2iπ³/3 (their H^v_w = H_w(1−v); the transcription of 15 October read H₂^v as Li₂(v), which produced a
spurious v-dependent term), fixes the finite part of the Glauber region to zero. The helicity-dependent term of the
gluon case comes from the rational terms of the one-loop splitting amplitudes: the helicity flip for ν = 2, a
term x(1−x)(C_A − 2T_F n_f)/3 times the tree, when both daughters keep the helicity of the parent, for ν = 0.

| file | what it does |
|---|---|
| `lc_engine.py` | Dirac projection on the good components, light-cone energy integrals by residues, transverse integrals (g → qq̄) |
| `vec_engine.py`, `gg_driver.py` | Lorentz tensors in light-cone components and Dirac traces; residues and transverse integrals of scalar products (g → gg) |
| `gqq_diagrams.py`, `gg_diagrams.py` | the graphs: Dirac strings or tensors, light-cone denominators, factors, colour (SU(3), explicitly) |
| `run_gqq.py`, `run_gg.py`, `driver.py` | the parametric integrands per interval of the loop momentum fraction (`*_integrands_x*.pkl`) |
| `factors.py`, `secdec.py` | factorisation into irreducible polynomials, sector decomposition, subtraction, Gauss–Legendre |
| `assemble_gqq.py`, `assemble_gg.py` | Laurent coefficients graph by graph (`gqq_laurent_x*.json`, `gg_laurent_x*.json`) |
| `combine_gqq.py`, `combine_gg.py`, `r0_theory_quarks.py`, `r0_theory_gluons.py` | the kernel and the comparison with the closed forms (`r0_theory_quarks.json`, `r0_theory_gluons.json`) |
| `glauber.py` → `glauber.json` | the imaginary part of the kernel (task 5): the b-space result, the end point of the eikonal integral, the off-shell eikonal pair in closed form against quadrature, its four η-regulated regions (part 2b, `glauber_regions.py`), the explicit two-loop graph of `../kernel` against its Glauber region, Im r₀ from the stored fits of `../data` and `../../../common/dps`, and planar N = 4 in 2 → 4 kinematics (about ten seconds) |
| `glauber_regions.py` | the off-shell exchange expanded by regions with the rapidity regulator at each region's scaling, all four in closed form (the soft region in light-cone polar coordinates): the poles cancel, the sum is the exact exchange, the Glauber region alone is its −iπ term (the δ-function piece plus the cross term of the two principal values, which the symmetric regulator couples), and for an outgoing parent the collinear function cancels it; workbook VII.15 (c′) (seconds) |
| `three_loop.py` → `three_loop.json` | the conditional three-loop prediction (the Letter and Supplemental Material): the kernel evolved to O(a_s³) with the two-loop cusp and non-cusp anomalous dimensions, β₀, β₁, r₀ and the boundary terms r₀⁽³⁾, r₁⁽³⁾; the transforms of eq. (5.5) and the hard functions; the coefficients of G⁵, G⁴, G³ (ν = 0) and G⁴, G³, G² (ν = 2), which contain no open quantity within the boundary ansatz described below, in general form and for gg → HH at μ = m; checks the two-loop terms against eqs. (5.3), (5.6), (7.1) and the conformal limit (seconds) |
| `test_theory.py` | recombines the stored Laurent coefficients and checks the stored results and the statements of the paper, reruns `glauber.py` and `three_loop.py` |
| `test_secdec.py`, `test_gg_tree.py`, `check_*.py`, `compare_stored.py` | the checks of the Supplemental Material: sector decomposition on known integrals, the g → gg tree and F₀, F₂ (gate A1), the normal form against the integrands, charge conjugation, sectors against direct quadrature; `check_helicity_terms.py` → `helicity_terms.json`, where the helicity-dependent term of r₀ comes from (Sec. 7.2); `check_formulas.py`, the analytic steps of Sec. 7 and the Supplemental Material one by one, as described in `derivation/README.md` |
| `*_nq96-72-40.json`, `*_neps16.json` | the precision runs at x = 1/2: quadrature orders raised by half, 16 instead of 32 points on the ε circle |

The three-loop results assume no independent quadratic or higher-power boundary logarithm at
that order. One-loop ultrasoft cancellation does not prove this restriction on the full kernel.
An additional r₂⁽³⁾X² term would change the quoted G³ coefficient for ν = 0 and G² coefficient for
ν = 2; the μ renormalisation-group equation alone does not exclude it.

### Planar N = 4: κ₁ = ζ₂ from the regions (Sec. 7.3, the corresponding planar table)

The same kernel with the fields of N = 4 super Yang–Mills in the loops: the fermion graphs with T_F n_f → 2C_A (four
adjoint Weyl fermions), the scalar triangle of a complex adjoint scalar (`n4_diagrams.py`; the seagull graphs vanish by
colour) with multiplicity 3 + ε (the three complex scalars and the 2ε ε-scalars of dimensional reduction, which make
the scheme the four-dimensional-helicity scheme of the BDS normalisation), and the self-energies, which add up to −2C_A
exactly. With independent momentum fractions x₁, x₂ the kernel has the poles 2C_A/ε² + 2C_A L_*/ε (β₀ = 0) and the finite
part C_A[−5ζ₂ − 2 ln x₁ ln(1−x₁) − 2 ln x₂ ln(1−x₂)], the same for both components (1 − 4 + 3 = 0). Together with the
one-loop four-gluon amplitudes of the two hard scatterings and the one-loop BDS-like normalisation of Dixon–Esterlis,
in 3 → 3 kinematics, the poles cancel, the finite part is −X²/2 + ln x ln(1−x) + ln y ln(1−y) + 7ζ₂/2 + r₀^{N=4}/(2C_A),
and κ₁ = ζ₂ to 2 × 10⁻⁷ for the fifteen pairs of fractions and both components.

| file | what it does |
|---|---|
| `n4_diagrams.py` | the scalar graphs (triangle, two orientations; seagull bubbles, zero by colour), their colour factors per complex adjoint scalar, and the scalar self-energy coefficient −1/(3−2ε) |
| `run_n4.py`, `assemble_n4.py`, `assemble_lib.py` | the integrands and Laurent coefficients of the scalar graphs (`n4_integrands_x*.pkl`, `n4_laurent_x*.json`), as `run_gg.py`/`assemble_gg.py`; with `--fdh`, the gluon-loop graphs with g^μ_μ = 4 (`fdh_laurent_x1o2.json`) |
| `check_fdh.py` → `fdh_check_x1o2.json` | the scheme: gluon-loop graphs with four-dimensional metric traces = conventional dimensional regularisation + ε × one complex adjoint scalar, to 2 × 10⁻¹³ |
| `check_bdslike.py` | eqs. (7.5) and (7.6): the one-loop four-point amplitude and the BDS-like function M̂₆ from the box functions of Bern–Dixon–Dunbar–Kosower and Y of Dixon–Esterlis (three-particle invariants cancel; constants 4ζ₂ and 9ζ₂) |
| `combine_n4.py` | the N = 4 kernel at (x₁, x₂): the three colour weights, the ε-scalars, the self-energies |
| `n4_hexagon.py` → `kappa1_n4.json` | the symbolic assembly (poles cancel, finite part −X²/2 + … + 7ζ₂/2, no dependence on the scattering angles) and κ₁ for all pairs of fractions |
| `check_formulas.py`, check 13 | the self-energies and the assembly, as the workbook (VII.14) states them |

    cd theory && make n4        # about ten minutes on one core: the scalar graphs at the five x, the scheme check, kappa1_n4.json

    cd theory && make test      # about fifteen seconds (glauber.py included)
    cd theory && make glauber   # about ten seconds: the imaginary part of the kernel, glauber.json and glauber.log
    cd theory && make threeloop # seconds: conditional three-loop coefficients, three_loop.json and three_loop.log
    cd theory && make checks    # about three minutes
    cd theory && make quarks    # about half an hour: g -> q qbar from scratch, compared with the stored results
    cd theory && make gluons    # about two hours: g -> gg from scratch
    cd theory && make n4        # about ten minutes: the scalar graphs of planar N = 4 and kappa_1 (Sec. 7.3)
