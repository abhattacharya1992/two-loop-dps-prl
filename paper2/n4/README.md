# The b-space structure against planar N = 4 (Dixon–Esterlis)

Paper 2, Sec. 5.1, Tables 5 and 6. (Stage D of the roadmap of task K; in the P1 package this
folder is `paper2/n4/`.)

    make test       reproduce results_stage_d.json and compare seven constants with published expressions (about 10 s)
    make analyse    the comparison: results_stage_d.json and stage_d.log
    make blind      the two independent reimplementations (about 6 minutes)

Only sympy is needed, no input files: the data are the printed coefficients of Dixon and
Esterlis, transcribed into `de_data.json`.

Test of the structure of Paper 2, Sec. 5, in the conformal theory, where exact results exist:
the self-crossing limit of the six-gluon amplitudes of planar N = 4 super Yang–Mills
(L. J. Dixon and I. Esterlis, JHEP 07 (2016) 116, arXiv:1602.02107, and the erratum JHEP 08 (2016)
131 for the NMHV amplitude), which mimics double parton scattering. Their variable is
δ = p_T²/Q_*² with Q_*² = s_36 x y (1−x)(1−y), the scale of Paper 2.

Structure tested, in b-space: a K₀(b) (1 + κ₁a + κ₂a² + …) exp[−(γ_K(a)/8) X² + B(a) X], with
γ_K = 4a − 4ζ₂a² + 22ζ₄a³ + γ_K⁽⁴⁾a⁴ (γ_K⁽⁴⁾ left free), transformed with Φ_ν of Paper 2, eq. (5.5):
ν = 0 for the MHV amplitude (|ΔJ_z| = 0), where G = ln|δ|/2 and C₁ = 2πi in 3 → 3 kinematics and
G = ln δ/2 + iπ/2, C₁ = −2πi in 2 → 4; ν = 2 for the (1+w)/(1+w*) part of the NMHV amplitude
(|ΔJ_z| = 2) in the v → 0⁻ self-crossing/MRK overlap, one-loop term K. Both generating functions are exp[−2Gα − Σ_{k odd ≥ 3}(2ζ_k/k)αᵏ],
times 1/(2α) for ν = 0 and 1/(1−α) for ν = 2 (checked against the Γ functions in the script).

| file | what |
|---|---|
| `de_data.json` | the coefficients of 2πi lnᵏ|δ| in E⁽ᴸ⁾, L = 1–5, eqs. (3.17)–(3.21) of Dixon–Esterlis, and of πi (1+w)/(1+w*) lnᵏ|δ| in ρ⁽ᴸ⁾, L = 1–4, k = 2L−2, …, 0, eqs. (D.14)–(D.17) of the erratum; Y and R₆⁽²⁾ for 2 → 4, eqs. (2.32), (3.10). Transcribed from the PDF texts two (MHV) and three (NMHV) times; see `_source` and `_transcription` inside |
| `stage_d.py` | the comparison: the one-loop terms define the normalisations C₁ and K, B⁽²⁾ and B⁽³⁾ are taken from independently published virtual anomalous dimensions, not fitted to the 2016 MHV coefficients. Six other parameters are determined from those MHV data, and the remaining coefficients, MHV and NMHV, are checks; also what the three-loop ln³ tests, the exact relation between the two components, the five-loop prediction for the |ΔJ_z| = 2 component and the 2 → 4 two-loop coefficients (about 10 s) |
| `results_stage_d.json`, `stage_d.log` | its output; `results_stage_d.json` is what the manuscript reads (Tables 5 and 6, Sec. 5.1) |
| `test_n4.py`, `Makefile` | the test and the targets above |
| `test_literature.py` | an independent symbolic comparison of κ₁–κ₄ and B⁽²⁾–B⁽⁴⁾ with the published expressions of arXiv:1903.10890, without fitting the MHV data |
| `stage_d_blind.py`, `.log` | separate reimplementation of the MHV comparison from the written definitions, written with the AI assistant in a fresh session (same results; about 30 s) |
| `nmhv_blind.py`, `.log` | the same for the NMHV comparison, written in another fresh session from `stage_d_blind.py` and the data only (same results; about 5 minutes, spent in the sympy series of the Γ functions) |

Result, MHV (3 → 3, 25 coefficients through five loops): one fixes C₁, six
fit parameters, and B⁽²⁾, B⁽³⁾ enter independently from the published
2019 virtual anomalous dimension. The other 18 are tests and all agree.
The external inputs are B⁽²⁾ = 3ζ₃ and B⁽³⁾ = −2ζ₂ζ₃ − 10ζ₅ in the
2016 BDS-like convention; test_literature.py checks them symbolically. In particular: the leading logarithms at L = 1–5
(C₁G∫₀¹dt e^{−2aG²t²}, the transform of the Gaussian); no next-to-leading logarithms at the tested orders
(B⁽¹⁾ = 0 and the scale Q_*); the three-loop ln³|δ| term, −(2/3)(κ₁ + γ_K⁽²⁾/4)C₁G³, vanishes after
κ₁ = ζ₂ (the first test of the two-loop cusp −4ζ₂); the ζ₂ tower at L = 4, 5; the ζ₃ tower at
L = 4, 5 after B⁽²⁾ = 3ζ₃; the ζ₄ tower at L = 4, 5 after κ₂ = ζ₄/2 (tests the three-loop cusp 22ζ₄);
the five-loop ln⁴ term after B⁽³⁾. The five-loop ln³|δ| coefficient then determines
γ_K⁽⁴⁾ = −(219/2)ζ₆ − 4ζ₃², the four-loop cusp anomalous dimension of Beisert–Eden–Staudacher. The
constants have uniform weight (2n for κ_n, 2n − 1 for B⁽ⁿ⁾). In 2 → 4 at two loops the ln³δ, ln²δ and
ln δ coefficients agree (the second, −π²/2, from the iπ/2 of G; the third, −2iπ³/3, independent of v).
Until 20 October 2026 the transcription read H₂^v of eq. (3.10) as Li₂(v); Dixon–Esterlis write H^v_w for
H_w(1−v) (their eq. (3.9) is the value at v = 1, where every H^v vanishes), so H₂^v = Li₂(1−v), and the
v-dependent excess iπ[Li₂(1−v) − Li₂(v)] reported earlier was an artefact of the transcription. The 2 → 4
single logarithm fixes the finite part of the Glauber region of the kernel to zero (`r0/theory/glauber.py`).

Result, NMHV (3 → 3, v → 0⁻, 16 coefficients through four loops, the constants k = 0 included): one
defines K, the other 15 are predictions with no new constant, and all agree, the constants
ζ₂ − 1, ζ₃ + ζ₄/2 + 3 and the four-loop one among them (the constants of the polynomials in
ln|δ| as grouped in the erratum; the ln|v| terms are outside the comparison). The leading
logarithms sum to K e^{−2aG²}. Within the ansatz the two components are related exactly, M₂(G)/K =
∫₀^∞ds e^{−s} (∂M₁/∂G)(G − s/2)/C₁ (checked algebraically through five loops, with NMHV data through four), and the structure predicts the
five-loop |ΔJ_z| = 2 coefficients (in `results_stage_d.json`, `nmhv_five_loop_prediction`, and
eq. (5.10) of the paper). The two singularities of Dixon–Esterlis are the two components of one
kernel.

History: the first two transcriptions of eq. (D.17) had distributed the minus sign in front of
the bracket of the four-loop ln²|δ| coefficient, −(7ζ₄/4 − ζ₃/2 + 3ζ₂/2 + 15/2), onto its first
term only, and that one coefficient appeared to differ from the prediction; a third reading,
term by term, corrected the transcription and added the constants. Nothing in the definitions
was changed.

The MHV representation and constants also follow from Caron-Huot et al.,
[JHEP 08 (2019) 016](https://doi.org/10.1007/JHEP08(2019)016),
[arXiv:1903.10890v2](https://arxiv.org/abs/1903.10890v2), Sec. 5. In this repository's 2016 BDS-like
convention, a = 2g², γ_K = 2Γ_cusp, B = −Γ_virt and
1 + Σ κ_n aⁿ = exp(ζ₂γ_K/4 + 2Γ₃), with the additional 2019 cosmic-normalization factor removed.
`test_literature.py` compares all seven constants with eqs. (2.20), (5.4), (5.10) and (5.11).

Consequences for Paper 2: the comparisons support the absence of an independent quadratic
boundary term in the exponent and of a one-loop linear term at the tested orders. The b-space
ansatz agrees with MHV data through five loops and with the specified NMHV component through
five loops: the nine five-loop coefficients have also been checked against an independent
endpoint expansion of the primary multi-Regge ancillary amplitude, including all odd-zeta
terms. See [PRIMARY_NMHV_CHECK.md](PRIMARY_NMHV_CHECK.md), `make primary-inputs` and
`make primary-check`. These comparisons do not establish
an all-order restriction on QCD boundary logarithms.
