# Equations used in stage2a_bookkeeping.py

Transcribed from the sources. Equation numbers of Banerjee et al. are those of the arXiv v2
version; check them against the journal version before quoting.

## Banerjee, Borowka, Dhani, Gehrmann, Ravindran, JHEP 11 (2018) 130 [arXiv:1809.05388]

d = 4 − 2ε, a_s = α_s(μ_R²)/(4π).

Coupling, Eq. (2.15):

    â_s μ^(2ε) S_ε = a_s μ_R^(2ε) [ 1 − a_s β₀/ε + a_s² (β₀²/ε² − β₁/(2ε)) + O(a_s³) ],
    S_ε = exp[(ln 4π − γ)ε],
    β₀ = (11/3) C_A − (4/3) T_F n_f,   β₁ = (34/3) C_A² − 4 T_F n_f C_F − (20/3) T_F n_f C_A.

Operator renormalisation: M^B_i = Z_O² M̂^B_i + Z^L_11 M̂^(A,(0))_i |_(λ=0), with
Z_O = 1 − a_s r_(O1;1)/ε + a_s² (r_(O2;2)/ε² − r_(O2;1)/ε), r_(O1;1) = β₀, r_(O2;2) = β₀²,
r_(O2;1) = β₁, and Z^L_11 = a_s² β₁/ε.

Renormalised amplitudes, Eq. (2.23):

    M^(B,(1))_i = μ_R^(2ε) [ M̂^(B,(1))_i − μ_R^(−2ε) (2 r_(O1;1)/ε) M̂^(B,(0))_i ]
    M^(B,(2))_i = μ_R^(4ε) [ M̂^(B,(2))_i − μ_R^(−2ε) ((2 r_(O1;1) + β₀)/ε) M̂^(B,(1))_i
                 + μ_R^(−4ε) ((r_(O1;1)² + 2 r_(O2;2))/ε² − 2 r_(O2;1)/ε) M̂^(B,(0))_i
                 + μ_R^(−4ε) (β₁/ε) M̂^(A,(0))_i |_(λ=0) ]

so the term acting on the one-loop amplitude is −3β₀/ε.

IR factorisation, Eq. (2.24):

    M^(B,(1))_i = 2 I^(1)_g(ε) M^(B,(0))_i + M^(B,(1),fin)_i
    M^(B,(2))_i = 4 I^(2)_g(ε) M^(B,(0))_i + 2 I^(1)_g(ε) M^(B,(1))_i + M^(B,(2),fin)_i

Subtraction operators, Eq. (2.25):

    I^(1)_g(ε) = − e^(εγ)/Γ(1−ε) · ( C_A/ε² + β₀/(2ε) ) · ( −μ_R²/s )^ε
    I^(2)_g(ε) = −(1/2) I^(1)_g(ε) [ I^(1)_g(ε) + β₀/ε ]
                 + e^(−εγ) Γ(1−2ε)/Γ(1−ε) · [ β₀/(2ε) + K ] · I^(1)_g(2ε) + 2 H^(2)_g(ε)
    K = (67/18 − ζ₂) C_A − (10/9) T_F n_f

Only I^(1) matters for the DPS logarithms; I^(2) multiplies M^(0).

## P1 (Bhattacharya, gg → AA threshold paper)

d = 4 + ε, following Ref. [8] (Bhattacharya, Mahakhud, Mathews, Ravindran, JHEP 02 (2020) 121).

- Eq. (2.13): M_i = M_i^(0) + a_s M_i^(1) + a_s² M_i^(2) + …, a_s renormalised at μ_R.
- Eq. (2.14): M^(1) = 2 I^(1) M^(0) + M^(1),fin;  M^(2) = 4 I^(2) M^(0) + 2 I^(1) M^(1) + M^(2),fin.
- Eq. (2.16): M_i^(1) = M_(i,publ)^(1) + 2β₀ ln(μ_R²/GeV²) M_i^(0);  M_i^(2) = M_(i,publ)^(2)/64.
- App. B: "The original calculation renormalises the amplitude with Z_GG², where Z_GG = Z_as
  and C_G ∝ a_s."
- Eq. (B.1): ∂M^(1)/∂ln μ_R² = 2β₀ M^(0).
- Eq. (B.3): ∂M^(2)/∂ln μ_R² = 3β₀ M^(1) + (2β₁ − 8H_g^(2) − 12 C_F n_f) M^(0).

In d = 4 + ε the same operator reads
I^(1) = − e^(−εγ/2)/Γ(1+ε/2) · (4C_A/ε² − β₀/ε) · (−s/μ_R²)^(ε/2), and the UV term acting on
M^(1) is +6β₀/ε. Substituting ε → −2ε gives the forms above; the script checks that the
O(ε⁰) terms agree.
