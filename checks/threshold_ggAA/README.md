# Threshold validation of the gluon-fusion DPS formula

This is a targeted symbolic comparison with A. Bhattacharya,
*Two-loop gg -> AA amplitudes near threshold and the double parton scattering
singularity*, arXiv:2610.19024v2.

Run from the repository root:

```bash
python3 checks/threshold_ggAA/check.py
```

The checker transcribes the independently published coefficients of
Eqs. (4.5), (4.6) and (4.7) and compares them with the PRL formula at
`beta_A -> 0`, `x = 1/2`, `s = 4*m_A^2`, `mu_R = m_A`, with coefficients
normalised to `m_A^2`. It verifies the one-loop equal-helicity logarithm
`C1`, the opposite-helicity constant `K`, four two-loop leading and
subleading logarithms, and the two-loop equal-helicity **single logarithm**.
The latter check includes the finite subtraction-scheme conversion and
the pseudoscalar hard form factor. The scalar/pseudoscalar weight ratio
`F_S/F_A=9/7` is also checked.

All equalities hold symbolically for arbitrary `N_c` and `n_f`.

**Scope:** These are checks between two kinematic analyses of the *same*
underlying two-loop amplitude, not independent raw-amplitude evaluations.
Agreement at `x=1/2` does not prove the complete matching function at
arbitrary `x` or the completeness of the factorisation theorem.
The earlier complete `beta_A^6` series and its validation package belong to
the earlier paper and need not be copied into this PRL repository.
