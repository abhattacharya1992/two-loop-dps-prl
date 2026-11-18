# Links of the DPS term to other vector-boson pairs and to DPS cross sections

Two exact checks behind statements of the Letter and the Supplemental Material. Both use the Python
standard library only and run in a few seconds.

```bash
python3 checks/dps_links/check_couplings.py
python3 checks/dps_links/check_dpd_kernels.py
```

## `check_couplings.py`

Numerator of the massless quark loop in `g g -> V1 V2`, with the vector bosons coupled through
`gamma^mu (L P_L + R P_R)`, for arbitrary loop momentum and arbitrary vectors at the four vertices
(explicit Dirac matrices, Gaussian rationals). Mixed-chirality terms vanish, the term with one
`gamma_5` changes sign when the fermion flow is reversed, and the terms with two axial-vector couplings
equal those with two vector couplings. The sum of the two orientations is `(L1 L2 + R1 R2)` times a
function that does not depend on the couplings. The same holds with a gluon attached twice to the loop.

This is the reason why the coefficient functions of von Manteuffel and Tancredi (arXiv:1503.08835), which
`paper2/zz` and `paper2/unequal_virtuality` evaluate without electroweak couplings, describe
`gamma* gamma*`, `Z gamma*`, `ZZ` and `W+ W-` production through massless quark loops.

## `check_dpd_kernels.py`

The one-loop DPS functions `F_S`, `F_A` and `W` of the Supplement are linear combinations of the
lowest-order kernels for the splitting of a gluon into two unpolarised, two longitudinally polarised
and two linearly polarised gluons in double parton distributions (Diehl, Kasemets and Keane,
arXiv:1401.1233, Appendix B). The script derives the kernels from the tree splitting amplitudes in the
helicity basis, compares them with the printed forms and with three values quoted by Diehl, Gaunt and
Schoenwald (arXiv:1702.06486), and verifies the relations of the Supplement, the helicity average of
the squared one-loop DPS term, the corresponding relations for quark splittings, and the statement
that the phase of the kernel enters its modulus squared only at relative order `a_s^2`.
As a control, the same helicity sum for two scalars coupled to massless quarks gives `(1 + beta^4)/4`,
the dependence on the velocity of the subtraction term in Eq. (9.22) of Diehl, Gaunt and Schoenwald.

## Scope

The checks are algebraic. They do not compute a cross section, its normalisation or the colour
channels of the double parton distributions, and they contain no real emission. No relation between
the finite matching function `r_0` and the two-loop splitting kernels of double parton distributions
is derived or assumed.
