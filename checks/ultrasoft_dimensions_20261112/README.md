# The ultrasoft identities in d dimensions (12 November 2026)

Records of `checks/ultrasoft_numerators.py` with the vectors, the metric contractions and the Dirac
matrices in more than four dimensions, and of `checks/ultrasoft_dimension_dependence.py`, which compares
them. The derivation note is `paper2/derivations/ultrasoft_and_three_loop_boundary.md`, Sec. 8, item 4.
These JSON files preserve the original higher-dimensional checks. Their
original source scripts are not included in this repository, so the
calculations cannot currently be rerun from the retained records alone.
They are historical numerical and symbolic evidence rather than an
executable regression test.

| File | Run | Result |
|---|---|---|
| `n6_d6.json` | algebra in 6 dimensions, vectors in 6; sections A to F; 3, 2, 2, 2 points | every identity of the four-dimensional run holds; the same suppression powers |
| `n6_d7.json` | algebra in 7 dimensions, the same six-dimensional points; sections A, D, F | the g → gg identities hold |
| `dependence_n6.json` | comparison of the two | every one-gluon Wilson-line term and 24 of the 32 two-gluon terms change with the dimension |
| `n5_d5.json`, `n5_d6.json`, `n5_d7.json` | algebra in 5, 6, 7 dimensions, one five-dimensional point; sections A, D | the identity holds in each |
| `dependence_n5.json` | comparison of the three | the nine Wilson-line terms are linear in the dimension |

What follows. For fixed vectors a g → gg numerator is a polynomial of degree one in the dimension (one
closed loop of Lorentz indices), so the identity at two dimensions and the same points holds for every
dimension; with six-dimensional vectors the transverse vectors of a point obey no linear relation, so it
holds as an identity in their scalar products. The g → q q̄ numerators have no closed index loop; with open
spinor indices the identity is one of the Clifford algebra of the vectors, for which the 8 × 8 matrices of
six dimensions are a faithful representation. γ5 is that of 't Hooft and Veltman.

Not covered, as in four dimensions: ultrasoft gluons that interact with each other or carry a
self-energy, and the regulated loop integrals themselves.

Python 3.13, SymPy 1.13.3. A mutation in which one vertex of the prediction loses its components in the
additional dimensions passes in four dimensions and fails in six.
