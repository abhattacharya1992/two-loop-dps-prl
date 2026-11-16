# Independent five-loop NMHV comparison

The selected component proportional to `pi i (1+w)/(1+w*)`, in the
self-crossing limit's multi-Regge overlap, agrees through five loops with an
independent expansion of `hexMRKL1-7.m`, the ancillary data to arXiv:1903.10890v2.
All nine five-loop coefficients, including the constant and odd-zeta terms,
agree exactly. The construction also reproduces all 16 coefficients through
four loops from the Dixon--Esterlis erratum. The one-loop entry fixes the
normalization, so the higher-loop comparisons comprise 24 coefficients.

## Reproduction

Run from this directory:

```bash
make primary-inputs
make primary-check
```

The downloader verifies the archive and extracted file against
`inputs/SHA256SUMS`. The calculation needs Python and SymPy; it does not need
Mathematica or a numerical special-function package. The result is
`results_primary_nmhv.json`. The primary amplitude data are downloaded rather
than redistributed.

## Endpoint construction

Use the Brown single-valued HPL construction described in section 3.3 of
arXiv:1207.0186. Let `Z` denote the shuffle-regularized Drinfeld associator and
let a tilde reverse its noncommutative words. Determine the second alphabet by

```
y0 = x0,
tilde(Z)(y) y1 tilde(Z)(y)^(-1) = Z(x)^(-1) x1 Z(x).
```

The script solves this equation by weight through eight, and verifies its
residual and the independently published weight-four terms of `y1`.
Multiple-zeta constants through weight eight are derived using shuffle,
stuffle, duality and the single-zeta(1) regularization relation. The
weight-eight basis constant zeta(5,3) is kept symbolic; it cancels from the
reported answer. No constants from the fitted DPS kernel enter these steps.

Set `t=1-z`, `tbar=1-zbar`, `u=log(t*tbar)` and
`A=Z(x) tilde(Z)(y)`. The endpoint constant series is

```
L0 = exp(-u x1) A.
```

Write the single-valued generating series as `L0 + t L10 + tbar L01 + ...`.
Its holomorphic differential equation implies

```
(1 + d/du) L10 = -x1 L10 - x0 L0.
```

It determines the polynomial coefficient of every word recursively in its
weight. On a polynomial, `(1+d/du)^(-1)` is the finite sum of successive
derivatives with alternating signs. Thus the first endpoint term at weight
nine needs endpoint constants only through weight eight; setting odd zetas
to zero is unnecessary.

For each primary `FMtres[L,0]/pi`, retain its coefficient of `t` with this
expansion. The NMHV rational prefactors and their target-projectile partner
select the phase `(1+w)/(1+w*)`. In the ancillary conventions the coefficient
polynomial is minus this holomorphic endpoint coefficient. Convert to the
2016 BDS-like convention with `exp(gamma_K*zeta(2)/2)` after extracting the
polynomial independent of `log|v|`, using
`a=g^2 N_c/(8 pi^2)` and

```
gamma_K = 4a - 4 zeta(2) a^2 + 22 zeta(4) a^3
          - (219 zeta(6)/2 + 4 zeta(3)^2) a^4.
```

The full lower-loop coefficients provide normalization and branch controls
before the five-loop coefficients are compared. No coefficient is fitted to
the five-loop target, and a missing target coefficient fails the command.

## Scope

This tests the displayed five-loop component in the overlap limit. It does
not test terms proportional to log|v|, determine the general NMHV
self-crossing amplitude away from that overlap, or prove the common-kernel
relation to all orders. It places no all-order constraint on QCD boundary
logarithms.
