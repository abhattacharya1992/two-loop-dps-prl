#!/usr/bin/env python3
"""Symbolic threshold consistency checks for the new two-loop DPS PRL.

The published threshold inputs are transcribed independently from
A. Bhattacharya, arXiv:2610.19024v2, Eqs. (2.12), (4.5)--(4.7).

The predictions are the x=1/2 limits of the new PRL Letter's
one-loop weights, two-loop evolution, and finite r0 for the AA parity.
The single-log comparison includes the AA hard form factor and the
finite Catani-to-minimal-subtraction conversion of the PRL Supplement.

These are symbolic consistency checks of two kinematic analyses of
*the same underlying two-loop amplitude*. They do not regenerate its
full original expressions or prove factorisation at arbitrary x.
"""

import sympy as sp

N, nf = sp.symbols('N_c n_f', integer=True, positive=True)
pi, I, L2 = sp.pi, sp.I, sp.log(2)
zeta2 = pi**2 / 6
beta0 = (11*N - 2*nf) / 3

# All dimensional coefficients are divided by m_A^2.
# At threshold: x=1/2, m_A=m_H=m, mu_R=m, s=4m^2.
x = sp.Rational(1, 2)
sigma = x*(1-x)
F_A = 1/sigma - 2*sigma
F_S = 1/sigma + 2*sigma
W = (x**4 + (1-x)**4)/sigma

# Independent reference: P1 threshold expansion M_i^(2)/m^2 = sum A_ik L^k.
# These are the explicit coefficients in Eqs. (4.5),(4.6) of arXiv:2610.19024v2.
old = {
    'C1': 7*pi*N,                        # P1 Eq. (4.7), M_1^(1)/m^2
    'K': pi*N/2,                         # P1 Eq. (4.7), M_2^(1)/m^2
    'A13': -sp.Rational(28, 3)*pi*N**2,
    'A12': -sp.Rational(14, 3)*pi*N*(11*N-2*nf),
    'A22': -2*pi*N**2,
    'A21': -sp.Rational(4, 3)*pi*N*(4*N-nf),
}
old['A11'] = (pi*N/sp.Integer(9)) * (
    N*(1454 + 462*L2 + 252*L2**2 - 42*pi**2 + 693*I*pi)
    - nf*(152 + 84*L2 + 126*I*pi)
)  # P1 Eq. (4.6) at L_mu = ln(m^2/mu_R^2) = 0.

# PRL's C1 and K: the AA tree has M_1^(0)/m^2=2i.
C1 = -4*pi*I*N*sigma*F_A*(2*I)
K = 2*pi*N*sigma*W*2  # positive with the P1 helicity convention
B = 2*beta0
pred = {
    'C1': C1,
    'K': K,
    'A13': C1*(-sp.Rational(4, 3)*N),
    'A12': C1*(-B),
    'A22': K*(-4*N),
    'A21': K*(4*N - 2*B),
}

# At the threshold Q_*^2=m^4/s=m^2/4, and G = ln(2 p_T/m)+i*pi/2,
# equal to the P1 leading-power L variable. The r0 term of the PRL
# for the AA equal-helicity projection, evaluated at x=1/2:
r0_AA = (
    N*(sp.Rational(134, 9) - 5*zeta2 - 4*L2**2)
    - sp.Rational(40, 9)*sp.Rational(1, 2)*nf
    + 2*(N-nf)/(3*F_A)
)

# AA one-loop hard-subprocess correction, m_A = mu_R = m.
# From the PRL Supplement: h_AA = C_A[-L_m^2 + 4 + zeta(2)], L_m=-i*pi.
Lm = -I*pi
h_AA = N*(-Lm**2 + 4 + zeta2)
Lstar = 2*L2 + I*pi      # ln[mu^2 / (Q_*^2 e^{-i*pi})]
Ls = -2*L2 + I*pi        # ln(mu^2 / s) + i*pi

# Minimal-subtraction c11 from the PRL, followed by the finite
# Catani subtraction used for the old AA remainder.
finite_2I = -N*Ls**2 + N*zeta2 - beta0*Ls
c11_MS = 2*h_AA + N*Lstar**2 + B*Lstar + r0_AA
pred['A11'] = C1*(c11_MS - finite_2I)

# A separate check of the published CP-even/odd one-loop weight ratio.
assert sp.simplify(F_A-sp.Rational(7, 2)) == 0
assert sp.simplify(F_S-sp.Rational(9, 2)) == 0
assert sp.simplify(W-sp.Rational(1, 2)) == 0
assert sp.simplify(F_S/F_A-sp.Rational(9, 7)) == 0

for name in ('C1', 'K', 'A13', 'A12', 'A22', 'A21', 'A11'):
    delta = sp.expand(pred[name] - old[name]).simplify()
    if delta != 0:
        raise AssertionError(f'{name}: the two calculations do not agree: {delta}')
    print(f'PASS {name}: exact symbolic equality for arbitrary N_c and n_f')

# Negative controls: a wrong gluon-cusp sign or forgetting the
# Catani finite conversion must not pass the benchmark.
assert sp.simplify(pred['A13'] + old['A13']) != 0
assert sp.simplify(C1*c11_MS - old['A11']) != 0

print('PASS scalar/pseudoscalar one-loop weight ratio F_S(1/2)/F_A(1/2) = 9/7')
print('PASS negative controls: wrong sign or missing scheme conversion rejected')
print('SCOPE: threshold beta_A -> 0, x=1/2 and mu_R=m only;')
print('       compares published coefficients, not new evaluations of raw amplitudes.')
