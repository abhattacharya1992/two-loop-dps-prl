"""The one-loop BDS-like function of the six-gluon amplitude, eq. (7.6) of Paper 2, from the box functions.

Planar N = 4, MHV, BDS normalisation (a = g^2 N_c/(8 pi^2) (4 pi e^-gamma)^eps): the one-loop amplitude is
M_6^(1) = (1/2) V_6 (1 - zeta_2 eps^2/2 + ...), V_6 the sum of the six one-mass and three two-mass-easy box functions of
Bern, Dixon, Dunbar and Kosower (c_Gamma normalisation; the factor converts c_Gamma e^{gamma eps} to 1):
    F^1m(s,t;P^2)      = -(1/eps^2)[(-s)^-eps + (-t)^-eps - (-P^2)^-eps] + Li2(1-P^2/s) + Li2(1-P^2/t) + ln^2(s/t)/2 + zeta_2,
    F^2me(s,t;P^2,Q^2) = -(1/eps^2)[(-s)^-eps + (-t)^-eps - (-P^2)^-eps - (-Q^2)^-eps] + Li2(1-P^2/s) + Li2(1-P^2/t)
                         + Li2(1-Q^2/s) + Li2(1-Q^2/t) - Li2(1-P^2 Q^2/(s t)) + ln^2(s/t)/2,
and M6hat = 2 M_6^(1) + Y with Y = sum_i [Li2(1-u_i) + ln^2(u_i)/2] (Dixon-Esterlis eq. (2.10)).  Checked here, at random
Euclidean points (all invariants negative, S = -s > 0) to 25 digits:
  (1) M6hat does not depend on the three-particle invariants;
  (2) M6hat = sum_i [-(1/eps^2)(mu^2/S_i)^eps] + Q(l) + 9 zeta_2 through eps^0, with l_i = ln(S_i/mu^2) and
      Q = (1/2) sum l_i^2 - sum l_i l_{i+1} + (1/2) sum_{i=1}^6 l_i l_{i+3}, which is invariant under l_i -> l_i + c;
  (3) the four-point analogue: F^0m = -(1/eps^2)[(-s)^-eps + (-t)^-eps] + ln^2(s/t)/2 + 3 zeta_2 and the same conversion
      give M_4^(1) = ... + ln^2(s/t)/2 + 4 zeta_2, eq. (7.5), and F^1m -> F^0m as P^2 -> 0, F^2me -> F^1m as Q^2 -> 0.

    python3 check_bdslike.py
"""
import random
import sympy as sp

S12, S23, S34, S45, S56, S61, S123, S234, S345 = sp.symbols('S12 S23 S34 S45 S56 S61 S123 S234 S345', positive=True)
z2 = sp.zeta(2)
ln = sp.log
Li2 = lambda x: sp.polylog(2, x)
two = [S12, S23, S34, S45, S56, S61]
three = [S123, S234, S345, S123, S234, S345]           # s_{i,i+1,i+2}


def F1m_fin(s, t, P2):
    return Li2(1 - P2/s) + Li2(1 - P2/t) + ln(s/t)**2/2 + z2


def F2me_fin(s, t, P2, Q2):
    return Li2(1 - P2/s) + Li2(1 - P2/t) + Li2(1 - Q2/s) + Li2(1 - Q2/t) - Li2(1 - P2*Q2/(s*t)) + ln(s/t)**2/2


V_fin = sum(F1m_fin(two[i], two[(i + 1) % 6], three[i]) for i in range(6))
V_fin += F2me_fin(S123, S345, S12, S45) + F2me_fin(S123, S234, S23, S56) + F2me_fin(S234, S345, S34, S61)
# poles of V: -(1/eps^2) sum_i S_i^-eps (each two-particle invariant once; the three-particle ones cancel between the
# one-mass and two-mass-easy boxes); their eps^0 part is -(1/2) sum ln^2 S_i; the conversion adds (3/2) zeta_2 to M_6
M6_eps0 = sp.Rational(1, 2)*(V_fin - sum(ln(s)**2/2 for s in two)) + sp.Rational(3, 2)*z2
u = S12*S45/(S123*S345); v = S23*S56/(S234*S123); w = S34*S61/(S345*S234)
Y = Li2(1 - u) + Li2(1 - v) + Li2(1 - w) + (ln(u)**2 + ln(v)**2 + ln(w)**2)/2
M6hat_eps0 = 2*M6_eps0 + Y
l = [ln(s) for s in two]                                 # mu = 1
Q = sum(x**2 for x in l)/2 - sum(l[i]*l[(i + 1) % 6] for i in range(6)) + sum(l[i]*l[(i + 3) % 6] for i in range(6))/2
form_eps0 = -sum(x**2 for x in l)/2 + Q + 9*z2            # eps^0 of eq. (7.6)

random.seed(7)
worst_three = worst_form = 0
for trial in range(6):
    vals = {s: sp.Rational(random.randint(2, 60), random.randint(1, 13)) for s in two}
    a = {**vals, S123: sp.Rational(7, 3), S234: sp.Rational(11, 5), S345: sp.Rational(5, 2)}
    b = {**vals, S123: sp.Rational(13, 4), S234: sp.Rational(3, 7), S345: sp.Rational(9, 11)}
    ma, mb = sp.N(M6hat_eps0.subs(a), 30), sp.N(M6hat_eps0.subs(b), 30)
    worst_three = max(worst_three, abs(ma - mb))
    worst_form = max(worst_form, abs(ma - sp.N(form_eps0.subs(a), 30)))
assert worst_three < 1e-25 and worst_form < 1e-25, (worst_three, worst_form)
c = sp.Symbol('c')
assert sp.expand(Q.subs({x: x + c for x in l}, simultaneous=True) - Q) == 0        # scale invariance of Q
# the four-point analogue and the massless limits of the box functions
s_, t_, P2, Q2 = sp.symbols('s t P2 Q2', positive=True)
F0m_fin = ln(s_/t_)**2/2 + 3*z2
assert sp.simplify(sp.limit(F1m_fin(s_, t_, P2), P2, 0) - F0m_fin) == 0
assert sp.simplify(sp.limit(F2me_fin(s_, t_, P2, Q2), Q2, 0) - F1m_fin(s_, t_, P2)) == 0
M4_const = 3*z2 + z2            # (1/2) V_4 = F^0m, double pole -2/eps^2, conversion +zeta_2
assert M4_const == 4*z2
print(f'M6hat independent of the three-particle invariants to {float(worst_three):.0e}; '
      f'equal to the poles + Q + 9 zeta_2 to {float(worst_form):.0e}; Q scale invariant; '
      f'F^1m -> F^0m, F^2me -> F^1m; M_4^(1) constant 4 zeta_2')
print('passed')
