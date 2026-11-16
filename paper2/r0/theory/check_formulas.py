"""The analytic steps behind Sec. 7 and appendix F of Paper 2, one by one, as the re-derivation workbook
(Part VII, and step VI.5) states them.  Each check prints one line and the script stops at the first failure.

    python3 check_formulas.py        (about 10 s; sympy and mpmath; 12 checks)

 1  the bracket of eq. (5.3) at a general scale: X_mu = X + L_*, the X-independent part Dc L_*^2 + B L_*
 2  the transforms of eq. (5.5): finite parts for nu = 0 (no G^1 from the X and X^2 insertions), the
    nu = 2 constant B - 2 Dc of -Dc X^2 + B X
 3  quark splittings in the gluon variable: X = X_g + kappa, L_* = L_g - kappa, Dc kappa = 2 (C_A - C_F) L_x
 4  [2 I^(1)]_{eps^0} at mu^2 = s and the conversions Catani -> MS, qT -> MS
 5  the constants of eq. (7.2): 2 gamma_cusp^(1)/gamma_cusp^(0) - C_A zeta_2
 6  F_eps of eq. (F.2) from the D-dimensional Fourier transform with the MS-bar factor (mu^2 e^gamma)^eps
 7  the soft function, eq. (F.3): the transform of [2/eta + ln(nu^2/l^2)]/l^2; no eps^0 in its 1/eta term
    and in the transform of 1/l^2 (step VI.5)
 8  the T_F n_f single pole: +4/3 per collinear function, -4/3 per daughter self-energy, total -8/3
 9  the transform of the cubic structure: its spin-3 part changes sign, T4 -> (S1 + S2 + S3)/2 - T4
10  the helicity projections of S1, S2, S3, T4, of the tree and of R = (c/6)(S1 - S2 - S3), and the T4
    coefficient of J1, 2c/3, against the rational terms of the one-loop splitting amplitudes
11  the poles of 2 I_K: real parts plus beta0/eps are those of the regions; the imaginary pole 2 pi i Dc/eps
12  the hard functions of gg -> H and gg -> A from the renormalised one-loop form factors (d = 4 + eps,
    finite parts C_A zeta_2 and C_A(4 + zeta_2) at Q^2 = mu^2): restoring (Q^2/mu^2)^{eps/2} and removing the
    poles gives h = C_A(-L^2 + c + zeta_2) with no beta0 L term, c = 0 (scalar), 4 (pseudoscalar)
"""
import json
import numpy as np
import sympy as sp
import mpmath as mp

eps, al, de, eta = sp.symbols('epsilon alpha delta eta')
G, X, Xg, L, Lg, Lx, kap = sp.symbols('G X X_g L_* L_g L_x kappa', real=True)
Dc, B, CA, CF, TF, nf, b0s = sp.symbols('Delta_C B C_A C_F T_F n_f beta0', real=True)
z2, z3 = sp.zeta(2), sp.zeta(3)
ok = lambda e: sp.simplify(sp.expand(e)) == 0


def report(n, text):
    print(f'{n:2d}  {text}')


# 1 ------------------------------------------------------------------------------------------------------------
Xmu = X + L
br = Dc*Xmu**2 - 2*Dc*X*Xmu + B*Xmu
assert ok(br - (-Dc*X**2 + B*X + Dc*L**2 + B*L))
report(1, 'Dc X_mu^2 - 2 Dc X X_mu + B X_mu = -Dc X^2 + B X + Dc L_*^2 + B L_*')

# 2 ------------------------------------------------------------------------------------------------------------
def phi(nu):
    return sp.exp(2*(sp.EulerGamma - G)*al)*sp.gamma(sp.Rational(nu, 2) + al)/sp.gamma(1 + sp.Rational(nu, 2) - al)/2


s0 = sp.series(phi(0), al, 0, 4).removeO()
fp = [sp.nsimplify(sp.expand(sp.factorial(n)*s0.coeff(al, n))) for n in range(3)]
fp = [sp.simplify(f) for f in fp]
assert ok(fp[0] + G) and ok(fp[1] - G**2) and ok(fp[2] - (-sp.Rational(4, 3)*G**3 - sp.Rational(2, 3)*z3))
assert all(sp.expand(f).coeff(G, 1) == 0 for f in fp[1:])
s2 = sp.series(phi(2), al, 0, 3).removeO()
r1 = sp.simplify(s2.coeff(al, 1)/s2.coeff(al, 0))
r2 = sp.simplify(2*s2.coeff(al, 2)/s2.coeff(al, 0))
assert ok(r1 - (1 - 2*G)) and ok(r2 - (4*G**2 - 4*G + 2))
assert ok(sp.expand(-Dc*r2 + B*r1).subs(G, 0) - (B - 2*Dc))
report(2, 'FP Phi_0, Phi_0\', Phi_0\'\' = -G, G^2, -4/3 G^3 - 2/3 zeta3 (no G^1 from X, X^2); '
          'Phi_2\'/Phi_2 = 1 - 2G, Phi_2\'\'/Phi_2 = 4G^2 - 4G + 2; nu = 2 constant B - 2 Dc')

# 3 ------------------------------------------------------------------------------------------------------------
e3 = (-Dc*X**2 + B*X + Dc*L**2 + B*L).subs({X: Xg + kap, L: Lg - kap}, simultaneous=True)
assert ok(e3 - (-Dc*Xg**2 + (B - 2*Dc*kap)*Xg + Dc*Lg**2 - 2*Dc*kap*Lg + B*Lg))
lnm2, lns = sp.symbols('lnm2 lns')
dcq = 2*CF - CA
lnQ = (CF*lnm2 + CF*lnm2 - CA*lns)/dcq                     # Delta_C ln Q_*^2 for q qbar -> Z, Z
lnQg = 2*lnm2 - lns                                          # m^4/s
kappa_q = sp.simplify((lnQ - lnQg).subs(lnm2, Lx + lns))     # m^2 = x(1-x) s
assert ok(dcq*kappa_q - 2*(CA - CF)*Lx)
report(3, 'gluon variable: -Dc X_g^2 + (B - 2 Dc kappa) X_g + Dc L_g^2 - 2 Dc kappa L_g + B L_g, '
          'Dc kappa = 2 (C_A - C_F) L_x')

# 4 ------------------------------------------------------------------------------------------------------------
def catani(C, g, Q2_over_mu2=1):
    """I^(1) of eq. (4.1), (-mu^2/Q^2)^eps = e^{i pi eps} (mu^2/Q^2)^eps"""
    return -sp.exp(eps*sp.EulerGamma)/sp.gamma(1 - eps)*(C/eps**2 + g/eps)*sp.exp(sp.I*sp.pi*eps)*Q2_over_mu2**(-eps)


twoI0 = sp.simplify(sp.series(2*catani(CA, b0s/2), eps, 0, 1).removeO().coeff(eps, 0))
assert ok(twoI0 - (7*z2*CA - sp.I*sp.pi*b0s))
# ggvvamp: c^Catani = c^qT + (-pi^2 C_A + i pi beta0) c1 (a_s units); Catani = MS - [2I]_0 c1
assert ok((-sp.pi**2*CA + sp.I*sp.pi*b0s) + twoI0 - CA*z2)
report(4, '[2 I^(1)]_{eps^0} = 7 zeta2 C_A - i pi beta0 at mu^2 = s; remainders = MS - [2I]_0 c1 (Catani), '
          'MS - C_A zeta2 c1 (qT)')

# 5 ------------------------------------------------------------------------------------------------------------
CR = sp.Symbol('C_R')
g0, g1 = 4*CR, 4*CR*((sp.Rational(67, 9) - 2*z2)*CA - sp.Rational(20, 9)*TF*nf)
assert ok(2*g1/g0 - CA*z2 - (CA*(sp.Rational(134, 9) - 5*z2) - sp.Rational(40, 9)*TF*nf))
report(5, 'C_A(134/9 - 5 zeta2) - (40/9) T_F n_f = 2 gamma^(1)/gamma^(0) - C_A zeta2')

# 6 ------------------------------------------------------------------------------------------------------------
# int d^Dk e^{ik.b} k^i (k^2)^-a  ~  b^i Gamma(D/2 + 1 - a)/(4^a Gamma(a)) (b^2)^(a - D/2 - 1),  D = 2 - 2 eps
D = 2 - 2*eps
ft = lambda a: sp.gamma(D/2 + 1 - a)/(4**a*sp.gamma(a))
mub2 = sp.Symbol('mub2')                                    # mu^2 b^2
ratio = ft(1 + eps)/ft(1)*mub2**eps*sp.exp(eps*sp.EulerGamma)
Feps = sp.exp(-eps*sp.EulerGamma)*sp.gamma(1 - 2*eps)/(sp.gamma(1 + eps)*sp.gamma(1 - eps))
val = sp.simplify(ratio.subs(mub2, 4*sp.exp(-2*sp.EulerGamma)))   # mu b = b0 = 2 e^{-gamma}
for ev in (sp.Rational(1, 7), -sp.Rational(2, 9), sp.Rational(3, 10)):
    assert abs(sp.N((val - Feps).subs(eps, ev), 30)) < 1e-25
report(6, 'F_eps = e^{-eps gamma} Gamma(1-2eps)/(Gamma(1+eps) Gamma(1-eps)): transform of k^i (mu^2 e^gamma/k^2)^eps '
          'relative to k^i, at mu b = b0')

# 7 ------------------------------------------------------------------------------------------------------------
# rapidity integral 2/eta + ln(nu^2/l^2); transform of (l^2)^(-1-delta): Gamma(-eps-delta)/Gamma(1+delta) (b^2/4)^(eps+delta)
mp.mp.dps = 30
for ev in (mp.mpf('0.13'), mp.mpf('-0.21'), mp.mpf('0.37')):
    b = mp.mpf(1); nu = 2*mp.exp(-mp.euler)/b; mu = nu
    f = lambda d: mp.gamma(-ev - d)/mp.gamma(1 + d)*(b**2/4)**(ev + d)
    lhs = (mp.log(nu**2)*f(0) + mp.diff(f, 0))*mp.exp(ev*mp.euler)*mu**(2*ev)
    rhs = mp.exp(-ev*mp.euler)*mp.gamma(-ev)*(mp.log(nu**2*b**2/(4*mp.exp(-2*mp.euler))) - mp.euler - mp.digamma(-ev))
    assert abs(lhs - rhs) < mp.mpf(10)**-20, (ev, lhs, rhs)
ser = sp.series(sp.exp(-eps*sp.EulerGamma)*sp.gamma(-eps), eps, 0, 2).removeO()
assert ok(ser.coeff(eps, -1) + 1) and ok(ser.coeff(eps, 0)) and ok(ser.coeff(eps, 1) + z2/2)
report(7, 'S = 2 Dc e^{-eps gamma} Gamma(-eps)[2/eta + ln(nu^2 b^2/b0^2) - gamma - psi(-eps)]; '
          'e^{-eps gamma} Gamma(-eps) = -1/eps - zeta2 eps/2 + ..., no eps^0 term')

# 8 ------------------------------------------------------------------------------------------------------------
se = (-4*(1 - eps)/(3 - 2*eps))*sp.gamma(eps)*sp.gamma(1 - eps)**2/sp.gamma(2 - 2*eps)
assert ok(sp.series(se, eps, 0, 1).removeO().coeff(eps, -1) + sp.Rational(4, 3))
import combine_gg as CG
worst = 0
for tag in ('1o10', '3o10', '1o2', '7o10', '9o10'):
    J, tau = CG.load(tag, 'nf')
    worst = max(worst, abs(J['S1'].get('0,-1', 0) - sp.Rational(4, 3)), abs(J['S1'].get('0,-2', 0)))
assert worst < 1e-8, worst
report(8, f'T_F n_f: 1/eps of J1 = +4/3 times the tree (to {worst:.0e}), of each daughter self-energy -4/3; '
          'kernel 2(4/3) - 4(4/3) = -8/3')

# 9 ------------------------------------------------------------------------------------------------------------
kx, ky = sp.symbols('k_x k_y', real=True)
k = sp.Matrix([kx, ky]); k2 = kx**2 + ky**2
dlt = sp.eye(2)
T = lambda a, i, j: k[a]*k[i]*k[j]/k2
spin3 = lambda a, i, j: T(a, i, j) - (dlt[i, j]*k[a] + dlt[a, i]*k[j] + dlt[a, j]*k[i])/4
for a in range(2):
    assert sp.simplify(sum(spin3(a, i, i) for i in range(2))) == 0 and sp.simplify(sum(spin3(i, i, a) for i in range(2))) == 0
# 2D: int d^2k/(2pi)^2 e^{ik.b} e^{i m phi} g(k) = i^m e^{i m phi_b}/(2pi) int k dk J_m(kb) g(k); int_0^oo J_m = 1
assert abs(mp.quadosc(lambda t: mp.besselj(1, t), [0, mp.inf], omega=1) - 1) < 1e-10
assert abs(mp.quadosc(lambda t: mp.besselj(3, t), [0, mp.inf], omega=1) - 1) < 1e-10
assert sp.I**1 == -sp.I**3
report(9, 'k-space T4 = spin 3 + (S1+S2+S3)/4; spin 1 -> i^1, spin 3 -> i^3 = -i^1: '
          'T4 -> (S1 + S2 + S3)/2 - T4 in b space (eps^0)')

# 10 -----------------------------------------------------------------------------------------------------------
ep, em = np.array([1, 1j])/np.sqrt(2), np.array([1, -1j])/np.sqrt(2)
rng = np.random.default_rng(3)
kv = rng.normal(size=2)
ksq = kv @ kv
def structs(e, i, j):
    return {'S1': (e @ kv)*(i @ j), 'S2': (e @ i)*(j @ kv), 'S3': (e @ j)*(i @ kv),
            'T4': (e @ kv)*(i @ kv)*(j @ kv)/ksq}
for xv in (0.1, 0.3, 0.5):
    for c in (1, -2):
        rat = {}
        for name, (i, j) in {'++': (em, em), '+-': (em, ep), '-+': (ep, em), '--': (ep, ep)}.items():
            s = structs(ep, i, j)                           # incoming +; outgoing daughter +: eps_+^* = eps_-
            tree = s['S1'] - s['S2']/(1 - xv) - s['S3']/xv
            R = c/6*(s['S1'] - s['S2'] - s['S3'])
            rat[name] = (tree, R, s['T4'])
        tr, R, _ = rat['++']
        assert abs(R/tr - c*xv*(1 - xv)/3) < 1e-12                        # BDKS (4.33) at eps = 0: 2z(1-z)/6
        assert abs(rat['+-'][1]) < 1e-12 and abs(rat['-+'][1]) < 1e-12
        assert abs(rat['--'][0]) < 1e-12 and abs(rat['--'][2]) > 1e-3     # the flip: no tree, T4 only
worstT4 = 0
for tag in ('1o10', '3o10', '1o2', '7o10', '9o10'):
    for col, c in (('CA', 1), ('nf', -2)):
        J, tau = CG.load(tag, col)
        worstT4 = max(worstT4, abs(J['T4'].get('0,0', 0) - 2*c/3), abs(J['T4'].get('0,-1', 0)))
assert worstT4 < 1e-6, worstT4
report(10, f'R/tree = x(1-x)c/3 in (+;++), 0 in (+;+-), (+;-+); the flip (+;--) has T4 only; '
           f'T4 of J1 = 2c/3, finite, at the five x (to {worstT4:.0e})')

# 11 -----------------------------------------------------------------------------------------------------------
def poles(e):
    s = sp.series(sp.expand(e), eps, 0, 0).removeO()
    return [sp.simplify(s.coeff(eps, -2)), sp.simplify(s.coeff(eps, -1))]


m2 = sp.Symbol('m2', positive=True)                          # m_1^2 = m_2^2 = x(1-x) in units s = mu = 1
beta0 = sp.Rational(11, 3)*CA - sp.Rational(4, 3)*TF*nf
gl = 2*catani(CA, beta0/2) - 2*2*catani(CA, beta0/2, m2)
qu = 2*catani(CA, beta0/2) - 2*2*catani(CF, sp.Rational(3, 2)*CF, m2)
pg, pq = poles(gl), poles(qu)
pg = [sp.expand(p.subs(sp.log(m2), Lx)) for p in pg]
pq = [sp.expand(p.subs(sp.log(m2), Lx)) for p in pq]
assert ok(pg[0] - 2*CA) and ok(sp.re(pg[1]) + beta0 - (CA*(sp.Rational(22, 3) - 4*Lx) - sp.Rational(8, 3)*TF*nf))
assert ok(sp.im(pg[1]) - 2*sp.pi*CA)
assert ok(pq[0] - (4*CF - 2*CA)) and ok(sp.re(pq[1]) + beta0 - CF*(6 - 4*Lx)) and ok(sp.im(pq[1]) - 2*sp.pi*(2*CF - CA))
report(11, 'poles of 2 I_K + beta0/eps: C_A[2/eps^2 + (22/3 - 4L_x)/eps] - (8/3)T_F n_f/eps (gluons), '
           'C_F[4/eps^2 + (6 - 4L_x)/eps] - 2C_A/eps^2 (quarks); imaginary pole 2 pi i Dc/eps')

# 12 -----------------------------------------------------------------------------------------------------------
# Renormalised one-loop form factors of G G (scalar) and G Gtilde (pseudoscalar) in d = 4 + eps', at mu^2 = Q^2:
# C_A{-8/eps'^2 + 22/(3 eps') + c + zeta_2} - 8 T_F n_f/(3 eps'), c = 0 and 4 (Ahmed et al. 2015, eq. 3.51 for c = 4).
# The pole 2 beta0/eps' is the operator renormalisation Z_GG = Z_as and does not carry the scale factor; the rest
# does: a_s(mu) [(Q^2/mu^2)^{eps'/2} F_bare + 2 beta0/eps'] with Q^2 = -m^2 - i0, (Q^2/mu^2)^{eps'/2} = e^{eps' L/2},
# L = ln(m^2/mu^2) - i pi.  The eps'^0 term is the minimally subtracted hard function h.
epsp, Lh, c = sp.symbols('epsilonp L c')
ff = sp.exp(epsp*Lh/2)*CA*(-8/epsp**2 + c + z2) + 2*beta0/epsp
ser = sp.series(ff, epsp, 0, 1).removeO()
assert ok(ser.coeff(epsp, 0) - CA*(-Lh**2 + c + z2))                      # h, with no beta0 L term
assert ok(ser.coeff(epsp, -2) + 8*CA) and ok(ser.coeff(epsp, -1) - (2*beta0 - 4*CA*Lh))
assert ok(ser.coeff(epsp, 0).subs(c, 4) - ser.coeff(epsp, 0).subs(c, 0) - 4*CA)
report(12, 'h = C_A(-L^2 + c + zeta_2) from the renormalised form factors, c = 0 (gg -> H), 4 (gg -> A); '
           'poles -8C_A/eps^2 + (2 beta0 - 4 C_A L)/eps; the two hard functions differ by 4 C_A')
# 13 -----------------------------------------------------------------------------------------------------------
# Planar N = 4 (Sec. 7.3): (a) the field content in the daughter self-energies, four-dimensional-helicity scheme:
# gluon + ghost (5-3eps)/(3-2eps), four adjoint Weyl fermions = 2 C_A in T_F n_f units, -8(1-eps)/(3-2eps), and
# (3 + eps) complex adjoint scalars (the 2 eps epsilon-scalars included), -1/(3-2eps) each: the sum is -2 exactly.
# (b) The two-loop singular terms of E = A_6/A_6^{BDS-like} in units of the one-loop DPS term, h_1 + h_2 - M6hat/2
# + K^(1)/K_0 in 3 -> 3 kinematics: the poles cancel and the finite part is -X^2/2 + ln x ln(1-x) + ln y ln(1-y)
# + 7 zeta_2/2, independent of the scattering angles (n4_hexagon.symbolic).
import n4_diagrams, n4_hexagon
cs = n4_diagrams.self_energy_coefficient()
e_ = sp.Symbol('eps')
assert sp.simplify((5 - 3*e_)/(3 - 2*e_) - 8*(1 - e_)/(3 - 2*e_) + (3 + e_)*cs + 2) == 0
sym = n4_hexagon.symbolic()
assert sym['poles'] == 'cancel' and sym['X^2'] == '-1/2' and sym['X^1'] == '0' and not sym['depends on t1, t2']
report(13, 'N = 4, FDH: daughter self-energies (5-3eps)/(3-2eps) - 8(1-eps)/(3-2eps) - (3+eps)/(3-2eps) = -2; '
           'hexagon: h_1 + h_2 - M6hat/2 + kernel poles finite, = -X^2/2 + ln x ln(1-x) + ln y ln(1-y) + 7 zeta_2/2')
print('passed')
