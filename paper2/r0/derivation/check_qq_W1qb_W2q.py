"""Wilson line couplings in J1(g -> q qbar) with colour factor C_F - C_A/2: diagrams W2q and W1qb.
The gluon from the Wilson line of one operator attaches to the daughter created by the other operator.
These two diagrams carry the rapidity divergence of the quark channel.

Units k^2 = mu^2 = nu_r = 1, a_s = alpha_s/(4 pi). Each diagram is (C_F - C_A/2) [A T + (B/2) (e.k)] with the tree
splitting T = (x/2) eslash kslash - ((1-x)/2) kslash eslash. For W2q, with xi = z/(1-x) the plus momentum of the
exchanged gluon in units of that of the antiquark, the residues and the transverse integral leave
    A =  2 F' (1-x)^(-eta) Int_0^1 dxi xi^(-1-eta) (1-xi)^(1-eps) I(xi),
    B = -4 (1-x) F'        Int_0^1 dxi xi^(-eta)   (1-xi)^(1-eps) I(xi),
    I(xi) = Int_0^1 dt t^(-eps) [(1-t)(1-xi) + xi/x]^(-1-eps),      F' = e^(-eps gamma_E) Gamma(1-2eps)/Gamma(1-eps),
in impact parameter space at b = b0. Expanding in eta and then in eps, with L = ln(x(1-x)),
    W2q:   A = 2/(eta eps) - 2/eps^2 - 2L/eps + 2 Li2(1-x) - 3 zeta2,    B = 4 x ln x,
    W1qb:  A with x <-> 1-x,                                             B = -4 (1-x) ln(1-x)   (charge conjugation).

Steps checked here:
  1  B at eps = eta = 0
  2  the coefficient of eta^0 in A: the subtracted integral at finite eps, its Laurent expansion, and the assembly of A
  3  the results against the stored Laurent coefficients
Usage: python3 check_qq_W1qb_W2q.py <repository root>"""
import json, sys, pathlib
import numpy as np
from scipy.special import gamma, spence, beta
from scipy.integrate import quad
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from report import check, require, finish
repo = pathlib.Path(sys.argv[1])
Z2 = np.pi**2/6
Li2 = lambda t: spence(1 - t)
Q = dict(epsabs=1e-13, epsrel=1e-13, limit=200)


def I(xi, x, e):
    return quad(lambda t: t**(-e)*((1 - t)*(1 - xi) + xi/x)**(-1 - e), 0, 1, **Q)[0]


I_at_0 = lambda e: gamma(1 - e)*gamma(-e)/gamma(1 - 2*e)                 # I(0)

print('1  B at eps = eta = 0: -4 (1-x) Int_0^1 dxi (1-xi) I(xi) = 4 x ln x')
d_ = max(abs(-4*(1 - x)*quad(lambda xi: (1 - xi)*I(xi, x, 0.0), 0, 1, epsabs=1e-12, epsrel=1e-12)[0] - 4*x*np.log(x))
         for x in (0.1, 0.3, 0.5, 0.7, 0.9))
check('x = 0.1, 0.3, 0.5, 0.7, 0.9', d_, 1e-10)

print('2  coefficient of eta^0: R(eps) = Int_0^1 dxi/xi [(1-xi)^(1-eps) I(xi) - I(0)]')


def R(x, e):
    """split used in the derivation: I(xi) with (1-t)^(-eps) replaced by 1 is elementary; the remainder is O(eps)"""
    a = (1/e)*x**e*beta(-e, 1 - e)
    b = (1/e)*quad(lambda xi: (1 - ((1 - xi)*(1 - xi + xi/x))**(-e))/xi, 0, 1, **Q)[0]
    rest = lambda xi: quad(lambda s: ((1 - s)**(-e) - 1)*(s*(1 - xi) + xi/x)**(-1 - e), 0, 1, **Q)[0]
    c = quad(lambda xi: ((1 - xi)**(1 - e)*rest(xi) - (I_at_0(e) + 1/e))/xi, 0, 1, epsabs=1e-11, epsrel=1e-11, limit=200)[0]
    return a + b + c


R_series = lambda x, e: -1/e**2 - np.log(x)/e - 0.5*np.log(x)**2 - Li2(1 - 1/x)
x, e = 0.3, -0.4
direct = quad(lambda xi: ((1 - xi)**(1 - e)*I(xi, x, e) - I_at_0(e))/xi, 0, 1, epsabs=1e-10, epsrel=1e-10, limit=400)[0]
check('x = 0.3, eps = -0.4: split form against direct quadrature', R(x, e) - direct, 1e-8)
for x in (0.3, 0.7):
    d_ = [R(x, e) - R_series(x, e) for e in (-0.04, -0.02, -0.01)]
    require(f'x = {x}: R approaches -1/eps^2 - ln x/eps - ln^2 x/2 - Li2(1 - 1/x) linearly in eps '
            f'(ratios {d_[0]/d_[1]:.3f}, {d_[1]/d_[2]:.3f})', abs(d_[2]) < 0.05 and 1.85 < d_[0]/d_[1] < 2.1 and 1.9 < d_[1]/d_[2] < 2.1)
# A = 2 F' [(-1/eta + ln(1-x)) I(0) + R], with F' I(0) = e^(-eps gamma_E) Gamma(-eps) = -1/eps + O(eps) and F' = 1 + (3/2) zeta2 eps^2
xs = np.array([0.1, 0.3, 0.5, 0.7, 0.9]); lx, lb = np.log(xs), np.log(1 - xs)
finite = 2*(-lx**2/2 - Li2(1 - 1/xs)) - 3*Z2            # eps^0; the poles are 2/(eta eps), -2/eps^2 and -2 (ln x + ln(1-x))/eps
check('finite part of A equals 2 Li2(1-x) - 3 zeta2', np.max(np.abs(finite - (2*Li2(1 - xs) - 3*Z2))), 1e-12)

print('3  results against the stored Laurent coefficients (five x, two helicity structures, six orders)')
worst = {'W2q': 0.0, 'W1qb': 0.0}
for tag, x in {'1o10': 0.1, '3o10': 0.3, '1o2': 0.5, '7o10': 0.7, '9o10': 0.9}.items():
    st = json.load(open(repo/f'paper2/r0/theory/gqq_laurent_x{tag}.json'))['diagrams']
    y = 1 - x; L = np.log(x*y)
    forms = {'W2q': (+1, {'-1,-1': 2, '0,-2': -2, '0,-1': -2*L, '0,0': 2*Li2(y) - 3*Z2}, 4*x*np.log(x)),
             'W1qb': (-1, {'-1,-1': 2, '0,-2': -2, '0,-1': -2*L, '0,0': 2*Li2(x) - 3*Z2}, -4*y*np.log(y))}   # stored prefactor of W1qb: -(C_F - C_A/2)
    for name, (sign, A, B) in forms.items():
        for key in ('-1,-2', '-1,-1', '-1,0', '0,-2', '0,-1', '0,0'):
            a = A.get(key, 0); b = B if key == '0,0' else 0
            for h, pred in (('ek', a - b/(2*x)), ('ke', a + b/(2*y))):
                worst[name] = max(worst[name], abs(sign*st[name]['laurent'][h][key][0] - pred))
for name, v in worst.items(): check(name, v, 1e-8)
finish()
