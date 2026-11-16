"""The collinear function J1(g -> gg) at one loop: all 23 entries of paper2/r0/theory/gg_diagrams.py, their colour
sum, and r0 against Eq. (7) of the Letter.

The diagrams are the gluon, ghost and quark triangles (Vtri, Vgh1, Vgh2, Vq1, Vq2), the four-gluon vertex bubbles on
either daughter (Vbuba1-3, Vbubb1-3), and the Wilson line couplings to the parent gluon (W1g, W2g), to the daughter
created by the same operator (W1k, W2kb) and to the other daughter (W1kb, W2k). The six entries in which the gluon
from the Wilson line enters a four-gluon vertex (W1v41-3, W2v41-3) vanish.

Helicity structures, with e the polarisation of the parent gluon and e_1, e_2 those of the daughters:
    S1 = (e.k)(e_1.e_2),  S2 = (e.e_1)(e_2.k),  S3 = (e.e_2)(e_1.k),  T4 = (e.k)(e_1.k)(e_2.k)/k^2.
The tree splitting is proportional to S1 - S2/(1-x) - S3/x. Below, c1, c2, c3 are the coefficients of S1, S2, S3
normalised to the tree splitting, and c4 is the coefficient of T4 in units of the tree coefficient of S1.
Units k^2 = mu^2 = nu_r = 1, a_s = alpha_s/(4 pi), impact parameter space at b = b0.

The momentum fraction x is kept as a symbol in steps 1, 4, 6 and 8, so the equalities there hold for all x.
Steps:
  1  numerators (exact tensor algebra)
  2  residues in the minus components, against quadrature in two variables
  3  parameter integrals, against quadrature
  4  the nine independent diagrams in closed form
  5  all entries against the stored Laurent coefficients
  6  colour sum
  7  daughter self-energies and soft function
  8  transform of T4, contraction of the two splittings, and r0 against Eq. (7) of the Letter for nu = 0 and 2,
     for scalar (gg -> HH) and pseudoscalar (gg -> AA) hard vertices
Usage: python3 check_gg.py <repository root>"""
import sys, json, pathlib, math, warnings
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from fractions import Fraction as F
import numpy as np
from scipy.integrate import quad, dblquad
from scipy.special import beta as Beta, spence, jv
from algebra import Ser, V_, X
from gluon_diagrams import build, to_r
from integrals import ParameterIntegrals, rt, LB, LX, LI, Z2, STRUCT, evaluate, average, ser_structures, structures, show_series
from gluon_channel import diagram, numerator, split_ll, at_one_minus_x
from supplement_formulas import soft, sigma, laurent
from report import check, require, finish
warnings.filterwarnings('ignore')
repo = pathlib.Path(sys.argv[1]); data = repo/'paper2/r0/theory'
TREE = {'S1': 1 + 0*X, 'S2': -1/(1 - X), 'S3': -1/X}          # tree splitting
L = LX + LB


def S(d): return Ser({n: (v if hasattr(v, 't') else rt(v)) for n, v in d.items()})


def closed(name):
    """(coefficient of 1/eta, coefficient of eta^0) of the nine independent diagrams: {helicity structure: series in eps}"""
    x = X; t2, t3 = TREE['S2'], TREE['S3']
    if name == 'W1k':
        U = S({-2: 2, -1: 1 + 2*LB, 0: 2 + Z2 + LB + LB*LB})
        return {}, {'S1': U, 'S2': U*t2, 'S3': U*t3}
    if name.startswith('Vbuba'):
        U = S({-1: 1, 0: 2 + LB})*{'1': F(3, 4), '2': F(3, 8), '3': F(-3, 8)}[name[-1]]
        return {}, {'S1': U, 'S2': U*t2}
    if name == 'W1g':
        return {}, {'S1': S({-2: -2, -1: -3*LB, 0: -3*Z2 + 2*LI - LB*F(1, 2) - LB*LB*F(3, 2)}),
                    'S2': S({-2: -2, -1: -(2*LB + x), 0: -3*Z2 + 2*LI - LB - LB*LB - 2*x})*t2,
                    'S3': S({-2: -2, -1: -(LB*F(5, 2) + x/2), 0: -3*Z2 + 2*LI + LB*(F(3, 4)*(1 - x)) - LB*LB*F(5, 4) - F(5, 4)*x})*t3}
    if name == 'W1kb':
        P = S({-1: -2})
        return ({'S1': P, 'S2': P*t2, 'S3': P*t3},
                {'S1': S({-2: 2, -1: 2*L, 0: 3*Z2 - 2*LI + LB}),
                 'S2': S({-2: 2, -1: 2*L + LB + x/2, 0: 3*Z2 - 2*LI + LB*F(1, 2) + LB*LB*F(1, 2) + x})*t2,
                 'S3': S({-2: 2, -1: 2*L, 0: 3*Z2 - 2*LI + LB*(F(7, 4) - F(3, 4)*x) - x/4})*t3})
    if name == 'Vtri':
        return {}, {'S1': S({-1: F(5, 4) - L, 0: F(7, 3) - L*F(3, 2) - (LX*LX + LB*LB)*F(1, 2)}),
                    'S2': S({-1: x - F(1, 2) - LX*F(3, 2) - LB,
                             0: F(5, 4)*x - F(5, 12) - LX*F(17, 4) - LX*LX*F(3, 4) - LB*F(3, 2) - LB*LB*F(1, 2)})*t2,
                    'S3': S({-1: F(1, 2) - x - LB*F(3, 2) - LX,
                             0: F(5, 6) - F(5, 4)*x - LB*F(17, 4) - LB*LB*F(3, 4) - LX*F(3, 2) - LX*LX*F(1, 2)})*t3,
                    'T4': S({0: F(3, 2)})}
    if name == 'Vgh1':
        U = S({-1: F(-1, 24), 0: L*F(-1, 24) - F(1, 9)})
        return {}, {'S1': U, 'S2': (U - S({0: (1 - x)/24}))*t2, 'S3': (U - S({0: x/24}))*t3, 'T4': S({0: F(-1, 12)})}
    if name == 'Vq1':        # coefficient of -i
        U = S({-1: F(-4, 3), 0: L*F(-4, 3) - F(20, 9)})
        return {}, {'S1': U, 'S2': (U + S({0: F(2, 3)*(1 - x)}))*t2, 'S3': (U + S({0: F(2, 3)*x}))*t3, 'T4': S({0: F(4, 3)})}
    raise KeyError(name)


def J1_closed():
    """colour sum: {'CA': .., 'nf': ..} for the coefficient of 1/eta and for the coefficient of eta^0"""
    rho = {'S1': 0*X, 'S2': 1 - X, 'S3': X}
    pole = {'CA': {s: S({-1: 2})*TREE[s] for s in TREE}, 'nf': {}}
    fin = {'CA': {s: S({-1: F(1, 3) - 2*L, 0: F(5, 9) - 2*Z2 - L*F(5, 3) - 2*LX*LB + rho[s]/3})*TREE[s] for s in TREE},
           'nf': {s: S({-1: F(4, 3), 0: F(20, 9) + L*F(4, 3) - rho[s]*F(2, 3)})*TREE[s] for s in TREE}}
    fin['CA']['T4'] = S({0: F(2, 3)}); fin['nf']['T4'] = S({0: F(-4, 3)})
    return pole, fin


def same(a, b):
    """equality of two {helicity structure: series} at the orders eps^-2, eps^-1, eps^0, as functions of x"""
    for s in set(a) | set(b):
        d = a.get(s, Ser()) - b.get(s, Ser())
        if any(n <= 0 and c for n, c in d.c.items()): return False
    return True


print('1  numerators, with x a symbol')
G = build(X)
tree = structures(numerator('tree', X).subs('u', -V_('kk')/(1 - X)))
t1 = tree['S1'].const_value()
require(f'tree splitting on the pole of (p1 - K)^2: {t1} [S1 - S2/(1-x) - S3/x]',
        t1 == -2 and tree['S2'].const_value() == t1*TREE['S2'] and tree['S3'].const_value() == t1*TREE['S3'] and 'T4' not in tree)
require('gluon from the Wilson line entering a four-gluon vertex: the six numerators vanish', all(not numerator(n, X) for n in G if 'v4' in n))
require('the two orientations of the ghost loop have opposite numerators', not (numerator('Vgh1', X) + numerator('Vgh2', X)))
require('triangles: numerators linear in K^- and in l^-, with no term K^- l^-',
        all(sorted(numerator(n, X).collect(['u', 'w'])) == [(0, 0), (0, 1), (1, 0)] for n in ('Vtri', 'Vq1', 'Vq2')))
require('all other numerators are independent of K^- and l^-',
        all(numerator(n, X).degree('u') <= 0 and numerator(n, X).degree('w') <= 0 for n in G if n not in ('Vtri', 'Vq1', 'Vq2', 'tree')))

print('2  residues in the minus components, against quadrature in two variables with finite imaginary parts (x = 0.3, z = 0.12)')
QK = dict(complex_func=True, epsabs=1e-11, epsrel=1e-10, limit=400)
xf, zz = 0.3, 0.12; xbf = 1 - xf
# bubble and W1k: denominators l^2, (K-l)^2, K^2, (p1-K)^2; K^- on the pole of (p1-K)^2, l^- on the pole of l^2
m = np.array([0.45 - 0.35j, 0.7 - 0.30j, 1.0 - 0.40j, 1.0 - 0.50j])
inner = lambda u: quad(lambda w: 1/((zz*w - m[0])*((xf - zz)*(u - w) - m[1])), -np.inf, np.inf, **QK)[0]
numv = quad(lambda u: inner(u)/((xf*u - m[2])*(-xbf*u - m[3])), -np.inf, np.inf, **QK)[0]
ub, wd = -m[3]/xbf, m[0]/zz
res = (2j*np.pi)*(-1/xbf)/(xf*ub - m[2])*(-2j*np.pi/zz)/((xf - zz)*(ub - wd) - m[1])
check('bubble: quadrature against residues', abs(numv - res), 1e-9)
kk, qq, ll = 1.0, 0.7, 0.45; ze = zz/xf; D = (1 - ze)*ll + ze*qq - ze*(1 - ze)*kk + ze*(1 - ze)*kk/xbf
ub, wd = -kk/xbf, ll/zz
res0 = ((2j*np.pi)*(-1/xbf)/(xf*ub - kk)*(-2j*np.pi/zz)/((xf - zz)*(ub - wd) - qq)/(2*np.pi)**2).real
check('bubble: residues/(2 pi)^2 = -1/(x k^2 D)', res0 + 1/(xf*kk*D), 1e-12)
# W1kb: denominators (K-l)^2, (p1-K+l)^2, l^2, (p1-K)^2; K^- on the pole of (K-l)^2, l^- on the pole of l^2
m = np.array([0.7 - 0.30j, 0.7 - 0.45j, 0.45 - 0.35j, 1.0 - 0.50j])
inner_num = lambda w: quad(lambda u: 1/(((xf - zz)*(u - w) - m[0])*((xbf + zz)*(w - u) - m[1])*(-xbf*u - m[3])), -np.inf, np.inf, **QK)[0]


def inner_res(w, mm):
    u1 = w + mm[0]/(xf - zz)
    return (-2j*np.pi/(xf - zz))/(((xbf + zz)*(w - u1) - mm[1])*(-xbf*u1 - mm[3]))


both_res = lambda mm: inner_res(mm[2]/zz, mm)*(-2j*np.pi/zz)
check('W1kb: K^- integral against its residue at three values of l^-', max(abs(inner_num(w_) - inner_res(w_, m)) for w_ in (-3.0, 0.4, 5.0)), 1e-9)
check('W1kb: l^- integral against its residue', abs(quad(lambda w: inner_res(w, m)/(zz*w - m[2]), -np.inf, np.inf, **QK)[0] - both_res(m)), 1e-9)
check('W1kb: residues/(2 pi)^2 = -(1-zeta)/((1-x) q_perp^2 D), zeta = z/x',
      (both_res(np.array([qq, qq, ll, kk]))/(2*np.pi)**2).real + (1 - ze)/(xbf*qq*D), 1e-12)

print('3  parameter integrals against quadrature (x = 0.3), V = 1 - x zeta - (1-x) s')
M = ParameterIntegrals(X); Vf = lambda z_, s_: 1 - xf*z_ - xbf*s_; worst = 0
for a in (-1, 0, 1, 2):
    for b in (0, 1, 2):
        m0 = evaluate(M.M0(a, b), xf)
        n0 = dblquad(lambda s_, z_: z_**a*s_**b, 0, 1, 0, lambda z_: z_)[0]
        n1 = -dblquad(lambda s_, z_: z_**a*s_**b*(np.log(s_) + np.log(Vf(z_, s_))), 0, 1, 0, lambda z_: z_, epsabs=1e-11, epsrel=1e-11)[0]
        worst = max(worst, abs(m0.get(0, 0) - n0), abs(m0.get(1, 0) - n1))
        if b >= 1:
            q_ = dblquad(lambda s_, z_: z_**a*s_**(b - 1)/Vf(z_, s_), 0, 1, 0, lambda z_: z_, epsabs=1e-11, epsrel=1e-11)[0]
            worst = max(worst, abs(evaluate(M.M1(a, b), xf)[0] - q_))
for a in (0, 1, 2):
    m1 = evaluate(M.M1(a, 0), xf)
    pole = -quad(lambda z_: z_**a/(1 - xf*z_), 0, 1)[0]
    fin = quad(lambda z_: z_**a*(np.log(z_) + 2*np.log(1 - xf*z_) - np.log(1 - z_))/(1 - xf*z_), 0, 1, epsabs=1e-12, epsrel=1e-12)[0]
    worst = max(worst, abs(m1[-1] - pole), abs(m1[0] - fin))
check('Int dzeta zeta^a Int ds s^b (s V)^(-eps) to order eps, and Int dzeta zeta^a Int ds s^(b-1-eps) V^(-1-eps) to order eps^0', worst, 1e-10)


def double_pole(ep):
    """Int_0^1 dzeta/zeta Int_0^zeta ds s^(-1-eps) V^(-1-eps) at eps < 0, with s = e^(-y)"""
    def g(s_):
        a_ = 1 - xbf*s_
        return quad(lambda z_: ((1 - xf*z_ - xbf*s_)**(-1 - ep) - a_**(-1 - ep))/z_, s_, 1, epsabs=1e-13, epsrel=1e-12)[0], a_**(-1 - ep)
    f = lambda y: np.exp(ep*y)*(g(np.exp(-y))[0] + y*g(np.exp(-y))[1])
    return quad(f, 0, np.inf, epsabs=1e-10, epsrel=1e-10, limit=400)[0]


mX = evaluate(M.M1(-1, 0), xf)
d_ = [double_pole(ep) - (mX[-2]/ep**2 + mX[-1]/ep + mX[0]) for ep in (-0.04, -0.02, -0.01)]
require(f'a = -1, b = 0: 1/eps^2 + ln(1-x)/eps + zeta2 - 2 Li2(x) - ln^2(1-x)/2, approached linearly in eps (ratios {d_[0]/d_[1]:.3f}, {d_[1]/d_[2]:.3f})',
        abs(d_[2]) < 0.05 and 1.85 < d_[0]/d_[1] < 2.1 and 1.9 < d_[1]/d_[2] < 2.1)
ep, et = -0.3, -0.2; worst = 0
for b in (0, 1):
    lhs = quad(lambda z_: (1 - z_)**(-1 - et)*quad(lambda s_: s_**(b - 1 - ep)*Vf(z_, s_)**(-1 - ep), 0, z_, epsabs=1e-12, epsrel=1e-11)[0], 0, 1,
               epsabs=1e-10, epsrel=1e-10, limit=300)[0]
    rhs = Beta(-ep - et, b - ep)*quad(lambda r: r**(-1 - et)*(xbf + xf*r)**(-1 - ep), 0, 1, epsabs=1e-12, epsrel=1e-12)[0]
    worst = max(worst, abs(lhs/rhs - 1))
check('W1kb: Int dzeta (1-zeta)^(-1-eta) Int ds s^(b-1-eps) V^(-1-eps) = B(-eps-eta, b-eps) Int drho rho^(-1-eta) [(1-x) + x rho]^(-1-eps), '
      'relative, at eps = -0.3, eta = -0.2', worst, 1e-6)
g1 = lambda e_: quad(lambda r: ((1 + xf/xbf*r)**(-1 - e_) - 1)/r, 0, 1, epsabs=1e-13, epsrel=1e-13)[0]
check('W1kb: Int drho/rho {[1 + x rho/(1-x)]^(-1-eps) - 1} = ln(1-x) - eps Li2(x) + O(eps^2)',
      max(abs(g1(0) - math.log(xbf)), abs((g1(1e-4) - g1(-1e-4))/2e-4 + spence(1 - xf))), 1e-8)

print('4  the nine independent diagrams in closed form, with x a symbol (c4 is absent where not listed)')
for name in ('W1k', 'Vbuba1', 'Vbuba2', 'Vbuba3', 'W1g', 'W1kb', 'Vtri', 'Vgh1', 'Vq1'):
    got = diagram(name, X); want = closed(name)
    require(f'{name}: computed result equals the closed form' + (' (factor -i left out)' if name == 'Vq1' else ''),
            same(got[0], want[0]) and same(got[1], want[1]))
    for label, part in (('1/eta', got[0]), ('eta^0', got[1])):
        for st, cn in (('S1', 'c1'), ('S2', 'c2'), ('S3', 'c3'), ('T4', 'c4')):
            if st not in part: continue
            norm = part[st]*(1/TREE[st]) if st in TREE else part[st]
            if any(norm.c.get(n) for n in (-2, -1, 0)): print(f'       {label}  {cn}:  {show_series(norm)}')
# rapidity divergence of W1kb: at zeta = 1 the numerator is s times the tree splitting
N0, N1 = split_ll(to_r(numerator('W1kb', X).subs('z', X*(1 - V_('ze')))))
st = ser_structures(average(N0, V_('s'), 2)[-1])
at1 = {s_: {k_: sum(c.const_value() for (a_, b_), c in p.c[0].collect(['ze', 's']).items() if b_ == k_) for k_ in (0, 1)} for s_, p in st.items()}
require('W1kb: at zeta = 1 the numerator is s times the tree splitting, so the coefficient of 1/eta is 2 e^(-eps gamma_E) Gamma(-eps) to all orders in eps',
        all(not at1[s_][0] for s_ in at1) and at1['S2'][1] == at1['S1'][1]*TREE['S2'] and at1['S3'][1] == at1['S1'][1]*TREE['S3'])

print('5  all entries against the stored Laurent coefficients (five x, four helicity structures, eta^-1 and eta^0, eps^-2 to eps^0)')
RESULT = {n: diagram(n, X) for n in G if n != 'tree'}


def compare(tag, x, names=None):
    stored = json.load(open(data/f'gg_laurent_x{tag}.json'))['diagrams']; worst = {}
    for name in (names or stored):
        pole, fin = RESULT[name]; imaginary = name.startswith('Vq')
        for st in STRUCT:
            vp = evaluate(pole[st], x) if st in pole else {}; vf = evaluate(fin[st], x) if st in fin else {}
            for n in (-2, -1, 0):
                sp_ = stored[name]['laurent'][st].get(f'-1,{n}', [0, 0]); sf = stored[name]['laurent'][st].get(f'0,{n}', [0, 0])
                sp_, sf = (-sp_[1], -sf[1]) if imaginary else (sp_[0], sf[0])
                worst[name] = max(worst.get(name, 0), abs(vp.get(n, 0) - sp_), abs(vf.get(n, 0) - sf))
    return worst


total = {}
for tag, x in {'1o10': 0.1, '3o10': 0.3, '1o2': 0.5, '7o10': 0.7, '9o10': 0.9}.items():
    for k_, v in compare(tag, x).items(): total[k_] = max(total.get(k_, 0), v)
require(f'the stored files hold the same {len(RESULT)} entries', set(total) == set(RESULT))
for k_, v in total.items():
    if 'v4' in k_: require(f'{k_}: stored coefficients vanish', v == 0)
    else: check(k_, v, 2e-6 if k_ == 'Vtri' else 5e-8)
check('Vtri at x = 1/2, stored with quadrature orders 96, 72, 40', compare('1o2_nq96-72-40', 0.5, ['Vtri'])['Vtri'], 5e-8)

print('6  colour sum, with x a symbol')
colour = {n: (v['colour']['CA'][0], v['colour']['nf'][1]) for n, v in json.load(open(data/'gg_laurent_x3o10.json'))['diagrams'].items()}
print('       colour factors in units of C_A (the quark triangles in units of T_F n_f, with the factor -i absorbed):')
items = [f'{n} {F(c[0] if c[0] else c[1]).limit_denominator(100)}' for n, c in colour.items() if 'v4' not in n]
for i in range(0, len(items), 9): print('         ' + ', '.join(items[i:i + 9]))
poleJ, finJ = J1_closed()
for ci, cname, text in ((0, 'CA', 'C_A:      tree splitting times 2/(eta eps) + (1/3 - 2L)/eps + 5/9 - 2 zeta2 - (5/3) L - 2 ln x ln(1-x) + rho/3,  c4 = 2/3'),
                        (1, 'nf', 'T_F n_f:  tree splitting times 4/(3 eps) + 20/9 + (4/3) L - (2/3) rho,  c4 = -4/3')):
    sum_pole, sum_fin = {}, {}
    for name, c in colour.items():
        if c[ci] == 0: continue
        p, f = RESULT[name]; w_ = F(c[ci]).limit_denominator(100)
        for s_, v in p.items(): sum_pole[s_] = sum_pole.get(s_, Ser()) + v*w_
        for s_, v in f.items(): sum_fin[s_] = sum_fin.get(s_, Ser()) + v*w_
    require(text, same(sum_pole, poleJ[cname]) and same(sum_fin, finJ[cname]))
print('       with L = ln(x(1-x)) and rho = 0, 1-x, x for S1, S2, S3')

print('7  daughter self-energies, Eq. (eq:sigma), and soft function, Eq. (eq:softeta), of the Supplement')
z2 = math.pi**2/6; worst = 0
for xv in (0.1, 0.3, 0.5, 0.8):
    Lv = math.log(xv*(1 - xv))
    for i, (p_, f_) in enumerate(((10/3, 62/9 + 5/3*Lv), (-8/3, -40/9 - 4/3*Lv))):
        c = laurent(lambda e: sigma(e, xv, 'gluon')[i] + sigma(e, 1 - xv, 'gluon')[i])
        worst = max(worst, abs(c[-2]), abs(c[-1] - p_), abs(c[0] - f_))
check('the two daughters: C_A [10/(3 eps) + 62/9 + (5/3) L] + T_F n_f [-8/(3 eps) - 40/9 - (4/3) L]', worst, 1e-10)
c1, c0 = laurent(lambda e: soft(e)[0]), laurent(lambda e: soft(e)[1])
check('S/Delta_C = -4/(eta eps) + 2/eps^2 - zeta2', max(abs(c1[-2]), abs(c1[-1] + 4), abs(c1[0]), abs(c0[-2] - 2), abs(c0[-1]), abs(c0[0] + z2)), 1e-10)
eta = 0.3
check('rapidity integral Int_0^inf dt/t |1/t - t|^(-eta) = B(eta/2, 1-eta) at eta = 0.3',
      2*quad(lambda t: t**(eta - 1)*(1 - t*t)**(-eta), 0, 1, epsabs=1e-13, epsrel=1e-13)[0] - Beta(eta/2, 1 - eta), 1e-8)
d_ = [Beta(e_/2, 1 - e_) - 2/e_ for e_ in (0.02, 0.01, 0.005)]
require(f'B(eta/2, 1-eta) = 2/eta + O(eta) (ratios {d_[0]/d_[1]:.3f}, {d_[1]/d_[2]:.3f})', 1.9 < d_[0]/d_[1] < 2.1 and 1.9 < d_[1]/d_[2] < 2.1)
require('rapidity divergences cancel: J1 + J2 gives 4 C_A/(eta eps), the soft function -4 Delta_C/(eta eps) with Delta_C = C_A',
        same(poleJ['CA'], {s: S({-1: 2})*TREE[s] for s in TREE}) and not poleJ['nf'] and abs(c1[-1] + 4) < 1e-10)
require('single pole of J1 + J2 + the four daughter self-energies: (B - 4 C_A L)/eps times the tree splitting, with B = 2 beta_0 as in the Letter',
        all(not (finJ['CA'][s].c[-1]*(2/TREE[s]) + F(20, 3) - (F(22, 3) - 4*L)) and not (finJ['nf'][s].c[-1]*(2/TREE[s]) - F(16, 3) + F(8, 3)) for s in TREE))

print('8  r0 from J1, J2, the soft function and the daughter self-energies, against Eq. (7) of the Letter, with x a symbol')
# transform of T4 to impact parameter space at eps = 0
def third_harmonic(k):
    """T4 - (S1 + S2 + S3)/4 as a tensor in the polarisations of the parent gluon and the two daughters"""
    k2 = k[0]*k[0] + k[1]*k[1]; dl = lambda a, b: 1 if a == b else 0
    return {(a, i, j): F(k[a]*k[i]*k[j], k2) - F(k[a]*dl(i, j) + dl(a, i)*k[j] + dl(a, j)*k[i], 4)
            for a in (0, 1) for i in (0, 1) for j in (0, 1)}


ok = True
for kv in ((1, 0), (2, 3), (-1, 5)):
    R = third_harmonic(kv)
    ok &= all(R[(a, i, j)] == R[(i, a, j)] == R[(a, j, i)] for (a, i, j) in R) and all(R[(0, 0, a)] + R[(1, 1, a)] == 0 for a in (0, 1))
require('T4 - (S1 + S2 + S3)/4 is symmetric and traceless in two transverse dimensions: the third harmonic in the direction of k', ok)
a_ = 0.1; worst = 0
for m_ in (1, 3):
    num_ = quad(lambda t: np.exp(-a_*t)*jv(m_, t), 0, 600, limit=4000, epsabs=1e-12, epsrel=1e-12)[0]
    worst = max(worst, abs(num_ - (math.sqrt(1 + a_*a_) - a_)**m_/math.sqrt(1 + a_*a_)))
check('Int_0^inf dt e^(-a t) J_m(t) = [sqrt(1+a^2) - a]^m/sqrt(1+a^2) for m = 1, 3 at a = 0.1; both tend to 1 as a -> 0', worst, 1e-8)
print('       The transform of the harmonic m with radial dependence 1/|k| is i^m/b, so the third harmonic changes sign')
print('       relative to the tree splitting: T4 -> (S1 + S2 + S3)/4 - [T4 - (S1 + S2 + S3)/4] = (S1 + S2 + S3)/2 - T4.')


# contraction of the two splittings over the polarisations of the daughters
def matrix(s_, e):
    """a helicity structure as a matrix in the polarisations of the two daughters, for b = (1, 0) and parent polarisation e"""
    if s_ == 'S1': return [[e[0], 0], [0, e[0]]]
    if s_ == 'S2': return [[e[0], 0], [e[1], 0]]
    if s_ == 'S3': return [[e[0], e[1]], [0, 0]]
    return [[e[0], 0], [0, 0]]


ANTI = [[0, 1], [-1, 0]]


def overlaps(process):
    """structure s of the first splitting times structure t of the second, summed over the polarisations of the daughters.
    Each hard vertex joins one daughter of each splitting: with delta_ij for gg -> HH, with epsilon_ij for gg -> AA.
    The sum is alpha (e.e') + beta (e.b)(e'.b)/b^2, with harmonics alpha + beta/2 for nu = 0 and beta for nu = 2."""
    def Q(s_, t_, e, e2):
        A, B = matrix(s_, e), matrix(t_, e2)
        if process == 'HH': return sum(F(A[i][j]*B[i][j]) for i in (0, 1) for j in (0, 1))
        return sum(F(ANTI[i][a]*ANTI[j][b]*A[i][j]*B[a][b]) for i in (0, 1) for j in (0, 1) for a in (0, 1) for b in (0, 1))
    ex, ey = (1, 0), (0, 1); out = {0: {}, 2: {}}; diagonal = True
    for s_ in STRUCT:
        for t_ in STRUCT:
            diagonal &= Q(s_, t_, ex, ey) == 0 and Q(s_, t_, ey, ex) == 0
            out[0][(s_, t_)] = (Q(s_, t_, ex, ex) + Q(s_, t_, ey, ey))/2; out[2][(s_, t_)] = Q(s_, t_, ex, ex) - Q(s_, t_, ey, ey)
    return out, diagonal


swap = {'S1': 'S1', 'S2': 'S3', 'S3': 'S2', 'T4': 'T4'}
require('the colour sum is symmetric under the exchange of the two daughters: x <-> 1-x with S2 <-> S3',
        all(not (at_one_minus_x(finJ[c_][swap[s_]].c[n]) - finJ[c_][s_].c[n]) for c_ in ('CA', 'nf') for s_ in STRUCT for n in finJ[c_][s_].c))
T1 = dict(TREE, T4=0*X); T2 = {s_: T1[s_].reflect() for s_ in STRUCT}        # in the second splitting the daughter that
J1b, J2b = {}, {}                                                             # meets the one with fraction x carries 1 - x
for c_ in ('CA', 'nf'):
    a = {s_: finJ[c_][s_].c[0] for s_ in STRUCT}
    J1b[c_] = {s_: a[s_] + a['T4']*F(1, 2) for s_ in TREE}; J1b[c_]['T4'] = -a['T4']          # impact parameter space
    J2b[c_] = {s_: at_one_minus_x(v) for s_, v in J1b[c_].items()}
inv = 1/(X*(1 - X))
F0, FA, F2 = inv + 2*X*(1 - X), inv - 2*X*(1 - X), (X**4 + (1 - X)**4)*inv
# name of the overlap, the overlap, tree overlap divided by it, sign of the last term of Eq. (7)
CASES = {('HH', 0): ('F_0', F0, inv*F(1, 2), 1, 'F_0/(2x(1-x))'), ('HH', 2): ('F_2', F2, inv, 1, 'F_2/(x(1-x))'),
         ('AA', 0): ('F_A', FA, inv*F(-1, 2), 1, '-F_A/(2x(1-x))'), ('AA', 2): ('F_2', F2, inv, -1, 'F_2/(x(1-x))')}
x_dependent = {'CA': F(134, 9) - 5*Z2 - 4*LX*LB, 'nf': rt(F(-40, 9))}                       # Eq. (7) without its last term
last_term = {'CA': F(2, 3), 'nf': F(-4, 3)}                                                 # 2 (C_A - 2 T_F n_f)/3
self_energies = {'CA': 2*(F(62, 9) + L*F(5, 3)), 'nf': 2*(F(-40, 9) - L*F(4, 3))}           # four daughters
soft_finite = {'CA': -Z2, 'nf': rt(0)}                                                      # -zeta2 Delta_C, Delta_C = C_A
for process in ('HH', 'AA'):
    O, diagonal = overlaps(process)
    require(f"gg -> {process}: every overlap has the form alpha (e.e') + beta (e.b)(e'.b)/b^2", diagonal)
    for nu in (0, 2):
        name, Fnu, ratio, sign, text = CASES[(process, nu)]
        tree_overlap = 0*X
        for s_ in STRUCT:
            for t_ in STRUCT: tree_overlap = tree_overlap + T1[s_]*T2[t_]*O[nu][(s_, t_)]
        require(f'gg -> {process}, nu = {nu}: tree overlap = {text}', tree_overlap == ratio*Fnu)
        for c_, cname in (('CA', 'C_A'), ('nf', 'T_F n_f')):
            one_loop = rt(0)
            for s_ in STRUCT:
                for t_ in STRUCT:
                    if O[nu][(s_, t_)]: one_loop = one_loop + (J1b[c_][s_]*T2[t_] + J2b[c_][t_]*T1[s_])*O[nu][(s_, t_)]
            # r0 = one_loop/tree_overlap + self-energies + soft function; tree_overlap/F_nu = ratio
            residual = one_loop + (self_energies[c_] + soft_finite[c_] - x_dependent[c_])*tree_overlap - rt(ratio*(sign*last_term[c_]))
            form = '134/9 - 5 zeta2 - 4 ln x ln(1-x)' if c_ == 'CA' else '-40/9'
            require(f'gg -> {process}, nu = {nu}, {cname} part: {form} {"+" if sign*last_term[c_] > 0 else "-"} ({abs(last_term[c_])})/{name}', not residual)
print('       Eq. (7) holds for nu = 0 and 2 as an equality in x. For gg -> AA, F_A replaces F_0 and the last term changes sign for nu = 2.')
finish()
