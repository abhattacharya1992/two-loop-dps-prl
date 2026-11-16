"""The collinear function J1(g -> q qbar) at one loop: sum of the eight diagrams, and r0 against Eq. (8) of the Letter.

Each diagram is (colour factor) [A T + (B/2) (e.k)], with the tree splitting T = (x/2) eslash kslash - ((1-x)/2) kslash eslash.
The coefficients of the two quark helicity structures, normalised to the tree splitting, are A - B/(2x) for
eslash kslash and A + B/(2(1-x)) for kslash eslash. Units k^2 = mu^2 = nu_r = 1, a_s = alpha_s/(4 pi), impact parameter
space at b = b0. The closed forms of A and B are derived in the scripts check_qq_*.py of the individual diagrams.

The momentum fraction x is kept as a symbol in steps 2 and 4, so the equalities there hold for all x.
Steps:
  1  the eight closed forms against the stored Laurent coefficients
  2  sum of the eight diagrams
  3  daughter self-energies and soft function
  4  r0 against Eq. (8) of the Letter
Usage: python3 check_qq_sum.py <repository root>"""
import sys, json, pathlib, math
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from fractions import Fraction as F
from algebra import Ser, X
from integrals import rt, LB, LX, LI, Z2, evaluate, show, show_series
from supplement_formulas import soft, sigma, laurent
from report import check, require, finish
repo = pathlib.Path(sys.argv[1])
L = LX + LB
LI1 = Z2 - LX*LB - LI                  # Li2(1-x)


def S(d): return Ser({n: (v if hasattr(v, 't') else rt(v)) for n, v in d.items()})


x = X
# name: (colour factor as coefficients of C_F and C_A, coefficient of 1/eta in A, A, B)
DIAGRAMS = {
    'W1q':  ((1, 0), {}, {-2: 2, -1: 2 + 2*LB, 0: 4 + Z2 + 2*LB + LB*LB}, {}),
    'W2qb': ((1, 0), {}, {-2: 2, -1: 2 + 2*LX, 0: 4 + Z2 + 2*LX + LX*LX}, {}),
    'W2q':  ((1, F(-1, 2)), {-1: 2}, {-2: -2, -1: -2*L, 0: 2*LI1 - 3*Z2}, {0: 4*x*LX}),
    'W1qb': ((1, F(-1, 2)), {-1: 2}, {-2: -2, -1: -2*L, 0: 2*LI - 3*Z2}, {0: -4*(1 - x)*LB}),
    'Vqq':  ((1, F(-1, 2)), {}, {-1: 1, 0: -L}, {0: -4*x*LX + 4*(1 - x)*LB}),
    'Vggq': ((0, F(1, 2)), {}, {-1: 2 - 2*L, 0: 4 - LX*LX - LB*LB}, {-1: -2*x*LX + 2*(1 - x)*LB, 0: -x*LX*LX + (1 - x)*LB*LB}),
    'W1g':  ((0, F(1, 2)), {}, {-2: -2, -1: -(2*LB + x), 0: -3*Z2 - LB - LB*LB + 2*LI - 2*x},
             {-1: 2*x*(1 - x) - 2*(1 - x)*LB, 0: 4*x*(1 - x) - (1 - x)*LB*LB}),
    'W2g':  ((0, F(1, 2)), {}, {-2: -2, -1: -(2*LX + (1 - x)), 0: -3*Z2 - LX - LX*LX + 2*LI1 - 2*(1 - x)},
             {-1: -2*x*(1 - x) + 2*x*LX, 0: -4*x*(1 - x) + x*LX*LX}),
}
DIAGRAMS = {n: (c, S(p), S(a), S(b)) for n, (c, p, a, b) in DIAGRAMS.items()}

print('1  the eight closed forms against the stored Laurent coefficients (five x, two helicity structures, six orders)')
worst = {n: 0.0 for n in DIAGRAMS}; signs = {}
for tag, xv in {'1o10': 0.1, '3o10': 0.3, '1o2': 0.5, '7o10': 0.7, '9o10': 0.9}.items():
    stored = json.load(open(repo/f'paper2/r0/theory/gqq_laurent_x{tag}.json'))['diagrams']
    for name, (colour, pole, A, B) in DIAGRAMS.items():
        pref = [F(stored[name]['pref'][k_]) for k_ in ('CF', 'CA')]
        sign = 1 if pref == [F(c) for c in colour] else -1 if pref == [-F(c) for c in colour] else 0     # stored prefactor = sign times colour factor
        signs[name] = sign
        vp, va, vb = evaluate(pole, xv), evaluate(A, xv), evaluate(B, xv)
        for n in (-2, -1, 0):
            a, b = va.get(n, 0), vb.get(n, 0)
            for h, pred in (('ek', a - b/(2*xv)), ('ke', a + b/(2*(1 - xv)))):
                worst[name] = max(worst[name], abs(sign*stored[name]['laurent'][h][f'0,{n}'][0] - pred),
                                  abs(sign*stored[name]['laurent'][h][f'-1,{n}'][0] - vp.get(n, 0)))
require('colour factors agree with the stored prefactors up to the sign convention of each entry', all(signs.values()))
for name, v in worst.items(): check(name, v, 5e-7 if name == 'Vqq' else 1e-8)

print('2  sum of the eight diagrams, with x a symbol')
total = {}
for i, col in enumerate(('C_F', 'C_A')):
    for j, part in enumerate(('pole', 'A', 'B')):
        tot = Ser()
        for name, entry in DIAGRAMS.items():
            if entry[0][i]: tot = tot + entry[1 + j]*F(entry[0][i])
        total[(col, part)] = tot
zero = lambda ser: not any(n <= 0 and c for n, c in ser.c.items())
for col in ('C_F', 'C_A'):
    print(f'       {col}:  1/eta  {show_series(total[(col, "pole")])}')
    print(f'             eta^0  {show_series(total[(col, "A")])}')
require('C_F part: 4/(eta eps) + (5 - 2L)/eps + 8 - 2 zeta2 + L + ln^2(x/(1-x)),  L = ln(x(1-x))',
        zero(total[('C_F', 'pole')] - S({-1: 4})) and zero(total[('C_F', 'A')] - S({-1: 5 - 2*L, 0: 8 - 2*Z2 + L + (LX - LB)*(LX - LB)})))
require('C_A part: -2/(eta eps) + 1 - ln^2 x - ln^2(1-x)',
        zero(total[('C_A', 'pole')] - S({-1: -2})) and zero(total[('C_A', 'A')] - S({0: 1 - LX*LX - LB*LB})))
require('B cancels in the C_F part and in the C_A part, poles and finite part', zero(total[('C_F', 'B')]) and zero(total[('C_A', 'B')]))

print('3  daughter self-energies, Eq. (eq:sigma), and soft function, Eq. (eq:softeta), of the Supplement')
z2 = math.pi**2/6; w_ = 0
for xv in (0.1, 0.3, 0.5, 0.8):
    c = laurent(lambda e: sigma(e, xv, 'quark') + sigma(e, 1 - xv, 'quark'))
    w_ = max(w_, abs(c[-2]), abs(c[-1] + 2), abs(c[0] - (-2 - math.log(xv*(1 - xv)))))
check('quark and antiquark: C_F [-2/eps - 2 - L]', w_, 1e-10)
c1, c0 = laurent(lambda e: soft(e)[0]), laurent(lambda e: soft(e)[1])
check('S/Delta_C = -4/(eta eps) + 2/eps^2 - zeta2', max(abs(c1[-2]), abs(c1[-1] + 4), abs(c1[0]), abs(c0[-2] - 2), abs(c0[-1]), abs(c0[0] + z2)), 1e-10)
require('rapidity divergences cancel: J1 + J2 gives 4 Delta_C/(eta eps) with Delta_C = 2 C_F - C_A, the soft function -4 Delta_C/(eta eps)',
        zero(total[('C_F', 'pole')]*2 - S({-1: 8})) and zero(total[('C_A', 'pole')]*2 - S({-1: -4})) and abs(c1[-1] + 4) < 1e-10)
require('single pole of J1 + J2 + the four daughter self-energies: (B - 4 C_F L)/eps with B = 6 C_F as in the Letter',
        zero(S({-1: total[('C_F', 'A')].c[-1]})*2 + S({-1: -4}) - S({-1: 6 - 4*L})) and not total[('C_A', 'A')].c.get(-1))

print('4  r0 = finite part of J1 + J2 + S + the four daughter self-energies, with x a symbol')
# J2 is J1 at 1 - x with the two quark helicity structures interchanged; the sum above is symmetric under x <-> 1 - x.
finite = lambda col: S({0: total[(col, 'A')].c[0]})
r0_F = finite('C_F')*2 + S({0: -2*Z2}) + S({0: 2*(-2 - L)})       # soft: -zeta2 Delta_C; self-energies: 2 C_F (-2 - L)
r0_A = finite('C_A')*2 + S({0: Z2})
print(f'       C_F: {show(r0_F.c[0])}')
print(f'       C_A: {show(r0_A.c[0])}')
require('Eq. (8), C_F part: 12 - 6 zeta2 + 2 ln^2(x/(1-x))', zero(r0_F - S({0: 12 - 6*Z2 + 2*(LX - LB)*(LX - LB)})))
require('Eq. (8), C_A part: 2 + zeta2 - 2 ln^2 x - 2 ln^2(1-x)', zero(r0_A - S({0: 2 + Z2 - 2*LX*LX - 2*LB*LB})))
finish()
