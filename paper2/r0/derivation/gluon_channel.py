"""The one-loop diagrams of J1(g -> gg) as Laurent series in eta and eps.

For each diagram: numerator by exact tensor contraction, residues in the minus components of the two loop momenta,
transverse integral in 2 - 2 eps dimensions, parameter integrals. The transverse denominators are l_perp^2 or
r_perp^2 = (k - l_perp)^2 and
    D = (l_perp - zeta k)^2 + zeta (1-zeta) k^2/(1-x),
with zeta the plus momentum of the on-shell loop gluon in units of x. After the residues,
    triangles:                                          zeta/(k^2 r_perp^2 D)
    four-gluon vertex bubble on a daughter:             -1/(x k^2 D)
    Wilson line coupling to the same daughter (W1k):    -1/(x k^2 D), eikonal denominator x zeta
    Wilson line coupling to the parent gluon (W1g):     -zeta/((1-x) l_perp^2 D), eikonal denominator x zeta
    Wilson line coupling to the other daughter (W1kb):  -zeta/((1-x) r_perp^2 D), eikonal factor [x (1-zeta)]^(-1-eta)
Only W1kb and its counterpart W2k have a rapidity divergence.

Units k^2 = mu^2 = nu_r = 1. The momentum fraction x is a rational number or the symbol X of algebra.py."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from fractions import Fraction as Fr
from math import comb
from algebra import C, V_, contract, Ser, coef, Poly, RatX
from gluon_diagrams import build, on_pole, to_r
from integrals import (average, integrate_a, integrate_b, add, scale, CGAMMA, FOURIER, STRUCT, evaluate, RT, rt,
                       LB, LX, LI, Z2, ParameterIntegrals, ser_structures, ser_exp)

_NUM = {}


def numerator(name, x):
    key = (name, coef(x))
    if key not in _NUM: _NUM[key] = contract(build(x)[name]['num'])
    return _NUM[key]


def split_ll(N):
    """N = N0 + l_perp^2 N1"""
    N0 = C(0); N1 = C(0)
    for (p,), c in N.collect(['ll']).items():
        if p == 0: N0 = N0 + c
        else: N1 = N1 + c.shift('ll', p - 1)
    return N0, N1


def unit(name, x):
    """factor of the diagram divided by the tree coefficient of S1, which is -2i: (rational number, power of i)"""
    f, p = build(x)[name]['factor']
    return f*Fr(-1, 2), (p - 1) % 4


def region(name, x):
    """the region of the loop plus momentum in which the residues are taken as listed above, in momentum space and
    without the power of i: {helicity structure: series in eps}"""
    x = coef(x); xb = 1 - x; ze, s = V_('ze'), V_('s')
    kind = build(x)[name]['kind']; num = numerator(name, x)
    if kind == 'buba':
        N = on_pole(num, 'l', x)
        res = integrate_a(average(N, ze, 1), x, 0, -1)
    elif kind == 'w1k':
        N = on_pole(num, 'l', x)
        res = integrate_a(average(N, ze, 1), x, -1, -1/x)
    elif kind == 'tri':
        N0, N1 = split_ll(on_pole(num, 'r', x))
        res = add(integrate_a(average(N1, ze, 1), x, 1, x), integrate_b(average(N0, s, 2), x, 0, x))
    elif kind == 'w1g':
        N0, N1 = split_ll(on_pole(num, 'l', x))
        res = add(integrate_a(average(N1, ze, 1), x, 0, -1/xb), integrate_b(average(N0, s, 2), x, -1, -1/xb))
    else:
        raise ValueError(kind)
    f, p = unit(name, x)
    return {st: v*CGAMMA*f for st, v in res.items()}


def relabel(p):
    """x -> 1 - x in a coefficient: ln(1-x) <-> ln x, Li2(x) -> zeta2 - ln x ln(1-x) - Li2(x)"""
    if not hasattr(p, 't'): return p
    out = rt(0)
    for (a, b, c, d, e1, e2), v in p.t.items():
        out = out + (LX**a)*(LB**b)*((Z2 - LX*LB - LI)**c)*(Z2**d)*v
        assert e1 == 0 and e2 == 0
    return out


def at_one_minus_x(c):
    """a coefficient written at x, evaluated at 1 - x"""
    if isinstance(c, RatX): return c.reflect()
    if not isinstance(c, Poly): return c
    return relabel(Poly(c.R, {k: (v.reflect() if isinstance(v, RatX) else v) for k, v in c.t.items()}))


def exchange_daughters(d):
    """a result computed at 1 - x, rewritten at x, with S2 <-> S3"""
    m = {'S1': 'S1', 'S2': 'S3', 'S3': 'S2', 'T4': 'T4'}
    return {m[s]: Ser({n: relabel(c) for n, c in v.c.items()}) for s, v in d.items()}


# a triangle has two regions; the second is the first region of the diagram with the two daughters exchanged
OTHER_ORIENTATION = {'Vtri': ('Vtri', 1), 'Vgh1': ('Vgh2', -1), 'Vgh2': ('Vgh1', -1), 'Vq1': ('Vq2', -1), 'Vq2': ('Vq1', -1)}


def triangle(name, x):
    other, sign = OTHER_ORIENTATION[name]
    return add(region(name, x), scale(exchange_daughters(region(other, 1 - coef(x))), sign))


def other_daughter(name, x):
    """W1kb: the gluon from the Wilson line of the first operator attaches to the daughter created by the second.
    r = K - l is on shell with r^+ = x zeta. The parameter integrals factorise,
        Int_0^1 dzeta (1-zeta)^(-1-eta) Int_0^zeta ds s^(b-1-eps) V^(-1-eps)
          = B(-eps-eta, b-eps) Int_0^1 drho rho^(-1-eta) [(1-x) + x rho]^(-1-eps)
          = (1-x)^(-1-eps) B(-eps, b-eps) [-1/eta + psi(-eps) - psi(b-2eps) + ln(1-x) - eps Li2(x)] + O(eta, eps^2).
    Returns (coefficient of 1/eta, coefficient of eta^0), each {helicity structure: series in eps}, in momentum space."""
    x = coef(x); xb = 1 - x; ze, s = V_('ze'), V_('s'); M = ParameterIntegrals(x)
    N0, N1 = split_ll(to_r(numerator(name, x).subs('z', x*(1 - ze))))
    assert not N1 and numerator(name, x).degree('u') <= 0 and numerator(name, x).degree('w') <= 0
    av = average(N0, s, 2); assert set(av) == {-1}
    g = Ser({0: LB, 1: -LI})
    psi = {0: Ser({-1: Fr(1, 2), 1: Z2}), 1: Ser({-1: Fr(1), 1: Z2})}           # psi(-eps) - psi(b - 2 eps)
    pole, fin = {}, {}
    for st, sr in ser_structures(av[-1]).items():
        for n, poly in sr.c.items():
            for (a, b), c in poly.collect(['ze', 's']).items():
                cv = c.const_value()
                Bb = M.beta(-1, b - 1)*ser_exp(-LB)*(1/xb)                       # (1-x)^(-1-eps) B(-eps, b-eps)
                p = Bb*(-1)
                f = Bb*(psi[b] + g + Ser({0: LX}))                               # x^(-eta) gives + ln x
                for j in range(1, a + 1):                                        # zeta^a = sum_j C(a,j) (-1)^j (1-zeta)^j
                    for i in range(j):
                        f = f + M.M1(i, b)*(comb(a, j)*(-1)**j*comb(j - 1, i)*(-1)**i)
                o = ser_exp(LB)*(-cv)                                            # -(1/(1-x)) Gamma(1+eps) [s V/(1-x)]^(-1-eps)
                pole[st] = pole.get(st, Ser()) + (p*o).shift(n)
                fin[st] = fin.get(st, Ser()) + (f*o).shift(n)
    f_, p_ = unit(name, x)
    return ({st: v*CGAMMA*f_ for st, v in pole.items()}, {st: v*CGAMMA*f_ for st, v in fin.items()})


SIGN = {0: 1, 2: -1, 3: 1}        # i^p for p = 0, 2; for the quark triangles (p = 3) the factor -i is kept with the colour factor

# diagrams obtained from another one by exchanging the two daughters, and the ratio of the two stored colour factors
EXCHANGED = {'W2g': ('W1g', 1), 'W2kb': ('W1k', 1), 'W2k': ('W1kb', -1),
             'Vbubb1': ('Vbuba1', -1), 'Vbubb2': ('Vbuba2', -1), 'Vbubb3': ('Vbuba3', -1)}


def diagram(name, x):
    """a diagram as stored in paper2/r0/theory/gg_laurent_x*.json: divided by the tree coefficient of S1, in impact
    parameter space at b = b0. Returns (coefficient of 1/eta, coefficient of eta^0), each {helicity structure: series
    in eps}. The coefficient of T4 is that of the momentum space structure; its transform is taken when the two
    splittings are contracted. For the quark triangles the factor -i is left out."""
    x = coef(x)
    if name in EXCHANGED:
        source, sign = EXCHANGED[name]
        p, f = diagram(source, 1 - x)
        return scale(exchange_daughters(p), sign), scale(exchange_daughters(f), sign)
    kind = build(x)[name]['kind']
    if kind in ('w1v4', 'w2v4'):                  # the gluon from the Wilson line enters a four-gluon vertex
        assert not numerator(name, x)
        return {}, {}
    if kind == 'w1kb':
        pole, fin = other_daughter(name, x)
    elif kind == 'tri':
        pole, fin = {}, triangle(name, x)
    else:
        pole, fin = {}, region(name, x)
    sg = SIGN[unit(name, x)[1]]
    to_b = lambda d: {st: v*(FOURIER if st != 'T4' else Ser({0: rt(1)}))*sg for st, v in d.items()}
    return to_b(pole), to_b(fin)
