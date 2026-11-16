"""Parametric integrands of the one-loop collinear function J1 for g -> gg (and of its tree), from diagrams
given as Lorentz tensors (vec_engine.py) and light-cone denominators (gg_diagrams.py).

Steps per diagram and z-interval: contract the tensor numerator; integrate over K^- (u) and l^- (w) by residues
(lc_engine.lc_integrate_multi, int dv/(2 pi)); integrate over l_perp in n = 2 - 2 eps dimensions with Feynman
parameters (transverse_scalar); read off the four structures
    S1 = (e.k)(i.j),  S2 = (e.i)(j.k),  S3 = (e.j)(i.k),  T4 = (e.k)(i.k)(j.k)/k^2     (Euclidean products),
at k^2 = 1.  Normalisation: the result is divided by the tree-level coefficient of S1 (computed with the same
rules), so that, as for g -> q qbar, the one-loop coefficient of structure s is
    a_s (mu^2 e^gamma/k^2)^eps int dz [eik] int dF I_s(z, F)   (in units of the tree coefficient of S1),
the (4 pi)^(-n/2), 1/2 and 1/(2 pi) of the loop measure being absorbed in a_s (mu^2 e^gamma/k^2)^eps.
"""
import itertools
import sympy as sp
import lc_engine as E
from vec_engine import (eps, kk, ll, lk, ek, el, ik, il, jk, jl, ei, ej, ij, tdot, contract)

x, z, u, w = sp.symbols('x z u w')
STRUCT = ('S1', 'S2', 'S3', 'T4')
_LDOT = {el: 'e', il: 'i', jl: 'j'}            # l . a for the external directions a


def _pairings(items):
    if not items:
        yield []
        return
    a = items[0]
    for n in range(1, len(items)):
        rest = items[1:n] + items[n + 1:]
        for m in _pairings(rest):
            yield [(a, items[n])] + m


def structures(expr):
    """polynomial in ek, ik, jk, ei, ej, ij (kk = 1) -> {S1: c, S2: c, S3: c, T4: c}"""
    p = sp.Poly(sp.expand(expr), ek, ik, jk, ei, ej, ij)
    out = {s: 0 for s in STRUCT}
    for mon, c in p.terms():
        if c == 0:
            continue
        dmon = dict(zip(('ek', 'ik', 'jk', 'ei', 'ej', 'ij'), mon))
        key = tuple(sorted((k, v) for k, v in dmon.items() if v))
        if key == (('ek', 1), ('ij', 1)):
            out['S1'] += c
        elif key == (('ei', 1), ('jk', 1)):
            out['S2'] += c
        elif key == (('ej', 1), ('ik', 1)):
            out['S3'] += c
        elif key == (('ek', 1), ('ik', 1), ('jk', 1)):
            out['T4'] += c
        else:
            raise ValueError(f'unexpected structure {key}')
    return out


def transverse_scalar(num, dens, test):
    """int d^n l/(2 pi)^n num/prod(dens), n = 2 - 2 eps, num polynomial in the l-dependent dots
    (ll, lk, el, il, jl); dens linear in (ll, lk, kk), each of one sign at the test point (flipped to positive).
    Returns ({structure: integrand over the Feynman simplex at kk = 1}, params)."""
    m = len(dens)
    params = sp.symbols(f'a1:{m}')
    alphas = list(params) + [1 - sum(params)]
    signs = []
    tv = {k_: v for k_, v in test.items() if k_ not in (ll, lk, kk)}
    for dd in dens:
        dd = sp.expand(sp.cancel(sp.together(dd)))
        A_, B_, C_ = dd.coeff(ll), dd.coeff(lk), dd.coeff(kk)
        An, Bn, Cn = (sp.N(c.subs(tv)) for c in (A_, B_, C_))
        if An != 0:
            assert Bn**2 - 4*An*Cn < 1e-12*abs(An*Cn) + 1e-30 and An*Cn >= 0, ('not definite', dd)
            signs.append(1 if An > 0 else -1)
        else:
            assert Bn == 0, ('linear in lk', dd)
            signs.append(1 if Cn > 0 else -1)
    dens = [sg*sp.expand(sp.cancel(sp.together(dd))) for sg, dd in zip(signs, dens)]
    sgn_total = 1
    for sg in signs:
        sgn_total *= sg
    comb = sp.expand(sum(a*dd for a, dd in zip(alphas, dens)))
    A = sp.expand(comb.coeff(ll))
    B = sp.expand(comb.coeff(lk))
    C = sp.expand(comb.coeff(kk))
    assert sp.expand(comb - A*ll - B*lk - C*kk) == 0
    beta = sp.cancel(-B/(2*A))
    M = sp.cancel(C - B**2/(4*A))
    n = 2 - 2*eps
    lp2, lpk, elp, ilp, jlp = sp.symbols('lp2 lpk elp ilp jlp')
    sub = {ll: lp2 + 2*beta*lpk + beta**2*kk, lk: lpk + beta*kk, el: elp + beta*ek, il: ilp + beta*ik,
           jl: jlp + beta*jk}
    cs = sp.expand(sp.sympify(num).subs(sub, simultaneous=True))
    poly = sp.Poly(cs, lp2, lpk, elp, ilp, jlp)
    out = 0
    for (p2, pk, pe, pi, pj), cc in poly.terms():
        items = ['k']*pk + ['e']*pe + ['i']*pi + ['j']*pj
        nl = len(items)
        if nl % 2:
            continue
        P = p2 + nl//2
        denom = sp.Integer(1)
        for jj in range(nl//2):
            denom *= (n + 2*jj)
        ang = 0
        for match in _pairings(items):
            pr = sp.Integer(1)
            for a_, b_ in match:
                pr *= tdot(a_, b_)
            ang += pr
        rad = (4*sp.pi)**(-n/2)*sp.gamma(P + n/2)*sp.gamma(m - P - n/2)/sp.gamma(n/2) \
            * A**(-P - n/2)*M**(P + n/2 - m)
        out += sgn_total*cc*ang/denom*rad
    res = structures(sp.expand(out).subs(kk, 1))
    return res, params


def residues(num, dens, at, order=(u, w)):
    """double residue integral of num/prod(dens); returns {tuple(dens): numerator}"""
    out = {}
    d0 = [(dd, 1) for dd in dens]
    for n1, d1 in E.lc_integrate_multi({'N': num}, d0, order[0], at):
        if len(order) == 1:
            key = tuple(sp.factor(e) for e, s in d1)
            out[key] = sp.cancel(out.get(key, 0) + n1['N'])
            continue
        for n2, d2 in E.lc_integrate_multi(n1, d1, order[1], at):
            key = tuple(sp.factor(e) for e, s in d2)
            out[key] = sp.cancel(out.get(key, 0) + n2['N'])
    return out


def tree_coefficients(diag, xval):
    """the tree: numerator, K^- residue; returns {structure: coefficient} at kk = 1 (rational in x)"""
    num = sp.expand(contract(diag['numerator']).subs(x, xval))
    dens = [sp.expand(dd.subs(x, xval)) for dd in diag['dens']]
    at = {kk: 1, ll: sp.Rational(7, 3), lk: sp.Rational(1, 5)}
    rr = residues(num, dens, at, order=(u,))
    tot = 0
    for dl, nm in rr.items():
        tot += nm/sp.Mul(*dl)
    tot = sp.expand(sp.cancel(tot*diag['factor']).subs(kk, 1))
    return structures(tot)


def diagram_integrands(diag, xval, intervals, norm, order=(u, w)):
    """returns [(z-interval, [(params, {structure: I})])]; norm = tree coefficient of S1 (with its factors)"""
    num = sp.expand(contract(diag['numerator']).subs(x, xval))
    dens = [sp.expand(dd.subs(x, xval)) for dd in diag['dens']]
    res = []
    strip = (4*sp.pi)**((2 - 2*eps)/2)
    for (za, zb) in intervals:
        zt = (sp.Rational(za) + sp.Rational(zb))/2
        at = {z: zt, kk: 1, ll: sp.Rational(7, 3), lk: sp.Rational(1, 5)}
        rr = residues(num, dens, at, order)
        pieces = []
        for dl, nm in rr.items():
            dep = [dd for dd in dl if sp.diff(dd, ll) != 0 or sp.diff(dd, lk) != 0]
            ind = [dd for dd in dl if dd not in dep]
            fac = 1/sp.Mul(*ind) if ind else 1
            if not dep:
                raise ValueError('no transverse denominator')
            r, prm = transverse_scalar(sp.expand(nm*fac), dep, at)
            out = {s: (sp.sympify(r[s])*strip*diag['factor']/norm if r[s] != 0 else 0) for s in STRUCT}
            pieces.append((prm, out))
        res.append(((za, zb), pieces))
    return res
