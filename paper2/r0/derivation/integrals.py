"""Transverse integration, tensor reduction and parameter integrals for the diagrams of J1(g -> gg).

Transverse integral in n = 2 - 2 eps dimensions, after the shift l = l' + beta k:
    Int d^n l'/pi^(n/2)  (l'^2)^p (l'.a_1)...(l'.a_2m) / (l'^2 + Delta)^M
      = [sum over pairings of the a's] 2^(-m) prod_{j=m}^{P-1} (1 - eps + j)  Gamma(M-P-1+eps)/Gamma(M)  Delta^(P+1-M-eps),   P = p + m.
Parameter integrals over 0 < s < zeta < 1, with V = 1 - x zeta - (1-x) s. The results are Laurent series in eps whose
coefficients are polynomials in lb = ln(1-x), lx = ln x, Li = Li2(x) and z2 = zeta2. Their coefficients are rational
numbers, or rational functions of x when x is the symbol X of algebra.py."""
from fractions import Fraction as Fr
from math import comb
import itertools
from algebra import Ring, Poly, Ser, ser_lin, ser_inv_lin, ser_prod, gamma_ratio, ser_exp, C, V_, RN, tdot, EMAX, RatX, coef

RT = Ring(['lb', 'lx', 'Li', 'z2', 'DIV', 'DIV2'])
LB, LX, LI, Z2, DIV, DIV2 = (RT.var(n) for n in RT.names)
rt = RT.const


# ----------------------------------------------------------------------------- transverse integral and tensor reduction
def _pairings(items):
    if not items:
        yield []
        return
    a = items[0]
    for n in range(1, len(items)):
        for m in _pairings(items[1:n] + items[n + 1:]):
            yield [(a, items[n])] + m


def poly_to_ser(poly):
    """polynomial in eps and other variables -> series in eps with polynomial coefficients"""
    return Ser({p: c for (p,), c in poly.collect(['e']).items()})


def average(N, beta, m):
    """N is a polynomial in the scalar products of l. Shift l = l' + beta k and integrate over l' with (l'^2 + Delta)^-m,
    m = 1 or 2. Returns {q: series with polynomial coefficients}, to be multiplied by Gamma(1+eps) Delta^(q - eps)."""
    kk = V_('kk')
    n = N.subs('ll', V_('p2') + 2*beta*V_('pk') + beta*beta*kk).subs('lk', V_('pk') + beta*kk)
    n = n.subs('el', V_('pe') + beta*V_('ek')).subs('il', V_('pi') + beta*V_('ik')).subs('jl', V_('pj') + beta*V_('jk'))
    out = {}
    for (p2, pk, pe, pi, pj), c in n.collect(['p2', 'pk', 'pe', 'pi', 'pj']).items():
        items = ['k']*pk + ['e']*pe + ['i']*pi + ['j']*pj
        nl = len(items)
        if nl % 2: continue
        P = p2 + nl//2
        ang = C(0)
        for match in _pairings(items):
            pr = C(1)
            for a_, b_ in match: pr = pr*tdot(a_, b_)
            ang = ang + pr
        fac = ser_prod([ser_lin(1 + j, -1) for j in range(nl//2, P)])*Fr(1, 2**(nl//2))
        Ng = m - P - 1                                    # Gamma(Ng + eps)/Gamma(1 + eps)
        fac = fac*gamma_ratio(Ng, 1)
        q = P + 1 - m
        term = poly_to_ser(c*ang)*fac
        out[q] = out[q] + term if q in out else term
    return out


STRUCT = ('S1', 'S2', 'S3', 'T4')
_SKEY = {(('ek', 1), ('ij', 1)): 'S1', (('ei', 1), ('jk', 1)): 'S2', (('ej', 1), ('ik', 1)): 'S3', (('ek', 1), ('ik', 1), ('jk', 1)): 'T4'}


def structures(poly):
    """polynomial in the scalar products of k, e, e_1, e_2 (with k^2 = 1) -> {helicity structure: polynomial in the other variables}"""
    names = ('ek', 'ik', 'jk', 'ei', 'ej', 'ij')
    out = {}
    for mon, c in poly.subs('kk', 1).collect(names).items():
        key = tuple(sorted((k, v) for k, v in zip(names, mon) if v))
        if key not in _SKEY: raise ValueError(f'not one of S1, S2, S3, T4: {key}')
        s = _SKEY[key]
        out[s] = out[s] + c if s in out else c
    return out


def ser_structures(ser):
    """series with polynomial coefficients -> {helicity structure: series}"""
    out = {}
    for n, c in ser.c.items():
        for s, p in structures(c).items():
            out.setdefault(s, {})[n] = p
    return {s: Ser(d) for s, d in out.items()}


# ----------------------------------------------------------------------------- Laurent polynomials in zeta
def pmul(a, b):
    out = {}
    for i, u in a.items():
        for j, v in b.items():
            out[i + j] = out.get(i + j, 0) + u*v
    return {k: v for k, v in out.items() if v}


def ppow(a, n):
    out = {0: Fr(1)}
    for _ in range(n): out = pmul(out, a)
    return out


def H(n): return sum(Fr(1, j) for j in range(1, n + 1))


class ParameterIntegrals:
    """integrals over 0 < s < zeta < 1 at a given x"""
    def __init__(self, x):
        self.x = coef(x); self.xb = 1 - self.x
        self.c = {0: Fr(1), 1: -self.x}              # c = 1 - x zeta
        self.y = {0: Fr(1), 1: Fr(-1)}               # 1 - zeta
        self._m0 = {}; self._m1 = {}

    def R(self, m):
        """Int_0^1 dzeta zeta^m/(1 - x zeta), m >= 0"""
        x = self.x
        return (LB*(-1/x) - sum(x**j*Fr(1, j + 1) for j in range(m)))*x**(-m)

    def i_poly(self, P):
        tot = rt(0)
        for n, v in P.items():
            tot = tot + (DIV*v if n == -1 else rt(v*Fr(1, n + 1)))
            assert n >= -1
        return tot

    def i_lnz(self, P):
        tot = rt(0)
        for n, v in P.items():
            assert n >= -1
            tot = tot + (DIV2*v if n == -1 else rt(-v*Fr(1, (n + 1)**2)))
        return tot

    def i_lnc(self, P):
        """Int_0^1 dzeta P(zeta) ln(1 - x zeta)"""
        tot = rt(0); x = self.x
        for n, v in P.items():
            assert n >= -1
            if n == -1: tot = tot - LI*v
            else: tot = tot + (LB*Fr(1, n + 1) + self.R(n + 1)*(x*Fr(1, n + 1)))*v
        return tot

    def i_ln1(self, P):
        """Int_0^1 dzeta P(zeta) ln(1 - zeta)"""
        tot = rt(0)
        for n, v in P.items():
            assert n >= -1
            tot = tot + (Z2*(-v) if n == -1 else rt(-v*H(n + 1)*Fr(1, n + 1)))
        return tot

    def M0(self, a, b):
        """Int_0^1 dzeta zeta^a Int_0^zeta ds s^b (s V)^-eps  through eps^1;  a >= -1, b >= 0"""
        key = (a, b)
        if key in self._m0: return self._m0[key]
        za = {a: Fr(1)}; xb = self.xb
        inner0 = {b + 1: Fr(1, b + 1)}
        I00 = self.i_poly(pmul(za, inner0))
        Ils = self.i_lnz(pmul(za, inner0)) - self.i_poly(pmul(za, {b + 1: Fr(1, (b + 1)**2)}))
        cb1 = pmul(za, ppow(self.c, b + 1))
        S1 = sum(comb(b, j)*(-1)**j*Fr(1, j + 1) for j in range(b + 1))
        S2 = sum(comb(b, j)*(-1)**j*Fr(1, (j + 1)**2) for j in range(b + 1))
        tot = self.i_lnc(cb1)*S1 - self.i_poly(cb1)*S2
        for j in range(b + 1):
            co = comb(b, j)*(-1)**j
            pj = pmul(za, pmul(ppow(self.c, b - j), ppow(self.y, j + 1)))
            tot = tot - self.i_ln1(pj)*(co*Fr(1, j + 1)) + self.i_poly(pj)*(co*Fr(1, (j + 1)**2))
        IlV = tot*xb**(-(b + 1))
        res = Ser({0: I00, 1: -(Ils + IlV)})
        self._m0[key] = res
        return res

    def M1(self, a, b):
        """Int_0^1 dzeta zeta^a Int_0^zeta ds s^(b-1-eps) V^(-1-eps)  through eps^0;  a >= -1, b >= 0"""
        key = (a, b)
        if key in self._m1: return self._m1[key]
        za = {a: Fr(1)}; x, xb = self.x, self.xb
        if b >= 1:
            cb = pmul(za, ppow(self.c, b - 1))
            tot = self.i_lnc(cb) - self.i_ln1(cb)
            for j in range(1, b):
                co = comb(b - 1, j)*(-1)**j*Fr(1, j)
                tot = tot + (self.i_poly(cb) - self.i_poly(pmul(za, pmul(ppow(self.c, b - 1 - j), ppow(self.y, j)))))*co
            res = Ser({0: tot*xb**(-b)})
        elif a >= 0:
            Q = {j: x**(j - a) for j in range(a)}
            fin = LB*LB*(-x**(-a)/(2*x)) - (self.i_lnz(Q) + self.i_lnc(Q)*2 - self.i_ln1(Q))
            res = Ser({-1: -self.R(a), 0: fin})
        else:
            # a = -1, b = 0: 1/eps^2 + ln(1-x)/eps + zeta2 - 2 Li2(x) - ln^2(1-x)/2
            res = Ser({-2: rt(1), -1: LB, 0: Z2 - LI*2 - LB*LB*Fr(1, 2)})
        self._m1[key] = res
        return res

    def beta(self, A, B, cA=1, cB=1):
        """Int_0^1 dzeta zeta^(A - cA eps) (1-zeta)^(B - cB eps) = Gamma(A+1-cA eps) Gamma(B+1-cB eps)/Gamma(A+B+2-(cA+cB) eps)"""
        num = gamma_ratio(A + 1, -cA)*gamma_ratio(B + 1, -cB)
        N3 = A + B + 2; c3 = -(cA + cB)
        if N3 >= 1: den = ser_prod([ser_inv_lin(j, c3) for j in range(1, N3)])
        else: den = ser_prod([ser_lin(j, c3) for j in range(N3, 1)])
        # Gamma(1 - cA eps) Gamma(1 - cB eps)/Gamma(1 - (cA+cB) eps) = 1 - cA cB zeta2 eps^2 + O(eps^3)
        g = Ser({0: rt(1), 2: Z2*(-cA*cB)})
        return num*den*g


CGAMMA = Ser({0: rt(1), 2: Z2*Fr(1, 2)})          # e^(eps gamma_E) Gamma(1+eps), from the transverse integral
FOURIER = Ser({0: rt(1), 2: Z2})                  # e^(-eps gamma_E) F_eps, from the transform to impact parameter space at b = b0


def integrate_a(av, x, zshift=0, pref=1):
    """one transverse denominator: av = {q: series with coefficients polynomial in zeta} multiplies
    Gamma(1+eps) [zeta(1-zeta) k^2/(1-x)]^(q-eps). Returns {helicity structure: series} for pref Int_0^1 dzeta zeta^zshift ..."""
    M = ParameterIntegrals(x); xb = 1 - coef(x); out = {}
    for q, ser in av.items():
        for s, sr in ser_structures(ser).items():
            for n, poly in sr.c.items():
                for (a,), c in poly.collect(['ze']).items():
                    val = M.beta(a + zshift + q, q)*ser_exp(LB)*(xb**(-q)*c.const_value()*coef(pref))
                    out[s] = out[s] + val.shift(n) if s in out else val.shift(n)
    return out


def integrate_b(av, x, zshift=0, pref=1):
    """two transverse denominators: av = {q: series with coefficients polynomial in zeta and s} multiplies
    Gamma(1+eps) [s V k^2/(1-x)]^(q-eps). Returns {helicity structure: series} for pref Int_0^1 dzeta zeta^zshift Int_0^zeta ds ..."""
    M = ParameterIntegrals(x); x = coef(x); xb = 1 - x; out = {}
    Vp = {(0, 0): Fr(1), (1, 0): -x, (0, 1): -xb}          # V as a polynomial in (zeta, s)
    for q, ser in av.items():
        if q >= 0:
            sv = {(0, 0): Fr(1)}
            for _ in range(q):                              # (s V)^q
                new = {}
                for (i, j), v in sv.items():
                    for (a, b), w_ in Vp.items():
                        key = (i + a, j + b + 1); new[key] = new.get(key, 0) + v*w_
                sv = new
        for s, sr in ser_structures(ser).items():
            for n, poly in sr.c.items():
                for (a, b), c in poly.collect(['ze', 's']).items():
                    cv = c.const_value()*coef(pref)*xb**(-q)
                    if q == -1:
                        val = M.M1(a + zshift, b)
                    else:
                        val = Ser()
                        for (i, j), v in sv.items():
                            val = val + M.M0(a + zshift + i, b + j)*v
                    val = val*ser_exp(LB)*cv
                    out[s] = out[s] + val.shift(n) if s in out else val.shift(n)
    return out


def add(d1, d2):
    out = dict(d1)
    for s, v in d2.items(): out[s] = out[s] + v if s in out else v
    return out


def scale(d, c): return {s: v*c for s, v in d.items()}


def evaluate(ser, x):
    """Laurent coefficients {order in eps: number} at a number x"""
    import math
    from scipy.special import spence
    x = float(x)
    vals = {'x': x, 'lb': math.log(1 - x), 'lx': math.log(x), 'Li': float(spence(1 - x)), 'z2': math.pi**2/6}
    num = lambda c: float(c.ev(vals)) if isinstance(c, Poly) else float(c.at(x)) if isinstance(c, RatX) else float(c)
    return {n: num(c) for n, c in ser.c.items()}


_NAMES = {'lb': ('ln(1-x)', 'ln^{}(1-x)'), 'lx': ('ln x', 'ln^{} x'), 'Li': ('Li2(x)', 'Li2(x)^{}'), 'z2': ('zeta2', 'zeta2^{}')}


def show(c):
    """a coefficient of a series as text"""
    if not isinstance(c, Poly): return str(c)
    parts = []
    for k, v in sorted(c.t.items(), key=lambda kv: (sum(kv[0]), kv[0])):
        mon = ' '.join(_NAMES[nm][0] if k[c.R.idx[nm]] == 1 else _NAMES[nm][1].format(k[c.R.idx[nm]]) for nm in ('lx', 'lb', 'Li', 'z2') if k[c.R.idx[nm]])
        t = str(v)
        neg = t.startswith('-') and not (' + ' in t or ' - ' in t[1:])
        if neg: t = t[1:]
        if mon:
            if ' + ' in t or ' - ' in t or '/' in t: t = f'({t})'
            t = mon if t == '1' else f'{t} {mon}'
        parts.append(('-' if neg else '+', t))
    if not parts: return '0'
    return ('-' if parts[0][0] == '-' else '') + parts[0][1] + ''.join(f' {sg} {t}' for sg, t in parts[1:])


def show_series(ser, orders=(-2, -1, 0)):
    """the coefficients of a series in eps as text"""
    return ';  '.join(f'eps^{n}: {show(ser.c[n])}' for n in orders if ser.c.get(n)) or '0'
