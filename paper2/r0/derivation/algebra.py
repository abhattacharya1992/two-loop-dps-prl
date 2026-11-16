"""Exact algebra for the collinear function J1(g -> gg).

Polynomials in several variables (negative powers allowed), truncated Laurent series in eps, and Lorentz tensors
in light-cone components. Coefficients are rational numbers or rational functions of the momentum fraction x with
denominators x^a (1-x)^b, so that x can be kept as a symbol. Only the standard library is used."""
from fractions import Fraction as Fr
import itertools


def _pmul(p, q):
    if not p or not q: return []
    out = [Fr(0)]*(len(p) + len(q) - 1)
    for i, u in enumerate(p):
        if u:
            for j, v in enumerate(q): out[i + j] += u*v
    return out


def _xpow(m, n):
    """coefficients of x^m (1-x)^n"""
    p = [Fr(1)]
    for _ in range(n): p = _pmul(p, [Fr(1), Fr(-1)])
    return [Fr(0)]*m + p


class RatX:
    """P(x)/(x^a (1-x)^b), P a polynomial with rational coefficients; common factors x and 1-x are cancelled,
    which makes the representation unique"""
    __slots__ = ('p', 'a', 'b')

    def __init__(self, p, a=0, b=0):
        p = [Fr(c) for c in p]
        while p and p[-1] == 0: p.pop()
        if not p:
            self.p, self.a, self.b = (), 0, 0
            return
        while a > 0 and p[0] == 0:
            p.pop(0); a -= 1
        while b > 0 and sum(p) == 0:                 # P(1) = 0: P = (1-x) Q with Q_k = p_0 + ... + p_k
            q, acc = [], Fr(0)
            for c in p[:-1]:
                acc += c; q.append(acc)
            p = q; b -= 1
        self.p, self.a, self.b = tuple(p), a, b

    @staticmethod
    def co(o):
        if isinstance(o, RatX): return o
        if isinstance(o, (int, Fr)): return RatX([o])
        return None

    def __add__(self, o):
        o = RatX.co(o)
        if o is None: return NotImplemented
        A, B = max(self.a, o.a), max(self.b, o.b)
        p1 = _pmul(list(self.p), _xpow(A - self.a, B - self.b)); p2 = _pmul(list(o.p), _xpow(A - o.a, B - o.b))
        n = max(len(p1), len(p2))
        return RatX([(p1[i] if i < len(p1) else 0) + (p2[i] if i < len(p2) else 0) for i in range(n)], A, B)
    __radd__ = __add__

    def __neg__(self): return RatX([-c for c in self.p], self.a, self.b)

    def __sub__(self, o):
        o = RatX.co(o)
        return NotImplemented if o is None else self + (-o)

    def __rsub__(self, o):
        o = RatX.co(o)
        return NotImplemented if o is None else o + (-self)

    def __mul__(self, o):
        o = RatX.co(o)
        if o is None: return NotImplemented
        return RatX(_pmul(list(self.p), list(o.p)), self.a + o.a, self.b + o.b)
    __rmul__ = __mul__

    def inv(self):
        """inverse; the numerator must be a constant times powers of x and 1-x"""
        p = list(self.p); m = n = 0
        while p and p[0] == 0:
            p.pop(0); m += 1
        while len(p) > 1:
            if sum(p) != 0: raise ZeroDivisionError('the inverse is not of the form P/(x^a (1-x)^b)')
            q, acc = [], Fr(0)
            for c in p[:-1]:
                acc += c; q.append(acc)
            p = q; n += 1
        return RatX(_pmul([1/p[0]], _xpow(self.a, self.b)), m, n)

    def __truediv__(self, o): return self*RatX.co(o).inv()
    def __rtruediv__(self, o): return RatX.co(o)*self.inv()

    def __pow__(self, n):
        base = self if n >= 0 else self.inv()
        out = RatX([1])
        for _ in range(abs(n)): out = out*base
        return out

    def __bool__(self): return bool(self.p)

    def __eq__(self, o):
        o = RatX.co(o)
        return o is not None and (self.p, self.a, self.b) == (o.p, o.a, o.b)

    def __hash__(self): return hash((self.p, self.a, self.b))

    def reflect(self):
        """the same function at 1 - x"""
        q = [Fr(0)]*len(self.p)
        for i, c in enumerate(self.p):
            for j, v in enumerate(_xpow(0, i)): q[j] += c*v
        return RatX(q, self.b, self.a)

    def at(self, x):
        """value at a number x"""
        return sum(c*x**i for i, c in enumerate(self.p))/(x**self.a*(1 - x)**self.b)

    def __str__(self):
        terms = []
        for i, c in enumerate(self.p):
            if c == 0: continue
            mon = '' if i == 0 else 'x' if i == 1 else f'x^{i}'
            mag = abs(c); num = '' if (mag == 1 and mon) else str(mag)
            terms.append(('-' if c < 0 else '+', (num + ' ' + mon).strip()))
        if not terms: return '0'
        s = ('-' if terms[0][0] == '-' else '') + terms[0][1] + ''.join(f' {sg} {t}' for sg, t in terms[1:])
        den = ' '.join(([('x' if self.a == 1 else f'x^{self.a}')] if self.a else []) + ([('(1-x)' if self.b == 1 else f'(1-x)^{self.b}')] if self.b else []))
        if den: s = (f'({s})' if len(terms) > 1 else s) + '/' + (f'({den})' if self.a and self.b else den)
        return s
    __repr__ = __str__


X = RatX([0, 1])          # the momentum fraction x as a symbol


def coef(c):
    """a coefficient: a rational function of x is kept, anything else becomes a Fraction"""
    return c if isinstance(c, RatX) else Fr(c)


class Ring:
    def __init__(self, names):
        self.names = list(names); self.idx = {n: i for i, n in enumerate(self.names)}; self.n = len(self.names)
        self.zero = (0,)*self.n

    def const(self, c):
        c = coef(c)
        return Poly(self, {self.zero: c} if c else {})

    def var(self, name, p=1):
        e = [0]*self.n; e[self.idx[name]] = p
        return Poly(self, {tuple(e): Fr(1)})


class Poly:
    __slots__ = ('R', 't')

    def __init__(self, R, t=None):
        self.R = R; self.t = t if t is not None else {}

    def _co(self, o):
        return o if isinstance(o, Poly) else self.R.const(o)

    def __bool__(self):
        return bool(self.t)

    def __add__(self, o):
        o = self._co(o); out = dict(self.t)
        for k, v in o.t.items():
            nv = out.get(k, 0) + v
            if nv: out[k] = nv
            else: out.pop(k, None)
        return Poly(self.R, out)
    __radd__ = __add__

    def __neg__(self):
        return Poly(self.R, {k: -v for k, v in self.t.items()})

    def __sub__(self, o):
        return self + (-self._co(o))

    def __rsub__(self, o):
        return self._co(o) - self

    def __mul__(self, o):
        if not isinstance(o, Poly):
            o = coef(o)
            return Poly(self.R, {k: v*o for k, v in self.t.items()} if o else {})
        a, b = (self.t, o.t) if len(self.t) <= len(o.t) else (o.t, self.t)
        out = {}
        for k1, v1 in a.items():
            for k2, v2 in b.items():
                k = tuple([p + q for p, q in zip(k1, k2)])
                out[k] = out.get(k, 0) + v1*v2
        return Poly(self.R, {k: v for k, v in out.items() if v})
    __rmul__ = __mul__

    def __truediv__(self, c):
        c = coef(c)
        return self*(c.inv() if isinstance(c, RatX) else 1/c)

    def __pow__(self, n):
        out = self.R.const(1)
        for _ in range(n): out = out*self
        return out

    def subs(self, name, val):
        """replace the variable name (non-negative powers only) by a polynomial or number"""
        i = self.R.idx[name]; val = self._co(val); groups = {}
        for k, v in self.t.items():
            assert k[i] >= 0, (name, k[i])
            groups.setdefault(k[i], {})[k[:i] + (0,) + k[i + 1:]] = v
        out = Poly(self.R); pw = self.R.const(1); p = 0
        for q in sorted(groups):
            while p < q: pw = pw*val; p += 1
            out = out + Poly(self.R, groups[q])*pw
        return out

    def shift(self, name, p):
        """multiply by name^p (p may be negative)"""
        i = self.R.idx[name]
        return Poly(self.R, {k[:i] + (k[i] + p,) + k[i + 1:]: v for k, v in self.t.items()})

    def collect(self, names):
        """{exponents of names: polynomial in the other variables}"""
        ii = [self.R.idx[n] for n in names]; out = {}
        for k, v in self.t.items():
            key = tuple(k[i] for i in ii); kk = list(k)
            for i in ii: kk[i] = 0
            out.setdefault(key, {})[tuple(kk)] = v
        return {key: Poly(self.R, t) for key, t in out.items()}

    def degree(self, name):
        i = self.R.idx[name]
        return max((k[i] for k in self.t), default=-99)

    def mindeg(self, name):
        i = self.R.idx[name]
        return min((k[i] for k in self.t), default=99)

    def ev(self, vals):
        """numerical value; vals: {name: number}; coefficients that depend on x are evaluated at vals['x']"""
        tot = 0
        for k, v in self.t.items():
            term = v.at(vals['x']) if isinstance(v, RatX) else v
            for i, p in enumerate(k):
                if p:
                    term = term*vals[self.R.names[i]]**p
            tot = tot + term
        return tot

    def const_value(self):
        assert all(k == self.R.zero for k in self.t), self
        return self.t.get(self.R.zero, Fr(0))

    def __repr__(self):
        if not self.t: return '0'
        parts = []
        for k, v in sorted(self.t.items()):
            mon = '*'.join(f'{self.R.names[i]}^{p}' if p != 1 else self.R.names[i] for i, p in enumerate(k) if p)
            parts.append(f'{v}' + ('*' + mon if mon else ''))
        return ' + '.join(parts)


# ----------------------------------------------------------------------------- truncated Laurent series in eps
EMAX = 2          # keep powers eps^n with n <= EMAX


class Ser:
    """sum_n c_n eps^n, n <= EMAX, coefficients in a commutative ring (Fraction or Poly)"""
    __slots__ = ('c',)

    def __init__(self, c=None):
        self.c = {n: v for n, v in (c or {}).items() if n <= EMAX and v}

    @staticmethod
    def one(): return Ser({0: Fr(1)})

    def __add__(self, o):
        if not isinstance(o, Ser): o = Ser({0: o})
        out = dict(self.c)
        for n, v in o.c.items():
            out[n] = out[n] + v if n in out else v
        return Ser(out)
    __radd__ = __add__

    def __neg__(self): return Ser({n: -v for n, v in self.c.items()})
    def __sub__(self, o): return self + (-(o if isinstance(o, Ser) else Ser({0: o})))

    def __mul__(self, o):
        if not isinstance(o, Ser):
            return Ser({n: v*o for n, v in self.c.items()})
        out = {}
        for n1, v1 in self.c.items():
            for n2, v2 in o.c.items():
                if n1 + n2 <= EMAX:
                    out[n1 + n2] = out[n1 + n2] + v1*v2 if n1 + n2 in out else v1*v2
        return Ser(out)
    __rmul__ = __mul__

    def shift(self, p): return Ser({n + p: v for n, v in self.c.items()})

    def get(self, n, default=0): return self.c.get(n, default)

    def __repr__(self): return ' + '.join(f'({v})*eps^{n}' for n, v in sorted(self.c.items())) or '0'


def ser_lin(a, b):
    """a + b eps"""
    return Ser({0: Fr(a), 1: Fr(b)})


def ser_inv_lin(a, b):
    """1/(a + b eps), a != 0, or 1/(b eps)"""
    a, b = Fr(a), Fr(b)
    if a == 0: return Ser({-1: 1/b})
    return Ser({n: (1/a)*(-b/a)**n for n in range(0, EMAX + 4)})


def ser_prod(factors):
    out = Ser.one()
    for f in factors: out = out*f
    return out


def gamma_ratio(N, c=1):
    """Gamma(N + c eps)/Gamma(1 + c eps) for integer N, as a series"""
    if N >= 1:
        return ser_prod([ser_lin(j, c) for j in range(1, N)])
    return ser_prod([ser_inv_lin(j, c) for j in range(N, 1)])


def ser_exp(coeff):
    """exp(eps*coeff), coeff in the coefficient ring"""
    out = {0: Fr(1)}; term = Fr(1)
    for n in range(1, EMAX + 4):
        term = term*coeff*Fr(1, n); out[n] = term
    return Ser(out)


# ----------------------------------------------------------------------------- Lorentz tensors in light-cone components
NUMVARS = ['u', 'w', 'z', 'ze', 's', 'e', 'kk', 'll', 'lk', 'ek', 'el', 'ik', 'il', 'jk', 'jl', 'ei', 'ej', 'ij', 'ee', 'ii', 'jj',
           'p2', 'pk', 'pe', 'pi', 'pj']
RN = Ring(NUMVARS)
TN = ('k', 'l', 'e', 'i', 'j')
_DOTNAME = {}
for a_, b_ in itertools.combinations_with_replacement(TN, 2):
    nm = a_ + b_ if a_ == b_ else ''.join(sorted((a_, b_), key=TN.index))
    nm = {'kl': 'lk', 'ke': 'ek', 'ki': 'ik', 'kj': 'jk', 'le': 'el', 'li': 'il', 'lj': 'jl'}.get(nm, nm)
    _DOTNAME[(a_, b_)] = nm; _DOTNAME[(b_, a_)] = nm
C = RN.const; V_ = RN.var
DTRACE = C(4) - 2*V_('e')          # g^mu_mu = d = 4 - 2 eps


def tdot(a, b): return V_(_DOTNAME[(a, b)])


def vec(plus=0, minus=0, perp=None):
    co = lambda v: v if isinstance(v, Poly) else C(v)
    return (co(plus), co(minus), {k: co(v) for k, v in (perp or {}).items()})


def vadd(*vs, coeffs=None):
    coeffs = coeffs or [1]*len(vs)
    p = C(0); m = C(0); perp = {}
    for c, v in zip(coeffs, vs):
        p = p + v[0]*c; m = m + v[1]*c
        for k, xx in v[2].items(): perp[k] = perp.get(k, C(0)) + xx*c
    return (p, m, {k: v for k, v in perp.items() if v})


def vneg(v): return vadd(v, coeffs=[-1])


_MD = {}


def mdot(a, b):
    key = (id(a), id(b))
    if key in _MD and _MD[key][0] is a and _MD[key][1] is b: return _MD[key][2]
    s = (a[0]*b[1] + a[1]*b[0])*Fr(1, 2)
    for p, xx in a[2].items():
        for q, yy in b[2].items():
            s = s - xx*yy*tdot(p, q)
    _MD[key] = (a, b, s)
    return s


_cnt = itertools.count()


def new_index(tag='x'): return f'{tag}{next(_cnt)}'


class T:
    """sum of terms coeff * prod(factors); factors ('V', vector, index) or ('G', i1, i2)"""
    def __init__(self, terms=None): self.terms = terms or []
    @staticmethod
    def scalar(c): return T([(C(c) if not isinstance(c, Poly) else c, [])])
    @staticmethod
    def V(v, idx): return T([(C(1), [('V', v, idx)])])
    @staticmethod
    def G(a, b): return T([(C(1), [('G', a, b)])])
    def __add__(self, o): return T(self.terms + o.terms)
    def __sub__(self, o): return T(self.terms + [(-c, f) for c, f in o.terms])
    def __mul__(self, o):
        if not isinstance(o, T): return T([(c*o, f) for c, f in self.terms])
        return T([(c1*c2, f1 + f2) for c1, f1 in self.terms for c2, f2 in o.terms])
    __rmul__ = __mul__


def _contract_term(factors):
    fs = list(factors); coeff = C(1); changed = True
    while changed:
        changed = False
        for n_, f in enumerate(fs):
            if f[0] != 'G': continue
            a, b = f[1], f[2]
            if a == b:
                coeff = coeff*DTRACE; fs.pop(n_); changed = True; break
            for m_, g in enumerate(fs):
                if m_ == n_: continue
                if g[0] == 'G' and a in (g[1], g[2]):
                    other = g[2] if g[1] == a else g[1]
                    fs[m_] = ('G', b, other); fs.pop(n_); changed = True; break
                if g[0] == 'V' and g[2] == a:
                    fs[m_] = ('V', g[1], b); fs.pop(n_); changed = True; break
            if changed: break
            raise ValueError(f'free index {a}')
    byidx = {}
    for f in fs: byidx.setdefault(f[2], []).append(f[1])
    for idx, vs in byidx.items():
        assert len(vs) == 2, (idx, len(vs))
        coeff = coeff*mdot(vs[0], vs[1])
        if not coeff: return coeff
    return coeff


def contract(t):
    tot = C(0)
    for c, f in t.terms:
        if c: tot = tot + c*_contract_term(f)
    return tot


def V3(k1, i1, k2, i2, k3, i3):
    return (T.G(i1, i2)*T.V(vadd(k1, vneg(k2)), i3) + T.G(i2, i3)*T.V(vadd(k2, vneg(k3)), i1)
            + T.G(i3, i1)*T.V(vadd(k3, vneg(k1)), i2))


def V4_parts(i1, i2, i3, i4):
    return [T.G(i1, i3)*T.G(i2, i4) - T.G(i1, i4)*T.G(i2, i3),
            T.G(i1, i2)*T.G(i3, i4) - T.G(i1, i4)*T.G(i2, i3),
            T.G(i1, i2)*T.G(i3, i4) - T.G(i1, i3)*T.G(i2, i4)]


def trace(slots):
    idxs = []; pre = T.scalar(1)
    for s in slots:
        if s[0] == 'idx': idxs.append(s[1])
        else:
            m = new_index('t'); idxs.append(m); pre = pre*T.V(s[1], m)

    def tr(ix):
        if not ix: return T.scalar(4)
        if len(ix) % 2: return T([])
        out = T([]); a = ix[0]
        for kpos in range(1, len(ix)):
            out = out + T.G(a, ix[kpos])*tr(ix[1:kpos] + ix[kpos + 1:])*((-1)**(kpos - 1))
        return out
    return pre*tr(idxs)


NVEC = vec(2, 0); NBVEC = vec(0, 2)
