"""Lorentz-tensor algebra for the gluon collinear functions (g -> gg) in light-cone components.

Vectors.  A 4-vector is (plus, minus, perp) with v^+ = nbar.v, v^- = n.v and perp a dict {name: coefficient}
over Euclidean transverse basis vectors 'k', 'l' (in D = d - 2 dimensions) and 'e', 'i', 'j' (the
polarisation of the incoming gluon and the transverse directions of the two operators).  Minkowski product
    a.b = (a^+ b^- + a^- b^+)/2 - a_perp . b_perp,
with the Euclidean products as symbols kk, ll, lk, ek, el, ik, il, jk, jl, ei, ej, ij, ee, ii, jj.
n = (2, 0, 0) and nbar = (0, 2, 0).

Tensors.  A tensor is a list of terms (coefficient, factors); a factor is ('V', vector, index) or
('G', index1, index2) (metric).  contract() sums over repeated indices (each index occurs exactly twice) and
returns a sympy expression; g^mu_mu = d.

Dirac traces of gamma matrices with vector or index slots are expanded into metrics (trace of the identity 4).
"""
import itertools
import sympy as sp

eps = sp.Symbol('eps')
d = 4 - 2*eps
METRIC_TRACE = d          # g^mu_mu; check_fdh.py sets it to 4 (four-dimensional-helicity scheme)
TNAMES = ('k', 'l', 'e', 'i', 'j')
_DOT = {}
for a, b in itertools.combinations_with_replacement(TNAMES, 2):
    s = sp.Symbol(''.join(sorted((a, b), key=TNAMES.index)) if a != b else a + a)
    _DOT[(a, b)] = s
    _DOT[(b, a)] = s
kk, ll, lk = _DOT[('k', 'k')], _DOT[('l', 'l')], _DOT[('k', 'l')]
ek, el, ik, il, jk, jl = (_DOT[('e', 'k')], _DOT[('e', 'l')], _DOT[('i', 'k')], _DOT[('i', 'l')],
                          _DOT[('j', 'k')], _DOT[('j', 'l')])
ei, ej, ij = _DOT[('e', 'i')], _DOT[('e', 'j')], _DOT[('i', 'j')]


def tdot(a, b):
    return _DOT[(a, b)]


def vec(plus=0, minus=0, perp=None):
    return (sp.sympify(plus), sp.sympify(minus), {k: sp.sympify(v) for k, v in (perp or {}).items() if v != 0})


N = vec(2, 0)
NB = vec(0, 2)


def vadd(*vs, coeffs=None):
    coeffs = coeffs or [1]*len(vs)
    p = sum(c*v[0] for c, v in zip(coeffs, vs))
    m = sum(c*v[1] for c, v in zip(coeffs, vs))
    perp = {}
    for c, v in zip(coeffs, vs):
        for k, x in v[2].items():
            perp[k] = perp.get(k, 0) + c*x
    return vec(p, m, perp)


def vneg(v):
    return vadd(v, coeffs=[-1])


def mdot(a, b):
    s = (a[0]*b[1] + a[1]*b[0])/2
    for p, x in a[2].items():
        for q, y in b[2].items():
            s -= x*y*tdot(p, q)
    return sp.expand(s)


# ----------------------------------------------------------------------------- tensors
_counter = itertools.count()


def new_index(tag='x'):
    return f'{tag}{next(_counter)}'


class T:
    """a sum of terms coeff * prod(factors)"""
    def __init__(self, terms=None):
        self.terms = terms or []

    @staticmethod
    def scalar(c):
        return T([(sp.sympify(c), [])])

    @staticmethod
    def V(v, idx):
        return T([(sp.Integer(1), [('V', v, idx)])])

    @staticmethod
    def G(a, b):
        return T([(sp.Integer(1), [('G', a, b)])])

    def __add__(self, other):
        return T(self.terms + other.terms)

    def __sub__(self, other):
        return T(self.terms + [(-c, f) for c, f in other.terms])

    def __mul__(self, other):
        if not isinstance(other, T):
            return T([(c*other, f) for c, f in self.terms])
        return T([(c1*c2, f1 + f2) for c1, f1 in self.terms for c2, f2 in other.terms])

    __rmul__ = __mul__


def _contract_term(factors):
    fs = list(factors)
    coeff = sp.Integer(1)
    # remove metrics
    changed = True
    while changed:
        changed = False
        for n_, f in enumerate(fs):
            if f[0] != 'G':
                continue
            a, b = f[1], f[2]
            if a == b:
                coeff *= METRIC_TRACE
                fs.pop(n_)
                changed = True
                break
            # find the other occurrence of a
            for m_, g in enumerate(fs):
                if m_ == n_:
                    continue
                if g[0] == 'G' and a in (g[1], g[2]):
                    other = g[2] if g[1] == a else g[1]
                    fs[m_] = ('G', b, other)
                    fs.pop(n_)
                    changed = True
                    break
                if g[0] == 'V' and g[2] == a:
                    fs[m_] = ('V', g[1], b)
                    fs.pop(n_)
                    changed = True
                    break
            if changed:
                break
            raise ValueError(f'free index {a}')
    # pair vectors
    byidx = {}
    for f in fs:
        assert f[0] == 'V'
        byidx.setdefault(f[2], []).append(f[1])
    for idx, vs in byidx.items():
        assert len(vs) == 2, (idx, len(vs))
        coeff *= mdot(vs[0], vs[1])
    return coeff


def contract(t):
    tot = sp.Integer(0)
    for c, f in t.terms:
        if c == 0:
            continue
        tot += c*_contract_term(f)
    return sp.expand(tot)


# ----------------------------------------------------------------------------- vertices (D = d + i g A)
def V3(k1, i1, k2, i2, k3, i3):
    """three-gluon vertex without the factor -g f^{a1 a2 a3}: all momenta incoming
    g^{i1 i2}(k1 - k2)^{i3} + g^{i2 i3}(k2 - k3)^{i1} + g^{i3 i1}(k3 - k1)^{i2}"""
    return (T.G(i1, i2)*T.V(vadd(k1, vneg(k2)), i3) + T.G(i2, i3)*T.V(vadd(k2, vneg(k3)), i1)
            + T.G(i3, i1)*T.V(vadd(k3, vneg(k1)), i2))


def V4_parts(i1, i2, i3, i4):
    """four-gluon vertex without -i g^2: the three Lorentz structures, multiplying in turn
    f^{a1 a2 e} f^{a3 a4 e}, f^{a1 a3 e} f^{a2 a4 e}, f^{a1 a4 e} f^{a2 a3 e}"""
    return [T.G(i1, i3)*T.G(i2, i4) - T.G(i1, i4)*T.G(i2, i3),
            T.G(i1, i2)*T.G(i3, i4) - T.G(i1, i4)*T.G(i2, i3),
            T.G(i1, i2)*T.G(i3, i4) - T.G(i1, i3)*T.G(i2, i4)]


def trace(slots):
    """Dirac trace of a product of slots, each ('idx', mu) for gamma^mu or ('vec', v) for vslash, in d dims"""
    idxs = []
    pre = T.scalar(1)
    for s in slots:
        if s[0] == 'idx':
            idxs.append(s[1])
        else:
            m = new_index('t')
            idxs.append(m)
            pre = pre*T.V(s[1], m)

    def tr(ix):
        if not ix:
            return T.scalar(4)
        if len(ix) % 2:
            return T([])
        out = T([])
        a = ix[0]
        for kpos in range(1, len(ix)):
            sign = (-1)**(kpos - 1)
            rest = ix[1:kpos] + ix[kpos + 1:]
            out = out + T.G(a, ix[kpos])*tr(rest)*sign
        return out
    return pre*tr(idxs)
