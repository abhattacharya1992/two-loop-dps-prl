"""Colour decomposition in one pass: Laurent polynomials in N and n_f as the values of the symbols
N and Nf.

The published expressions are polynomials in N, 1/N and Nf with coefficients that depend on the
kinematics only.  Evaluating them with N and Nf replaced by the elementary polynomials of this
class, through the evaluator of the P1 package, gives the coefficient of every colour monomial
N^a n_f^b at once, exactly as the expression contains it; no choice of a monomial basis is made
beforehand, and a division by anything but a monomial or a number raises an error.  Function
arguments (Log, PolyLog, Li22) are always kinematic, so a polynomial reaching a function is also
an error.
"""
import mpmath as mp


class Laurent:
    """sum c_{a,b} N^a Nf^b with mpmath coefficients; keys (a, b) are integers."""
    __slots__ = ('terms',)

    def __init__(self, terms=None):
        self.terms = {k: v for k, v in (terms or {}).items() if v != 0}

    @classmethod
    def const(cls, c):
        return cls({(0, 0): mp.mpmathify(c)})

    @staticmethod
    def lift(other):
        return other if isinstance(other, Laurent) else Laurent.const(other)

    def is_monomial(self):
        return len(self.terms) == 1

    # arithmetic --------------------------------------------------------------------------------
    def __add__(self, other):
        o = Laurent.lift(other)
        t = dict(self.terms)
        for k, v in o.terms.items():
            t[k] = t.get(k, 0) + v
        return Laurent(t)

    __radd__ = __add__

    def __neg__(self):
        return Laurent({k: -v for k, v in self.terms.items()})

    def __sub__(self, other):
        return self + (-Laurent.lift(other))

    def __rsub__(self, other):
        return Laurent.lift(other) + (-self)

    def __mul__(self, other):
        if not isinstance(other, Laurent):
            return Laurent({k: v * other for k, v in self.terms.items()})
        t = {}
        for (a, b), v in self.terms.items():
            for (c, d), w in other.terms.items():
                k = (a + c, b + d)
                t[k] = t.get(k, 0) + v * w
        return Laurent(t)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if not isinstance(other, Laurent):
            return Laurent({k: v / other for k, v in self.terms.items()})
        if not other.is_monomial():
            raise ValueError('division by a polynomial in N, Nf')
        ((c, d), w), = other.terms.items()
        return Laurent({(a - c, b - d): v / w for (a, b), v in self.terms.items()})

    def __rtruediv__(self, other):
        if not self.is_monomial():
            raise ValueError('division by a polynomial in N, Nf')
        ((c, d), w), = self.terms.items()
        return Laurent({(-c, -d): mp.mpmathify(other) / w})

    def __pow__(self, n):
        if n != int(n):
            raise ValueError('Laurent polynomial exponent must be an integer')
        n = int(n)
        if n < 0:
            return 1 / (self ** (-n))
        out = Laurent.const(1)
        base = self
        while n:
            if n & 1:
                out = out * base
            base = base * base
            n >>= 1
        return out

    def __eq__(self, other):
        return isinstance(other, Laurent) and self.terms == other.terms

    def __repr__(self):
        return ' + '.join(f'({mp.nstr(v, 8)}) N^{a} Nf^{b}' for (a, b), v in sorted(self.terms.items()))

    # conversion --------------------------------------------------------------------------------
    def monomials(self):
        return sorted(self.terms)

    def at(self, N, nf):
        """Numerical value at given N and n_f."""
        N, nf = mp.mpmathify(N), mp.mpmathify(nf)
        return sum((v * N ** a * nf ** b for (a, b), v in self.terms.items()), mp.mpc(0))


N = Laurent({(1, 0): mp.mpf(1)})
NF = Laurent({(0, 1): mp.mpf(1)})


def monomial_name(key):
    """Names used in the point files: N2, N_Nf, Nf2, Nf_over_N, and N^a_Nf^b otherwise."""
    a, b = key
    names = {(2, 0): 'N2', (1, 1): 'N_Nf', (0, 2): 'Nf2', (-1, 1): 'Nf_over_N', (1, 0): 'N', (0, 1): 'Nf',
             (0, 0): '1'}
    return names.get(key, f'N^{a}_Nf^{b}')
