"""Sector decomposition and numerical Laurent expansion of the parametric integrals of the collinear functions.

An integral is  C(eps, eta) * int_{[0,1]^n} dv  prod_j P_j(v)^(e_j)  with P_j polynomials (positive inside the
cube) and e_j = n_j + m_j eps + l_j eta.  The cube is split at 1/2 in every variable and each half reflected,
so that all zeros on the boundary sit on the faces through the origin; iterated sector decomposition at the
origin then factorises every polynomial into a monomial times a polynomial that does not vanish on the cube.
Each sector integral  int dv prod v_i^(a_i + b_i eps + c_i eta) H(v)  is continued to the neighbourhood of
eps = eta = 0 by subtracting H at v_i = 0 for every exponent with integer part -1, and evaluated numerically for
complex (eps, eta); the Laurent coefficients are read off by FFT on circles (eta first, |eta| << |eps|).
"""
import itertools
import numpy as np
import sympy as sp
from lc_engine import eps
from factors import eta


# ----------------------------------------------------------------------------- decomposition
def _monomial_factor(P, vs):
    """write P = prod v^k * Q with Q not divisible by any v: returns (dict v->k, Q)"""
    ks = {}
    Q = sp.expand(P)
    for v in vs:
        k = 0
        while Q != 0 and sp.expand(Q.subs(v, 0)) == 0:
            Q = sp.cancel(Q/v)
            Q = sp.expand(Q)
            k += 1
        if k:
            ks[v] = k
    return ks, Q


def _vanishes_at_origin(Q, vs):
    return sp.expand(Q.subs({v: 0 for v in vs})) == 0


def _positive_on_cube(Q, vs, n=5):
    pts = np.linspace(0, 1, n)
    f = sp.lambdify(vs, Q.subs({eps: sp.Rational(1, 100), eta: sp.Rational(1, 1000)}), 'numpy')
    for p in itertools.product(pts, repeat=len(vs)):
        val = complex(f(*p)) if len(vs) else complex(Q)
        if not (val.real > 0):
            return False
    return True


def decompose(facs, vs):
    """facs: {P: exponent}, vs: variables on [0,1]^n.  Returns a list of sectors
    (jac, mono, rest, vmap) with jac a sympy Jacobian factor (regular), mono {new var: exponent},
    rest {Q: exponent} (Q non-vanishing on the cube) and vmap the substitution old var -> expression."""
    vs = list(vs)
    out = []
    # 1. split every variable at 1/2 and reflect the upper half
    for halves in itertools.product((0, 1), repeat=len(vs)):
        sub = {}
        jac = sp.Integer(1)
        for v, h in zip(vs, halves):
            sub[v] = v/2 if h == 0 else 1 - v/2
            jac *= sp.Rational(1, 2)
        f2 = {}
        for P, e in facs.items():
            P2 = sp.expand(P.subs(sub, simultaneous=True))
            f2[P2] = sp.expand(f2.get(P2, 0) + e)
        out.extend(_iterate(f2, vs, jac, dict(sub)))
    return out


def _iterate(facs, vs, jac, vmap, depth=0):
    """iterated sector decomposition at the origin"""
    assert depth < 6, 'sector decomposition does not terminate'
    mono = {v: 0 for v in vs}
    rest = {}
    bad = False
    for P, e in facs.items():
        ks, Q = _monomial_factor(P, vs)
        for v, k in ks.items():
            mono[v] = sp.expand(mono[v] + k*e)
        if Q.is_number:
            jac *= Q**e if sp.sympify(e).is_integer else sp.Abs(Q)**e
            continue
        harmless = sp.sympify(e).is_integer and sp.sympify(e) > 0
        if not harmless and _vanishes_at_origin(Q, vs):
            bad = True
        rest[Q] = sp.expand(rest.get(Q, 0) + e)

    def _harmless(e):
        e = sp.sympify(e)
        return bool(e.is_integer and e > 0)

    if not bad:
        for Q, e in rest.items():
            if not _harmless(e):
                assert _positive_on_cube(Q, vs), ('polynomial not positive on the cube', Q, vmap)
        return [(jac, mono, rest, vmap)]
    # the minimal set of variables whose vanishing makes a bad polynomial vanish identically
    badQ = [Q for Q, e in rest.items() if not _harmless(e) and _vanishes_at_origin(Q, vs)][0]
    S = None
    for size in range(1, len(vs) + 1):
        for comb in itertools.combinations(vs, size):
            if sp.expand(badQ.subs({v: 0 for v in comb})) == 0:
                S = list(comb)
                break
        if S:
            break
    assert S and len(S) >= 2, ('singular set not resolvable', badQ, S)
    # primary sectors in the variables of S: v_k = t, v_i = t u_i (i in S, i != k)
    res = []
    for vk in S:
        sub = {v: (vk if v == vk else (v*vk if v in S else v)) for v in vs}
        f2 = {}
        # the monomials already extracted are re-expressed
        for v, ex in mono.items():
            if ex != 0:
                key = sp.expand(sub[v])
                f2[key] = sp.expand(f2.get(key, 0) + ex)
        for Q, e in rest.items():
            Q2 = sp.expand(Q.subs(sub, simultaneous=True))
            f2[Q2] = sp.expand(f2.get(Q2, 0) + e)
        # Jacobian of the substitution: vk^(|S|-1)
        f2[vk] = sp.expand(f2.get(vk, 0) + len(S) - 1)
        vmap2 = {old: sp.expand(expr.subs(sub, simultaneous=True)) for old, expr in vmap.items()}
        # monomial pieces like (u*t) must be split into factors
        f3 = {}
        jac2 = jac
        centre = {v: sp.Rational(1, 2) for v in vs}
        for P, e in f2.items():
            for fac in sp.Mul.make_args(sp.factor(P)):
                pw = 1
                if isinstance(fac, sp.Pow):
                    fac, pw = fac.args
                ee = sp.expand(pw*e)
                if fac.is_number:
                    jac2 *= fac**ee if sp.sympify(ee).is_integer else sp.Abs(fac)**ee
                    continue
                if sp.N(fac.subs(centre).subs({eps: sp.Rational(1, 100), eta: sp.Rational(1, 1000)})) < 0:
                    fac = sp.expand(-fac)
                    if sp.sympify(ee).is_integer:
                        jac2 *= (-1)**ee
                f3[fac] = sp.expand(f3.get(fac, 0) + ee)
        res.extend(_iterate(f3, vs, jac2, vmap2, depth + 1))
    return res


# ----------------------------------------------------------------------------- numerics
_GL = {}


def _gl(n):
    if n not in _GL:
        x, w = np.polynomial.legendre.leggauss(n)
        _GL[n] = ((x + 1)/2, w/2)
    return _GL[n]


class Sector:
    def __init__(self, jac, mono, rest, vs, const):
        self.vs = list(vs)
        self.mono = mono
        self.const = const
        args = list(vs) + [eps, eta]
        H = jac
        for Q, e in rest.items():
            H = H*Q**e
        import scipy.special as _ss
        self.Hf = sp.lambdify(args, H, [{'gamma': _ss.gamma, 'loggamma': _ss.loggamma}, 'numpy'])
        self.cf = sp.lambdify([eps, eta], const, 'mpmath')
        self.expo = []
        for v in vs:
            e = sp.expand(mono.get(v, 0))
            a = e.subs({eps: 0, eta: 0})
            b = e.coeff(eps)
            c = e.coeff(eta)
            assert sp.expand(e - a - b*eps - c*eta) == 0
            self.expo.append((int(a), complex(b), complex(c)))
            assert int(a) >= -1, ('power divergence', v, e)

    def value(self, ep, et, n=48, p=3):
        """the sector integral at complex (eps, eta)"""
        dim = len(self.vs)
        kap = [complex(a + b*ep + c*et) for a, b, c in self.expo]
        sub = [a == -1 for a, b, c in self.expo]
        t, wt = _gl(n)
        tot = 0j
        # loop over subsets of subtracted variables set to zero (the h(0)/kappa terms)
        for mask in itertools.product((0, 1), repeat=dim):
            # mask[i] = 1: variable i replaced by 0 (term h(0)/(kappa_i + 1)) -- only for subtracted variables
            if any(m and not s for m, s in zip(mask, sub)):
                continue
            free = [i for i in range(dim) if not mask[i]]
            coef = 1.0 + 0j
            for i in range(dim):
                if mask[i]:
                    coef /= (kap[i] + 1)
            # integrand over the free variables: prod v^kap * [H with subtractions over the free subtracted vars]
            if not free:
                vals = [0.0]*dim
                tot += coef*complex(self.Hf(*vals, ep, et))
                continue
            grids = np.meshgrid(*([t]*len(free)), indexing='ij')
            # substitution v = s^p on each free variable
            V = [g**p for g in grids]
            J = np.ones_like(grids[0])
            for g in grids:
                J = J*p*g**(p - 1)
            W = np.ones_like(grids[0])
            for k in range(len(free)):
                sh = [1]*len(free); sh[k] = n
                W = W*wt.reshape(sh)
            # subtracted H: sum over subsets of the free subtracted variables
            fsub = [i for i in free if sub[i]]
            Hs = np.zeros_like(grids[0], dtype=complex)
            for smask in itertools.product((0, 1), repeat=len(fsub)):
                args = []
                vi = {}
                for k, i in enumerate(free):
                    vi[i] = V[k]
                for i in range(dim):
                    if mask[i]:
                        vi[i] = np.zeros_like(grids[0])
                sign = 1
                for bit, i in zip(smask, fsub):
                    if bit:
                        vi[i] = np.zeros_like(grids[0]); sign = -sign
                args = [vi[i] for i in range(dim)]
                Hs = Hs + sign*np.asarray(self.Hf(*args, ep, et), dtype=complex)*np.ones_like(grids[0])
            mon = np.ones_like(grids[0], dtype=complex)
            for k, i in enumerate(free):
                mon = mon*V[k]**kap[i]
            tot += coef*np.sum(W*J*mon*Hs)
        return complex(self.cf(ep, et))*tot


def laurent(fun, r_eps=0.2, n_eps=16, r_eta_frac=0.05, n_eta=8, eta_orders=(-1, 0), eps_orders=(-2, -1, 0),
            has_eta=True):
    """Laurent coefficients c[(i, j)] of fun(eps, eta) ~ sum eta^i eps^j, eta expanded first.
    fun must be analytic in eta on |eta| <= r_eta except at eta = 0, and in eps on |eps| = r_eps except at 0."""
    th_e = 2*np.pi*(np.arange(n_eps) + 0.5)/n_eps
    th_h = 2*np.pi*(np.arange(n_eta) + 0.5)/n_eta
    res = {}
    # inner: eta-Laurent coefficients at each eps on the circle
    inner = {i: np.zeros(n_eps, dtype=complex) for i in eta_orders}
    for a, te in enumerate(th_e):
        ep = r_eps*np.exp(1j*te)
        if not has_eta:
            v0 = fun(ep, 0.0)
            for i in eta_orders:
                inner[i][a] = v0 if i == 0 else 0.0
            continue
        rh = r_eta_frac*abs(ep)
        vals = np.array([fun(ep, rh*np.exp(1j*th)) for th in th_h])
        for i in eta_orders:
            # c_i = (1/N) sum f(eta_k) eta_k^(-i)
            inner[i][a] = np.mean(vals*(rh*np.exp(1j*th_h))**(-i))
    for i in eta_orders:
        epsv = r_eps*np.exp(1j*th_e)
        for j in eps_orders:
            res[(i, j)] = np.mean(inner[i]*epsv**(-j))
    return res
