"""One-loop collinear function J1 for g -> q qbar (n-collinear sector, Feynman gauge, Wilson lines W1 at x1
for the incoming antiquark r and W2 at x2 for the incoming quark rbar, eta regulator on the Wilson lines).

Momenta (p1^+ = 1): K = (x, u, k) flows into operator 1 (quark field at x1), p1 - K into operator 2.
Loop momentum l = (z, w, l_perp).  Every diagram is written as a Dirac string in matrix order from operator 1
to operator 2, with
    S(q) = i qslash/(q^2 + i0),   quark-gluon vertex -i g gamma^alpha t^c,   gluon -i g_{ab}/(l^2 + i0),
    W1 vertex (incoming antiquark): -g t^c nbar^b/(l^+ + i0),   W2 vertex (incoming quark): +g t^c nbar^b/(l^+ + i0)
for l flowing into the Wilson line, and the three-gluon vertex -g f^{abc}[...] (D = d + i g A).
'pref' is the product of all factors of i and g divided by g^3, times the colour factor relative to the tree
tr(t^a t^a'), as {'CF': c, 'CA': c}.  Tree: i g * t^a, n-coefficient (1/2)(x e k - (1-x) k e), K^- integral i/k2.
"""
import sympy as sp

x, z, u, w = sp.symbols('x z u w')
k2, l2, lk = sp.symbols('k2 l2 lk')


def mom(plus, minus, perp):
    return ('slash', plus, minus, perp)


K = mom(x, u, {'k': 1})
Kmp = mom(x - 1, u, {'k': 1})                         # K - p1
p1 = mom(1, 0, {})


def q2(c_k, c_l):
    """(c_k k + c_l l)^2 in terms of k2, l2, lk"""
    return c_k**2*k2 + 2*c_k*c_l*lk + c_l**2*l2


def den(plus, minus, ck, cl):
    """q^2 = q^+ q^- - qperp^2"""
    return sp.expand(plus*minus - q2(ck, cl))


L2 = den(z, w, 0, 1)                                  # l^2
DIAGRAMS = {}

# (8) W2-q: [K, nbar, K+l, eps, K+l-p1], dens K^2, (K+l)^2, (K+l-p1)^2, l^2; prefactor i^3 (-ig)^2 (-i) (+g) = g^3
DIAGRAMS['W2q'] = dict(
    strings=[(1, [K, ('nbar',), mom(x + z, u + w, {'k': 1, 'l': 1}), ('eps',), mom(x + z - 1, u + w, {'k': 1, 'l': 1})])],
    dens=[den(x, u, 1, 0), den(x + z, u + w, 1, 1), den(x + z - 1, u + w, 1, 1), L2],
    pref={'CF': 1, 'CA': -sp.Rational(1, 2)}, eik=True)

# (6) W2-qbar: [K, eps, K-p1, nbar, K-p1+l]; colour tr(t^a t^c t^c t^a') = C_F
DIAGRAMS['W2qb'] = dict(
    strings=[(1, [K, ('eps',), Kmp, ('nbar',), mom(x - 1 + z, u + w, {'k': 1, 'l': 1})])],
    dens=[den(x, u, 1, 0), den(x - 1, u, 1, 0), den(x - 1 + z, u + w, 1, 1), L2],
    pref={'CF': 1, 'CA': 0}, eik=True)

# (5) W1-q: [K-l, nbar, K, eps, K-p1]; prefactor -g^3, colour C_F
DIAGRAMS['W1q'] = dict(
    strings=[(1, [mom(x - z, u - w, {'k': 1, 'l': -1}), ('nbar',), K, ('eps',), Kmp])],
    dens=[den(x - z, u - w, 1, -1), den(x, u, 1, 0), den(x - 1, u, 1, 0), L2],
    pref={'CF': -1, 'CA': 0}, eik=True)

# (7) W1-qbar: [K-l, eps, K-l-p1, nbar, K-p1]; prefactor -g^3, colour C_F - C_A/2
DIAGRAMS['W1qb'] = dict(
    strings=[(1, [mom(x - z, u - w, {'k': 1, 'l': -1}), ('eps',), mom(x - z - 1, u - w, {'k': 1, 'l': -1}), ('nbar',), Kmp])],
    dens=[den(x - z, u - w, 1, -1), den(x - z - 1, u - w, 1, -1), den(x - 1, u, 1, 0), L2],
    pref={'CF': -1, 'CA': sp.Rational(1, 2)}, eik=True)

# (3) V_qqbar: [K, gam, K-l, eps, K-l-p1, gam, K-p1]; prefactor i^4 (-ig)^3 (-i) = g^3, colour C_F - C_A/2
DIAGRAMS['Vqq'] = dict(
    strings=[(1, [K, ('gam', 1), mom(x - z, u - w, {'k': 1, 'l': -1}), ('eps',), mom(x - z - 1, u - w, {'k': 1, 'l': -1}),
                  ('gam', 1), Kmp])],
    dens=[den(x, u, 1, 0), den(x - z, u - w, 1, -1), den(x - z - 1, u - w, 1, -1), den(x - 1, u, 1, 0), L2],
    pref={'CF': 1, 'CA': -sp.Rational(1, 2)}, eik=False)

# (4) V_ggq: p1 -> gluon b (l, to the quark side) + gluon c (p1 - l, to the antiquark side).
# [K, gamma^b, K-l, gamma^c, K-p1] with the vertex -g f^{abc}[eps^b (p1+l)^c + g^{bc} eps.(p1-2l) + eps^c (l-2p1)^b];
# prefactor x colour: (-i)(-g^2)(-1)(-g)(i C_A/2) = -(C_A/2) g^3.  eps.(p1 - 2l) = 2 eps.l_perp(Eucl.) = -(e l + l e)
Kml = mom(x - z, u - w, {'k': 1, 'l': -1})
DIAGRAMS['Vggq'] = dict(
    strings=[(1, [K, ('eps',), Kml, mom(1 + z, w, {'l': 1}), Kmp]),
             ('-(el+le)', [K, ('gam', 1), Kml, ('gam', 1), Kmp]),
             (1, [K, mom(z - 2, w, {'l': 1}), Kml, ('eps',), Kmp])],
    dens=[den(x, u, 1, 0), den(x - z, u - w, 1, -1), den(x - 1, u, 1, 0), L2, den(1 - z, -w, 0, -1)],
    pref={'CF': 0, 'CA': -sp.Rational(1, 2)}, eik=False)

# (9) W1-g: p1 -> gluon c (l, into W1) + gluon b (p1 - l, to the fermion line): [K-l, gamma^b, K-p1];
# vertex terms: (2 - z) eps^b and nbar^b x eps.(2l - p1) = nbar^b (-2 eps.l_E) = nbar^b (e l + l e);
# prefactor x colour: (-i g^3)(-(i/2) C_A) = -(C_A/2) g^3
DIAGRAMS['W1g'] = dict(
    strings=[(2 - z, [mom(x - z, u - w, {'k': 1, 'l': -1}), ('eps',), Kmp]),
             ('el+le', [mom(x - z, u - w, {'k': 1, 'l': -1}), ('nbar',), Kmp])],
    dens=[den(x - z, u - w, 1, -1), den(x - 1, u, 1, 0), den(1 - z, -w, 0, -1), L2],
    pref={'CF': 0, 'CA': -sp.Rational(1, 2)}, eik=True)

# (10) W2-g: [K, gamma^b, K-p1+l]; prefactor x colour: (i g^3)(i/2 C_A) = -(C_A/2) g^3
DIAGRAMS['W2g'] = dict(
    strings=[(2 - z, [K, ('eps',), mom(x - 1 + z, u + w, {'k': 1, 'l': 1})]),
             ('el+le', [K, ('nbar',), mom(x - 1 + z, u + w, {'k': 1, 'l': 1})])],
    dens=[den(x, u, 1, 0), den(x - 1 + z, u + w, 1, 1), den(1 - z, -w, 0, -1), L2],
    pref={'CF': 0, 'CA': -sp.Rational(1, 2)}, eik=True)
