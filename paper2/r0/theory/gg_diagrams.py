"""The tree and the one-loop graphs of the collinear function J1 for g -> gg (n-collinear sector, Feynman gauge).

Matrix element <0| B^(i)c(0) B^(j)d(b) |g^a0(p1, eps)>, p1^+ = 1: K = (x, u, k) flows into operator 1, p1 - K into
operator 2; loop momentum l = (z, w, l_perp).  Rules (D = d + i g A; all momenta into the vertex; the factors -g,
-i g^2 and the colour structures are taken out):
    gluon propagator -i g/q^2, ghost propagator i/q^2, quark propagator i qslash/q^2,
    three-gluon vertex -g f^{a1 a2 a3} [g12 (k1-k2)_3 + g23 (k2-k3)_1 + g31 (k3-k1)_2],
    four-gluon vertex -i g^2 sum_r f_r L_r (vec_engine.V4_parts),
    ghost vertex +g f^{a b c} p^mu (a outgoing ghost, b gluon, c incoming ghost, p outgoing ghost momentum),
    quark vertex -i g gamma^mu t^a, closed loops -1;
    operator B^(i) = -i . B_perp (i a transverse unit vector), from B_perp = -(i/K^+) [W^dagger nbar_nu G^{nu mu} W]:
      one gluon (momentum K):            + O^sigma(K, K^+) A_sigma,
      two gluons (q1 into the Wilson line, q2 into the field strength, q1 + q2 = K):
                                          - i g f^{a b c} nbar^rho/(q1^+ + i0) O^sigma(q2, K^+) A^b_rho A^c_sigma,
      O^sigma(q, K^+) = i^sigma - (i.q) nbar^sigma/K^+   (i.q Minkowski).
The eikonal 1/(q1^+) = 1/z of the Wilson-line graphs is left out here ('eik': True) and supplied, with the
eta regulator nu^eta |z|^-eta, by the assembly, as for g -> q qbar.
'factor' collects the factors of i and the signs (no g); 'colour' is computed in colour(), relative to the tree
colour f^{a0 c d}, as {'CA': ., 'nf': .} (n_f in units of T_F n_f).
"""
import itertools
import numpy as np
import sympy as sp
from vec_engine import vec, vadd, vneg, mdot, T, V3, V4_parts, trace, N, NB, kk, ll, lk, new_index
from gg_driver import x, z, u, w

P1 = vec(1, 0)
K = vec(x, u, {'k': 1})
L = vec(z, w, {'l': 1})
E = vec(0, 0, {'e': 1})
I_ = vec(0, 0, {'i': 1})
J_ = vec(0, 0, {'j': 1})


def mom(plus, minus, ck, cl):
    return vec(plus, minus, {'k': ck, 'l': cl})


def den(plus, minus, ck, cl):
    return sp.expand(plus*minus - (ck**2*kk + 2*ck*cl*lk + cl**2*ll))


def O(dirv, q, Kplus):
    """operator vector O^sigma(q, K^+) = a^sigma - (a.q) nbar^sigma/K^+ for the transverse direction a"""
    return vadd(dirv, NB, coeffs=[1, -mdot(dirv, q)/Kplus])


# momenta used below
KmL = mom(x - z, u - w, 1, -1)              # K - l
QmK = mom(1 - x, -u, -1, 0)                 # p1 - K
QmKpL = mom(1 - x + z, w - u, -1, 1)        # p1 - K + l
QmKmL = mom(1 - x - z, -u - w, -1, -1)      # p1 - K - l
QmL = mom(1 - z, -w, 0, -1)                 # p1 - l
KpL = mom(x + z, u + w, 1, 1)               # K + l

D_K = den(x, u, 1, 0)
D_QmK = den(1 - x, -u, -1, 0)
D_L = den(z, w, 0, 1)
D_KmL = den(x - z, u - w, 1, -1)
D_QmKpL = den(1 - x + z, w - u, -1, 1)
D_QmKmL = den(1 - x - z, -u - w, -1, -1)
D_QmL = den(1 - z, -w, 0, -1)
D_KpL = den(x + z, u + w, 1, 1)


def ix(*names):
    return [new_index(n) for n in names]


# ----------------------------------------------------------------------------- tree
def tree():
    al, mu, nu = ix('al', 'mu', 'nu')
    num = T.V(E, al)*V3(P1, al, vneg(K), mu, vneg(QmK), nu)*T.V(O(I_, K, x), mu)*T.V(O(J_, QmK, 1 - x), nu)
    # (-g) (-i)^2 -> +1
    return dict(numerator=num, dens=[D_K, D_QmK], factor=1, colour=('tree',))


# ----------------------------------------------------------------------------- one loop
DIAGRAMS = {}


def _reg(name, **kw):
    DIAGRAMS[name] = kw


# gluon triangle: A (p1 -> K-l, p1-K+l), B (K-l, l -> K), C (p1-K+l -> p1-K, l)
def _tri():
    al, be, ga, de, mu, nu = ix('al', 'be', 'ga', 'de', 'mu', 'nu')
    num = (T.V(E, al)*V3(P1, al, vneg(KmL), be, vneg(QmKpL), ga)*V3(KmL, be, L, de, vneg(K), mu)
           * V3(QmKpL, ga, vneg(L), de, vneg(QmK), nu)*T.V(O(I_, K, x), mu)*T.V(O(J_, QmK, 1 - x), nu))
    # (-g)^3 (-i)^5 = i
    return num
_reg('Vtri', numerator=_tri(), dens=[D_K, D_KmL, D_QmKpL, D_L, D_QmK], factor=sp.I, eik=False,
     colour=('fff', 'a0 e f', 'e h c', 'f h d'))


# ghost triangles (two orientations), same lines as the gluon triangle
def _ghost(orient):
    al, mu, nu = ix('al', 'mu', 'nu')
    if orient == 1:      # A -> B -> C -> A: q_AB = K - l, q_BC = -l, q_CA = K - l - p1
        qAB, qBC, qCA = KmL, vneg(L), vneg(QmKpL)
        num = T.V(qAB, al)*T.V(qBC, mu)*T.V(qCA, nu)
    else:                # A -> C -> B -> A: q_AC = p1 - K + l, q_CB = l, q_BA = l - K
        qAC, qCB, qBA = QmKpL, L, vneg(KmL)
        num = T.V(qAC, al)*T.V(qBA, mu)*T.V(qCB, nu)
    num = num*T.V(E, al)*T.V(O(I_, K, x), mu)*T.V(O(J_, QmK, 1 - x), nu)
    return num
# (-1) loop, i^3 ghost propagators, (-i)^2 gluon propagators -> -i
_reg('Vgh1', numerator=_ghost(1), dens=[D_K, D_KmL, D_QmKpL, D_L, D_QmK], factor=-sp.I, eik=False,
     colour=('fff', 'gAB a0 gCA', 'gBC c gAB', 'gCA d gBC'))
_reg('Vgh2', numerator=_ghost(2), dens=[D_K, D_KmL, D_QmKpL, D_L, D_QmK], factor=-sp.I, eik=False,
     colour=('fff', 'gAC a0 gBA', 'gCB d gAC', 'gBA c gCB'))


# quark triangles (two orientations of the fermion flow)
def _quark(orient):
    al, mu, nu = ix('al', 'mu', 'nu')
    if orient == 1:      # flow A -> B -> C -> A: Tr[gA S(qCA) gC S(qBC) gB S(qAB)]
        qAB, qBC, qCA = KmL, vneg(L), vneg(QmKpL)
        tr = trace([('idx', al), ('vec', qCA), ('idx', nu), ('vec', qBC), ('idx', mu), ('vec', qAB)])
    else:                # flow A -> C -> B -> A: Tr[gA S(qBA) gB S(qCB) gC S(qAC)]
        qAC, qCB, qBA = QmKpL, L, vneg(KmL)
        tr = trace([('idx', al), ('vec', qBA), ('idx', mu), ('vec', qCB), ('idx', nu), ('vec', qAC)])
    return tr*T.V(E, al)*T.V(O(I_, K, x), mu)*T.V(O(J_, QmK, 1 - x), nu)
# (-1) (-i)^3 i^3 (-i)^2 = +1;   colour Tr[t^a0 t^d t^c] (orientation 1), Tr[t^a0 t^c t^d] (orientation 2)
_reg('Vq1', numerator=_quark(1), dens=[D_K, D_KmL, D_QmKpL, D_L, D_QmK], factor=1, eik=False,
     colour=('ttt', 'a0 d c'))
_reg('Vq2', numerator=_quark(2), dens=[D_K, D_KmL, D_QmKpL, D_L, D_QmK], factor=1, eik=False,
     colour=('ttt', 'a0 c d'))


# four-gluon bubbles: (a) on the K line, (b) on the p1 - K line; symmetry factor 1/2
def _bub(which, r):
    al, mu, nu, be, ga, be2, ga2 = ix('al', 'mu', 'nu', 'be', 'ga', 'be2', 'ga2')
    if which == 'a':     # 4-vertex (p1, -(p1-K), -l, -(K-l)); 3-vertex (l, K-l, -K)
        lor = V4_parts(al, nu, be, ga)[r]
        num = T.V(E, al)*lor*V3(L, be, KmL, ga, vneg(K), mu)*T.V(O(I_, K, x), mu)*T.V(O(J_, QmK, 1 - x), nu)
    else:                # 4-vertex (p1, -K, -l, -(p1-K-l)); 3-vertex (l, p1-K-l, -(p1-K))
        lor = V4_parts(al, mu, be, ga)[r]
        num = T.V(E, al)*lor*V3(L, be, QmKmL, ga, vneg(QmK), nu)*T.V(O(I_, K, x), mu)*T.V(O(J_, QmK, 1 - x), nu)
    return num
for r in range(3):
    # (-i)(-1)(-i)^4 (1/2) = i/2
    _reg(f'Vbuba{r+1}', numerator=_bub('a', r), dens=[D_L, D_KmL, D_K, D_QmK], factor=sp.I/2, eik=False,
         colour=('V4', r, 'a0 d e f', 'e f c'))
    _reg(f'Vbubb{r+1}', numerator=_bub('b', r), dens=[D_L, D_QmKmL, D_K, D_QmK], factor=sp.I/2, eik=False,
         colour=('V4', r, 'a0 c e f', 'e f d'))


# Wilson lines of operator 1 (q1 = l into W, q2 = K - l into the field strength): -i f^{c b e} nbar_rho O1_sigma(K-l)
def _w1(which, r=None):
    al, rho, sig, nu, be = ix('al', 'rho', 'sig', 'nu', 'be')
    op1 = T.V(NB, rho)*T.V(O(I_, KmL, x), sig)
    op2 = T.V(O(J_, QmK, 1 - x), nu)
    if which == 'g':      # A: (p1, -l, -(p1-l)); B: (p1-l, -(K-l), -(p1-K))
        num = T.V(E, al)*V3(P1, al, vneg(L), rho, vneg(QmL), be)*V3(QmL, be, vneg(KmL), sig, vneg(QmK), nu)
    elif which == 'kb':   # A: (p1, -(K-l), -(p1-K+l)); C: (p1-K+l, -l, -(p1-K))
        num = T.V(E, al)*V3(P1, al, vneg(KmL), sig, vneg(QmKpL), be)*V3(QmKpL, be, vneg(L), rho, vneg(QmK), nu)
    elif which == 'k':    # A: (p1, -K, -(p1-K)); B: (K, -l, -(K-l))
        num = T.V(E, al)*V3(P1, al, vneg(K), be, vneg(QmK), nu)*V3(K, be, vneg(L), rho, vneg(KmL), sig)
    else:                 # 4-vertex (p1, -l, -(K-l), -(p1-K))
        num = T.V(E, al)*V4_parts(al, rho, sig, nu)[r]
    return num*op1*op2


# factors: g, kb, k: (-1)^2 (-i)^4 (-i) = -i;  4g: (-i)(-i)(-i)^3 = -i
_reg('W1g', numerator=_w1('g'), dens=[D_L, D_QmL, D_KmL, D_QmK], factor=-sp.I, eik=True,
     colour=('fff', 'a0 b h', 'h e d', 'c b e'))
_reg('W1kb', numerator=_w1('kb'), dens=[D_KmL, D_QmKpL, D_L, D_QmK], factor=-sp.I, eik=True,
     colour=('fff', 'a0 e h', 'h b d', 'c b e'))
_reg('W1k', numerator=_w1('k'), dens=[D_K, D_L, D_KmL, D_QmK], factor=-sp.I, eik=True,
     colour=('fff', 'a0 h d', 'h b e', 'c b e'))
for r in range(3):
    _reg(f'W1v4{r+1}', numerator=_w1('4', r), dens=[D_L, D_KmL, D_QmK], factor=-sp.I, eik=True,
         colour=('V4W', r, 'a0 b e d', 'c b e'))


# Wilson lines of operator 2 (q1 = l into W, q2 = p1 - K - l): -i f^{d b e} nbar_rho O2_sigma(p1-K-l)
def _w2(which, r=None):
    al, rho, sig, mu, be = ix('al', 'rho', 'sig', 'mu', 'be')
    op2 = T.V(NB, rho)*T.V(O(J_, QmKmL, 1 - x), sig)
    op1 = T.V(O(I_, K, x), mu)
    if which == 'g':      # A: (p1, -l, -(p1-l)); (p1-l, -K, -(p1-K-l))
        num = T.V(E, al)*V3(P1, al, vneg(L), rho, vneg(QmL), be)*V3(QmL, be, vneg(K), mu, vneg(QmKmL), sig)
    elif which == 'k':    # A: (p1, -(K+l), -(p1-K-l)); B: (K+l, -K, -l)
        num = T.V(E, al)*V3(P1, al, vneg(KpL), be, vneg(QmKmL), sig)*V3(KpL, be, vneg(K), mu, vneg(L), rho)
    elif which == 'kb':   # A: (p1, -K, -(p1-K)); C: (p1-K, -l, -(p1-K-l))
        num = T.V(E, al)*V3(P1, al, vneg(K), mu, vneg(QmK), be)*V3(QmK, be, vneg(L), rho, vneg(QmKmL), sig)
    else:                 # 4-vertex (p1, -K, -l, -(p1-K-l))
        num = T.V(E, al)*V4_parts(al, mu, rho, sig)[r]
    return num*op1*op2


_reg('W2g', numerator=_w2('g'), dens=[D_L, D_QmL, D_K, D_QmKmL], factor=-sp.I, eik=True,
     colour=('fff', 'a0 b h', 'h c e', 'd b e'))
_reg('W2k', numerator=_w2('k'), dens=[D_KpL, D_K, D_L, D_QmKmL], factor=-sp.I, eik=True,
     colour=('fff', 'a0 h e', 'h c b', 'd b e'))
_reg('W2kb', numerator=_w2('kb'), dens=[D_K, D_QmK, D_L, D_QmKmL], factor=-sp.I, eik=True,
     colour=('fff', 'a0 c h', 'h b e', 'd b e'))
for r in range(3):
    _reg(f'W2v4{r+1}', numerator=_w2('4', r), dens=[D_K, D_L, D_QmKmL], factor=-sp.I, eik=True,
         colour=('V4W', r, 'a0 c b e', 'd b e'))


# ----------------------------------------------------------------------------- colour (SU(3), numerically)
def _su3():
    lam = np.zeros((8, 3, 3), dtype=complex)
    lam[0][0, 1] = lam[0][1, 0] = 1
    lam[1][0, 1] = -1j; lam[1][1, 0] = 1j
    lam[2][0, 0] = 1; lam[2][1, 1] = -1
    lam[3][0, 2] = lam[3][2, 0] = 1
    lam[4][0, 2] = -1j; lam[4][2, 0] = 1j
    lam[5][1, 2] = lam[5][2, 1] = 1
    lam[6][1, 2] = -1j; lam[6][2, 1] = 1j
    lam[7][0, 0] = lam[7][1, 1] = 1/np.sqrt(3); lam[7][2, 2] = -2/np.sqrt(3)
    t = lam/2
    f = np.zeros((8, 8, 8))
    for a, b, c in itertools.product(range(8), repeat=3):
        f[a, b, c] = np.real(-2j*np.trace((t[a] @ t[b] - t[b] @ t[a]) @ t[c]))
    return t, f


_T, _F = _su3()


def colour(spec):
    """colour factor relative to the tree f^{a0 c d}, as {'CA': ., 'nf': .} (SU(3): C_A = 3, T_F = 1/2)"""
    kind = spec[0]
    letters = 'abcdefghijklmnopqrstuvwxyz'
    names = {}

    def sym(n):
        if n not in names:
            names[n] = letters[len(names)]
        return names[n]
    ext = [sym('a0'), sym('c'), sym('d')]
    out = ''.join(ext)
    if kind == 'tree':
        return {'CA': 0, 'nf': 0}
    if kind == 'fff':
        subs = [''.join(sym(n) for n in s.split()) for s in spec[1:]]
        val = np.einsum(','.join(subs) + '->' + out, _F, _F, _F)
        ratio = val/np.where(np.abs(_F) > 1e-9, _F, 1)
        r = _ratio(ratio)
        return {'CA': r/3, 'nf': 0}
    if kind == 'V4':
        r_, legs, v3 = spec[1], spec[2].split(), spec[3].split()
        a1, a2, a3, a4 = (sym(n) for n in legs)
        xx = sym('x4')
        pairs = [((a1, a2), (a3, a4)), ((a1, a3), (a2, a4)), ((a1, a4), (a2, a3))][r_]
        s1 = pairs[0][0] + pairs[0][1] + xx
        s2 = pairs[1][0] + pairs[1][1] + xx
        s3 = ''.join(sym(n) for n in v3)
        val = np.einsum(f'{s1},{s2},{s3}->' + out, _F, _F, _F)
        return {'CA': _ratio(val/np.where(np.abs(_F) > 1e-9, _F, 1))/3, 'nf': 0}
    if kind == 'V4W':
        r_, legs, opc = spec[1], spec[2].split(), spec[3].split()
        a1, a2, a3, a4 = (sym(n) for n in legs)
        xx = sym('x4')
        pairs = [((a1, a2), (a3, a4)), ((a1, a3), (a2, a4)), ((a1, a4), (a2, a3))][r_]
        s1 = pairs[0][0] + pairs[0][1] + xx
        s2 = pairs[1][0] + pairs[1][1] + xx
        s3 = ''.join(sym(n) for n in opc)
        val = np.einsum(f'{s1},{s2},{s3}->' + out, _F, _F, _F)
        return {'CA': _ratio(val/np.where(np.abs(_F) > 1e-9, _F, 1))/3, 'nf': 0}
    if kind == 'ttt':
        order = spec[1].split()            # the generators in the order of the trace, names among a0, c, d
        tr = np.einsum('aij,bjk,cki->abc', _T, _T, _T)
        V = np.transpose(tr, [order.index(n) for n in ('a0', 'c', 'd')])     # V[a0, c, d]
        fpart = (V - np.transpose(V, (0, 2, 1)))/2                         # the d^{abc} part cancels (Furry)
        return {'CA': 0, 'nf': _ratio(fpart/np.where(np.abs(_F) > 1e-9, _F, 1))/0.5}   # per T_F n_f
    raise ValueError(spec)


def _ratio(arr):
    arr = np.asarray(arr)
    if arr.shape != _F.shape or not np.all(np.isfinite(arr)):
        raise ValueError('nonfinite or malformed colour tensor')
    mask = np.abs(_F) > 1e-9
    vals = arr[mask]
    if not np.allclose(arr[~mask], 0, atol=1e-9, rtol=0):
        raise ValueError('nonzero colour component where the tree vanishes')
    if not np.allclose(vals, vals[0], atol=1e-9, rtol=0):
        raise ValueError('colour structure not proportional to the tree')
    v = complex(vals[0])
    return v.real if abs(v.imag) < 1e-12 else v
