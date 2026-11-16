"""The scalar-loop graphs of the one-loop collinear function J1 for g -> gg: a complex scalar in the adjoint
representation (the six real scalars of N = 4 super Yang-Mills are three of them; the 2 eps epsilon-scalars of
dimensional reduction are eps of them).  Same conventions as gg_diagrams.py (D = d + i g A, all momenta into the
vertex, the factors g^n and the colour structures taken out):
    scalar propagator i/q^2,
    scalar vertex      +g f^{o b i} (p_in + p_out)^mu   (o outgoing scalar, b gluon, i incoming scalar;
                       the ghost vertex of gg_diagrams.py is the same with p_out alone),
    seagull            +i g^2 g^{mu nu} [f^{x o b1} f^{x i b2} + f^{x o b2} f^{x i b1}]   (gluons b1 mu, b2 nu),
    no sign for the closed loop, no symmetry factor (complex field; a real adjoint scalar is half of this).
Both are read off L = (D_mu phi)^a* (D^mu phi)^a with (D_mu phi)^a = d_mu phi^a + g f^{cab} A^c_mu phi^b.
'colour' is computed in colour() relative to the tree f^{a0 c d}, as {'ns': .} in units of C_A, per complex adjoint
scalar.  The identity checked by check_fdh.py: gluon-loop graphs in the four-dimensional-helicity scheme (metric
traces g^mu_mu = 4, loop momenta in d dimensions) = the same graphs in conventional dimensional regularisation
+ eps times these scalar graphs (2 eps real epsilon-scalars = eps complex ones).
"""
import numpy as np
import sympy as sp
from vec_engine import vec, vadd, vneg, T, N, NB, new_index
import gg_diagrams as G
from gg_diagrams import (P1, K, L, E, I_, J_, O, KmL, QmK, QmKpL, QmKmL, QmL, KpL,
                         D_K, D_QmK, D_L, D_KmL, D_QmKpL, D_QmKmL, D_QmL, D_KpL, ix, x)

DIAGRAMS = {}


def _reg(name, **kw):
    DIAGRAMS[name] = kw


# scalar triangles (two orientations of the charge flow), same lines as the ghost triangles of gg_diagrams.py:
# vertex A (incoming gluon, index al), B (operator 1, mu), C (operator 2, nu); factor (p_in + p_out) at each vertex
def _scalar(orient):
    al, mu, nu = ix('al', 'mu', 'nu')
    if orient == 1:      # A -> B -> C -> A: q_AB = K - l, q_BC = -l, q_CA = K - l - p1
        qAB, qBC, qCA = KmL, vneg(L), vneg(QmKpL)
        num = T.V(vadd(qCA, qAB), al)*T.V(vadd(qAB, qBC), mu)*T.V(vadd(qBC, qCA), nu)
    else:                # A -> C -> B -> A: q_AC = p1 - K + l, q_CB = l, q_BA = l - K
        qAC, qCB, qBA = QmKpL, L, vneg(KmL)
        num = T.V(vadd(qBA, qAC), al)*T.V(vadd(qCB, qBA), mu)*T.V(vadd(qAC, qCB), nu)
    return num*T.V(E, al)*T.V(O(I_, K, x), mu)*T.V(O(J_, QmK, 1 - x), nu)


# i^3 scalar propagators, (-i)^2 gluon propagators, vertices +1 each, no loop sign -> (-i)(-1) = +i
_reg('Vs1', numerator=_scalar(1), dens=[D_K, D_KmL, D_QmKpL, D_L, D_QmK], factor=sp.I, eik=False,
     colour=('fffs', 'gAB a0 gCA', 'gBC c gAB', 'gCA d gBC'))
_reg('Vs2', numerator=_scalar(2), dens=[D_K, D_KmL, D_QmKpL, D_L, D_QmK], factor=sp.I, eik=False,
     colour=('fffs', 'gAC a0 gBA', 'gCB d gAC', 'gBA c gCB'))


# seagull bubbles: (a) on the K line, (b) on the p1 - K line.  Scalar flow: from the seagull along l to the cubic
# vertex, back along the other line.  (a): seagull (p1 al, -(p1-K) nu; scalars l out, K-l in), cubic vertex
# (scalar l in, scalar l-K out, gluon -K mu): (p_in + p_out) = 2l - K.  (b): seagull (p1 al, -K mu; scalars l out,
# p1-K-l in), cubic vertex (l in, l-(p1-K) out, gluon -(p1-K) nu): 2l - (p1 - K).
def _sbub(which):
    al, mu, nu = ix('al', 'mu', 'nu')
    if which == 'a':
        num = T.G(al, nu)*T.V(vadd(L, L, vneg(K)), mu)
    else:
        num = T.G(al, mu)*T.V(vadd(L, L, vneg(QmK)), nu)
    return num*T.V(E, al)*T.V(O(I_, K, x), mu)*T.V(O(J_, QmK, 1 - x), nu)


# seagull +i, cubic vertex +1, i^2 scalar propagators, (-i)^2 gluon propagators -> i (-1)(-1) = +i
_reg('Vsbuba', numerator=_sbub('a'), dens=[D_L, D_KmL, D_K, D_QmK], factor=sp.I, eik=False,
     colour=('V4s', 'a0 d e f', 'f c e'))
_reg('Vsbubb', numerator=_sbub('b'), dens=[D_L, D_QmKmL, D_K, D_QmK], factor=sp.I, eik=False,
     colour=('V4s', 'a0 c e f', 'f d e'))


def colour(spec):
    """colour factor relative to the tree f^{a0 c d}, as {'ns': .} in units of C_A (SU(3): C_A = 3), per complex
    adjoint scalar; 'fffs' as 'fff' of gg_diagrams.colour; 'V4s': seagull (gluon1 gluon2 scalar_out scalar_in)
    and cubic vertex (out gluon in)"""
    kind = spec[0]
    if kind == 'fffs':
        r = G.colour(('fff',) + tuple(spec[1:]))
        return {'ns': r['CA']}
    if kind == 'V4s':
        letters = 'abcdefghijklmnopqrstuvwxyz'
        names = {}

        def sym(n):
            if n not in names:
                names[n] = letters[len(names)]
            return names[n]
        out = ''.join(sym(n) for n in ('a0', 'c', 'd'))
        b1, b2, o, i = (sym(n) for n in spec[1].split())
        xx = sym('x4')
        s3 = ''.join(sym(n) for n in spec[2].split())
        F = G._F
        val = (np.einsum(f'{xx}{o}{b1},{xx}{i}{b2},{s3}->' + out, F, F, F)
               + np.einsum(f'{xx}{o}{b2},{xx}{i}{b1},{s3}->' + out, F, F, F))
        return {'ns': G._ratio(val/np.where(np.abs(F) > 1e-9, F, 1))/3}
    raise ValueError(spec)


def self_energy_coefficient():
    """the gluon self-energy from a complex scalar loop, relative to the structure Gamma(eps) Gamma(1-eps)^2/Gamma(2-2eps)
    of combine_gg.self_energy: -1/(3 - 2 eps) per T_R (the fermion loop of combine_gg.py is -4(1-eps)/(3-2eps) per
    T_F n_f, the gluon and ghost loops (5-3eps)/(3-2eps) C_A).  Derived here from the Feynman-parameter integral
    int_0^1 dx (1-2x)^2 [x(1-x)]^(-eps) = B(1-eps,1-eps) - 4 B(2-eps,2-eps) and the transversality of the sum with the
    g^{mu nu} term, both checked symbolically."""
    e = sp.Symbol('eps')
    B = lambda a, b: sp.gamma(a)*sp.gamma(b)/sp.gamma(a + b)
    qq = B(1 - e, 1 - e) - 4*B(2 - e, 2 - e)                 # coefficient of q^mu q^nu, times Gamma(eps) (-q^2)^-eps
    gg = -2*B(2 - e, 2 - e)/(1 - e)                          # coefficient of g^{mu nu} q^2 from 4 l^mu l^nu/d
    assert sp.simplify(sp.gammasimp(qq + gg)) == 0, 'not transverse'
    ratio = sp.gammasimp(qq/(sp.gamma(1 - e)**2/sp.gamma(2 - 2*e)))
    assert sp.simplify(ratio - 1/(3 - 2*e)) == 0, ratio
    # the fermion loop with the same normalisation: -8 int x(1-x) [x(1-x)]^-eps ... = -4(1-eps)/(3-2eps) per T_F
    ferm = sp.gammasimp(-8*B(2 - e, 2 - e)/(sp.gamma(1 - e)**2/sp.gamma(2 - 2*e)))
    assert sp.simplify(ferm + 4*(1 - e)/(3 - 2*e)) == 0, ferm
    return -1/(3 - 2*e)


if __name__ == '__main__':
    for name, dd in DIAGRAMS.items():
        print(name, colour(dd['colour']))
    print('scalar self-energy coefficient per T_R:', self_energy_coefficient())
    e = sp.Symbol('eps')
    tot = sp.simplify((5 - 3*e)/(3 - 2*e) - 2*4*(1 - e)/(3 - 2*e) + (3 + e)*self_energy_coefficient())
    print('N = 4, FDH: gluon+ghost + 4 adjoint Weyl fermions (T_F n_f -> 2 C_A) + (3 + eps) complex adjoint scalars:', tot)
