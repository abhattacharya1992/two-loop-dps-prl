#!/usr/bin/env python3
# Checks of the one-loop DPS crossed-box coefficient of gg -> AA and gg -> HH
# (paper Sec. 2). Labels in brackets are those of the author's derivation notes.
# Everything is built from the Feynman rules: stripped three-gluon vertex
# (3gvertex), stripped ggA / ggH vertices (HVvertices), polarisations (pols),
# routing (boxroute). Nothing is read from the full-amplitude files; the
# 7 pi N_c m^2, 9/7 and K_1 numbers are compared against, not fed in.
#
# Frame: p1=(1,0,0,1), p2=(1,0,0,-1), so E=1, s=4, m^2 = 4 x(1-x) at q_T=0.
#
# Not here: the completeness of K_1 (hh_spin2_numerator.py, landau_boxes.py and
# ../hh_two_loop/forward_approach.py) and the quark-loop orientation factor between
# Gaunt-Stirling Eqs. (2.6)-(2.7) and (3.17)-(3.18) (../orientation/). This script
# reproduces their Section-3 numbers, i.e. a single crossed box.
#
# python3 verify_dps_one_loop.py        (sympy, numpy, scipy; ~20 s)

import itertools
import numpy as np
import sympy as S
from scipy.integrate import quad


def check(name, lhs, rhs, tol=1e-9):
    if isinstance(lhs, (S.Basic, S.MatrixBase)) or isinstance(rhs, (S.Basic, S.MatrixBase)):
        d = S.simplify(S.expand(lhs - rhs))
        ok = all(S.simplify(c) == 0 for c in d) if isinstance(d, S.MatrixBase) else d == 0
    else:
        ok = np.allclose(lhs, rhs, rtol=tol, atol=tol)
    print(('ok   ' if ok else 'FAIL ') + name)
    if not ok:
        raise AssertionError(f'{name}: lhs={lhs}, rhs={rhs}')


x, kx, ky, qx, qy = S.symbols('x k_x k_y q_x q_y', real=True)
Nc = S.symbols('N_c', positive=True)
D = x * (1 - x)
dot = kx * (kx - qx) + ky * (ky - qy)          # k.(k-q)
wedge = ky * qx - kx * qy                       # antisymmetric part, integrates to 0
kp = kx + S.I * ky
qp = qx + S.I * qy

# ---------------------------------------------------------------- kinematics
def kinematics():
    s, mA2, mB2 = S.symbols('s m_A2 m_B2', positive=True)
    lam = (s - mA2 - mB2) ** 2 - 4 * mA2 * mB2
    X = (s + mA2 - mB2 + S.sqrt(lam)) / (2 * s)
    Y = (s + mA2 - mB2 - S.sqrt(lam)) / (2 * s)
    check('(xycond) s x y = m_A^2 for the (unequal) solution', s * X * Y, mA2)
    check('(xycond) s(1-x)(1-y) = m_B^2', s * (1 - X) * (1 - Y), mB2)
    beta = S.symbols('beta', positive=True)
    m2 = s * (1 - beta ** 2) / 4
    check('(xEqual) equal masses: x = (1+beta)/2', X.subs({mA2: m2, mB2: m2}), (1 + beta) / 2)


# ------------------------------------------------------------ splitting weights
def splitting():
    z = S.symbols('z', real=True)
    V = z * (1 - z)
    # |w|^2 for the three channels, Eqs. (ggweights), (gqqweights), (qgweights);
    # z is the first daughter's fraction (outgoing quark for q -> qg).
    gg = (1 / V, z ** 4 / V, (1 - z) ** 4 / V)
    gqq = (z ** 2, (1 - z) ** 2)
    qg = (1 / (1 - z), z ** 2 / (1 - z))
    check('(Pgg)  sum|w|^2 = P_gg/C_A', sum(gg), 2 * (z / (1 - z) + (1 - z) / z + V))
    check('(Pqg)  sum|w|^2 = P_qg/T_R', sum(gqq), z ** 2 + (1 - z) ** 2)
    check('(Pqq)  sum|w|^2 = P_qq/C_F', sum(qg), (1 + z ** 2) / (1 - z))
    # The weights come from the BDDK collinear amplitudes (II.4)-(II.6) with
    # |<12>|^2 = |[12]|^2 = k_T^2/[z(1-z)] (bracketmagnitude); |Split|^2 k_T^2
    # is listed, and the light-front conversion (convphi) divides by z(1-z).
    bddk = {'gg': (S.Integer(1), z ** 4, (1 - z) ** 4),
            'gqq': (z ** 2 * V, (1 - z) ** 2 * V),
            'qg': (z, z ** 3)}
    for lab, w in (('gg', gg), ('gqq', gqq), ('qg', qg)):
        check(f'(convphi) {lab}: |Split|^2 k_T^2 / [z(1-z)] = |w|^2',
              S.Matrix([a / V for a in bddk[lab]]), S.Matrix(w))
    # Pseudoscalar/scalar weights from the helicity products, Sec. "helicity
    # wavefunctions": ++ pairing minus (A) or plus (S) the two mixed pairings.
    w = (1 / S.sqrt(V), z ** 2 / S.sqrt(V), (1 - z) ** 2 / S.sqrt(V))
    pair = [wi * wi.subs(z, 1 - z) for wi in w]
    FA = 1 / V - 2 * V
    FS = 1 / V + 2 * V
    check('F_A = w++w++ - w+-w+- - w-+w-+', pair[0] - pair[1] - pair[2], FA)
    check('F_S = w++w++ + w+-w+- + w-+w-+', pair[0] + pair[1] + pair[2], FS)
    half = S.Rational(1, 2)
    check('F_A(1/2) = 7/2', FA.subs(z, half), S.Rational(7, 2))
    check('F_S(1/2) = 9/2', FS.subs(z, half), S.Rational(9, 2))
    check('(ratio97) F_S/F_A = 9/7 at threshold', FS.subs(z, half) / FA.subs(z, half), S.Rational(9, 7))


# ------------------------------------------------------------------- colour
def colour():
    J = S.I
    r3 = S.sqrt(3)
    gm = [[[0, 1, 0], [1, 0, 0], [0, 0, 0]],
          [[0, -J, 0], [J, 0, 0], [0, 0, 0]],
          [[1, 0, 0], [0, -1, 0], [0, 0, 0]],
          [[0, 0, 1], [0, 0, 0], [1, 0, 0]],
          [[0, 0, -J], [0, 0, 0], [J, 0, 0]],
          [[0, 0, 0], [0, 0, 1], [0, 1, 0]],
          [[0, 0, 0], [0, 0, -J], [0, J, 0]],
          [[1 / r3, 0, 0], [0, 1 / r3, 0], [0, 0, -2 / r3]]]
    T = [S.Matrix(m) / 2 for m in gm]
    check('(colordef) tr(T^a T^b) = delta/2, all 64 entries',
          S.Matrix(8, 8, lambda a, b: S.trace(T[a] * T[b])), S.eye(8) / 2)
    # f^{abc} from tr([T^a,T^b]T^c) = i f^{abc}/2
    F = np.zeros((8, 8, 8))
    for a, b, c in itertools.product(range(8), repeat=3):
        F[a, b, c] = float(S.re(S.simplify(2 * S.trace((T[a] * T[b] - T[b] * T[a]) * T[c]) / S.I)))
    check('(CA) f^{acd} f^{bcd} = 3 delta_ab', np.einsum('acd,bcd->ab', F, F), 3 * np.eye(8))
    check('(CF) sum_a T^a T^a = 4/3', sum((t * t for t in T), S.zeros(3)), S.Rational(4, 3) * S.eye(3))
    check('(TR) tr(T^1 T^1) = 1/2', S.trace(T[0] * T[0]), S.Rational(1, 2))


# ------------------------------------------- the gluonic crossed-box numerator
metric = S.diag(1, -1, -1, -1)
lower = lambda v: metric * v


def eps4(a, b, c, d):
    if len({a, b, c, d}) < 4:
        return 0
    seq = (int(a), int(b), int(c), int(d))     # Matrix(...) hands over sympy Integers
    return (-1) ** sum(seq[i] > seq[j] for i in range(4) for j in range(i + 1, 4))


def triple(P, Q, R):
    # (3gvertex), all momenta outgoing, indices down
    A, B, C = lower(P - Q), lower(Q - R), lower(R - P)
    return [[[metric[u, i] * A[j] + metric[i, j] * B[u] + metric[j, u] * C[i]
              for j in range(4)] for i in range(4)] for u in range(4)]


def VA(u, v):   # eps^{mu alpha nu beta} u_mu v_nu, indices up
    a, b = lower(u), lower(v)
    return S.Matrix(4, 4, lambda i, j: sum(eps4(r, i, n, j) * a[r] * b[n] for r in range(4) for n in range(4)))


def VH(u, v):   # (u.v) g^{ab} - v^a u^b
    return (u.T * metric * v)[0] * metric - v * u.T


p1 = S.Matrix([1, 0, 0, 1])
p2 = S.Matrix([1, 0, 0, -1])
e1p = S.Matrix([0, 1, S.I, 0]) / S.sqrt(2)
e2p = S.Matrix([0, 1, -S.I, 0]) / S.sqrt(2)
e2m = e1p

# Transverse momenta scaled by lam; the coefficient of lam^2 is the leading
# numerator.  O(lam^2) corrections to the longitudinal components cannot feed
# into it because each splitting vanishes linearly in k_T.
lam = S.symbols('lambda')
QA = S.Matrix([1, lam * qx, lam * qy, 2 * x - 1])
K = S.Matrix([x, lam * kx, lam * ky, x])
R = QA - K
A = p1 - K
B = p2 - R


def coeffs(M):      # [M_0, M_1, M_2] in powers of lam
    return [M.applyfunc(lambda f: S.expand(f).coeff(lam, n)) for n in range(3)]


def polyprod(P, Q):  # product of two truncated matrix polynomials
    return [sum([P[i] * Q[n - i] for i in range(n + 1)], S.zeros(4)) for n in range(3)]


def box_numerator(vertex, pol2):
    # (fullnumerator): e1.V3(-p1,k,a)  e2.V3(-p2,r,b)  V_Y(k,r)  V_Y(a,b),
    # contracted as tr(L^T V1 R V2^T) with the index pattern (alpha gamma)(beta delta).
    T1 = triple(-p1, K, A)
    T2 = triple(-p2, R, B)
    L = coeffs(S.Matrix(4, 4, lambda i, j: sum(e1p[u] * T1[u][i][j] for u in range(4))))
    Rm = coeffs(S.Matrix(4, 4, lambda i, j: sum(pol2[u] * T2[u][i][j] for u in range(4))))
    V1 = coeffs(vertex(K, R))
    V2 = coeffs(vertex(A, B))
    prod = polyprod(polyprod(polyprod([a.T for a in L], V1), Rm), [a.T for a in V2])
    return S.factor(S.trace(prod[2]))


def numerators():
    w0 = 1 - 2 * D ** 2
    wS = 1 + 2 * D ** 2
    w2 = x ** 4 + (1 - x) ** 4
    check('(AA00) N_++^AA  = 8[(1-2D^2) k.(k-q) + i(1+2D^2) k^q]', box_numerator(VA, e2p), 8 * (w0 * dot + S.I * wS * wedge))
    check('(HH00) N_++^HH  = -8[(1+2D^2) k.(k-q) + i(1-2D^2) k^q]', box_numerator(VH, e2p), -8 * (wS * dot + S.I * w0 * wedge))
    check('(AA20) N_+-^AA  = -8[x^4+(1-x)^4] k+(k+-q+)', box_numerator(VA, e2m), -8 * w2 * kp * (kp - qp))


# ------------------------------------------------- normalisation and targets
def normalisation():
    E = S.symbols('E', positive=True)
    s = 4 * E ** 2
    m2 = s * D
    a_s = S.symbols('a_s', positive=True)
    g2 = 16 * S.pi ** 2 * a_s
    # (overallpref): i^2 (vertices) (-i)^4 (propagators) g^2 N_c (2 pi i)^2/(2s) (2pi)^-4
    pref = (-1) * (+1) * g2 * Nc * (2 * S.pi * S.I) ** 2 / (2 * s) / (2 * S.pi) ** 4
    check('(overallpref) prefactor = 2 a_s N_c / s', pref, 2 * a_s * Nc / s)
    check('(T1project) e1+ . e2+ = -1', (e1p.T * metric * e2p)[0], -1)
    epq = -(qx + S.I * qy) / S.sqrt(2)          # e1+ . Q_A = e2- . Q_A
    check('(T2project) 2(e1+.Q)(e2-.Q)/q^2 = e^{2 i phi_q}', 2 * epq ** 2 / (qx ** 2 + qy ** 2), qp ** 2 / (qx ** 2 + qy ** 2))

    # Tree, numerically: t+u stripped contraction, then x i C^2 and / (-1).
    g = np.diag([1., -1., -1., -1.])
    p1n, p2n = np.array([1., 0, 0, 1]), np.array([1., 0, 0, -1])
    e1, e2 = np.array([0, 1, 1j, 0]) / np.sqrt(2), np.array([0, 1, -1j, 0]) / np.sqrt(2)
    levi = np.zeros((4,) * 4)
    for idx in itertools.permutations(range(4)):
        levi[idx] = (-1.) ** sum(idx[i] > idx[j] for i in range(4) for j in range(i + 1, 4))
    mdot = lambda a, b: a @ (g @ b)
    va = lambda u, v: np.einsum('rasb,r,s->ab', levi, g @ u, g @ v)
    vh = lambda u, v: mdot(u, v) * g - np.outer(v, u)

    def tree(vertex, el, er):
        q3 = np.array([1., 0, 0, -.28])      # any on-shell-pair split of p1+p2 along the beam
        tot = 0j
        for Q in (q3, p1n + p2n - q3):
            r = p1n - Q
            tot += np.einsum('i,j,iz,jz,z->', el, er, vertex(p1n, -r), vertex(p2n, r), np.diag(g)) / mdot(r, r)
        return tot
    check('(treeAA) AA t+u stripped contraction = -s/2 = -2', tree(va, e1, e2), -2)
    check('        HH t+u stripped contraction = +2', tree(vh, e1, e2), 2)
    for V, lab in ((va, 'AA'), (vh, 'HH')):
        check(f'        {lab} Ward identity, e1 -> p1', tree(V, p1n, e2), 0)
        check(f'        {lab} Ward identity, e2 -> p2', tree(V, e1, p2n), 0)
    # tree:  i C^2 (-s/2) / (-1) = i s/2
    check('(treeAA) M1^(0) = i s/2 with s=4', S.I * (-2) / (-1), 2 * S.I)

    # Chain (aaallfactors): pref x T1 projection x k.(k-q) part of (AA00) x I_0,
    # with I_0 = pi ln(R^2/pT^2), so the coefficient of ln pT^2 flips sign.
    w0 = 1 - 2 * D ** 2
    wS = 1 + 2 * D ** 2
    w2 = x ** 4 + (1 - x) ** 4
    FA = 1 / D - 2 * D
    FS = 1 / D + 2 * D
    cA = (2 * Nc / s) * (-1) * (8 * E ** 4 * w0) * S.pi * (-1)
    check('(AAfinalgeneral) coeff of ln pT^2 = pi N_c m_A^2 F_A(x)', cA, S.pi * Nc * m2 * FA)
    # L = ln(2pT/m) + i pi/2, so coefficient of L is twice that of ln pT^2
    check('(7piresult) threshold: 2 x coeff = 7 pi N_c m_A^2', (2 * cA / (S.pi * Nc * m2)).subs(x, S.Rational(1, 2)), 7)
    cH = (2 * Nc / s) * (-1) * (-8 * E ** 4 * wS) * S.pi * (-1)
    check('(HHfinal) coeff of ln pT^2 = -pi N_c m_H^2 F_S(x)', cH, -S.pi * Nc * m2 * FS)
    # (K1final): pref x (AA20) x I_2 = (2Nc/s)(-8E^4 W2)(-pi e^{2i phi}); strip e^{2i phi}
    c2 = (2 * Nc / s) * (-8 * E ** 4 * w2) * (-S.pi)
    check('(K1final) M2^(1)|DPS = pi N_c m^2 [x^4+(1-x)^4]/[x(1-x)]', c2, S.pi * Nc * m2 * w2 / D)
    check('(K1threshold) K_1^DPS(0) = pi N_c / 2', (c2 / m2).subs(x, S.Rational(1, 2)), S.pi * Nc / 2)


# ---------------------------------------------- massless quark boxes (GS Sec. 3)
def quark_boxes():
    J = S.I
    pauli = [S.Matrix([[0, 1], [1, 0]]), S.Matrix([[0, -J], [J, 0]]), S.diag(1, -1)]
    zero = S.zeros(2)
    gam = [S.diag(1, 1, -1, -1)] + [S.BlockMatrix([[zero, P], [-P, zero]]).as_explicit() for P in pauli]
    slash = lambda v: sum((metric[i, i] * v[i] * gam[i] for i in range(4)), S.zeros(4))

    def trace2(vecs):   # lam^2 coefficient of tr(prod slash(v))
        cur = [S.eye(4), S.zeros(4), S.zeros(4)]
        for v in vecs:
            cur = polyprod(cur, coeffs(slash(v)))
        return S.factor(S.trace(cur[2]))
    # Gaunt-Stirling routing: K ~ x p2, a=p2-K, b=p1+K-Q2, c=Q2-K  (HHtrace, ZZtrace)
    Kq = S.Matrix([x, lam * kx, lam * ky, -x])
    Q2 = S.Matrix([1, lam * qx, lam * qy, 1 - 2 * x])
    a, b, c = p2 - Kq, p1 + Kq - Q2, Q2 - Kq
    mix = kp * (kx - S.I * ky - qx + S.I * qy)      # k+ (k- - q-)
    check('(quarknum) N_HH   = -16 x(1-x) k+(k- - q-)', trace2([Kq, e1p, a, b, e2p, c]), -16 * D * mix)
    check('(quarknum) N_ZZ++++ = 16 x^2 k+(k- - q-)', trace2([e2p, Kq, e1p, a, e1p, b, e2p, c]), 16 * x ** 2 * mix)
    check('(quarknum) N_ZZ++-- = 16 (1-x)^2 k+(k- - q-)', trace2([e1p, Kq, e1p, a, e2p, b, e2p, c]), 16 * (1 - x) ** 2 * mix)
    # Real parts with s=4, m^2=4x(1-x): -4 m^2, 4 s x^2, 4 s (1-x)^2 times k.(k-q) (quarkshort);
    # then L = (2 pi i)^2/(2s) N I_0 with I_0 = pi ln(R^2/pT^2) (pinchmaster, I0result).
    s, m2 = S.symbols('s m2', positive=True)
    xb = (1 + S.sqrt(1 - 4 * m2 / s)) / 2
    pinch = (2 * S.pi * S.I) ** 2 / (2 * s) * S.pi * (-1)     # coefficient of ln pT^2
    check('(masslessHH) L_HH = -8 pi^3 m^2/s ln pT^2', pinch * (-4 * m2), -8 * S.pi ** 3 * m2 / s)
    check('(ZZpp) L_ZZ++++ = 4 pi^3/s [s - 2m^2 + s beta] ln pT^2',
          pinch * 4 * s * xb ** 2, 4 * S.pi ** 3 / s * (s - 2 * m2 + s * S.sqrt(1 - 4 * m2 / s)))
    check('(ZZmm) L_ZZ++-- = 4 pi^3/s [s - 2m^2 - s beta] ln pT^2',
          pinch * 4 * s * (1 - xb) ** 2, 4 * S.pi ** 3 / s * (s - 2 * m2 - s * S.sqrt(1 - 4 * m2 / s)))
    # GS Eqs. (2.6)-(2.7) print twice these ZZ coefficients; not resolved here.


# ------------------------------------------------ the transverse tensor integral
def I_tensor(qv, R):
    # I_ij = int_{|k|<R} d^2k k_i (k-q)_j / (k^2 (k-q)^2), by direct quadrature
    # in polar coordinates; radial breakpoint at |q|, angular at phi_q.
    q = np.hypot(*qv)
    phq = np.arctan2(qv[1], qv[0])

    def ang(rho, i, j):
        def f(phi):
            k = rho * np.array([np.cos(phi), np.sin(phi)])
            kq = k - qv
            return rho * k[i] * kq[j] / (rho * rho * (kq @ kq))
        return quad(f, phq - np.pi, phq + np.pi, points=[phq], limit=400)[0]
    out = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            out[i, j] = quad(lambda r: ang(r, i, j), 0, q, limit=200)[0] + quad(lambda r: ang(r, i, j), q, R, limit=200)[0]
    return out


def tensor_integral():
    # Angular densities from the contour, (I0angular) and the residue below
    # (I2contour): scalar 2 pi/rho for rho>q, zero inside; spin-2 -2 pi rho/q^2
    # inside, zero outside (the pole z = rho/q is inside the unit circle only
    # for rho<q).  Checked by an angular Fourier sum with q=(1,0).
    phi = np.linspace(0, 2 * np.pi, 16384, endpoint=False)
    for r in (.15, .6, .88, 1.12, 2.1):
        kpl = r * np.exp(1j * phi)
        den = np.abs(kpl - 1) ** 2
        check(f'(I2contour) spin-2 angular density, rho={r}', 2 * np.pi * np.mean(r * kpl * (kpl - 1) / (r * r * den)),
              -2 * np.pi * r if r < 1 else 0)
        check(f'(I0angular) scalar angular density, rho={r}', 2 * np.pi * np.mean(r * (r * r - r * np.cos(phi)) / (r * r * den)),
              0 if r < 1 else 2 * np.pi / r)
    # Full 2-d tensor integral against (universalIij) at several q (direction
    # and magnitude) and cutoffs.  For R>q the formula is exact: the scalar
    # density is 2 pi/rho beyond q and the spin-2 one vanishes there.
    for qv, R in (((1., 0.), 10.), ((.6, .8), 5.), ((-1.3, .4), 20.), ((.3, -2.), 8.)):
        lab = f'q={qv}'
        qv = np.array(qv)
        q2 = qv @ qv
        I = I_tensor(qv, R)
        want = np.eye(2) / 2 * np.pi * np.log(R ** 2 / q2) - np.pi * (np.outer(qv, qv) / q2 - np.eye(2) / 2)
        check(f'(universalIij) I_ij by quadrature, {lab}, R={R}', I, want, tol=1e-7)
        phq = np.arctan2(qv[1], qv[0])
        check(f'(I2result) I_xx - I_yy + 2i I_xy = -pi e^{{2 i phi_q}}, {lab}',
              I[0, 0] - I[1, 1] + 2j * I[0, 1], -np.pi * np.exp(2j * phq), tol=1e-7)
        check(f'(I0result) I_xx + I_yy = pi ln(R^2/q^2), {lab}', np.trace(I), np.pi * np.log(R ** 2 / q2), tol=1e-7)
    # Feynman-parameter route (fpidentity): I_2 = -q+^2 int dz z(1-z) int d^2l /(l^2 + z(1-z) q^2)^2,
    # with the l integral done numerically rather than as pi/Delta.
    for q in (0.7, 1.0, 2.5):
        inner = lambda z: quad(lambda r: 2 * np.pi * r / (r * r + z * (1 - z) * q * q) ** 2, 0, np.inf)[0]
        val = -quad(lambda z: z * (1 - z) * inner(z), 0, 1)[0] * q * q      # q along x: q+^2 = q^2
        check(f'(fpidentity) Feynman-parameter I_2 = -pi, q={q}', val, -np.pi)


if __name__ == '__main__':
    kinematics()
    splitting()
    colour()
    numerators()
    normalisation()
    quark_boxes()
    tensor_integral()
