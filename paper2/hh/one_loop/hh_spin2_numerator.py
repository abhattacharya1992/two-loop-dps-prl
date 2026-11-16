#!/usr/bin/env python3
# The DPS value of the opposite-helicity one-loop amplitude for gg -> HH.
#
# verify_dps_one_loop.py derives the +- crossed-box numerator only for CP-odd vertices,
# (AA20): N_+-^AA = -8 E^4 [x^4 + (1-x)^4] k+(k+ - q+). The completeness argument
# (paper Sec. 2.5) is tested numerically on gg -> HH, so the HH numerator is needed as
# well. This script repeats the (AA20) calculation of verify_dps_one_loop.py with
# the CP-even vertex V_H^{ab}(u,v) = (u.v) g^{ab} - v^a u^b and the same routing,
# polarisations and three-gluon vertex (functions copied from that script), and
# then the normalisation chain (K1final) of the same script.
#
# Result: N_+-^HH = N_+-^AA, so M_2^(1)|DPS = m^2 K_1(x) for HH as for AA, with
# K_1 = pi N_c [x^4 + (1-x)^4]/[x(1-x)]. Per unit N_c and with x = w = (1+beta)/2
# this is the pi*Ktilde that forward_approach.py compares with.
#
#   python3 hh_spin2_numerator.py          (sympy; about twenty seconds)
import sympy as S

x, kx, ky, qx, qy = S.symbols('x k_x k_y q_x q_y', real=True)
Nc = S.symbols('N_c', positive=True)
D = x*(1-x)
kp = kx + S.I*ky
qp = qx + S.I*qy
metric = S.diag(1, -1, -1, -1)
lower = lambda v: metric*v


def eps4(a, b, c, d):
    if len({a, b, c, d}) < 4:
        return 0
    seq = (int(a), int(b), int(c), int(d))
    return (-1)**sum(seq[i] > seq[j] for i in range(4) for j in range(i+1, 4))


def triple(P, Q, R):
    # three-gluon vertex, all momenta outgoing, indices down
    A, B, C = lower(P-Q), lower(Q-R), lower(R-P)
    return [[[metric[u, i]*A[j] + metric[i, j]*B[u] + metric[j, u]*C[i]
              for j in range(4)] for i in range(4)] for u in range(4)]


def VA(u, v):   # eps^{mu alpha nu beta} u_mu v_nu, indices up
    a, b = lower(u), lower(v)
    return S.Matrix(4, 4, lambda i, j: sum(eps4(r, i, n, j)*a[r]*b[n] for r in range(4) for n in range(4)))


def VH(u, v):   # (u.v) g^{ab} - v^a u^b
    return (u.T*metric*v)[0]*metric - v*u.T


# frame and routing of verify_dps_one_loop.py: E = 1, s = 4, transverse momenta scaled by lam
p1 = S.Matrix([1, 0, 0, 1])
p2 = S.Matrix([1, 0, 0, -1])
e1p = S.Matrix([0, 1, S.I, 0])/S.sqrt(2)
e2m = e1p                                  # helicity - for the gluon moving along -z
lam = S.symbols('lambda')
QA = S.Matrix([1, lam*qx, lam*qy, 2*x-1])
K = S.Matrix([x, lam*kx, lam*ky, x])
R = QA-K
A = p1-K
B = p2-R


def coeffs(M):
    return [M.applyfunc(lambda f: S.expand(f).coeff(lam, n)) for n in range(3)]


def polyprod(P, Q):
    return [sum([P[i]*Q[n-i] for i in range(n+1)], S.zeros(4)) for n in range(3)]


def box_numerator(vertex, pol2):
    # e1.V3(-p1,k,a) e2.V3(-p2,r,b) V_Y(k,r) V_Y(a,b); the lam^2 term is the leading numerator
    T1 = triple(-p1, K, A)
    T2 = triple(-p2, R, B)
    L = coeffs(S.Matrix(4, 4, lambda i, j: sum(e1p[u]*T1[u][i][j] for u in range(4))))
    Rm = coeffs(S.Matrix(4, 4, lambda i, j: sum(pol2[u]*T2[u][i][j] for u in range(4))))
    V1 = coeffs(vertex(K, R))
    V2 = coeffs(vertex(A, B))
    prod = polyprod(polyprod(polyprod([a.T for a in L], V1), Rm), [a.T for a in V2])
    return S.factor(S.trace(prod[2]))


w2 = x**4+(1-x)**4
NA = box_numerator(VA, e2m)
NH = box_numerator(VH, e2m)
target = -8*w2*kp*(kp-qp)
okA = S.simplify(S.expand(NA-target)) == 0
okH = S.simplify(S.expand(NH-target)) == 0
print('N_+-^AA =', NA)
print('N_+-^HH =', NH)
print(('ok   ' if okA else 'FAIL ') + '(AA20) N_+-^AA = -8[x^4+(1-x)^4] k+(k+-q+)   (verify_dps_one_loop.py, reproduced)')
print(('ok   ' if okH else 'FAIL ') + '       N_+-^HH = -8[x^4+(1-x)^4] k+(k+-q+)   (same numerator for CP-even vertices)')

# (K1final): prefactor 2 a_s N_c/s, numerator coefficient, I_2 = -pi e^{2 i phi_q};
# strip e^{2 i phi_q} and a_s. E = 1 here, so s = 4 and m^2 = s x(1-x).
s = 4
m2 = s*D
c2 = (2*Nc/s)*(-8*w2)*(-S.pi)
okK = S.simplify(c2-S.pi*Nc*m2*w2/D) == 0
print(('ok   ' if okK else 'FAIL ') + 'M_2^(1)|DPS = m^2 K_1(x), K_1 = pi N_c [x^4+(1-x)^4]/[x(1-x)], for HH as for AA')
if not (okA and okH and okK):
    raise SystemExit(1)
