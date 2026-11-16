#!/usr/bin/env python3
# The two fermion-flow orientations of the massless-quark crossed box in
# gg -> ZZ (and gg -> HH with scalar vertices).
#
# Vertices around the loop: g(p1), Z(p3), g(p2), Z(p4); the gluons sit on
# opposite corners, which is the crossed box of Gaunt-Stirling Fig. 1.
# Orientation A: fermion flows V1 -> V3 -> V2 -> V4 -> V1.
# Orientation B: fermion flows V1 -> V4 -> V2 -> V3 -> V1.
# With loop momentum l on the line V4 -> V1 in A, the same four propagators
# appear in B with reversed momenta if l' = -(l + p1). The denominators are
# therefore identical and only the Dirac traces can differ.
#
# Charge conjugation predicts trace_B = trace_A for an even number of
# vector-type vertices: gamma^mu picks up a sign under C, gamma^mu gamma5 and 1
# do not, and the four propagators give (-1)^4. So:
#   v^2 and a^2 parts of the Z couplings: B = A   (the two orientations add)
#   v a cross terms:                      B = -A  (they cancel)
#   scalar (Yukawa) H vertices:           B = A
# The script checks this numerically at random loop momenta and at loop
# momenta collinear to p1 (the DPS region), with physical helicity states for
# the gluons and all three polarisations of each Z.
#
#   python3 orientation_check.py
import numpy as np

rng = np.random.default_rng(20261012)

# Dirac representation, metric (+,-,-,-)
I2 = np.eye(2); Z2 = np.zeros((2, 2))
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)
g0 = np.block([[I2, Z2], [Z2, -I2]]).astype(complex)
gi = [np.block([[Z2, s], [-s, Z2]]) for s in (sx, sy, sz)]
gam = [g0]+gi
g5 = 1j*gam[0]@gam[1]@gam[2]@gam[3]
metric = np.diag([1., -1., -1., -1.])


def slash(p):
    return sum(metric[m, m]*p[m]*gam[m] for m in range(4))


def gluon_pol(direction, helicity):
    # massless gluon along +z or -z; epsilon^mu with lower-case components
    # epsilon(+z, h) = (0, 1, i h, 0)/sqrt2 ; for -z rotate by pi about y.
    if direction > 0:
        return np.array([0, 1, 1j*helicity, 0])/np.sqrt(2)
    return np.array([0, -1, 1j*helicity, 0])/np.sqrt(2)


def z_pol(p, m, kind):
    # three polarisation vectors orthogonal to a massive momentum p (p.p = m^2)
    pvec = p[1:]
    n = pvec/np.linalg.norm(pvec)
    a = np.cross(n, [0, 0, 1.]) if abs(n[2]) < 0.9 else np.cross(n, [1., 0, 0])
    a /= np.linalg.norm(a)
    b = np.cross(n, a)
    if kind == 0:
        return np.array([0, *a])
    if kind == 1:
        return np.array([0, *b])
    E = p[0]; P = np.linalg.norm(pvec)
    return np.array([P/m, *(E/m*n)])          # longitudinal


def vertex_Z(eps, v, a):
    # gamma^mu (v - a gamma5) eps*_mu ; eps is taken real or conjugated by the caller
    return slash(eps)@(v*np.eye(4)-a*g5)


def traces(l, p1, p2, p3, p4, V1, V2, V3, V4):
    # orientation A: tr[ V4 S(l+p4) V2 S(l+p1-p3) V3 S(l+p1) V1 S(l) ]
    A = np.trace(V4@slash(l+p4)@V2@slash(l+p1-p3)@V3@slash(l+p1)@V1@slash(l))
    # orientation B with l' = -(l+p1): tr[ V3 S(l'+p3) V2 S(l'+p1-p4) V4 S(l'+p1) V1 S(l') ]
    lp = -(l+p1)
    B = np.trace(V3@slash(lp+p3)@V2@slash(lp+p1-p4)@V4@slash(lp+p1)@V1@slash(lp))
    return A, B


def kinematics(beta, cos_th, m=1.0):
    E = m/np.sqrt(1-beta**2)              # energy of each boson, sqrt(s) = 2E
    P = E*beta
    st = np.sqrt(1-cos_th**2)
    p1 = np.array([E, 0, 0, E]); p2 = np.array([E, 0, 0, -E])
    p3 = np.array([E, P*st, 0, P*cos_th]); p4 = p1+p2-p3
    return p1, p2, p3, p4, m


worst = 0.0
for beta, cth in ((0.5, 0.3), (0.3, 0.9), (0.7, -0.4)):
    p1, p2, p3, p4, m = kinematics(beta, cth)
    for h1, h2 in ((+1, +1), (+1, -1), (-1, +1), (-1, -1)):
        V1 = slash(gluon_pol(+1, h1)); V2 = slash(gluon_pol(-1, h2))
        for e3 in range(3):
            for e4 in range(3):
                eps3 = z_pol(p3, m, e3).conj(); eps4 = z_pol(p4, m, e4).conj()
                res = {'vv': [], 'aa': [], 'va': [], 'scalar': []}
                for trial in range(4):
                    if trial < 2:   # generic off-shell loop momentum
                        l = rng.normal(size=4)
                    else:           # collinear to p1, as in the DPS region: l + p1 = x p1 + k_T with small k_T
                        x = (1+beta)/2
                        kT = np.array([0, *rng.normal(size=2)*1e-3, 0])
                        l = (x-1)*p1+kT
                    # Z couplings split into v^2, a^2 and va pieces by linearity
                    A_vv, B_vv = traces(l, p1, p2, p3, p4, V1, V2, vertex_Z(eps3, 1, 0), vertex_Z(eps4, 1, 0))
                    A_aa, B_aa = traces(l, p1, p2, p3, p4, V1, V2, vertex_Z(eps3, 0, 1), vertex_Z(eps4, 0, 1))
                    A_va, B_va = traces(l, p1, p2, p3, p4, V1, V2, vertex_Z(eps3, 1, 0), vertex_Z(eps4, 0, 1))
                    A_s, B_s = traces(l, p1, p2, p3, p4, V1, V2, np.eye(4), np.eye(4))
                    for key, (Ax, Bx) in (('vv', (A_vv, B_vv)), ('aa', (A_aa, B_aa)),
                                          ('va', (A_va, B_va)), ('scalar', (A_s, B_s))):
                        if abs(Ax) > 1e-9:
                            res[key].append(Bx/Ax)
                expected = {'vv': 1, 'aa': 1, 'va': -1, 'scalar': 1}
                for key, r in res.items():
                    if r:
                        worst = max(worst, np.max(np.abs(np.array(r)-expected[key])))

print('Crossed box, orientations A and B, 3 kinematic points x 4 gluon helicities x 9 Z polarisations,')
print('generic loop momenta and loop momenta collinear to p1:')
for key, exp in (('v^2 Z couplings', 1), ('a^2 Z couplings', 1), ('v a cross terms', -1), ('scalar (Yukawa) vertices', 1)):
    print(f'   {key:28s} expected B/A = {exp:+d}')
print(f'largest deviation of B/A from the expected value: {worst:.1e}')
assert worst < 1e-9
print('So the full amplitude contains the crossed box twice for the v^2 + a^2 and scalar parts:')
print('a coefficient extracted from the full gg -> ZZ or gg -> HH amplitude is 2 x the single-box value.')
