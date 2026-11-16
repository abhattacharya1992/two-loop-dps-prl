#!/usr/bin/env python3
"""Identify the 20 tensor structures T_j of arXiv:1503.08835 from the published projectors.

Usage: check_tensors.py INPUTS

The projectors satisfy P_j . T_k = delta_jk with P_j = sum_i B_ji T_i^dagger, so B is the
inverse of M_ik = T_i^dagger . T_k, the contraction with the d-dimensional polarisation sums
    gluons:  -g + (p1 p2 + p2 p1)/(p1.p2)    (gauge e1.p2 = e2.p1 = 0)
    bosons:  -g + p3 p3/p3^2,  -g + p4 p4/p4^2.
M is computed for a candidate basis with explicit vectors in D = 6 dimensions and inverted;
the result is compared with B read from gg_projectors-1.0.inc at d = 6.
Index order of each tensor: (mu, nu, rho, sigma) <-> (gluon 1, gluon 2, boson 3, boson 4).
"""
import itertools
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import mpmath as mp
from formeval import Evaluator

mp.mp.dps = 40
D = 6
G = [1] + [-1] * (D - 1)                      # metric diag


def vec(*c):
    return [mp.mpf(v) for v in c]


def dot(a, b):
    return sum(G[i] * a[i] * b[i] for i in range(D))


# kinematics in the centre-of-mass frame, beam along the last axis
s, ma2, mb2 = mp.mpf(10), mp.mpf("1.3"), mp.mpf("0.7")
rs = mp.sqrt(s)
E = rs / 2
p1 = vec(E, 0, 0, 0, 0, E)
p2 = vec(E, 0, 0, 0, 0, -E)
E3 = (s + ma2 - mb2) / (2 * rs)
k = mp.sqrt(E3**2 - ma2)
th, ph = mp.mpf("0.7"), mp.mpf("0.4")
n = [mp.sin(th) * mp.cos(ph), mp.sin(th) * mp.sin(ph), 0, 0, mp.cos(th)]
p3 = [E3] + [k * c for c in n]
p4 = [p1[i] + p2[i] - p3[i] for i in range(D)]
t = dot([p1[i] - p3[i] for i in range(D)], [p1[i] - p3[i] for i in range(D)])
u = dot([p2[i] - p3[i] for i in range(D)], [p2[i] - p3[i] for i in range(D)])
assert abs(dot(p4, p4) - mb2) < 1e-30

g = [[(G[i] if i == j else 0) for j in range(D)] for i in range(D)]   # g^{ij} (upper)
P = {"1": p1, "2": p2, "3": p3}


def basis():
    """Candidate tensors as functions T[mu][nu][rho][sigma] (all indices upper)."""
    def gg(a, b, c, d):         # g^{ab} g^{cd} with slots named 0..3
        return lambda I: g[I[a]][I[b]] * g[I[c]][I[d]]

    def gpp(a, b, c, pc, d, pd):  # g^{ab} p_pc^{c} p_pd^{d}
        return lambda I: g[I[a]][I[b]] * P[pc][I[c]] * P[pd][I[d]]

    def pppp(q):                  # p3^mu p3^nu p_q0^rho p_q1^sigma
        return lambda I: P["3"][I[0]] * P["3"][I[1]] * P[q[0]][I[2]] * P[q[1]][I[3]]
    mu, nu, rho, sig = 0, 1, 2, 3
    return [
        gg(mu, nu, rho, sig), gg(mu, rho, nu, sig), gg(mu, sig, nu, rho),
        gpp(mu, nu, rho, "1", sig, "1"), gpp(mu, nu, rho, "1", sig, "2"),
        gpp(mu, nu, rho, "2", sig, "1"), gpp(mu, nu, rho, "2", sig, "2"),
        gpp(mu, rho, nu, "3", sig, "1"), gpp(mu, rho, nu, "3", sig, "2"),
        gpp(mu, sig, nu, "3", rho, "1"), gpp(mu, sig, nu, "3", rho, "2"),
        gpp(nu, rho, mu, "3", sig, "1"), gpp(nu, rho, mu, "3", sig, "2"),
        gpp(nu, sig, mu, "3", rho, "1"), gpp(nu, sig, mu, "3", rho, "2"),
        gpp(rho, sig, mu, "3", nu, "3"),
        pppp("11"), pppp("12"), pppp("21"), pppp("22")]


def lower(i):
    return G[i]


def polsum_gluon(pa, pb):
    pab = dot(pa, pb)
    # -g_{ij} + (pa_i pb_j + pb_i pa_j)/(pa.pb), lower indices
    return [[-(G[i] if i == j else 0) + (G[i] * pa[i] * G[j] * pb[j] + G[i] * pb[i] * G[j] * pa[j]) / pab
             for j in range(D)] for i in range(D)]


def polsum_massive(p, m2):
    return [[-(G[i] if i == j else 0) + G[i] * p[i] * G[j] * p[j] / m2 for j in range(D)] for i in range(D)]


S1, S2, S3, S4 = polsum_gluon(p1, p2), polsum_gluon(p2, p1), polsum_massive(p3, ma2), polsum_massive(p4, mb2)

T = basis()
R = range(D)
IDX = list(itertools.product(R, R, R, R))
dense = [[f(I) for I in IDX] for f in T]       # flat arrays, index ((mu*D+nu)*D+rho)*D+sig


def apply(S, a, slot):
    """Contract index `slot` of the flat tensor a with S (S_{ij} a^{..j..})."""
    out = [mp.mpf(0)] * len(a)
    stride = D ** (3 - slot)
    for n, I in enumerate(IDX):
        i = I[slot]
        base = n - i * stride
        acc = mp.mpf(0)
        for j in range(D):
            if S[i][j]:
                v = a[base + j * stride]
                if v:
                    acc += S[i][j] * v
        out[n] = acc
    return out


lowered = []
for a in dense:
    b = apply(S1, a, 0)
    b = apply(S2, b, 1)
    b = apply(S3, b, 2)
    b = apply(S4, b, 3)
    lowered.append(b)
M = mp.matrix(20, 20)
for i in range(20):
    nz = [(n, v) for n, v in enumerate(dense[i]) if v]
    for k in range(20):
        M[i, k] = sum(v * lowered[k][n] for n, v in nz)
asym = max(abs(M[i, k] - M[k, i]) for i in range(20) for k in range(20))
print("asymmetry of M:", mp.nstr(asym, 3))
B_mine = M ** -1

# B from the file at d = 6
text = pathlib.Path(sys.argv[1], "gg_projectors-1.0.inc").read_text()
B_file = mp.matrix(20, 20)
blocks = re.findall(r"id P(\d+) =(.*?);", text, re.S)
fn = lambda name, args: (1 / args[0]) if name == "Den" else args[0]
for j, body in blocks:
    j = int(j)
    for i in range(1, 21):
        syms = {"d": mp.mpf(D), "t": t, "u": u, "ma2": ma2, "mb2": mb2}
        for k in range(1, 21):
            syms[f"T{k}dag"] = mp.mpf(1 if k == i else 0)
        B_file[j - 1, i - 1] = Evaluator(syms, fn).evaluate(body)

scale = max(abs(B_file[a, b]) for a in range(20) for b in range(20))
diff = max(abs(B_file[a, b] - B_mine[a, b]) for a in range(20) for b in range(20))
print(f"max |B_file - B_candidate| / max |B| = {mp.nstr(diff / scale, 3)}")
if not all(mp.isfinite(v) for v in [*B_file,*B_mine,asym,diff,scale]) or scale<=0:
    raise ValueError('Nonfinite or empty tensor comparison')
if diff / scale > mp.mpf('1e-20'):
    for a in range(20):
        row = [mp.nstr(B_file[a, b] / B_mine[a, b], 4) if B_mine[a, b] != 0 else "-" for b in range(20)]
        print(a + 1, row)
    raise ValueError('Tensor projectors disagree beyond 1e-20')
if len(blocks)!=20 or {int(j) for j,_ in blocks}!=set(range(1,21)):
    raise ValueError('Incomplete tensor projector set')
