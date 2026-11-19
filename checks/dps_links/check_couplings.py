"""Coupling structure of the massless quark loop in g g -> V1 V2 (integrand level, exact).

Supports the statement of the Letter ("Up to electroweak couplings, these amplitudes describe
gamma* gamma*, Z gamma*, ZZ and W+ W- pairs") and of the Supplement (Sec. "Quark splittings and two-loop
gg -> ZZ", last paragraph of "Prediction"). Standard library only; run: python3 check_couplings.py

The closed quark loop has the four vertices g(p1), V1(p3), g(p2), V2(p4) in this cyclic order (the
order of the DPS box). The vector bosons couple with gamma^mu (L P_L + R P_R). The check is on the
numerator of the integrand, for arbitrary loop momentum and arbitrary vectors at the four vertices,
with exact rational arithmetic (Gaussian rationals, explicit 4x4 Dirac matrices).

Statements tested:
 (1) mixed-chirality terms vanish:               T_LR = T_RL = 0;
 (2) T_LL + T_RR = T_VV  and  T_RR - T_LL = T_5  (T_5: one gamma_5), so that
     T = 1/2 (L1 L2 + R1 R2) T_VV + 1/2 (R1 R2 - L1 L2) T_5;
 (3) axial-axial equals vector-vector:            T_AA = T_VV;
 (4) reversing the fermion flow:                  T_VV(rev) = T_VV,  T_5(rev) = -T_5,
     so the sum of the two orientations is (L1 L2 + R1 R2) T_VV: one coupling factor times a
     function that does not depend on the couplings;
 (5) (1)-(3) with a gluon attached twice to the loop (a two-loop numerator).

Not tested here: the integration, the reduction to master integrals, and diagrams in which the two
vector bosons couple to different quark loops (they vanish for massless quarks, arXiv:1503.08835).
"""
from fractions import Fraction as Fr
import random

class C:
    """Gaussian rational"""
    __slots__ = ('r', 'i')
    def __init__(s, r=0, i=0): s.r, s.i = Fr(r), Fr(i)
    def __add__(s, o): o = o if isinstance(o, C) else C(o); return C(s.r + o.r, s.i + o.i)
    __radd__ = __add__
    def __sub__(s, o): o = o if isinstance(o, C) else C(o); return C(s.r - o.r, s.i - o.i)
    def __neg__(s): return C(-s.r, -s.i)
    def __mul__(s, o): o = o if isinstance(o, C) else C(o); return C(s.r*o.r - s.i*o.i, s.r*o.i + s.i*o.r)
    __rmul__ = __mul__
    def __eq__(s, o): o = o if isinstance(o, C) else C(o); return s.r == o.r and s.i == o.i
    def iszero(s): return s.r == 0 and s.i == 0
    def __repr__(s): return f'({s.r}+{s.i}i)'

Z, O, I = C(0), C(1), C(0, 1)
def mat(rows): return [[x if isinstance(x, C) else C(x) for x in row] for row in rows]
def mm(a, b): return [[sum((a[i][k]*b[k][j] for k in range(4)), C(0)) for j in range(4)] for i in range(4)]
def madd(a, b, f=1): return [[a[i][j] + f*b[i][j] for j in range(4)] for i in range(4)]
def msc(a, f): return [[f*a[i][j] for j in range(4)] for i in range(4)]
def tr(a): return sum((a[i][i] for i in range(4)), C(0))
def prod(*ms):
    out = ms[0]
    for m in ms[1:]: out = mm(out, m)
    return out
ID = mat([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
# Dirac representation
g0 = mat([[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,-1]])
g1 = mat([[0,0,0,1],[0,0,1,0],[0,-1,0,0],[-1,0,0,0]])
g2 = mat([[0,0,0,-I],[0,0,I,0],[0,I,0,0],[-I,0,0,0]])
g3 = mat([[0,0,1,0],[0,0,0,-1],[-1,0,0,0],[0,1,0,0]])
G = [g0, g1, g2, g3]
g5 = msc(prod(g0, g1, g2, g3), I)
PL, PR = msc(madd(ID, g5, -1), Fr(1, 2)), msc(madd(ID, g5, 1), Fr(1, 2))
# Clifford algebra and gamma_5 properties (sanity)
metric = [1, -1, -1, -1]
for m in range(4):
    for n in range(4):
        anti = madd(mm(G[m], G[n]), mm(G[n], G[m]))
        assert anti == msc(ID, 2*metric[m] if m == n else 0)
    assert madd(mm(G[m], g5), mm(g5, G[m])) == msc(ID, 0)
assert mm(g5, g5) == ID

def slash(p): return madd(madd(madd(msc(g0, p[0]), msc(g1, -p[1])), msc(g2, -p[2])), msc(g3, -p[3]))
def vec(rng): return [Fr(rng.randint(-9, 9), rng.randint(1, 7)) for _ in range(4)]

def run(seed):
    rng = random.Random(seed)
    e1, e2, e3, e4 = (vec(rng) for _ in range(4))          # arbitrary vectors at the four vertices
    p1, p2, p3 = vec(rng), vec(rng), vec(rng)              # arbitrary external momenta (no on-shell condition needed)
    k = vec(rng)                                            # arbitrary loop momentum
    sub = lambda a, b: [x - y for x, y in zip(a, b)]
    add = lambda a, b: [x + y for x, y in zip(a, b)]
    q1 = k; q2 = sub(k, p3); q3 = add(q2, p2); q4 = sub(k, p1)   # massless quark momenta around the loop
    S = [slash(q) for q in (q1, q2, q3, q4)]
    E1, E2, E3, E4 = (slash(e) for e in (e1, e2, e3, e4))
    def T(v1, v2):          # orientation 1: g1 -> V1 -> g2 -> V2 (read against the fermion flow)
        return tr(prod(E1, S[3], mm(E4, v2), S[2], E2, S[1], mm(E3, v1), S[0]))
    def Trev(v1, v2):       # orientation 2: opposite fermion flow, same denominators
        m = lambda a: msc(a, -1)
        return tr(prod(E1, m(S[0]), mm(E3, v1), m(S[1]), E2, m(S[2]), mm(E4, v2), m(S[3])))
    res = {}
    for name, f in (('or1', T), ('or2', Trev)):
        res[name] = dict(VV=f(ID, ID), V5=f(g5, ID), FV=f(ID, g5), AA=f(g5, g5),
                         LL=f(PL, PL), RR=f(PR, PR), LR=f(PL, PR), RL=f(PR, PL))
    a = res['or1']; b = res['or2']
    checks = {
        '(1) T_LR = T_RL = 0 (both orientations)': a['LR'].iszero() and a['RL'].iszero() and b['LR'].iszero() and b['RL'].iszero(),
        '(2) T_LL + T_RR = T_VV': a['LL'] + a['RR'] == a['VV'],
        '(2) T_RR - T_LL = T_5 (gamma_5 at V1) = T_5 (gamma_5 at V2)': a['RR'] - a['LL'] == a['V5'] and a['V5'] == a['FV'],
        '(3) T_AA = T_VV': a['AA'] == a['VV'],
        '(4) T_VV(reversed) = T_VV': b['VV'] == a['VV'],
        '(4) T_5(reversed) = -T_5': (b['V5'] + a['V5']).iszero(),
        'T_VV and T_5 are both nonzero at this point': not a['VV'].iszero() and not a['V5'].iszero(),
    }
    # general couplings: sum of both orientations equals (L1 L2 + R1 R2) T_VV
    L1, R1, L2, R2 = (Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(4))
    c1 = madd(msc(PL, L1), msc(PR, R1)); c2 = madd(msc(PL, L2), msc(PR, R2))
    checks['(4) sum of orientations = (L1 L2 + R1 R2) T_VV for random couplings'] = T(c1, c2) + Trev(c1, c2) == (L1*L2 + R1*R2)*a['VV']
    # (5) two-loop numerator: a gluon exchanged inside the loop between the lines q1 and q3
    l = vec(rng)
    def T2(v1, v2):
        out = C(0)
        for mu in range(4):   # gamma^mu ... gamma_mu, metric sign
            out = out + metric[mu]*tr(prod(E1, S[3], mm(E4, v2), slash(add(q3, l)), G[mu], S[2], E2, S[1], mm(E3, v1),
                                           S[0], G[mu], slash(add(q1, l))))
        return out
    checks['(5) two-loop numerator: T_LR = T_RL = 0'] = T2(PL, PR).iszero() and T2(PR, PL).iszero()
    checks['(5) two-loop numerator: T_AA = T_VV'] = T2(g5, g5) == T2(ID, ID)
    checks['(5) two-loop numerator: T_LL + T_RR = T_VV'] = T2(PL, PL) + T2(PR, PR) == T2(ID, ID)
    return checks

if __name__ == '__main__':
    total = {}
    n = 6
    for seed in range(1, n + 1):
        for key, ok in run(seed).items():
            total[key] = total.get(key, 0) + (1 if ok else 0)
    for key, v in total.items():
        print(f'{"PASS" if v == n else "FAIL"}  {v}/{n}  {key}')
    assert all(v == n for v in total.values())
    print('all exact (rational arithmetic), 6 random kinematic points')
