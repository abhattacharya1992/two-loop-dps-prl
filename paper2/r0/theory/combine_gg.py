"""The one-loop kernel for gluon splittings (gg -> HH) at mu = mu_b, from J1 at x and at 1 - x, the self-energies of
the daughters and the soft function, for the two components nu = 0 and nu = 2; against the closed form of the data.

    python3 combine_gg.py 3/10 [CA|nf]        (needs gg_laurent_x3o10.json and gg_laurent_x7o10.json)

Units p1^+ = p2^- = mu = nu = 1, b = b0, s = 1, m1^2 = m2^2 = x(1-x).  J1 is stored per structure (S1, S2, S3, T4 of
gg_driver.py) in units of the tree coefficient of S1, which is the same at x and 1 - x.  b space: S1-S3 keep their
form; T4 -> (S1 + S2 + S3)/2 - T4 at eps^0 (its spin-3 part transforms with the opposite sign; T4 has no poles).
The n-bar collinear function is J1 at 1 - x (same tensor form, with eps_2 and b_r = -b; the overall sign cancels).
Kernel: sum_ij J1^ij(eps1, b) J2^ij(eps2, b) = alpha (eps1.eps2) b^2 + beta (eps1.b)(eps2.b): nu = 0 is alpha +
beta/2, nu = 2 is beta.  Coupling renormalisation: -beta0/eps, pure pole, not included (so the expected single pole
contains +beta0 relative to 2 I_K).
"""
import sys, json
import numpy as np
import sympy as sp
import mpmath as mp
from laurent_inputs import load_laurent

STRUCT = ('S1', 'S2', 'S3', 'T4')
Z2 = float(mp.zeta(2))


def series(f, orders=(-2, -1, 0), r=0.2, n=32):
    th = 2*np.pi*(np.arange(n) + 0.5)/n
    e = r*np.exp(1j*th)
    v = np.array([complex(f(mp.mpc(ee))) for ee in e])
    return {j: np.mean(v*e**(-j)) for j in orders}


def ftf(ep):
    return mp.exp(-mp.euler*ep)*mp.gamma(1 - 2*ep)/(mp.gamma(1 + ep)*mp.gamma(1 - ep))


def self_energy(xx, col):
    """the two daughters of J1(xx), b space at L_b = 0"""
    if col == 'CA':
        c = lambda ep: (5 - 3*ep)/(3 - 2*ep)
    else:
        c = lambda ep: -4*(1 - ep)/(3 - 2*ep)
    return series(lambda ep: ftf(ep)*c(ep)*mp.gamma(ep)*mp.gamma(1 - ep)**2/mp.gamma(2 - 2*ep)
                  * ((1 - xx)**ep + xx**ep))


SOFT_ETA = series(lambda ep: 2*mp.exp(-mp.euler*ep)*mp.gamma(-ep)*2)
SOFT_0 = series(lambda ep: 2*mp.exp(-mp.euler*ep)*mp.gamma(-ep)*(-mp.euler - mp.digamma(-ep)))


def contraction_table():
    """C[nu][(s, t)] for the b-space structures of the two sides (2 dimensions)"""
    e1 = sp.symbols('e1x e1y'); e2 = sp.symbols('e2x e2y'); bb = sp.symbols('bx by')
    dot = lambda a, c: a[0]*c[0] + a[1]*c[1]
    b2 = dot(bb, bb)

    def tens(s, e):
        out = [[0, 0], [0, 0]]
        for i in range(2):
            for j in range(2):
                if s == 'S1':
                    out[i][j] = dot(e, bb)*(1 if i == j else 0)
                elif s == 'S2':
                    out[i][j] = e[i]*bb[j]
                elif s == 'S3':
                    out[i][j] = e[j]*bb[i]
                else:
                    out[i][j] = dot(e, bb)*bb[i]*bb[j]/b2
        return out
    al, be = sp.symbols('alpha beta')
    C = {0: {}, 2: {}}
    for s in STRUCT:
        for t in STRUCT:
            A, B = tens(s, e1), tens(t, e2)
            val = sp.expand(sp.cancel(sum(A[i][j]*B[i][j] for i in range(2) for j in range(2))*b2**2))
            # match val = [alpha (e1.e2) b^2 + beta (e1.b)(e2.b)] b^4
            ansatz = sp.expand((al*dot(e1, e2)*b2 + be*dot(e1, bb)*dot(e2, bb))*b2**2)
            sol = sp.solve(sp.Poly(val - ansatz, *e1, *e2, *bb).coeffs(), [al, be], dict=True)
            assert len(sol) == 1, (s, t)
            a_, b_ = sol[0][al], sol[0][be]
            C[0][(s, t)] = float(a_ + b_/2)
            C[2][(s, t)] = float(b_)
    return C


CT = contraction_table()


def load(tag, col):
    d = load_laurent(f'gg_laurent_x{tag}.json', 'gg', expected_x=tag.split('_', 1)[0].replace('o', '/'))
    tree = {s: complex(*d['tree'][s]) for s in STRUCT}
    tot = {s: {} for s in STRUCT}
    for name, dd in d['diagrams'].items():
        w = complex(*dd['colour'][col])
        if w == 0:
            continue
        for s in STRUCT:
            for key, (re, im) in dd['laurent'][s].items():
                tot[s][key] = tot[s].get(key, 0) + w*complex(re, im)
    norm = tree['S1']
    tau = {s: tree[s]/norm for s in STRUCT}
    return tot, tau


def to_b(a):
    """k-space coefficients (S1, S2, S3, T4) -> b-space structures (T4 at eps^0)"""
    return {'S1': a['S1'] + a['T4']/2, 'S2': a['S2'] + a['T4']/2, 'S3': a['S3'] + a['T4']/2, 'T4': -a['T4']}


def ratio(nu, A1, tau1, A2, tau2):
    """[<A1, tau2> + <tau1, A2>]/<tau1, tau2> for coefficient dicts"""
    num = sum(A1[s]*tau2[t]*CT[nu][(s, t)] + tau1[s]*A2[t]*CT[nu][(s, t)] for s in STRUCT for t in STRUCT)
    den = sum(tau1[s]*tau2[t]*CT[nu][(s, t)] for s in STRUCT for t in STRUCT)
    return num/den


def r0_closed(xx, col, nu):
    F = 1/(xx*(1 - xx)) + 2*xx*(1 - xx) if nu == 0 else (xx**4 + (1 - xx)**4)/(xx*(1 - xx))
    if col == 'CA':
        return 134/9 - 5*Z2 - 4*np.log(xx)*np.log(1 - xx) + (2/3)/F
    return -40/9 - (4/3)/F


def kernel(num, den, col, suffix=''):
    """{nu: {'poles': (c-2, c-1), 'eps^0': ., '1/eta': ., 'imag': .}} for one colour part"""
    xx = num/den
    tx, ty = f'{num}o{den}{suffix}', f'{den - num}o{den}{suffix}'
    J1, tau1 = load(tx, col)
    J2, tau2 = load(ty, col)
    se = {j: self_energy(xx, col)[j] + self_energy(1 - xx, col)[j] for j in (-2, -1, 0)}
    dc = 1 if col == 'CA' else 0
    out = {}
    for nu in (0, 2):
        res, eta = {}, {}
        for j in (-2, -1, 0):
            A1 = to_b({s: J1[s].get(f'0,{j}', 0) for s in STRUCT})
            A2 = to_b({s: J2[s].get(f'0,{j}', 0) for s in STRUCT})
            res[j] = ratio(nu, A1, tau1, A2, tau2) + dc*SOFT_0[j] + se[j]
            E1 = {s: J1[s].get(f'-1,{j}', 0) for s in STRUCT}
            E2 = {s: J2[s].get(f'-1,{j}', 0) for s in STRUCT}
            eta[j] = ratio(nu, E1, tau1, E2, tau2) + dc*SOFT_ETA[j]
        out[nu] = {'eps^-2': res[-2].real, 'eps^-1': res[-1].real, 'eps^0': res[0].real,
                   'max |1/eta coefficient|': max(abs(v) for v in eta.values()),
                   'max |imag|': max(abs(v.imag) for v in res.values())}
    return out


def expected_poles(xx, col):
    Lx = np.log(xx*(1 - xx))
    return (2.0, -4*Lx + 22/3) if col == 'CA' else (0.0, -8/3)


def main():
    num, den = map(int, (sys.argv[1] if len(sys.argv) > 1 else '3/10').split('/'))
    cols = [sys.argv[2]] if len(sys.argv) > 2 else ['CA', 'nf']
    xx = num/den
    for col in cols:
        out = kernel(num, den, col)
        ep = expected_poles(xx, col)
        for nu, r in out.items():
            print(f'{col} nu={nu}: 1/eta {r["max |1/eta coefficient|"]:.1e} | eps^-2 {r["eps^-2"]:+.9f} '
                  f'(exp {ep[0]:+.4f})  eps^-1 {r["eps^-1"]:+.9f} (exp {ep[1]:+.9f})  eps^0 {r["eps^0"]:+.9f}  '
                  f'data {r0_closed(xx, col, nu):+.9f}  diff {r["eps^0"] - r0_closed(xx, col, nu):+.1e}  '
                  f'|imag| {r["max |imag|"]:.1e}')


if __name__ == '__main__':
    main()
