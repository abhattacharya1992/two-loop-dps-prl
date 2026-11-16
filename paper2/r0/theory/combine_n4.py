"""The one-loop DPS kernel of planar N = 4 super Yang-Mills at mu = mu_b from the regions, for gluon splittings with the
two momentum fractions x1, x2 (the fractions of the two daughters that enter hard vertex 1), in the four-dimensional-
helicity scheme, in units a_s C_A = a/2 (a = g^2 N_c/(8 pi^2), the coupling of Dixon and Esterlis):

    K^(1)/K_0 = J1(x1) + J2(x2) + S + Sigma(x1) + Sigma(x2)      (Paper 2, eq. (F.1), combine_gg.py)

with the N = 4 field content in the loops: the gluon and ghost graphs of gg_diagrams.py (colour C_A), the fermion graphs
with T_F n_f -> 2 C_A (four adjoint Weyl fermions: an adjoint Weyl fermion is N_c fundamental Dirac flavours, T_F N_c =
C_A/2 each), and the scalar graphs of n4_diagrams.py with (3 + eps) complex adjoint scalars (the six real scalars, and
the 2 eps epsilon-scalars of dimensional reduction, which turn conventional dimensional regularisation into the
four-dimensional-helicity scheme of the BDS normalisation; check_fdh.py).  The daughter self-energies: (5-3eps)/(3-2eps)
- 8(1-eps)/(3-2eps) - (3+eps)/(3-2eps) = -2 exactly.  Expected poles 2/eps^2 + 2 L_*/eps with L_* =
-ln[x1 (1-x1) x2 (1-x2)] (beta_0 = 0); the finite part is r_0^{N=4}(x1, x2)/C_A, the same for nu = 0 and nu = 2
(no rational helicity term: 1 - 4 + 3 = 0).

    python3 combine_n4.py 3/10 7/10        (needs gg_laurent_x3o10.json, gg_laurent_x7o10.json, n4_laurent_x*.json)
"""
import sys, json, os
import numpy as np
import mpmath as mp
import combine_gg as C
from laurent_inputs import load_laurent

STRUCT = C.STRUCT
Z2 = C.Z2
WEIGHT = {'CA': 1.0, 'nf': 2.0}          # T_F n_f -> 2 C_A
NS = 3.0                                 # complex adjoint scalars at eps^0; the epsilon-scalars add eps x (the same graphs)


def load_n4(tag, suffix=''):
    """J1 at x = tag: {structure: {eps order: value}} (eta^0 terms) and {structure: {eps order}} of the 1/eta terms,
    weighted with the N = 4 field content; tau = tree structure"""
    x = tag.split('_', 1)[0].replace('o', '/')
    d = load_laurent(f'gg_laurent_x{tag}{suffix}.json', 'gg', expected_x=x)
    n = load_laurent(f'n4_laurent_x{tag}{suffix}.json', 'n4', expected_x=x)
    tree = {s: complex(*d['tree'][s]) for s in STRUCT}
    norm = tree['S1']
    tau = {s: tree[s]/norm for s in STRUCT}
    tot = {s: {} for s in STRUCT}

    def add(dd, w, shift):
        for name, g in dd.items():
            for s in STRUCT:
                for key, (re, im) in g['laurent'][s].items():
                    i, j = map(int, key.split(','))
                    tot[s][(i, j + shift)] = tot[s].get((i, j + shift), 0) + w*complex(re, im)
    for name, g in d['diagrams'].items():
        for col, w in WEIGHT.items():
            c = complex(*g['colour'][col])
            if c != 0:
                add({name: g}, w*c, 0)
    for name, g in n['diagrams'].items():
        c = complex(*g['colour']['ns'])
        if c != 0:
            add({name: g}, NS*c, 0)            # the three scalars
            add({name: g}, c, 1)               # eps x the same graphs: the epsilon-scalars
    return tot, tau


def self_energy_n4(xx):
    return C.series(lambda ep: C.ftf(ep)*(-2)*mp.gamma(ep)*mp.gamma(1 - ep)**2/mp.gamma(2 - 2*ep)*((1 - xx)**ep + xx**ep))


def kernel_n4(x1, x2, suffix=''):
    """x1, x2 as sympy-like fractions 'p/q'; returns {nu: {...}} as combine_gg.kernel"""
    t1 = x1.replace('/', 'o'); t2 = x2.replace('/', 'o')
    J1, tau1 = load_n4(t1, suffix)
    J2, tau2 = load_n4(t2, suffix)
    xx1 = eval(x1); xx2 = eval(x2)
    se1, se2 = self_energy_n4(xx1), self_energy_n4(xx2)
    out = {}
    for nu in (0, 2):
        res, eta = {}, {}
        for j in (-2, -1, 0):
            A1 = C.to_b({s: J1[s].get((0, j), 0) for s in STRUCT})
            A2 = C.to_b({s: J2[s].get((0, j), 0) for s in STRUCT})
            res[j] = C.ratio(nu, A1, tau1, A2, tau2) + C.SOFT_0[j] + se1[j] + se2[j]
            E1 = {s: J1[s].get((-1, j), 0) for s in STRUCT}
            E2 = {s: J2[s].get((-1, j), 0) for s in STRUCT}
            eta[j] = C.ratio(nu, E1, tau1, E2, tau2) + C.SOFT_ETA[j]
        out[nu] = {'eps^-2': res[-2].real, 'eps^-1': res[-1].real, 'eps^0': res[0].real,
                   'max |1/eta coefficient|': max(abs(v) for v in eta.values()),
                   'max |imag|': max(abs(v.imag) for v in res.values())}
    return out


def expected_poles(xx1, xx2):
    return (2.0, -2*np.log(xx1*(1 - xx1)) - 2*np.log(xx2*(1 - xx2)))


def r0_guess(xx1, xx2):
    """the structure the hexagon requires (n4_hexagon.py): -5 zeta_2 - 2 ln x1 ln(1-x1) - 2 ln x2 ln(1-x2)"""
    return -5*Z2 - 2*np.log(xx1)*np.log(1 - xx1) - 2*np.log(xx2)*np.log(1 - xx2)


def main():
    x1 = sys.argv[1] if len(sys.argv) > 1 else '3/10'
    x2 = sys.argv[2] if len(sys.argv) > 2 else '7/10'
    out = kernel_n4(x1, x2)
    ep = expected_poles(eval(x1), eval(x2))
    for nu, r in out.items():
        print(f'x1 = {x1}, x2 = {x2}, nu={nu}: 1/eta {r["max |1/eta coefficient|"]:.1e} | eps^-2 {r["eps^-2"]:+.9f} (exp {ep[0]:+.4f})'
              f'  eps^-1 {r["eps^-1"]:+.9f} (exp {ep[1]:+.9f})  eps^0 {r["eps^0"]:+.9f}  '
              f'[-5 z2 - 2 ln x ln(1-x) terms: {r0_guess(eval(x1), eval(x2)):+.9f}]  |imag| {r["max |imag|"]:.1e}')


if __name__ == '__main__':
    main()
