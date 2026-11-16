#!/usr/bin/env python3
# Paper Sec. 3.4: the one-loop class-B gg -> HH amplitude, read off from the mu_R
# dependence of the two-loop remainders, so that the two-loop logarithms can be
# divided by one-loop numbers taken from the same files.
#
# Renormalisation-group invariance of C_H^2 (M^(0) + a_s M^(1) + a_s^2 M^(2)) with
# C_H ~ a_s gives, as in P1 Eqs. (B.1), (B.3) and Eq. (2.23) of Banerjee et al.
# (arXiv:1809.05388v2),
#     dM^(1)/dln mu^2 = 2 beta0 M^(0),   dM^(2)/dln mu^2 = 3 beta0 M^(1) + c M^(0),
# with c a constant (beta1, H_g^(2), the C_HH contact counterterm). With
# L = ln(m_H^2/mu_R^2) the two-loop remainder is therefore A + B L + C L^2 with
#     B = -[3 beta0 M^(1)(mu_R = m_H) + c M^(0)],     C = 3 beta0^2 M^(0).
# In colour terms (3 beta0 = 11N - 2n_f, M^(1) = N m_A + n_f m_F):
#     B_{N^2}  = -11 m_A  + (tree),   B_{N n_f} = 2 m_A - 11 m_F + (tree),
#     B_{n_f^2} = 2 m_F   + (tree),   B_{n_f/N} = (tree),
#     C_{N^2} : C_{N n_f} : C_{n_f^2} : C_{n_f/N} = 121/3 : -44/3 : 4/3 : 0  (times M^(0)).
# "(tree)" terms are proportional to M^(0): constant in G for spin 0 and O(p_T^2)
# for spin 2, so they do not touch a G slope or a spin-two forward limit.
#
# Read off:
#   spin 0: the G slope of m_A (one-loop DPS coefficient per unit N), DPS value -2 pi F_S;
#   spin 2: the forward limit of m_A2, DPS value pi Kt (K_1 per unit N, paper Sec. 2);
#   m_F has neither (the n_f one-loop diagrams have no DPS pinch);
#   M^(0) from C, compared point by point with the class-B tree from Feynman rules.
# Then c13/c_1, c12/c_1 and c22/K use only numbers from the HH files.
#
#   python3 analyse_one_loop.py          (mpmath; a few seconds)
import json, sys
from fractions import Fraction
import mpmath as mp
from scan_tools import load, fit_model, L_parts, weights, rel, fmt, COLOURS

pi = mp.pi


def tree_hh(x, y):
    """Class-B tree M_1^(0)/m^2 from the t- and u-channel gluon exchanges, m_H = 1.

    Effective vertex (i C_H stripped) V^{ab}(r, t) = (r.t) g^{ab} - t^a r^b for gluons with
    incoming momenta r, t (paper Sec. 2.2); Feynman-gauge propagator -i g/k^2.
    The ++ helicity amplitude is i C_H^2 (S_t + S_u) and the T_1 contraction of the ++
    polarisations is -1, so M_1 = -i (S_t + S_u). Exact kinematics of the point (x, y), evaluated
    in mpmath at the working precision, so that the comparison with the coefficient of L^2 does
    not depend on the floating-point arithmetic of the machine."""
    x, y = mp.mpf(x), mp.mpf(y)
    s = (1+x)**2/x
    t = -y
    E = mp.sqrt(s)/2
    P3 = mp.sqrt(E**2-1)
    cth = (t-1+2*E**2)/(2*E*P3)             # t = m^2 - 2E^2 + 2 E P3 cos(theta)
    sth = mp.sqrt(max(mp.mpf(0), 1-cth**2))
    g = (1, -1, -1, -1)                     # metric, diagonal
    R = range(4)
    dot = lambda a, b: sum(g[i]*a[i]*b[i] for i in R)
    p1 = (E, 0, 0, E); p2 = (E, 0, 0, -E)
    p3 = (E, P3*sth, 0, P3*cth); p4 = tuple(p1[i]+p2[i]-p3[i] for i in R)
    h = 1/mp.sqrt(2)
    e1 = (0, h, mp.mpc(0, h), 0); e2 = (0, h, mp.mpc(0, -h), 0)

    def V(r, tt):   # upper indices; g^{ab} = g_ab numerically
        rt = dot(r, tt)
        return [[(rt*g[a] if a == b else 0)-tt[a]*r[b] for b in R] for a in R]

    S = 0
    for q in (p3, p4):           # the boson produced at the p1 vertex
        k = tuple(q[i]-p1[i] for i in R)     # exchanged gluon, incoming at that vertex
        A, B = V(p1, k), V(p2, tuple(-c for c in k))
        # e1_mu V^{mu a} g_ab V^{nu b} e2_nu / k^2
        S += sum(e1[m]*g[m]*A[m][a]*g[a]*B[n][a]*g[n]*e2[n] for m in R for a in R for n in R)/dot(k, k)
    return mp.mpc(0, -1)*S


def two_point_limit(recs, colour, part):
    """c from c + d p_T^2 through the two smallest deltas (spin-two forward limit)."""
    r1, r2 = recs[-2:]
    v1, v2 = L_parts(r1, colour)[part], L_parts(r2, colour)[part]
    p1, p2 = mp.mpf(r1['pT2']), mp.mpf(r2['pT2'])
    return (v2*p1-v1*p2)/(p1-p2)


def main(points='points', outjson='results_one_loop.json'):
    reg, _ = load(points)
    betas = sorted({b for (_, b) in reg}, key=Fraction)
    report = {}
    for beta in betas:
        if (1, beta) not in reg or (2, beta) not in reg:
            continue
        w, FS, Kt = weights(beta)
        r1, r2 = reg[(1, beta)], reg[(2, beta)]
        # free fits of the L^1 part: G^3, G^2 (spin 0) and G^2, G (spin 2) must be absent
        f1 = {c: fit_model(r1, c, 1, 3, 2) for c in COLOURS}
        f2 = {c: fit_model(r2, c, 1, 2, 2) for c in COLOURS}
        # with those absent: linear in G (spin 0) and constant (spin 2), plus p_T^2 corrections
        s0 = {c: (fit_model(r1, c, 1, 1, 2, 0), fit_model(r1, c, 1, 1, 2, 1)) for c in COLOURS}
        sA = -s0['N2'][0][1]/11
        sA_shift = abs(s0['N2'][0][1]-s0['N2'][1][1])/11
        sF = s0['Nf2'][0][1]/2
        kA = -two_point_limit(r2, 'N2', 1)/11
        kF = (2*kA-two_point_limit(r2, 'N_Nf', 1))/11
        kF2 = two_point_limit(r2, 'Nf2', 1)/2
        print(f'\nbeta = {beta}')
        print('  one loop, read from the L^1 part of the two-loop remainder (per unit N or n_f):')
        print(f'    spin 0, G slope of m_A  = {fmt(mp.re(sA), 16):>22s}   DPS: -2 pi F_S = {fmt(-2*pi*FS, 16)}'
              f'   rel {mp.nstr(rel(sA, -2*pi*FS), 2)} (shift {mp.nstr(sA_shift/abs(sA), 2)})')
        print(f'    spin 0, G slope of m_F  = {mp.nstr(abs(sF), 2):>22s}   (no DPS pinch: 0)')
        print(f'    spin 2, forward m_A2    = {fmt(mp.re(kA), 16):>22s}   DPS: pi Kt     = {fmt(pi*Kt, 16)}'
              f'   rel {mp.nstr(rel(kA, pi*Kt), 2)}')
        print(f'    spin 2, forward m_F2    = {mp.nstr(abs(kF), 2):>22s} (from N n_f), {mp.nstr(abs(kF2), 2)} (from n_f^2)   (0)')
        lead0 = abs(2*pi*FS)
        print('    largest |coefficient| of the powers that must be absent, relative to |c_1|:',
              mp.nstr(max(abs(f1[c][k]) for c in COLOURS for k in (2, 3))/lead0, 2), '(spin 0),',
              mp.nstr(max(abs(f2[c][k]) for c in COLOURS for k in (1, 2))/(pi*Kt), 2), '(spin 2)')
        # tree from the L^2 part, point by point, against the Feynman-rule tree
        worst_tree, worst_ratio = 0, 0
        for rec in r1:
            C = {c: L_parts(rec, c)[2] for c in COLOURS}
            M0 = 3*C['N2']/121
            fr = tree_hh(rec['x'], rec['y'])
            worst_tree = max(worst_tree, abs(M0-fr)/abs(fr))
            worst_ratio = max(worst_ratio, float(abs(C['N_Nf']/C['N2']+mp.mpf(4)/11)),
                              float(abs(C['Nf2']/C['N2']-mp.mpf(4)/121)), float(abs(C['Nf_over_N']/C['N2'])))
        print(f'  tree: 3 C_N2/121 vs Feynman rules (-i s/2 here), largest relative difference over the deltas {float(worst_tree):.1e};')
        print(f'        colour ratios of the L^2 part against -4/11, 4/121, 0: largest deviation {worst_ratio:.1e}')
        # ratio tests using only these files; leading two-loop coefficients from analyse_colour's model
        A13 = fit_model(r1, 'N2', 0, 3, 2)[3]
        A12 = fit_model(r1, 'N2', 0, 3, 2)[2]
        A12f = fit_model(r1, 'N_Nf', 0, 3, 2)[2]
        A22 = fit_model(r2, 'N2', 0, 2, 2)[2]
        A21 = fit_model(r2, 'N2', 0, 2, 2)[1]
        A21f = fit_model(r2, 'N_Nf', 0, 2, 2)[1]
        ratios = {'c13(N^2)/slope': (A13/sA, mp.mpf(-4)/3), 'c12(N^2)/slope': (A12/sA, mp.mpf(-22)/3),
                  'c12(N n_f)/slope': (A12f/sA, mp.mpf(4)/3), 'c22(N^2)/K': (A22/kA, mp.mpf(-4)),
                  'c21(N^2)/K': (A21/kA, mp.mpf(-32)/3), 'c21(N n_f)/K': (A21f/kA, mp.mpf(8)/3)}
        print('  ratios, every number from the HH files (expect -(4/3) C_A, -2 beta0 = -(22/3) N + (4/3) n_f, -4 C_A,')
        print('  4 C_A - 4 beta0 = -(32/3) N + (8/3) n_f):')
        for k, (v, e) in ratios.items():
            print(f'    {k:17s} = {fmt(mp.re(v), 14):>18s}   expect {fmt(e, 6):>9s}   rel {mp.nstr(rel(v, e), 2)}')
        report[beta] = {
            'one_loop': {'spin0_slope_mA': str(mp.re(sA)), 'DPS_minus_2piFS': str(-2*pi*FS),
                         'spin0_slope_rel_dev': str(rel(sA, -2*pi*FS)), 'spin0_slope_shift': str(sA_shift/abs(sA)),
                         'spin0_slope_mF_abs': str(abs(sF)),
                         'spin2_forward_mA2': str(mp.re(kA)), 'DPS_piKt': str(pi*Kt),
                         'spin2_forward_rel_dev': str(rel(kA, pi*Kt)),
                         'spin2_forward_mF2_abs': [str(abs(kF)), str(abs(kF2))]},
            'tree_max_rel_diff': float(worst_tree), 'L2_colour_ratio_max_dev': worst_ratio,
            'ratios': {k: {'value': str(mp.re(v)), 'expected': str(e), 'rel_dev': str(rel(v, e))} for k, (v, e) in ratios.items()}}
    json.dump(report, open(outjson, 'w'), indent=1)


if __name__ == '__main__':
    main(*sys.argv[1:])
