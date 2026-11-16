"""kappa_1 of planar N = 4 super Yang-Mills from the regions (Paper 2, Sec. 5.1: the one-loop constant of the DPS kernel
at its own scale in the self-crossing limit of the six-gluon amplitude of Dixon and Esterlis, fitted there to the
two-loop single logarithm as kappa_1 = zeta_2).

The singular terms of E = A_6/A_6^{BDS-like} (DE eq. (2.12)) in b space, at two loops and in units of the one-loop DPS
term a K_0(b):  E^(2)/(a K_0) = h_1 + h_2 - M6hat/2 + K^(1)/K_0,  with
    h_k = M_4^(1)(s_k, t_k) = -(1/eps^2)[(mu^2/-s_k)^eps + (mu^2/-t_k)^eps] + ln^2(s_k/t_k)/2 + 4 zeta_2   (BDS normalisation,
          a = g^2 N_c/(8 pi^2) (4 pi e^-gamma)^eps, FDH scheme), the one-loop four-gluon amplitudes of the two hard
          2 -> 2 subprocesses, with their adjacent invariants in the planar ordering (6',1,2,3') and (3'',4,5,6''),
    M6hat = 2 (M_6^{1-loop} + Y) = sum_i [-(1/eps^2)(mu^2/S_i)^eps] + Q(ln S_i) + 9 zeta_2,  S_i = -s_{i,i+1},
          Q = (1/2) sum_i l_i^2 - sum_i l_i l_{i+1} + (1/2) sum_i l_i l_{i+3}  (bdslike_check.py: from the one-mass and
          two-mass-easy box functions of the one-loop MHV amplitude and Y = sum Li_2(1-u_i) + (1/2) sum ln^2 u_i, the
          three-particle invariants cancel; the BDS-like normalisation depends on two-particle invariants only),
    K^(1)/K_0 = a_s C_A [2/eps^2 + 2 L_*/eps + r_0/C_A] = a [1/eps^2 + L_*/eps + r_0/(2 C_A)] at mu = mu_b, L_* = -X,
          r_0 = r_0^{N=4}(x1, x2) from the regions (combine_n4.py).
In 3 -> 3 kinematics (DE: particles 1, 3, 5 incoming) every invariant above is negative, so everything is real and
X = ln(|Q_*^2| b^2/b_0^2).  The self-crossing kinematics (DE eqs. (2.17)-(2.22)): s_12 = (1-x)(1-y) s, s_45 = x y s,
s_23 = t_1/(1-x), s_61 = t_1/(1-y), s_34 = t_2/x, s_56 = t_2/y, with t_1 = (p_2 + (1-x)p_3)^2 and t_2 = (p_4 + x p_3)^2
the t channels of the two subprocesses, Q_*^2 = s x y (1-x)(1-y) (= m_1^2 m_2^2/s, Paper 2).

Part 1 (symbolic): the poles cancel, and the finite part is
    B(X; x, y) = -X^2/2 + ln x ln(1-x) + ln y ln(1-y) + 7 zeta_2/2,
independent of t_1, t_2 (the scattering angles) and of the relative azimuth (the cross ratio v of DE): -X^2/2 is the
cusp term -gamma_K^(1)/8 X^2, there is no term linear in X (B^(1) = 0), and the x, y dependence has the form of the
kernel's -2 ln x ln(1-x) terms.
Part 2 (numerical): kappa_1 = B|_{X^0} + r_0^{N=4}(x1, x2)/(2 C_A) for all pairs (x1, x2) of the computed fractions,
nu = 0 and nu = 2; must be the same number for all of them: zeta_2 if the structure is right.

    python3 n4_hexagon.py          (writes kappa1_n4.json and prints the table)
"""
import json, os, itertools
import numpy as np
import sympy as sp
import combine_n4 as N4

Z2 = float(sp.zeta(2))


def symbolic():
    eps, X, lx, l1x, ly, l1y, lT1, lT2 = sp.symbols('eps X lnx ln1mx lny ln1my lamT1 lamT2', real=True)
    z2 = sp.zeta(2)
    # l(A) = ln(A/mu_b^2) = X + ln(A/|Q_*^2|)
    lS = {1: X - lx - ly, 2: X + lT1 - l1x, 3: X + lT2 - lx, 4: X - l1x - l1y, 5: X + lT2 - ly, 6: X + lT1 - l1y}
    lT = {1: X + lT1, 2: X + lT2}
    pole = lambda l: -sp.exp(-eps*l)/eps**2                     # -(1/eps^2)(mu^2/A)^eps with ln(A/mu^2) = l
    h = {k: pole(lS[s]) + pole(lT[k]) + sp.Rational(1, 2)*(lS[s] - lT[k])**2 + 4*z2 for k, s in ((1, 1), (2, 4))}
    Q = (sp.Rational(1, 2)*sum(lS[i]**2 for i in lS) - sum(lS[i]*lS[i % 6 + 1] for i in lS)
         + sp.Rational(1, 2)*sum(lS[i]*lS[(i + 2) % 6 + 1] for i in lS))
    M6hat = sum(pole(lS[i]) for i in lS) + Q + 9*z2
    Lstar = -X
    kernel_poles = 1/eps**2 + Lstar/eps
    total = h[1] + h[2] - M6hat/2 + kernel_poles
    ser = sp.series(total, eps, 0, 1).removeO()
    ser = sp.expand(ser)
    poles = {j: sp.simplify(ser.coeff(eps, j)) for j in (-2, -1)}
    B = sp.expand(ser.coeff(eps, 0))
    B_expected = -X**2/2 + lx*l1x + ly*l1y + sp.Rational(7, 2)*z2
    assert poles[-2] == 0 and poles[-1] == 0, poles
    assert sp.simplify(B - B_expected) == 0, sp.simplify(B - B_expected)
    return {'poles': 'cancel', 'B': str(B_expected), 'X^2': str(sp.Poly(B, X).coeff_monomial(X**2)),
            'X^1': str(sp.Poly(B, X).coeff_monomial(X)), 'depends on t1, t2': bool(B.has(lT1) or B.has(lT2))}


def numerical():
    tags = ['1/10', '3/10', '1/2', '7/10', '9/10']
    avail = [t for t in tags if os.path.exists(f'gg_laurent_x{t.replace("/", "o")}.json')
             and os.path.exists(f'n4_laurent_x{t.replace("/", "o")}.json')]
    rows = []
    for x1, x2 in itertools.combinations_with_replacement(avail, 2):
        res = N4.kernel_n4(x1, x2)
        xx1, xx2 = eval(x1), eval(x2)
        ep = N4.expected_poles(xx1, xx2)
        rec = {'x1': x1, 'x2': x2, 'expected poles (a_s C_A units)': {'eps^-2': ep[0], 'eps^-1': ep[1]}}
        for nu, r in res.items():
            r0 = r['eps^0']                                        # r_0^{N=4}/C_A in a_s units
            kappa = 3.5*Z2 + np.log(xx1)*np.log(1 - xx1) + np.log(xx2)*np.log(1 - xx2) + r0/2
            rec[f'nu={nu}'] = dict(r, kappa_1=kappa, **{'kappa_1 - zeta_2': kappa - Z2,
                                                        'r0 - [-5 z2 - 2 ln x ln(1-x) terms]': r0 - N4.r0_guess(xx1, xx2)})
        rows.append(rec)
    return rows


def main():
    sym = symbolic()
    print('symbolic: poles', sym['poles'], '; B =', sym['B'], '; X^2 coefficient', sym['X^2'], '; X coefficient', sym['X^1'],
          '; depends on the scattering angles:', sym['depends on t1, t2'])
    rows = numerical()
    print(f'{"x1":>5} {"x2":>5} {"nu":>2} {"eps^-2":>11} {"exp":>6} {"eps^-1":>12} {"exp":>12} {"r0/C_A":>13} '
          f'{"kappa_1":>12} {"- zeta_2":>10} {"1/eta":>8} {"imag":>8}')
    devs = []
    for rec in rows:
        for nu in (0, 2):
            r = rec[f'nu={nu}']
            ep = rec['expected poles (a_s C_A units)']
            print(f'{rec["x1"]:>5} {rec["x2"]:>5} {nu:2d} {r["eps^-2"]:11.7f} {ep["eps^-2"]:6.2f} {r["eps^-1"]:12.8f} '
                  f'{ep["eps^-1"]:12.8f} {r["eps^0"]:13.9f} {r["kappa_1"]:12.9f} {r["kappa_1 - zeta_2"]:+10.1e} '
                  f'{r["max |1/eta coefficient|"]:8.1e} {r["max |imag|"]:8.1e}')
            devs.append(abs(r['kappa_1 - zeta_2']))
    out = {'symbolic': sym, 'zeta_2': Z2, 'pairs': rows,
           'largest |kappa_1 - zeta_2|': max(devs) if devs else None,
           'largest spread of kappa_1 between pairs and components':
               (max(rec[f'nu={nu}']['kappa_1'] for rec in rows for nu in (0, 2))
                - min(rec[f'nu={nu}']['kappa_1'] for rec in rows for nu in (0, 2))) if rows else None}
    # precision runs at x = 1/2
    for label, suffix in (('quadrature orders 96, 72, 40', '_nq96-72-40'), ('16 points on the eps circle', '_neps16')):
        if os.path.exists(f'n4_laurent_x1o2{suffix}.json') and os.path.exists(f'gg_laurent_x1o2{suffix}.json'):
            base, new = N4.kernel_n4('1/2', '1/2'), N4.kernel_n4('1/2', '1/2', suffix)
            out[f'change of r0/C_A at x = 1/2, {label}'] = {f'nu={nu}': new[nu]['eps^0'] - base[nu]['eps^0'] for nu in (0, 2)}
            print(f'x1 = x2 = 1/2, {label}: change of r0/C_A', {k: f'{v:+.1e}' for k, v in out[f'change of r0/C_A at x = 1/2, {label}'].items()})
    json.dump(out, open('kappa1_n4.json', 'w'), indent=1)
    if devs:
        print(f'largest |kappa_1 - zeta_2| = {max(devs):.1e}; spread {out["largest spread of kappa_1 between pairs and components"]:.1e}')
    print('wrote kappa1_n4.json')


if __name__ == '__main__':
    main()
