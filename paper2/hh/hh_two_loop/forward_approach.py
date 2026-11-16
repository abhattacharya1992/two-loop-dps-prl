#!/usr/bin/env python3
# Paper Sec. 2.5: numerical side of the completeness argument for the finite
# opposite-helicity one-loop term, using gg -> HH.
#
# The one-loop class-B amplitude is not in the ancillary files of Banerjee et al.,
# but it is contained in the mu_R dependence of the two-loop remainder
# (analyse_one_loop.py): with L = ln(m_H^2/mu_R^2) the N^2 colour term of the
# L^1 coefficient of the spin-two remainder is -11 m_A2 + (tree), where m_A2 is the
# one-loop spin-two amplitude per unit N and the tree term is O(p_T^2).
#
# The argument of Sec. 2.5 says that m_A2 tends to the DPS value pi*Ktilde with no
# logarithm, and that everything else is analytic in the transverse momentum, so the
# difference vanishes like p_T^2. This script prints, point by point,
#   r = (m_A2 - pi Kt)/p_T^2      which must tend to a constant, and
#   dr/dp_T^2                     between neighbouring deltas, which must also be
#                                 constant if r = r0 + r1 p_T^2 + ..., with no log.
# A p_T^2 ln p_T^2 term in m_A2 would make dr/dp_T^2 grow like 1/p_T^2, and a
# p_T^4 ln p_T^2 term would make it drift like ln p_T^2.
#
# The point-by-point values (m_A2 - pi Kt)/(pi Kt) are stored as well, for the figure of the
# approach to the forward limit.
#
#   python3 forward_approach.py
import json
import mpmath as mp
from scan_tools import load, L_parts, weights, fit_model, fmt

mp.mp.dps = 40
reg, _ = load('points')
print(f'{"beta":5s} {"delta":>6s} {"p_T^2/m^2":>9s} {"m_A2 (one loop, per unit N)":>28s} '
      f'{"r = (m_A2 - pi Kt)/p_T^2":>25s} {"dr/dp_T^2":>13s} {"Im m_A2/p_T^2":>14s}')
out = {}
for beta in ('1/4', '1/2', '2/3'):
    if (2, beta) not in reg:
        continue
    w, FS, Kt = weights(beta)
    rs, slopes, ims, curve = [], [], [], []
    prev = None
    for rec in reg[(2, beta)]:
        mA2 = -L_parts(rec, 'N2')[1]/11
        pT2 = mp.mpf(rec['pT2'])
        # (m_A2 - pi Kt)/(pi Kt) point by point, for the figure of the approach to the forward limit
        curve.append({'delta_exp': rec['delta_exp'], 'pT2': str(pT2),
                      'rel_re': str(mp.re(mA2-mp.pi*Kt)/(mp.pi*Kt)), 'rel_im': str(mp.im(mA2)/(mp.pi*Kt))})
        r = mp.re(mA2-mp.pi*Kt)/pT2
        im = mp.im(mA2)/pT2
        slope = None if prev is None else (r-prev[0])/(pT2-prev[1])
        print(f'{beta:5s} {"1e-%d" % rec["delta_exp"]:>6s} {mp.nstr(pT2, 3):>9s} {fmt(mp.re(mA2), 16):>28s} '
              f'{mp.nstr(r, 13):>25s} {"" if slope is None else mp.nstr(slope, 7):>13s} {mp.nstr(im, 6):>14s}')
        if rec['delta_exp'] >= 6:          # p_T^2 <= 7e-7 here: the p_T^4 term in r is below 1e-5 of r
            rs.append(r)
            ims.append(im)
            if slope is not None and rec['delta_exp'] >= 7:
                slopes.append(slope)
        prev = (r, pT2)
    best = fit_model(reg[(2, beta)], 'N2', 1, 2, 2)
    print(f'      free fit of the L^1 part in powers of G: |G^2| = {mp.nstr(abs(best[2]), 2)}, '
          f'|G| = {mp.nstr(abs(best[1]), 2)} (absolute; the constant is 11 pi Kt = {fmt(11*mp.pi*Kt, 6)})')
    # constant + p_T^2 through the two smallest deltas
    r1, r2 = reg[(2, beta)][-2:]
    m1, m2 = -L_parts(r1, 'N2')[1]/11, -L_parts(r2, 'N2')[1]/11
    p1, p2 = mp.mpf(r1['pT2']), mp.mpf(r2['pT2'])
    K = (m2*p1-m1*p2)/(p1-p2)
    print(f'      forward limit from c + d p_T^2 at the two smallest deltas: {fmt(mp.re(K), 15)}'
          f'   pi Kt = {fmt(mp.pi*Kt, 15)}   rel {mp.nstr(abs(K-mp.pi*Kt)/(mp.pi*Kt), 2)}')
    spread = lambda v: (max(v)-min(v))/abs(sum(v)/len(v))
    print(f'      delta <= 1e-6: r constant to {mp.nstr(spread(rs), 2)}, Im m_A2/p_T^2 constant to '
          f'{mp.nstr(spread(ims), 2)}; delta <= 1e-7: dr/dp_T^2 constant to {mp.nstr(spread(slopes), 2)}')
    out[beta] = {'forward_limit': str(mp.re(K)), 'pi_Kt': str(mp.pi*Kt),
                 'rel_dev': str(abs(K-mp.pi*Kt)/(mp.pi*Kt)),
                 'r_values': [str(x) for x in rs], 'r_spread': str(spread(rs)),
                 'constant_digits': int(mp.floor(-mp.log10(spread(rs)))),
                 'dr_dpT2_values': [str(x) for x in slopes], 'dr_dpT2_spread': str(spread(slopes)),
                 'Im_over_pT2_spread': str(spread(ims)),
                 'free_fit_abs_G2_G1': [str(abs(best[2])), str(abs(best[1]))],
                 'curve': curve}
json.dump(out, open('results_forward.json', 'w'), indent=1)
