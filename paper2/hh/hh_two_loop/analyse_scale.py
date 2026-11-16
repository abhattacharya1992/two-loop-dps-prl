#!/usr/bin/env python3
# Paper Sec. 3.5: renormalisation-scale dependence of the two-loop DPS logarithms in gg -> HH.
#
# The remainder is A + B L + C L^2 with L = ln(m_H^2/mu_R^2) (see analyse_one_loop.py).
# The leading logarithms are mu_R independent if and only if
#   spin 0: B has no G^3 or G^2 term and C has no G^2 or G term,
#   spin 2: B has no G^2 term and C has no G^2 or G term,
# for every colour term separately. RG invariance predicts more: B is linear in G
# (spin 0) or constant (spin 2), and C is 3 beta0^2 times the tree.
#
# Fits include p_T^2 (1, G, G^2) corrections (scan_tools.fit_model).
#
# The script also refits the total amplitude (N_c = 3, n_f = 5) at three scales,
# mu_R^2 = m_H^2/2 (the scale of the stored scans and of the published plots),
# m_H^2 and 4 m_H^2, and compares the leading coefficients; and it rebuilds the
# amplitude values in stored_forward_results.jsonl (an earlier, direct evaluation of the
# full expressions at N_c = 3, mu_R^2 = m_H^2/2) from the colour/scale decomposition.
# 28 of the 38 stored points lie on the grid of this scan; they check the pipeline
# point by point.
#
#   python3 analyse_scale.py              (mpmath; a few seconds)
import json, sys
from fractions import Fraction
import mpmath as mp
from scan_tools import load, fit_model, fit_values, L_parts, weights, rel, fmt, COLOURS, COLOUR_WEIGHT

pi = mp.pi


def total_at(rec, L, N=3, nf=5):
    return sum(COLOUR_WEIGHT[c](N, nf)*(lambda f: f[0]+f[1]*L+f[2]*L**2)(L_parts(rec, c)) for c in COLOURS)


def leading_at(records, L, gauge, N=3, nf=5):
    # leading-power polynomial plus p_T^2 (1, G, G^2), smallest-delta window
    vals = [total_at(r, L, N, nf) for r in records]
    return fit_values(records, vals, 3 if gauge == 1 else 2, 2)


def main(points='points', stored='stored_forward_results.jsonl', outjson='results_scale.json'):
    reg, control = load(points)
    betas = sorted({b for (_, b) in reg}, key=Fraction)
    report = {'absent_terms': {}, 'three_scales': {}, 'stored_points': {}, 'precision_control': []}

    print('Terms that must vanish if the leading logs are mu_R independent (|coefficient|,')
    print('divided by the size of the corresponding leading coefficient of the N^2 term at L^0):')
    for beta in betas:
        for gauge in (1, 2):
            if (gauge, beta) not in reg:
                continue
            recs = reg[(gauge, beta)]
            lead = fit_model(recs, 'N2', 0, 3 if gauge == 1 else 2, 2)
            scale = abs(lead[3]) if gauge == 1 else abs(lead[2])
            row = {}
            for col in COLOURS:
                Bc = fit_model(recs, col, 1, 3 if gauge == 1 else 2, 2)
                Cc = fit_model(recs, col, 2, 2, 2)
                if gauge == 1:
                    row[col] = {'B G^3': abs(Bc[3])/scale, 'B G^2': abs(Bc[2])/scale,
                                'C G^2': abs(Cc[2])/scale, 'C G^1': abs(Cc[1])/scale}
                else:
                    row[col] = {'B G^2': abs(Bc[2])/scale, 'B G^1': abs(Bc[1])/scale,
                                'C G^2': abs(Cc[2])/scale, 'C G^1': abs(Cc[1])/scale, 'C G^0': abs(Cc[0])/scale}
            worst = max(v for d in row.values() for v in d.values())
            print(f'  beta = {beta:4s} spin {0 if gauge == 1 else 2}: largest relative size {mp.nstr(worst, 2)}')
            report['absent_terms'][f'{beta} gauge{gauge}'] = {c: {k: str(v) for k, v in d.items()} for c, d in row.items()}

    print('\nLeading coefficients of the full remainder (N_c = 3, n_f = 5) at three scales:')
    print(f'  {"beta":5s} {"coefficient":12s} {"mu^2 = m^2/2":>20s} {"mu^2 = m^2":>20s} {"mu^2 = 4 m^2":>20s} {"max spread":>11s}')
    for beta in betas:
        if (1, beta) not in reg or (2, beta) not in reg:
            continue
        Ls = (mp.log(2), mp.mpf(0), -mp.log(4))
        f1 = [leading_at(reg[(1, beta)], L, 1) for L in Ls]
        f2 = [leading_at(reg[(2, beta)], L, 2) for L in Ls]
        entries = {'spin0 G^3': [f[3] for f in f1], 'spin0 G^2': [f[2] for f in f1],
                   'spin2 G^2': [f[2] for f in f2], 'spin2 G^1': [f[1] for f in f2],
                   'spin0 G^1': [f[1] for f in f1]}
        report['three_scales'][beta] = {}
        for k, v in entries.items():
            spread = max(abs(a-v[1]) for a in v)/abs(v[1])
            note = '' if k != 'spin0 G^1' else '   <- moves with mu_R, see next line'
            print(f'  {beta:5s} {k:12s} ' + ' '.join(f'{fmt(mp.re(a), 12):>20s}' for a in v) + f' {mp.nstr(spread, 2):>11s}{note}')
            report['three_scales'][beta][k] = {'values_re': [str(mp.re(a)) for a in v], 'spread': str(spread)}
        # the spin-0 single log must move by -3 beta0 c_1 per unit of L, with c_1 = -2 pi N F_S (one loop)
        w, FS, Kt = weights(beta)
        c1 = -2*pi*3*FS
        expected = -(11*3-2*5)*c1                       # 3 beta0 = 11 N - 2 n_f
        slope = (entries['spin0 G^1'][0]-entries['spin0 G^1'][2])/(Ls[0]-Ls[2])
        print(f'  {beta:5s} G^1 of M_1 moves by {fmt(mp.re(slope), 14)} per unit L; -3 beta0 c_1 = {fmt(expected, 14)}'
              f'   rel {mp.nstr(rel(slope, expected), 2)}')
        report['three_scales'][beta]['spin0 G^1 slope in L'] = {'value': str(mp.re(slope)), 'expected': str(expected),
                                                               'rel_dev': str(rel(slope, expected))}

    # pointwise reconstruction of the stored N_c = 3 scans at mu^2 = m^2/2
    stored_rows = [json.loads(l) for l in open(stored)]
    worst, count = mp.mpf(0), 0
    for s in stored_rows:
        b = str(Fraction(s['beta']).limit_denominator(100))
        k = round(-mp.log10(mp.mpf(s['delta'])))
        match = [r for r in reg.get((s['gauge'], b), []) if r['delta_exp'] == k]
        if not match or abs(mp.mpf(s['delta'])-mp.mpf(10)**(-k)) > mp.mpf(10)**(-k-20):
            continue
        nf = s.get('nf', 5) if s.get('nf') is not None else 5
        rebuilt = total_at(match[0], mp.log(2), 3, nf)
        ref = mp.mpc(mp.mpf(s['physical_re']), mp.mpf(s['physical_im']))
        d = abs(rebuilt-ref)/abs(ref)
        worst = max(worst, d)
        count += 1
    print(f'\nStored N_c = 3 scans rebuilt from the colour/scale split at mu^2 = m^2/2: {count} points,'
          f' largest relative difference {mp.nstr(worst, 3)}')
    report['stored_points'] = {'compared': count, 'max_rel_diff': str(worst)}

    # extra-digit control points
    for c in control:
        twin = [r for r in reg[(c['gauge'], c['beta'])] if r['delta_exp'] == c['delta_exp']][0]
        d = max(abs(total_at(c, L)-total_at(twin, L))/abs(total_at(c, L)) for L in (-1, 0, 1))
        print(f'Precision control, gauge {c["gauge"]}, beta {c["beta"]}, delta 1e-{c["delta_exp"]}: '
              f'{twin["dps"]} vs {c["dps"]} digits agree to {mp.nstr(d, 3)}')
        report['precision_control'].append({'gauge': c['gauge'], 'beta': c['beta'], 'delta_exp': c['delta_exp'],
                                            'dps': [twin['dps'], c['dps']], 'rel_diff': str(d)})
    json.dump(report, open(outjson, 'w'), indent=1)


if __name__ == '__main__':
    main(*sys.argv[1:])
