#!/usr/bin/env python3
# Paper Sec. 3.3, Table 2, Figure 1: colour structure of the two-loop DPS logarithms in gg -> HH.
#
# The class-B remainders are A = n_f^2 a + (n_f/N) b + N n_f c + N^2 d. The DPS
# relations (one-loop coefficients c_1 = -2 pi N F_S and K_H = pi N Kt, with
# c13 = -(4/3) C_A c_1, c12 = -2 beta0 c_1, c22 = -4 C_A K_H) fix how each colour
# term contributes to the leading logarithms at mu_R = m_H:
#
#                 G^3 (spin 0)       G^2 (spin 0)        G^2 (spin 2)     G (spin 2)
#   N^2           (8/3) pi F_S       (44/3) pi F_S       -4 pi Kt         -(32/3) pi Kt
#   N n_f         0                  -(8/3) pi F_S       0                 (8/3) pi Kt
#   n_f/N         0                  0                   0                 0
#   n_f^2         0                  0                   0                 0
#
# F_S = 1/(w(1-w)) + 2w(1-w), Kt = (w^4 + (1-w)^4)/(w(1-w)), w = (1+beta)/2.
# The last column is the RG prediction (4 C_A - 4 beta0) K for the spin-two single
# logarithm (paper Sec. 5), which P1 found in AA at threshold.
#
# So the test is that the leading logarithms are exactly N^2 and N n_f, with n_f
# entering only through beta0 and no C_F n_f (1/N) piece. Each colour term is
# extracted separately, so the result holds for every N_c; the script also prints
# the totals for N_c = 2..5 and n_f = 0, 5.
#
# Fit: leading power plus p_T^2 (1, G, G^2) corrections through consecutive deltas
# (scan_tools.fit_model); "shift" is the change when the window moves one step up in
# delta, used as the uncertainty. The plain leading-power fit is shown too.
#
#   python3 analyse_colour.py            (mpmath; a few seconds)
import json, sys
from fractions import Fraction
import mpmath as mp
from scan_tools import load, fits, fit_model, weights, rel, fmt, COLOURS, COLOUR_WEIGHT

pi = mp.pi
SPEC = {'spin0 G^3': (1, 3, 3), 'spin0 G^2': (1, 3, 2), 'spin0 G^1': (1, 3, 1),
        'spin2 G^2': (2, 2, 2), 'spin2 G^1': (2, 2, 1)}       # gauge, degree, power


def predictions(beta):
    w, FS, Kt = weights(beta)
    z = mp.mpf(0)
    return {
        'spin0 G^3': {'N2': mp.mpf(8)/3*pi*FS, 'N_Nf': z, 'Nf_over_N': z, 'Nf2': z},
        'spin0 G^2': {'N2': mp.mpf(44)/3*pi*FS, 'N_Nf': -mp.mpf(8)/3*pi*FS, 'Nf_over_N': z, 'Nf2': z},
        'spin2 G^2': {'N2': -4*pi*Kt, 'N_Nf': z, 'Nf_over_N': z, 'Nf2': z},
        'spin2 G^1': {'N2': -mp.mpf(32)/3*pi*Kt, 'N_Nf': mp.mpf(8)/3*pi*Kt, 'Nf_over_N': z, 'Nf2': z},
    }


def extract(reg, beta):
    out = {}
    for col in COLOURS:
        out[col] = {}
        for coeff, (gauge, deg, power) in SPEC.items():
            recs = reg[(gauge, beta)]
            w0 = fit_model(recs, col, 0, deg, 2, 0)
            w1 = fit_model(recs, col, 0, deg, 2, 1)
            plain = fits(recs, col, 0, deg)[0]
            out[col][coeff] = (w0[power], w1[power] if w1 else None, plain[power])
    # one power beyond the expected leading-power polynomial, as a diagnostic
    r1, r2 = reg[(1, beta)], reg[(2, beta)]
    e1, e2 = fit_model(r1, 'N2', 0, 4, 2), fit_model(r2, 'N2', 0, 3, 2)
    out['extra'] = {}
    if e1: out['extra']['spin0 G^4 (N^2)'] = e1[4]
    if e2: out['extra']['spin2 G^3 (N^2)'] = e2[3]
    return out


def main(points='points', outjson='results_colour.json'):
    reg, _ = load(points)
    betas = sorted({b for (_, b) in reg}, key=Fraction)
    report = {}
    for beta in betas:
        if (1, beta) not in reg or (2, beta) not in reg:
            continue
        P = predictions(beta)
        X = extract(reg, beta)
        w, FS, Kt = weights(beta)
        print(f'\nbeta = {beta}   w = {fmt(w, 6)}, F_S = {fmt(FS, 10)}, Kt = {fmt(Kt, 10)}, mu_R = m_H, '
              f'{len(reg[(1, beta)])} deltas')
        print(f'  {"coefficient":10s} {"colour":9s} {"extracted (Re)":>22s} {"prediction":>18s} {"rel. dev.":>10s} {"shift":>8s} {"plain fit dev.":>15s}')
        rows = []
        for coeff in ('spin0 G^3', 'spin0 G^2', 'spin2 G^2', 'spin2 G^1'):
            scale = abs(P[coeff]['N2'])
            for col in COLOURS:
                val, nxt, plain = X[col][coeff]
                pred = P[coeff][col]
                dev = abs(val-pred)/(abs(pred) if pred != 0 else scale)
                shift = abs(val-nxt)/scale if nxt is not None else mp.nan
                pdev = abs(plain-pred)/(abs(pred) if pred != 0 else scale)
                tag = '' if pred != 0 else ' (0)'
                print(f'  {coeff:10s} {col:9s} {fmt(mp.re(val), 15):>22s} {fmt(pred, 13)+tag:>18s} '
                      f'{mp.nstr(dev, 2):>10s} {mp.nstr(shift, 2):>8s} {mp.nstr(pdev, 2):>15s}')
                rows.append({'coefficient': coeff, 'colour': col, 're': str(mp.re(val)), 'im': str(mp.im(val)),
                             'prediction': str(pred), 'deviation': str(dev), 'window_shift': str(shift),
                             'plain_fit_deviation': str(pdev)})
        print('  (deviation of a vanishing term is quoted relative to the N^2 coefficient)')
        print('  free extra powers, relative to the leading N^2 coefficient:',
              {k: mp.nstr(abs(v)/abs(P['spin0 G^3' if 'spin0' in k else 'spin2 G^2']['N2']), 2) for k, v in X['extra'].items()})
        print('  spin-0 single log G^1 by colour (no prediction):',
              {c: fmt(mp.re(X[c]['spin0 G^1'][0]), 12) for c in COLOURS})
        print(f'  totals / prediction:   {"N_c":>3s} {"n_f":>3s} {"spin0 G^3":>16s} {"spin0 G^2":>16s} {"spin2 G^2":>16s} {"spin2 G^1":>16s}')
        tot = []
        for N in (2, 3, 4, 5):
            for nf in (0, 5):
                vals = {}
                for coeff in ('spin0 G^3', 'spin0 G^2', 'spin2 G^2', 'spin2 G^1'):
                    v = sum(COLOUR_WEIGHT[c](N, nf)*X[c][coeff][0] for c in COLOURS)
                    p = sum(COLOUR_WEIGHT[c](N, nf)*P[coeff][c] for c in COLOURS)
                    vals[coeff] = mp.re(v)/p
                print(f'                         {N:3d} {nf:3d} ' + ' '.join(f'{mp.nstr(r, 13):>16s}' for r in vals.values()))
                tot.append({'N_c': N, 'n_f': nf, **{k: str(r) for k, r in vals.items()}})
        report[beta] = {'rows': rows, 'totals_ratio_to_prediction': tot,
                        'extra_powers': {k: str(abs(v)) for k, v in X['extra'].items()},
                        'spin0_single_log': {c: [str(mp.re(X[c]['spin0 G^1'][0])), str(mp.im(X[c]['spin0 G^1'][0]))]
                                             for c in COLOURS}}
    json.dump(report, open(outjson, 'w'), indent=1)


if __name__ == '__main__':
    main(*sys.argv[1:])
