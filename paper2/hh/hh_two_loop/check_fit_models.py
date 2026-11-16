#!/usr/bin/env python3
# Paper App. B: the choice of the fit model without the predictions, and least-squares fits.
#
# For the N^2 terms of the four leading two-loop gg -> HH coefficients (c_13, c_12 of M_1;
# c_22, c_21 of M_2) at the three velocities, each fit model of analyse_windows.py (the
# polynomial in G alone, and the polynomial plus p_T^2 times a polynomial of degree P = 0..3 in
# G) is fitted through the window ending at the smallest delta (window 0) and through the window
# moved one step up (window 1). Two criteria choose a model:
#   window shift  |c(window 0) - c(window 1)|, which does not use the predictions;
#   deviation     |c(window 0) - prediction|, the criterion of the text.
# The script reports, for every coefficient, the model each criterion selects. For the model of
# the paper (P = 2) it also fits by least squares through one and two more points than unknowns
# (the added points lie further from the forward limit) and reports the change of the
# coefficient in units of the window shift.
#
# Before anything else it reproduces the window-0 coefficients and window shifts of
# results_colour.json (the fits of the paper) and stops if they differ.
#
#   python3 check_fit_models.py            (writes results_fit_models.json; seconds)
import json
from fractions import Fraction
import mpmath as mp
from scan_tools import load, fit_values, window_fit, L_parts, G_of
from analyse_colour import predictions, SPEC

mp.mp.dps = 50
MODELS = {'plain': None, 'P=0': 0, 'P=1': 1, 'P=2': 2, 'P=3': 3}
PAPER = 'P=2'
COEFFS = ('spin0 G^3', 'spin0 G^2', 'spin2 G^2', 'spin2 G^1')


def model_fit(recs, vals, deg, pdeg, window):
    """window fit of one model; None if the scan has too few points for that window"""
    if pdeg is None:
        n = deg+1
        end = len(recs)-window
        if end-n < 0:
            return None
        return window_fit([G_of(r) for r in recs[end-n:end]], vals[end-n:end], deg)
    return fit_values(recs, vals, deg, pdeg, window)


def lsq_fit(recs, vals, deg, pdeg, npts):
    """least-squares fit of sum_j c_j G^j + p_T^2 sum_k d_k G^k through the npts smallest deltas"""
    n = deg+pdeg+2
    if npts > len(recs):
        return None
    with mp.workdps(2*mp.mp.dps):
        A = mp.matrix(npts, n)
        b = mp.matrix(npts, 1)
        for i, r in enumerate(recs[-npts:]):
            G, p = G_of(r), mp.mpf(r['pT2'])
            for j in range(deg+1):
                A[i, j] = G**j
            for k in range(pdeg+1):
                A[i, deg+1+k] = p*G**k
            b[i] = vals[-npts:][i]
        AH = A.H
        sol = mp.lu_solve(AH*A, AH*b)       # normal equations at twice the working precision
        return [+sol[j] for j in range(deg+1)]


def main(points='points', outjson='results_fit_models.json'):
    reg, _ = load(points)
    stored = json.load(open('results_colour.json'))
    betas = sorted({b for (_, b) in reg}, key=Fraction)
    out, chosen_shift, chosen_dev, lsq_ratios = {}, {}, {}, {1: [], 2: []}
    lsq_devs, paper_devs = {1: [], 2: []}, []
    for beta in betas:
        P = predictions(beta)
        out[beta] = {}
        for coeff in COEFFS:
            gauge, deg, power = SPEC[coeff]
            recs = reg[(gauge, beta)]
            vals = [L_parts(r, 'N2')[0] for r in recs]
            pred = P[coeff]['N2']
            rows = {}
            for name, pdeg in MODELS.items():
                c0 = model_fit(recs, vals, deg, pdeg, 0)
                c1 = model_fit(recs, vals, deg, pdeg, 1)
                if c0 is None or c1 is None:
                    continue
                rows[name] = {'shift': abs(c0[power]-c1[power])/abs(pred), 'deviation': abs(c0[power]-pred)/abs(pred),
                              'c0': c0[power], 'c1': c1[power]}
            # the paper's fit, as stored by analyse_colour.py
            ref = [r for r in stored[beta]['rows'] if r['coefficient'] == coeff and r['colour'] == 'N2'][0]
            mine = rows[PAPER]
            assert abs(mp.re(mine['c0'])-mp.mpf(ref['re'])) < mp.mpf(10)**-30*abs(pred), (beta, coeff)
            assert abs(mine['shift']-mp.mpf(ref['window_shift'])) < mp.mpf(10)**-25, (beta, coeff)
            by_shift = min(rows, key=lambda m: rows[m]['shift'])
            by_dev = min(rows, key=lambda m: rows[m]['deviation'])
            chosen_shift[by_shift] = chosen_shift.get(by_shift, 0)+1
            chosen_dev[by_dev] = chosen_dev.get(by_dev, 0)+1
            n = deg+MODELS[PAPER]+2
            shift_abs = abs(mine['c0']-mine['c1'])
            paper_devs.append(abs(mine['c0']-pred)/shift_abs)
            lsq = {}
            for extra in (1, 2):
                cl = lsq_fit(recs, vals, deg, MODELS[PAPER], n+extra)
                if cl is not None:
                    lsq[extra] = abs(cl[power]-mine['c0'])/shift_abs
                    lsq_ratios[extra].append(lsq[extra])
                    lsq_devs[extra].append(abs(cl[power]-pred)/shift_abs)
            out[beta][coeff] = {
                'models': {m: {'window_shift': mp.nstr(v['shift'], 6), 'deviation': mp.nstr(v['deviation'], 6)}
                           for m, v in rows.items()},
                'selected_by_window_shift': by_shift, 'selected_by_deviation': by_dev,
                'lsq_change_in_window_shifts': {f'{n+e} points': mp.nstr(v, 4) for e, v in lsq.items()}}
            print(f"beta = {beta:3s}  {coeff:9s} (N^2): smallest window shift {by_shift}, smallest deviation {by_dev}; "
                  + '; '.join(f"{m}: shift {mp.nstr(v['shift'], 2)}, dev. {mp.nstr(v['deviation'], 2)}" for m, v in rows.items())
                  + '; least squares (P=2): ' + ', '.join(f"{n+e} points {mp.nstr(v, 2)} shifts" for e, v in lsq.items()))
    rng = lambda v: [mp.nstr(min(v), 3), mp.nstr(max(v), 3)] if v else None
    summary = {'selected_by_window_shift': chosen_shift, 'selected_by_deviation': chosen_dev,
               'paper_model_deviation_in_window_shifts': rng(paper_devs),
               'lsq_one_more_point_deviation_in_window_shifts': rng(lsq_devs[1]),
               'lsq_two_more_points_deviation_in_window_shifts': rng(lsq_devs[2]),
               'lsq_one_more_point_in_window_shifts': [mp.nstr(min(lsq_ratios[1]), 3), mp.nstr(max(lsq_ratios[1]), 3)],
               'lsq_two_more_points_in_window_shifts': ([mp.nstr(min(lsq_ratios[2]), 3), mp.nstr(max(lsq_ratios[2]), 3)]
                                                        if lsq_ratios[2] else None)}
    print('models selected by the window shift:', chosen_shift, ' by the deviation:', chosen_dev)
    print('paper model (P=2), deviation from the prediction in units of the window shift:', summary['paper_model_deviation_in_window_shifts'])
    print('least squares with the paper model, change in units of the window shift: one more point',
          summary['lsq_one_more_point_in_window_shifts'], '; two more points', summary['lsq_two_more_points_in_window_shifts'])
    print('least squares with the paper model, deviation from the prediction in units of the window shift: one more point',
          summary['lsq_one_more_point_deviation_in_window_shifts'], '; two more points', summary['lsq_two_more_points_deviation_in_window_shifts'])
    json.dump({'summary': summary, 'fits': out}, open(outjson, 'w'), indent=1)


if __name__ == '__main__':
    main()
