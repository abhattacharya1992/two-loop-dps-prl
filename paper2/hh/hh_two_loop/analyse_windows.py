#!/usr/bin/env python3
# Paper App. B, figure of the fit stability: how the leading two-loop gg -> HH coefficients
# depend on the fit window and on the fit model.
#
# For the N^2 terms of c_13, c_12 (spin 0) and c_22, c_21 (spin 2) at mu_R = m_H, the
# coefficients are refitted through windows of consecutive deltas that end at every delta of
# the scan, with four models: the polynomial in G alone (deg+1 points), and the polynomial
# plus p_T^2 times a polynomial of degree 0, 1 or 2 in G (deg+pdeg+2 points; pdeg = 2 is the
# model of the paper). The relative deviation from the DPS prediction is stored for every
# window and model, so that the figure shows the deviations shrinking as the window moves to
# smaller delta and as power corrections are included.
#
#   python3 analyse_windows.py            (writes results_windows.json)
import json
from fractions import Fraction
import mpmath as mp
from scan_tools import load, fit_values, window_fit, L_parts, G_of
from analyse_colour import predictions, SPEC

mp.mp.dps = 50
MODELS = {'plain': None, 'pT2': 0, 'pT2 G': 1, 'pT2 G^2': 2}      # name -> pdeg


def main(points='points', outjson='results_windows.json'):
    reg, _ = load(points)
    betas = sorted({b for (_, b) in reg}, key=Fraction)
    out = {}
    for beta in betas:
        P = predictions(beta)
        out[beta] = {}
        for coeff, (gauge, deg, power) in SPEC.items():
            if coeff not in P:            # the spin-0 single logarithm has no prediction
                continue
            recs = reg[(gauge, beta)]
            vals = [L_parts(r, 'N2')[0] for r in recs]
            pred = P[coeff]['N2']
            rows = []
            for name, pdeg in MODELS.items():
                n = deg+1 if pdeg is None else deg+pdeg+2
                for end in range(n, len(recs)+1):           # window = recs[end-n:end]
                    if pdeg is None:
                        Gs = [G_of(r) for r in recs[end-n:end]]
                        c = window_fit(Gs, vals[end-n:end], deg)
                    else:
                        c = fit_values(recs[:end], vals[:end], deg, pdeg, 0)
                    rows.append({'model': name, 'points': n, 'delta_exp_min': recs[end-1]['delta_exp'],
                                 'delta_exp_max': recs[end-n]['delta_exp'],
                                 'rel_dev': str(abs(c[power]-pred)/abs(pred))})
            out[beta][coeff] = {'prediction': str(pred), 'rows': rows}
            best = [r for r in rows if r['model'] == 'pT2 G^2'][-1]
            print(f"beta = {beta}  {coeff} (N^2): paper model, smallest window: rel. dev. {mp.nstr(mp.mpf(best['rel_dev']), 2)}; "
                  + '; '.join(f"{m}: {mp.nstr(mp.mpf([r for r in rows if r['model'] == m][-1]['rel_dev']), 2)}" for m in MODELS))
    json.dump(out, open(outjson, 'w'), indent=1)


if __name__ == '__main__':
    main()
