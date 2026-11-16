#!/usr/bin/env python3
"""Fit stability of the two-loop gg -> ZZ coefficients (paper App. B, figure of the fit windows).

Usage: analyse_windows.py POINTS_DIR [OUT_JSON]

For the C_F and C_A parts of the ++++ amplitude (coefficients of G^3 and G^2, in units of c_1)
and of the +-++ amplitude (G^2 and G, in units of K), the coefficients are refitted through
windows of consecutive deltas ending at every delta of the scan, with the models of
analyse_zz.py: the polynomial in G alone (4 points) and the polynomial plus p_T^2 times a
polynomial of degree 0, 1, 2 or 3 in G (5-8 points; degree 3 is the model of the paper). The
relative deviation from the prediction of Sec. 6 (units of the source, alpha_s/2pi per loop)
is stored for every window and model.
"""
import json
import pathlib
import sys
from fractions import Fraction

import mpmath as mp

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from analyse_zz import load, fit, cplx, bnum

mp.mp.dps = 50
MODELS = {"plain": None, "pT2": 0, "pT2 G": 1, "pT2 G^2": 2, "pT2 G^3": 3}


def main():
    regular, _ = load(sys.argv[1])
    out = {}
    for beta, recs in sorted(regular.items(), key=lambda kv: float(Fraction(kv[0]))):
        x = (1 + bnum(beta)) / 2
        L = mp.log(x * (1 - x))
        # predictions in units of c1 or K, source normalisation (analyse_zz.py)
        pred = {("++++", "2l_CF"): {3: mp.mpf(-4) / 3, 2: -2 * L - 3}, ("++++", "2l_CA"): {3: mp.mpf(2) / 3, 2: 2 * L},
                ("+-++", "2l_CF"): {2: mp.mpf(-4), 1: -4 * L - 2}, ("+-++", "2l_CA"): {2: mp.mpf(2), 1: 4 * L - 2}}
        res = {}
        for (h, part), pp in pred.items():
            one = fit(recs, [cplx(r["M"][h]["1l"]) for r in recs], 1, 2)
            unit = one[1] if h == "++++" else one[0]              # c1 or K from the one-loop fit
            vals = [cplx(r["M"][h][part]) for r in recs]
            rows = []
            for name, pdeg in MODELS.items():
                n = 4 if pdeg is None else 3 + pdeg + 2
                for end in range(n, len(recs) + 1):
                    sub, v = recs[:end], vals[:end]
                    if pdeg is None:
                        G = [cplx(r["G"]) for r in sub[-n:]]
                        M = mp.matrix([[g ** j for j in range(4)] for g in G])
                        c = list(mp.lu_solve(M, mp.matrix(v[-n:])))
                    else:
                        c = fit(sub, v, 3, pdeg)
                    for power, p in pp.items():
                        rows.append({"model": name, "power": power, "points": n,
                                     "delta_exp_min": sub[-1]["delta_exp"], "delta_exp_max": sub[-n]["delta_exp"],
                                     "rel_dev": mp.nstr(abs(c[power] / unit - p) / abs(p), 15)})
            res[f"{h} {part}"] = {"predictions": {str(k): mp.nstr(v, 20) for k, v in pp.items()}, "rows": rows}
            top = max(pp)
            print(f"beta = {beta}  {h} {part}, G^{top}: " + "; ".join(
                f"{m}: {mp.nstr(mp.mpf([r for r in rows if r['model'] == m and r['power'] == top][-1]['rel_dev']), 2)}"
                for m in MODELS))
        out[beta] = res
    if len(sys.argv) > 2:
        pathlib.Path(sys.argv[2]).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
