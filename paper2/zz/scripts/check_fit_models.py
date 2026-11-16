#!/usr/bin/env python3
"""The choice of the fit model without the predictions, and least-squares fits (paper App. B).

Usage: check_fit_models.py POINTS_DIR [OUT_JSON]

For the C_F and C_A parts of the eight amplitudes of the two groups of Sec. 6 (first group
++++, ++--, ----, --++, coefficients of G^3 and G^2 in units of c_1; second group +-++, +---,
-+--, -+++, coefficients of G^2 and G in units of K), at the three velocities, every fit model
of analyse_windows.py (the polynomial of degree 3 in G alone, and the polynomial plus p_T^2 times
a polynomial of degree P = 0..3 in G) is fitted through the window ending at the smallest delta
(window 0) and through the window moved one step up (window 1). Two criteria choose a model:
    window shift  |c(window 0) - c(window 1)|, which does not use the predictions;
    deviation     |c(window 0) - prediction|.
The script reports the model each criterion selects, the largest ratio of the deviation to the
window shift for the model of the paper (P = 3 in both groups), and least-squares fits of that
model through one more point than unknowns, in units of the window shift.

Before anything else it reproduces the window-0 coefficients of results_zz.json (the fits of the
paper) and stops if they differ.
"""
import json
import pathlib
import sys
from fractions import Fraction

import mpmath as mp

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from analyse_zz import load, fit, cplx, bnum

mp.mp.dps = 50
MODELS = {"plain": None, "P=0": 0, "P=1": 1, "P=2": 2, "P=3": 3}
PAPER = "P=3"
GROUPS = {1: ("++++", "++--", "----", "--++"), 2: ("+-++", "+---", "-+--", "-+++")}


def model_fit(recs, vals, pdeg, window):
    if pdeg is None:
        end = len(recs) - window
        if end - 4 < 0:
            return None
        G = [cplx(r["G"]) for r in recs[end - 4:end]]
        M = mp.matrix([[g ** j for j in range(4)] for g in G])
        return list(mp.lu_solve(M, mp.matrix(vals[end - 4:end])))
    return fit(recs, vals, 3, pdeg, window)


def lsq_fit(recs, vals, pdeg, npts):
    n = 3 + pdeg + 2
    if npts > len(recs):
        return None
    with mp.workdps(2 * mp.mp.dps):
        A = mp.matrix(npts, n)
        b = mp.matrix(npts, 1)
        for i, r in enumerate(recs[-npts:]):
            G, p = cplx(r["G"]), mp.mpf(r["pT2"])
            for j in range(4):
                A[i, j] = G ** j
            for k in range(pdeg + 1):
                A[i, 4 + k] = p * G ** k
            b[i] = vals[-npts:][i]
        AH = A.H
        sol = mp.lu_solve(AH * A, AH * b)
        return [+sol[j] for j in range(4)]


def main():
    regular, _ = load(sys.argv[1])
    stored = json.load(open(pathlib.Path(__file__).parent.parent / "results_zz.json"))["betas"]
    out, sel_shift, sel_dev = {}, {}, {}
    worst_ratio, lsq_ratios, lsq_devs, paper_devs = 0, [], [], []
    for beta, recs in sorted(regular.items(), key=lambda kv: float(Fraction(kv[0]))):
        x = (1 + bnum(beta)) / 2
        L = mp.log(x * (1 - x))
        pred = {1: {"2l_CF": {3: mp.mpf(-4) / 3, 2: -2 * L - 3}, "2l_CA": {3: mp.mpf(2) / 3, 2: 2 * L}},
                2: {"2l_CF": {2: mp.mpf(-4), 1: -4 * L - 2}, "2l_CA": {2: mp.mpf(2), 1: 4 * L - 2}}}
        out[beta] = {}
        for grp, hels in GROUPS.items():
            for h in hels:
                one = fit(recs, [cplx(r["M"][h]["1l"]) for r in recs], 1, 2)
                unit = one[1] if grp == 1 else one[0]
                for part, pp in pred[grp].items():
                    vals = [cplx(r["M"][h][part]) for r in recs]
                    ref = stored[beta]["hel"][h][part]["f"]
                    for power, p in pp.items():
                        rows = {}
                        for name, pdeg in MODELS.items():
                            c0, c1 = model_fit(recs, vals, pdeg, 0), model_fit(recs, vals, pdeg, 1)
                            if c0 is None or c1 is None:
                                continue
                            rows[name] = {"shift": abs((c0[power] - c1[power]) / unit),
                                          "deviation": abs(c0[power] / unit - p), "c0": c0[power], "c1": c1[power]}
                        mine = rows[PAPER]["c0"]
                        assert abs(mine - cplx(ref[power])) <= mp.mpf(10) ** -25 * max(1, abs(mine)), (beta, h, part, power)
                        by_shift = min(rows, key=lambda m: rows[m]["shift"])
                        by_dev = min(rows, key=lambda m: rows[m]["deviation"])
                        sel_shift[(grp, by_shift)] = sel_shift.get((grp, by_shift), 0) + 1
                        sel_dev[(grp, by_dev)] = sel_dev.get((grp, by_dev), 0) + 1
                        ratio = rows[PAPER]["deviation"] / rows[PAPER]["shift"]
                        worst_ratio = max(worst_ratio, ratio)
                        paper_devs.append(ratio)
                        cl = lsq_fit(recs, vals, MODELS[PAPER], 3 + MODELS[PAPER] + 3)
                        lsq = None
                        if cl is not None:
                            shift_abs = abs(rows[PAPER]["c0"] - rows[PAPER]["c1"])
                            lsq = abs(cl[power] - rows[PAPER]["c0"]) / shift_abs
                            lsq_ratios.append(lsq)
                            lsq_devs.append(abs(cl[power] / unit - p) / rows[PAPER]["shift"])
                        out[beta][f"{h} {part} G^{power}"] = {
                            "group": grp,
                            "models": {m: {"window_shift": mp.nstr(v["shift"], 6), "deviation": mp.nstr(v["deviation"], 6)}
                                       for m, v in rows.items()},
                            "selected_by_window_shift": by_shift, "selected_by_deviation": by_dev,
                            "paper_model_deviation_over_shift": mp.nstr(ratio, 4),
                            "lsq_one_more_point_change_in_window_shifts": mp.nstr(lsq, 4) if lsq is not None else None}
        print(f"beta = {beta}: done", flush=True)
    summary = {"selected_by_window_shift": {f"group {g}: {m}": n for (g, m), n in sorted(sel_shift.items())},
               "selected_by_deviation": {f"group {g}: {m}": n for (g, m), n in sorted(sel_dev.items())},
               "paper_model_largest_deviation_over_shift": mp.nstr(worst_ratio, 3),
               "lsq_one_more_point_change_in_window_shifts": [mp.nstr(min(lsq_ratios), 3), mp.nstr(max(lsq_ratios), 3)] if lsq_ratios else None,
               "lsq_one_more_point_deviation_in_window_shifts": [mp.nstr(min(lsq_devs), 3), mp.nstr(max(lsq_devs), 3)] if lsq_devs else None,
               "paper_model_deviation_in_window_shifts": [mp.nstr(min(paper_devs), 3), mp.nstr(max(paper_devs), 3)]}
    for k, v in summary.items():
        print(f"{k}: {v}")
    if len(sys.argv) > 2:
        pathlib.Path(sys.argv[2]).write_text(json.dumps({"summary": summary, "fits": out}, indent=1))


if __name__ == "__main__":
    main()
