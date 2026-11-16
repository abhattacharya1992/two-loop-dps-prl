#!/usr/bin/env python3
"""The one-loop gg -> ZZ amplitudes along the lines of fixed beta from 90 degrees to the forward
limit, in units of their DPS terms (paper Fig. 2).

Usage: wide_curve.py POINTS_DIR WIDE_POINTS_DIR [OUT_JSON]

Combines the points of the forward scan (zz_points.py, delta = 10^-4 ... 10^-12) with those of
zz_wide.py (from 90 degrees down to sin^2(theta) = 10^-3.5) and stores, for every point,
    M_++++^(1)/C_1 - G      with C_1 = -8 pi i x^2 (source normalisation): tends to a constant,
    M_+-++^(1)/K            with K = -4 pi i x (1 - x): tends to 1,
with x = (1 + beta)/2, ordered from large to small p_T^2. The ratios do not depend on the
normalisation, so the closed forms are given in that of the source (alpha_s/2pi per loop).
"""
import json
import pathlib
import sys
from fractions import Fraction

import mpmath as mp

mp.mp.dps = 30


def C(v):
    return mp.mpc(mp.mpf(v[0]), mp.mpf(v[1]))


def main():
    files = sorted(pathlib.Path(sys.argv[1]).glob("zz_*.json")) + sorted(pathlib.Path(sys.argv[2]).glob("zzw_*.json"))
    recs = [json.loads(f.read_text()) for f in files]
    for f, r in zip(files, recs):        # every Li_{2,2} series converged (as for the scan, summary.json)
        assert r["li22_bad_order"] == 0 and r["li22_max_rate"] < 1, f.name
    out = {}
    for beta in sorted({r["beta"] for r in recs}, key=Fraction):
        q = Fraction(beta)
        x = (1 + mp.mpf(q.numerator) / q.denominator) / 2
        c1, K = -8j * mp.pi * x ** 2, -4j * mp.pi * x * (1 - x)
        sel = {}
        for r in recs:
            if r["beta"] != beta:
                continue
            key = mp.nstr(mp.mpf(r["pT2"]), 12)          # the precision control repeats a scan point
            if key not in sel or r["dps"] > sel[key]["dps"]:
                sel[key] = r
        rows = []
        for r in sorted(sel.values(), key=lambda r: -mp.mpf(r["pT2"])):
            G = C(r["G"])
            r1 = C(r["M"]["++++"]["1l"]) / c1 - G
            r2 = C(r["M"]["+-++"]["1l"]) / K
            rows.append({"pT2": mp.nstr(mp.mpf(r["pT2"]), 30), "G": [mp.nstr(mp.re(G), 30), mp.nstr(mp.im(G), 30)],
                         "M1_over_C1_minus_G": [mp.nstr(mp.re(r1), 30), mp.nstr(mp.im(r1), 30)],
                         "M2_over_K": [mp.nstr(mp.re(r2), 30), mp.nstr(mp.im(r2), 30)],
                         "source": "wide" if "j" in r else "scan"})
        out[beta] = {"C1": [mp.nstr(mp.re(c1), 30), mp.nstr(mp.im(c1), 30)], "K": [mp.nstr(mp.re(K), 30), mp.nstr(mp.im(K), 30)],
                     "rows": rows}
        print(f"beta = {beta}: {len(rows)} points, p_T^2/m^2 from {rows[0]['pT2'][:8]} to {mp.nstr(mp.mpf(rows[-1]['pT2']), 3)};"
              f" M_+-++/K from {rows[0]['M2_over_K'][0][:8]} to {rows[-1]['M2_over_K'][0][:16]} (real part)")
    if len(sys.argv) > 3:
        pathlib.Path(sys.argv[3]).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
