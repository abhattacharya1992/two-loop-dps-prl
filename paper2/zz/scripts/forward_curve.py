#!/usr/bin/env python3
"""Approach of the one-loop gg -> ZZ amplitudes with |Delta J_z| = 2 to their forward limit
(paper Sec. 6.2 and the figure of the approach to the forward limit).

Usage: forward_curve.py POINTS_DIR RESULTS_JSON [OUT_JSON]

For +-++ at every point of the scan: p_T^2 and (M^(1) - K)/K, real and imaginary parts, with K
the forward limit in its closed form, K = -4 pi i x(1 - x) in the normalisation of the source
(x = (1 + beta)/2; the fitted limit of results_zz.json agrees with it to 10^-17, which is also
recorded). The closed form is used because at delta = 10^-12 the difference is only about
10^-12 of K, so that a fitted K would limit the ratio (M - K)/(K p_T^2) to about five digits. Like the
gg -> HH amplitude of Sec. 2.5, the difference should vanish like p_T^2 with no logarithm.
"""
import json
import pathlib
import sys
from fractions import Fraction

import mpmath as mp

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from analyse_zz import load, cplx, bnum

mp.mp.dps = 50
H = "+-++"
FROM = 6          # the spread of (M - K)/(K p_T^2) is quoted for delta <= 10^-FROM


def main():
    regular, _ = load(sys.argv[1])
    res = json.loads(pathlib.Path(sys.argv[2]).read_text())
    out = {}
    for beta, recs in sorted(regular.items(), key=lambda kv: float(Fraction(kv[0]))):
        x = (1 + bnum(beta)) / 2
        K = -4j * mp.pi * x * (1 - x)
        K_fit = cplx(res["betas"][beta]["hel"][H]["K"])
        curve = []
        for r in recs:
            d = (cplx(r["M"][H]["1l"]) - K) / K
            curve.append({"delta_exp": r["delta_exp"], "pT2": r["pT2"], "rel_re": mp.nstr(mp.re(d), 20),
                          "rel_im": mp.nstr(mp.im(d), 20)})
        ratios = [mp.mpf(c["rel_re"]) / mp.mpf(c["pT2"]) for c in curve if c["delta_exp"] >= FROM]
        spread = float((max(ratios) - min(ratios)) / abs(sum(ratios) / len(ratios)))
        out[beta] = {"helicity": H, "K": [mp.nstr(mp.re(K), 30), mp.nstr(mp.im(K), 30)],
                     "K_fit_rel_diff": mp.nstr(abs(K_fit / K - 1), 3), "curve": curve,
                     "re_over_pT2_spread": spread, "spread_from_delta_exp": FROM}
        print(f"beta = {beta}: Re (M - K)/K at delta = 10^-{curve[0]['delta_exp']}: {mp.nstr(mp.mpf(curve[0]['rel_re']), 3)}, "
              f"at 10^-{curve[-1]['delta_exp']}: {mp.nstr(mp.mpf(curve[-1]['rel_re']), 3)}; "
              f"Re/p_T^2 constant to {spread:.1e} for delta <= 10^-{FROM}")
    if len(sys.argv) > 3:
        pathlib.Path(sys.argv[3]).write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
