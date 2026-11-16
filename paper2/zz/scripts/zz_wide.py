#!/usr/bin/env python3
"""gg -> ZZ helicity amplitudes away from the forward limit, from 90 degrees down to
sin^2(theta) = 10^-3.5, next to the first point of the scan of zz_points.py, delta = 10^-4 (paper
Fig. 2).

Usage: zz_wide.py BUILD OUTDIR BETA J [J ...] [--recompute]

Existing output is reused only when its provenance matches. Use a fresh OUTDIR
to preserve an older run, or --recompute to replace its points after evaluation.

The point J = 0, 1, 2, ... has sin^2(theta) = 10^(-J/2), so that p_T^2 = p_max^2 10^(-J/2) with
p_max^2/m^2 = beta^2/(1 - beta^2) at 90 degrees; J = 0 is the 90-degree point. The kinematics
are those of zz_points.py with delta = -t - x = (s/2) beta (1 - cos theta), and the working
precision is 30 + 5k digits for delta = 10^-k (k rounded up, at least 0), as in the scan.
Each point file holds M for all 36 helicities and the parts 1l, 2l_CF, 2l_CA.
"""
import argparse
import pathlib
import sys
import time
from fractions import Fraction

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import mpmath as mp
from ggvv import GGVV, Functions
from zz_cache import cache_provenance, translated_provenance, validate_point
from zz_points import PARTS, DIRECT, all_coefficients, pol, dot
from itertools import product


def kinematics(beta, j):
    fb = Fraction(beta)
    b = mp.mpf(fb.numerator) / fb.denominator
    x = (1 - b) / (1 + b)
    s = (1 + x) ** 2 / x
    sin2 = mp.mpf(10) ** (-mp.mpf(j) / 2)
    delta = s / 2 * b * (1 - mp.sqrt(1 - sin2))
    y = x + delta
    t = -y
    u = 2 - s - t
    z = -u
    G = mp.log(1 - y * z) / 2
    pT2 = (y * z - 1) * x / (1 + x) ** 2
    return dict(beta=b, x=x, s=s, t=t, u=u, G=G, pT2=pT2, delta=delta, sin2=sin2)


def helicity_amplitudes(A, kin):
    """As in zz_points.py, but with cos(theta) = 1 - delta/(s beta/2) taken from the kinematics:
    zz_points.py computes it as sqrt(1 - sin^2), which loses half the digits near 90 degrees."""
    s = kin["s"]
    E = mp.sqrt(s) / 2
    kmod = mp.sqrt(E ** 2 - 1)
    cth = 1 - kin["delta"] / (s * kin["beta"] / 2)
    sth = mp.sqrt(kin["pT2"]) / kmod
    assert abs(cth ** 2 + sth ** 2 - 1) < mp.mpf(10) ** (10 - mp.mp.dps)
    theta = mp.atan2(sth, cth)
    p1 = [E, 0, 0, E]
    p2 = [E, 0, 0, -E]
    p3 = [E, kmod * sth, 0, kmod * cth]
    p4 = [E, -kmod * sth, 0, -kmod * cth]
    tchk = dot([p1[i] - p3[i] for i in range(4)], [p1[i] - p3[i] for i in range(4)])
    assert abs(tchk - kin["t"]) < mp.mpf(10) ** (10 - mp.mp.dps), "t mismatch"
    out = {}
    for l1, l2, l3, l4 in product((1, -1), (1, -1), (1, 0, -1), (1, 0, -1)):
        e1 = pol(0, 0, l1, E, E)
        e2 = pol(mp.pi, 0, l2, E, E)
        e3 = [mp.conj(c) for c in pol(theta, 0, l3, E, kmod)]
        e4 = [mp.conj(c) for c in pol(mp.pi - theta, mp.pi, l4, E, kmod)]
        e12, e13, e14, e23, e24, e34 = dot(e1, e2), dot(e1, e3), dot(e1, e4), dot(e2, e3), dot(e2, e4), dot(e3, e4)
        e1p3, e2p3 = dot(e1, p3), dot(e2, p3)
        e3p1, e3p2, e4p1, e4p2 = dot(e3, p1), dot(e3, p2), dot(e4, p1), dot(e4, p2)
        T = {1: e12 * e34, 2: e13 * e24, 3: e14 * e23,
             4: e12 * e3p1 * e4p1, 5: e12 * e3p1 * e4p2, 6: e12 * e3p2 * e4p1, 7: e12 * e3p2 * e4p2,
             8: e13 * e2p3 * e4p1, 9: e13 * e2p3 * e4p2, 10: e14 * e2p3 * e3p1, 11: e14 * e2p3 * e3p2,
             12: e23 * e1p3 * e4p1, 13: e23 * e1p3 * e4p2, 14: e24 * e1p3 * e3p1, 15: e24 * e1p3 * e3p2,
             16: e34 * e1p3 * e2p3,
             17: e1p3 * e2p3 * e3p1 * e4p1, 18: e1p3 * e2p3 * e3p1 * e4p2,
             19: e1p3 * e2p3 * e3p2 * e4p1, 20: e1p3 * e2p3 * e3p2 * e4p2}
        key = "".join({1: "+", -1: "-", 0: "0"}[l] for l in (l1, l2, l3, l4))
        out[key] = {part: sum(A[part][i] * T[i] for i in range(1, 21)) for part in PARTS}
    return out


def run_point(g, beta, j, outdir, *, recompute=False):
    mp.mp.dps = 30
    k = max(0, int(mp.ceil(-mp.log10(kinematics(beta, j)["delta"]))))
    dps = 30 + 5 * k
    path = pathlib.Path(outdir) / f"zzw_b{beta.replace('/', 'o')}_j{j}_p{dps}.json"
    metadata = translated_provenance(g, "zz-wide-amplitude",
        {"beta": beta, "j": j, "dps": dps, "parts": PARTS,
         "precision_rule": "30+5*max(0,ceil(-log10(delta))); initial dps=30",
         "ma2": 1, "mb2": 1, "scheme": "qT", "mu2": "s"}, wide=True)
    validator = lambda record: validate_point(record, wide=True)
    if cache_provenance.load_result(path, metadata, recompute=recompute, validator=validator) is not None:
        return path, "cached"
    mp.mp.dps = dps
    kin = kinematics(beta, j)
    assert abs(kin["pT2"] / (kin["beta"] ** 2 / (1 - kin["beta"] ** 2) * kin["sin2"]) - 1) < mp.mpf(10) ** (5 - dps)
    t0 = time.time()
    Ft = Fu = Functions()          # shared cache, as in zz_points.py
    At = g.compute(kin["s"], kin["t"], 1, 1, Ft)
    Au = g.compute(kin["s"], kin["u"], 1, 1, Fu)
    A = {part: all_coefficients({i: At[i][part] for i in DIRECT}, {i: Au[i][part] for i in DIRECT})
         for part in PARTS}
    H = helicity_amplitudes(A, kin)
    c = lambda v: [mp.nstr(mp.re(v), dps), mp.nstr(mp.im(v), dps)]
    rec = {"beta": beta, "j": j, "sin2theta": mp.nstr(kin["sin2"], dps), "delta": mp.nstr(kin["delta"], dps),
           "dps": dps, "x": mp.nstr(kin["x"], dps), "s": mp.nstr(kin["s"], dps), "t": mp.nstr(kin["t"], dps),
           "pT2": mp.nstr(kin["pT2"], dps), "G": c(kin["G"]),
           "li22_max_rate": Ft.max_rate, "li22_bad_order": Ft.order_bad,
           "M": {h: {p: c(v) for p, v in d.items()} for h, d in H.items()},
           "sec": round(time.time() - t0, 1), "provenance": metadata}
    cache_provenance.validate_result(rec, validator, path)
    cache_provenance.atomic_write(path, rec)
    return path, rec["sec"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build")
    parser.add_argument("outdir")
    parser.add_argument("beta")
    parser.add_argument("j", nargs="+", type=int)
    parser.add_argument("--recompute", action="store_true",
                        help="evaluate and replace existing output, including legacy records")
    args = parser.parse_args()
    pathlib.Path(args.outdir).mkdir(parents=True, exist_ok=True)
    g = GGVV(args.build)
    for j in args.j:
        p, sec = run_point(g, args.beta, j, args.outdir, recompute=args.recompute)
        print(time.strftime("%H:%M:%S"), p.name, sec, flush=True)
