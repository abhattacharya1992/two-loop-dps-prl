#!/usr/bin/env python3
"""Forward points of the one- and two-loop gg -> ZZ helicity amplitudes (massless quark loops,
class A of von Manteuffel and Tancredi, qT scheme at mu^2 = s), from the translated ggvvamp code.

Usage: zz_points.py BUILD OUTDIR BETA K DPS [--recompute]          one point
       zz_points.py BUILD OUTDIR --jobs FILE [--recompute]

Existing output is reused only when its provenance matches. Use a fresh OUTDIR
to preserve an older run, or --recompute to replace its points after evaluation.

Units m_Z = 1. Kinematics as for gg -> HH in Paper 2:
    x = (1 - beta)/(1 + beta),  s = (1 + x)^2/x,  -t = y = x + delta,  delta = 10^-K,
    G = (1/2) ln(1 - y z) with z = -u,  p_T^2 = (y z - 1) x/(1 + x)^2.
The amplitude is M = sum_j A_j T_j(e1, e2, e3*, e4*) with the tensors of eq. (2.14) of
arXiv:1503.08835 (identified from the published projectors, check_tensors.py), helicity
polarisation vectors in the centre-of-mass frame, gluon 1 along +z, Z_3 at polar angle theta.
Each point file holds M for all 36 helicities and for the parts 1l, 2l_CF, 2l_CA.
"""
import argparse
import pathlib
import sys
import time
from fractions import Fraction
from itertools import product

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import mpmath as mp
from ggvv import GGVV, Functions
from zz_cache import cache_provenance, translated_provenance, validate_point

PARTS = ("1l", "2l_CF", "2l_CA")
DIRECT = (1, 2, 4, 5, 8, 9, 16, 17, 18)


def kinematics(beta, k):
    fb = Fraction(beta)
    b = mp.mpf(fb.numerator) / fb.denominator
    x = (1 - b) / (1 + b)
    delta = mp.mpf(10) ** (-k)
    y = x + delta
    s = (1 + x) ** 2 / x
    t = -y
    u = 2 - s - t
    z = -u
    G = mp.log(1 - y * z) / 2
    pT2 = (y * z - 1) * x / (1 + x) ** 2
    return dict(beta=b, x=x, s=s, t=t, u=u, G=G, pT2=pT2)


def all_coefficients(At, Au):
    """The 20 A_j from the nine computed coefficients at (s,t) and at (s,u), equal masses
    (crossing relations of ggvvampprec.cpp; for ma2 = mb2, x34 = x12 and x12x34 = identity)."""
    A = {}
    for j in DIRECT:
        A[j] = At[j]
    A[3], A[6], A[7] = Au[2], Au[5], Au[4]
    A[10], A[11] = -Au[8], -Au[9]
    A[12], A[13] = Au[9], Au[8]
    A[14], A[15] = -At[9], -At[8]
    A[19], A[20] = Au[18], Au[17]
    return A


def pol(theta, phi, lam, E, kmod, m=1):
    ct, st, cp, sp = mp.cos(theta), mp.sin(theta), mp.cos(phi), mp.sin(phi)
    if lam == 0:
        return [kmod / m, E / m * st * cp, E / m * st * sp, E / m * ct]
    r2 = 1 / mp.sqrt(2)
    return [mp.mpc(0), r2 * (-lam * ct * cp + 1j * sp), r2 * (-lam * ct * sp - 1j * cp), r2 * lam * st]


def dot(a, b):
    return a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]


def helicity_amplitudes(A, kin):
    s = kin["s"]
    E = mp.sqrt(s) / 2
    kmod = mp.sqrt(E ** 2 - 1)
    sth = mp.sqrt(kin["pT2"]) / kmod
    cth = mp.sqrt(1 - sth ** 2)
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
        out[key] = {part: sum(A[part][j] * T[j] for j in range(1, 21)) for part in PARTS}
    return out


def run_point(g, beta, k, dps, outdir, *, recompute=False):
    path = pathlib.Path(outdir) / f"zz_b{beta.replace('/', 'o')}_d{k}_p{dps}.json"
    metadata = translated_provenance(g, "zz-forward-amplitude",
        {"beta": beta, "delta_exp": k, "dps": dps, "parts": PARTS,
         "ma2": 1, "mb2": 1, "scheme": "qT", "mu2": "s"})
    if cache_provenance.load_result(path, metadata, recompute=recompute, validator=validate_point) is not None:
        return path, "cached"
    mp.mp.dps = dps
    kin = kinematics(beta, k)
    t0 = time.time()
    Ft = Fu = Functions()          # shared cache: most Li22 arguments do not depend on z
    At = g.compute(kin["s"], kin["t"], 1, 1, Ft)
    Au = g.compute(kin["s"], kin["u"], 1, 1, Fu)
    A = {part: all_coefficients({j: At[j][part] for j in DIRECT}, {j: Au[j][part] for j in DIRECT})
         for part in PARTS}
    H = helicity_amplitudes(A, kin)
    c = lambda z: [mp.nstr(mp.re(z), dps), mp.nstr(mp.im(z), dps)]
    rec = {"beta": beta, "delta_exp": k, "dps": dps, "x": mp.nstr(kin["x"], dps),
           "s": mp.nstr(kin["s"], dps), "t": mp.nstr(kin["t"], dps), "pT2": mp.nstr(kin["pT2"], dps),
           "G": c(kin["G"]), "li22_max_rate": Ft.max_rate,
           "li22_bad_order": Ft.order_bad,
           "M": {h: {p: c(v) for p, v in d.items()} for h, d in H.items()},
           "A": {p: {str(j): c(v) for j, v in A[p].items()} for p in PARTS},
           "sec": round(time.time() - t0, 1), "provenance": metadata}
    cache_provenance.validate_result(rec, validate_point, path)
    cache_provenance.atomic_write(path, rec)
    return path, rec["sec"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build")
    parser.add_argument("outdir")
    parser.add_argument("point", nargs="*", metavar="BETA_K_DPS")
    parser.add_argument("--jobs", type=pathlib.Path)
    parser.add_argument("--recompute", action="store_true",
                        help="evaluate and replace existing output, including legacy records")
    args = parser.parse_args()
    if args.jobs:
        if args.point:
            parser.error("choose a single point or --jobs, not both")
        jobs = [line.split() for line in args.jobs.read_text().splitlines() if line.strip()]
    else:
        jobs = [args.point]
    if any(len(job) != 3 for job in jobs):
        parser.error("each point requires BETA K DPS")
    pathlib.Path(args.outdir).mkdir(parents=True, exist_ok=True)
    g = GGVV(args.build)
    for beta, k, dps in jobs:
        p, sec = run_point(g, beta, int(k), int(dps), args.outdir, recompute=args.recompute)
        print(time.strftime("%H:%M:%S"), p.name, sec, flush=True)


if __name__ == "__main__":
    main()
