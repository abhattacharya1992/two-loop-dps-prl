#!/usr/bin/env python3
"""Li_{2,2} of ggvv.py at the working precision of the scan, for every argument of the integral branch.

Usage: check_li22_scan.py BUILD_GGV [--workers N] [POINT ...]
       POINT  BETA:K[:DIGITS] (forward line, delta = 10^-K, equal masses m = 1; DIGITS defaults to
              30 + 5K, the precision of the scan); default: every point of the scan (beta = 1/4, 1/2, 2/3,
              K = 4 ... 12), its precision control (beta = 1/2, K = 10, 95 digits) and the blind test of
              ../../r0/data/points_zz (beta = 4/5, K = 4 ... 12)

ggvv.py sums the defining series of Li22(a, b) for max(|a|, |ab|) < 0.95 and otherwise integrates
-a Int_0^1 dt ln(t) Li2(a b t)/(1 - a t) with mpmath, split at t = 1/2, 9/10, 99/100. Towards the
forward limit a b approaches 1 from below for some arguments, with b > 1 and 1 - a b of the order
of delta, and the branch point t = 1/(ab) lies just beyond the end of the path, where a fixed split
can let the tanh-sinh quadrature report convergence before it resolves the branch point (it did
for gg -> HH, ../../hh/hh_two_loop/check_li22_values.py). check_li22_values.py here compares the
evaluation with an independent representation at 30 digits; this script compares it, at the
precision of each scan point, with the same integral split at points 1 - d, d = 1/2, 1/8, ...,
accumulating at t = 1 until d is no larger than the distance of the nearest singular point, and
evaluated with 15 more digits (the evaluation of ../../hh/hh_two_loop/evaluate_hh.li22, which agrees
with two independent representations to the working precision for 1 - ab and 1 - a down to
10^-14 at 42-90 digits). The series branch is not compared.
"""
import argparse
import pathlib
import sys
import time
from multiprocessing import Pool

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import mpmath as mp
from ggvv import GGVV, li22
from zz_points import kinematics
from check_li22_values import Probe
from numeric_checks import require_close

GUARD, RATIO = 15, 4
SCAN = [(b, k, 30 + 5 * k) for b in ("1/4", "1/2", "2/3") for k in range(4, 13)] + [("1/2", 10, 95)]
BLIND = [("4/5", k, 30 + 5 * k) for k in range(4, 13)]
G = None


def points_of(a, w):
    near = mp.mpf(1)
    for v in (a, w):
        if 0 < v < 1:
            near = min(near, 1 / v - 1)
    pts, d = [mp.mpf(0), mp.mpf(1) / 2], mp.mpf(1) / 2
    while d > near:
        d /= RATIO
        pts.append(1 - d)
    pts.append(mp.mpf(1))
    return pts


def accurate(a, b):
    with mp.workdps(mp.mp.dps + GUARD):
        a, b = mp.mpf(a), mp.mpf(b)
        w = a * b
        if a == 0 or w == 0:
            return mp.mpf(0)
        f = lambda t: -mp.log(t) * (mp.polylog(2, w * t) / w) / (1 - a * t) if t else mp.mpf(0)
        v = a * w * mp.quad(f, points_of(a, w))
    return +v


def run_point(point):
    beta, k, dps = point
    mp.mp.dps = dps
    t0 = time.time()
    kin = kinematics(beta, k)
    s, ma2, mb2 = kin["s"], mp.mpf(1), mp.mpf(1)
    Q = [mp.mpf(int(q.split("/")[0])) / int(q.split("/")[1]) for _, q in G.rationals]
    P = Probe()
    for t in (kin["t"], kin["u"]):          # the scan evaluates the code at (s, t) and at (s, u)
        u, x, y, z = G.variables(s, t, ma2, mb2)
        r = [None] * 2000
        r[1], r[2], r[3], r[4], r[5], r[6], r[7] = ma2, t / ma2, u / ma2, mb2 / ma2, x, y, z
        G.p0002(r, Q, P)
    pairs = sorted(set(P.pairs), key=lambda p: max(abs(p[0]), abs(p[0] * p[1])))
    integral = [(a, b) for a, b in pairs if max(abs(a), abs(a * b)) >= mp.mpf("0.95")]
    real = [(mp.re(a), mp.re(b)) for a, b in integral if mp.im(a) == 0 and mp.im(b) == 0]
    if not real or len(real)!=len(integral):
        raise ValueError('Integral-branch Li22 comparison is empty or contains unsupported complex arguments')
    # the relative difference, and the absolute one in units of 10^-digits: the values are at most of
    # order one, and for small values (a b of order 10^-12) mpmath's absolute error target limits the
    # relative accuracy of either evaluation, which is not the effect looked for here
    worst, where, worst_abs, near = mp.mpf(0), None, mp.mpf(0), mp.mpf(1)
    for a, b in real:
        v, ref = li22(a, b), accurate(a, b)
        # The paper allows the last three digits of this absolute comparison
        # to differ. Enforce that budget instead of only printing it.
        require_close(v,ref,rtol=0,atol=1000*mp.mpf(10)**(-dps),label='Li22 scan')
        dev = abs(v - ref) / max(abs(ref), mp.mpf(10) ** -dps)
        worst_abs = max(worst_abs, abs(v - ref) * mp.mpf(10) ** dps)
        near = min([near] + [1 / q - 1 for q in (a, a * b) if 0 < q < 1])
        if where is None or dev > worst:
            worst, where = dev, (a, b)
    text = (f"beta = {beta}, delta = 10^-{k}, {dps} digits: {len(pairs)} distinct arguments, {len(real)} real ones "
            f"on the integral branch" + (f" ({len(integral) - len(real)} complex, not compared)" if len(real) < len(integral) else "")
            + f", nearest singular point at t = 1 + {mp.nstr(near, 3)}; largest relative difference {mp.nstr(worst, 3)}"
            + (f" at a = {mp.nstr(where[0], 12)}, b = {mp.nstr(where[1], 12)}" if where else "")
            + f"; largest absolute difference {mp.nstr(worst_abs, 3)} x 10^-{dps} ({time.time() - t0:.0f} s)")
    return text, float(worst), float(worst_abs)


def main():
    global G
    ap = argparse.ArgumentParser()
    ap.add_argument("build")
    ap.add_argument("points", nargs="*")
    ap.add_argument("--workers", type=int, default=1)
    a = ap.parse_args()
    G = GGVV(a.build)
    pts = []
    for spec in a.points:
        beta, k, *dps = spec.split(":")
        pts.append((beta, int(k), int(dps[0]) if dps else 30 + 5 * int(k)))
    pts = pts or SCAN + BLIND
    worst, ratio = 0.0, 0.0
    with Pool(a.workers) as pool:
        for text, w, rr in pool.imap(run_point, pts):
            print(text, flush=True)
            worst, ratio = max(worst, w), max(ratio, rr)
    print(f"largest relative difference over all points: {worst:.3g}; largest absolute difference in units of "
          f"10^-digits: {ratio:.3g}")


if __name__ == "__main__":
    main()
