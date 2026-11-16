#!/usr/bin/env python3
"""Check of the Li_{2,2} evaluation of ggvv.py against an independent representation, for the
arguments with which the code calls it.

Usage: check_li22_values.py BUILD_GGV POINT...
       POINT  'bench' (the benchmark point of evalcmdlinegg.cpp) or BETA:K (forward line,
              delta = 10^-K, equal masses m = 1)

ggvv.py evaluates Li22(a, b) = sum_{i>j>=1} a^i b^j/(i^2 j^2) from this series or, near the edge
of convergence, from -a Int_0^1 dt ln(t) Li2(a b t)/(1 - a t). The reference here writes
1/j^2 = Int_0^1 du (-ln u) u^(j-1) and sums over j first:

    Li22(a, b) = Int_0^1 du (-ln u) [b Li2(a) - Li2(a b u)/u] / (1 - b u),

which shares no code with ggvv.py apart from mpmath's Li2 and quadrature. For every call of the
generated code (p0002) at the given points, ggvv.py's value at 30 digits is compared with the
reference computed with 60 digits (for arguments such as a ~ 10^-12, b ~ 10^12 the quadrature of
the reference needs more than 30 digits).
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import mpmath as mp
from ggvv import GGVV, li22
from zz_points import kinematics
from numeric_checks import require_close

mp.mp.dps = 30


class Probe:
    """Stands in for the special functions of p0002 and records the Li22 arguments."""

    def __init__(self):
        self.pairs = []

    def li22(self, a, b):
        self.pairs.append((a, b))
        return mp.mpf(0)

    def lin(self, n, x):
        return mp.mpf(0)

    def log(self, x):
        return mp.mpf(0)

    def zeta(self, n):
        return mp.mpf(0)


def reference(a, b):
    with mp.workdps(2 * mp.mp.dps):
        return +_reference(mp.mpf(a), mp.mpf(b))


def _reference(a, b):
    l2a = mp.polylog(2, a)
    f = lambda u: -mp.log(u) * (b * l2a - mp.polylog(2, a * b * u) / u) / (1 - b * u)
    pts = [0, mp.mpf(1) / 2, mp.mpf(9) / 10, mp.mpf(99) / 100, 1]
    if 0 < 1 / b < 1:
        pts = sorted(set(pts + [1 / b]))
    return mp.quad(f, pts)


def main():
    if len(sys.argv)<3:raise SystemExit('Provide a build and at least one comparison point')
    g = GGVV(sys.argv[1])
    worst = 0
    for spec in sys.argv[2:]:
        if spec == "bench":
            label, args = "benchmark point", ("15625", "-2384.89727896833", "6467.215561", "625")
            s, t, ma2, mb2 = map(mp.mpf, args)
        else:
            beta, k = spec.split(":")
            kin = kinematics(beta, int(k))
            label, (s, t, ma2, mb2) = f"beta = {beta}, delta = 10^-{k}", (kin["s"], kin["t"], mp.mpf(1), mp.mpf(1))
        u, x, y, z = g.variables(s, t, ma2, mb2)
        Q = [mp.mpf(int(q.split("/")[0])) / int(q.split("/")[1]) for _, q in g.rationals]
        r = [None] * 2000
        r[1], r[2], r[3], r[4], r[5], r[6], r[7] = ma2, t / ma2, u / ma2, mb2 / ma2, x, y, z
        P = Probe()
        g.p0002(r, Q, P)
        pairs = sorted(set(P.pairs), key=lambda p: max(abs(p[0]), abs(p[0] * p[1])))
        if not pairs:raise ValueError('No Li22 arguments were compared')
        dev = 0
        for a, b in pairs:
            v, ref = li22(a, b), reference(a, b)
            scale=max(abs(ref),mp.mpf('1e-20'))
            # At 30 working digits, require the independently tested 24-digit
            # accuracy (including the small-value floor used in the report).
            require_close(v/scale,ref/scale,rtol=0,atol='1e-24',label='Li22 reference')
            dev = max(dev, float(abs(v - ref) / scale))
        rate = max(max(abs(a), abs(a * b)) for a, b in pairs)
        print(f"{label}: {len(pairs)} distinct arguments, max(|a|, |ab|) up to {mp.nstr(rate, 15)}; "
              f"largest relative difference from the reference: {mp.nstr(dev, 3)}", flush=True)
        worst = max(worst, dev)
    print(f"largest relative difference over all points: {mp.nstr(worst, 3)}")


if __name__ == "__main__":
    main()
