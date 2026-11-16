#!/usr/bin/env python3
# Paper App. A: the evaluation of Li_{2,2} in the gg -> HH expressions against an independent
# representation, for every argument with which the published expressions call it on the forward
# lines of the scan and of the blind test of r0. Needs the ancillary files (make fetch in ..).
#
# evaluate_hh.li22 computes Li_{2,2}(a, b) = sum_{m>n>=1} a^m b^n/(m^2 n^2) from
#     -a Int_0^1 dt ln(t) Li_2(a b t)/(1 - a t).
# The reference writes 1/n^2 = Int_0^1 du (-ln u) u^(n-1) and sums over the smaller index first,
#     Li_{2,2}(a, b) = Int_0^1 du (-ln u) [b Li_2(a) - Li_2(a b u)/u]/(1 - b u),
# the representation that ../../zz/scripts/check_li22_values.py uses for gg -> ZZ (at u = 1/b
# the bracket vanishes, so b > 1 is harmless; the quadrature is split there). It shares no code
# with evaluate_hh.li22 apart from mpmath's Li_2 and quadrature, and is evaluated with twice the
# digits.
#
# Each point is evaluated as the scan evaluates it (scan_colour_scale.run_job: the four colour
# terms of Gauge1Simp and Gauge2Simp separately, at L = -1, 0, 1, N = n_f = 1, with 30 + 3k digits
# at delta = 10^-k), with Li22new replaced by a probe that records its arguments; then
# evaluate_hh.li22 at that precision is compared with the reference for every distinct argument.
#
# This check found that the earlier evaluate_hh.li22, which split the integral at t = 1/2 only,
# lost up to 21 digits for arguments with a b just below 1 and b > 1 (a relative error of
# 10^-29 at 51 digits for a = 1/4 + 6 x 10^-9, b = 4 - 9 x 10^-7): there the branch point
# t = 1/(ab) of Li_2 lies just beyond t = 1, and mpmath's quadrature reported convergence before
# resolving it. All points were recomputed with the present li22 on 3 November 2026.
#
#   python3 check_li22_values.py [--workers N]          every point: the scan (beta = 1/2, 1/4, 2/3,
#                                                       delta = 10^-4 ... 10^-11), its precision control
#                                                       and the blind-test velocities beta = 1/10, 4/5
#                                                       (several CPU hours)
#   python3 check_li22_values.py 1/2:4 1/2:10:75        selected points, BETA:K[:DIGITS] for delta = 10^-K
import argparse, sys, time
from fractions import Fraction
from multiprocessing import Pool
import mpmath as mp
import evaluate_hh as ev
from scan_colour_scale import BETAS, DELTA_EXPONENTS, L_VALUES, working_digits, split_terms, env

BLIND_BETAS = ('1/10', '4/5')       # ../../r0/data/run_hh_new_velocities.py
CONTROL = ('1/2', 10, working_digits(10) + 15)     # scan_colour_scale.jobs()


def reference(a, b):
    """Li_{2,2}(a, b) from the sum over the smaller index first, at twice the working digits. The
    integral is split at u = 1/b, where the bracket vanishes, and at points 1 - d, d = 1/2, 1/8, ...,
    accumulating at u = 1 until d is no larger than the distance of the branch point u = 1/(ab) of
    Li_2(a b u) from u = 1 (with fixed points 1/2, 9/10, 99/100 only, the quadrature can report
    convergence before it resolves a branch point just beyond u = 1)."""
    with mp.workdps(2*mp.mp.dps):
        a, b = mp.mpf(a), mp.mpf(b)
        w = a*b
        l2a = mp.polylog(2, a)
        f = lambda u: -mp.log(u)*(b*l2a - mp.polylog(2, w*u)/u)/(1 - b*u) if u else mp.mpf(0)
        s = 1/w - 1 if 0 < w < 1 else mp.mpf(1)
        pts, d = [mp.mpf(0), mp.mpf(1)/2], mp.mpf(1)/2
        while d > s:
            d /= 4
            pts.append(1 - d)
        pts.append(mp.mpf(1))
        if b > 1:
            pts = sorted(set(pts + [1/b]))
        return +mp.quad(f, pts)


def arguments(codes, x, y):
    """the distinct (a, b) with which the compiled colour terms call Li22new at (x, y), as in the scan"""
    seen = []

    def probe(a, b):
        if not any(a == c and b == d for c, d in seen):
            seen.append((a, b))
        return mp.mpf(0)
    for L in L_VALUES:
        muR2 = mp.e**(-L)
        for c in codes:
            e = env(x, y, muR2, 1, 1)
            e['Li22new'] = probe
            eval(c, {'__builtins__': {}}, e)
    return seen


def compare(pairs, cache=None):
    """largest relative difference between evaluate_hh.li22 and the reference over the pairs, and the
    pair at which it occurs; cache, if given, keeps the reference values of a point for both gauges"""
    worst, where = mp.mpf(0), None
    cache = {} if cache is None else cache
    for a, b in pairs:
        if mp.im(a) != 0 or mp.im(b) != 0:
            raise ValueError(f'complex argument {a}, {b}: the reference here is for real arguments')
        a, b = mp.re(a), mp.re(b)
        if (a, b) not in cache:
            cache[(a, b)] = reference(a, b)
        v, ref = ev.li22(a, b), cache[(a, b)]
        dev = abs(v - ref)/max(abs(ref), mp.mpf(10)**(-mp.mp.dps))
        if where is None or dev > worst:
            worst, where = dev, (a, b)
    return worst, where


def points(specs):
    """(beta, K, digits) of the selected points, or of every point of the scan, its precision control
    and the blind test"""
    if specs:
        for s in specs:
            beta, k, *dps = s.split(':')
            yield beta, int(k), int(dps[0]) if dps else working_digits(int(k))
    else:
        for beta in BETAS + BLIND_BETAS:
            for k in DELTA_EXPONENTS:
                yield beta, k, working_digits(k)
        yield CONTROL


CODES = {}       # compiled colour terms, set before the worker processes are forked


def run_point(point):
    beta, k, dps = point
    mp.mp.dps = dps
    for f in (ev.li22, ev.LiClassical, ev.Log, ev.Zeta):
        f.cache_clear()
    fb = Fraction(beta)
    bv = mp.mpf(fb.numerator)/fb.denominator
    x = (1 - bv)/(1 + bv)
    y = x + mp.mpf(10)**(-k)
    lines, overall, cache = [], 0.0, {}
    for n, cs in CODES.items():
        t = time.time()
        pairs = arguments(cs, x, y)
        dev, (a, b) = compare(pairs, cache)
        overall = max(overall, float(dev))
        # distance from t = 1 of the nearest singular point of the integrand of evaluate_hh.li22
        s = min(min([1/v - 1 for v in (mp.re(p), mp.re(p*q)) if 0 < v < 1] or [mp.mpf(1)]) for p, q in pairs)
        lines.append(f'Gauge{n}Simp, beta = {beta}, delta = 10^-{k}, {dps} digits: {len(pairs)} distinct arguments, '
                     f'nearest singular point of the integrand at t = 1 + {mp.nstr(s, 3)}, '
                     f'largest relative difference {mp.nstr(dev, 3)} at a = {mp.nstr(a, 12)}, b = {mp.nstr(b, 12)} '
                     f'({time.time()-t:.0f} s)')
    return '\n'.join(lines), overall


def main():
    parser = argparse.ArgumentParser(description='Li_{2,2} of evaluate_hh.py against the reference')
    parser.add_argument('points', nargs='*', help='BETA:K[:DIGITS]; default: every point')
    parser.add_argument('--workers', type=int, default=1)
    args = parser.parse_args()
    CODES.update({n: split_terms(n) for n in (1, 2)})
    overall = 0.0
    with Pool(args.workers) as pool:
        for text, worst in pool.imap(run_point, list(points(args.points))):
            print(text, flush=True)
            overall = max(overall, worst)
    print(f'largest relative difference over all points: {overall:.3g}')


if __name__ == '__main__':
    main()
