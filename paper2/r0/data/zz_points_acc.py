#!/usr/bin/env python3
"""gg -> ZZ forward points with the gg -> HH evaluation of Li_{2,2} on its integral branch.

    python3 zz_points_acc.py REPO BUILD OUTDIR BETA K DPS     one point (as zz_points.py BUILD OUTDIR BETA K DPS)
    python3 zz_points_acc.py REPO --selftest [BUILD]          seconds (about a minute with BUILD)

REPO is the clone of two-loop-dps-prl. Nothing in it is changed or written.

paper2/zz/scripts/zz_points.run_point is run unchanged except for one function. ggvv.py sums the
defining series of Li22(a, b) for max(|a|, |ab|) < 0.95 and otherwise integrates
-a Int_0^1 dt ln(t) Li2(abt)/(1 - at) split at the fixed points t = 1/2, 9/10, 99/100. At beta = 19/20
(x = 0.975) paper2/zz/scripts/check_li22_scan.py stops at delta = 1e-4: this integral differs by
2.7e-32 from the accurate value at 50 digits. Here the integral branch is evaluated instead by
paper2/hh/hh_two_loop/evaluate_hh.li22, the routine of every gg -> HH point (integration points
accumulating at t = 1 up to the nearest singular point, 15 more digits). The series branch is unchanged.

For every point, OUTDIR/li22_<point>.json records, over the distinct integral-branch arguments,
  max_abs_diff_fixed_split       the largest |fixed split - evaluate_hh.li22|;
  abs_diff_to_reference_there    for the pair where that is largest, |evaluate_hh.li22 - reference|,
                                 the reference being the same integral with 30 more digits and
                                 points 1 - d, d = 1/2, 1/4, ..., down to a quarter of the distance
                                 of the nearest singular point (and the same at t = 0).
The point file is that of zz_points.py; its provenance lists this script and evaluate_hh.py among the
sources and has the extra parameter li22_integral_branch.
"""
import json
import pathlib
import sys
import tempfile
import time

REPO = pathlib.Path(sys.argv[1]).resolve()
sys.path.insert(0, str(REPO / 'paper2' / 'zz' / 'scripts'))
sys.path.append(str(REPO / 'paper2' / 'hh' / 'hh_two_loop'))

import mpmath as mp                                                   # noqa: E402
import ggvv                                                           # noqa: E402
import zz_points                                                      # noqa: E402
import evaluate_hh as ev                                              # noqa: E402
from zz_cache import cache_provenance, translated_inputs, source_files  # noqa: E402

HERE = pathlib.Path(__file__).resolve()
LABEL = ('paper2/hh/hh_two_loop/evaluate_hh.li22 for max(|a|, |ab|) >= 0.95 '
         '(fixed-split integral of ggvv.py replaced); series of ggvv.py below')
STATS = {'n': 0, 'worst': None, 'pair': None, 'value': None}


def li22(x, y):
    """Replaces ggvv.li22, which ggvv.Functions.li22 calls once per distinct argument pair."""
    x, y = mp.mpmathify(x), mp.mpmathify(y)
    if max(abs(x), abs(x * y)) < mp.mpf('0.95'):
        return ggvv.li22_series(x, y)
    new = ev.li22(x, y)
    old = ggvv.li22_integral(x, y)
    diff = abs(old - new)
    STATS['n'] += 1
    if STATS['worst'] is None or diff > STATS['worst']:
        STATS.update(worst=diff, pair=(x, y), value=new)
    return new


def provenance(g, stage, parameters, *, wide=False):
    """zz_cache.translated_provenance with this script, evaluate_hh.py and the label added."""
    sources = source_files('ggvv.py', 'translate_ggvvamp.py', 'zz_points.py')
    sources['evaluate_hh.py'] = pathlib.Path(ev.__file__).resolve()
    sources['zz_points_acc.py'] = HERE
    metadata = cache_provenance.fingerprint(
        stage, inputs=translated_inputs(g.build), sources=sources,
        parameters=dict(parameters, li22_integral_branch=LABEL), libraries=('mpmath',))
    if metadata['inputs'] != g.input_hashes:
        raise cache_provenance.CacheMismatch('Translated ZZ inputs changed after GGVV was loaded')
    return metadata


ggvv.li22 = li22
zz_points.translated_provenance = provenance


def reference(a, b, guard=30):
    """The integral of evaluate_hh.li22 with `guard` more digits and twice as many points per decade."""
    if mp.im(a) != 0 or mp.im(b) != 0:
        raise ValueError('complex Li22 argument')
    a, b = mp.re(a), mp.re(b)
    with mp.workdps(mp.mp.dps + guard):
        w = a * b
        if a >= 1 or w >= 1:
            raise ValueError('Li22 argument on the integration cut')
        if a == 0 or w == 0:
            return mp.mpf(0)
        near1, near0 = mp.mpf(1), mp.mpf(1)
        for v in (a, w):
            if 0 < v < 1:
                near1 = min(near1, 1 / v - 1)
            elif v < -1:
                near0 = min(near0, -1 / v)
        pts, d = [mp.mpf(1) / 2], mp.mpf(1) / 2
        while d > near0 / 4:
            d /= 2
            pts.insert(0, d)
        pts.insert(0, mp.mpf(0))
        d = mp.mpf(1) / 2
        while d > near1 / 4:
            d /= 2
            pts.append(1 - d)
        pts.append(mp.mpf(1))
        f = lambda t: -mp.log(t) * (mp.polylog(2, w * t) / w) / (1 - a * t) if t else mp.mpf(0)
        value = a * w * mp.quad(f, pts)
    return +value


def write_sidecar(outdir, beta, k, dps, point_name, sec):
    rec = {'beta': beta, 'delta_exp': k, 'dps': dps, 'point_file': point_name, 'sec': sec,
           'li22_integral_branch': LABEL, 'integral_branch_pairs': STATS['n']}
    if STATS['pair'] is not None:
        a, b = STATS['pair']
        rec['max_abs_diff_fixed_split'] = mp.nstr(STATS['worst'], 6)
        rec['at'] = [mp.nstr(a, 25), mp.nstr(b, 25)]
        rec['abs_diff_to_reference_there'] = mp.nstr(abs(STATS['value'] - reference(a, b)), 6)
    side = pathlib.Path(outdir) / f"li22_b{beta.replace('/', 'o')}_d{k}_p{dps}.json"
    side.write_text(json.dumps(rec, indent=1) + '\n')
    return rec


def run(build, outdir, beta, k, dps):
    g = ggvv.GGVV(build)
    path, sec = zz_points.run_point(g, beta, k, dps, outdir)
    if sec == 'cached':
        print(time.strftime('%H:%M:%S'), path.name, 'cached (matching provenance), nothing computed', flush=True)
        return
    rec = write_sidecar(outdir, beta, k, dps, path.name, sec)
    print(time.strftime('%H:%M:%S'), path.name, sec, 's; integral-branch pairs', rec['integral_branch_pairs'],
          '; fixed split off by', rec.get('max_abs_diff_fixed_split'),
          '; evaluate_hh.li22 vs reference', rec.get('abs_diff_to_reference_there'), flush=True)


def selftest(build=None):
    cases = [('series branch', mp.mpf('0.94'), mp.mpf('0.9')),
             ('a = 1 - 1e-6', 1 - mp.mpf(10) ** -6, mp.mpf('0.987654321')),
             ('ab = 1 - 1e-12, b > 1', mp.mpf('0.5'), 2 * (1 - mp.mpf(10) ** -12))]
    for dps in (50, 90):
        mp.mp.dps = dps
        ev.li22.cache_clear()
        for name, a, b in cases:
            hh, ref, fixed = ev.li22(a, b), reference(a, b), ggvv.li22_integral(a, b)
            line = (f'{dps} digits, {name}: |evaluate_hh - reference| = {mp.nstr(abs(hh - ref), 3)}, '
                    f'|fixed split - reference| = {mp.nstr(abs(fixed - ref), 3)}')
            if name == 'series branch':
                line += f', |series - reference| = {mp.nstr(abs(ggvv.li22_series(a, b) - ref), 3)}'
            print(line)
    mp.mp.dps = 50
    ev.li22.cache_clear()
    STATS.update(n=0, worst=None, pair=None, value=None)
    for _, a, b in cases:
        li22(a, b)
    with tempfile.TemporaryDirectory() as tmp:
        rec = write_sidecar(tmp, '19/20', 4, 50, 'test', 0.0)
    print('record:', json.dumps(rec))
    if build:
        g = ggvv.GGVV(build)
        md = provenance(g, 'zz-forward-amplitude', {'beta': '19/20', 'delta_exp': 4, 'dps': 50})
        print('provenance sources:', sorted(md['sources']), '; label present:',
              md['parameters']['li22_integral_branch'] == LABEL)
    print('selftest finished')


if __name__ == '__main__':
    args = sys.argv[2:]
    if args and args[0] == '--selftest':
        selftest(args[1] if len(args) > 1 else None)
    elif len(args) == 5:
        run(args[0], args[1], args[2], int(args[3]), int(args[4]))
    else:
        sys.exit(__doc__)
