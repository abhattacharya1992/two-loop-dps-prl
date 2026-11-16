#!/usr/bin/env python3
"""gg -> ZZ on the forward line: the point driver of zz/scripts (zz_points.py), called unchanged.
Needs the translated VVamp code in zz/build/ggv (`make inputs` and `make build` in zz/).

    python3 -m dps.backends.zz --beta 1/2 --k 4 [--dps D] [--out DIR]

The record written is that of the stored points (zz/points); dps.points.load_zz reads it.  The
working precision of the stored scan is 30 + 5k.
"""
import argparse
import pathlib
import sys

from .. import layout

ZZ = layout.path('zz_points').parent


def evaluate(beta, k, dps=None, out_dir=None, build=None, recompute=False):
    """One point through zz_points.run_point; returns (path, seconds) or (path, 'cached')."""
    sys.path.insert(0, str(ZZ / 'scripts'))
    import zz_points                           # noqa: E402
    from ggvv import GGVV                      # noqa: E402
    build = pathlib.Path(build) if build else ZZ / 'build' / 'ggv'
    if not build.is_dir():
        sys.exit(f'{build} not found: run `make inputs` and `make build` in {ZZ}')
    out_dir = pathlib.Path(out_dir) if out_dir else ZZ / 'points'
    out_dir.mkdir(exist_ok=True)
    return zz_points.run_point(GGVV(str(build)), beta, k, dps or 30 + 5 * k, out_dir,
                               recompute=recompute)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--beta', default='1/2')
    p.add_argument('--k', type=int, default=4)
    p.add_argument('--dps', type=int)
    p.add_argument('--out')
    p.add_argument('--build')
    p.add_argument('--recompute', action='store_true', help='explicitly replace an existing point result')
    a = p.parse_args()
    print(evaluate(a.beta, a.k, a.dps, a.out, a.build, recompute=a.recompute))
