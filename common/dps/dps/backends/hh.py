#!/usr/bin/env python3
"""gg -> HH on the forward line: the scan driver of hh/hh_two_loop (scan_colour_scale.py), called
unchanged.  Needs the two ancillary files of arXiv:1809.05388 in hh/hh_two_loop/anc (`make fetch` there).

    python3 -m dps.backends.hh --gauge 1 --beta 1/2 --k 4 [--dps D] [--out DIR]

The record written is that of the stored points (hh/hh_two_loop/points); dps.points.load_hh reads it.
"""
import argparse
import pathlib
import sys

from .. import layout

HH = layout.path('hh_points').parent


def evaluate(gauge, beta, k, dps=None, out_dir=None, recompute=False):
    """One point through scan_colour_scale.run_job; returns (path, seconds) or (path, 'cached')."""
    sys.path.insert(0, str(HH))
    import scan_colour_scale as scan          # noqa: E402
    if not (HH / 'anc' / f'Gauge{gauge}Simp.in').exists():
        sys.exit(f'{HH}/anc/Gauge{gauge}Simp.in not found: run `make fetch` in {HH.parent}')
    return scan.run_job((gauge, beta, k, dps or scan.working_digits(k)),
                        recompute=recompute, out_dir=out_dir)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--gauge', type=int, default=1)
    p.add_argument('--beta', default='1/2')
    p.add_argument('--k', type=int, default=4)
    p.add_argument('--dps', type=int)
    p.add_argument('--out')
    p.add_argument('--recompute', action='store_true', help='explicitly replace an existing point result')
    a = p.parse_args()
    print(evaluate(a.gauge, a.beta, a.k, a.dps, a.out, recompute=a.recompute))
