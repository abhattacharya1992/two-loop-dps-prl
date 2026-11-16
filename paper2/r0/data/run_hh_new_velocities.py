#!/usr/bin/env python3
"""Blind test of the r0 closed form: gg -> HH points at two velocities not used to find it.

Runs the task-B scan (scan_colour_scale.run_job, unchanged) at beta = 1/10 and 4/5, writing to
points_hh/ here so that the points of task B are not touched. Run isolated:
    python3 -I run_hh_new_velocities.py [--workers 2] [--recompute]
Existing points are kept if their provenance (input files, code, parameters, runtime) matches;
--recompute replaces them, as needed after a change of the evaluation code.
"""
import sys, pathlib, time, argparse
from multiprocessing import Pool

HERE = pathlib.Path(__file__).resolve().parent
BDIR = HERE.parents[1] / 'hh' / 'hh_two_loop'
sys.path.insert(0, str(BDIR))
import scan_colour_scale as scan          # noqa: E402

BETAS = ('1/10', '4/5')
OUT = HERE / 'points_hh'


def job(j, recompute=False):
    scan.OUT = OUT
    return scan.run_job(j, recompute=recompute, out_dir=OUT)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--workers', type=int, default=2)
    p.add_argument('--recompute', action='store_true', help='replace existing points')
    a = p.parse_args()
    OUT.mkdir(exist_ok=True)
    js = [(g, b, k, scan.working_digits(k)) for b in BETAS for g in (2, 1) for k in scan.DELTA_EXPONENTS]
    with Pool(a.workers) as pool:
        from functools import partial
        for path, sec in pool.imap_unordered(partial(job, recompute=a.recompute), js):
            print(time.strftime('%H:%M:%S'), path.split('/')[-1], sec, flush=True)
