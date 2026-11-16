#!/usr/bin/env python3
"""gg -> AA on the forward line: the one- and two-loop Type II finite remainders of
arXiv:1909.08993, evaluated through the P1 package (ggaa.twoloop), colour term by colour term.

The published files are single expressions in x, z, Y = 1 + y, Z = 1 + z, M2, muR2, N and Nf.  They are
evaluated here with N and Nf replaced by the Laurent polynomials of dps.colour, so that one pass gives
the coefficient of every colour monomial (N^2, N n_f, n_f^2, n_f/N at two loops; N, n_f at one loop);
the cost is that of one numerical evaluation.  Units m_A = mu_R = 1, for which the published one-loop
files are the remainders at mu_R = m_A (P1, eq. 2.18) and the two-loop files are 64 times the
remainders in the normalisation of P1 (ggaa.normalization.FILE_FACTOR); the records hold the
remainders in that normalisation, M^(2) = file/64.  Kinematics as for gg -> HH (dps.kinematics):
xk = (1 - beta)/(1 + beta), y = xk + 10^-k.

Precision: the two-loop files cancel strongly, by 57 to 98 digits at 1 - cos(theta) = 10^-10 in P1's
points (results/m1_forward_points.jsonl, 180 against 220 digits), more at small beta; here the N^2
term loses 46 digits at delta = 10^-4 and 58 at 10^-8 for beta = 1/2.  The scan uses 50 + 9k digits
at beta >= 1/2, 15 more below 1/2 and another 15 below 1/4, and a control point with 15 more digits;
`--check` evaluates a point at two precisions and prints the agreeing digits.

    python3 -m dps.backends.aa --gauge 1 --beta 1/2 --k 4                  one point (default digits)
    python3 -m dps.backends.aa --jobs [--workers W]                          the scan, resumable
    python3 -m dps.backends.aa --check --gauge 1 --beta 1/2 --k 4           digit check

The AA evaluator is under common/aa_source of the master folder (legacy: paper1/), or the directory with paper2/ next to
src/ and inputs/), or the directory named by the environment variable P1_ROOT.  Its inputs must be in
place (inputs/fetch_inputs.sh there).
"""
import argparse
import json
import os
import pathlib
import sys
import time
from fractions import Fraction
from multiprocessing import Pool

import mpmath as mp

HERE = pathlib.Path(__file__).resolve().parent
PACKAGE_DIR = HERE.parents[1]                         # paper2/dps
sys.path.insert(0, str(PACKAGE_DIR))
from dps import kinematics                            # noqa: E402
from dps.colour import N, NF, monomial_name           # noqa: E402
from dps.layout import p1_root                        # noqa: E402
from dps.provenance import (atomic_write, fingerprint, load_result, sha256_file,
                            require_finite, require_keys, require_pair)  # noqa: E402

P1_ROOT = p1_root()
ANC = P1_ROOT / 'inputs' / 'arXiv-1909.08993' / 'anc'
OUT = PACKAGE_DIR / 'points_aa'
FILE_FACTOR = 64

BETAS = ('1/2', '1/4', '2/3')
DELTA_EXPONENTS = (4, 5, 6, 7, 8, 9, 10, 11)
BLIND_BETAS = ('1/10', '4/5')


def p1_package():
    if not (P1_ROOT / 'src' / 'python' / 'ggaa').is_dir():
        sys.exit(f'the P1 package was not found at {P1_ROOT} (expected common/aa_source or a legacy Paper I tree); '
                 'set P1_ROOT to its top directory')
    if not ANC.is_dir():
        sys.exit(f'{ANC} not found; download the inputs of the P1 package first (inputs/fetch_inputs.sh)')
    sys.path.insert(0, str(P1_ROOT / 'src' / 'python'))
    from ggaa import twoloop                         # noqa: E402
    return twoloop


def digits(beta, k):
    """Working precision of the scan at delta = 10^-k: the N^2 term of the two-loop files loses about
    35 + 3k digits at beta = 1/2 (check at k = 4 and 8), more at smaller velocities (P1's points:
    34 more at beta = 0.05 than at 0.5), so 50 + 9k digits leave at least 30 for the fits, with 15 more
    below beta = 1/2 and another 15 below 1/4."""
    b = Fraction(beta)
    return 50 + 9 * k + (15 if b < Fraction(1, 2) else 0) + (15 if b < Fraction(1, 4) else 0)


def file_name(gauge, beta, k, dps):
    return f'aa_g{gauge}_b{beta.replace("/", "o")}_d{k}_p{dps}.json'


def symbols(kin):
    """The symbol values of the published files at this point, N and Nf as polynomials."""
    return {'x': kin['xk'], 'z': kin['z'], 'Y': 1 + kin['y'], 'Z': 1 + kin['z'],
            'M2': mp.mpf(1), 'muR2': mp.mpf(1), 'N': N, 'Nf': NF, 'Pi': +mp.pi, 'I': mp.mpc(0, 1)}


def strings(poly, dps, factor=1):
    return {monomial_name(key): [mp.nstr(mp.re(v / factor), dps), mp.nstr(mp.im(v / factor), dps)]
            for key, v in sorted(poly.terms.items())}


_TEXT = {}


def expression(tl, gauge, loops, digest=None):
    path = ANC / tl.FILES[gauge, loops]
    key = (gauge, loops, digest or sha256_file(path))
    if key not in _TEXT:
        _TEXT[key] = tl.read_expression(path)
    return _TEXT[key]


def provenance(tl, gauge, beta, k, dps, regular):
    sources = {f'dps/{name}': PACKAGE_DIR / 'dps' / name
               for name in ('backends/aa.py', 'colour.py', 'kinematics.py', 'provenance.py')}
    sources.update({f'ggaa/{name}': P1_ROOT / 'src/python/ggaa' / name
                    for name in ('twoloop.py', 'li22.py', 'cache.py')})
    metadata = fingerprint('aa-amplitude',
        inputs={tl.FILES[gauge, loops]: ANC / tl.FILES[gauge, loops] for loops in (1, 2)},
        sources=sources, libraries=('mpmath',),
        parameters={'gauge': gauge, 'beta': beta, 'delta_exp': k, 'dps': dps,
                    'regular': regular, 'file_factor': FILE_FACTOR, 'M2': 1, 'muR2': 1})
    expected = {}
    for line in (P1_ROOT / 'inputs/SHA256SUMS').read_text().splitlines():
        fields = line.split()
        if len(fields) == 2:
            expected[pathlib.Path(fields[1].lstrip('*')).name] = fields[0]
    for name, digest in metadata['inputs'].items():
        if name not in expected or digest != expected[name]:
            raise ValueError(f'{ANC / name}: actual SHA-256 does not match inputs/SHA256SUMS')
    return metadata


def validate_point(record):
    parameters = record['provenance']['parameters']
    for key in ('gauge', 'beta', 'delta_exp', 'dps', 'regular'):
        if record[key] != parameters[key]:
            raise ValueError(f'{key} differs from provenance')
    if record['inputs'] != record['provenance']['inputs'] or record['process'] != 'aa':
        raise ValueError('input hashes or process differs from provenance')
    for key in ('xk', 'y', 'costh', 'pT2', 'sec', 'sec_one_loop'):
        require_finite(record[key])
    require_pair(record['G'])
    for key, names in (('A', ('Nf_over_N', 'Nf2', 'N_Nf', 'N2')), ('A1', ('Nf', 'N'))):
        require_keys(record[key], names, exact=True)
        for value in record[key].values():
            require_pair(value)


def evaluate(gauge, beta, k, dps, out_dir=OUT, regular=True, recompute=False):
    """One point: writes out_dir/aa_g{gauge}_b{beta}_d{k}_p{dps}.json and returns (path, seconds)."""
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / file_name(gauge, beta, k, dps)
    tl = p1_package()
    metadata = provenance(tl, gauge, beta, k, dps, regular)
    if load_result(path, metadata, recompute=recompute, validator=validate_point) is not None:
        return str(path), 'cached'
    mp.mp.dps = dps
    kin = kinematics.point(beta, k)
    sym = symbols(kin)
    t0 = time.time()
    one = tl.Evaluator(sym).evaluate(expression(tl, gauge, 1, metadata['inputs'][tl.FILES[gauge, 1]]))
    t1 = time.time()
    ev = tl.Evaluator(sym)
    two = ev.evaluate(expression(tl, gauge, 2, metadata['inputs'][tl.FILES[gauge, 2]]))
    t2 = time.time()
    rec = {'process': 'aa', 'gauge': gauge, 'beta': beta, 'delta_exp': k, 'dps': dps, 'regular': regular,
           'xk': mp.nstr(kin['xk'], dps), 'y': mp.nstr(kin['y'], dps), 'costh': mp.nstr(kin['costh'], dps),
           'pT2': mp.nstr(kin['pT2'], dps),
           'G': [mp.nstr(mp.re(kin['G']), dps), mp.nstr(mp.im(kin['G']), dps)],
           'A': strings(two, dps, FILE_FACTOR), 'A1': strings(one, dps),
           'normalisation': 'M^(2) = file/64, M^(1) = file, m_A = mu_R = 1 (P1 eq. 2.18)',
           'inputs': metadata['inputs'], 'provenance': metadata,
           'function_calls': len(ev._cache), 'sec': round(t2 - t0, 1), 'sec_one_loop': round(t1 - t0, 1)}
    validate_point(rec)
    atomic_write(path, rec)
    return str(path), rec['sec']


def jobs(blind=False):
    """The scan: both components at the three velocities and eight deltas, plus the control point
    at delta = 10^-10 of beta = 1/2 with 15 more digits; with blind=True the two velocities not
    used to find the closed forms of r_0."""
    betas = BLIND_BETAS if blind else BETAS
    js = [(g, b, k, digits(b, k), True) for b in betas for k in DELTA_EXPONENTS for g in (1, 2)]
    if not blind:
        js += [(g, '1/2', 10, digits('1/2', 10) + 15, False) for g in (1, 2)]
    return js


def run_job(job):
    gauge, beta, k, dps, regular, *options = job
    out_dir, recompute = options or (OUT, False)
    return evaluate(gauge, beta, k, dps, out_dir=out_dir, regular=regular, recompute=recompute)


CHECKS = PACKAGE_DIR / 'aa_digit_check.json'


def check(gauge, beta, k, dps, extra=20):
    """Evaluate the two-loop file at dps and dps + extra digits, print the agreeing digits of every
    colour term and record them in aa_digit_check.json."""
    tl = p1_package()
    metadata = provenance(tl, gauge, beta, k, [dps, dps + extra], False)
    vals, secs = {}, {}
    for d in (dps, dps + extra):
        mp.mp.dps = d
        kin = kinematics.point(beta, k)
        t0 = time.time()
        two = tl.Evaluator(symbols(kin)).evaluate(expression(tl, gauge, 2))
        vals[d] = {monomial_name(key): v for key, v in two.terms.items()}
        secs[d] = round(time.time() - t0)
        print(f'dps {d}: {secs[d]} s', flush=True)
    mp.mp.dps = dps + extra + 10
    rec = {'gauge': gauge, 'beta': beta, 'delta_exp': k, 'dps': [dps, dps + extra],
           'provenance': metadata, 'seconds': secs, 'loss': {}}
    for name in sorted(vals[dps]):
        a, b = vals[dps][name], vals[dps + extra][name]
        agree = -mp.log10(abs(a - b) / abs(b)) if a != b else dps
        rec['loss'][name] = round(dps - float(agree), 1)
        print(f'  {name:10s} agreeing digits {float(agree):6.1f}   loss {dps - float(agree):5.1f}   |value| {mp.nstr(abs(b), 6)}')
    all_checks = json.loads(CHECKS.read_text()) if CHECKS.exists() else {}
    all_checks[f'g{gauge}_b{beta.replace("/", "o")}_k{k}'] = rec
    CHECKS.write_text(json.dumps(all_checks, indent=1))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--gauge', type=int, default=1)
    p.add_argument('--beta', default='1/2')
    p.add_argument('--k', type=int, default=4)
    p.add_argument('--dps', type=int)
    p.add_argument('--out', default=str(OUT))
    p.add_argument('--jobs', action='store_true')
    p.add_argument('--blind', action='store_true')
    p.add_argument('--workers', type=int, default=2)
    p.add_argument('--check', action='store_true')
    p.add_argument('--recompute', action='store_true', help='explicitly replace existing point results')
    a = p.parse_args()
    dps = a.dps or digits(a.beta, a.k)
    if a.check:
        check(a.gauge, a.beta, a.k, dps)
    elif a.jobs or a.blind:
        with Pool(a.workers) as pool:
            tasks = [job + (a.out, a.recompute) for job in jobs(blind=a.blind)]
            for path, sec in pool.imap_unordered(run_job, tasks):
                print(time.strftime('%H:%M:%S'), path.split('/')[-1], sec, flush=True)
    else:
        print(evaluate(a.gauge, a.beta, a.k, dps, a.out, recompute=a.recompute))
