#!/usr/bin/env python3
# Forward scan of the two-loop class-B gg -> HH finite remainders of Banerjee,
# Borowka, Dhani, Gehrmann, Ravindran (arXiv:1809.05388), split by colour
# structure and by power of the renormalisation-scale logarithm.
#
# Each ancillary expression has exactly four colour terms,
#
#     Gauge_i = Nf^2 a_i + (Nf/N) b_i + N Nf c_i + N^2 d_i ,
#
# and mu_R appears only through L = Log[M2] - Log[muR2] = ln(m_H^2/mu_R^2),
# at most squared. At every forward point we evaluate the four colour
# functions a_i, b_i, c_i, d_i at L = -1, 0, +1. A quadratic is fixed by three
# values, so this gives each colour function at any mu_R exactly:
#
#     f(L) = f0 + f1 L + f2 L^2,  f1 = (f(+1) - f(-1))/2,  f2 = (f(+1) + f(-1))/2 - f0 .
#
# The polylogarithms are cached per point (evaluate_hh.py), so the twelve
# evaluations cost little more than one. Results are stored as
# A = i Gauge/(4 M2) with M2 = 1, the authors' plotting normalisation (paper Sec. 3.1).
#
# Kinematics: x = (1-beta)/(1+beta), y = x + delta, z = x + 1/x - y;
# G = (1/2) log(1 - y z) on the principal branch (1 - yz < 0, so Im G = pi/2).
#
#   python3 scan_colour_scale.py            # full scan, 2 worker processes
#   python3 scan_colour_scale.py --test     # one cheap point, prints timings
import ast, json, sys, time, argparse, pathlib
from fractions import Fraction
from multiprocessing import Pool
import mpmath as mp
import evaluate_hh as ev

HERE = pathlib.Path(__file__).resolve().parent
for parent in HERE.parents:
    for package in (parent / 'common/dps', parent / 'dps'):
        if (package / 'dps/provenance.py').is_file():
            sys.path.insert(0, str(package))
            break
    else:
        continue
    break
from dps import provenance as cache_provenance

OUT = HERE/'points'
MONOMIALS = ('Nf2', 'Nf_over_N', 'N_Nf', 'N2')      # names used in the output files
L_VALUES = (-1, 0, 1)

BETAS = ('1/2', '1/4', '2/3')
DELTA_EXPONENTS = (4, 5, 6, 7, 8, 9, 10, 11)


def working_digits(k):
    # The rational prefactors carry inverse powers of the Gram determinant
    # (yz - 1 ~ delta), so roughly 3k digits cancel at delta = 10^-k. These are
    # the precisions of the stored scans (45 at 1e-5, 57 at 1e-9) continued;
    # the control points at delta = 1e-10 with 15 more digits test the choice.
    return 30+3*k


def split_terms(n):
    """Compile the four top-level colour terms of Gauge{n}Simp.in separately."""
    sys.setrecursionlimit(100000)
    raw = (ev.DIR/'anc'/f'Gauge{n}Simp.in').read_text(encoding='utf8').split('=', 1)[1]
    raw = raw.strip().rstrip(';').replace('^', '**').replace('[', '(').replace(']', ')')
    body = ast.parse(raw, mode='eval').body
    # same syntax check as evaluate_hh.compile_amp
    allowed = ev.VALID_FUNCTIONS | ev.VALID_PARAMS
    assert all(isinstance(z, ev.VALID_OPS) for z in ast.walk(body))
    assert {z.id for z in ast.walk(body) if isinstance(z, ast.Name)} <= allowed
    terms = []

    def flatten(node):
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
            flatten(node.left); flatten(node.right)
        else:
            terms.append(node)
    flatten(body)
    assert len(terms) == 4, f'expected four colour terms, found {len(terms)}'
    codes = []
    for t in terms:
        expr = ast.Expression(body=t)
        ast.fix_missing_locations(expr)
        codes.append(compile(expr, f'Gauge{n}Simp.in', 'eval'))
    return codes


def env(x, y, muR2, N, Nf):
    return dict(x=x, y=y, M2=mp.mpf(1), muR2=muR2, N=mp.mpf(N), Nf=mp.mpf(Nf), I=mp.j, Pi=mp.pi,
                Log=ev.Log, PolyLog=ev.LiClassical, Li22new=ev.li22, Zeta=ev.Zeta)


def identify(codes):
    """Which colour monomial multiplies each term: ratio of the term at (N,Nf)=(2,3) to (1,1).
    Nf^2 -> 9, Nf/N -> 3/2, N Nf -> 6, N^2 -> 4. Done once, at a generic non-forward point."""
    mp.mp.dps = 20
    for f in (ev.li22, ev.LiClassical, ev.Log, ev.Zeta): f.cache_clear()
    x, y = mp.mpf(1)/3, mp.mpf(1)/2
    expected = {mp.mpf(9): 'Nf2', mp.mpf(3)/2: 'Nf_over_N', mp.mpf(6): 'N_Nf', mp.mpf(4): 'N2'}
    names = []
    for c in codes:
        r = eval(c, {'__builtins__': {}}, env(x, y, mp.mpf('0.5'), 2, 3)) / \
            eval(c, {'__builtins__': {}}, env(x, y, mp.mpf('0.5'), 1, 1))
        match = [v for k, v in expected.items() if abs(r-k) < mp.mpf(10)**-10]
        assert len(match) == 1, r
        names.append(match[0])
    assert sorted(names) == sorted(MONOMIALS)
    return names


_CODES = {}


def validate_point(record):
    parameters = record['provenance']['parameters']
    for key in ('gauge', 'beta', 'delta_exp', 'dps'):
        if record[key] != parameters[key]:
            raise ValueError(f'{key} differs from provenance')
    if record['inputs'] != record['provenance']['inputs']:
        raise ValueError('input hashes differ from provenance')
    for key in ('x', 'y', 'pT2', 'sec'):
        cache_provenance.require_finite(record[key])
    cache_provenance.require_pair(record['G'])
    cache_provenance.require_keys(record['A'], MONOMIALS, exact=True)
    for values in record['A'].values():
        cache_provenance.require_keys(values, [str(v) for v in L_VALUES], exact=True)
        for value in values.values():
            cache_provenance.require_pair(value)


def run_job(job, recompute=False, out_dir=None):
    gauge, beta, k, dps = job
    path = pathlib.Path(out_dir if out_dir is not None else OUT)/f'g{gauge}_b{beta.replace("/", "o")}_d{k}_p{dps}.json'
    input_path = ev.DIR/'anc'/f'Gauge{gauge}Simp.in'
    metadata = cache_provenance.fingerprint('hh-amplitude',
        inputs={input_path.name: input_path},
        sources={'scan_colour_scale.py': __file__, 'evaluate_hh.py': ev.__file__,
                 'dps/provenance.py': cache_provenance.__file__},
        parameters={'gauge': gauge, 'beta': beta, 'delta_exp': k, 'dps': dps,
                    'L_values': L_VALUES, 'M2': 1, 'normalisation': 'i*Gauge/4'},
        libraries=('mpmath',))
    if cache_provenance.load_result(path, metadata, recompute=recompute, validator=validate_point) is not None:
        return str(path), 'cached'
    code_key = (gauge, metadata['inputs'][input_path.name])
    if code_key not in _CODES:
        codes = split_terms(gauge)
        _CODES[code_key] = (codes, identify(codes))
    codes, names = _CODES[code_key]
    mp.mp.dps = dps
    for f in (ev.li22, ev.LiClassical, ev.Log, ev.Zeta): f.cache_clear()
    fb = Fraction(beta)
    b = mp.mpf(fb.numerator)/fb.denominator
    x = (1-b)/(1+b)
    delta = mp.mpf(10)**(-k)
    y = x+delta
    z = x+1/x-y
    G = mp.log(1-y*z)/2
    t0 = time.time()
    vals = {}
    for L in L_VALUES:
        muR2 = mp.e**(-L)          # M2 = 1, so L = -log(muR2)
        for name, c in zip(names, codes):
            raw = eval(c, {'__builtins__': {}}, env(x, y, muR2, 1, 1))
            A = mp.j*raw/4
            vals.setdefault(name, {})[str(L)] = [mp.nstr(mp.re(A), dps), mp.nstr(mp.im(A), dps)]
    rec = {'gauge': gauge, 'beta': beta, 'delta_exp': k, 'dps': dps,
           'x': mp.nstr(x, dps), 'y': mp.nstr(y, dps),
           'pT2': mp.nstr((y*z-1)*x/(1+x)**2, dps),
           'G': [mp.nstr(mp.re(G), dps), mp.nstr(mp.im(G), dps)],
           'A': vals, 'sec': round(time.time()-t0, 1),
           'inputs': metadata['inputs'], 'provenance': metadata}
    validate_point(rec)
    cache_provenance.atomic_write(path, rec)
    return str(path), rec['sec']


def jobs():
    js = [(g, b, k, working_digits(k)) for g in (2, 1) for b in BETAS for k in DELTA_EXPONENTS]
    # precision control: delta = 10^-10 again with 15 more digits (the scan extends to 10^-11)
    js += [(g, '1/2', 10, working_digits(10)+15) for g in (2, 1)]
    return js


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--test', action='store_true')
    p.add_argument('--workers', type=int, default=2)
    p.add_argument('--out', type=pathlib.Path, help='directory for new or resumable results')
    p.add_argument('--recompute', action='store_true', help='explicitly replace existing point results')
    a = p.parse_args()
    if a.out is not None:
        OUT = a.out
    OUT.mkdir(parents=True, exist_ok=True)
    if a.test:
        OUT_TEST = a.out or HERE/'points_test'
        OUT_TEST.mkdir(parents=True, exist_ok=True)
        OUT = OUT_TEST
        for g in (1,):
            print(run_job((g, '1/2', 3, 30), recompute=a.recompute), flush=True)
        sys.exit()
    with Pool(a.workers) as pool:
        from functools import partial
        for path, sec in pool.imap_unordered(partial(run_job, recompute=a.recompute, out_dir=OUT), jobs()):
            print(time.strftime('%H:%M:%S'), path.split('/')[-1], sec, flush=True)
