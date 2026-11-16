#!/usr/bin/env python3
# Tests for the Paper 2 package, run from the top directory:
#
#   python3 tests/test_release.py           checks that need no third-party input (~1 min)
#   python3 tests/test_release.py --full    also the published expressions themselves: checksums,
#                                           one direct evaluation, the colour split (2-3 min more;
#                                           run fetch_inputs.py first), and, if the P1 package and
#                                           its inputs are found, the gg -> AA curves of Fig. 2
#                                           (a minute more; the P1 package is ../paper1 of the master
#                                           folder, the one this folder lives in, or the one in the
#                                           environment variable P1)
#
# The functions are named test_* so that pytest finds them too.
#
# What is tested:
#   - the stored scan points are complete;
#   - rerunning every analysis on them reproduces the stored result files;
#   - the stored results satisfy the statements made in the paper, with margins;
#   - the tables and the numbers quoted in the paper (../manuscript) are regenerated exactly
#     from them and from the gg -> ZZ results in ../zz;
#   - the one-loop, subtraction, RG and orientation scripts run and pass their own checks.
import hashlib, json, os, pathlib, shutil, subprocess, sys, tempfile
from decimal import Decimal, InvalidOperation, localcontext

ROOT = pathlib.Path(__file__).resolve().parents[1]
HH = ROOT/'hh_two_loop'
PY = sys.executable
CHECKSUMS = {'Gauge1Simp.in': '3b2f6edede722228e8526a73b1aed1feaca6ec5bb834c5632c1995a1a444b7a3',
             'Gauge2Simp.in': '6fc90c7cc0bcb2fd8025aae8ca5279a4249519873a94168ff7373e0bc5ee7a44'}


class Skipped(Exception):
    pass


def skip(msg):
    try:
        import pytest
        pytest.skip(msg)
    except ImportError:
        raise Skipped(msg)


def run(script, cwd, *args):
    r = subprocess.run([PY, script, *args], cwd=cwd, capture_output=True, text=True)
    assert r.returncode == 0, f'{script} failed:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}'
    return r.stdout


def same(a, b, rel=1e-20, path=''):
    """Compare two JSON trees; strings that parse as numbers are compared numerically."""
    if isinstance(a, dict):
        assert isinstance(b, dict) and set(a) == set(b), f'{path}: keys differ'
        for k in a:
            same(a[k], b[k], rel, f'{path}/{k}')
    elif isinstance(a, list):
        assert isinstance(b, list) and len(a) == len(b), f'{path}: lengths differ'
        for i, (x, y) in enumerate(zip(a, b)):
            same(x, y, rel, f'{path}[{i}]')
    else:
        try:
            x, y = Decimal(str(a)), Decimal(str(b))
        except (InvalidOperation, TypeError, ValueError):
            assert a == b, f'{path}: {a!r} != {b!r}'
            return
        if x.is_nan() and y.is_nan():
            return
        assert x.is_finite() and y.is_finite(), f'{path}: non-finite comparison {a} vs {b}'
        tolerance = Decimal(str(rel))
        with localcontext() as ctx:
            ctx.prec = max(50, len(x.as_tuple().digits), len(y.as_tuple().digits),
                           len(tolerance.as_tuple().digits)) + 10
            assert abs(x-y) <= tolerance*max(abs(x), abs(y), Decimal('1e-300')), f'{path}: {a} vs {b}'


def test_comparison_precision():
    # Both pairs round to the same binary64 value, but differ by more than 1e-20 relatively.
    for a, b in (('1.000000000000000001', '1'), ('1.000000000000000001e-100', '1e-100')):
        try:
            same(a, b)
        except AssertionError:
            pass
        else:
            raise AssertionError(f'high-precision discrepancy was ignored: {a} vs {b}')
    same('1.000000000000000000001', '1')
    same({'value': ['1e-100', 'unchanged']}, {'value': ['1e-100', 'unchanged']})


def test_scan_reports_failures():
    if not shutil.which('make'):
        skip('make is needed to check the scan recipe')
    with tempfile.TemporaryDirectory() as tmp:
        (pathlib.Path(tmp)/'hh_two_loop'/'logs').mkdir(parents=True)
        result = subprocess.run(['make', '--no-print-directory', '-f', str(ROOT/'Makefile'),
                                 'scan', 'PY=false'], cwd=tmp, capture_output=True, text=True)
        assert result.returncode != 0, 'make scan hid a failing evaluator behind tee'


def test_points_complete():
    pts = sorted((HH/'points').glob('*.json'))
    assert len(pts) == 50, len(pts)
    for p in pts:
        r = json.loads(p.read_text())
        assert set(r['A']) == {'N2', 'N_Nf', 'Nf_over_N', 'Nf2'}, p.name
        for v in r['A'].values():
            assert set(v) == {'-1', '0', '1'}, p.name


def test_analyses_reproduce_stored_results():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp)
        for f in ('scan_tools.py', 'analyse_colour.py', 'analyse_one_loop.py', 'analyse_scale.py',
                  'forward_approach.py', 'analyse_windows.py', 'stored_forward_results.jsonl'):
            shutil.copy(HH/f, tmp/f)
        (tmp/'points').symlink_to(HH/'points')
        for script in ('analyse_colour.py', 'analyse_one_loop.py', 'analyse_scale.py', 'forward_approach.py',
                       'analyse_windows.py'):
            run(script, tmp)
        for name in ('results_colour.json', 'results_one_loop.json', 'results_scale.json', 'results_forward.json',
                     'results_windows.json'):
            same(json.loads((tmp/name).read_text()), json.loads((HH/name).read_text()), path=name)


def test_statements_of_the_paper():
    colour = json.loads((HH/'results_colour.json').read_text())
    one = json.loads((HH/'results_one_loop.json').read_text())
    scale = json.loads((HH/'results_scale.json').read_text())
    fwd = json.loads((HH/'results_forward.json').read_text())
    for b, v in colour.items():
        for r in v['rows']:
            limit = 1e-10 if float(r['prediction']) != 0 else 1e-11     # Sec. 3.3, Table 2
            assert float(r['deviation']) < limit, (b, r['coefficient'], r['colour'])
        for t in v['totals_ratio_to_prediction']:
            assert all(abs(float(t[k])-1) < 1e-10 for k in t if k not in ('N_c', 'n_f'))
    for b, v in one.items():                                           # Sec. 3.4, Table 3
        assert float(v['one_loop']['spin0_slope_rel_dev']) < 1e-15
        assert float(v['one_loop']['spin2_forward_rel_dev']) < 1e-17
        assert float(v['one_loop']['spin0_slope_mF_abs']) < 1e-14
        assert float(v['tree_max_rel_diff']) < 1e-14
        assert all(float(x['rel_dev']) < 1e-10 for x in v['ratios'].values())
    for d in scale['absent_terms'].values():                           # Sec. 3.5
        assert all(float(x) < 1e-10 for c in d.values() for x in c.values())
    for b, d in scale['three_scales'].items():
        assert all(float(x['spread']) < 1e-10 for k, x in d.items() if 'spread' in x and k != 'spin0 G^1')
        assert float(d['spin0 G^1 slope in L']['rel_dev']) < 1e-10
    assert scale['stored_points']['compared'] == 28 and float(scale['stored_points']['max_rel_diff']) < 1e-15
    assert all(float(c['rel_diff']) < 1e-20 for c in scale['precision_control'])
    for b, v in fwd.items():                                           # Sec. 2.5
        assert float(v['rel_dev']) < 1e-17
        assert float(v['dr_dpT2_spread']) < 1e-4


def test_one_loop_scripts():
    run('verify_dps_one_loop.py', ROOT/'one_loop')
    out = run('hh_spin2_numerator.py', ROOT/'one_loop')
    assert 'FAIL' not in out
    run('landau_boxes.py', ROOT/'one_loop')


def test_stage2a():
    stored = (ROOT/'stage2a'/'results_stage2a.json').read_text()
    out = run('stage2a_bookkeeping.py', ROOT/'stage2a')
    assert out.count('give the same O(eps^0) terms: True') == 2
    new = json.loads((ROOT/'stage2a'/'results_stage2a.json').read_text())
    assert new['rows'] == json.loads(stored)['rows'], 'Table 4 entries changed'


def test_tables_reproduced():
    # Regenerate the standalone PRL tables and require exactly the stored values.
    paper = ROOT.parent/'manuscript_prl'
    files = [paper/'numbers.tex'] + sorted((paper/'tables').glob('*.tex'))
    stored = {f: f.read_text() for f in files}
    run('make_tables.py', paper)
    for f, text in stored.items():
        assert f.read_text() == text, f'{f.relative_to(ROOT.parent)} changed'


def test_rg_argument():
    import sympy as sp
    run('src/verify_integrated_setup.py', ROOT/'rg_argument')
    new = json.loads((ROOT/'rg_argument'/'results'/'verification.json').read_text())
    ref = json.loads((ROOT/'rg_argument'/'results'/'verification_reference.json').read_text())
    equal = lambda a, b: sp.simplify(sp.sympify(a)-sp.sympify(b)) == 0
    for key in ('A0_over_C1', 'A2_over_K1'):
        assert equal(new[key], ref[key]), key
    for key in ('gluon_channel', 'boundary_shift'):
        for k, v in ref[key].items():
            assert v == new[key][k] if k == 'note' else equal(new[key][k], v), (key, k)
    for row, v in ref['representation_table'].items():
        assert all(equal(new['representation_table'][row][k], v[k]) for k in v), row


def test_orientation():
    run('orientation_check.py', ROOT/'orientation')


def test_aa_curve():
    cands = [ROOT.parents[1]/'common'/'aa_source', ROOT.parents[1]/'paper1', ROOT.parents[1]]             # the master folder, or the P1 folder paper2/ sits in
    p1 = pathlib.Path(os.environ.get('P1') or next((c for c in cands if (c/'src'/'python'/'ggaa').is_dir()), cands[0])).resolve()
    if not (p1/'src'/'python'/'ggaa').is_dir() or not (p1/'inputs'/'arXiv-1909.08993'/'anc').is_dir():
        skip(f'P1 package with its inputs not found at {p1}; set P1 to its top directory')
    with tempfile.TemporaryDirectory() as tmp:
        new = pathlib.Path(tmp)/'results_aa_curve.json'
        run('aa_full_range.py', ROOT/'one_loop', str(p1), str(new))
        same(json.loads(new.read_text()), json.loads((ROOT/'one_loop'/'results_aa_curve.json').read_text()),
             path='results_aa_curve.json')


def test_published_expressions():
    anc = HH/'anc'
    if not all((anc/n).exists() for n in CHECKSUMS):
        skip('ancillary files of arXiv:1809.05388 not present; run fetch_inputs.py')
    for n, h in CHECKSUMS.items():
        assert hashlib.sha256((anc/n).read_bytes()).hexdigest() == h, n
    code = r'''
import sys, mpmath as mp
sys.path.insert(0, ".")
import evaluate_hh as ev, scan_colour_scale as sc
mp.mp.dps = 36
x, y = mp.mpf(1)/3, mp.mpf(1)/2
# raw Gauge_i at x = 1/3, y = 1/2, mu_R^2 = m_H^2/2, N = 3, n_f = 5, as stored by an earlier
# run of this evaluator at 36 digits (bose_symmetry.json; its at_z entries are the same
# amplitudes with t and u exchanged, which is the Bose-symmetry check quoted in the paper)
import json
bose = json.load(open("bose_symmetry.json"))
ref = {g: mp.mpmathify(bose[str(g)]["at_y"].replace(" ", "")) for g in (1, 2)}
w = {"N2": 9, "N_Nf": 15, "Nf_over_N": mp.mpf(5)/3, "Nf2": 25}
def flush():   # the polylog caches are keyed on arguments, not on precision
    for f in (ev.li22, ev.LiClassical, ev.Log, ev.Zeta):
        f.cache_clear()
for g in (1, 2):
    flush()
    a = ev.eval_amp(ev.compile_amp(g), x, y, 1, mp.mpf("0.5"), 3, 5)
    assert abs(a-ref[g]) < mp.mpf(10)**-30*abs(ref[g]), (g, a)
    codes = sc.split_terms(g)
    names = sc.identify(codes)          # works at 20 digits
    mp.mp.dps = 36
    flush()
    env = sc.env(x, y, mp.mpf("0.5"), 1, 1)
    split = sum(w[n]*eval(c, {"__builtins__": {}}, env) for n, c in zip(names, codes))
    assert abs(split-a) < mp.mpf(10)**-30*abs(a), (g, split, a)
print("ok")
'''
    out = subprocess.run([PY, '-c', code], cwd=HH, capture_output=True, text=True)
    assert out.returncode == 0 and out.stdout.strip().endswith('ok'), out.stdout + out.stderr


if __name__ == '__main__':
    full = '--full' in sys.argv
    tests = [test_comparison_precision, test_scan_reports_failures,
             test_points_complete, test_analyses_reproduce_stored_results, test_statements_of_the_paper,
             test_one_loop_scripts, test_stage2a, test_tables_reproduced, test_rg_argument, test_orientation]
    if full:
        tests += [test_published_expressions, test_aa_curve]
    failed = 0
    for t in tests:
        try:
            t()
            print(f'ok      {t.__name__}', flush=True)
        except Skipped as e:
            print(f'skipped {t.__name__}: {e}', flush=True)
        except BaseException as e:
            if type(e).__name__ == 'Skipped':
                print(f'skipped {t.__name__}: {e}', flush=True)
                continue
            failed += 1
            print(f'FAIL    {t.__name__}: {e}', flush=True)
    if not full:
        print('(the published expressions themselves are tested with --full, after fetch_inputs.py)')
    sys.exit(1 if failed else 0)
