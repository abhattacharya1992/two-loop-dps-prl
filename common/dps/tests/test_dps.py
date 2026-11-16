#!/usr/bin/env python3
"""The shared code against the stored results of Paper 2.

1. kinematics: the x, y, p_T^2 and G stored in every gg -> HH and gg -> ZZ point file are reproduced
   to the last digit at the precision of the file;
2. gg -> HH leading logarithms (hh/hh_two_loop/results_colour.json) to 1e-45 and r_0
   (r0/data/r0_from_data.json) to the last digit;
3. gg -> ZZ leading logarithms (zz/results_zz.json, 30-digit strings) to 1e-28 and r_0
   (r0/data/r0_from_zz.json) to 1e-13 (that script took x from a float);
4. the Laurent polynomials of the colour decomposition against numerical evaluation;
5. all 80 regular gg -> AA points and the two precision controls are present; the fitted one-loop
   coefficients C_1 and K, the two-loop logarithms, the absent colour terms and r_0 agree with the
   analytic formulae. These checks fit the point records directly, without reading results_aa.json.
6. the HH and ZZ backend paths resolve to their evaluation scripts in the current layout.

    python3 tests/test_dps.py        (a few seconds)
"""
import json
from fractions import Fraction
import os
import pathlib
import sys

import mpmath as mp

HERE = pathlib.Path(__file__).resolve().parent
PKG = HERE.parent
sys.path.insert(0, str(PKG))
from dps import analysis, fit, kinematics, layout, points, theory      # noqa: E402
from dps.colour import N, NF                                   # noqa: E402

PAPER2 = PKG.parent
mp.mp.dps = 50


def C(pair):
    return mp.mpc(mp.mpf(pair[0]), mp.mpf(pair[1]))


def test_backend_paths():
    from dps.backends import hh, zz
    assert (hh.HH / 'scan_colour_scale.py').is_file(), hh.HH
    assert (zz.ZZ / 'scripts' / 'zz_points.py').is_file(), zz.ZZ
    assert (zz.ZZ / 'scripts' / 'ggvv.py').is_file(), zz.ZZ
    print('backends: HH and ZZ evaluation scripts found in the current layout')


def test_kinematics():
    n = 0
    for p in sorted(layout.path('hh_points').glob('*.json')):
        r = json.loads(p.read_text())
        mp.mp.dps = r['dps']
        k = kinematics.point(r['beta'], r['delta_exp'])
        assert mp.nstr(k['xk'], r['dps']) == r['x'] and mp.nstr(k['y'], r['dps']) == r['y'], p.name
        assert mp.nstr(k['pT2'], r['dps']) == r['pT2'], p.name
        assert [mp.nstr(mp.re(k['G']), r['dps']), mp.nstr(mp.im(k['G']), r['dps'])] == r['G'], p.name
        n += 1
    for p in sorted(layout.path('zz_points').glob('zz_*.json')):
        r = json.loads(p.read_text())
        mp.mp.dps = r['dps']
        k = kinematics.point(r['beta'], r['delta_exp'])
        assert mp.nstr(k['xk'], r['dps']) == r['x'] and mp.nstr(k['t'], r['dps']) == r['t'], p.name
        # zz_points.py computes z = s - 2 - y instead of xk + 1/xk - y: equal up to rounding at the
        # last digits of the working precision, so p_T^2 and G are compared to 10^-(dps - 12)
        tol = mp.mpf(10) ** (12 - r['dps'])
        assert abs(k['pT2'] / mp.mpf(r['pT2']) - 1) < tol, p.name
        assert abs(k['G'] - mp.mpc(mp.mpf(r['G'][0]), mp.mpf(r['G'][1]))) < tol, p.name
        n += 1
    mp.mp.dps = 50
    print(f'kinematics: {n} point files reproduced (gg -> HH to the last digit, gg -> ZZ to the rounding of z)')


def test_hh():
    reg, _ = points.load_hh(layout.path('hh_points'))
    stored = json.load(open(layout.path('hh_colour')))
    r0s = json.load(open(layout.path('r0_hh')))
    name = {(1, 3): 'spin0 G^3', (1, 2): 'spin0 G^2', (2, 2): 'spin2 G^2', (2, 1): 'spin2 G^1'}
    worst, n = mp.mpf(0), 0
    for beta in ('1/4', '1/2', '2/3'):
        for row in analysis.leading_logs_gluon(theory.HH, reg, beta):
            s, = [q for q in stored[beta]['rows']
                  if q['coefficient'] == name[row['component'], row['power']] and q['colour'] == row['colour']]
            worst = max(worst, abs(row['value'] - C((s['re'], s['im']))) / max(abs(row['value']), 1))
            n += 1
        r0 = analysis.r0_gluon(theory.HH, reg, beta)
        for nu in ('nu0', 'nu2'):
            for key in ('r0_N', 'r0_nf'):
                a, s = r0[nu][key], r0s[beta][nu][key]
                assert (str(mp.re(a)), str(mp.im(a))) == (s[0], s[1]), (beta, nu, key)
    assert worst < mp.mpf(10) ** -45, worst
    print(f'gg -> HH: {n} leading-log coefficients to {mp.nstr(worst, 2)}, r_0 identical to r0_from_data.json')


def test_zz():
    reg, _ = points.load_zz(layout.path('zz_points'))
    st = json.load(open(layout.path('zz_results')))['betas']
    r0z = json.load(open(layout.path('r0_zz')))
    worst, worst_r0, n = mp.mpf(0), mp.mpf(0), 0
    for beta in ('1/4', '1/2', '2/3'):
        for row in analysis.leading_logs_quark(reg, beta):
            e = st[beta]['hel'][row['helicity']]
            f = C(e[row['part']]['f'][row['power']])
            norm = C(e['c1']) if row['kind'] == 'dJz0' else C(e['K'])
            worst = max(worst, abs(row['ratio'] - f / norm) / abs(f / norm))
            n += 1
        r0 = analysis.r0_quark(reg, beta)
        for nu in ('nu0', 'nu2'):
            for h, row in r0[nu].items():
                for part in ('2l_CF', '2l_CA'):
                    worst_r0 = max(worst_r0, abs(row[part]['r0'] - C(r0z[beta][nu][h][part])))
    assert worst < mp.mpf(10) ** -28 and worst_r0 < mp.mpf(10) ** -13, (worst, worst_r0)
    print(f'gg -> ZZ: {n} leading-log ratios to {mp.nstr(worst, 2)}, r_0 to {mp.nstr(worst_r0, 2)} of r0_from_zz.json')


def test_laurent():
    sys.path.insert(0, str(layout.p1_root() / 'src' / 'python'))
    try:
        from ggaa import twoloop as tl
    except ImportError as exc:
        raise AssertionError(f'P1 evaluator unavailable at {layout.p1_root()}; set P1_ROOT if needed') from exc
    mp.mp.dps = 30
    base = {'x': mp.mpf('0.3'), 'z': mp.mpf('2.1'), 'Y': mp.mpf('1.4'), 'Z': mp.mpf('3.1'), 'M2': mp.mpf(1),
            'muR2': mp.mpf(1), 'Pi': +mp.pi, 'I': mp.mpc(0, 1)}
    tests = ['(N^2-1)/N*Nf/(2*N)', '-16*(-1 + N^2)*(6*M2*Nf*(1 + x)^4*(-I + Pi + I*Log[M2]))/(x*(1 + z))/N',
             '(4*M2*(1 + x)^2*((-27680*I)*N^3 - (1080*I)*Nf + (6200*I)*N^2*Nf + 15640*N^3*Pi))/N',
             'N^2*Nf^2/N^3 + 1/N^2*N - N^(-1)', '(2*N - Nf)^2/N', '-N*Nf*(-x)^2/(-N)',
             'N*PolyLog[2, 1 - z]*Nf - Log[Y]*N^2/Nf*Nf', '(N*(1+x) - Nf*(1-x))*(N*(1+x) + Nf*(1-x))/(N*Nf)']
    for t in tests:
        for Nc, nf in ((3, 5), (4, 2)):
            num = tl.Evaluator({**base, 'N': mp.mpf(Nc), 'Nf': mp.mpf(nf)}).evaluate(t)
            lau = tl.Evaluator({**base, 'N': N, 'Nf': NF}).evaluate(t)
            assert abs(num - lau.at(Nc, nf)) < mp.mpf(10) ** -24, t
    mp.mp.dps = 50
    print(f'Laurent polynomials: {len(tests)} expressions at two (N, n_f) agree with numerical evaluation')


def test_aa():
    d = layout.AA_POINTS
    assert d.is_dir() and any(d.glob('aa_*.json')), f'gg -> AA point records missing from {d}'
    reg, ctrl = points.load_aa(d)
    expected = {(comp, beta) for comp in (1, 2) for beta in ('1/10', '1/4', '1/2', '2/3', '4/5')}
    assert set(reg) == expected, f'gg -> AA scans differ: missing {expected - set(reg)}, extra {set(reg) - expected}'
    for key, recs in reg.items():
        assert [r['delta_exp'] for r in recs] == list(range(4, 12)), (key, 'incomplete or duplicate scan points')
    assert sorted((r['component'], r['beta'], r['delta_exp']) for r in ctrl) == [(1, '1/2', 10), (2, '1/2', 10)]
    for control in ctrl:
        reference, = [r for r in reg[control['component'], control['beta']]
                      if r['delta_exp'] == control['delta_exp']]
        assert control['dps'] > reference['dps'], control['file']
        assert points.compare_stored(d, control['file'], reference['file']) < mp.mpf('1e-20'), control['file']
    worst, fitted = mp.mpf(0), 0
    for (comp, beta), recs in sorted(reg.items()):
        x = kinematics.momentum_fraction(beta)
        fitted += 1
        one = fit.leading(recs, points.series(recs, '1l:N'), 1, 2)
        if comp == 1:
            worst = max(worst, abs(one[1] / (2 * mp.pi * theory.F_A(x)) - 1))
        else:
            worst = max(worst, abs(one[0] / (mp.pi * theory.W(x)) - 1))
        onef = fit.leading(recs, points.series(recs, '1l:Nf'), 1, 2)
        worst = max(worst, abs(onef[1 if comp == 1 else 0]) / abs(one[1 if comp == 1 else 0]))
    assert worst < mp.mpf(10) ** -8, worst
    print(f'gg -> AA: one-loop C_1 = 2 pi N F_A and K = pi N W, no n_f term, to {mp.nstr(worst, 2)} ({fitted} scans)')
    # Fit the two-loop point data afresh: cached comparison metrics must not certify changed inputs.
    worst_logs, worst_r0, worst_absent, n = mp.mpf(0), mp.mpf(0), mp.mpf(0), 0
    for beta in analysis.betas_of(reg):
        rows = analysis.leading_logs_gluon(theory.AA, reg, beta)
        assert len(rows) == 16, (beta, len(rows))
        worst_logs = max(worst_logs, *(r['deviation'] for r in rows))
        r0 = analysis.r0_gluon(theory.AA, reg, beta)
        closed = analysis.r0_closed_gluon(theory.AA, beta)
        for nu in ('nu0', 'nu2'):
            worst_r0 = max(worst_r0, abs(r0[nu]['r0_N'] - closed[nu]['N']),
                           2 * abs(r0[nu]['r0_nf'] - closed[nu]['nf']))
            worst_absent = max(worst_absent, abs(r0[nu]['Nf2']['over_oneloop']),
                               abs(r0[nu]['Nf_over_N']['over_oneloop']))
        n += 1
    assert worst_logs < mp.mpf(10) ** -8, worst_logs
    assert worst_r0 < mp.mpf(10) ** -8, worst_r0
    assert worst_absent < mp.mpf(10) ** -8, worst_absent
    print(f'gg -> AA: fresh two-loop fits, logarithm deviations {mp.nstr(worst_logs, 2)}, '
          f'r0 absolute error {mp.nstr(worst_r0, 2)}, absent colour terms {mp.nstr(worst_absent, 2)} '
          f'({n} velocities; results_aa.json is not used)')


def test_closed_forms():
    """eq. (7.2) for gg -> HH and eq. (7.4) for gg -> AA differ by the helicity term alone: through F_A in place
    of F_S for nu = 0, through the sign of the flip term for nu = 2"""
    for beta in ('1/10', '1/4', '1/2', '2/3', '4/5'):
        x = kinematics.momentum_fraction(Fraction(beta))
        for N, nf in ((1, 0), (0, 1)):
            c = mp.mpf(2) / 3 * (N - 2 * theory.TF * nf)
            d0 = theory.r0_closed(theory.AA, 0, x, N, nf) - theory.r0_closed(theory.HH, 0, x, N, nf)
            d2 = theory.r0_closed(theory.AA, 2, x, N, nf) - theory.r0_closed(theory.HH, 2, x, N, nf)
            assert abs(d0 - c * (1 / theory.F_A(x) - 1 / theory.F_S(x))) < mp.mpf(10) ** -25
            assert abs(d2 + 2 * c / theory.W(x)) < mp.mpf(10) ** -25
    print('closed forms: gg -> AA minus gg -> HH is the helicity term alone in both components (F_A for nu = 0, the sign for nu = 2)')


if __name__ == '__main__':
    test_backend_paths()
    test_kinematics()
    test_closed_forms()
    test_hh()
    test_zz()
    test_laurent()
    test_aa()
    print('passed')
