#!/usr/bin/env python3
"""Fits and comparisons for every process whose points are present; writes results_PROCESS.json
next to this package and prints a summary.

For gluon splittings (hh, aa) each velocity gets
    leading_logs   the G^3, G^2 (M_1) and G^2, G^1 (M_2) coefficients of the colour monomials at
                   mu = m, with the window shift, the prediction of eq. (5.5) and the deviation
    r0             r_0 of both components for the C_A and T_F n_f parts: value, shift, closed form
                   (7.2), deviation; and the n_f^2 and n_f/N terms of the same coefficients over
                   the one-loop term, which must vanish
    one_loop       (aa only) C_1 and K fitted from the one-loop files against 2 pi N F_A and pi N W
For quark splittings (zz) the rows of analysis.leading_logs_quark and r0_quark, with the closed
form (7.3).  Numbers are strings with 30 significant digits; deviations are floats.

    python3 -m dps.report [hh|aa|zz ...]
"""
import json
import pathlib
import sys

import mpmath as mp

HERE = pathlib.Path(__file__).resolve().parent
PKG = HERE.parent
sys.path.insert(0, str(PKG))
from dps import analysis, fit, kinematics, layout, points, theory      # noqa: E402

POINTS = {'hh': layout.path('hh_points'), 'zz': layout.path('zz_points'), 'aa': layout.AA_POINTS}
mp.mp.dps = 50


def S(z, n=30):
    """Complex number as two 30-digit strings, real number as one."""
    if z is None:
        return None
    if isinstance(z, mp.mpc):
        return [mp.nstr(mp.re(z), n), mp.nstr(mp.im(z), n)]
    return mp.nstr(z, n)


def F(z):
    return float(z) if z is not None else None


def gluon(process, regular, control):
    out = {}
    for beta in analysis.betas_of(regular):
        if not all((c, beta) in regular and len(regular[(c, beta)]) >= 7 for c in (1, 2)):
            continue
        x = kinematics.momentum_fraction(beta)
        res = {'x': S(x), 'n_points': {c: len(regular[(c, beta)]) for c in (1, 2)},
               'deltas': {c: [r['delta_exp'] for r in regular[(c, beta)]] for c in (1, 2)}}
        res['leading_logs'] = [{'component': r['component'], 'power': r['power'], 'colour': r['colour'],
                                'value': S(r['value']), 'shift': S(r['shift']), 'prediction': S(r['prediction']),
                                'deviation': F(r['deviation'])} for r in analysis.leading_logs_gluon(process, regular, beta)]
        r0 = analysis.r0_gluon(process, regular, beta)
        closed = analysis.r0_closed_gluon(process, beta)
        res['r0'] = {}
        for nu in ('nu0', 'nu2'):
            e = r0[nu]
            res['r0'][nu] = {
                'N': {'value': S(e['r0_N']), 'shift': S(e['N2']['shift']), 'closed': S(closed[nu]['N']),
                      'deviation': F(abs(e['r0_N'] - closed[nu]['N'])), 'imag': F(abs(mp.im(e['r0_N'])))},
                # the N n_f term over N gives the coefficient of n_f; eq. (7.2) and Table 9 use T_F n_f = n_f/2
                'TFnf': {'value': S(2 * e['r0_nf']), 'shift': S(2 * e['N_Nf']['shift']), 'closed': S(2 * closed[nu]['nf']),
                         'deviation': F(2 * abs(e['r0_nf'] - closed[nu]['nf'])), 'imag': F(2 * abs(mp.im(e['r0_nf'])))},
                'Nf2_over_oneloop': F(abs(e['Nf2']['over_oneloop'])),
                'Nf_over_N_over_oneloop': F(abs(e['Nf_over_N']['over_oneloop'])),
                'known_N': S(e['known_N']), 'known_nf': S(e['known_nf'])}
        if process.name == 'aa':
            # the two components against gg -> HH (Sec. 7.1): nu = 0 differs by the helicity term with F_A in
            # place of F_S, nu = 2 by the sign of the helicity-flip term, -(4/3)(C_A - 2 T_F n_f)/W
            hhc = analysis.r0_closed_gluon(theory.HH, beta)
            res['r0']['against_hh'] = {
                'nu0': {'N': F(abs(r0['nu0']['r0_N'] - hhc['nu0']['N']
                                   - mp.mpf(2) / 3 * (1 / theory.F_A(x) - 1 / theory.F_S(x)))),
                        'TFnf': F(2 * abs(r0['nu0']['r0_nf'] - hhc['nu0']['nf']
                                          + mp.mpf(2) / 3 * (1 / theory.F_A(x) - 1 / theory.F_S(x))))},
                'nu2': {'N': F(abs(r0['nu2']['r0_N'] - hhc['nu2']['N'] + mp.mpf(4) / 3 / theory.W(x))),
                        'TFnf': F(2 * abs(r0['nu2']['r0_nf'] - hhc['nu2']['nf'] - mp.mpf(4) / 3 / theory.W(x)))}}
            res['one_loop'] = {}
            for comp, idx, target, name in ((1, 1, 2 * mp.pi * theory.F_A(x), 'C1/N'), (2, 0, mp.pi * theory.W(x), 'K/N')):
                recs = regular[(comp, beta)]
                oneN, sN = fit.with_shift(recs, points.series(recs, '1l:N'), 1, 2)
                onef, sf = fit.with_shift(recs, points.series(recs, '1l:Nf'), 1, 2)
                res['one_loop'][name] = {'fit': S(oneN[idx]), 'shift': S(sN[idx]), 'analytic': S(target),
                                         'deviation': F(abs(oneN[idx] / target - 1)),
                                         'nf_term_over_fit': F(abs(onef[idx] / oneN[idx]))}
                res['one_loop'][f'{"c10" if comp == 1 else "K_slope"}/N'] = S(oneN[1 - idx])
        out[beta] = res
    ctrl = []
    for r in control:
        base = [q for q in regular.get((r['component'], r['beta']), []) if q['delta_exp'] == r['delta_exp']]
        if base:
            q = base[0]
            # compared at the full precision of the stored digits, not at the 50 digits of the analysis,
            # at which the two evaluations coincide
            worst = points.compare_stored(POINTS[process.name], r['file'], q['file'])
            ctrl.append({'component': r['component'], 'beta': r['beta'], 'delta_exp': r['delta_exp'], 'dps': r['dps'],
                         'base_dps': q['dps'], 'max_rel_diff': F(worst)})
    return {'betas': out, 'control': ctrl}


def quark(regular, control):
    out = {}
    for beta in analysis.betas_of(regular):
        x = kinematics.momentum_fraction(beta)
        rows = analysis.leading_logs_quark(regular, beta)
        if not rows:
            continue
        res = {'x': S(x), 'leading_logs': [{k: (S(v) if isinstance(v, (mp.mpc, mp.mpf)) and k not in ('deviation',) else F(v) if k == 'deviation' else v)
                                            for k, v in r.items()} for r in rows]}
        r0 = analysis.r0_quark(regular, beta)
        CFc, CAc = theory.r0_closed(theory.ZZ, 0, x)
        res['r0'] = {}
        for nu in ('nu0', 'nu2'):
            res['r0'][nu] = {}
            for h, row in r0[nu].items():
                res['r0'][nu][h] = {part: {'value': S(row[part]['r0']), 'shift': S(row[part]['shift']),
                                           'closed': S(c), 'deviation': F(abs(row[part]['r0'] - c))}
                                    for part, c in (('2l_CF', CFc), ('2l_CA', CAc))}
        out[beta] = res
    return {'betas': out}


def summary(name, res):
    print(f'== {name}')
    for beta, r in res['betas'].items():
        if 'leading_logs' in r and name != 'zz':
            dev = max(e['deviation'] for e in r['leading_logs'] if e['prediction'] not in ('0.0', None) and mp.mpf(e['prediction']) != 0) \
                if any(mp.mpf(e['prediction']) != 0 for e in r['leading_logs']) else 0
            print(f'  beta = {beta}: leading logs within {dev:.1e} of eq. (5.5)')
            for nu in ('nu0', 'nu2'):
                e = r['r0'][nu]
                print(f'     r0({nu}): C_A part {e["N"]["value"][0][:14]} (closed {e["N"]["closed"][:14]}, dev {e["N"]["deviation"]:.1e}, Im {e["N"]["imag"]:.1e}); '
                      f'T_F n_f part {e["TFnf"]["value"][0][:14]} (closed {e["TFnf"]["closed"][:14]}, dev {e["TFnf"]["deviation"]:.1e})')
                if 'against_hh' in r['r0']:
                    a = r['r0']['against_hh'][nu]
                    rel = '(2/3)(C_A - 2 T_F n_f)(1/F_A - 1/F_S)' if nu == 'nu0' else '-(4/3)(C_A - 2 T_F n_f)/W'
                    print(f'            minus gg -> HH: {rel} to {max(a["N"], a["TFnf"]):.1e}')
            if 'one_loop' in r:
                o = r['one_loop']
                print(f'     one loop: C1/N dev {o["C1/N"]["deviation"]:.1e}, K/N dev {o["K/N"]["deviation"]:.1e}')
        elif name == 'zz':
            dev = max(e['deviation'] for e in r['leading_logs'])
            worst = max(v[p]['deviation'] for nu in r['r0'].values() for v in nu.values() for p in v)
            print(f'  beta = {beta}: leading logs within {dev:.1e}; r0 within {worst:.1e} of eq. (7.3)')


def main(names):
    for name in names:
        d = POINTS[name]
        if not d.is_dir() or not any(d.glob('*.json')):
            print(f'== {name}: no points')
            continue
        regular, control = points.LOADERS[name](d)
        res = gluon(theory.PROCESSES[name], regular, control) if name != 'zz' else quark(regular, control)
        if not res['betas']:
            print(f'== {name}: not enough points for a fit yet')
            continue
        (PKG / f'results_{name}.json').write_text(json.dumps(res, indent=1))
        summary(name, res)


if __name__ == '__main__':
    main(sys.argv[1:] or ['hh', 'zz', 'aa'])
