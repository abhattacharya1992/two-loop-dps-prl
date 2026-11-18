#!/usr/bin/env python3
"""Endpoint test of the r0 closed forms: gg -> HH and gg -> ZZ at beta = 19/20 (x = 0.975).

At x = 0.975 the x-dependent part of the quark-splitting r0 (Eq. (8) of the Letter, N_c = 3) is 2.4 times its
value at x = 0.9, the largest momentum fraction of the other tests. The points in points_endpoint/ were computed on
17 November 2026:
  hh/              scan_colour_scale.run_job of ../../hh/hh_two_loop, unchanged: delta = 1e-4 ... 1e-11 with
                   30 + 3k digits, and delta = 1e-10 with 15 more digits (precision control);
  zz/              zz_points_acc.py (this folder): ../../zz/scripts/zz_points.py unchanged, except that
                   Li_{2,2} on its integral branch is evaluated by ../../hh/hh_two_loop/evaluate_hh.li22;
                   delta = 1e-4 ... 1e-12 with 30 + 5k digits and delta = 1e-10 with 95 digits; li22_*.json
                   record, for each point, how far the fixed subdivision of ggvv.py is from that evaluation;
  zz_fixed_split/  the same ZZ points with zz_points.py unchanged.
The points were made, from the repository root, with
  gg -> HH  scan_colour_scale.run_job((g, '19/20', k, 30 + 3*k), out_dir=OUT) for g = 1, 2 and k = 4 ... 11, and
            run_job((g, '19/20', 10, 75), out_dir=OUT) for the controls (paper2/hh/hh_two_loop on sys.path;
            the input files of `make -C paper2/hh fetch`);
  gg -> ZZ  python3 paper2/r0/data/zz_points_acc.py $PWD BUILD OUT 19/20 K DPS for (K, DPS) = (k, 30 + 5k),
            k = 4 ... 12, and (10, 95); BUILD is the translated ggvvamp code, made by
            paper2/zz/scripts/translate_ggvvamp.py from the unpacked ggvvamp-1.0.tar.gz (see paper2/zz/README.md).
            zz_points.py with the same arguments, without REPO, gives zz_fixed_split/.
r0 is extracted as in blind_test_hh.py and blind_test_zz.py. For ZZ the amplitudes with Delta J_z = 0 enter
when |c_1| exceeds 1e-6 (not 1e-3) of the largest one, so that ++-- and --++, whose one-loop logarithm is
(1-x)^2/x^2 = 1/1521 of that of ++++, are included. The precision controls are compared as in
analyse_scale.py (N_c = 3, n_f = 5) and analyse_zz.py.
    python3 endpoint_test.py        (a few seconds) -> endpoint_test.json
"""
import json, math, pathlib, sys
import mpmath as mp

HERE = pathlib.Path(__file__).resolve().parent
TASKS = HERE.parents[1]
sys.path.insert(0, str(TASKS / 'zz' / 'scripts'))
import r0_from_data as RD                    # noqa: E402  (puts ../../hh/hh_two_loop on the path)
import r0_closed_forms as CF                 # noqa: E402
import r0_from_zz as RZ                      # noqa: E402
import analyse_zz as AZ                      # noqa: E402
from scan_tools import load as load_hh, cplx, COLOUR_WEIGHT   # noqa: E402

POINTS = HERE / 'points_endpoint'


def hh_part(folder):
    mp.mp.dps = 50
    rows = []
    for beta, r in RD.extract(folder).items():
        x = mp.mpf(r['x'])
        for comp, nu in (('nu0', 0), ('nu2', 2)):
            pa, pf = CF.r0_gluon(x, nu)
            a = mp.mpc(mp.mpf(r[comp]['r0_N'][0]), mp.mpf(r[comp]['r0_N'][1]))
            f = mp.mpc(mp.mpf(r[comp]['r0_nf'][0]), mp.mpf(r[comp]['r0_nf'][1]))
            rows.append({'beta': beta, 'x': float(x), 'component': comp,
                         'C_A': [float(mp.re(a)), float(mp.im(a))], 'C_A_formula': float(pa),
                         'n_f': [float(mp.re(f)), float(mp.im(f))], 'n_f_formula': float(pf),
                         'dev_C_A': float(abs(a - pa)), 'dev_n_f': float(abs(f - pf)),
                         'shift_C_A': float(mp.mpf(r[comp]['r0_N'][2])), 'shift_n_f': float(mp.mpf(r[comp]['r0_nf'][2]))})
    reg, control = load_hh(folder)
    total = lambda rec, L: sum(COLOUR_WEIGHT[c](3, 5) * cplx(rec['A'][c][L]) for c in COLOUR_WEIGHT)
    ctrl = []
    for c in control:
        twin = [q for q in reg[(c['gauge'], c['beta'])] if q['delta_exp'] == c['delta_exp']][0]
        d = max(abs(total(c, L) - total(twin, L)) / abs(total(c, L)) for L in ('-1', '0', '1'))
        ctrl.append({'gauge': c['gauge'], 'beta': c['beta'], 'delta_exp': c['delta_exp'],
                     'dps': [twin['dps'], c['dps']], 'rel_diff': float(d)})
    return rows, ctrl


def zz_part(folder, fixed):
    mp.mp.dps = 40                                    # as blind_test_zz.py
    regular, control = AZ.load(folder)
    res = AZ.analyse(regular)
    rows = []
    for beta, r in res.items():
        x = mp.mpf(r['x'])
        pF, pA = CF.r0_quark(x)
        cmax = max(abs(RZ.C(e['c1'])) for e in r['hel'].values())
        for h, e in r['hel'].items():
            if '2l_CF' not in e:
                continue
            c1, K = RZ.C(e['c1']), RZ.C(e['K'])
            if e['dJz'] == 0 and abs(c1) > 1e-6 * cmax:
                comp, norm, j = 'nu0', c1, 1
            elif abs(e['dJz']) == 2 and abs(K) > 1e-6:
                comp, norm, j = 'nu2', K, 0
            else:
                continue
            row = {'beta': beta, 'x': float(x), 'hel': h, 'component': comp}
            for part, (cf, ca), pred in (('2l_CF', (1, 0), pF), ('2l_CA', (0, 1), pA)):
                f = RZ.C(e[part]['f'][j])
                fw = RZ.C(e[part]['f_w1'][j])
                kn = RZ.known(x, cf, ca)[0 if comp == 'nu0' else 1]
                r0 = 2 * f / norm - kn
                row[part] = {'r0': [float(mp.re(r0)), float(mp.im(r0))], 'formula': float(pred),
                             'dev': float(abs(r0 - pred)), 'shift': float(2 * abs(f - fw) / abs(norm))}
            rows.append(row)
    with mp.workdps(50):                              # as analyse_zz.py
        ctrl = AZ.controls(regular, control)
    side = [json.loads(p.read_text()) for p in sorted(folder.glob('li22_*.json'))]
    lost = max(side, key=lambda s: float(s['max_abs_diff_fixed_split']) * 10 ** s['dps'])
    li22 = {'points': len(side),
            'max_digits_lost_fixed_split': round(math.log10(float(lost['max_abs_diff_fixed_split'])) + lost['dps']),
            'at': {'delta_exp': lost['delta_exp'], 'dps': lost['dps']},
            'max_reference_diff_in_units_of_last_digit': max(float(s['abs_diff_to_reference_there']) * 10 ** s['dps']
                                                              for s in side)}
    # largest change of an amplitude, relative to the largest helicity amplitude of the same point and part
    change = 0.0
    with mp.workdps(120):
        for p in sorted(folder.glob('zz_*.json')):
            a, b = json.loads(p.read_text()), json.loads((fixed / p.name).read_text())
            for part in AZ.PARTS:
                big = max(abs(AZ.cplx(v[part])) for v in a['M'].values())
                d = max(abs(AZ.cplx(a['M'][h][part]) - AZ.cplx(b['M'][h][part])) for h in a['M'])
                change = max(change, float(d / big))
    return rows, ctrl, li22, change


def quark_x_part(x):
    lx, ly = mp.log(x), mp.log(1 - x)
    return mp.mpf(4) / 3 * 2 * (lx - ly) ** 2 + 3 * (-2 * lx ** 2 - 2 * ly ** 2)


def main():
    hh_rows, hh_ctrl = hh_part(POINTS / 'hh')
    zz_rows, zz_ctrl, li22, change = zz_part(POINTS / 'zz', POINTS / 'zz_fixed_split')
    mp.mp.dps = 30
    ratio = quark_x_part(mp.mpf('0.975')) / quark_x_part(mp.mpf('0.9'))
    out = {'beta': '19/20', 'x': 0.975, 'quark_x_dependence_ratio_to_x_0.9': float(ratio),
           'hh': {'rows': hh_rows, 'control': hh_ctrl},
           'zz': {'rows': zz_rows, 'control': zz_ctrl, 'li22': li22, 'max_amplitude_change_fixed_split': change}}
    json.dump(out, open(HERE / 'endpoint_test.json', 'w'), indent=1)
    for r in hh_rows:
        print(f"HH {r['component']}: C_A |diff| {r['dev_C_A']:.1e} (shift {r['shift_C_A']:.1e}),"
              f" n_f |diff| {r['dev_n_f']:.1e} (shift {r['shift_n_f']:.1e})")
    for r in zz_rows:
        print(f"ZZ {r['hel']} {r['component']}: " + ', '.join(
            f"{p[3:]} |diff| {r[p]['dev']:.1e} (shift {r[p]['shift']:.1e})" for p in ('2l_CF', '2l_CA')))
    print('controls:', [c['rel_diff'] for c in hh_ctrl], [c['max_rel_diff'] for c in zz_ctrl])
    print('Li22:', li22, '; largest amplitude change from the fixed subdivision:', f'{change:.1e}')
    print(f'x-dependent part of the quark r0, x = 0.975 against x = 0.9: {float(ratio):.3f}')


if __name__ == '__main__':
    main()
