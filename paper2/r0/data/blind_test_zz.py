#!/usr/bin/env python3
"""Blind test for quark splittings: r0 at beta = 4/5 (x = 0.9), a velocity not used to find the closed form.
The points (points_zz/) were made with the unchanged task-J code (zz_points.py); the fit is that of
analyse_zz.py, and r0 is extracted as in r0_from_zz.py and compared with r0_closed_forms.r0_quark.
    python3 blind_test_zz.py
"""
import json, sys, pathlib
import mpmath as mp

HERE = pathlib.Path(__file__).resolve().parent
TASKS = HERE.parents[1]
sys.path.insert(0, str(TASKS / 'zz' / 'scripts'))
import analyse_zz as AZ                      # noqa: E402
import r0_from_zz as RZ                      # noqa: E402
import r0_closed_forms as CF                 # noqa: E402

mp.mp.dps = 40


def main():
    regular, control = AZ.load(HERE / 'points_zz')
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
            if e['dJz'] == 0 and abs(c1) > 1e-3*cmax:
                comp, norm, j = 'nu0', c1, 1
            elif abs(e['dJz']) == 2 and abs(K) > 1e-6:
                comp, norm, j = 'nu2', K, 0
            else:
                continue
            row = {'beta': beta, 'x': float(x), 'hel': h, 'component': comp}
            for part, (cf, ca), pred in (('2l_CF', (1, 0), pF), ('2l_CA', (0, 1), pA)):
                f = RZ.C(e[part]['f'][j])
                fw = RZ.C(e[part]['f_w1'][j]) if e[part]['f_w1'] else None
                kn = RZ.known(x, cf, ca)[0 if comp == 'nu0' else 1]
                r0 = 2*f/norm - kn
                row[part] = {'r0': [float(mp.re(r0)), float(mp.im(r0))], 'formula': float(pred),
                             'dev': float(abs(r0 - pred)),
                             'shift': float(2*abs(f - fw)/abs(norm)) if fw is not None else None}
            rows.append(row)
            print(f"beta={beta} x={mp.nstr(x, 4)} {comp} {h}: C_F {row['2l_CF']['r0'][0]:.10f} (formula {row['2l_CF']['formula']:.10f},"
                  f" |diff| {row['2l_CF']['dev']:.1e}, shift {row['2l_CF']['shift']})   C_A {row['2l_CA']['r0'][0]:.10f}"
                  f" (formula {row['2l_CA']['formula']:.10f}, |diff| {row['2l_CA']['dev']:.1e}, shift {row['2l_CA']['shift']})")
    json.dump(rows, open(HERE / 'blind_test_zz.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
