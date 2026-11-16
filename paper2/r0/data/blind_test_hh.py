#!/usr/bin/env python3
"""Blind test: r0 at velocities that were not used to find the closed form (points_hh/, made by
run_hh_new_velocities.py), extracted exactly as in r0_from_data.py and compared with r0_closed_forms.py.
    python3 blind_test_hh.py
"""
import json, pathlib
import mpmath as mp
import r0_from_data as RD
import r0_closed_forms as CF

HERE = pathlib.Path(__file__).resolve().parent


def main():
    out = RD.extract(HERE / 'points_hh')
    rows = []
    for beta, r in out.items():
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
            print(f"beta={beta} x={mp.nstr(x, 6)} {comp}: C_A part data {mp.nstr(a, 12)}  formula {mp.nstr(pa, 12)}"
                  f"  |diff| {mp.nstr(abs(a - pa), 2)} (shift {mp.nstr(mp.mpf(r[comp]['r0_N'][2]), 2)})")
            print(f"                 n_f part data {mp.nstr(f, 12)}  formula {mp.nstr(pf, 12)}"
                  f"  |diff| {mp.nstr(abs(f - pf), 2)} (shift {mp.nstr(mp.mpf(r[comp]['r0_nf'][2]), 2)})")
    json.dump(rows, open(HERE / 'blind_test_hh.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
