#!/usr/bin/env python3
"""Gate C, step 2: closed forms for r0 read off from the data, checked against every extracted value.

Gluon splittings (gg -> HH, r0_from_data.json; 'N' = C_A part, 'nf' = 2 T_F n_f part):
    r0 = C_A [134/9 - 5 zeta2 - 4 ln x ln(1-x)] - (40/9) T_F n_f + (2/3)(C_A - 2 T_F n_f)/F_nu(x),
    F_0 = F_S = 1/(x(1-x)) + 2x(1-x),  F_2 = W = (x^4 + (1-x)^4)/(x(1-x)).
Quark splittings (gg -> ZZ, r0_from_zz.json), the same for all eight helicity amplitudes:
    r0 = C_F [12 - 6 zeta2 + 2 ln^2(x/(1-x))] + C_A [2 + zeta2 - 2 (ln^2 x + ln^2(1-x))].
The constants 134/9 - 5 zeta2 and 12 - 6 zeta2 were identified with PSLQ on the basis {1, zeta2} (integer
coefficients up to 200, tolerance 1e-8; the constants are known to 3e-9 or better). At that precision the basis
cannot be enlarged: with ln 2 and zeta3 added, PSLQ gives them coefficient zero for the n_f and quark-splitting
constants, but for 134/9 - 5 zeta2 it finds relations with coefficients of order 100, and one exists among
1, zeta2, ln 2, zeta3 themselves (-44 + 6 zeta2 + 103 ln 2 - 31 zeta3 = 8e-10). The basis rests on the weight of
the one-loop calculation, which gives constants of weight at most two that do not depend on x, and the constants
on the regions calculation of Paper I. The coefficients of the logarithms and of 1/F_nu are the integers and fractions
that the three velocities give (exact solves, residuals ~1e-10).
    python3 r0_closed_forms.py
"""
import json, pathlib
import mpmath as mp

HERE = pathlib.Path(__file__).resolve().parent
mp.mp.dps = 30
z2 = mp.zeta(2)


def r0_gluon(x, nu):
    """returns (C_A coefficient, n_f coefficient) with n_f = 2 T_F n_f"""
    P = mp.log(x)*mp.log(1-x)
    F = 1/(x*(1-x)) + 2*x*(1-x) if nu == 0 else (x**4 + (1-x)**4)/(x*(1-x))
    return (mp.mpf(134)/9 - 5*z2 - 4*P + mp.mpf(2)/3/F, -mp.mpf(20)/9 - mp.mpf(2)/3/F)


def r0_quark(x):
    """returns (C_F coefficient, C_A coefficient)"""
    lx, ly = mp.log(x), mp.log(1-x)
    return (12 - 6*z2 + 2*(lx - ly)**2, 2 + z2 - 2*(lx**2 + ly**2))


def main():
    worst = {'HH': 0, 'ZZ': 0}
    hh = json.load(open(HERE / 'r0_from_data.json'))
    for beta, r in hh.items():
        x = mp.mpf(r['x'])
        for comp, nu in (('nu0', 0), ('nu2', 2)):
            pa, pf = r0_gluon(x, nu)
            da = abs(mp.mpc(mp.mpf(r[comp]['r0_N'][0]) - pa, mp.mpf(r[comp]['r0_N'][1])))
            df = abs(mp.mpc(mp.mpf(r[comp]['r0_nf'][0]) - pf, mp.mpf(r[comp]['r0_nf'][1])))
            worst['HH'] = max(worst['HH'], da, df)
            print(f'HH beta={beta} {comp}: C_A part {mp.nstr(pa, 12)} (data - formula {mp.nstr(da, 2)}), '
                  f'n_f part {mp.nstr(pf, 12)} ({mp.nstr(df, 2)})')
    zz = json.load(open(HERE / 'r0_from_zz.json'))
    for beta, r in zz.items():
        x = mp.mpf(r['x'])
        pF, pA = r0_quark(x)
        for comp in ('nu0', 'nu2'):
            for h, row in r[comp].items():
                dF = abs(mp.mpc(mp.mpf(row['2l_CF'][0]) - pF, mp.mpf(row['2l_CF'][1])))
                dA = abs(mp.mpc(mp.mpf(row['2l_CA'][0]) - pA, mp.mpf(row['2l_CA'][1])))
                worst['ZZ'] = max(worst['ZZ'], dF, dA)
        print(f'ZZ beta={beta}: C_F part {mp.nstr(pF, 12)}, C_A part {mp.nstr(pA, 12)}')
    print('largest complex |data - formula|:', {k: mp.nstr(v, 2) for k, v in worst.items()})
    assert worst['HH'] < 1e-8 and worst['ZZ'] < 1e-7
    json.dump({k: float(v) for k, v in worst.items()}, open(HERE / 'r0_closed_forms.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
