#!/usr/bin/env python3
"""Gate C, step 1 for quark splittings: r0 from the two-loop gg -> ZZ amplitudes (task J fits).

The fits of task J (results_zz.json) are in VVamp units (alpha_s/(2 pi) per loop) and in the qT scheme at
mu^2 = s.  In the a_s = alpha_s/(4 pi) expansion of Paper 2 the two-loop/one-loop ratios are twice the
VVamp ratios, and ggvvamp's conversion to Catani's scheme, A2 -> A2 + (-pi^2 C_A/2 + i pi b0) A1 with
b0 = 11 C_A/6 - n_f/3 (ggvvamp.cpp, convert_A_to_Catani), becomes c^Catani = c^qT + (-pi^2 C_A + i pi beta0) c1.
With delta0 = [2 I^(1)]_{eps^0} = C_A pi^2 + C_A zeta2 - i pi beta0 at mu^2 = s, minimal subtraction is
    c^MS = c^qT + C_A zeta2 c1.
In the variable G of the data (gluon scale Q_g^2 = m^4/s) the bracket of Paper 2 eq. (5.3) is
    -Dc X_g^2 + (B - 2 Dc kappa) X_g + Dc l^2 - 2 Dc kappa l + B l,
with Dc = 2C_F - C_A, B = 6 C_F, Dc kappa = 2 (C_A - C_F) L_x, L_x = ln[x(1-x)], l = ln(mu^2/(Q_g^2 e^{-i pi})) = -2 L_x + i pi.
Its constant and the hard functions h = C_F(-L^2 + 3L - 8 + zeta2), L = ln(m^2/mu^2) - i pi = L_x - i pi, give
    c11/c1 = 2h + Dc l^2 - 2 Dc kappa l + B l - C_A zeta2 + r0(nu=0),
    c20/K  = (B - 2 Dc - 2 Dc kappa) + 2h + Dc l^2 - 2 Dc kappa l + B l - C_A zeta2 + r0(nu=2),
for each colour factor (C_F, C_A) separately.
    python3 r0_from_zz.py
"""
import json, pathlib
import mpmath as mp

HERE = pathlib.Path(__file__).resolve().parent
TASKS = HERE.parents[1]
mp.mp.dps = 40
z2 = mp.zeta(2)


def C(v):
    return mp.mpc(mp.mpf(v[0]), mp.mpf(v[1]))


def known(x, CF, CA):
    Lx = mp.log(x*(1-x))
    l = -2*Lx + 1j*mp.pi
    Dc = 2*CF - CA
    Dk = 2*(CA - CF)*Lx
    B = 6*CF
    L = Lx - 1j*mp.pi
    h = CF*(-L**2 + 3*L - 8 + z2)
    nu0 = 2*h + Dc*l**2 - 2*Dk*l + B*l - CA*z2
    nu2 = nu0 + (B - 2*Dc - 2*Dk)
    return nu0, nu2


def main():
    d = json.load(open(TASKS / 'zz' / 'results_zz.json'))['betas']
    out = {}
    for beta, r in d.items():
        x = mp.mpf(r['x'])
        res = {'x': r['x'], 'nu0': {}, 'nu2': {}}
        cmax = max(abs(C(e['c1'])) for e in r['hel'].values())
        for h, e in r['hel'].items():
            if '2l_CF' not in e:
                continue
            c1, K = C(e['c1']), C(e['K'])
            if e['dJz'] == 0 and abs(c1) > 1e-3*cmax:
                comp, norm, j = 'nu0', c1, 1
            elif abs(e['dJz']) == 2 and abs(K) > 1e-6:
                comp, norm, j = 'nu2', K, 0
            else:
                continue
            row = {}
            for part, (CF, CA) in (('2l_CF', (1, 0)), ('2l_CA', (0, 1))):
                f = C(e[part]['f'][j]); fw = C(e[part]['f_w1'][j])
                kn = known(x, CF, CA)[0 if comp == 'nu0' else 1]
                r0 = 2*f/norm - kn
                row[part] = [str(mp.re(r0)), str(mp.im(r0)), str(2*abs(f - fw)/abs(norm))]
            res[comp][h] = row
        out[beta] = res
        print(f'beta = {beta}, x = {mp.nstr(x, 6)}')
        for comp in ('nu0', 'nu2'):
            for h, row in res[comp].items():
                print(f'  {comp} {h}: ' + '   '.join(f"{p[3:]}: {mp.nstr(mp.mpf(v[0]), 12)} + {mp.nstr(mp.mpf(v[1]), 3)} i (shift {mp.nstr(mp.mpf(v[2]), 2)})"
                                                    for p, v in row.items()))
    json.dump(out, open(HERE / 'r0_from_zz.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
