"""Combine the one-loop collinear functions J1(x) and J2(y = 1 - x), the self-energies, the coupling counterterm
and the soft function into the one-loop kernel for quark splittings at mu = mu_b, and compare with the data.

    python3 combine_gqq.py 3/10        (needs gqq_laurent_x3o10.json and gqq_laurent_x7o10.json)

J2(y) is taken as J1 at y (charge conjugation maps the n-bar sector with an antiquark at x1 onto J1); both
helicity structures are printed, so that their equality (helicity independence) is visible.
Units p1^+ = p2^- = mu = nu = 1, b = b0: L_b = L_nu = 0, s = 1, m^2 = x(1-x).
"""
import sys, json
import numpy as np
import mpmath as mp
from laurent_inputs import load_laurent

x = mp.mpf(mp.fraction(*map(int, sys.argv[1].split('/')))) if len(sys.argv) > 1 else mp.mpf(3)/10
tagx = sys.argv[1].replace('/', 'o') if len(sys.argv) > 1 else '3o10'
num, den = map(int, (sys.argv[1] if len(sys.argv) > 1 else '3/10').split('/'))
tagy = f'{den - num}o{den}'
z2 = float(mp.zeta(2))


def series(f, orders=(-2, -1, 0), r=0.2, n=32):
    th = 2*np.pi*(np.arange(n) + 0.5)/n
    e = r*np.exp(1j*th)
    v = np.array([complex(f(mp.mpc(ee))) for ee in e])
    return {j: np.mean(v*e**(-j)) for j in orders}


def ftf(ep):
    return mp.exp(-mp.euler*ep)*mp.gamma(1 - 2*ep)/(mp.gamma(1 + ep)*mp.gamma(1 - ep))


def analytic(xx):
    """self-energies of the quark and the antiquark of J1(xx), b-space at L_b = 0: C_F coefficient"""
    se = lambda ep: -ftf(ep)*mp.gamma(ep)*mp.gamma(1 - ep)**2/mp.gamma(1 - 2*ep)*(1 - ep)/(1 - 2*ep)*((1 - xx)**ep + xx**ep)
    return series(se)


def load(tag):
    d = load_laurent(f'gqq_laurent_x{tag}.json', 'gqq', expected_x=tag.split('_', 1)[0].replace('o', '/'))
    tot = {st: {col: {} for col in ('CF', 'CA')} for st in ('ek', 'ke')}
    per = {}
    for name, dd in d['diagrams'].items():
        pref = {k: float(eval(v)) for k, v in dd['pref'].items()}
        for st in ('ek', 'ke'):
            for key, (re, im) in dd['laurent'][st].items():
                for col in ('CF', 'CA'):
                    tot[st][col][key] = tot[st][col].get(key, 0) + pref[col]*complex(re, im)
                per.setdefault(name, {}).setdefault(st, {})[key] = complex(re, im)
    return tot, per


def main():
    J1, per1 = load(tagx)
    J2, per2 = load(tagy)
    se1, se2 = analytic(x), analytic(1 - x)
    # soft function: 2 Delta_C e^{-g eps} Gamma(-eps) [2/eta - gamma - psi(-eps)], Delta_C = 2 C_F - C_A
    soft_eta = series(lambda ep: 2*mp.exp(-mp.euler*ep)*mp.gamma(-ep)*2)                 # coefficient of 1/eta
    soft_0 = series(lambda ep: 2*mp.exp(-mp.euler*ep)*mp.gamma(-ep)*(-mp.euler - mp.digamma(-ep)))
    Lx = float(mp.log(x*(1 - x)))
    for st in ('ek', 'ke'):
        print(f'--- helicity structure {st}')
        for col, dc in (('CF', 2), ('CA', -1)):
            tot = {}
            for j in (-2, -1, 0):
                key0 = f'0,{j}'
                v = J1[st][col].get(key0, 0) + J2[st][col].get(key0, 0) + dc*soft_0[j]
                if col == 'CF':
                    v += se1[j] + se2[j]
                tot[j] = v
            eta_pole = {j: J1[st][col].get(f'-1,{j}', 0) + J2[st][col].get(f'-1,{j}', 0) + dc*soft_eta[j] for j in (-2, -1, 0)}
            print(f'  {col}: 1/eta coefficients (must vanish):', {j: f'{eta_pole[j].real:.2e}' for j in eta_pole})
            print(f'  {col}: eps^-2 {tot[-2].real:.10f}  eps^-1 {tot[-1].real:.10f}  eps^0 {tot[0].real:.10f}'
                  f'   (imag. parts {abs(tot[-2].imag):.1e} {abs(tot[-1].imag):.1e} {abs(tot[0].imag):.1e})')
            if col == 'CF':
                print(f'      expected poles: eps^-2 4, eps^-1 {-4*Lx + 6:.10f}')
                lx, ly = float(mp.log(x)), float(mp.log(1 - x))
                print(f'      data r0: {12 - 6*z2 + 2*(lx - ly)**2:.10f}')
            else:
                print('      expected poles: eps^-2 -2, eps^-1 0')
                lx, ly = float(mp.log(x)), float(mp.log(1 - x))
                print(f'      data r0: {2 + z2 - 2*(lx**2 + ly**2):.10f}')


if __name__ == '__main__':
    main()
