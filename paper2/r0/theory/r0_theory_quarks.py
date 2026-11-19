"""r0 for quark splittings from the regions, against the closed form found in the gg -> ZZ data.

The one-loop kernel at mu = mu_b, in units p1^+ = p2^- = mu = nu = 1, b = b0 (L_b = L_nu = 0), s = 1:
    K1/K0 = J1(x) + J2(x) + S + (self-energies of the four daughters),
J2(x) = J1(1 - x) with the helicity structures exchanged (charge conjugation), S the soft function of
section C3 of ../README.md.  The renormalisation of the coupling and of the external gluon field adds pure poles and no finite
part.  The Glauber region and the Coulomb phase of the soft function are not computed: the regions computed are real
and lack the pole 2 pi i Delta_C/eps of 2 I_K, so the comparison is for the real parts (Im r0 = 0 in the data).
Reads gqq_laurent_x*.json (assemble_gqq.py) and writes r0_theory_quarks.json.

    python3 r0_theory_quarks.py
"""
import json, os
from fractions import Fraction
import numpy as np
import mpmath as mp
from laurent_inputs import load_laurent

Z2 = float(mp.zeta(2))
PAIRS = [(1, 10), (3, 10), (1, 2)]          # x = num/den; J1 is needed at x and at 1 - x


def series(f, orders=(-2, -1, 0), r=0.2, n=32):
    th = 2*np.pi*(np.arange(n) + 0.5)/n
    e = r*np.exp(1j*th)
    v = np.array([complex(f(mp.mpc(ee))) for ee in e])
    return {j: np.mean(v*e**(-j)) for j in orders}


def ftf(ep):
    return mp.exp(-mp.euler*ep)*mp.gamma(1 - 2*ep)/(mp.gamma(1 + ep)*mp.gamma(1 - ep))


def self_energies(xx):
    """quark and antiquark of J1(xx), b space at L_b = 0: coefficient of C_F"""
    return series(lambda ep: -ftf(ep)*mp.gamma(ep)*mp.gamma(1 - ep)**2/mp.gamma(1 - 2*ep)*(1 - ep)/(1 - 2*ep)
                  * ((1 - xx)**ep + xx**ep))


SOFT_ETA = series(lambda ep: 2*mp.exp(-mp.euler*ep)*mp.gamma(-ep)*2)          # x Delta_C, coefficient of 1/eta
SOFT_0 = series(lambda ep: 2*mp.exp(-mp.euler*ep)*mp.gamma(-ep)*(-mp.euler - mp.digamma(-ep)))
DC = {'CF': 2, 'CA': -1}                                                       # Delta_C = 2 C_F - C_A


def load(tag, suffix=''):
    d = load_laurent(f'gqq_laurent_x{tag}{suffix}.json', 'gqq', expected_x=tag.split('_', 1)[0].replace('o', '/'))
    tot = {st: {c: {} for c in DC} for st in ('ek', 'ke')}
    for name, dd in d['diagrams'].items():
        pref = {k: float(Fraction(v)) for k, v in dd['pref'].items()}
        for st in ('ek', 'ke'):
            for key, (re, im) in dd['laurent'][st].items():
                for c in DC:
                    tot[st][c][key] = tot[st][c].get(key, 0) + pref[c]*complex(re, im)
    return tot


def r0_closed(xx):
    lx, ly = np.log(xx), np.log(1 - xx)
    return {'CF': 12 - 6*Z2 + 2*(lx - ly)**2, 'CA': 2 + Z2 - 2*(lx**2 + ly**2)}


def kernel(num, den, suffix=''):
    xx = num/den
    J1 = load(f'{num}o{den}', suffix)
    J2s = load(f'{den - num}o{den}', suffix)
    se_x, se_y = self_energies(xx), self_energies(1 - xx)
    se = {k: se_x[k] + se_y[k] for k in (-2, -1, 0)}
    out = {}
    flip = {'ek': 'ke', 'ke': 'ek'}
    # the four helicity configurations: structure s1 of J1 and s2 of J2, where J2 in structure s2 is J1(1 - x) in
    # the other structure (charge conjugation exchanges quark and antiquark, so that the parton carrying the
    # gluon helicity changes side)
    for s1, s2 in (('ek', 'ek'), ('ek', 'ke'), ('ke', 'ek'), ('ke', 'ke')):
        st1, st2 = s1, flip[s2]
        res = {}
        for c, dc in DC.items():
            fin = {j: J1[st1][c].get(f'0,{j}', 0) + J2s[st2][c].get(f'0,{j}', 0) + dc*SOFT_0[j]
                   + (se[j] if c == 'CF' else 0) for j in (-2, -1, 0)}
            eta = {j: J1[st1][c].get(f'-1,{j}', 0) + J2s[st2][c].get(f'-1,{j}', 0) + dc*SOFT_ETA[j] for j in (-2, -1, 0)}
            res[c] = {'eps^-2': fin[-2].real, 'eps^-1': fin[-1].real, 'eps^0': fin[0].real,
                      'max |imag|': max(abs(fin[j].imag) for j in fin),
                      'max |1/eta coefficient|': max(abs(eta[j]) for j in eta)}
        out[(s1, s2)] = res
    # J1 alone, at x and at 1 - x: the finite part of the C_F and C_A coefficients per helicity structure
    j1 = {st: {c: J1[st][c].get('0,0', 0).real for c in DC} for st in ('ek', 'ke')}
    j1y = {st: {c: J2s[st][c].get('0,0', 0).real for c in DC} for st in ('ek', 'ke')}
    return out, (j1, j1y)


def main():
    results = {'x': {}}
    print(f'{"x":>5} {"col":>3} {"eps^-2":>13} {"expected":>9} {"eps^-1":>14} {"expected":>14} {"r0 (regions)":>14} '
          f'{"r0 (data)":>14} {"diff":>9}')
    for num, den in PAIRS:
        xx = num/den
        out, (j1, j1y) = kernel(num, den)
        data = r0_closed(xx)
        lx = np.log(xx*(1 - xx))
        expect = {'CF': (4.0, -4*lx + 6), 'CA': (-2.0, 0.0)}
        rec = {}
        for c in DC:
            th = out[('ek', 'ek')][c]
            spread = max(out[h][c]['eps^0'] for h in out) - min(out[h][c]['eps^0'] for h in out)
            print(f'{xx:5.2f} {c:>3} {th["eps^-2"]:13.9f} {expect[c][0]:9.4f} {th["eps^-1"]:14.9f} {expect[c][1]:14.9f}'
                  f' {th["eps^0"]:14.9f} {data[c]:14.9f} {th["eps^0"] - data[c]:9.1e}   helicity spread {spread:.1e}')
            rec[c] = {'theory, J1 structure ek, J2 structure ek': th,
                      'max spread of eps^0 over the four helicity configurations': spread,
                      'expected poles': {'eps^-2': expect[c][0], 'eps^-1': expect[c][1]},
                      'r0 closed form (data)': data[c], 'difference': th['eps^0'] - data[c],
                      'J1 finite part, ek': j1['ek'][c], 'J1 finite part, ke': j1['ke'][c],
                      'J1 at 1 - x, finite part, ek': j1y['ek'][c], 'J1 at 1 - x, finite part, ke': j1y['ke'][c]}
        results['x'][f'{num}/{den}'] = rec
    # precision: (a) higher quadrature orders, (b) 16 instead of 32 points on the eps circle (aliasing), at x = 1/2
    base, _ = kernel(1, 2)
    checks = {}
    for label, suffix in (('quadrature orders 96, 72, 40 instead of 64, 48, 28', '_nq96-72-40'),
                          ('16 instead of 32 points on the eps circle', '_neps16')):
        if os.path.exists(f'gqq_laurent_x1o2{suffix}.json'):
            out, _ = kernel(1, 2, suffix)
            checks[label] = {c: out[('ek', 'ek')][c]['eps^0'] - base[('ek', 'ek')][c]['eps^0'] for c in DC}
            print(f'x = 1/2, {label}: change of r0', {c: f'{v:+.1e}' for c, v in checks[label].items()})
    results['precision checks at x = 1/2 (change of r0)'] = checks
    json.dump(results, open('r0_theory_quarks.json', 'w'), indent=1)
    print('wrote r0_theory_quarks.json')


if __name__ == '__main__':
    main()
