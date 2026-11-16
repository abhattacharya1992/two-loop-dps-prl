#!/usr/bin/env python3
"""Gate C, step 1: the constant r0 of the one-loop kernel as the two-loop gg -> HH remainders give it.

At mu_R = m_H the DPS terms of the two-loop remainders (Catani's subtraction, Paper 2 eq. (4.1)) are,
with H1 H2 K factorisation in minimal subtraction and the Catani reading of section 5,
    c11 / c1 = 2h + Dc l^2 + B l - delta0 + r0(nu=0),
    c20 / K1 = (B - 2Dc) + 2h + Dc l^2 + B l - delta0 + r0(nu=2),
where c11 is the coefficient of G in A1^(2), c20 the forward limit of A2^(2), and
    l      = ln(mu^2/(Q*^2 e^{-i pi})) = ln(s/m_H^2) + i pi           (Q*^2 = m_H^4/s),
    h      = C_A(-L_H^2 + zeta2), L_H = ln(m_H^2/mu^2) - i pi           (MS-subtracted gg -> H, HEFT),
    delta0 = [2 I^(1)]_{eps^0} = -C_A L_s^2 + C_A zeta2 - beta0 L_s,  L_s = ln(mu^2/s) + i pi,
and (B - 2Dc) is the constant of the nu = 2 transform of -Dc X^2 + B X (Paper 2, below eq. (5.5)).
Dc = C_A = N, B = 2 beta0, beta0 = 11N/3 - 2 n_f/3.  The colour terms of the data are separated as in
Paper 2 eq. (3.1); c1 = -2 pi N F_S and K1 = pi N W are proportional to N, so the N^2 term of the data
gives the N part of r0 and the N n_f term its n_f part.

    python3 r0_from_data.py      (reads the point files of task B; a few seconds)
"""
import json, sys, pathlib
from fractions import Fraction
import mpmath as mp

HERE = pathlib.Path(__file__).resolve().parent
TASKS = HERE.parents[1]
sys.path.insert(0, str(TASKS / 'hh' / 'hh_two_loop'))
from scan_tools import load, fit_model, weights, beta_number, COLOURS   # noqa: E402

mp.mp.dps = 50
z2 = mp.zeta(2)


def known_parts(beta):
    """N- and n_f-coefficients of 2h + Dc l^2 + B l - delta0 at mu = m_H (relative to the one-loop term)."""
    b = beta_number(beta)
    v = (1-b)/(1+b)
    s = (1+v)**2/v                      # s/m_H^2
    Ls = mp.mpf(s)
    l = mp.log(Ls) + 1j*mp.pi           # ln(mu^2/(Q*^2 e^{-i pi}))
    LH = -1j*mp.pi                      # ln(m_H^2/mu^2) - i pi
    Lsub = -mp.log(Ls) + 1j*mp.pi       # ln(mu^2/s) + i pi
    # coefficients of N (C_A and the N part of beta0) and of n_f (the n_f part of beta0)
    def total(CA, b0):
        h = CA*(-LH**2 + z2)
        delta0 = -CA*Lsub**2 + CA*z2 - b0*Lsub
        return 2*h + CA*l**2 + 2*b0*l - delta0
    N_part = total(1, mp.mpf(11)/3)
    nf_part = total(0, -mp.mpf(2)/3)
    nu2_N = 2*mp.mpf(11)/3 - 2          # (B - 2Dc) for the N part
    nu2_nf = 2*(-mp.mpf(2)/3)
    return {'nu0': {'N': N_part, 'nf': nf_part}, 'nu2': {'N': N_part + nu2_N, 'nf': nf_part + nu2_nf},
            'l': l, 's_over_mH2': s}


def extract(points):
    reg, _ = load(points)
    betas = sorted({b for (_, b) in reg}, key=Fraction)
    out = {}
    for beta in betas:
        if len(reg.get((1, beta), [])) < 8 or len(reg.get((2, beta), [])) < 8:
            continue                                      # scan of this velocity not complete
        w, FS, Kt = weights(beta)
        c1_over_N = -2*mp.pi*FS
        K1_over_N = mp.pi*Kt
        kp = known_parts(beta)
        res = {'x': str(w), 'F_S': str(FS), 'W': str(Kt)}
        for comp, gauge, deg, power, norm in (('nu0', 1, 3, 1, c1_over_N), ('nu2', 2, 2, 0, K1_over_N)):
            recs = reg[(gauge, beta)]
            vals = {}
            for col in COLOURS:
                for part in (0, 1, 2):                       # L^0, L^1, L^2 of A^(2)
                    f0 = fit_model(recs, col, part, deg, 2, 0)[power]
                    f1 = fit_model(recs, col, part, deg, 2, 1)[power]
                    vals[(col, part)] = (f0/norm, abs(f1 - f0)/abs(norm))
            r0N = vals[('N2', 0)][0] - kp[comp]['N']
            r0nf = vals[('N_Nf', 0)][0] - kp[comp]['nf']
            res[comp] = {
                'data_over_oneloop': {col: [str(mp.re(vals[(col, 0)][0])), str(mp.im(vals[(col, 0)][0])),
                                            str(vals[(col, 0)][1])] for col in COLOURS},
                'L1_over_oneloop': {col: [str(mp.re(vals[(col, 1)][0])), str(mp.im(vals[(col, 1)][0]))] for col in COLOURS},
                'known_N': [str(mp.re(kp[comp]['N'])), str(mp.im(kp[comp]['N']))],
                'known_nf': [str(mp.re(kp[comp]['nf'])), str(mp.im(kp[comp]['nf']))],
                'r0_N': [str(mp.re(r0N)), str(mp.im(r0N)), str(vals[('N2', 0)][1])],
                'r0_nf': [str(mp.re(r0nf)), str(mp.im(r0nf)), str(vals[('N_Nf', 0)][1])],
            }
        out[str(beta)] = res
    return out


def main():
    out = extract(TASKS / 'hh' / 'hh_two_loop' / 'points')
    for beta, r in out.items():
        print(f"beta = {beta}, x = {mp.nstr(mp.mpf(r['x']), 6)}")
        for comp in ('nu0', 'nu2'):
            c = r[comp]
            print(f"  {comp}: r0 N-part  = {mp.nstr(mp.mpf(c['r0_N'][0]), 12)} + {mp.nstr(mp.mpf(c['r0_N'][1]), 3)} i"
                  f"   (shift {mp.nstr(mp.mpf(c['r0_N'][2]), 2)})")
            print(f"       r0 nf-part = {mp.nstr(mp.mpf(c['r0_nf'][0]), 12)} + {mp.nstr(mp.mpf(c['r0_nf'][1]), 3)} i"
                  f"   (shift {mp.nstr(mp.mpf(c['r0_nf'][2]), 2)})")
            print('       n_f/N, n_f^2 terms of the data / one-loop:',
                  [mp.nstr(abs(mp.mpc(mp.mpf(c['data_over_oneloop'][k][0]), mp.mpf(c['data_over_oneloop'][k][1]))), 2)
                   for k in ('Nf_over_N', 'Nf2')])
    json.dump(out, open(HERE / 'r0_from_data.json', 'w'), indent=1)
    print('written r0_from_data.json')


if __name__ == '__main__':
    main()
