#!/usr/bin/env python3
# Paper Sec. 2 and Fig. 2: the exact one-loop gg -> AA amplitudes of P1 along lines of fixed
# beta, from 90 degrees to the forward limit, compared with the DPS terms of Sec. 2.3.
#
# The amplitudes are the one-loop finite remainders of arXiv:1909.08993 (ancillary files
# Finite_1Loop_P1 and Finite_1Loop_P2), evaluated with the module ggaa.oneloop of the P1
# package: by default common/aa_source of the master folder or a legacy Paper I source tree (paper2/
# next to src/ and inputs/), otherwise the one whose top directory is given as the first argument. Its inputs
# must have been downloaded (`make inputs` there). With m_A = mu_R = 1 the files need no scale
# correction. For every point the script stores
#   M_1^(1)/C_1 - G     with C_1 = 2 pi N_c m^2 F_A(x) and G of eq. (2.8): tends to a constant,
#   M_2^(1)/K           with K = pi N_c m^2 W(x): tends to 1,
# at sin^2(theta) = 10^(-j/8), j = 0, 1, ..., 88 (p_T^2 = p_max^2 sin^2(theta), p_max^2 =
# beta^2/(1 - beta^2) m^2), for beta = 1/4, 1/2, 2/3, with 60 working digits.
#
#   python3 aa_full_range.py [P1_ROOT] [OUT_JSON]       (about a minute)
import json
import pathlib
import sys
from fractions import Fraction

import mpmath as mp

HERE = pathlib.Path(__file__).resolve().parent
def default_root():
    # common/aa_source of the master folder, or a legacy Paper I source tree
    for cand in (HERE.parents[2] / 'common' / 'aa_source', HERE.parents[2] / 'paper1', HERE.parents[2]):
        if (cand / 'src' / 'python' / 'ggaa').is_dir():
            return cand
    return HERE.parents[2] / 'common' / 'aa_source'


ROOT = pathlib.Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else default_root()
OUT = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else HERE/'results_aa_curve.json'
ANC = ROOT/'inputs'/'arXiv-1909.08993'/'anc'
if not (ROOT/'src'/'python'/'ggaa').is_dir():
    sys.exit(f'the P1 package was not found at {ROOT}; give its top directory as the first argument')
if not ANC.is_dir():
    sys.exit(f'{ANC} not found; download the inputs of the P1 package first (make inputs there)')
sys.path.insert(0, str(ROOT/'src'/'python'))
from ggaa.oneloop import OneLoop          # noqa: E402  (the P1 package)

mp.mp.dps = 60
NC = 3
BETAS = ('1/4', '1/2', '2/3')
JMAX = 88


def main():
    ol = OneLoop(ANC)
    out = {}
    for beta in BETAS:
        q = Fraction(beta)
        b = mp.mpf(q.numerator)/q.denominator
        x = (1+b)/2
        FA = 1/(x*(1-x))-2*x*(1-x)
        W = (x**4+(1-x)**4)/(x*(1-x))
        C1, K = 2*mp.pi*NC*FA, mp.pi*NC*W
        s = 4/(1-b**2)
        rows = []
        for j in range(JMAX+1):
            sin2 = mp.mpf(10)**(-mp.mpf(j)/8)
            c = mp.sqrt(1-sin2)                   # 0 <= theta <= 90 degrees
            pT2 = b**2/(1-b**2)*sin2
            G = mp.log(s*pT2)/2+mp.mpc(0, 1)*mp.pi/2
            r1 = ol(1, b, c, mA=1, muR=1)/C1-G
            r2 = ol(2, b, c, mA=1, muR=1)/K
            rows.append({'j': j, 'pT2': mp.nstr(pT2, 30), 'G': [mp.nstr(mp.re(G), 30), mp.nstr(mp.im(G), 30)],
                         'M1_over_C1_minus_G': [mp.nstr(mp.re(r1), 30), mp.nstr(mp.im(r1), 30)],
                         'M2_over_K': [mp.nstr(mp.re(r2), 30), mp.nstr(mp.im(r2), 30)]})
        out[beta] = {'C1': mp.nstr(C1, 30), 'K': mp.nstr(K, 30), 'rows': rows}
        first, last = rows[0], rows[-1]
        print(f"beta = {beta}: p_T^2/m^2 from {mp.nstr(mp.mpf(first['pT2']), 3)} to {mp.nstr(mp.mpf(last['pT2']), 3)};"
              f" M_1/C_1 - G from {first['M1_over_C1_minus_G'][0][:8]} to {last['M1_over_C1_minus_G'][0][:12]} (real part),"
              f" M_2/K from {first['M2_over_K'][0][:8]} to {last['M2_over_K'][0][:14]} (real part)")
    OUT.write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
