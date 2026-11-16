#!/usr/bin/env python3
"""Check two-loop boundary logarithms against retained ZZ amplitudes.

All coefficients come from unconstrained fits. This verifies the stored
amplitudes in the tested kinematics, not a general factorisation theorem.
"""
import json
from pathlib import Path
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
LOG=('++++','++--','----','--++')
ANG=('+-++','+---','-+++','-+--')

def pair(v):
    return mp.mpc(mp.mpf(v[0]),mp.mpf(v[1]))

def main():
    data=json.loads((ROOT/'paper2/zz/results_zz.json').read_text())['betas']
    with mp.workdps(60):
        ca,cf=mp.mpf(3),mp.mpf(4)/3
        delta,B=2*cf-ca,6*cf
        shift_factor=(ca-cf)/delta
        largest1=largest2=window=mp.mpf(0)
        count=0
        for beta,rec in sorted(data.items()):
            x=mp.mpf(rec['x'])
            shift=shift_factor*mp.log(x*(1-x))
            for h in LOG+ANG:
                row=rec['hel'][h]
                one=pair(row['c1'] if h in LOG else row['K'])
                if abs(one)<mp.mpf('1e-15'):
                    continue
                def coefficient(j,key):
                    return 2*(cf*pair(row['2l_CF'][key][j])+ca*pair(row['2l_CA'][key][j]))/one
                def boundary(key):
                    if h in LOG:
                        a3=coefficient(3,key)
                        a2=coefficient(2,key)+3*a3*shift
                        return -B-a2,delta+3*a3/4
                    b2=coefficient(2,key)
                    b1=coefficient(1,key)+2*b2*shift
                    return -B-(b1+b2)/2,delta+b2/4
                r1,r2=boundary('f')
                w1,w2=boundary('f_w1')
                largest1=max(largest1,abs(r1))
                largest2=max(largest2,abs(r2))
                window=max(window,abs(r1-w1),abs(r2-w2))
                count+=1
        assert count==24,count
        assert largest1<mp.mpf('1.5e-9'),largest1
        assert largest2<mp.mpf('2.3e-11'),largest2
        assert window>max(largest1,largest2),window
        print('PASS: {} ZZ configurations; max |r1|={}, |r2|={}, window shift={}'.format(
            count,mp.nstr(largest1,5),mp.nstr(largest2,5),mp.nstr(window,5)))

if __name__=='__main__':
    main()
