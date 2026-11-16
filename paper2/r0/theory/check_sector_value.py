"""Sector-decomposed value against direct quadrature (tanh-sinh in every variable, vectorised), for one
diagram at one real eps where the parametric integral converges (no eta).

    python3 check_sector_value.py 3/10 Vqq ek -3/10
"""
import sys, warnings
warnings.filterwarnings('ignore')
import numpy as np
import sympy as sp
import scipy.special as ss
from lc_engine import eps, D
from gqq_diagrams import z
import assemble_gqq as A
from numeric_checks import close,interior_values

args = [a for a in sys.argv[1:] if a != '--allow-legacy-input']
name, st = args[1], args[2]
epsv = sp.Rational(args[3])
H = float(sp.Rational(args[4])) if len(args) > 4 else 1/16
A.setup(sp.Rational(args[0]), allow_legacy='--allow-legacy-input' in sys.argv)


def ts_nodes(h=H, tmax=3.6):
    """tanh-sinh nodes on (0, 1): x, 1 - x and weights"""
    t = np.arange(-tmax, tmax + h/2, h)
    u = 0.5*np.pi*np.sinh(t)
    x = 1/(1 + np.exp(-2*u))              # (1 + tanh u)/2
    xc = 1/(1 + np.exp(2*u))              # 1 - x, without cancellation
    wgt = h*0.5*np.pi*np.cosh(t)/(2*np.cosh(u)**2)
    interior=(x>0)&(x<1)&(xc>0)&(xc<1)
    return x[interior], xc[interior], wgt[interior]


xs, xcs, ws = ts_nodes()


tot_sec, tot_dir = 0, 0
checked=0
for iv, pcs in A.data['results'][name]:
    for prm, o in pcs:
        if o[st] == 0:
            continue
        secs = A.piece_sectors(o[st], iv, prm, False)
        vs = sum(S.value(complex(epsv), 0.0, n={1: 64, 2: 48, 3: 28}[len(S.vs)]) for S in secs)
        e = o[st].subs(D, 2 - 2*eps).subs(eps, epsv)
        za, zb = sp.Rational(iv[0]), sp.Rational(iv[1])
        zz = float(za) + float(zb - za)*xs
        jz = float(zb - za)
        if len(prm) == 0:
            vd = np.sum(ws*jz*interior_values(e,[z],[zz+0j]))
        elif len(prm) == 1:
            Z, Aa = np.meshgrid(zz, xs, indexing='ij')
            W = np.outer(ws*jz, ws)
            vd = np.sum(W*interior_values(e,[z]+list(prm),[Z+0j,Aa+0j]))
        else:
            a1, a2 = prm
            Z, S1, T1 = np.meshgrid(zz, xs, xs, indexing='ij')
            _, S1c, _ = np.meshgrid(zz, xcs, xs, indexing='ij')
            W = np.einsum('i,j,k->ijk', ws*jz, ws, ws)
            vd = np.sum(W*S1c*interior_values(e,[z,a1,a2],[Z+0j,S1+0j,S1c*T1+0j]))
        print(f'  interval {iv} params {len(prm)}: sectors {vs.real:+.10f} (imag {vs.imag:+.1e})   direct {vd.real:+.10f} (imag {vd.imag:+.1e})', flush=True)
        close(vs,vd,label=f'{iv} sector integral')
        tot_sec += vs; tot_dir += vd;checked+=1
if not checked:raise ValueError('No nonzero integrands were compared')
close(tot_sec,tot_dir,label='total sector integral')
print(f'{name} {st} eps = {epsv}: sectors {tot_sec.real:+.10f}   direct {tot_dir.real:+.10f}')
