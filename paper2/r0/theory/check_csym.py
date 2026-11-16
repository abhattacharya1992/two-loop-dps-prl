"""Charge conjugation at the level of the integrands: J1 at x with the gluon fraction z, structure ek, must equal
J1 at 1 - x with -z, structure ke, diagram by diagram for the C-even diagrams (Vqq, Vggq).  The Feynman-parameter
integral at fixed z and fixed eps is done with mpmath (adaptive).

    python3 check_csym.py Vqq 1/10 -0.3
"""
import sys, pickle
import sympy as sp
import mpmath as mp
import numpy as np
from scipy import integrate
import scipy.special as ss
from lc_engine import eps, D
from gqq_diagrams import z

name = sys.argv[1] if len(sys.argv) > 1 else 'Vqq'
zabs = sp.Rational(sys.argv[2]) if len(sys.argv) > 2 else sp.Rational(1, 10)
epsv = sp.Rational(sys.argv[3]) if len(sys.argv) > 3 else sp.Rational(-3, 10)
mp.mp.dps = 20


def density(tag, zval, st):
    d = pickle.load(open(f'gqq_integrands_x{tag}.pkl', 'rb'))
    tot = 0j
    for iv, pcs in d['results'][name]:
        if not (sp.Rational(iv[0]) < zval < sp.Rational(iv[1])):
            continue
        for prm, o in pcs:
            e = o[st]
            if e == 0:
                continue
            e = e.subs(D, 2 - 2*eps).subs({z: zval, eps: epsv})
            if len(prm) == 0:
                tot += complex(sp.N(e, 20))
            elif len(prm) == 1:
                f = sp.lambdify(prm, e, [{'gamma': ss.gamma}, 'numpy'])
                re_ = integrate.quad(lambda a: complex(f(a)).real, 0, 1, limit=400, epsabs=1e-11, epsrel=1e-10)[0]
                im_ = integrate.quad(lambda a: complex(f(a)).imag, 0, 1, limit=400, epsabs=1e-11, epsrel=1e-10)[0]
                tot += complex(re_, im_)
            else:
                a1, a2 = prm
                f = sp.lambdify([a1, a2], e, [{'gamma': ss.gamma}, 'numpy'])
                g = lambda t, s_: (1 - s_)*complex(f(s_, (1 - s_)*t)).real
                tot += integrate.dblquad(g, 0, 1, 0, 1, epsabs=1e-10, epsrel=1e-9)[0]
    return tot


for zz in (zabs, -zabs):
    a = density('3o10', zz, 'ek')
    b = density('7o10', -zz, 'ke')
    c = density('3o10', zz, 'ke')
    dd = density('7o10', -zz, 'ek')
    print(f'z = {float(zz):+.3f}: x=3/10 ek {a}   x=7/10 (-z) ke {b}   |   '
          f'x=3/10 ke {c}   x=7/10 (-z) ek {dd}')
