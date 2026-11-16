#!/usr/bin/env python3
"""Stage B2, step 1c: the cross-pair contribution to the one-loop DPS kernel in momentum space,
  K_cross(k, q) = (1/4) int dl^+ d^2l_perp  [ int dl^- N R^- R^+ / D7 ]   (model_lminus.py, closed form),
evaluated numerically (double precision; the closed form was validated against mpmath), as a
function of lambda with k = lambda*khat, q = lambda*qhat fixed in direction and ratio, to read off
the logarithmic structure in lambda.  Units s = 1.
    python3 model_exact.py [lambda ...]
"""
import math, cmath, sys
import numpy as np
from scipy import integrate

PI = math.pi
I = 1j

def lminus_closed(lp, ltx, lty, kx, ky, qx, qy, x):
    """closed form of the l^- integral (model_lminus.py, double precision)"""
    k2 = kx*kx + ky*ky
    kl2 = (kx+ltx)**2 + (ky+lty)**2
    a = (kx+ltx-qx)**2 + (ky+lty-qy)**2
    a6 = (kx-qx)**2 + (ky-qy)**2
    N = kx*(kx+ltx-qx) + ky*(ky+lty-qy)
    lt2 = ltx*ltx + lty*lty
    kp = x - lp
    rminus = []
    km3 = -k2/(1-kp)
    D1_3 = kp*km3 - k2
    rminus.append((2*PI*I/(-(1-kp))/D1_3/x, (kl2 - x*km3)/x))
    if kp < 0:
        km1 = k2/kp
        D3_1 = -(1-kp)*km1 - k2
        rminus.append((2*PI*I/kp/D3_1/x, (kl2 - x*km1)/x))
    B = x*lp + a6
    z7 = lt2/lp
    s7 = 1 if lp > 0 else -1
    tot = 0
    for pref, z2 in rminus:
        z6 = (1-x) + x*a/B
        z6p = -x + x*a/B
        for cut, pre, numf, poles in ((1-x, 2*PI*I*pref*N/(a*B*lp), lambda lm: (1-x) - lm, [(z6, -1), (z7, s7), (z2, +1)]),
                                      (-x, -2*PI*I*pref*N/(a*(-B)*lp), lambda lm: lm + x, [(z6p, +1), (z7, s7), (z2, +1)])):
            acc = 0
            for j, (zj, sj) in enumerate(poles):
                resid = numf(zj)
                for i, (zi, si) in enumerate(poles):
                    if i != j:
                        resid /= (zj - zi)
                acc += resid*(math.log(abs(cut - zj)) - I*PI*sj*(1 if cut > zj else 0))
            tot += pre*acc
    return tot

def lplus_integral(ltx, lty, kx, ky, qx, qy, x):
    """int_{-(1-x)}^{oo} dl^+ of the l^- closed form; breakpoints at the branch points"""
    a6 = (kx-qx)**2 + (ky-qy)**2
    f_re = lambda lp: lminus_closed(lp, ltx, lty, kx, ky, qx, qy, x).real
    f_im = lambda lp: lminus_closed(lp, ltx, lty, kx, ky, qx, qy, x).imag
    lo = -(1-x)
    brk = sorted({lo, -a6/x, 0.0, x, 1.0, 5.0})
    pieces = list(zip(brk[:-1], brk[1:])) + [(brk[-1], np.inf)]
    tot = 0
    for aa, bb in pieces:
        if bb == np.inf:
            r = integrate.quad(f_re, aa, bb, limit=200)[0]; im = integrate.quad(f_im, aa, bb, limit=200)[0]
        else:
            r = integrate.quad(f_re, aa, bb, limit=200, points=None)[0]; im = integrate.quad(f_im, aa, bb, limit=200)[0]
        tot += r + I*im
    return tot

def kernel_cross(lam, khat=(0.6, 0.3), qhat=(1.0, 0.0), x=0.75, rmax_factor=400.0):
    """(1/4) int d^2 l_perp of the l^+ integral, in polar coordinates about t = l_perp + k - q = 0"""
    kx, ky = lam*khat[0], lam*khat[1]; qx, qy = lam*qhat[0], lam*qhat[1]
    cx, cy = qx - kx, qy - ky            # l_perp = c + t
    def integrand(r, phi):
        ltx, lty = cx + r*math.cos(phi), cy + r*math.sin(phi)
        return lplus_integral(ltx, lty, kx, ky, qx, qy, x)*r
    def radial(phi, part):
        f = lambda r: getattr(integrand(r, phi), part)
        brk = [0.0, 0.5*lam, lam, 2*lam, 5*lam, 20*lam, 100*lam, rmax_factor*lam]
        tot = 0
        for aa, bb in zip(brk[:-1], brk[1:]):
            tot += integrate.quad(f, aa, bb, limit=100)[0]
        tot += integrate.quad(f, brk[-1], np.inf, limit=100)[0]
        return tot
    re = integrate.quad(lambda phi: radial(phi, 'real'), 0, 2*PI, limit=50)[0]
    im = integrate.quad(lambda phi: radial(phi, 'imag'), 0, 2*PI, limit=50)[0]
    return 0.25*(re + I*im)

if __name__ == '__main__':
    import time
    lams = [float(a) for a in sys.argv[1:]] or [0.03]
    for lam in lams:
        t0 = time.time()
        val = kernel_cross(lam)
        print(f'lambda = {lam}: K_cross = {val!r}   ({time.time()-t0:.0f} s)', flush=True)
