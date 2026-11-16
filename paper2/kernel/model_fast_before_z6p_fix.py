#!/usr/bin/env python3
"""Vectorised double-precision version of model_exact2.py (same formulas; validated against it by
model_ref_check.py).  F(l_perp) = int dl^+ F_lminus(l^+, l_perp; eps0), F_lminus = closed-form l^- integral.
Gauss-Legendre (64 nodes) on segments: log-spaced around l^+ = 0 and around B = 0, plus windows around
the pinch points (width 0.1 ... 1e4 eps0).  The transverse integral is in model_2d.py.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
PI = np.pi

def pinch_points(ltx, lty, kx, ky, qx, qy, x):
    """real l^+ in (-(1-x), oo) where an upper and a lower l^- pole coincide.  Sides: z6 above, z6p below,
    z2 (and z2' of the D1 branch) below, z7 below for l^+ > 0 and above for l^+ < 0.  Pinches:
    z6 = z2 (any l^+), z6 = z2' (l^+ > x), z6 = z7 (l^+ > 0), z7 = z2 (l^+ < 0), z6p = z7 (l^+ < 0)."""
    k2 = kx*kx + ky*ky; kl2 = (kx+ltx)**2 + (ky+lty)**2
    a = (kx+ltx-qx)**2 + (ky+lty-qy)**2; a6 = (kx-qx)**2 + (ky-qy)**2; lt2 = ltx*ltx + lty*lty
    lo = -(1-x)
    def real_roots(poly, cond):
        return [v.real for v in np.roots(poly) if abs(v.imag) < 1e-10*max(1.0, abs(v.real)) and cond(v.real)]
    pts = []
    # z6 = z7: [(1-x)(x l + a6) + x a] l = lt2 (x l + a6)
    pts += real_roots([(1-x)*x, (1-x)*a6 + x*a - lt2*x, -lt2*a6], lambda v: v > 0)
    # z6 = z2, z2 = kl2/x + k2/(1-x+l):  [(1-x)(x l + a6) + x a](1-x+l) = (kl2/x)(x l + a6)(1-x+l) + k2 (x l + a6)
    p1 = np.polymul([(1-x)*x, (1-x)*a6 + x*a], [1.0, 1-x]); p2 = np.polymul([kl2, kl2*a6/x], [1.0, 1-x]); p3 = [k2*x, k2*a6]
    pts += real_roots(np.polysub(np.polysub(p1, p2), p3), lambda v: v > lo)
    # z6 = z2' (D1 branch, l^+ > x), z2' = kl2/x + k2/(x - l) [km1 = k2/kp, kp = x - l]:
    #   [(1-x)(x l + a6) + x a](x - l) = (kl2/x)(x l + a6)(x - l) + k2 (x l + a6)
    p1 = np.polymul([(1-x)*x, (1-x)*a6 + x*a], [-1.0, x]); p2 = np.polymul([kl2, kl2*a6/x], [-1.0, x]); p3 = [k2*x, k2*a6]
    pts += real_roots(np.polysub(np.polysub(p1, p2), p3), lambda v: v > x)
    # z7 = z2 for l^+ < 0: lt2 (1-x+l) = (kl2/x) l (1-x+l) + k2 l
    p1 = np.polymul([lt2], [1.0, 1-x]); p2 = np.polymul([kl2/x], [1.0, 1-x, 0.0]); p3 = [k2, 0.0]
    pts += real_roots(np.polysub(np.polysub(p1, p2), p3), lambda v: lo < v < 0)
    # z6p = z7 for l^+ < 0: [-x(x l + a6) + x a] l = lt2 (x l + a6)
    pts += real_roots([-x*x, -x*a6 + x*a - lt2*x, -lt2*a6], lambda v: lo < v < 0)
    return sorted(set(pts))

def halfline(cut, pre, num_is_cut_minus, zs, sides):
    """sum_j c_j [ln(cut - z_j) - i pi s_j] for num = (cut - lm) or (lm + x) [num_is_cut_minus=True -> (cut - lm)]"""
    acc = 0
    for j in range(3):
        zj = zs[j]
        resid = (cut - zj) if num_is_cut_minus else (zj - cut)      # note: (lm + x) with cut = -x is (lm - cut)
        for i in range(3):
            if i != j:
                resid = resid/(zj - zs[i])
        acc = acc + resid*(np.log(cut - zj) - 1j*PI*sides[j])
    return pre*acc

def F_lminus(lp, ltx, lty, kx, ky, qx, qy, x, eps0):
    """vectorised over the array lp"""
    lp = np.asarray(lp, dtype=float)
    k2 = kx*kx + ky*ky; kl2 = (kx+ltx)**2 + (ky+lty)**2
    a = (kx+ltx-qx)**2 + (ky+lty-qy)**2; a6 = (kx-qx)**2 + (ky-qy)**2
    N = kx*(kx+ltx-qx) + ky*(ky+lty-qy); lt2 = ltx*ltx + lty*lty
    kp = x - lp
    B = x*lp + a6
    z7 = lt2/lp; s7 = np.where(lp > 0, 1.0, -1.0)
    tot = np.zeros_like(lp, dtype=complex)
    branches = []
    km3 = -k2/(1-kp); D1_3 = kp*km3 - k2
    branches.append((2*PI*1j/(-(1-kp))/D1_3/x, (kl2 - x*km3)/x, np.ones_like(lp, dtype=bool)))
    km1 = np.where(kp != 0, k2/kp, 1.0); D3_1 = -(1-kp)*km1 - k2
    branches.append((np.where(kp != 0, 2*PI*1j/np.where(kp != 0, kp, 1.0)/D3_1/x, 0), (kl2 - x*km1)/x, kp < 0))
    for pref, z2, mask in branches:
        z6 = (1-x) + x*a/B; z6p = -x + x*a/B
        A = halfline(1-x, 2*PI*1j*pref*N/(a*B*lp), True, [z6 + 1j*eps0, z7 - 1j*s7*eps0, z2 - 1j*eps0], [-1.0, s7, 1.0])
        Bp = halfline(-x, -2*PI*1j*pref*N/(a*(-B)*lp), False, [z6p - 1j*eps0, z7 - 1j*s7*eps0, z2 - 1j*eps0], [1.0, s7, 1.0])
        tot += np.where(mask, A + Bp, 0)
    return tot

def segments(ltx, lty, kx, ky, qx, qy, x, eps0, L=1000.0):
    a6 = (kx-qx)**2 + (ky-qy)**2
    pts = {-(1-x), -a6/x, 0.0, x, 1.0, L}
    v = 5.0
    while v < L:                                  # 5, 50, 500, ... below the cap L
        pts.add(v); v *= 10.0
    for e in range(1, 9):
        for sgn in (1, -1):
            v = sgn*10.0**(-e)
            if -(1-x) < v < L: pts.add(v)
    for e in range(1, 7):
        for sgn in (1, -1):
            pts.add(-a6/x*(1 + sgn*10.0**(-e)))
    for p0 in pinch_points(ltx, lty, kx, ky, qx, qy, x):
        for m in (0.1, 0.3, 1, 3, 10, 30, 100, 1000, 10000):
            for sgn in (1, -1):
                v = p0 + sgn*m*eps0
                if -(1-x) < v < L: pts.add(v)
        pts.add(p0)
    return sorted(pts)

GL_N = 64
GL_X, GL_W = leggauss(GL_N)

def lplus_integral(ltx, lty, kx, ky, qx, qy, x, eps0, L=1000.0):
    """F(l_perp) = int_{-(1-x)}^{L} dl^+ F_lminus.  The integrand falls like 1/l^{+2} beyond l^+ ~ max(1, l_perp^2),
    so L should be >> l_perp^2 (the 2D integrator scales it)."""
    seg = segments(ltx, lty, kx, ky, qx, qy, x, eps0, L=L)
    a = np.array(seg[:-1]); b = np.array(seg[1:])
    mid = 0.5*(a + b); half = 0.5*(b - a)
    nodes = (mid[:, None] + half[:, None]*GL_X[None, :]).ravel()
    w = (half[:, None]*GL_W[None, :]).ravel()
    vals = F_lminus(nodes, ltx, lty, kx, ky, qx, qy, x, eps0)
    return np.sum(w*vals)

if __name__ == '__main__':
    import time, sys
    # validation against model_exact2 (mpmath) at the two test points
    lam = 0.03; x = 0.75; kx, ky, qx, qy = 0.6*lam, 0.3*lam, lam, 0.0
    for (ltx, lty), ref in [((0.2*lam, -0.5*lam), (-3112528.12281+66868258264.1j)), ((-3*lam, 2*lam), (-353726.71151+16876270.7777j))]:
        t0 = time.time()
        v = lplus_integral(ltx, lty, kx, ky, qx, qy, x, 1e-4*lam*lam)
        print('l_perp/lam =', (ltx/lam, lty/lam), 'numpy', v, ' mpmath', ref, ' rel', abs(v-ref)/abs(ref), f'{time.time()-t0:.2f} s')
