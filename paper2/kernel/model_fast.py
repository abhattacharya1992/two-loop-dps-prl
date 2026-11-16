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
    """real l^+ where an upper and a lower l^- pole coincide INSIDE the integration range of their piece.
    Sides: z6 and z6p above, z2 (and z2' of the D1 branch) below, z7 below for l^+ > 0 and above for l^+ < 0.
    Piece A (l^- < 1-x): z6 is inside the range only for B = x l^+ + a6 < 0; there z7 is above like z6, so the only
    pinch is z6 = z2 (needs z2 < 1-x), with l^+ in (-(1-x), -a6/x).  Piece B (l^- < -x): z2 > 0 is outside,
    z7 (inside only for l^+ < 0) is above like z6p: no pinch.  The D1 branch (l^+ > x) has B > 0: no pinch.
    [Until 16 Oct 2026 z6p was taken below, which produced a fake pinch z6p = z7 on the line l_perp || k - q.]"""
    k2 = kx*kx + ky*ky; kl2 = (kx+ltx)**2 + (ky+lty)**2
    a = (kx+ltx-qx)**2 + (ky+lty-qy)**2; a6 = (kx-qx)**2 + (ky-qy)**2
    lo = -(1-x); hi = -a6/x
    # z6 = z2, z2 = kl2/x + k2/(1-x+l):  [(1-x)(x l + a6) + x a](1-x+l) = (kl2/x)(x l + a6)(1-x+l) + k2 (x l + a6)
    p1 = np.polymul([(1-x)*x, (1-x)*a6 + x*a], [1.0, 1-x]); p2 = np.polymul([kl2, kl2*a6/x], [1.0, 1-x]); p3 = [k2*x, k2*a6]
    roots = np.roots(np.polysub(np.polysub(p1, p2), p3))
    pts = []
    for v in roots:
        if abs(v.imag) < 1e-10*max(1.0, abs(v.real)) and lo < v.real < hi:
            z2 = kl2/x + k2/(1-x+v.real)
            if z2 < 1-x:
                pts.append(v.real)
    return sorted(set(pts))

def kinks(ltx, lty, kx, ky, qx, qy, x):
    """l^+ where a pole crosses the end point of its range ((c - z) ln(c - z) behaviour): z7 = 1-x, z7 = -x,
    z2 = 1-x, z2' = 1-x.  Used as plain breakpoints."""
    k2 = kx*kx + ky*ky; kl2 = (kx+ltx)**2 + (ky+lty)**2; lt2 = ltx*ltx + lty*lty
    pts = [lt2/(1-x), -lt2/x]
    den = (1-x) - kl2/x
    if den > 0:
        pts += [-(1-x) + k2/den, x + k2/den]
    return pts

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
        Bp = halfline(-x, -2*PI*1j*pref*N/(a*(-B)*lp), False, [z6p + 1j*eps0, z7 - 1j*s7*eps0, z2 - 1j*eps0], [-1.0, s7, 1.0])   # z6p above
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
    # windows around the pinches and around B = 0 (l^+ = -a6/x, where z6 and z6p run to infinity): three per
    # decade from 1e-2 eps0 to 1e5 eps0.  Near t = 0 the Glauber pinch sits at -a6/x - a/(1-x), within a few
    # eps0 of B = 0, and the coarser windows used before (0.1 ... 1e4 eps0, one per half decade) lost the real
    # part of F there (checked against tanh-sinh, model_ref_check.py).
    # ... and around the kinks, where (c - z_j) ln(c - z_j + i s eps0) has large derivatives on the scale eps0.  Near
    # t = 0 the kink of z7 at the cut -x (l^+ = -l_perp^2/x), the point B = 0 and the Glauber pinch all lie within
    # O(|t| lam) of each other, and without these windows the real part of F was off by several per cent there.
    centres = list(pinch_points(ltx, lty, kx, ky, qx, qy, x)) + [-a6/x] + kinks(ltx, lty, kx, ky, qx, qy, x)
    for p0 in centres:
        for m in WINDOW_M:
            for sgn in (1, -1):
                v = p0 + sgn*m*eps0
                if -(1-x) < v < L: pts.add(v)
        if -(1-x) < p0 < L: pts.add(p0)
    return sorted(pts)

GL_N = 64
GL_X, GL_W = leggauss(GL_N)
WINDOW_M = np.geomspace(1e-3, 1e5, 25)      # window offsets in units of eps0 (3 per decade)

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

def check_reference(value, reference, rtol=1e-8):
    """Require a finite numerical comparison within the stated tolerance."""
    if not np.isfinite(rtol) or rtol < 0 or not np.isfinite(value) or not np.isfinite(reference):
        raise ValueError('Nonfinite comparison or invalid tolerance')
    error = abs(value-reference)
    if error > rtol*abs(reference):
        raise ValueError(f'Reference comparison exceeds relative tolerance {rtol:g}')
    return error/abs(reference) if reference else error


if __name__ == '__main__':
    import time, sys
    # Standalone model_exact2 reference values at 30 decimal digits, L=1000.
    lam = 0.03; x = 0.75; kx, ky, qx, qy = 0.6*lam, 0.3*lam, lam, 0.0
    for (ltx, lty), ref in [((0.2*lam, -0.5*lam), (29315454560.610445+31675394086.539404j)), ((-3*lam, 2*lam), (21294354.473517448-1780295.6052552994j))]:
        t0 = time.time()
        v = lplus_integral(ltx, lty, kx, ky, qx, qy, x, 1e-4*lam*lam)
        rel = check_reference(v, ref)
        print('l_perp/lam =', (ltx/lam, lty/lam), 'numpy', v, ' mpmath', ref, ' rel', rel, f'{time.time()-t0:.2f} s')
