#!/usr/bin/env python3
"""Stage B2, step 1a: the cross-pair exchange graph of the two-loop DPS term (scalar model with the
DPS numerator N = k_perp.(k+l-q)_perp), and the reduction of the first loop (k) to its DPS residues.

Units s = p^+ p^- = 1, p1 = (1, 0, 0), p2 = (0, 1, 0), p3 = (x, 1-x, q) + O(lambda^2), equal masses.
Light-cone: l = (l^+, l^-, l_perp), l^2 = l^+ l^- - l_perp^2, d^4l = dl^+ dl^- d^2l_perp / 2.

Propagators (gluon l emitted by rbar = p2 - p3 + k + l, absorbed by k; vertex 1 gets k + l and
r = p3 - k - l, vertex 2 gets kbar = p1 - k and rbar - l = p2 - p3 + k):
  D1 = k^2, D2 = (k+l)^2, D3 = (p1-k)^2, D4 = (p3-k-l)^2, D5 = (p2-p3+k+l)^2, D6 = (p2-p3+k)^2, D7 = l^2.

Residue-first reduction (leading power in lambda = p_T/m, k near the DPS pinch): the k^- integral
sees only D1, D2, D3 (D4-D6 are k^- independent at leading power); closing in the upper half plane
picks D3 (always) and D1 (for l^+ > x, where k^+ = x - l^+ < 0). The k^+ integral sees D4, D5, D6
at the pinch k^+ + l^+ = x; closing in the upper half plane picks D4 (for l^- < 1-x) and D5 (for
l^- < -x). This script checks the residue formulas against direct numerical quadrature of the
k^- and k^+ integrals of the leading-power integrand, with a finite i0.
    python3 model_kres.py
"""
import mpmath as mp
mp.mp.dps = 20
I = 1j

def kminus_integrand(km, kp, lp, lm, kt, lt, i0):
    """1/(D1 D2 D3) at fixed k^+, l; kt, lt: transverse 2-vectors (mpmath matrices or tuples)"""
    kt2 = kt[0]**2 + kt[1]**2
    klt2 = (kt[0]+lt[0])**2 + (kt[1]+lt[1])**2
    D1 = kp*km - kt2 + I*i0
    D2 = (kp+lp)*(km+lm) - klt2 + I*i0
    D3 = -(1-kp)*km - kt2 + I*i0
    return 1/(D1*D2*D3)

def kminus_residues(kp, lp, lm, kt, lt):
    """2 pi i x (sum of upper-half-plane residues) of 1/(D1 D2 D3): D3 always; D1 if k^+ < 0"""
    kt2 = kt[0]**2 + kt[1]**2
    klt2 = (kt[0]+lt[0])**2 + (kt[1]+lt[1])**2
    tot = 0
    # D3 pole: km3 = -kt2/(1-kp); dD3/dkm = -(1-kp)
    km3 = -kt2/(1-kp)
    D1 = kp*km3 - kt2; D2 = (kp+lp)*(km3+lm) - klt2
    tot += 1/(-(1-kp)) * 1/(D1*D2)
    if kp < 0:   # D1 pole in the upper half plane: km1 = kt2/kp, dD1/dkm = kp
        km1 = kt2/kp
        D2 = (kp+lp)*(km1+lm) - klt2; D3 = -(1-kp)*km1 - kt2
        tot += (1/kp) * 1/(D2*D3)
    return 2*mp.pi*I*tot

def kplus_integrand(u, lp, lm, kt, lt, q, x, i0):
    """1/(D4 D5 D6) with u = k^+ - x (the small variable), k^- dropped"""
    a = (kt[0]+lt[0]-q[0])**2 + (kt[1]+lt[1]-q[1])**2
    a6 = (kt[0]-q[0])**2 + (kt[1]-q[1])**2
    D4 = -(u+lp)*((1-x)-lm) - a + I*i0
    D5 = (u+lp)*(x+lm) - a + I*i0
    D6 = u*x - a6 + I*i0
    return 1/(D4*D5*D6)

def kplus_residues(lp, lm, kt, lt, q, x):
    """2 pi i x (sum of upper-half-plane residues) of 1/(D4 D5 D6): D4 if l^- < 1-x; D5 if l^- < -x"""
    a = (kt[0]+lt[0]-q[0])**2 + (kt[1]+lt[1]-q[1])**2
    a6 = (kt[0]-q[0])**2 + (kt[1]-q[1])**2
    tot = 0
    if lm < 1-x:   # D4 pole: u4 = -lp - a/((1-x)-lm); dD4/du = -((1-x)-lm)
        u4 = -lp - a/((1-x)-lm)
        D5 = (u4+lp)*(x+lm) - a; D6 = u4*x - a6
        tot += 1/(-((1-x)-lm)) * 1/(D5*D6)
    if lm < -x:    # D5 pole: u5 = -lp + a/(x+lm) (x + lm < 0); dD5/du = (x+lm)
        u5 = -lp + a/(x+lm)
        D4 = -(u5+lp)*((1-x)-lm) - a; D6 = u5*x - a6
        tot += 1/(x+lm) * 1/(D4*D6)
    return 2*mp.pi*I*tot


def kminus_residues_i0(kp, lp, lm, kt, lt, i0):
    """the same residues with the i0 kept in the remaining denominators (for finite-i0 checks)"""
    kt2 = kt[0]**2 + kt[1]**2
    klt2 = (kt[0]+lt[0])**2 + (kt[1]+lt[1])**2
    tot = 0
    km3 = -(kt2 - I*i0)/(1-kp)
    D1 = kp*km3 - kt2 + I*i0; D2 = (kp+lp)*(km3+lm) - klt2 + I*i0
    tot += 1/(-(1-kp)) * 1/(D1*D2)
    if kp < 0:
        km1 = (kt2 - I*i0)/kp
        D2 = (kp+lp)*(km1+lm) - klt2 + I*i0; D3 = -(1-kp)*km1 - kt2 + I*i0
        tot += (1/kp) * 1/(D2*D3)
    return 2*mp.pi*I*tot

def kplus_residues_i0(lp, lm, kt, lt, q, x, i0):
    a = (kt[0]+lt[0]-q[0])**2 + (kt[1]+lt[1]-q[1])**2
    a6 = (kt[0]-q[0])**2 + (kt[1]-q[1])**2
    tot = 0
    if lm < 1-x:
        u4 = -lp - (a - I*i0)/((1-x)-lm)
        D5 = (u4+lp)*(x+lm) - a + I*i0; D6 = u4*x - a6 + I*i0
        tot += 1/(-((1-x)-lm)) * 1/(D5*D6)
    if lm < -x:
        u5 = -lp + (a - I*i0)/(x+lm)
        D4 = -(u5+lp)*((1-x)-lm) - a + I*i0; D6 = u5*x - a6 + I*i0
        tot += 1/(x+lm) * 1/(D4*D6)
    return 2*mp.pi*I*tot

def check(label, exact_fun, res_fun, pts, i0s):
    worst = 0
    for p in pts:
        vals = [exact_fun(i0, *p) for i0 in i0s]
        # extrapolate linearly in i0 to 0 from the two smallest
        ex = 2*vals[-1] - vals[-2]
        r = res_fun(*p)
        worst = max(worst, abs(ex - r)/abs(r))
    print(f'{label}: worst relative difference residues vs quadrature (i0 -> 0): {mp.nstr(worst, 3)}')
    return worst

if __name__ == '__main__':
    x = mp.mpf('0.75')
    lam = mp.mpf('0.02')
    kt = (lam*mp.mpf('0.6'), lam*mp.mpf('0.3')); q = (lam, mp.mpf(0))
    # ---- k^- integral at a few (k^+, l) including the branch k^+ < 0
    pts = []
    for kp, lp, lm, lt in [(mp.mpf('0.7'), mp.mpf('0.05'), lam**2*mp.mpf('0.4'), (lam*mp.mpf('0.2'), -lam*mp.mpf('0.5'))),
                           (mp.mpf('0.3'), mp.mpf('0.45'), lam*mp.mpf('0.3'), (lam*mp.mpf('-0.4'), lam*mp.mpf('0.1'))),
                           (mp.mpf('-0.1'), mp.mpf('0.85'), lam**2*mp.mpf('2'), (lam*mp.mpf('0.3'), lam*mp.mpf('0.3')))]:
        pts.append((kp, lp, lm, lt))
    def quad_with_poles(f, poles, i0):
        pts = sorted({-mp.inf, mp.inf, *poles, *[b - 10*i0 for b in poles], *[b + 10*i0 for b in poles]})
        return mp.quad(f, pts)

    def ex_km(i0, kp, lp, lm, lt):
        f = lambda km: kminus_integrand(km, kp, lp, lm, kt, lt, i0)
        kt2 = kt[0]**2 + kt[1]**2
        poles = [-kt2/(1-kp), kt2/kp, -lm + ((kt[0]+lt[0])**2+(kt[1]+lt[1])**2)/(kp+lp)]
        return quad_with_poles(f, poles, i0)
    w1 = check('k^- integral', ex_km, lambda kp, lp, lm, lt: kminus_residues(kp, lp, lm, kt, lt), pts, [lam**2*mp.mpf('1e-3'), lam**2*mp.mpf('5e-4')])
    # ---- k^+ integral at a few (l) including the branch l^- < -x
    pts2 = [(mp.mpf('0.05'), lam**2*mp.mpf('0.4'), (lam*mp.mpf('0.2'), -lam*mp.mpf('0.5'))),
            (lam*mp.mpf('0.3'), mp.mpf('0.1'), (lam*mp.mpf('-0.4'), lam*mp.mpf('0.1'))),
            (lam**2*mp.mpf('0.5'), mp.mpf('-0.9'), (lam*mp.mpf('0.3'), lam*mp.mpf('0.3')))]
    def ex_kp(i0, lp, lm, lt):
        f = lambda u: kplus_integrand(u, lp, lm, kt, lt, q, x, i0)
        a = (kt[0]+lt[0]-q[0])**2 + (kt[1]+lt[1]-q[1])**2
        a6 = (kt[0]-q[0])**2 + (kt[1]-q[1])**2
        poles = [-lp - a/((1-x)-lm), -lp + a/(x+lm), a6/x]
        return quad_with_poles(f, poles, i0)
    w2 = check('k^+ integral', ex_kp, lambda lp, lm, lt: kplus_residues(lp, lm, kt, lt, q, x), pts2, [lam**2*mp.mpf('1e-3'), lam**2*mp.mpf('5e-4')])
    assert w1 < 1e-2 and w2 < 1e-2      # finite-i0 residuals; the convergence with i0 is checked by hand in the log
    print('residue formulas of the first loop verified')
