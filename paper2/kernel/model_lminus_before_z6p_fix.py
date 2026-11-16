#!/usr/bin/env python3
"""Stage B2, step 1b: the cross-pair graph after the first loop's residues (model_kres.py):

  K_cross(k, q) = (1/4) int dl^+ dl^- d^2l_perp  N R^-(l) R^+(l) / D7,
  N = k.(k+l_perp-q),  D7 = l^+ l^- - l_perp^2 + i0,
  R^- = 2 pi i [Res_{D3} + theta(l^+ - x) Res_{D1}] of 1/(D1 D2 D3) at k^+ = x - l^+,
  R^+ = 2 pi i [theta(1-x-l^-) Res_{D4} + theta(-x-l^-) Res_{D5}] of 1/(D4 D5 D6).

The l^- integrand is a rational function of l^- (all poles simple, each with a definite side of
the real axis from the i0's), integrated over the half-line (-oo, 1-x) with the extra piece on
(-oo, -x). For a rational R(z) = sum_j c_j/(z - z_j) with sum_j c_j = 0 (R ~ 1/z^2),
  int_{-oo}^{c} R dz = sum_j c_j [ ln|c - r_j| - i pi s_j theta(c - r_j) ],   z_j = r_j - i s_j 0,
s_j = +1 for a pole below the real axis, -1 above. This closed form is checked against direct
quadrature with a finite i0, then used in model_exact.py.
    python3 model_lminus.py
"""
import mpmath as mp
mp.mp.dps = 20
I = 1j

def pieces(lp, lt, k, q, x):
    """the l^- rational integrand as a list of (poles, side, prefactor-function) pieces:
    returns a list of (c_cut, numerator_fn, denominators) where each piece is
    prefactor * num(lm) / prod (lm - z_j), integrated over (-oo, c_cut)"""
    k2 = k[0]**2 + k[1]**2
    kl2 = (k[0]+lt[0])**2 + (k[1]+lt[1])**2
    a = (k[0]+lt[0]-q[0])**2 + (k[1]+lt[1]-q[1])**2
    a6 = (k[0]-q[0])**2 + (k[1]-q[1])**2
    N = k[0]*(k[0]+lt[0]-q[0]) + k[1]*(k[1]+lt[1]-q[1])
    lt2 = lt[0]**2 + lt[1]**2
    kp = x - lp
    out = []
    # ---- R^- pieces: (jacobian, D1-or-D3 value, D2 as linear function of l^-): D2 = x (lm - z2)
    rminus = []
    km3 = -k2/(1-kp)
    D1_3 = kp*km3 - k2                                  # D1 at the D3 residue (l^- independent)
    z2_3 = (kl2 - x*km3)/x                              # D2 = x(lm - z2) with z2 = z2_3 - i0  (below)
    rminus.append((2*mp.pi*I/(-(1-kp)) / D1_3 / x, z2_3))
    if kp < 0:
        km1 = k2/kp
        D3_1 = -(1-kp)*km1 - k2
        z2_1 = (kl2 - x*km1)/x
        rminus.append((2*mp.pi*I/kp / D3_1 / x, z2_1))
    # ---- R^+ pieces: 1/(a D6(u4)) on (-oo, 1-x); -1/(a D6(u5)) on (-oo, -x)
    #   D6(u4) = -x lp - a6 - x a/((1-x) - lm) = -(x lp + a6) (z6 - lm)/((1-x) - lm),  z6 = (1-x) + x a/(x lp + a6) + i0 (above)
    #   D6(u5) = -x lp - a6 + x a/(x + lm)     = -(x lp + a6) (lm - z6')/(x + lm),      z6' = -x + x a/(x lp + a6) - i0 ?
    B = x*lp + a6
    # D7 = lp (lm - z7), z7 = lt2/lp - i0 for lp > 0 (below), above for lp < 0
    z7 = lt2/lp
    s7 = +1 if lp > 0 else -1
    for pref, z2 in rminus:
        # piece A: on (-oo, 1-x): pref * N * [((1-x) - lm)/(-B (z6 - lm))] / (a) / (lp (lm - z7)) / (lm - z2)
        #        = pref*N/(a*(-B)*lp) * ((1-x) - lm) / [(lm - z6)(lm - z7)(lm - z2)] * (-1)   [since (z6 - lm) = -(lm - z6)]
        z6 = (1-x) + x*a/B
        out.append(dict(cut=1-x, pref=2*mp.pi*I*pref*N/(a*B*lp), num=lambda lm, c=1-x: (c - lm),
                        poles=[(z6, -1), (z7, s7), (z2, +1)]))
        if True:   # piece B on (-oo, -x): -1/(a D6(u5)); 1/D6(u5) = (x + lm)/[-(B)(lm - z6')], z6' = -x + x a/B
            z6p = -x + x*a/B
            # sign of i0 of z6': D6(u5) = -x lp + x(a - i0)/(x+lm) - a6 + i0 ; near the zero, as a function of lm the
            # imaginary part is -i0 x/(x+lm) + i0 with x+lm < 0 on this branch: both positive -> D6 + i0 -> root below? see check
            out.append(dict(cut=-x, pref=-2*mp.pi*I*pref*N/(a*(-B)*lp), num=lambda lm, c=-x: (lm + x),
                            poles=[(z6p, +1), (z7, s7), (z2, +1)]))
    return out

def halfline_integral(piece):
    """int_{-oo}^{cut} pref * num(lm) / prod (lm - z_j) dlm by partial fractions with the branch rule"""
    poles = piece['poles']; c = piece['cut']; num = piece['num']
    tot = 0
    for j, (zj, sj) in enumerate(poles):
        resid = num(zj)
        for i, (zi, si) in enumerate(poles):
            if i != j:
                resid /= (zj - zi)
        tot += resid*(mp.log(abs(c - zj)) - I*mp.pi*sj*(1 if c > zj else 0))
    return piece['pref']*tot

def lminus_closed(lp, lt, k, q, x):
    return sum(halfline_integral(p) for p in pieces(lp, lt, k, q, x))

def lminus_numeric(lp, lt, k, q, x, i0):
    """direct quadrature of the l^- integrand with finite i0, built from the residue formulas"""
    from model_kres import kminus_residues_i0, kplus_residues_i0
    N = k[0]*(k[0]+lt[0]-q[0]) + k[1]*(k[1]+lt[1]-q[1])
    lt2 = lt[0]**2 + lt[1]**2
    def f(lm):
        D7 = lp*lm - lt2 + I*i0
        return N*kminus_residues_i0(x-lp, lp, lm, k, lt, i0)*kplus_residues_i0(lp, lm, k, lt, q, x, i0)/D7
    k2 = k[0]**2 + k[1]**2; kl2 = (k[0]+lt[0])**2 + (k[1]+lt[1])**2
    a = (k[0]+lt[0]-q[0])**2 + (k[1]+lt[1]-q[1])**2; a6 = (k[0]-q[0])**2 + (k[1]-q[1])**2
    kp = x - lp; B = x*lp + a6
    pts = {lt2/lp, (kl2 + x*k2/(1-kp))/x, (1-x) + x*a/B, -x + x*a/B, 1-x, -x}
    if kp < 0: pts.add((kl2 - x*k2/kp)/x)
    pts = sorted(p for p in pts if p < 1-x)
    brk = sorted({-mp.inf, 1-x, *pts, *[p - 10*i0 for p in pts], *[p + 10*i0 for p in pts]})
    return mp.quad(f, brk)

if __name__ == '__main__':
    x = mp.mpf('0.75'); lam = mp.mpf('0.02')
    k = (lam*mp.mpf('0.6'), lam*mp.mpf('0.3')); q = (lam, mp.mpf(0))
    tests = [(mp.mpf('0.05'), (lam*mp.mpf('0.2'), -lam*mp.mpf('0.5'))),
             (lam*mp.mpf('0.3'), (lam*mp.mpf('-0.4'), lam*mp.mpf('0.1'))),
             (-mp.mpf('0.1'), (lam*mp.mpf('0.3'), lam*mp.mpf('0.3'))),
             (mp.mpf('0.85'), (lam*mp.mpf('0.1'), lam*mp.mpf('0.2'))),
             (lam**2*mp.mpf('0.7'), (lam*mp.mpf('0.5'), -lam*mp.mpf('0.2')))]
    worst = 0
    for lp, lt in tests:
        cl = lminus_closed(lp, lt, k, q, x)
        nums = [lminus_numeric(lp, lt, k, q, x, i0) for i0 in (lam**2*mp.mpf('1e-3'), lam**2*mp.mpf('5e-4'))]
        ex = 2*nums[1] - nums[0]
        rel = abs(cl - ex)/abs(cl)
        worst = max(worst, rel)
        print(f'l+={mp.nstr(lp,4)}: closed {mp.nstr(cl, 8)}  quad(i0->0) {mp.nstr(ex, 8)}  rel {mp.nstr(rel, 2)}')
    print('worst', mp.nstr(worst, 3))
