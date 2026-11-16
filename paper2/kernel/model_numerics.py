#!/usr/bin/env python3
"""Stage B2: what the explicit two-loop graph says (the cross-pair exchange of model_kres.py: gluon emitted by
rbar, absorbed by k; scalar lines, numerator N = k_perp.(k + l - q)_perp; x = 0.75, khat = (0.6, 0.3),
qhat = (1, 0); k_perp = lambda khat, q = lambda qhat).

F(l_perp) = int dl^+ dl^- [...] is evaluated by model_fast.py (closed-form l^- integral, Gauss-Legendre in
l^+ with windows at the pinches, finite i0 = 1e-4 lambda^2; checked against the mpmath evaluator in
model_ref_check.py).  Double precision is reliable down to lambda = 1e-3; at lambda = 3e-4 the windows of
width 1e-3 i0 near l^+ = -(1-x) reach the rounding level and the real part drifts, so lambda >= 1e-3 here.
Three statements are tested as lambda -> 0:

 (1) ultrasoft transverse momenta, l_perp = c lambda^2: lambda^6 F tends to a constant, so this region
     contributes at leading power to K = (1/4) int d^2 l_perp F ~ lambda^-2 in a single graph (it cancels
     only in the sum over attachments, us_telescoping.py);
 (2) l_perp = lambda lhat: lambda^4 F = f0(lhat) + f1(lhat) ln lambda + power corrections, a single
     logarithm of the hard scale (the rapidity logarithm), and the real part has none;
 (3) its coefficient is the eikonal one: from the soft-collinear limit of the c1 region (l^+ -> 0 at fixed
     l_perp, gluon on shell, rbar eikonal -1/(x l^+), k eikonal), int dl^+/l^+ over lambda^2 < |l^+| < 1:
       f1 = -2 (2 pi i)^3 N(l)/(x^2 k^2 a(l) l_perp^2) lambda^4,   a = |k + l - q|^2,
     i.e. the tree-level integrand at shifted momentum times 1/l_perp^2: in b space the Fourier transform of
     1/l_perp^2 times the lowest-order kernel, which in 2 - 2 eps dimensions and the b0 convention has no
     constant term (kernel_logs.py).  No Glauber, ultrasoft or pinch contribution modifies it.
    python3 model_numerics.py
"""
import numpy as np
import model_fast as MF

x = 0.75
KH, QH = (0.6, 0.3), (1.0, 0.0)

def F(lam, lt):
    kx, ky = lam*KH[0], lam*KH[1]; qx, qy = lam*QH[0], lam*QH[1]
    L = max(1000.0, 100.0*(lt[0]**2 + lt[1]**2))
    return MF.lplus_integral(lt[0], lt[1], kx, ky, qx, qy, x, 1e-4*lam*lam, L=L)

def part1():
    print('(1) l_perp = c lambda^2 (direction 0.7):  lambda^6 F')
    worst = 0
    for c in (0.1, 1.0, 3.0):
        vals = []
        for lam in (0.03, 0.01, 0.003, 0.001):
            r = c*lam*lam
            vals.append(lam**6*F(lam, (r*np.cos(0.7), r*np.sin(0.7))))
        print(f'    c = {c:4.1f}:  ' + '  '.join(f'{v.real:+.5e}{v.imag:+.5e}i' for v in vals))
        diffs = [abs(vals[i+1] - vals[i])/abs(vals[-1]) for i in range(len(vals) - 1)]
        print('             relative changes between successive lambdas: ' + ', '.join(f'{z:.1e}' for z in diffs))
        # converges to a non-zero constant: the changes shrink roughly like lambda
        assert all(diffs[i+1] < 0.6*diffs[i] for i in range(len(diffs) - 1)) and abs(vals[-1]) > 10
        worst = max(worst, diffs[-1])
    print(f'    largest change between the last two lambdas: {worst:.1e} (power corrections, shrinking like lambda)')
    return worst

def part2_3():
    print('(2), (3) l_perp = lambda lhat: lambda^4 F, the slope in ln lambda, and the eikonal prediction')
    worst = 0
    for lh in ((0.2, -0.5), (-3.0, 2.0), (0.15, 0.1), (0.9, 0.4), (-0.3, -1.1)):
        lams = (0.01, 0.003, 0.001)
        v = [lam**4*F(lam, (lam*lh[0], lam*lh[1])) for lam in lams]
        s1 = (v[1] - v[0])/np.log(lams[1]/lams[0]); s2 = (v[2] - v[1])/np.log(lams[2]/lams[1])
        k2 = KH[0]**2 + KH[1]**2
        a = (KH[0] + lh[0] - QH[0])**2 + (KH[1] + lh[1] - QH[1])**2
        N = KH[0]*(KH[0] + lh[0] - QH[0]) + KH[1]*(KH[1] + lh[1] - QH[1])
        l2 = lh[0]**2 + lh[1]**2
        pred = -2*(2j*np.pi)**3*N/(x**2*k2*a*l2)
        rel = abs(s2 - pred)/abs(pred)
        worst = max(worst, rel)
        print(f'    lhat = {lh}:  lambda^4 F = ' + ', '.join(f'{z:.5e}' for z in v))
        print(f'        slope d(lambda^4 F)/d ln lambda: {s1:.5e} -> {s2:.5e};  eikonal prediction {pred:.5e};  rel. diff {rel:.1e}')
        print(f'        real part changes by {abs(v[2].real - v[1].real)/abs(v[2]):.1e} of |lambda^4 F| (no logarithm in the real part)')
    print(f'    worst relative difference of the slope from the eikonal prediction: {worst:.1e}')
    return worst

if __name__ == '__main__':
    w1 = part1()
    w2 = part2_3()
    assert w1 < 5e-2 and w2 < 1e-2
    print('all checks passed')
