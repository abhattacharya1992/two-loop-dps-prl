#!/usr/bin/env python3
"""Stage B2, step 1c (second version): the cross-pair contribution to the one-loop DPS kernel,
  K_cross(k, q) = (1/4) int dl^+ d^2l_perp  [ int dl^- N R^- R^+ / D7 ],
with the l^- integral in closed form (model_lminus.py) but the i0's kept as finite imaginary parts
eps0 of the pole positions.  Reason: where an upper-half-plane pole and a lower-half-plane pole of
the l^- integrand pinch (z6 = z7, z6 = z2, or z7 = z2 for l^+ < 0), the l^- integral has a simple
pole in l^+ on the real axis whose i0 side is fixed by the pole positions; evaluating the closed
form with real pole positions loses that information and the l^+ quadrature becomes
prescription-dependent (observed: three quadrature rules gave three answers).  With finite eps0 the
l^+ integrand is a smooth function with Lorentzian peaks of width eps0 at the pinch points, which
are located by solving the pinch equations and resolved by breakpoints; the limit eps0 -> 0 is
taken by extrapolation.  mpmath throughout.
    python3 model_exact2.py
"""
import mpmath as mp
import numpy as np
mp.mp.dps = 20
I = 1j
PI = mp.pi

def lminus_closed(lp, lt, k, q, x, eps0):
    """int dl^- N R^- R^+ / D7 with the pole positions carrying -i s eps0 (s = side)"""
    k2 = k[0]**2 + k[1]**2
    kl2 = (k[0]+lt[0])**2 + (k[1]+lt[1])**2
    a = (k[0]+lt[0]-q[0])**2 + (k[1]+lt[1]-q[1])**2
    a6 = (k[0]-q[0])**2 + (k[1]-q[1])**2
    N = k[0]*(k[0]+lt[0]-q[0]) + k[1]*(k[1]+lt[1]-q[1])
    lt2 = lt[0]**2 + lt[1]**2
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
            zs = [mp.mpc(r, -s*eps0) for r, s in poles]
            acc = 0
            for j, zj in enumerate(zs):
                resid = numf(zj)
                for i, zi in enumerate(zs):
                    if i != j:
                        resid /= (zj - zi)
                sj = poles[j][1]
                acc += resid*(mp.log(cut - zj) - I*PI*sj)
            tot += pre*acc
    return tot

def pinch_points(lt, k, q, x):
    """real l^+ where two l^- poles of opposite sides coincide: z6 = z7, z6 = z2 (l^+ > 0 branch of z7 below),
    z7 = z2 (l^+ < 0), z6p = z7 (l^+ < 0: z6p above, z7 above -> no pinch; l^+ > 0: pinch), z6p = z2: same side, no pinch"""
    k2 = float(k[0]**2 + k[1]**2)
    kl2 = float((k[0]+lt[0])**2 + (k[1]+lt[1])**2)
    a = float((k[0]+lt[0]-q[0])**2 + (k[1]+lt[1]-q[1])**2)
    a6 = float((k[0]-q[0])**2 + (k[1]-q[1])**2)
    lt2 = float(lt[0]**2 + lt[1]**2)
    x = float(x)
    pts = []
    # z6 = z7: [(1-x)(x l + a6) + x a] l = lt2 (x l + a6)
    r = np.roots([(1-x)*x, (1-x)*a6 + x*a - lt2*x, -lt2*a6])
    pts += [float(v.real) for v in r if abs(v.imag) < 1e-12 and v.real > 0]
    # z6 = z2 with z2 = kl2/x + k2/(1-x+l): (1-x) + x a/(x l + a6) = kl2/x + k2/(1-x+l)
    #   multiply by (x l + a6)(1-x+l): [(1-x)(x l + a6) + x a](1-x+l) = kl2/x (x l + a6)(1-x+l) + k2 (x l + a6)
    c3 = [(1-x)*x, (1-x)*a6 + x*a, 0.0]           # (1-x)(x l + a6) + x a = (1-x)x l + (1-x)a6 + x a
    p1 = np.polymul([ (1-x)*x, (1-x)*a6 + x*a ], [1.0, 1-x])
    p2 = np.polymul([kl2, kl2*a6/x], [1.0, 1-x])
    p3 = [k2*x, k2*a6]
    r = np.roots(np.polysub(np.polysub(p1, p2), p3))
    pts += [float(v.real) for v in r if abs(v.imag) < 1e-12 and v.real > 0]
    # z7 = z2 for l^+ < 0: lt2/l = kl2/x + k2/(1-x+l)
    p1 = np.polymul([lt2], [1.0, 1-x]); p2 = np.polymul([kl2/x], np.polymul([1.0, 0.0], [1.0, 1-x])); p3 = np.polymul([k2], [1.0, 0.0])
    r = np.roots(np.polysub(np.polysub(p1, p2), p3))
    pts += [float(v.real) for v in r if abs(v.imag) < 1e-12 and -(1-x) < v.real < 0]
    # z6p = z7 for l^+ > 0: [-x(x l + a6) + x a] l = lt2 (x l + a6)
    r = np.roots([-x*x, -x*a6 + x*a - lt2*x, -lt2*a6])
    pts += [float(v.real) for v in r if abs(v.imag) < 1e-12 and v.real > 0]
    return sorted(set(pts))

def lplus_integral(lt, k, q, x, eps0, L=mp.mpf(1000)):
    a6 = (k[0]-q[0])**2 + (k[1]-q[1])**2
    f = lambda lp: lminus_closed(lp, lt, k, q, x, eps0)
    pts = {-(1-x), -a6/x, mp.mpf(0), x, mp.mpf(1), mp.mpf(5), mp.mpf(50), L}
    for e in range(1, 9):
        for sgn in (1, -1):
            v = sgn*mp.mpf(10)**(-e)
            if -(1-x) < v < L: pts.add(v)
    for e in range(1, 7):
        for sgn in (1, -1):
            pts.add(-a6/x*(1 + sgn*mp.mpf(10)**(-e)))
    for p0 in pinch_points(lt, k, q, x):
        for m in (1, 3, 10, 30, 100, 1000):
            for sgn in (1, -1):
                v = mp.mpf(p0) + sgn*m*eps0
                if -(1-x) < v < L: pts.add(v)
        pts.add(mp.mpf(p0))
    return mp.quad(f, sorted(pts))

if __name__ == '__main__':
    import time
    x = mp.mpf('0.75'); lam = mp.mpf('0.03')
    k = (mp.mpf('0.6')*lam, mp.mpf('0.3')*lam); q = (lam, mp.mpf(0))
    for lt in [(mp.mpf('0.2')*lam, -mp.mpf('0.5')*lam), (-3*lam, 2*lam)]:
        print('pinch points', pinch_points(lt, k, q, x))
        for eps0 in (lam**2*mp.mpf('1e-3'), lam**2*mp.mpf('3e-4'), lam**2*mp.mpf('1e-4')):
            t0 = time.time()
            v = lplus_integral(lt, k, q, x, eps0)
            print([float(c/lam) for c in lt], 'eps0', mp.nstr(eps0, 3), mp.nstr(v, 12), f'{time.time()-t0:.0f} s')
