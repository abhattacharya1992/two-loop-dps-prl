#!/usr/bin/env python3
"""Stage B2, step 1c (second version): the cross-pair contribution to the one-loop DPS kernel,
  K_cross(k, q) = (1/4) int dl^+ d^2l_perp  [ int dl^- N R^- R^+ / D7 ],
with the l^- integral in closed form (model_lminus.py) but the i0's kept as finite imaginary parts
eps0 of the pole positions.  Reason: where an upper-half-plane pole and a lower-half-plane pole of
the l^- integrand pinch inside its integration range (z6 = z2), the l^- integral has a simple
pole in l^+ on the real axis whose i0 side is fixed by the pole positions; evaluating the closed
form with real pole positions loses that information and the l^+ quadrature becomes
prescription-dependent (observed: three quadrature rules gave three answers).  With finite eps0 the
l^+ integrand is a smooth function with Lorentzian peaks of width eps0 at the pinch points, which
are located by solving the pinch equations and resolved by breakpoints; the limit eps0 -> 0 is
taken by extrapolation.  mpmath throughout.
    python3 model_exact2.py
"""
import mpmath as mp
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
                                      (-x, -2*PI*I*pref*N/(a*(-B)*lp), lambda lm: lm + x, [(z6p, -1), (z7, s7), (z2, +1)])):
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
    """Opposite-side poles that coincide inside the l^- integration range.

    Only z6 = z2 in the piece l^- < 1-x contributes. It requires
    -(1-x) < l^+ < -a6/x and z2 < 1-x. The other coincidences are
    outside their integration ranges or involve poles on the same side.
    """
    x = mp.mpf(x)
    k2 = k[0]**2 + k[1]**2
    kl2 = (k[0]+lt[0])**2 + (k[1]+lt[1])**2
    a = (k[0]+lt[0]-q[0])**2 + (k[1]+lt[1]-q[1])**2
    a6 = (k[0]-q[0])**2 + (k[1]-q[1])**2
    c = 1-x-kl2/x
    A = c*x
    B = c*(a6+x*(1-x))+x*a-k2*x
    C = c*a6*(1-x)+x*a*(1-x)-k2*a6
    if A == 0:
        roots = [-C/B] if B else []
    else:
        discriminant = B*B-4*A*C
        if discriminant < 0:
            return []
        square_root = mp.sqrt(discriminant)
        numerator = -(B + (square_root if B >= 0 else -square_root))/2
        roots = [numerator/A, C/numerator] if numerator else [mp.mpf(0)]
    return sorted({v for v in roots
                   if -(1-x) < v < -a6/x and kl2/x+k2/(1-x+v) < 1-x})


def segments(lt, k, q, x, eps0, L=mp.mpf(1000)):
    """Resolve pinches, endpoint crossings and B=0 at the working precision."""
    x, eps0, L = mp.mpf(x), mp.mpf(eps0), mp.mpf(L)
    if not all(mp.isfinite(v) for v in (x, eps0, L)) or not 0 < x < 1 or eps0 <= 0 or L <= 0:
        raise ValueError('Require 0 < x < 1, positive eps0 and a positive finite cap')
    lo = -(1-x)
    a6 = (k[0]-q[0])**2 + (k[1]-q[1])**2
    k2 = k[0]**2 + k[1]**2
    kl2 = (k[0]+lt[0])**2 + (k[1]+lt[1])**2
    lt2 = lt[0]**2 + lt[1]**2
    pts = {lo, L}

    def add(value):
        if lo < value < L:
            pts.add(value)

    for value in (-a6/x, mp.mpf(0), x, mp.mpf(1)):
        add(value)
    value = mp.mpf(5)
    while value < L:
        add(value)
        value *= 10
    for e in range(1, 9):
        for sign in (-1, 1):
            add(sign*mp.mpf(10)**(-e))
    for e in range(1, 7):
        for sign in (-1, 1):
            add(-a6/x*(1+sign*mp.mpf(10)**(-e)))
    centres = pinch_points(lt, k, q, x) + [-a6/x, lt2/(1-x), -lt2/x]
    denominator = 1-x-kl2/x
    if denominator > 0:
        centres += [lo+k2/denominator, x+k2/denominator]
    for centre in centres:
        add(centre)
        for j in range(25):
            width = mp.power(10, -3+mp.mpf(j)/3)*eps0
            add(centre-width)
            add(centre+width)
    return sorted(pts)


def lplus_integral(lt, k, q, x, eps0, L=mp.mpf(1000)):
    lt, k, q = [tuple(mp.mpf(v) for v in pair) for pair in (lt, k, q)]
    x, eps0 = mp.mpf(x), mp.mpf(eps0)
    return mp.quad(lambda lp: lminus_closed(lp, lt, k, q, x, eps0),
                   segments(lt, k, q, x, eps0, L))


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
