#!/usr/bin/env python3
"""Reference values for the l^+ integral F(l_perp) = int dl^+ F_lminus at test points, from the mpmath
standalone evaluator model_exact2.py (tanh-sinh, 30 digits), including its own pinch-point calculation,
against the double-precision Gauss-Legendre evaluator model_fast.py.  Points: the two original test
points at |l_perp| ~ lam, the soft neighbourhood (|l_perp| = 1e-2 lam and 1e-3 lam), the neighbourhood
of t = l_perp + k - q = 0 (|t| = 1e-2 lam), and a hard point |l_perp| = 30 lam with the l^+ cap raised.
Usage: python3 model_ref_check.py > model_ref_check.log  (a few minutes)."""
import time, math
import mpmath as mp
import model_exact2 as ME
import model_fast as MF

mp.mp.dps = 30      # at 20 digits tanh-sinh nodes in the narrowest windows next to l^+ = -(1-x) round onto the end point

lam = 0.03; x = 0.75
k = (0.6*lam, 0.3*lam); q = (lam, 0.0); eps0 = 1e-4*lam*lam
c = (q[0]-k[0], q[1]-k[1])
points = {
    'A (0.2,-0.5)lam': ((0.2*lam, -0.5*lam), 1000.0),
    'B (-3,2)lam': ((-3*lam, 2*lam), 1000.0),
    'soft 1e-2 lam': ((1e-2*lam*math.cos(0.7), 1e-2*lam*math.sin(0.7)), 1000.0),
    'soft 1e-3 lam': ((1e-3*lam*math.cos(0.7), 1e-3*lam*math.sin(0.7)), 1000.0),
    't=0 + 1e-2 lam': ((c[0] + 1e-2*lam*math.cos(0.7), c[1] + 1e-2*lam*math.sin(0.7)), 1000.0),
    'hard 30 lam': ((30*lam*math.cos(0.7), 30*lam*math.sin(0.7)), 1e6),
}
worst = 0.0
for name, (lt, L) in points.items():
    t0 = time.time()
    # Separate arithmetic, root calculation and quadrature in the standalone reference.
    ltm, km, qm = (mp.mpf(lt[0]), mp.mpf(lt[1])), (mp.mpf(k[0]), mp.mpf(k[1])), (mp.mpf(q[0]), mp.mpf(q[1]))
    ref = ME.lplus_integral(ltm, km, qm, mp.mpf(x), mp.mpf(eps0), L=mp.mpf(L))
    t1 = time.time()
    v = MF.lplus_integral(lt[0], lt[1], k[0], k[1], q[0], q[1], x, eps0, L=L)
    ref_c = complex(ref)
    tolerance = 1e-4 if name.startswith('t=0') else (1e-5 if name.startswith('hard') else 1e-8)
    rel = MF.check_reference(v, ref_c, rtol=tolerance)
    worst = max(worst, rel)
    print(f'{name:16s} L={L:g}  mpmath {ref_c:.10e}  ({t1-t0:.0f} s)\n{"":16s}          numpy  {v:.10e}  rel {rel:.1e}', flush=True)
print('worst relative difference', f'{worst:.1e}')
