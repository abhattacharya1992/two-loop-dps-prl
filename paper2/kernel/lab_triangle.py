#!/usr/bin/env python3
"""Stage B2.1, laboratory: expansion by regions of the one-loop massless triangle with two
off-shell legs in the Sudakov limit, p1^2 = -a, p2^2 = -b, |s| = c, a, b << c, against the exact
result. Normalisation: I3 = int d^dl/(i pi^{d/2}) 1/[(l^2+i0)((l+p1)^2+i0)((l-p2)^2+i0)], d = 4 - 2 eps.

Regions (derived in light-cone coordinates, p1 along n, p2 along nbar):
  hard       (l ~ 1):                 I_h  = -c_G/eps^2 (-s - i0)^{-1-eps}
  collinear-1 (l ~ (1, lam^2, lam)):  I_c1 = -c_G/eps^2 a^{-eps}/s        (eikonal 1/(-2 l.p2))
  collinear-2:                        I_c2 = -c_G/eps^2 b^{-eps}/s
  ultrasoft  (l ~ lam^2):             I_us = -Gamma(eps)^2 Gamma(1-eps) (a b)^{-eps} (-s - i0)^{eps-1}
with c_G = Gamma(1+eps) Gamma(1-eps)^2/Gamma(1-2eps). The soft region (lam, lam, lam) and the
Glauber region are scaleless here. The ultrasoft region keeps the off-shellness in the eikonal
denominators, (l + p1)^2 -> p1^2 + p1^+ l^-, and carries the Coulomb phase for timelike s
(the pole l^+ = b/p2^- lies on the contour). The sum is finite and must equal the exact
triangle at leading power; checked numerically for spacelike and timelike s.

    python3 lab_triangle.py
"""
import mpmath as mp
import sympy as sp

mp.mp.dps = 25

# ------------------------------------------------------------------ exact triangle (Feynman parameters, d = 4)
def exact_triangle(a, b, s, i0=None):
    """I3 = - int dy dz / (Delta - i0), Delta = x y a + x z b - y z s, x = 1 - y - z"""
    if s < 0:
        f = lambda y, z: 1/((1 - y - z)*y*a + (1 - y - z)*z*b - y*z*s)
        return -mp.quad(lambda y: mp.quad(lambda z: f(y, z), [0, 1 - y]), [0, 1])
    if i0 is None:
        # Integrate z analytically, keeping its principal value and physical
        # discontinuity. Use stable roots and boundary products near y=0,1.
        with mp.extradps(20):
            def inner(y):
                length=1-y
                linear=b*length-(a+s)*y
                constant=a*y*length
                disc=mp.sqrt(linear**2+4*b*constant)
                if linear>=0:
                    positive=(linear+disc)/(2*b)
                    negative=-constant/(b*positive)
                else:
                    negative=(linear-disc)/(2*b)
                    positive=-constant/(b*negative)
                ratio=s*y*length/(b*(length-negative)*positive)
                real=(mp.log(ratio)-mp.log((length-negative)/(-negative)))/disc
                return real-1j*mp.pi/disc
            endpoint=mp.mpf(10)**(-mp.mp.dps+10)
            value=mp.quad(inner,[endpoint,mp.mpf('.00001'),mp.mpf('.001'),
                                mp.mpf('.01'),mp.mpf('.1'),mp.mpf('.5'),
                                mp.mpf('.9'),mp.mpf('.999'),1-endpoint])
        return +value
    # timelike: Delta changes sign on a curve; 1/(Delta - i0) = PV + i pi delta(Delta).  Use a finite
    # i0 and extrapolate (the error is O(i0 ln i0)); good to ~1e-4 relative
    vals = []
    for e0 in i0:
        f = lambda y, z: 1/((1 - y - z)*y*a + (1 - y - z)*z*b - y*z*s - 1j*e0)
        vals.append(-mp.quad(lambda y: mp.quad(lambda z: f(y, z), mp.linspace(0, 1 - y, 6)), mp.linspace(0, 1, 6)))
    return vals

# ------------------------------------------------------------------ regions, expanded in eps
e = sp.symbols('epsilon')
a, b, c = sp.symbols('a b c', positive=True)
cG = sp.gamma(1 + e)*sp.gamma(1 - e)**2/sp.gamma(1 - 2*e)

def regions(s_sign):
    """s = s_sign * c; returns the eps-expansion (through eps^0) of the sum of regions times c"""
    if s_sign < 0:              # spacelike s = -c: (-s - i0) = c
        Ih = -cG/e**2*c**(-1 - e)
        Ic1 = -cG/e**2*a**(-e)/(-c)
        Ic2 = -cG/e**2*b**(-e)/(-c)
        Ius = -sp.gamma(e)**2*sp.gamma(1 - e)*(a*b)**(-e)*c**(e - 1)
    else:                       # timelike s = +c: (-s - i0) = c e^{-i pi}; the phase is written out
        Ih = -cG/e**2*c**(-1 - e)*sp.exp(sp.I*sp.pi*(1 + e))         # (c e^{-i pi})^{-1-eps}
        Ic1 = -cG/e**2*a**(-e)/c
        Ic2 = -cG/e**2*b**(-e)/c
        Ius = -sp.gamma(e)**2*sp.gamma(1 - e)*(a*b)**(-e)*c**(e - 1)*sp.exp(-sp.I*sp.pi*(e - 1))   # (c e^{-i pi})^{eps-1}
    tot = sp.expand(sp.series(sp.expand(c*(Ih + Ic1 + Ic2 + Ius)), e, 0, 1).removeO())
    return sp.simplify(tot)

for s_sign, label in ((-1, 'spacelike s = -c'), (+1, 'timelike s = +c')):
    tot = regions(s_sign)
    tot = sp.expand(sp.expand_log(tot, force=True))
    print(f'{label}:  c * sum of regions = {tot}')
    assert not tot.has(e), 'poles do not cancel'
    # numerical comparison at a = b = 1e-3 c, and a = 2e-3, b = 5e-4
    for av, bv in ((1e-3, 1e-3), (2e-3, 5e-4)):
        pred = complex(tot.subs({a: av, b: bv, c: 1}))
        if s_sign < 0:
            ex = exact_triangle(mp.mpf(av), mp.mpf(bv), mp.mpf(-1))
            print(f'   a={av}, b={bv}: regions {pred:.6f}   exact {float(ex):.6f}   rel.diff {abs(pred-float(ex))/abs(float(ex)):.1e}')
            assert abs(pred - float(ex)) < 3e-2*abs(float(ex))      # power corrections ~ a ln^2 a
        else:
            # First continue the leading-power formula, then compare separately
            # with the finite-virtuality Feynman-parameter integral.
            cont = (1/1.0)*((mp.log(av) + 1j*mp.pi)*(mp.log(bv) + 1j*mp.pi) + mp.pi**2/3)
            print(f'   a={av}, b={bv}: regions {pred:.5f}   continued leading-power result {complex(cont):.5f}')
            assert abs(pred - complex(cont)) < 1e-9
            ex=exact_triangle(mp.mpf(av),mp.mpf(bv),mp.mpf(1))
            print('   exact timelike integral',mp.nstr(ex,15),' relative power correction',mp.nstr(abs(pred-ex)/abs(ex),3))
            if not mp.isfinite(ex) or abs(pred-ex)>mp.mpf('.03')*abs(ex):
                raise ValueError('Timelike regions fail the finite-virtuality comparison')
print('laboratory passed: h + c1 + c2 + us reproduces the Sudakov limit, with the Coulomb phase from the ultrasoft region')
