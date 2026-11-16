#!/usr/bin/env python3
"""Stage B1: the one-loop soft function of the DPS kernel (Paper 2, Sec. 5; definitions note v3,
Definition 2): four incoming Wilson lines, n and nbar at transverse position 0 (lines k, r) and n
and nbar at position b (lines kbar, rbar). At one loop only the two cross pairs (k, rbar) and
(kbar, r) contribute (same-position n-nbar pairs are scaleless, n-n pairs vanish for light-like
lines), each with colour T_k.T_rbar = T_kbar.T_r = Delta_C/2.

The virtual n-nbar exchange at transverse separation b between two incoming lines, with the
delta-regulator on the eikonal propagators (1/(l^- + i delta1), 1/(-l^+ + i delta2)):

  J = int d^dl/(2pi)^d  e^{-i l_perp.b} / [(l^+ l^- - l_perp^2 + i0)(l^- + i delta1)(-l^+ + i delta2)]

l^+ by residues (gluon pole, l^- > 0 only), then l^- in closed form, then the transverse Fourier
transform in d - 2 = 2 - 2 eps dimensions (the delta1 delta2 in the transverse denominator is a
power correction):

  J = i pi^{2-eps}/(2pi)^d  Gamma(-eps) (b^2/4)^eps  [ psi(-eps) + gamma_E - ln(delta1 delta2 b^2/b0^2) ]

and, with -i g^2 (n.nbar) (Delta_C/2) J per pair, two pairs, MS-bar coupling a_s = alpha_s/(4 pi):

  S^(1) = 2 a_s Delta_C (mu^2 b^2/b0^2)^eps e^{-eps gamma_E} Gamma(-eps)
          [ psi(-eps) + gamma_E - ln(delta1 delta2 b^2/b0^2) ]
       = 2 a_s Delta_C [ -1/eps^2 - L_b/eps - L_b^2/2 + ln(delta1 delta2 b^2/b0^2)(1/eps + L_b) + zeta_2/2 ] + O(eps)

with L_b = ln(mu^2 b^2/b0^2), b0 = 2 e^{-gamma_E}. The coefficient of the rapidity logarithm,
-(mu^2 b^2/b0^2)^eps e^{-eps gamma_E} Gamma(-eps) = 1/eps + L_b + O(eps), has no constant term.

Every step is checked numerically below (mpmath):
  1. the l^+ residue and the l^- integral at fixed l_perp, against direct numerical integration;
  2. the transverse Fourier master formula in two dimensions (0 < alpha < 1, where it converges);
  3. the epsilon expansion.
    python3 soft_function.py
"""
import mpmath as mp
import sympy as sp

mp.mp.dps = 30

# ------------------------------------------------------------------ 1. the longitudinal integrals
def lminus_integral_closed(A, d1, d2):
    """int_0^oo dl^- / [(A - i d2 l^-)(l^- + i d1)] = ln(A/(d1 d2)) / (A - d1 d2)"""
    return mp.log(A/(d1*d2))/(A - d1*d2)

def lminus_integral_numeric(A, d1, d2):
    f = lambda lm: 1/((A - 1j*d2*lm)*(lm + 1j*d1))
    return mp.quad(f, [0, mp.inf])

def lplus_integral_numeric(lm, A, d2, eps0):
    """int dl^+ / [(l^+ l^- - A + i eps0)(-l^+ + i d2)], with a finite i0 = i eps0; the residue
    calculation is exact for any eps0 > 0: 2 pi i/(A - i eps0 - i d2 l^-) for l^- > 0, and 0 for l^- < 0"""
    f = lambda lp: 1/((lp*lm - A + 1j*eps0)*(-lp + 1j*d2))
    pts = sorted([-mp.inf, 0, A/lm - 1, A/lm, A/lm + 1, mp.inf])
    return mp.quad(f, pts)

A, d1, d2, eps0 = mp.mpf('0.7'), mp.mpf('0.013'), mp.mpf('0.021'), mp.mpf('0.05')
lm = mp.mpf('0.4')
closed_plus = 2j*mp.pi/(A - 1j*eps0 - 1j*d2*lm)
num_plus = lplus_integral_numeric(lm, A, d2, eps0)
print('l^+ integral (l^- > 0):   closed', mp.nstr(closed_plus, 12), '  numeric', mp.nstr(num_plus, 12))
assert abs(closed_plus - num_plus) < 1e-15
num_plus_neg = lplus_integral_numeric(-lm, A, d2, eps0)
print('l^+ integral (l^- < 0):   numeric', mp.nstr(num_plus_neg, 8), ' (vanishes: both poles in the upper half plane)')
assert abs(num_plus_neg) < 1e-15
print('l^- integral:  closed', mp.nstr(lminus_integral_closed(A, d1, d2), 15), '  numeric', mp.nstr(lminus_integral_numeric(A, d1, d2), 15))
if abs(lminus_integral_closed(A,d1,d2)-lminus_integral_numeric(A,d1,d2))>mp.mpf('1e-20'):
    raise ValueError('Longitudinal l-minus integral comparison failed')

# ------------------------------------------------------------------ 2. the transverse Fourier transform (two dimensions)
def fourier_master(alpha, b, n=2):
    """int d^n l e^{-i l.b} (l^2)^{alpha-1} = pi^{n/2} 2^{n-2(1-alpha)} Gamma(n/2-1+alpha)/Gamma(1-alpha) (b^2)^{1-alpha-n/2}"""
    a = 1 - alpha
    return mp.pi**(mp.mpf(n)/2)*2**(n - 2*a)*mp.gamma(mp.mpf(n)/2 - a)/mp.gamma(a)*(b*b)**(a - mp.mpf(n)/2)

b = mp.mpf('1.3')
# in two dimensions the master formula is the Hankel transform int_0^oo dl l^mu J_0(l b) =
# 2^mu Gamma((1+mu)/2)/(b^(1+mu) Gamma((1-mu)/2)), mu = 2 alpha - 1, times 2 pi: identical in closed form;
# numerically the oscillatory integral converges slowly, so the check is made at alpha = 0.15 to 1e-6
alpha = mp.mpf('0.15')
hankel = 2*mp.pi*2**(2*alpha-1)*mp.gamma(alpha)/(b**(2*alpha)*mp.gamma(1-alpha))
assert abs(hankel - fourier_master(alpha, b)) < 1e-25
# numerically: damp with e^{-kappa l} and extrapolate kappa -> 0 (the damped integral is not analytic
# at kappa = 0, so the extrapolation is approximate: agreement to 1e-4 is what it gives)
ks = [mp.mpf(k) for k in ('0.2', '0.1', '0.05', '0.025')]
vals = [mp.quad(lambda l: 2*mp.pi*l**(2*alpha-1)*mp.exp(-k*l)*mp.besselj(0, l*b), mp.linspace(0, 200, 400) + [mp.inf]) for k in ks]
num = mp.lu_solve(mp.matrix([[1, k, k**2, k**3] for k in ks]), mp.matrix(vals))[0]
print(f'Fourier master, alpha={alpha}: closed', mp.nstr(fourier_master(alpha, b), 12), ' numeric (damped, extrapolated)', mp.nstr(num, 12))
assert abs(num - fourier_master(alpha, b)) < 1e-4*abs(num)

# ------------------------------------------------------------------ 3. the epsilon expansion of S^(1)/(2 a_s Delta_C)
e, Lb, R = sp.symbols('epsilon L_b R')            # R = ln(delta1 delta2 b^2/b0^2)
# psi(-eps) = psi(1 - eps) + 1/eps
expr = sp.exp(e*Lb)*sp.exp(-e*sp.EulerGamma)*sp.gamma(-e)*(sp.polygamma(0, 1 - e) + 1/e + sp.EulerGamma - R)
ser = sp.expand(sp.series(expr, e, 0, 1).removeO())
print('S^(1)/(2 a_s Delta_C) =', ser)
target = -1/e**2 - Lb/e - Lb**2/2 + R*(1/e + Lb) + sp.zeta(2)/2
assert sp.simplify(sp.expand(ser - target).subs(sp.zeta(2), sp.pi**2/6)) == 0
coef_R = sp.series(-sp.exp(e*Lb)*sp.exp(-e*sp.EulerGamma)*sp.gamma(-e), e, 0, 2).removeO()
print('coefficient of the rapidity logarithm:', sp.expand(coef_R), '  (no constant term)')
assert sp.expand(coef_R).coeff(e, 0) == Lb
print('checks passed')
