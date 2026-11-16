#!/usr/bin/env python3
"""Stage B2: which regions can produce the boundary terms r2 X^2 + r1 X of the one-loop kernel.

The bare one-loop kernel, normalised to the d-dimensional lowest-order kernel, is a sum of regions.  At one
loop each region is homogeneous in its scale; after the rapidity regulator has cancelled, two kinds remain:

  (sc)  collinear-1, collinear-2, soft, Glauber: virtuality mu_b^2 = b0^2/b^2 (they resolve b), with the
        hard scale entering only through ONE rapidity logarithm L_Q = ln(Q^2/mu_b^2):
          (mu^2/mu_b^2)^eps [ c/eps^2 + (d + g(eps)... ) ],  g(eps) L_Q = (g/eps + g1) L_Q;
  (us)  ultrasoft (l ~ lambda^2): virtuality sigma_1 sigma_2 ~ mu_b^4/Q^2, no rapidity divergence:
          (mu^2 Q^2 kappa_u/mu_b^4)^eps [ cu/eps^2 + du/eps + eu ].

Minimal subtraction removes the poles (which must not depend on b); the finite part is then written in
X_mu = ln(mu^2/mu_b^2) and X = ln(Q^2/mu_b^2) and compared with the structure of Paper 2, eq. (5.3),
  D X_mu^2 - 2 D X X_mu + B X_mu + r2 X^2 + r1 X + r0.
Result (part 1):  r2 = cu/2,  r1 = g1 + du + cu ln(kappa_u):  only the ultrasoft regions can give r2, and
r1 comes from the constant g1 of the rapidity coefficient and from the ultrasoft regions.
Part 2 checks the bookkeeping on the laboratory triangle (lab_triangle.py), whose 'kernel' c1 + c2 + us
(the hard region removed) has an X^2 term equal to cu/2 because its ultrasoft region does not cancel.
Part 3 applies it to the DPS kernel: summed over the attachments to {p1, k, kbar} x {p2, r, rbar} the
ultrasoft region is scaleless (us_telescoping.py), so cu = du = 0, hence r2 = 0 and r1 = g1; the rapidity
coefficient is that of the soft function, 2 a_s Delta_C (mu^2 b^2/b0^2)^eps e^{-eps gamma_E} Gamma(-eps)
(soft_function.py), whose expansion has no eps^0 term, so g1 = 0 and r1 = 0.
    python3 kernel_logs.py
"""
import sympy as sp

eps = sp.symbols('epsilon')
Lmu, Lb, LQ = sp.symbols('L_mu L_b L_Q')            # ln mu^2, ln mu_b^2, ln Q^2 (Q^2 = Q_*^2 e^{-i pi})
c, d, e0, g, g1, cu, du, eu, ku, kr = sp.symbols('c d e0 g g1 c_u d_u e_u kappa_u kappa_r')
D, B, r0 = sp.symbols('Delta_C B r0')

def laurent(expr, n=0):
    return sp.expand(sp.series(expr, eps, 0, n + 1).removeO())

def part1():
    print('part 1: the general bookkeeping')
    sc = sp.exp(eps*(Lmu - Lb))*(c/eps**2 + d/eps + e0 + (g/eps + g1)*(LQ - Lb + kr))
    us = sp.exp(eps*(Lmu + LQ + ku - 2*Lb))*(cu/eps**2 + du/eps + eu)
    tot = laurent(sc + us)
    p2, p1, f = tot.coeff(eps, -2), tot.coeff(eps, -1), tot.coeff(eps, 0)
    # the poles must not depend on b (IR factorisation of the 2 -> 2 amplitude and of H1, H2)
    cond = sp.solve(sp.Eq(sp.diff(p1, Lb), 0), g)
    print('   b-independence of the single pole:  g =', cond[0])
    f = sp.expand(f.subs(g, cond[0]))
    X, Xmu = sp.symbols('X X_mu')
    # Lmu = Xmu + Lb, LQ = X + Lb: the finite part must then be independent of Lb
    fx = sp.expand(f.subs({Lmu: Xmu + Lb, LQ: X + Lb}))
    assert sp.simplify(sp.diff(fx, Lb)) == 0
    poly = sp.Poly(fx, X, Xmu)
    co = {m: sp.simplify(poly.coeff_monomial(m)) for m in (Xmu**2, X*Xmu, X**2, Xmu, X, 1)}
    print('   finite part:  X_mu^2:', co[Xmu**2], '  X X_mu:', co[X*Xmu], '  X^2:', co[X**2])
    print('                 X_mu:', co[Xmu], '  X:', co[X], '  1:', co[1])
    # compare with D X_mu^2 - 2 D X X_mu + B X_mu + r2 X^2 + r1 X + r0
    sol = sp.solve([sp.Eq(co[Xmu**2], D), sp.Eq(co[X*Xmu], -2*D)], [c, D], dict=True)[0]
    print('   matching X_mu^2 and X X_mu (both fixed by gamma_K):', sol, ' (consistent: one relation, c + cu = 2 Delta_C)')
    r2 = co[X**2]; r1 = co[X]
    print('   => r2 =', r2, ',  r1 =', r1)
    assert sp.simplify(r2 - cu/2) == 0 and sp.simplify(r1 - (g1 + du + cu*ku)) == 0
    return r2, r1

def part2():
    print('part 2: the bookkeeping on the laboratory triangle (spacelike s = -Q^2, legs a = b = mu_b^2)')
    # regions of lab_triangle.py, times Q^2, with a scale mu restored as (mu^2)^eps
    cG = sp.gamma(1 + eps)*sp.gamma(1 - eps)**2/sp.gamma(1 - 2*eps)
    mu2, mb2, Q2 = sp.symbols('mu2 mb2 Q2', positive=True)
    Ic1 = cG/eps**2*(mu2/mb2)**eps
    Ic2 = cG/eps**2*(mu2/mb2)**eps
    Ius = -sp.gamma(eps)**2*sp.gamma(1 - eps)*(mu2*Q2/mb2**2)**eps
    Ih = -cG/eps**2*(mu2/Q2)**eps
    kern = laurent(Ic1 + Ic2 + Ius)          # the 'kernel': the hard region removed
    full = laurent(Ic1 + Ic2 + Ius + Ih)
    assert not full.has(sp.Pow(eps, -1)) and full.coeff(eps, -2) == 0
    cu_tri = sp.limit(eps**2*(-sp.gamma(eps)**2*sp.gamma(1 - eps)), eps, 0)
    fin = kern.coeff(eps, 0)
    X, Xmu = sp.symbols('X X_mu')
    fin = sp.expand(sp.expand_log(fin, force=True))
    fx = sp.expand(fin.subs({sp.log(mu2): Xmu + sp.log(mb2), sp.log(Q2): X + sp.log(mb2)}))
    fx = sp.expand(sp.expand_log(fx, force=True))
    assert sp.simplify(sp.diff(fx, sp.Symbol('mb2'))) == 0 or not fx.has(sp.log(mb2))
    poly = sp.Poly(fx, X, Xmu)
    x2 = poly.coeff_monomial(X**2)
    print(f'   ultrasoft double pole c_u = {cu_tri};  X^2 coefficient of the finite part = {x2};  c_u/2 = {cu_tri/2}')
    assert sp.simplify(x2 - cu_tri/2) == 0
    print('   (with the hard region added the poles cancel and the X^2 term is that of ln(a/Q^2) ln(b/Q^2))')

def part3():
    print('part 3: the DPS kernel')
    print('   ultrasoft region summed over attachments: scaleless (us_telescoping.py)  =>  c_u = d_u = e_u = 0')
    print('   => r2 = 0,  r1 = g1 (the eps^0 term of the rapidity coefficient at mu = mu_b)')
    # rapidity coefficient of the soft function (soft_function.py, B1):
    #   coefficient of ln(delta1 delta2 b^2/b0^2) = 2 a_s Delta_C (mu^2/mu_b^2)^eps e^{-eps gamma_E} Gamma(-eps)
    ser = sp.series(sp.exp(-eps*sp.EulerGamma)*sp.gamma(-eps), eps, 0, 2).removeO()
    ser = sp.expand(sp.simplify(ser))
    print('   e^{-eps gamma_E} Gamma(-eps) =', ser, ' + O(eps^2)')
    assert sp.simplify(ser.coeff(eps, 0)) == 0
    print('   no eps^0 term: g1 = 0  =>  r1 = 0  (b0 = 2 e^{-gamma_E}; the same for the Fourier transform of')
    print('   1/l_perp^2 in 2 - 2 eps dimensions, which is what the rapidity divergence multiplies)')

if __name__ == '__main__':
    part1()
    part2()
    part3()
    print('all checks passed')
