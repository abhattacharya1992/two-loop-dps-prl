#!/usr/bin/env python3
"""Conditional three-loop DPS logarithms from the renormalisation group (Paper II, Sec. 7.2).

The kernel of eq. (5.3) is evolved from mu_b to mu to one order beyond the two-loop structure tested in the
paper.  With X = ln(Q_*^2 e^{-i pi} b^2/b_0^2), X_mu = ln(mu^2 b^2/b_0^2), a = a_s(mu) = alpha_s(mu)/(4 pi),

    d K~/d ln mu = gamma_K K~,   gamma_K(a, mu) = -Dc gamma_cusp(a) (X - X_mu) + 2 DG(a),
    gamma_cusp(a) = 4 a + g1 a^2 + ...  (per unit Casimir: g1 = 4[(67/9 - 2 zeta2) C_A - (20/9) T_F n_f]),
    DG(a) = DG0 a + DG1 a^2 + ...        (non-cusp: gamma_a - gamma_R1 - gamma_R2),
    da/d ln mu^2 = -beta0 a^2 - beta1 a^3,
    K~(mu_b) = a(mu_b) K0 [1 + a(mu_b) r0 + a(mu_b)^2 (r0_3 + r1_3 X)]        (the boundary condition at mu_b:
        r1 = r2 = 0 at two loops; at three loops the boundary ansatz assumes r2_3 = 0 and no higher powers of X,
        with r0_3 and r1_3 open; the one-loop ultrasoft argument does not prove this higher-order assumption),

so that K~(mu)/(a K0) = [a(mu_b)/a] [1 + ...] exp( (1/2) int_0^{X_mu} gamma_K(a(X'), X') dX' ).  Expanded to a^2
relative to the lowest order, the a^1 term must be Dc X_mu^2 - 2 Dc X X_mu + B X_mu with B = beta0 + DG0 (eq. (5.3)),
and the a^2 term is new.  The amplitude is H1 H2 K~ with H_i = H_i^(0) (1 + a h + a^2 h2); the Fourier transforms of
X^n are those of eq. (5.5): for nu = 0, X^n -> -C_1 F_n(G) with F_0 = -G, F_1 = G^2, F_2 = -(4/3) G^3 - (2/3) zeta3,
and for nu = 2, X^n -> K R_n(G) with R_0 = 1, R_1 = 1 - 2 G, R_2 = 4 G^2 - 4 G + 2.  At a general scale,
X_mu = X + L_*, L_* = ln[mu^2/(Q_*^2 e^{-i pi})], as in eq. (7.1).

Conditional result: in units of C_1 (nu = 0) and K (nu = 2), the three-loop amplitude is a polynomial in G whose coefficients of
G^5, G^4, G^3 (nu = 0) and G^4, G^3, G^2 (nu = 2) depend only on Dc, beta0, B, g1, r0, h and L_*, all known; the
next coefficient contains DG1, beta1, h2 and the boundary term r1_3, which this analysis does not fix. An independent
r2_3 X^2 term would shift the G^3 coefficient by 4 r2_3/3 (nu = 0) and the G^2 coefficient by 4 r2_3 (nu = 2).
The checks below verify the algebra within the boundary ansatz, not the absence of such a term. The script
prints the coefficients for gluon splittings (gg -> HH: Dc = C_A, B = 2 beta0, h = C_A(-L^2 + zeta2) with
L = ln(m^2/mu^2) - i pi, r0 of eq. (7.3)) and for quark splittings (gg -> ZZ: Dc = 2 C_F - C_A, B = 6 C_F,
h = C_F(-L^2 + 3L - 8 + zeta2), r0 of eq. (7.4)), checks the two-loop terms against eqs. (5.3), (5.6) and (7.1) and
the conformal limit against the exponentiated form of the planar N = 4 section, and writes three_loop.json.

    python3 three_loop.py          (a few seconds; sympy)
"""
import json, pathlib
import sympy as sp

HERE = pathlib.Path(__file__).resolve().parent
a, X, Xmu, Ls, G, al = sp.symbols('a X X_mu L_* G alpha')
Dc, B, b0, b1, g1, DG0, DG1, r0, r0_3, r1_3, h, h2 = sp.symbols('Delta_C B beta0 beta1 g1 DG0 DG1 r0 r0_3 r1_3 h h2')
CA, CF, nf, x = sp.symbols('C_A C_F n_f x')
z2, z3 = sp.zeta(2), sp.Symbol('zeta3')
TF = sp.Rational(1, 2)
NA = 2                                                   # orders in a beyond the lowest-order kernel: three loops


def running(t):
    """a(mu') in terms of a = a(mu), t = ln(mu'^2/mu^2), through a^3"""
    return a - (b0*a**2 + b1*a**3)*t + b0**2*a**3*t**2


def kernel_bracket():
    """K~(mu)/(a K0) as a polynomial in X and X_mu through a^2 (three loops)"""
    Xp = sp.Symbol('Xp')
    at = running(Xp - Xmu)                                                    # a at the scale X'
    gamma_cusp = 4*at + g1*at**2
    DG = DG0*at + DG1*at**2
    integrand = -Dc*gamma_cusp*(X - Xp) + 2*DG
    E = sp.Rational(1, 2)*sp.integrate(sp.expand(integrand), (Xp, 0, Xmu))
    ab = running(-Xmu)                                                        # a(mu_b)
    boundary = ab*(1 + ab*r0 + ab**2*(r0_3 + r1_3*X))/a
    full = sp.series(boundary*sp.exp(E), a, 0, NA + 1).removeO()
    return sp.expand(full)


def transforms(nmax):
    """F_n (nu = 0, units of -C_1) and R_n (nu = 2, units of K) of eq. (5.5), n = 0..nmax"""
    expo = sp.exp(-2*G*al - sp.Rational(2, 3)*z3*al**3 - sp.Rational(2, 5)*sp.Symbol('zeta5')*al**5)
    gen0 = sp.series(expo/(2*al) - 1/(2*al), al, 0, nmax + 1).removeO()
    gen2 = sp.series(expo/(1 - al), al, 0, nmax + 1).removeO()
    F = [sp.expand(sp.factorial(n)*gen0.coeff(al, n)) for n in range(nmax + 1)]
    R = [sp.expand(sp.factorial(n)*gen2.coeff(al, n)) for n in range(nmax + 1)]
    assert F[0] == -G and sp.expand(F[1] - G**2) == 0 and sp.expand(F[2] + sp.Rational(4, 3)*G**3 + sp.Rational(2, 3)*z3) == 0
    assert R[0] == 1 and sp.expand(R[1] - (1 - 2*G)) == 0 and sp.expand(R[2] - (4*G**2 - 4*G + 2)) == 0
    return F, R


def amplitude(bracket, F, R):
    """the DPS term H1 H2 K~ transformed: nu = 0 in units of C_1, nu = 2 in units of K, through a^2 beyond the
    one-loop term; X_mu = X + L_*"""
    br = sp.expand(bracket.subs(Xmu, X + Ls))
    hard = sp.series((1 + a*h + a**2*h2)**2, a, 0, NA + 1).removeO()
    p = sp.Poly(sp.expand(br*hard), X)
    nu0 = sp.expand(sum(-c*F[m[0]] for m, c in zip(p.monoms(), p.coeffs())))     # X^n -> -C_1 F_n: units of C_1
    nu2 = sp.expand(sum(c*R[m[0]] for m, c in zip(p.monoms(), p.coeffs())))
    nu0 = sp.expand(sp.series(nu0, a, 0, NA + 1).removeO())
    nu2 = sp.expand(sp.series(nu2, a, 0, NA + 1).removeO())
    return nu0, nu2


def coefficients(expr, order):
    """the G-polynomial at a^order as {power: coefficient}"""
    e = sp.expand(expr.coeff(a, order))
    p = sp.Poly(e, G)
    return {m[0]: sp.factor(c) for m, c in zip(p.monoms(), p.coeffs())}


def main():
    br = kernel_bracket()
    # ---- checks at two loops: eq. (5.3) with B = beta0 + DG0
    two = sp.expand(br.coeff(a, 1))
    assert sp.expand(two - (Dc*Xmu**2 - 2*Dc*X*Xmu + (b0 + DG0)*Xmu + r0)) == 0, two
    # ---- the conformal limit (no running, no one-loop non-cusp term): exp(-Dc gamma_cusp X^2/4 + DG1 a^2 X) at
    #      X_mu = X, times the boundary constants, as the exponentiated form of the planar N = 4 section
    conf = {b0: 0, b1: 0, DG0: 0, r1_3: 0}
    expo = sp.series(sp.exp(-Dc*(4*a + g1*a**2)*X**2/4 + DG1*a**2*X)*(1 + a*r0 + a**2*r0_3), a, 0, NA + 1).removeO()
    assert sp.expand(br.subs(conf).subs(Xmu, X) - sp.expand(expo)) == 0
    F, R = transforms(2*NA + 1)
    nu0, nu2 = amplitude(br, F, R)
    # ---- two-loop amplitude against eqs. (5.6) and (7.1)
    c = coefficients(nu0, 1)
    assert sp.expand(c[3] + sp.Rational(4, 3)*Dc) == 0 and sp.expand(c[2] + (b0 + DG0)) == 0
    assert sp.expand(c[1] - (2*h + Dc*Ls**2 + (b0 + DG0)*Ls + r0)) == 0, c[1]
    d = coefficients(nu2, 1)
    assert sp.expand(d[2] + 4*Dc) == 0 and sp.expand(d[1] - (4*Dc - 2*(b0 + DG0))) == 0
    assert sp.expand(d[0] - ((b0 + DG0) - 2*Dc + 2*h + Dc*Ls**2 + (b0 + DG0)*Ls + r0)) == 0, d[0]
    # ---- three loops
    c3 = coefficients(nu0.subs(DG0, B - b0), 2)
    d3 = coefficients(nu2.subs(DG0, B - b0), 2)
    open_syms = {DG1, b1, h2, r1_3, r0_3}
    out = {'nu=0': {}, 'nu=2': {}}
    print('three-loop coefficients, nu = 0 (units of C_1), general Delta_C, B = beta0 + DG0, L_* = ln[mu^2/(Q_*^2 e^{-i pi})]:')
    for k in sorted(c3, reverse=True):
        free = sorted(str(s) for s in c3[k].free_symbols & open_syms)
        tag = '  [open: ' + ', '.join(free) + ']' if free else ''
        print(f'  G^{k}: {sp.expand(c3[k])}{tag}')
        out['nu=0'][f'G^{k}'] = {'coefficient': str(sp.expand(c3[k])), 'open': free}
    print('nu = 2 (units of K):')
    for k in sorted(d3, reverse=True):
        free = sorted(str(s) for s in d3[k].free_symbols & open_syms)
        tag = '  [open: ' + ', '.join(free) + ']' if free else ''
        print(f'  G^{k}: {sp.expand(d3[k])}{tag}')
        out['nu=2'][f'G^{k}'] = {'coefficient': str(sp.expand(d3[k])), 'open': free}
    # the predicted coefficients contain no open quantity
    assert all(not (c3[k].free_symbols & open_syms) for k in (5, 4, 3))
    assert all(not (d3[k].free_symbols & open_syms) for k in (4, 3, 2))
    assert (c3[2].free_symbols & open_syms) and (d3[1].free_symbols & open_syms)
    # ---- the two splittings at mu^2 = m^2 (L = -i pi in h), expressed through C_A, n_f, C_F, x
    beta0 = sp.Rational(11, 3)*CA - sp.Rational(4, 3)*TF*nf
    g1val = 4*((sp.Rational(67, 9) - 2*z2)*CA - sp.Rational(20, 9)*TF*nf)
    lx, ly = sp.log(x), sp.log(1 - x)
    cases = {
        'gluon splittings (gg -> HH)': {
            Dc: CA, B: 2*beta0, b0: beta0, g1: g1val, h: CA*(-(-sp.I*sp.pi)**2 + z2),
            'r0': {0: CA*(sp.Rational(134, 9) - 5*z2 - 4*lx*ly) - sp.Rational(40, 9)*TF*nf
                      + sp.Rational(2, 3)*(CA - 2*TF*nf)/(1/(x*(1 - x)) + 2*x*(1 - x)),
                   2: CA*(sp.Rational(134, 9) - 5*z2 - 4*lx*ly) - sp.Rational(40, 9)*TF*nf
                      + sp.Rational(2, 3)*(CA - 2*TF*nf)/((x**4 + (1 - x)**4)/(x*(1 - x)))},
            'L_*': -sp.log(x*(1 - x)) + sp.I*sp.pi},                        # Q_*^2 = m^4/s, m^2 = x(1-x) s, mu = m
        'quark splittings (gg -> ZZ)': {
            Dc: 2*CF - CA, B: 6*CF, b0: beta0, g1: g1val, h: CF*(-(-sp.I*sp.pi)**2 + 3*(-sp.I*sp.pi) - 8 + z2),
            'r0': {0: CF*(12 - 6*z2 + 2*(lx - ly)**2) + CA*(2 + z2 - 2*(lx**2 + ly**2)),
                   2: CF*(12 - 6*z2 + 2*(lx - ly)**2) + CA*(2 + z2 - 2*(lx**2 + ly**2))},
            'L_*': sp.Symbol('L_*')},                                        # Q_* of eq. (5.2) for quarks: kept symbolic
    }
    out['cases'] = {}
    for name, sub in cases.items():
        subs = {k: v for k, v in sub.items() if isinstance(k, sp.Symbol)}
        rec = {}
        print(name + ':')
        for nu, cc, powers in (('nu=0', c3, (5, 4, 3)), ('nu=2', d3, (4, 3, 2))):
            for k in powers:
                e = sp.expand(cc[k].subs(subs).subs(r0, sub['r0'][int(nu[-1])]).subs(Ls, sub['L_*']))
                e = sp.collect(sp.expand(e), [CA, CF, nf])
                rec[f'{nu} G^{k}'] = str(e)
                print(f'  {nu} G^{k}: {e}')
        out['cases'][name] = rec
    # ---- numbers for gg -> HH at mu = m, x = 5/8, 3/4, 5/6, by colour structure (N^2, N n_f, n_f^2), real and imaginary parts
    sub = cases['gluon splittings (gg -> HH)']
    subs = {k: v for k, v in sub.items() if isinstance(k, sp.Symbol)}
    out['numbers gg -> HH at mu = m'] = {}
    for xv in (sp.Rational(5, 8), sp.Rational(3, 4), sp.Rational(5, 6)):
        rec = {}
        for nu, cc, powers in (('nu=0', c3, (5, 4, 3)), ('nu=2', d3, (4, 3, 2))):
            for k in powers:
                e = sp.expand(cc[k].subs(subs).subs(r0, sub['r0'][int(nu[-1])]).subs(Ls, sub['L_*']).subs(x, xv))
                e = sp.expand(e.subs(z3, sp.zeta(3)))
                parts = {}
                for name_c, mono in (('N^2', CA**2), ('N n_f', CA*nf), ('n_f^2', nf**2)):
                    coeff = sp.Poly(e, CA, nf).coeff_monomial(mono)
                    val = complex(sp.N(coeff, 20))
                    parts[name_c] = [val.real, val.imag]
                assert sp.expand(e - sum(sp.Poly(e, CA, nf).coeff_monomial(m)*m for m in (CA**2, CA*nf, nf**2))) == 0
                rec[f'{nu} G^{k}'] = parts
        out['numbers gg -> HH at mu = m'][str(xv)] = rec
    (HERE/'three_loop.json').write_text(json.dumps(out, indent=1))
    print('wrote three_loop.json')
    print('passed')


if __name__ == '__main__':
    main()
