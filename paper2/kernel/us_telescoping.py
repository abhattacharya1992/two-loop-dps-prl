#!/usr/bin/env python3
"""Stage B2: the ultrasoft gluons do not resolve the DPS splitting.

A gluon l enters the p1 system {p1, k, kbar} of the DPS loop; K = (x, K^-, K_perp) is the momentum that
reaches vertex 1 from that system (fixed by the p2-side residue), units p1^+ = 1.  The loop's minus
component u = k^- is integrated at leading power in the DPS loop (only the propagators of the p1 system
depend on it).  Three attachments:

  (a) on k (absorbed just before vertex 1):  1/(D1 D2 D3),  k = K - l,
        D1 = k^2,  D2 = (k+l)^2 = K^2,  D3 = (p1-k)^2;
  (b) on kbar (absorbed just before vertex 2):  1/(D1 D3 D3''),  k = K,  D3'' = (p1 - K + l)^2;
  (c) on the parent before the splitting:  1/(P D1 D3'),  P = (p1+l)^2,  D3' = (p1 + l - K)^2;
  tree:  1/(D1 D3) with k = K:  int du = 2 pi i/K_perp^2.

Gauge-theory eikonal numerators nbar.(2p + l) = 2 p^+ + l^+ give the weights w_a = x - l^+/2,
w_b = 1 - x + l^+/2, w_c = 1 + l^+/2 (an overall 2 dropped).  Claim, for l ~ lambda^2 in all components
(ultrasoft) and K_perp ~ lambda:

  w_a (a)/tree -> 1/(l^- - sigma + i0),   w_b (b)/tree -> 1/(l^- - sigma + i0),
  w_c (c)/tree -> 1/(l^- + i0) - 1/(l^- - sigma + i0),       sigma = K_perp^2/(x(1-x)),

so with colour conservation at the splitting, T_k + T_kbar = T_p1, the sum over attachments is
T_p1 tree/(l^- + i0) (checked component by component: the T_k part w_a(a) + w_c(c) and the T_kbar part
w_b(b) + w_c(c) each tend to tree/l^-): the eikonal line of the on-shell parent.  The same holds on the p2 side, and the
ultrasoft region of the one-loop kernel, summed over attachments, is T_p1.T_p2 times a scaleless integral,
i.e. zero in dimensional regularisation, although each attachment alone is leading power.  The script
checks (1) the residue formulas against direct quadrature in u with a finite i0, (2) the ultrasoft limits
and the cancellation, with corrections O(lambda), (3) that the cancellation needs the ultrasoft scaling
(Glauber l_perp ~ lambda: no cancellation; soft l ~ lambda: the parent attachment is suppressed by one
power of lambda), (4) the Glauber region: the l^- integral of the parent attachment vanishes (both poles on
the same side, 1/l^-^2 decay), and that of a daughter attachment, whose tails fall like 1/l^- at both ends,
tends with a cutoff symmetric in l_z to -i pi times the coefficient of the tail, independent of sigma (no
logarithm of the hard scale), and (5) the bookkeeping of the 3 x 3 attachments for the ultrasoft integral
of one pair: given that an integral with a zero off-shellness is scaleless, the sum vanishes identically,
while the daughters alone give T_p1.T_p2 J (the content is in (1) and (2); (5) only adds up the colours).
    python3 us_telescoping.py
"""
import mpmath as mp
mp.mp.dps = 30
I = mp.mpc(0, 1)
PI = mp.pi

def props(l, K, x, i0):
    """the propagators as functions of u (the loop minus component); l = (l+, l-, (lx, ly)), K = (Kx, Ky)"""
    lp, lm, (lx, ly) = l
    Kx, Ky = K
    K2 = Kx**2 + Ky**2
    kx, ky = Kx - lx, Ky - ly                       # k_perp in graph (a)
    k2 = kx**2 + ky**2
    kp = x - lp                                     # k^+ in graph (a)
    lK2 = (lx - Kx)**2 + (ly - Ky)**2
    l2 = lx**2 + ly**2
    a = lambda u: 1/((kp*u - k2 + I*i0)*(x*(u + lm) - K2 + I*i0)*(-(1 - kp)*u - k2 + I*i0))
    b = lambda u: 1/((x*u - K2 + I*i0)*(-(1 - x)*u - K2 + I*i0)*((1 - x + lp)*(lm - u) - lK2 + I*i0))
    P = (1 + lp)*lm - l2 + I*i0
    c = lambda u: 1/(P*(x*u - K2 + I*i0)*((1 + lp - x)*(lm - u) - lK2 + I*i0))
    tree = lambda u: 1/((x*u - K2 + I*i0)*(-(1 - x)*u - K2 + I*i0))
    poles = {'a': [k2/kp, -k2/(1 - kp), -lm + K2/x], 'b': [K2/x, -K2/(1 - x), lm - lK2/(1 - x + lp)],
             'c': [K2/x, lm - lK2/(1 + lp - x)], 'tree': [K2/x, -K2/(1 - x)]}
    return {'a': a, 'b': b, 'c': c, 'tree': tree}, poles

def residues(l, K, x):
    """the same integrals by residues (i0 -> 0 in the remaining factors; the side of each pole from its i0)"""
    lp, lm, (lx, ly) = l
    Kx, Ky = K
    K2 = Kx**2 + Ky**2
    kx, ky = Kx - lx, Ky - ly
    k2 = kx**2 + ky**2
    kp = x - lp
    lK2 = (lx - Kx)**2 + (ly - Ky)**2
    l2 = lx**2 + ly**2
    # (a): close above, pole of D3 at u3 = -k2/(1-kp)
    u3 = -k2/(1 - kp)
    ra = 2*PI*I/(-(1 - kp))/((kp*u3 - k2)*(x*(u3 + lm) - K2 + I*mp.mpf('1e-60')))
    # (b), (c), tree: close below, pole of D1 at u1 = K2/x
    u1 = K2/x
    rb = -2*PI*I/x/((-(1 - x)*u1 - K2)*((1 - x + lp)*(lm - u1) - lK2 + I*mp.mpf('1e-60')))
    P = (1 + lp)*lm - l2 + I*mp.mpf('1e-60')
    rc = -2*PI*I/x/(P*((1 + lp - x)*(lm - u1) - lK2 + I*mp.mpf('1e-60')))
    rt = -2*PI*I/x/(-(1 - x)*u1 - K2)
    return {'a': ra, 'b': rb, 'c': rc, 'tree': rt}

def quad_check():
    print('(1) residue formulas against quadrature in u with finite i0 (linear extrapolation i0 -> 0)')
    x = mp.mpf('0.75'); lam = mp.mpf('0.05')
    K = (mp.mpf('0.6')*lam, mp.mpf('0.3')*lam)
    worst = 0
    for lhat in [(mp.mpf('0.7'), mp.mpf('-1.3'), (mp.mpf('0.4'), mp.mpf('-0.2'))),
                 (mp.mpf('-0.5'), mp.mpf('2.1'), (mp.mpf('-0.3'), mp.mpf('0.6')))]:
        l = (lhat[0]*lam**2, lhat[1]*lam**2, (lhat[2][0]*lam**2, lhat[2][1]*lam**2))
        res = residues(l, K, x)
        vals = {}
        for i0 in (lam**2*mp.mpf('1e-4'), lam**2*mp.mpf('5e-5')):
            f, poles = props(l, K, x, i0)
            for g in ('a', 'b', 'c', 'tree'):
                pts = sorted({-mp.inf, mp.inf, *poles[g], *[p + s*m*i0 for p in poles[g] for s in (1, -1) for m in (1, 10, 100)]})
                vals.setdefault(g, []).append(mp.quad(f[g], pts))
        for g in ('a', 'b', 'c', 'tree'):
            ex = 2*vals[g][1] - vals[g][0]
            rel = abs(ex - res[g])/abs(res[g])
            worst = max(worst, rel)
            print(f'    l/lam^2 = {[float(lhat[0]), float(lhat[1])]}  graph {g:4s}  residue {mp.nstr(res[g], 10)}  quadrature {mp.nstr(ex, 10)}  rel {mp.nstr(rel, 2)}')
    print(f'    worst relative difference {mp.nstr(worst, 2)}')
    return worst

def us_limit():
    print('(2) ultrasoft scaling l = lambda^2 lhat: each attachment against its limit, and the colour sum')
    x = mp.mpf('0.75')
    lhat = (mp.mpf('0.7'), mp.mpf('-1.3'), (mp.mpf('0.4'), mp.mpf('-0.2')))
    out = []
    for lam in (mp.mpf('0.1'), mp.mpf('0.01'), mp.mpf('0.001'), mp.mpf('0.0001')):
        K = (mp.mpf('0.6')*lam, mp.mpf('0.3')*lam)
        l = (lhat[0]*lam**2, lhat[1]*lam**2, (lhat[2][0]*lam**2, lhat[2][1]*lam**2))
        r = residues(l, K, x)
        lp, lm = l[0], l[1]
        sig = (K[0]**2 + K[1]**2)/(x*(1 - x))
        wa, wb, wc = x - lp/2, 1 - x + lp/2, 1 + lp/2
        da = abs(wa*r['a']/r['tree']*(lm - sig) - 1)
        db = abs(wb*r['b']/r['tree']*(lm - sig) - 1)
        dc = abs(wc*r['c']/r['tree']/(1/lm - 1/(lm - sig)) - 1)
        # the colour sum T_k w_a(a) + T_kbar w_b(b) + (T_k + T_kbar) w_c(c), component by component:
        # T_k: w_a(a) + w_c(c),  T_kbar: w_b(b) + w_c(c); each must tend to tree/l^-
        dk = abs((wa*r['a'] + wc*r['c'])/r['tree']*lm - 1)
        dkb = abs((wb*r['b'] + wc*r['c'])/r['tree']*lm - 1)
        dsum = max(dk, dkb)
        single = abs(wa*r['a']/r['tree'])*abs(lm)         # one attachment alone, in units of the parent eikonal
        out.append((lam, da, db, dc, dsum))
        print(f'    lambda = {mp.nstr(lam, 2):7s}  |w_a a/tree (l^- - sigma) - 1| = {mp.nstr(da, 2):8s}  (b) {mp.nstr(db, 2):8s}  '
              f'(c) {mp.nstr(dc, 2):8s}  colour sum: |T_k part x l^-/tree - 1| = {mp.nstr(dk, 2):8s} |T_kbar part ...| = {mp.nstr(dkb, 2):8s}'
              f'  [one attachment alone: {mp.nstr(single, 3)}]')
    return out

def other_regions():
    print('(3) the cancellation needs the ultrasoft scaling')
    x = mp.mpf('0.75')
    for name, scal in (('Glauber (lambda^2, lambda^2, lambda)', (2, 2, 1)), ('soft (lambda, lambda, lambda)', (1, 1, 1))):
        for lam in (mp.mpf('0.01'), mp.mpf('0.001')):
            K = (mp.mpf('0.6')*lam, mp.mpf('0.3')*lam)
            l = (mp.mpf('0.7')*lam**scal[0], mp.mpf('-1.3')*lam**scal[1], (mp.mpf('0.4')*lam**scal[2], mp.mpf('-0.2')*lam**scal[2]))
            r = residues(l, K, x)
            lp, lm = l[0], l[1]
            wa, wb, wc = x - lp/2, 1 - x + lp/2, 1 + lp/2
            dk = abs((wa*r['a'] + wc*r['c'])/r['tree']*lm - 1)
            dkb = abs((wb*r['b'] + wc*r['c'])/r['tree']*lm - 1)
            print(f'    {name:38s} lambda = {mp.nstr(lam, 2):6s}  colour sum: |T_k part x l^-/tree - 1| = {mp.nstr(dk, 3):9s}'
                  f' |T_kbar part ...| = {mp.nstr(dkb, 3):9s}  |parent/daughter (a)| = {mp.nstr(abs(wc*r["c"]/(wa*r["a"])), 3)}')
    print('(4a) Glauber region, the l^- integral of the parent attachment (both poles below the axis, 1/l^-^2 decay): zero')
    lam = mp.mpf('0.01'); K = (mp.mpf('0.6')*lam, mp.mpf('0.3')*lam); i0 = mp.mpf('1e-6')*lam**2
    lp, lt = mp.mpf('0.7')*lam**2, (mp.mpf('0.4')*lam, mp.mpf('-0.2')*lam)
    K2 = K[0]**2 + K[1]**2; lt2 = lt[0]**2 + lt[1]**2; lK2 = (lt[0] - K[0])**2 + (lt[1] - K[1])**2
    c_of = lambda t: -2*PI*I/x/(((1 + lp)*t - lt2 + I*i0)*((1 + lp - x)*(t - K2/x) - lK2 + I*i0))
    p1, p2 = lt2/(1 + lp), K2/x + lK2/(1 + lp - x)
    brk = sorted({-mp.inf, mp.inf, p1, p2, *[p + s*m*i0 for p in (p1, p2) for s in (1, -1) for m in (10, 100)]})
    v = mp.quad(c_of, brk); a = mp.quad(lambda t: abs(c_of(t)), brk)
    print(f'    |int dl^- (c)| / int dl^- |(c)| = {mp.nstr(abs(v)/a, 3)}')
    glauber_parent_ratio = abs(v)/a
    print('(4b) Glauber region, daughter attachment: 1/l^- tails at both ends; with a cutoff |l^-| < L symmetric in l_z the')
    print('     integral tends to -i pi times the coefficient of the tail, whatever sigma (here two values of K_perp):')
    devs = []
    for Kh in ((mp.mpf('0.6'), mp.mpf('0.3')), (mp.mpf('-1.7'), mp.mpf('0.9'))):
        K = (Kh[0]*lam, Kh[1]*lam); K2 = K[0]**2 + K[1]**2
        k2 = (K[0] - lt[0])**2 + (K[1] - lt[1])**2; kp = x - lp
        u3 = -k2/(1 - kp)
        a_of = lambda t: 2*PI*I/(-(1 - kp))/((kp*u3 - k2)*(x*(u3 + t) - K2 + I*i0))
        z = K2/x - u3
        tail = 2*PI*I/(-(1 - kp))/((kp*u3 - k2)*x)         # (a) ~ tail/l^- at large |l^-|
        for L in (mp.mpf(10)**4*lam**2, mp.mpf(10)**8*lam**2):
            val = mp.quad(a_of, [-L, z - 100*i0, z - i0, z, z + i0, z + 100*i0, L])
            devs.append(abs(val/(-I*PI*tail) - 1))
            print(f'    K_perp/lam = {[float(Kh[0]), float(Kh[1])]}, sigma-type pole at l^- = {mp.nstr(z/lam**2, 4)} lam^2, L = {mp.nstr(L/lam**2, 2)} lam^2:'
                  f'  int/(-i pi tail) - 1 = {mp.nstr(abs(val/(-I*PI*tail) - 1), 3)}')
    print('(4c) Glauber region, parent attachment at large l^- (fixed l_perp ~ lambda, l^+ ~ lambda^2): (c)/tree ~ 1/l^-^2, (a)/tree ~ 1/l^-')
    lam = mp.mpf('0.01'); K = (mp.mpf('0.6')*lam, mp.mpf('0.3')*lam)
    for m in (10, 100, 1000, 10000):
        l = (mp.mpf('0.7')*lam**2, m*lam**2, (mp.mpf('0.4')*lam, mp.mpf('-0.2')*lam))
        r = residues(l, K, x)
        print(f'    l^- = {m:6d} lambda^2:  l^-^2 (c)/tree = {mp.nstr(abs(l[1]**2*r["c"]/r["tree"]), 6):12s}  l^- (a)/tree = {mp.nstr(abs(l[1]*r["a"]/r["tree"]), 6)}')
    return glauber_parent_ratio, max(devs)

def us_integrals():
    print('(5) ultrasoft integral of one pair and of the sum over the 3 x 3 attachments')
    import sympy as sp
    eps, s1, s2 = sp.symbols('epsilon sigma1 sigma2', positive=True)
    Tk, Tkb, Tr, Trb = sp.symbols('T_k T_kbar T_r T_rbar')
    # J(s1, s2) = int dl^+ dl^- d^{2-2eps} l_perp / [(l^2 + i0)(l^- - s1 + i0)(l^+ + s2 - i0)]
    #           = -2 i pi^2 pi^{-eps} Gamma(eps)^2 Gamma(1-eps) (s1 s2 e^{i pi})^{-eps}  (lab_triangle.py, model numerics)
    J = lambda a, b: -2*sp.I*sp.pi**2*sp.pi**(-eps)*sp.gamma(eps)**2*sp.gamma(1 - eps)*(a*b)**(-eps)*sp.exp(-sp.I*sp.pi*eps)
    # the parent factors 1/l^- - 1/(l^- - s1) and 1/l^+ - 1/(l^+ + s2) give J(0, .) - J(s1, .); J with a zero argument
    # is scaleless (zero in dimensional regularisation)
    Jz = lambda a, b: 0 if (a == 0 or b == 0) else J(a, b)
    p1side = [(Tk, s1, 1), (Tkb, s1, 1), (Tk + Tkb, 0, 1), (Tk + Tkb, s1, -1)]   # (colour, sigma, sign); last two: parent
    p2side = [(Tr, s2, 1), (Trb, s2, 1), (Tr + Trb, 0, 1), (Tr + Trb, s2, -1)]
    tot = sum(ca*cb*sa*sb*Jz(a, b) for ca, a, sa in p1side for cb, b, sb in p2side)
    ser = sp.series(J(s1, s2), eps, 0, 1).removeO()
    print('    one pair (k, rbar):  T_k.T_rbar J,  J = (', sp.simplify(ser.coeff(eps, -2)), ')/eps^2 + (',
          sp.simplify(sp.expand_log(ser.coeff(eps, -1), force=True)), ')/eps + O(1)')
    print('    sum over attachments to {p1, k, kbar} x {p2, r, rbar}:', sp.simplify(tot))
    print('    daughters only (no parent attachments):', sp.factor(sp.simplify(sum(ca*cb*Jz(a, b) for ca, a, sa in p1side[:2] for cb, b, sb in p2side[:2]))))
    return sp.simplify(tot)

if __name__ == '__main__':
    w = quad_check()
    assert w < 1e-3, 'residue formulas'
    out = us_limit()
    assert out[-1][4] < 1e-3 and out[-1][4] < out[0][4], 'ultrasoft cancellation'
    gp, gd = other_regions()
    assert gp < 1e-6 and gd < 1e-2
    tot = us_integrals()
    assert tot == 0
    print('all checks passed')
