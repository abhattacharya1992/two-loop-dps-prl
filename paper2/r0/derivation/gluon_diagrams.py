"""The tree splitting and the one-loop diagrams of J1(g -> gg), transcribed from paper2/r0/theory/gg_diagrams.py,
and their numerators after the residues in the minus components of the two loop momenta.

Momenta: p1 = (1, 0, 0) for the parent gluon, K = (x, K^-, k) and p1 - K for the daughters, l = (z, l^-, l_perp) for
the loop. Polarisations: e for the parent gluon, i and j for the daughters (e_1 and e_2 of the Supplement).
Each entry holds the numerator, the factor (rational number, power of the imaginary unit) and the kind of diagram.
The momentum fraction x is a rational number or the symbol X of algebra.py."""
from fractions import Fraction as Fr
from algebra import (C, V_, vec, vadd, vneg, mdot, T, V3, V4_parts, trace, contract, NBVEC, new_index, RN, Poly, coef)


def build(x):
    x = coef(x); xb = 1 - x
    u, w, z = V_('u'), V_('w'), V_('z')
    P1 = vec(1, 0); K = vec(x, u, {'k': 1}); L = vec(z, w, {'l': 1})
    E = vec(0, 0, {'e': 1}); I_ = vec(0, 0, {'i': 1}); J_ = vec(0, 0, {'j': 1})
    mom = lambda plus, minus, ck, cl: vec(plus, minus, {'k': ck, 'l': cl})
    O = lambda dirv, q, Kplus: vadd(dirv, NBVEC, coeffs=[1, -mdot(dirv, q)/Kplus])
    KmL = mom(x - z, u - w, 1, -1); QmK = mom(1 - x, -u, -1, 0); QmKpL = mom(1 - x + z, w - u, -1, 1)
    QmKmL = mom(1 - x - z, -u - w, -1, -1); QmL = mom(1 - z, -w, 0, -1); KpL = mom(x + z, u + w, 1, 1)
    ix = lambda *names: [new_index(n) for n in names]
    G = {}

    al, mu, nu = ix('al', 'mu', 'nu')
    G['tree'] = dict(num=T.V(E, al)*V3(P1, al, vneg(K), mu, vneg(QmK), nu)*T.V(O(I_, K, x), mu)*T.V(O(J_, QmK, 1 - x), nu),
                     factor=(Fr(1), 0), kind='tree')

    al, be, ga, de, mu, nu = ix('al', 'be', 'ga', 'de', 'mu', 'nu')
    G['Vtri'] = dict(num=(T.V(E, al)*V3(P1, al, vneg(KmL), be, vneg(QmKpL), ga)*V3(KmL, be, L, de, vneg(K), mu)
                          * V3(QmKpL, ga, vneg(L), de, vneg(QmK), nu)*T.V(O(I_, K, x), mu)*T.V(O(J_, QmK, 1 - x), nu)),
                     factor=(Fr(1), 1), kind='tri')

    def ghost(orient):
        al, mu, nu = ix('al', 'mu', 'nu')
        if orient == 1:
            qAB, qBC, qCA = KmL, vneg(L), vneg(QmKpL); num = T.V(qAB, al)*T.V(qBC, mu)*T.V(qCA, nu)
        else:
            qAC, qCB, qBA = QmKpL, L, vneg(KmL); num = T.V(qAC, al)*T.V(qBA, mu)*T.V(qCB, nu)
        return num*T.V(E, al)*T.V(O(I_, K, x), mu)*T.V(O(J_, QmK, 1 - x), nu)
    G['Vgh1'] = dict(num=ghost(1), factor=(Fr(1), 3), kind='tri')        # -i = i^3
    G['Vgh2'] = dict(num=ghost(2), factor=(Fr(1), 3), kind='tri')

    def quark(orient):
        al, mu, nu = ix('al', 'mu', 'nu')
        if orient == 1:
            qAB, qBC, qCA = KmL, vneg(L), vneg(QmKpL)
            tr = trace([('idx', al), ('vec', qCA), ('idx', nu), ('vec', qBC), ('idx', mu), ('vec', qAB)])
        else:
            qAC, qCB, qBA = QmKpL, L, vneg(KmL)
            tr = trace([('idx', al), ('vec', qBA), ('idx', mu), ('vec', qCB), ('idx', nu), ('vec', qAC)])
        return tr*T.V(E, al)*T.V(O(I_, K, x), mu)*T.V(O(J_, QmK, 1 - x), nu)
    G['Vq1'] = dict(num=quark(1), factor=(Fr(1), 0), kind='tri')
    G['Vq2'] = dict(num=quark(2), factor=(Fr(1), 0), kind='tri')

    def bub(which, r):
        al, mu, nu, be, ga = ix('al', 'mu', 'nu', 'be', 'ga')
        if which == 'a':
            lor = V4_parts(al, nu, be, ga)[r]
            return T.V(E, al)*lor*V3(L, be, KmL, ga, vneg(K), mu)*T.V(O(I_, K, x), mu)*T.V(O(J_, QmK, 1 - x), nu)
        lor = V4_parts(al, mu, be, ga)[r]
        return T.V(E, al)*lor*V3(L, be, QmKmL, ga, vneg(QmK), nu)*T.V(O(I_, K, x), mu)*T.V(O(J_, QmK, 1 - x), nu)
    for r in range(3):
        G[f'Vbuba{r+1}'] = dict(num=bub('a', r), factor=(Fr(1, 2), 1), kind='buba')
        G[f'Vbubb{r+1}'] = dict(num=bub('b', r), factor=(Fr(1, 2), 1), kind='bubb')

    def w1(which, r=None):
        al, rho, sig, nu, be = ix('al', 'rho', 'sig', 'nu', 'be')
        op1 = T.V(NBVEC, rho)*T.V(O(I_, KmL, x), sig); op2 = T.V(O(J_, QmK, 1 - x), nu)
        if which == 'g':
            num = T.V(E, al)*V3(P1, al, vneg(L), rho, vneg(QmL), be)*V3(QmL, be, vneg(KmL), sig, vneg(QmK), nu)
        elif which == 'kb':
            num = T.V(E, al)*V3(P1, al, vneg(KmL), sig, vneg(QmKpL), be)*V3(QmKpL, be, vneg(L), rho, vneg(QmK), nu)
        elif which == 'k':
            num = T.V(E, al)*V3(P1, al, vneg(K), be, vneg(QmK), nu)*V3(K, be, vneg(L), rho, vneg(KmL), sig)
        else:
            num = T.V(E, al)*V4_parts(al, rho, sig, nu)[r]
        return num*op1*op2
    G['W1g'] = dict(num=w1('g'), factor=(Fr(1), 3), kind='w1g')
    G['W1kb'] = dict(num=w1('kb'), factor=(Fr(1), 3), kind='w1kb')
    G['W1k'] = dict(num=w1('k'), factor=(Fr(1), 3), kind='w1k')
    for r in range(3):
        G[f'W1v4{r+1}'] = dict(num=w1('4', r), factor=(Fr(1), 3), kind='w1v4')

    def w2(which, r=None):
        al, rho, sig, mu, be = ix('al', 'rho', 'sig', 'mu', 'be')
        op2 = T.V(NBVEC, rho)*T.V(O(J_, QmKmL, 1 - x), sig); op1 = T.V(O(I_, K, x), mu)
        if which == 'g':
            num = T.V(E, al)*V3(P1, al, vneg(L), rho, vneg(QmL), be)*V3(QmL, be, vneg(K), mu, vneg(QmKmL), sig)
        elif which == 'k':
            num = T.V(E, al)*V3(P1, al, vneg(KpL), be, vneg(QmKmL), sig)*V3(KpL, be, vneg(K), mu, vneg(L), rho)
        elif which == 'kb':
            num = T.V(E, al)*V3(P1, al, vneg(K), mu, vneg(QmK), be)*V3(QmK, be, vneg(L), rho, vneg(QmKmL), sig)
        else:
            num = T.V(E, al)*V4_parts(al, mu, rho, sig)[r]
        return num*op1*op2
    G['W2g'] = dict(num=w2('g'), factor=(Fr(1), 3), kind='w2g')
    G['W2k'] = dict(num=w2('k'), factor=(Fr(1), 3), kind='w2k')
    G['W2kb'] = dict(num=w2('kb'), factor=(Fr(1), 3), kind='w2kb')
    for r in range(3):
        G[f'W2v4{r+1}'] = dict(num=w2('4', r), factor=(Fr(1), 3), kind='w2v4')
    return G


def sub_w(poly, expr_times_ze_inverse):
    """w -> expr/ze  (expr free of w)"""
    out = C(0)
    for (p,), c in poly.collect(['w']).items():
        out = out + (c*expr_times_ze_inverse**p).shift('ze', -p)
    return out


def to_r(poly):
    """transverse loop momentum l -> k - r; the scalar products of r take over the names of those of l"""
    kk = V_('kk')
    poly = poly.subs('ll', kk - 2*V_('pk') + V_('p2')).subs('lk', kk - V_('pk'))
    poly = poly.subs('el', V_('ek') - V_('pe')).subs('il', V_('ik') - V_('pi')).subs('jl', V_('jk') - V_('pj'))
    for a, b in (('p2', 'll'), ('pk', 'lk'), ('pe', 'el'), ('pi', 'il'), ('pj', 'jl')):
        poly = poly.subs(a, V_(b))
    return poly


def on_pole(num, on_shell, x):
    """numerator on the two poles in the minus components, as a polynomial in the transverse scalar products and a
    Laurent polynomial in ze (= zeta).
    on_shell = 'l':  z = x zeta, K^- = -k^2/(1-x), l^- = l_perp^2/(x zeta)
    on_shell = 'r':  r = K - l with r^+ = x zeta: z = x (1 - zeta), K^- = -k^2/(1-x), r^- = r_perp^2/(x zeta), l_perp = k - r_perp"""
    x = coef(x); xb = 1 - x; kk, ll, ze = V_('kk'), V_('ll'), V_('ze')
    if on_shell == 'l':
        n = num.subs('z', x*ze).subs('u', -kk/xb)
        return sub_w(n, ll/x)
    if on_shell == 'r':
        n = to_r(num.subs('z', x*(1 - ze)).subs('u', -kk/xb))
        # l^- = K^- - r^- = -k^2/(1-x) - r_perp^2/(x zeta)
        out = C(0)
        for (p,), c in n.collect(['w']).items():
            term = C(1)
            for _ in range(p):
                term = term*(-kk/xb) + (term*(-ll/x)).shift('ze', -1)
            out = out + c*term
        return out
    raise ValueError(on_shell)
