"""The one-loop DPS coefficients against the lowest-order 1 -> 2 splitting kernels of double parton
distributions (DPDs). Exact rational arithmetic, standard library only; run: python3 check_dpd_kernels.py

Supports Eqs. (eq:dpdkernels), (eq:dpdsquare) and (eq:Ksquared) of the Supplement and the sentence on
the DPD kernels in its section "Test in two-loop gg -> ZZ".

Kernels as printed in Diehl, Kasemets, Keane, JHEP 05 (2014) 118 [arXiv:1401.1233], Appendix B
(zb = 1 - z, z the momentum fraction of the first parton):
    T_{g->gg}        = 2 Nc [ zb/z + z/zb + z zb ]
    T_{g->Dg Dg}     = 2 Nc ( 2 - z zb )            two longitudinally polarised gluons
    T_{g->dg dg}     = 2 Nc z zb                    two linearly polarised gluons
    T_{g->g dg}      = 2 Nc z / zb                  first unpolarised, second linearly polarised
    T_{g->q qbar}    = 1/2 ( z^2 + zb^2 )
    T_{g->Dq Dqbar}  = - T_{g->q qbar}
    T_{g->dq dqbar}  = - z zb
Values quoted in Diehl, Gaunt, Schoenwald, JHEP 06 (2017) 083 [arXiv:1702.06486], Sec. 9:
    P_{g->gg}(1/2) = 27/2,  P_{g->q qbar}(1/2) = 1/4,  P_gg : P_DgDg : P_dgdg = 9 : 7 : 1 at 1/2.

Functions of the manuscript (Letter Eq. (eq:oneloop) and the text below it; Supplement Eqs. (eq:FA),
(eq:FS), (eq:K1)):
    F_S = F_0 = 1/(x(1-x)) + 2x(1-x)       gg -> HH, equal helicities
    F_A       = 1/(x(1-x)) - 2x(1-x)       gg -> AA, equal helicities
    W   = F_2 = (x^4 + (1-x)^4)/(x(1-x))   opposite helicities, both processes
    gg -> ZZ: C_1 proportional to x^2 (++++) or (1-x)^2 (++--), K proportional to x(1-x).

What is checked:
 (a) the kernels above follow from the tree splitting amplitudes in the helicity basis
     (probabilities and interference terms), including their normalisation;
 (b) the three values quoted by Diehl, Gaunt and Schoenwald;
 (c) 2 C_A F_S = T_gg + T_DgDg + 2 T_dgdg,  2 C_A F_A = T_gg + T_DgDg - 2 T_dgdg,
     2 C_A W = T_gg - T_DgDg;   x^2 + (1-x)^2 = 2 T_qq,  x(1-x) = -T_dqdq;
 (d) the helicity amplitudes A(lam1, lam2) built from the two splittings and two (pseudo)scalar
     vertices equal F_S or F_A (equal helicities) and W (opposite helicities);
 (e) (1/4) sum_helicities |A|^2 = [T T' + DT DT' + 2 dT dT' +- 2 T_gdg T_gdg'] / (2 C_A)^2;
 (f) unequal masses, quark channel: the products x1 (1 - x2) and (1 - x1) x2 of the one-loop
     coefficients in the Supplement's section on different invariant masses;
 (g) |1 + a_s(R + i I)|^2 = 1 + 2 a_s R + a_s^2 (R^2 + I^2): the phase 2 pi Delta_C X_mu enters the
     squared kernel only at relative order a_s^2;
 (h) control against a published subtraction term: for two scalars coupled to massless quarks, the
     same helicity sum gives (1 + beta^4)/4 = 2 T_qq^2 + 2 T_dqdq^2 with x = (1 + beta)/2, the
     dependence on beta of Eq. (9.22) of Diehl, Gaunt and Schoenwald beyond its flux factor, with
     their weights (unpolarised 1, longitudinal 1, transverse delta^{jj'} delta^{jj'} = 2);
 (i) vector vertices (gg -> ZZ): the helicity average of the squared one-loop DPS coefficients,
     summed over the helicities of the two bosons, is (1/2)(x^2 + (1-x)^2)^2 = 2 T_qq^2; the
     transversely polarised kernel does not contribute.

Not checked here: the absolute normalisation of the cross section, the colour channels of the DPDs,
real emission, and any relation of r_0 to the two-loop DPD splitting kernels.
"""
from fractions import Fraction as Fr
import math

Nc = 3
CA = Fr(Nc)
TF = Fr(1, 2)


# ---------- published kernels ---------------------------------------------------------------------
def T_gg(z):
    zb = 1 - z
    return 2*Nc*(zb/z + z/zb + z*zb)


def T_DgDg(z):
    return 2*Nc*(2 - z*(1 - z))


def T_dgdg(z):
    return 2*Nc*z*(1 - z)


def T_gdg(z):
    return 2*Nc*z/(1 - z)


def T_qq(z):
    return Fr(1, 2)*(z*z + (1 - z)**2)


def T_DqDq(z):
    return -T_qq(z)


def T_dqdq(z):
    return -z*(1 - z)


# ---------- functions of the manuscript -------------------------------------------------------------
def F_S(x):
    return 1/(x*(1 - x)) + 2*x*(1 - x)


def F_A(x):
    return 1/(x*(1 - x)) - 2*x*(1 - x)


def W(x):
    return (x**4 + (1 - x)**4)/(x*(1 - x))


# ---------- tree splitting amplitudes, helicity basis -----------------------------------------------
def n_gg(lam, h1, h2, z):
    """modulus of the g -> g g amplitude times sqrt(z(1-z)) |k_T|; z is the fraction of gluon 1"""
    if lam == -1:
        return n_gg(1, -h1, -h2, z)                      # parity
    return {(1, 1): Fr(1), (1, -1): z*z, (-1, 1): (1 - z)**2, (-1, -1): Fr(0)}[(h1, h2)]


def n_qq(lam, h, z):
    """modulus of the g -> q qbar amplitude (common factor removed); quark helicity h, fraction z"""
    return z if h == lam else 1 - z


def main():
    out = []

    def check(name, cond):
        out.append((name, bool(cond)))

    xs = [Fr(1, 2), Fr(1, 3), Fr(3, 5), Fr(5, 8), Fr(3, 4), Fr(5, 6), Fr(9, 10), Fr(7, 13), Fr(2, 11)]

    # (a)
    for z in xs:
        d = z*(1 - z)
        prob = {(h1, h2): n_gg(1, h1, h2, z)**2/d for h1 in (1, -1) for h2 in (1, -1)}
        unpol = CA*sum(prob.values())
        heli = CA*sum(h1*h2*p for (h1, h2), p in prob.items())
        lin = 2*CA*n_gg(1, 1, -1, z)*n_gg(1, -1, 1, z)/d
        mixed = 2*CA*sum(n_gg(1, h1, 1, z)*n_gg(1, h1, -1, z) for h1 in (1, -1))/d
        check('(a) gluon kernels from the helicity amplitudes',
              (unpol, heli, lin, mixed) == (T_gg(z), T_DgDg(z), T_dgdg(z), T_gdg(z)))
        pq = {h: n_qq(1, h, z)**2 for h in (1, -1)}
        check('(a) quark kernels from the helicity amplitudes',
              (TF*sum(pq.values()), -TF*sum(pq.values()), -2*TF*n_qq(1, 1, z)*n_qq(1, -1, z))
              == (T_qq(z), T_DqDq(z), T_dqdq(z)))

    # (b)
    h = Fr(1, 2)
    check('(b) P_gg(1/2) = 27/2 and P_qq(1/2) = 1/4', (T_gg(h), T_qq(h)) == (Fr(27, 2), Fr(1, 4)))
    check('(b) P_gg : P_DgDg : P_dgdg = 9 : 7 : 1 at 1/2', (T_gg(h)/T_dgdg(h), T_DgDg(h)/T_dgdg(h)) == (9, 7))

    # (c)
    for x in xs:
        check('(c) 2 C_A F_S = T_gg + T_DgDg + 2 T_dgdg', 2*CA*F_S(x) == T_gg(x) + T_DgDg(x) + 2*T_dgdg(x))
        check('(c) 2 C_A F_A = T_gg + T_DgDg - 2 T_dgdg', 2*CA*F_A(x) == T_gg(x) + T_DgDg(x) - 2*T_dgdg(x))
        check('(c) 2 C_A W = T_gg - T_DgDg', 2*CA*W(x) == T_gg(x) - T_DgDg(x))
        check('(c) x^2 + (1-x)^2 = 2 T_qq and x(1-x) = -T_dqdq',
              (x*x + (1 - x)**2, x*(1 - x)) == (2*T_qq(x), -T_dqdq(x)))

    # (d), (e): gluon 1 -> k (helicity h, fraction x) + kbar (h', 1-x); gluon 2 -> r (h, 1-x) + rbar (h', x);
    # (k, r) -> X_1 and (kbar, rbar) -> X_2 need equal helicities; vertex factors m_+ = 1, m_- = parity
    for x in xs:
        d = x*(1 - x)
        for parity, F_eq in ((1, F_S), (-1, F_A)):
            m = {1: 1, -1: parity}

            def A(l1, l2):
                return sum(n_gg(l1, a, b, x)*n_gg(l2, a, b, 1 - x)*m[a]*m[b] for a in (1, -1) for b in (1, -1))/d

            check('(d) A(++) = F_S or F_A, A(+-) = W up to the sign of the vertex',
                  (A(1, 1), A(1, -1), A(-1, -1), A(-1, 1)) == (F_eq(x), parity*W(x), F_eq(x), parity*W(x)))
            lhs = sum(A(l1, l2)**2 for l1 in (1, -1) for l2 in (1, -1))/4
            rhs = (T_gg(x)*T_gg(1 - x) + T_DgDg(x)*T_DgDg(1 - x) + 2*T_dgdg(x)*T_dgdg(1 - x)
                   + parity*2*T_gdg(x)*T_gdg(1 - x))/(2*CA)**2
            check('(e) helicity average of |A|^2 = sum of products of kernels',
                  lhs == rhs == (F_eq(x)**2 + W(x)**2)/2)

    # (f)
    x1, x2 = Fr(2, 3), Fr(1, 5)
    check('(f) unequal masses: x1 (1 - x2) and (1 - x1) x2',
          (n_qq(1, 1, x1)*n_qq(1, 1, 1 - x2), n_qq(1, -1, x1)*n_qq(1, -1, 1 - x2)) == (x1*(1 - x2), (1 - x1)*x2))

    # (g)
    DC, B, r0, a_s = 3.0, 2*(11 - 10/3), 2.817, 0.0094
    Xr, Xmu = -3.0, 0.7
    bracket = DC*Xmu**2 - 2*DC*complex(Xr, -math.pi)*Xmu + B*Xmu + r0
    exact = abs(1 + a_s*bracket)**2 - 1
    first = 2*a_s*(DC*Xmu**2 - 2*DC*Xr*Xmu + B*Xmu + r0)
    check('(g) the phase enters |K|^2 only at relative order a_s^2',
          abs(exact - first - a_s**2*abs(bracket)**2) < 1e-14 and abs(bracket.imag - 2*math.pi*DC*Xmu) < 1e-12)

    # (h) g -> q qbar twice, two scalar vertices (quark and antiquark of equal helicity at each vertex):
    # gluon 1 -> quark (h, fraction x) + antiquark (-h); gluon 2 -> quark (h', fraction x) + antiquark (-h');
    # quark of 1 with antiquark of 2 needs h = -h', and so does the other pair
    for beta in (Fr(1, 10), Fr(1, 4), Fr(1, 2), Fr(2, 3), Fr(4, 5), Fr(19, 20)):
        x = (1 + beta)/2

        def Aq(l1, l2):
            return sum(n_qq(l1, a, x)*n_qq(l2, -a, x) for a in (1, -1))

        lhs = sum(Aq(l1, l2)**2 for l1 in (1, -1) for l2 in (1, -1))/4
        check('(h) scalars coupled to quarks: helicity average = (1 + beta^4)/4 = 2 T_qq^2 + 2 T_dqdq^2',
              lhs == (1 + beta**4)/4 == 2*T_qq(x)**2 + 2*T_dqdq(x)**2
              and (Aq(1, 1), Aq(1, -1)) == (2*x*(1 - x), x*x + (1 - x)**2))

    # (i) vector vertices: quark and antiquark of opposite helicity at each vertex, so h = h';
    # the two terms h = +1 and h = -1 are different final helicity states (++ and --) and add in the square
    for x in xs:
        tot = Fr(0)
        for l1 in (1, -1):
            for l2 in (1, -1):
                tot += sum((n_qq(l1, a, x)*n_qq(l2, a, x))**2 for a in (1, -1))
        check('(i) gg -> ZZ: helicity average = (x^2 + (1-x)^2)^2/2 = 2 T_qq^2',
              tot/4 == (x*x + (1 - x)**2)**2/2 == 2*T_qq(x)**2
              and (n_qq(1, 1, x)*n_qq(1, 1, x), n_qq(1, -1, x)*n_qq(1, -1, x), n_qq(1, 1, x)*n_qq(-1, 1, x))
              == (x*x, (1 - x)**2, x*(1 - x)))

    groups = {}
    for name, ok in out:
        g = groups.setdefault(name, [0, 0])
        g[0] += ok
        g[1] += 1
    for name, (p, t) in groups.items():
        print(f'{"PASS" if p == t else "FAIL"}  {p}/{t}  {name}')
    assert all(ok for _, ok in out)
    print('all identities hold (exact rational arithmetic; (g) in floating point)')


if __name__ == '__main__':
    main()
