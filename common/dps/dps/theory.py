"""The factorised DPS formula of Paper 2, process by process.

Expansion in a_s = alpha_s/(4 pi).  Units of the mass m of the produced particles, mu = m unless
stated.  For every process the one-loop DPS terms are

    M_1^(1) = C_1 G + ...,      M_2^(1) -> K,        C_1 = 2 pi N_c m^2 F(x) sigma,   K = pi N_c m^2 W(x),

with x = (1 + beta)/2, sigma = +1 and F = F_A for gg -> AA, sigma = -1 and F = F_S for gg -> HH
(eqs. 2.4-2.6), and the two-loop logarithms fixed by the renormalisation group, eq. (5.5):

    M_1^(2) in C_1 [ -(4/3) Dc G^3 - B G^2 ],      M_2^(2) in K [ -4 Dc G^2 + (4 Dc - 2 B) G ],

Dc = C_A and B = 2 beta_0 for gluon splittings.  For quark splittings (gg -> ZZ) Dc = 2 C_F - C_A,
B = 6 C_F, and the logarithm is G_q = G + d_q, d_q = -(C_A - C_F) L_x/Dc, L_x = ln[x(1 - x)].

The constant of the boundary term enters the next coefficients, eq. (7.1): in units of C_1 and K
and in minimal subtraction,

    c_11 = 2h + Dc L*^2 + B L* + r_0(nu = 0),
    c_20 = B - 2 Dc + 2h + Dc L*^2 + B L* + r_0(nu = 2),      L* = ln[mu^2/(Q*^2 e^{-i pi})],

with h the minimally subtracted one-loop amplitude of each 2 -> 1 subprocess over its tree.  The
remainders of the three processes are not minimally subtracted: Catani's operator for gg -> HH and
gg -> AA (minus [2 I^(1)]_{eps^0} times the one-loop amplitude), the q_T scheme at mu^2 = s for
gg -> ZZ (minus C_A zeta_2 times it).  `known_parts` returns everything but r_0 for each process,
colour factor by colour factor, so that r_0 = (fitted coefficient)/(one-loop term) - known.
Closed forms of r_0, eqs. (7.2) and (7.3), are in `r0_closed`.  For gluon splittings the helicity term
(2/3)(C_A - 2 T_F n_f)/F_nu enters with the sign rho_nu of the configuration it comes from relative to
the tree kernel: rho_0 = 1 (the same-helicity rational term, in the configuration whose weight enters
F_S and F_A alike), rho_2 = M_{--}/M_{++} of the hard vertex, +1 for the scalar and -1 for the
pseudoscalar (the helicity flip, which puts both daughters of one parent into the vertices with the
helicity opposite to the tree configuration).  gg -> AA on the forward line (October 2026): the nu = 0
component gives eq. (7.2) with F_A, the nu = 2 component the gg -> HH values minus
(4/3)(C_A - 2 T_F n_f)/W, i.e. the flip term with the opposite sign, at all five velocities.
"""
import mpmath as mp

TF = mp.mpf(1) / 2


def Z2():
    """zeta_2 at the current precision (never cached: the precision is set by the caller)."""
    return mp.zeta(2)

# Colour monomials of the two-loop remainders and the (N, n_f) weights they carry.
GLUON_COLOURS = ('N2', 'N_Nf', 'Nf_over_N', 'Nf2')


def F_A(x):
    """Helicity overlap of two g -> gg splittings for a CP-odd vertex, eq. (2.4)."""
    return 1 / (x * (1 - x)) - 2 * x * (1 - x)


def F_S(x):
    """The same for a CP-even vertex, eq. (2.5)."""
    return 1 / (x * (1 - x)) + 2 * x * (1 - x)


def W(x):
    """The |J_z| = 2 overlap, eq. (2.6)."""
    return (x ** 4 + (1 - x) ** 4) / (x * (1 - x))


def beta0(N, nf):
    return (11 * N - 2 * nf) / mp.mpf(3)


class Process:
    """What the DPS formula needs to know about one process."""

    def __init__(self, name, splitting, F0, sigma, hard, scheme_shift, label, rho2=None):
        self.name, self.splitting, self.F0, self.sigma = name, splitting, F0, sigma
        self.hard, self.scheme_shift, self.label = hard, scheme_shift, label
        self.rho2 = rho2            # M_{--}/M_{++} of the hard vertex: the sign of the helicity-flip term for nu = 2

    def rho(self, nu):
        """the sign of the helicity term of eq. (7.2) in component nu (gluon splittings)"""
        return 1 if nu == 0 else self.rho2

    # one-loop DPS terms per N_c (gluon splittings) ------------------------------------------
    def c1_over_N(self, x):
        """C_1/(N_c m^2)."""
        return 2 * mp.pi * self.sigma * self.F0(x)

    def K_over_N(self, x):
        return mp.pi * W(x)

    def F(self, nu, x):
        """F_nu of eq. (7.2): the x dependence of the two components of the tree kernel."""
        return self.F0(x) if nu == 0 else W(x)


def h_scalar(L):
    """Minimally subtracted one-loop gg -> H amplitude over its tree (effective vertex),
    L = ln(m^2/mu^2) - i pi."""
    return -L ** 2 + Z2()          # times C_A


def h_pseudoscalar(L):
    """The same for gg -> A through the operator G Gtilde, with the Wilson coefficient factored out.
    The renormalised one-loop form factor of that operator is
    C_A{-8/eps^2 + 22/(3 eps) + 4 + zeta_2} - 8 n_f T_F/(3 eps) at mu_R^2 = Q^2 (Ahmed, Gehrmann,
    Mathews, Rana, Ravindran, JHEP 11 (2015) 169, eq. (3.50), d = 4 + eps); the scalar one has
    zeta_2 in place of 4 + zeta_2.  Restoring the scale with (Q^2/mu^2)^{eps/2}, Q^2 = -m^2 - i0, and
    removing the poles gives -L^2 + 4 + zeta_2 in units of C_A: the pseudoscalar hard function
    differs from the scalar one by the constant 4 C_A.  Both components of gg -> AA test it: a wrong
    constant would shift the C_A part of r_0 by a number independent of x, and the data at five
    velocities leave no such shift (the nu = 2 component differs from gg -> HH by the x-dependent
    -(4/3)(C_A - 2 T_F n_f)/W alone, the sign change of the helicity flip)."""
    return -L ** 2 + 4 + Z2()          # times C_A


def h_quark(L):
    """Minimally subtracted one-loop q qbar -> Z amplitude over its tree, times C_F."""
    return -L ** 2 + 3 * L - 8 + Z2()


HH = Process('hh', 'gluon', F_S, -1, h_scalar, 'catani', 'gg -> HH', rho2=+1)
AA = Process('aa', 'gluon', F_A, +1, h_pseudoscalar, 'catani', 'gg -> AA', rho2=-1)
ZZ = Process('zz', 'quark', None, None, h_quark, 'qT', 'gg -> ZZ')
PROCESSES = {p.name: p for p in (HH, AA, ZZ)}


def predictions_gluon(x, N, nf, sigma, F0):
    """Two-loop leading logarithms for gluon splittings, in units of N_c m^2 (eq. 5.5):
    {'M1': {'G3': ..., 'G2': ...}, 'M2': {'G2': ..., 'G1': ...}}, each a polynomial in N, nf
    given as the coefficients of the colour monomials N^2 and N nf (nothing else appears)."""
    c1 = 2 * mp.pi * sigma * F0(x)
    K = mp.pi * W(x)
    # Dc = N, B = 2 beta0 = (22 N - 4 nf)/3; per N: Dc -> N^2 (N2), B -> (22/3) N^2 - (4/3) N nf
    return {'M1': {'G3': {'N2': -mp.mpf(4) / 3 * c1, 'N_Nf': 0},
                   'G2': {'N2': -mp.mpf(22) / 3 * c1, 'N_Nf': mp.mpf(4) / 3 * c1}},
            'M2': {'G2': {'N2': -4 * K, 'N_Nf': 0},
                   'G1': {'N2': (4 - mp.mpf(44) / 3) * K, 'N_Nf': mp.mpf(8) / 3 * K}}}


def known_parts_gluon(beta, process):
    """Everything but r_0 in c_11/C_1 and c_20/K for gluon splittings, eq. (7.1), at mu = m
    (Catani-subtracted remainders), as {'nu0': {'N': ..., 'nf': ...}, 'nu2': {...}}: the N part
    collects C_A and the N part of beta_0, the nf part the n_f part of beta_0.  `beta` is the
    velocity as a fraction string; s/m^2 = (1 + xk)^2/xk = 1/[x(1 - x)] and Q*^2 = m^4/s."""
    from .kinematics import beta_number
    b = beta_number(beta)
    v = (1 - b) / (1 + b)
    s = (1 + v) ** 2 / v                               # s/m^2
    Lstar = mp.log(s) + 1j * mp.pi                     # ln(mu^2/(Q*^2 e^{-i pi}))
    L = -1j * mp.pi                                    # ln(m^2/mu^2) - i pi
    Lsub = -mp.log(s) + 1j * mp.pi                     # ln(mu^2/s) + i pi

    def total(CA, b0):
        h = CA * process.hard(L)
        delta0 = -CA * Lsub ** 2 + CA * Z2() - b0 * Lsub   # [2 I^(1)]_{eps^0}
        return 2 * h + CA * Lstar ** 2 + 2 * b0 * Lstar - delta0

    N_part, nf_part = total(1, mp.mpf(11) / 3), total(0, -mp.mpf(2) / 3)
    nu2_N, nu2_nf = 2 * mp.mpf(11) / 3 - 2, 2 * (-mp.mpf(2) / 3)        # B - 2 Dc
    return {'nu0': {'N': N_part, 'nf': nf_part},
            'nu2': {'N': N_part + nu2_N, 'nf': nf_part + nu2_nf}, 'Lstar': Lstar}


def known_parts_quark(x, CF, CA):
    """Everything but r_0 for quark splittings in the variable G of the gluon scale, q_T scheme at
    mu^2 = s (r0/data/r0_from_zz.py): (nu0, nu2) for the colour factor (CF, CA) in {(1,0), (0,1)}."""
    Lx = mp.log(x * (1 - x))
    l = -2 * Lx + 1j * mp.pi
    Dc = 2 * CF - CA
    Dk = 2 * (CA - CF) * Lx
    B = 6 * CF
    L = Lx - 1j * mp.pi
    h = CF * h_quark(L)
    nu0 = 2 * h + Dc * l ** 2 - 2 * Dk * l + B * l - CA * Z2()
    return nu0, nu0 + (B - 2 * Dc - 2 * Dk)


def r0_closed(process, nu, x, N=1, nf=0):
    """Eqs. (7.2) and (7.3).  Gluon splittings: the parts proportional to C_A (N = 1, nf = 0) and to
    T_F n_f (N = 0, nf = 1) separately when called that way; quark splittings: (C_F part, C_A part)."""
    if process.splitting == 'gluon':
        CA = N
        helicity = process.rho(nu) * mp.mpf(2) / 3 * (CA - 2 * TF * nf) / process.F(nu, x)
        return CA * (mp.mpf(134) / 9 - 5 * Z2() - 4 * mp.log(x) * mp.log(1 - x)) - mp.mpf(40) / 9 * TF * nf + helicity
    lx, l1 = mp.log(x), mp.log(1 - x)
    CF_part = 12 - 6 * Z2() + 2 * (lx - l1) ** 2
    CA_part = 2 + Z2() - 2 * lx ** 2 - 2 * l1 ** 2
    return CF_part, CA_part
