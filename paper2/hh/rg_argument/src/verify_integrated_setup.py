#!/usr/bin/env python3
# Checks of the renormalisation-group argument of paper Sec. 5 (labels in brackets are
# those of the author's notes on it).
#
#   1. RG constraint on the ansatz A_DPS = H1 H2 K:  gamma_K = gamma_A - gamma_H1 - gamma_H2
#      [gammafinal], its integral from mu_b [integratedgamma, fullcandidate, RGendpoint],
#      the Delta_C, B table, and Q*^2 = m1^2 m2^2/s [G] for the gluon channel.
#   2. Bessel-Mellin moments F_m(alpha) [Fm, F0a-c, F2moments] and the resulting
#      A_0^(2), A_2^(2) polynomials in G [scalarRG, spinRG, adjointcandidate].
#   3. The boundary freedom R^(1) = K0 (r2 X^2 + r1 X) [boundary] and exactly how it
#      shifts those coefficients [boundarychanges]. This is what makes the argument
#      conditional: nothing here fixes r1, r2.
#   4. Kinematics: G = (1/2) ln(-q_T^2/Q*^2) = (1/2) ln(1 - yz) for the HH/AA variables
#      s = m^2 (1+rho)^2/rho, t = -m^2 y, u = -m^2 z; and the n_f form of the spin-2
#      single log, 4C_A - 4beta_0 = -(8/3)(4N_c - n_f).
#   5. The four opposite-beam daughter residues with one common +i delta
#      [n1triplet-Rb4, commonresidue], the longitudinal kernel J [Janswer], the
#      rapidity-region subtraction [cancel] and the radial double log [radial].
#
# Nothing here evaluates a two-loop integrand. Output: results/verification.json.
from __future__ import annotations
import json
import platform
from pathlib import Path
import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
checks: list[dict[str, object]] = []


def exact(name, actual, expected=0):
    residual = sp.simplify(sp.factor(sp.together(actual - expected)))
    checks.append({"name": name, "kind": "symbolic", "residual": str(residual)})
    if residual != 0:
        raise AssertionError(f"{name}: residual {residual}")
    print(f"ok  {name}")


def close(name, a, b, tol):
    err = abs(a - b) / max(mp.mpf(1), abs(b))
    checks.append({"name": name, "kind": "numerical", "error": mp.nstr(err, 6), "tol": mp.nstr(tol, 3)})
    if not err < tol:
        raise AssertionError(f"{name}: error {err}")
    print(f"ok  {name}  (err {mp.nstr(err, 3)})")


# ---------------------------------------------------------------- 1. RG constraint
CA, CF, nf, TF = sp.symbols("C_A C_F n_f T_F", positive=True)
beta0 = sp.Symbol("beta0")
Cpar, C1, C2, gpar, g1, g2 = sp.symbols("C_a C_R1 C_R2 gamma_a gamma_R1 gamma_R2")
DC = C1 + C2 - Cpar            # Delta_C   [defDC]
Dgam = gpar - g1 - g2          # D_gamma   [defDC]
Bcoef = beta0 + Dgam           # B         [fullcandidate]

# gamma_K / a_s = -4 Delta_C ln(Q*^2/mu^2) + 2 D_gamma.  With t = ln(mu/Q*) this is
# 8 Delta_C t + 2 D_gamma; integrate from mu_b (t = -L_b/2) to Q* (t = 0).
t, Lb = sp.symbols("t L_b")
exact("integrated anomalous dimension [integratedgamma]",
      sp.integrate(8 * DC * t + 2 * Dgam, (t, -Lb / 2, 0)), -DC * Lb**2 + Dgam * Lb)
# a_s(mu_b) = a_s(Q*)(1 + beta0 a_s L_b): the Born coupling adds beta0 L_b  [running]
exact("endpoint with Born running coupling [RGendpoint]",
      -DC * Lb**2 + Dgam * Lb + beta0 * Lb, -DC * Lb**2 + Bcoef * Lb)

# General-scale form. X = ln(Q*^2 b^2/b0^2) is mu-independent, Y = ln(mu^2 b^2/b0^2).
X, Y = sp.symbols("X Y")
bracket = DC * Y**2 - 2 * DC * X * Y + Bcoef * Y               # [fullcandidate]
# d/dln(mu): dY = 2, and the Born a_s(mu) contributes -2 beta0.
exact("d/dln(mu) of the candidate reproduces gamma_K [fullcandidate]",
      -2 * beta0 + 2 * sp.diff(bracket, Y), 4 * DC * (Y - X) + 2 * Dgam)
exact("candidate at mu^2 = Q*^2 e^{-i pi} [RGendpoint]", bracket.subs(Y, X), -DC * X**2 + Bcoef * X)
# Real hard scale mu^2 = Q*^2 instead: differs by a b-independent term only.
real_endpoint = -DC * Y**2 + (Bcoef + 2 * sp.I * sp.pi * DC) * Y
complex_endpoint = (-DC * X**2 + Bcoef * X).subs(X, Y - sp.I * sp.pi)
exact("timelike continuation differs by b-independent term",
      sp.diff(complex_endpoint - real_endpoint, Y))

# Representation table.  gamma_g^(0) = -beta0, gamma_q^(0) = -3 C_F  [gammaqg]
channels = {
    "gluon parent, two gluon fusions": ({Cpar: CA, C1: CA, C2: CA, gpar: -beta0, g1: -beta0, g2: -beta0}, CA, 2 * beta0),
    "gluon parent, two quark fusions": ({Cpar: CA, C1: CF, C2: CF, gpar: -beta0, g1: -3 * CF, g2: -3 * CF}, 2 * CF - CA, 6 * CF),
    "quark parent, quark and gluon fusions": ({Cpar: CF, C1: CF, C2: CA, gpar: -3 * CF, g1: -3 * CF, g2: -beta0}, CA, 2 * beta0),
}
table = {}
for name, (sub, dc_exp, b_exp) in channels.items():
    exact(f"table row Delta_C, B: {name}", (DC.subs(sub) - dc_exp)**2 + (Bcoef.subs(sub) - b_exp)**2)
    table[name] = {"Delta_C": str(dc_exp), "B": str(b_exp)}
gluon = channels["gluon parent, two gluon fusions"][0]
quark = channels["gluon parent, two quark fusions"][0]
# The ordered eikonal exchange carries C_a; it coincides with Delta_C only for gluons.
exact("C_a - Delta_C = 2(C_A - C_F) for quark fusions [colourdiscriminator]",
      (Cpar - DC).subs(quark), 2 * (CA - CF))
exact("C_a = Delta_C for gluon fusions", (Cpar - DC).subs(gluon))

# ------------------------------------------------- 2. Bessel-Mellin moments in G
# F_m(alpha) = (i^m/2) e^{(h+2gamma_E) alpha} Gamma(m/2+alpha)/Gamma(m/2+1-alpha),
# h = ln(Q*^2/q_T^2) = -2G  [Fm].  Series taken from the Gamma functions directly.
al, h, G = sp.symbols("alpha h G")
pref = sp.exp((h + 2 * sp.EulerGamma) * al)
F0 = sp.expand(sp.series(pref * sp.gamma(al) / (2 * sp.gamma(1 - al)), al, 0, 3).removeO().subs(h, -2 * G))
F2 = sp.expand(sp.series(-pref * sp.gamma(1 + al) / (2 * sp.gamma(2 - al)), al, 0, 3).removeO().subs(h, -2 * G))
# finite parts of F_0, F_0', F_0'' at alpha = 0 (the 1/(2 alpha) pole is the contact term)
F00, F01, F02 = F0.coeff(al, 0), F0.coeff(al, 1), 2 * F0.coeff(al, 2)
F20, F21, F22 = F2.coeff(al, 0), F2.coeff(al, 1), 2 * F2.coeff(al, 2)
exact("FP F_0(0) = -G [F0a]", F00, -G)
exact("FP F_0'(0) = G^2 [F0b]", F01, G**2)
exact("FP F_0''(0) = -4G^3/3 - 2 zeta_3/3 [F0c]", F02, -sp.Rational(4, 3) * G**3 - 2 * sp.zeta(3) / 3)
exact("F_2'(0)/F_2(0) = 1 - 2G [F2moments]", F21 / F20, 1 - 2 * G)
exact("F_2''(0)/F_2(0) = 4G^2 - 4G + 2 [F2moments]", F22 / F20, 4 * G**2 - 4 * G + 2)

# Insert -Delta_C X^2 + B X: X -> d/dalpha.  Scalar normalised as -C1 F_0, spin-two as K1 F_2/F_2(0).
Dsym, Bsym = sp.symbols("Delta_C B")
A0 = sp.expand(-(-Dsym * F02 + Bsym * F01))         # / C1
A2 = sp.expand((-Dsym * F22 + Bsym * F21) / F20)     # / K1
# A0 also carries the zeta_3 contact constant from F0c; only G^3, G^2 are compared.
exact("A_0^(2)/C_1 = -(4/3) Delta_C G^3 - B G^2 + O(G) [scalarRG]",
      A0 - A0.coeff(G, 0), -sp.Rational(4, 3) * Dsym * G**3 - Bsym * G**2)
exact("A_2^(2)/K_1 = -4 Delta_C G^2 + (4 Delta_C - 2B) G + const [spinRG]",
      sp.expand(A2 - A2.coeff(G, 0)) - (-4 * Dsym * G**2 + (4 * Dsym - 2 * Bsym) * G), 0)
adjoint = {Dsym: CA, Bsym: 2 * beta0}
exact("gluon channel: -4C_A/3, -2 beta0, -4C_A, 4C_A - 4beta0 [adjointcandidate]",
      (A0.coeff(G, 3).subs(adjoint) + sp.Rational(4, 3) * CA)**2
      + (A0.coeff(G, 2).subs(adjoint) + 2 * beta0)**2
      + (A2.coeff(G, 2).subs(adjoint) + 4 * CA)**2
      + (A2.coeff(G, 1).subs(adjoint) - (4 * CA - 4 * beta0))**2)

# ------------------------------------------------- 3. the unfixed boundary term
# R^(1) contains K0 (r2 X^2 + r1 X).  X has no mu dependence, so the RG equation above is
# unchanged; the Fourier moments are not.  [boundary, boundarychanges]
r1, r2 = sp.symbols("r_1 r_2")
dA0 = sp.expand(-(r2 * F02 + r1 * F01))
dA2 = sp.expand((r2 * F22 + r1 * F21) / F20)
exact("delta A_0^(2)/C_1 = (4/3) r2 G^3 - r1 G^2 [boundarychanges]",
      dA0 - dA0.coeff(G, 0), sp.Rational(4, 3) * r2 * G**3 - r1 * G**2)
exact("delta A_2^(2)/K_1 = 4 r2 G^2 + (-4 r2 - 2 r1) G [boundarychanges]",
      dA2 - dA2.coeff(G, 0), 4 * r2 * G**2 + (-4 * r2 - 2 * r1) * G)
tot0, tot2 = sp.expand(A0 + dA0), sp.expand(A2 + dA2)
exact("shifted coefficients -(4/3)(Delta_C - r2), -(B + r1), -4(Delta_C - r2)",
      (tot0.coeff(G, 3) + sp.Rational(4, 3) * (Dsym - r2))**2
      + (tot0.coeff(G, 2) + (Bsym + r1))**2
      + (tot2.coeff(G, 2) + 4 * (Dsym - r2))**2
      + (tot2.coeff(G, 1) - (4 * (Dsym - r2) - 2 * (Bsym + r1)))**2)

# ------------------------------------------------- 4. kinematics and n_f form
# 2 -> 2 with massless beams and equal masses: tu - m^4 = s p_T^2 exactly.
s, pT, m = sp.symbols("s p_T m", positive=True)
E = sp.sqrt(s) / 2
pz = sp.sqrt(E**2 - pT**2 - m**2)
tt = m**2 - 2 * E * (E - pz)      # (p_A - q_1)^2, p_A = (E,0,0,E)
uu = m**2 - 2 * E * (E + pz)      # (p_B - q_1)^2
exact("tu - m^4 = s p_T^2", sp.expand(tt * uu) - m**4, s * pT**2)
# Amplitude variables y = -t/m^2, z = -u/m^2 (HH/AA files), Q*^2 = m^4/s for equal masses [G].
yv, zv = -tt / m**2, -uu / m**2
Qstar2 = m**4 / s
exact("-q_T^2/Q*^2 = 1 - yz, i.e. G = (1/2) ln(-q_T^2/Q*^2) = (1/2) ln(1 - yz)",
      sp.expand(-pT**2 / Qstar2), sp.expand(1 - yv * zv))
# Forward parametrisation used for the HH fits: s = m^2 (1+rho)^2/rho, y = rho + d, z = 1/rho - d.
rho, dlt, beta = sp.symbols("rho delta beta", positive=True)
s_rho = m**2 * (1 + rho)**2 / rho
exact("y + z = s/m^2 - 2 for y = rho + delta, z = 1/rho - delta",
      (rho + dlt) + (1 / rho - dlt), s_rho / m**2 - 2)
exact("1 - yz = -delta (1/rho - rho - delta)", 1 - (rho + dlt) * (1 / rho - dlt), -dlt * (1 / rho - rho - dlt))
exact("Q*^2 = m^2 (1 - beta^2)/4 for rho = (1 - beta)/(1 + beta)",
      (m**4 / s_rho).subs(rho, (1 - beta) / (1 + beta)), m**2 * (1 - beta**2) / 4)
# Same branch numerically: ln of a negative real takes +i pi in both forms.
with mp.workdps(30):
    rho_n, d_n, m_n = mp.mpf(1) / 3, mp.mpf("1e-3"), mp.mpf(1)
    s_n = m_n**2 * (1 + rho_n)**2 / rho_n
    yz_n = (rho_n + d_n) * (1 / rho_n - d_n)
    qT2_n = m_n**4 * (yz_n - 1) / s_n
    close("G from (1/2) log(1 - yz) equals (1/2) log(q_T^2/Q*^2) + i pi/2",
          mp.log(mp.mpc(1 - yz_n, 0)) / 2, mp.log(qT2_n / (m_n**4 / s_n)) / 2 + 1j * mp.pi / 2, mp.mpf("1e-25"))

Nc = sp.Symbol("N_c", positive=True)
beta0_qcd = (11 * CA - 4 * TF * nf) / 3
exact("4C_A - 4beta0 = -(8/3)(4N_c - n_f) for T_F = 1/2",
      (4 * CA - 4 * beta0_qcd).subs({CA: Nc, TF: sp.Rational(1, 2)}), -sp.Rational(8, 3) * (4 * Nc - nf))

# ------------------------------------------------- 5. causal residues and J
z, w = sp.symbols("z w", positive=True)
P, PB, A, Ap, Bt, Bp, d = sp.symbols("P PB A Ap B Bp delta", positive=True)
tau, lp, lm = sp.symbols("tau ell_plus ell_minus", real=True)
I = sp.I
Omega = A / ((1 - z) * P) + Ap / (z * P)          # [omegaxi]
Xi = Bt / ((1 - w) * PB) + Bp / (w * PB)
dn = d / (z * (1 - z) * P)                         # [shifteddelta]
db = d / (w * (1 - w) * PB)

# Every scalar denominator carries the same +i delta; the daughter eikonals get no
# independent i0 sign.  Close on the spectator pole in each case.
n1 = [z * P * tau - A + I * d, -(1 - z) * P * tau - A + I * d, z * P * (tau - lm) - Ap + I * d]
p1 = (-A + I * d) / ((1 - z) * P)                                       # upper half plane
R_n1 = I / (-(1 - z) * P * n1[0].subs(tau, p1) * n1[2].subs(tau, p1))
E_n1 = -I / (z * P**2 * (A - I * d) * (lm + Omega - I * dn))
n2 = [z * P * tau - Ap + I * d, -(1 - z) * P * tau - Ap + I * d, -(1 - z) * P * (tau + lm) - A + I * d]
p2 = (Ap - I * d) / (z * P)                                             # lower, clockwise
R_n2 = -I / (z * P * n2[1].subs(tau, p2) * n2[2].subs(tau, p2))
E_n2 = -I / ((1 - z) * P**2 * (Ap - I * d) * (lm + Omega - I * dn))
b3 = [-w * PB * tau - Bt + I * d, (1 - w) * PB * tau - Bt + I * d, -w * PB * (tau - lp) - Bp + I * d]
p3 = (Bt - I * d) / ((1 - w) * PB)                                      # lower
R_b3 = -I / ((1 - w) * PB * b3[0].subs(tau, p3) * b3[2].subs(tau, p3))
E_b3 = I / (w * PB**2 * (Bt - I * d) * (lp - Xi + I * db))
b4 = [-w * PB * tau - Bp + I * d, (1 - w) * PB * tau - Bp + I * d, (1 - w) * PB * (tau + lp) - Bt + I * d]
p4 = (-Bp + I * d) / (w * PB)                                           # upper
R_b4 = I / (-w * PB * b4[1].subs(tau, p4) * b4[2].subs(tau, p4))
E_b4 = I / ((1 - w) * PB**2 * (Bp - I * d) * (lp - Xi + I * db))
for name, got, want in [("R_n1 [Rn1]", R_n1, E_n1), ("R_n2 [Rn2]", R_n2, E_n2),
                        ("R_b3 [Rb3]", R_b3, E_b3), ("R_b4 [Rb4]", R_b4, E_b4)]:
    exact(name, got, want)

# Eikonal numerators 4 p_i.p_j = 2 z_i w_j s and the two 1/2 light-cone Jacobians
# remove all fraction prefactors: the four products share one causal denominator. [commonresidue]
common = 1 / ((lm + Omega - I * dn) * (lp - Xi + I * db))
mismatch = 0
for zi, rn, aa in [(z, E_n1, A), (1 - z, E_n2, Ap)]:
    for wj, rb, bb in [(w, E_b3, Bt), (1 - w, E_b4, Bp)]:
        mismatch += sp.simplify(sp.Rational(1, 4) * 2 * zi * wj * P * PB * rn * rb
                                - common / (2 * P * PB * (aa - I * d) * (bb - I * d)))
exact("all four pairs reduce to the common residue [commonresidue]", mismatch)
# Diagonal limit k' -> k: Omega Xi = A B / Q_kin^2 with Q_kin^2 = s z(1-z)w(1-w) = m1^2 m2^2/s. [diagscale]
exact("Omega Xi -> A B / Q_kin^2 in the diagonal limit [diagscale]",
      (Omega * Xi).subs({Ap: A, Bp: Bt}), A * Bt / (P * PB * z * (1 - z) * w * (1 - w)))

# Rapidity-region subtraction for the reduced integral [regionintegrals, cancel].
eta, nu, Om, Xs, T = sp.symbols("eta nu Omega Xi T", positive=True)
I_low = sp.pi * (nu / Xs)**eta * (sp.cot(sp.pi * eta) - sp.I)
I_high = -sp.pi * (nu / (T / Om))**eta / sp.sin(sp.pi * eta)
exact("I_low + I_high = ln(T/(Omega Xi)) - i pi + O(eta) [cancel]",
      sp.simplify(sp.series(I_low + I_high, eta, 0, 1).removeO()), sp.log(T / (Om * Xs)) - sp.I * sp.pi)

# ---- numerics
mp.mp.dps = 60
# Real-axis integration of each finite-delta triplet, independent of the residue algebra.
pars = {z: sp.Rational(3, 10), w: sp.Rational(3, 5), P: 5, PB: 7, A: sp.Rational(4, 5), Ap: sp.Rational(11, 10),
        Bt: sp.Rational(9, 10), Bp: sp.Rational(7, 10), d: sp.Rational(1, 5), lm: sp.Rational(1, 4), lp: sp.Rational(17, 100)}
for name, trip, want in [("n1", n1, E_n1), ("n2", n2, E_n2), ("b3", b3, E_b3), ("b4", b4, E_b4)]:
    f = sp.lambdify(tau, (1 / sp.prod(trip)).subs(pars), "mpmath")
    target = sp.lambdify((), want.subs(pars), "mpmath")()
    got = mp.quad(f, [-mp.inf, -2, -1, mp.mpf("-.3"), 0, mp.mpf(".3"), 1, 2, mp.inf]) / (2 * mp.pi)
    close(f"real-axis integral of triplet {name} at finite delta", got, target, mp.mpf("1e-40"))


def J_closed(om, xi, tau_):
    # J = -1/(2(T+Omega Xi)) - i ln(T/(Omega Xi)) / (2 pi (T+Omega Xi))   [Janswer]
    den = tau_ + om * xi
    return -1 / (2 * den) - 1j * mp.log(tau_ / (om * xi)) / (2 * mp.pi * den)


def J_pv(om, xi, tau_):
    # J = -(i/2pi) int_0^inf dx / ((x - Xi + i0)(T + Omega x))   [Jaftercontour], PV plus -i pi delta
    gap = xi * mp.mpf("1e-18")
    f = lambda x: 1 / ((x - xi) * (tau_ + om * x))
    pv = mp.quad(f, [0, xi - gap]) + mp.quad(f, [xi + gap, 2 * xi, 10 * xi, mp.inf])
    return -1j * (pv - 1j * mp.pi / (tau_ + om * xi)) / (2 * mp.pi)


for k, vals in enumerate([("0.2", "0.7", "0.6"), ("1.3", "0.4", "2.1"), ("0.03", "0.06", "0.9")], 1):
    om, xi, tau_ = map(mp.mpf, vals)
    close(f"J by principal-value quadrature, set {k} [Janswer]", J_pv(om, xi, tau_), J_closed(om, xi, tau_), mp.mpf("1e-16"))
    close(f"J with Omega, Xi interchanged, set {k}", J_pv(xi, om, tau_), J_closed(om, xi, tau_), mp.mpf("1e-16"))

# Radial integral int_0^R ln x/(1+x) dx = ln R ln(1+R) + Li2(-R) -> (1/2) ln^2 R - pi^2/6.  [radial]
for R in [mp.mpf("0.3"), mp.mpf("2"), mp.mpf("1000")]:
    quadv = mp.quad(lambda x: mp.log(x) / (1 + x), [0, min(1, R), R])
    close(f"radial log integral, R={R} [radial]", quadv, mp.log(R) * mp.log1p(R) + mp.polylog(2, -R), mp.mpf("1e-45"))

# Bessel-Mellin transform int_0^inf r^{2a-1} J_m(r) dr by quadrature, alpha = 0.2.  [Besselmellin]
with mp.workdps(30):
    av = mp.mpf("0.2")
    for mm in (0, 2):
        f = lambda r: r**(2 * av - 1) * mp.besselj(mm, r)
        val = mp.quad(f, [0, 1]) + mp.quadosc(f, [1, mp.inf], zeros=lambda n: mp.besseljzero(mm, n))
        ref = 2**(2 * av - 1) * mp.gamma(mp.mpf(mm) / 2 + av) / mp.gamma(mp.mpf(mm) / 2 + 1 - av)
        close(f"Bessel-Mellin integral m={mm} [Besselmellin]", val, ref, mp.mpf("1e-12"))

# ---- summary
summary = {
    "python": platform.python_version(), "sympy": sp.__version__, "mpmath": mp.__version__,
    "representation_table": table,
    "A0_over_C1": str(A0), "A2_over_K1": str(A2),
    "gluon_channel": {"A0_over_C1": str(sp.expand(A0.subs(adjoint))), "A2_over_K1": str(sp.expand(A2.subs(adjoint)))},
    "boundary_shift": {"dA0_over_C1": str(dA0), "dA2_over_K1": str(dA2),
                       "note": "r1, r2 are not determined by the RG argument"},
    "kinematics": "-q_T^2/Q*^2 = 1 - yz exactly; 4C_A - 4beta0 = -(8/3)(4N_c - n_f)",
    "checks": checks,
}
(OUT / "verification.json").write_text(json.dumps(summary, indent=2) + "\n")
print()
print("gluon channel: Delta_C = C_A, B = 2 beta0")
print("A0/C1 =", summary["gluon_channel"]["A0_over_C1"])
print("A2/K1 =", summary["gluon_channel"]["A2_over_K1"])
print("boundary shift dA0/C1 =", dA0, ";  dA2/K1 =", dA2)
