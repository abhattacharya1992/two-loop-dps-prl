"""Blind NMHV check: same impact-parameter kernel as stage_d_blind.py, but with the "nu = 2" transform,
compared with the 3->3 NMHV ratio-function data rho^(L) (Dixon-Esterlis erratum, eqs. D.14-D.17).

(i)   kernel   a*(1+r0 a+r0b a^2+r0c a^3+r0d a^4) exp[-(Gamma/8) X^2 + (B2 a^2+B3 a^3+B4 a^4) X],
      Gamma = 4a - 4 z2 a^2 + 22 z4 a^3 + G4 a^4, expanded through a^5; the eight constants are taken
      from the MHV solve of stage_d_blind.py (that script is re-run by importing it, nothing typed by hand)
(ii)  X^n -> K * R_n(G),  R_n = n! [alpha^n] { exp(2(EulerGamma - G) alpha) Gamma(1+alpha)/Gamma(2-alpha) }
(iii) G = ln|delta|/2 = Ld/2;  result = K * sum_L a^L P_L(G);  coefficient of a^L Ld^k  <->  c_k of
      pi i (1+w)/(1+w*) c_k ln^k|delta| in rho^(L);  the L=1, k=0 coefficient (=1) defines K.
Conventions as in stage_d_blind.py: even zeta values as exact powers of pi, odd ones as symbols z3, z5, z7.
"""
import json, sys, os, time, functools
print = functools.partial(print, flush=True)
import sympy as sp
from sympy import Rational as R, pi, zeta, gamma, EulerGamma, exp, symbols, expand, factorial

t0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "de_data.json")   # path: the only edit after the blind run

# ---------- re-run the MHV analysis (kernel, nu = 0 transform, solve for the constants) ----------------
print("=== re-running stage_d_blind.py (MHV data -> constants) ===")
import stage_d_blind as mhv
print("=== stage_d_blind.py finished (%.1fs) ===\n" % (time.time() - t0))

a, X, G, al, Ld = mhv.a, mhv.X, mhv.G, mhv.al, mhv.Ld
z, even_to_pi, odd_to_sym, unk = mhv.z, mhv.even_to_pi, mhv.odd_to_sym, mhv.unk
r0, r0b, r0c, r0d, B2, B3, B4, G4 = unk
solved = dict(mhv.solved)
missing = [u for u in unk if u not in solved]
print("constants solved from the MHV data:", ", ".join(str(u) for u in solved), "| not fixed:", missing or "none")

def pi_to_zeta(e):
    """readable form: first all even zetas -> pi powers (the kernel's Gamma carries z2, z4 symbols while the
    solved constants carry pi powers), then pi powers -> zeta symbols"""
    e = expand(sp.sympify(e).subs(even_to_pi))
    try:
        return mhv.pi_to_zeta(e)
    except Exception:
        return e

# ---------- (a): build the same kernel here (copy of the stage_d_blind.py construction) ---------------
LMAX, NMAX = 5, 8                                  # a^1..a^5;  X^8 is the highest power at a^5
Gam = 4*a - 4*z[2]*a**2 + 22*z[4]*a**3 + G4*a**4
Bf = B2*a**2 + B3*a**3 + B4*a**4
pref = 1 + r0*a + r0b*a**2 + r0c*a**3 + r0d*a**4
E = -Gam/8*X**2 + Bf*X
def trunc(e, N):
    p = sp.Poly(expand(e), a)
    return sum(c*a**m[0] for m, c in p.terms() if m[0] <= N)
expE, Em = sp.Integer(1), sp.Integer(1)
for m in range(1, LMAX):
    Em = trunc(Em*E, LMAX - 1)
    expE += Em / factorial(m)
S = expand(trunc(a * pref * expE, LMAX))
assert expand(S - mhv.S) == 0, "kernel differs from the one built in stage_d_blind.py"
Spoly = sp.Poly(S, X)
assert max(m[0] for m, _ in Spoly.terms()) == NMAX
print("kernel built and identical to stage_d_blind.S; %d X-monomials (%.1fs)" % (len(Spoly.terms()), time.time() - t0))

# ---------- (b): nu = 2 transform,  R_n(G) from the generating function -------------------------------
gen = exp(2*(EulerGamma - G)*al) * gamma(1 + al) / gamma(2 - al)
ser = sp.series(gen, al, 0, NMAX + 1).removeO()
print("generating-function series done (%.1fs)" % (time.time() - t0))
ser = expand(sp.expand_func(ser).rewrite(zeta))
assert not ser.has(EulerGamma), "EulerGamma did not cancel"
ser = expand(ser.subs(odd_to_sym))
leftover = ser.atoms(sp.Function)
assert not leftover, f"unconverted functions in the series: {leftover}"
Rn = [expand(factorial(n) * ser.coeff(al, n)) for n in range(NMAX + 1)]
print("R_0 =", Rn[0], "| R_1 =", Rn[1], "| R_2 =", Rn[2], "| R_3 =", Rn[3])
assert Rn[0] == 1
assert expand(Rn[1] - (1 - 2*G)) == 0
assert expand(Rn[2] - (4*G**2 - 4*G + 2)) == 0
print("R_0, R_1, R_2 sanity checks passed.")
# independent cross-check (printed only): Gamma(1+x)/Gamma(2-x) = exp(-2 EulerGamma x - 2 sum_{odd k>=3} zeta(k) x^k/k)/(1-x)
cf = sum((-2*G*al - 2*sum(z[k]*al**k/k for k in (3, 5, 7)))**m / factorial(m) for m in range(NMAX + 1))
cf = expand(cf * sum(al**j for j in range(NMAX + 1)))
cf_ok = all(expand(factorial(n) * cf.coeff(al, n) - Rn[n]) == 0 for n in range(NMAX + 1))
print("closed-form cross-check of R_0..R_8 (exp(-2G x - 2 sum zeta_odd x^k/k)/(1-x)):", "agree" if cf_ok else "DISAGREE")

# ---------- (c): transform, substitute the constants, G = Ld/2, compare with the NMHV data ------------
T = sum(c * Rn[m[0]] for m, c in Spoly.terms())           # = sum_L a^L P_L(G), overall K stripped
T = expand(T.subs(solved))
Tp = sp.Poly(T, a, G)
P = {L: expand(sum(Tp.coeff_monomial(a**L * G**k) * G**k for k in range(0, 2*L - 1))) for L in range(1, LMAX + 1)}
for L in range(1, LMAX + 1):                               # P_L has degree <= 2L-2 in G
    assert expand(Tp.as_expr().coeff(a, L) - P[L]) == 0
print("transform applied, constants substituted (%.1fs)\n" % (time.time() - t0))

nmhv = json.load(open(DATA))["nmhv_three_to_three"]
loc = {str(s): s for s in z.values()}
def pred(L, k):
    return sp.Poly(expand(P[L].subs(G, Ld/2)), Ld).coeff_monomial(Ld**k)

n_cmp = n_pass = 0
norm_ok = None
for L in range(1, 5):
    coeffs = nmhv[str(L)]["coefficients"]
    for k in range(2*L - 2, -1, -1):
        p = pred(L, k)
        if str(k) not in coeffs:
            print(f"L={L} k={k}  no datum   pred = {pi_to_zeta(p)}")
            continue
        datum = sp.sympify(coeffs[str(k)], locals=loc)
        diff = expand((p - datum).subs(even_to_pi))
        free = [u for u in unk if diff.has(u)]
        ok = (diff == 0) and not free
        tag = "PASS" if ok else ("UNDETERMINED (%s)" % free if free else "FAIL")
        if L == 1 and k == 0:
            norm_ok = ok and (expand(p) == 1)
            print(f"L={L} k={k}  {tag}  pred = {pi_to_zeta(p)}  data = {datum}   [normalisation: defines K; pred == 1 exactly: {norm_ok}]")
            continue
        n_cmp += 1
        n_pass += ok
        if ok:
            print(f"L={L} k={k}  PASS  {pi_to_zeta(p)}")
        else:
            print(f"L={L} k={k}  {tag}  pred = {pi_to_zeta(p)}   data = {datum}   diff = {pi_to_zeta(diff)}")
print(f"\nSUMMARY: {n_cmp} coefficients compared (L=2..4, excluding the L=1 normalisation), {n_pass} PASS, {n_cmp - n_pass} FAIL")
print(f"L=1 normalisation coefficient == 1: {norm_ok}")

# ---------- (d): leading logarithms vs K exp(-2 a G^2) ------------------------------------------------
print("\nLeading-log check: coefficient of G^(2L-2) in P_L vs (-2)^(L-1)/(L-1)!  [from K exp(-2 a G^2)]")
ll_ok = True
for L in range(1, LMAX + 1):
    k = 2*L - 2
    got = sp.Poly(P[L], G).coeff_monomial(G**k)
    want = R(-2)**(L - 1) / factorial(L - 1)
    ok = expand(got - want) == 0
    ll_ok &= ok
    extra = ""
    if L <= 4:                                             # also against the data (in Ld = 2G)
        datum = sp.sympify(nmhv[str(L)]["coefficients"][str(k)], locals=loc)
        extra = f"   | Ld^{k}: pred {pred(L, k)}  data {datum}  {'ok' if expand(pred(L, k) - datum) == 0 else 'MISMATCH'}"
    print(f"  L={L} k={k}: pred {got}  expected {want}  {'PASS' if ok else 'FAIL'}{extra}")
print("LL tower (L=1..5) equals K exp(-2 a G^2):", ll_ok)
if not norm_ok or n_cmp!=15 or n_pass!=n_cmp or not ll_ok:
    raise ValueError('Independent NMHV coefficient comparison failed')
print("done (%.1fs)" % (time.time() - t0))
