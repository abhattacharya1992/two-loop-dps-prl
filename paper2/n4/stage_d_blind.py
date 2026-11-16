"""Blind re-implementation: DPS-singularity structure vs Dixon-Esterlis 3->3 data.

(i)  impact-parameter kernel  a*(1+r0 a+r0b a^2+r0c a^3+r0d a^4) exp[-(Gamma/8) X^2 + (B2 a^2+B3 a^3+B4 a^4) X]
(ii) X^n -> -C1 * F_n(G),  F_n = n! [alpha^n] { Phi_0(alpha) - 1/(2 alpha) }
(iii) G = ln|delta|/2,  C1 = 2 pi i
"""
import json, sys, time, functools
print = functools.partial(print, flush=True)
import sympy as sp
from sympy import Rational as R, pi, zeta, gamma, EulerGamma, exp, symbols, expand, factorial, I

t0 = time.time()
import os
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "de_data.json")   # path: the only edit after the blind run
d = json.load(open(DATA))["three_to_three"]

a, X, G, al, Ld, C1 = symbols("a X G alpha Ld C1")          # Ld = ln|delta|
r0, r0b, r0c, r0d, B2, B3, B4, G4 = unk = symbols("r0 r0b r0c r0d B2 B3 B4 G4")
z = {k: symbols(f"z{k}") for k in range(2, 10)}               # symbols used in the data strings
even_to_pi = {z[2]: pi**2/6, z[4]: pi**4/90, z[6]: pi**6/945, z[8]: pi**8/9450}
odd_to_sym = {zeta(3): z[3], zeta(5): z[5], zeta(7): z[7], zeta(9): z[9]}

def pi_to_zeta(e):
    """rewrite pi^(2m) -> zeta(2m) symbols (weight-2m pure even zeta)"""
    e = expand(e)
    e = e.replace(lambda t: t.is_Pow and t.base == pi, lambda t: {2: 6*z[2], 4: 90*z[4], 6: 945*z[6], 8: 9450*z[8]}[int(t.exp)])
    return sp.factor(e) if e.is_polynomial() else e

# ---------- (ii): F_n from the Euler-gamma ratio -------------------------------------------
LMAX, NMAX = 5, 8                                  # X^8 is the highest power at a^5
Phi0 = R(1, 2) * exp(2*(EulerGamma - G)*al) * gamma(al) / gamma(1 - al)
rat = sp.series(gamma(al) / gamma(1 - al), al, 0, NMAX + 1).removeO()
print("gamma ratio series done (%.1fs)" % (time.time() - t0))
ex = sum((2*(EulerGamma - G)*al)**m / factorial(m) for m in range(NMAX + 2))
ser = expand(R(1, 2) * expand(rat) * ex)
ser = sum(ser.coeff(al, n) * al**n for n in range(-1, NMAX + 1))
ser = expand(sp.expand_func(ser).rewrite(zeta))
assert ser.coeff(al, -1) == R(1, 2), "pole is not 1/(2 alpha)"
reg = expand(ser - 1/(2*al))
assert not reg.has(EulerGamma)
F = [expand(factorial(n) * reg.coeff(al, n)).subs(odd_to_sym) for n in range(NMAX + 1)]
print("F0,F1,F2 =", F[0], "|", F[1], "|", F[2])
assert F[0] == -G and F[1] == G**2 and expand(F[2] - (-R(4, 3)*G**3 - R(2, 3)*z[3])) == 0
print("F_0..F_2 verified.  (%.1fs)" % (time.time() - t0))

# ---------- (i): expand the impact-parameter kernel through a^5 ------------------------------
Gam = 4*a - 4*z[2]*a**2 + 22*z[4]*a**3 + G4*a**4
Bf = B2*a**2 + B3*a**3 + B4*a**4
pref = 1 + r0*a + r0b*a**2 + r0c*a**3 + r0d*a**4
E = -Gam/8*X**2 + Bf*X                       # exponent, O(a)
def trunc(e, N):                              # keep a^0..a^N
    p = sp.Poly(expand(e), a)
    return sum(c*a**m[0] for m, c in p.terms() if m[0] <= N)
expE, Em = sp.Integer(1), sp.Integer(1)
for m in range(1, LMAX):                      # exp(E) through a^(LMAX-1)
    Em = trunc(Em*E, LMAX - 1)
    expE += Em / factorial(m)
print("exp expanded (%.1fs)" % (time.time() - t0))
S = trunc(a * pref * expE, LMAX)              # a^1 .. a^5
S = expand(S)

# ---------- apply (ii) and (iii) ---------------------------------------------------------------
Spoly = sp.Poly(S, X)
Emom = sum(c * (-C1) * F[m[0]] for m, c in Spoly.terms())
Emom = expand(Emom.subs(G, Ld/2))             # C1 kept as overall factor; coefficient of C1*Ld^k == coefficient of 2 pi i ln^k
Emom = expand(Emom / C1)
assert not Emom.has(C1)
Ep = sp.Poly(Emom, a, Ld)
def pred(L, k):
    return Ep.coeff_monomial(a**L * Ld**k)
print("S has %d X-monomials" % len(Spoly.terms()))
print("kernel expanded and transformed.  (%.1fs)" % (time.time() - t0))

# ---------- compare ------------------------------------------------------------------------------
solved = {}
results = []
for L in range(1, LMAX + 1):
    coeffs = d[str(L)]["coefficients"]
    for k in sorted((int(x) for x in coeffs), reverse=True):
        val = coeffs[str(k)]
        if val == "unknown":
            results.append((L, k, "no datum")); continue
        datum = sp.sympify(val, locals={str(s): s for s in z.values()})
        p = pred(L, k).subs(solved)
        diff = expand((p - datum).subs(even_to_pi))
        free = [u for u in unk if diff.has(u)]
        if len(free) == 0:
            ok = diff == 0
            results.append((L, k, "PASS" if ok else "FAIL"))
            if not ok:
                print(f"  FAIL L={L} k={k}: pred={pi_to_zeta(p.subs(even_to_pi))}  data={datum}")
        elif len(free) == 1:
            sol = sp.solve(diff, free[0])
            assert len(sol) == 1
            solved[free[0]] = expand(sol[0])
            results.append((L, k, f"fixes {free[0]} = {pi_to_zeta(sol[0])}"))
        else:
            results.append((L, k, f"UNDETERMINED: {free}"))
for r in results:
    print("L=%d k=%d  %s" % r)

print("\nSolved constants (pi powers -> zeta):")
for u in unk:
    print(f"  {u} = {pi_to_zeta(solved[u]) if u in solved else 'not fixed'}")
if G4 in solved:
    target = (-R(219, 2)*z[6] - 4*z[3]**2).subs(even_to_pi)
    print("G4 == -(219/2) z6 - 4 z3^2 ?", expand(solved[G4] - target) == 0)

# ---------- leading-log tower -------------------------------------------------------------------
t = symbols("t")
LLint = -G * sp.integrate(sp.series(exp(-2*a*G**2*t**2), a, 0, LMAX).removeO(), (t, 0, 1))
LLint = expand(LLint.subs(G, Ld/2))            # = -G int_0^1 dt exp(-2 a G^2 t^2), G = ln|delta|/2
ok_data = ok_pred = True
for L in range(1, LMAX + 1):
    k = 2*L - 1
    tower = -LLint.coeff(a, L - 1).coeff(Ld, k)   # E_LL = -C1 * a * [ -G int ... ]  => c_{2L-1} = - [a^(L-1) Ld^k] of the integral
    datum = sp.sympify(d[str(L)]["coefficients"][str(k)])
    ok_data &= (tower == datum)
    ok_pred &= (expand(pred(L, k).subs(solved)) == tower)
    print(f"  LL L={L} k={k}: integral {tower}  data {datum}  pred {pred(L, k)}")
print("LL tower: data match", ok_data, "| prediction match", ok_pred)
if (any(status.startswith(('FAIL','UNDETERMINED')) for _,_,status in results) or
        not ok_data or not ok_pred or G4 not in solved or
        expand(solved[G4]-(-R(219,2)*z[6]-4*z[3]**2).subs(even_to_pi))!=0):
    raise ValueError('Independent MHV coefficient comparison failed')
print("done (%.1fs)" % (time.time() - t0))
