#!/usr/bin/env python3
"""The b-space structure of Paper 2 (Sec. 5) against the self-crossing limit of the planar N = 4
six-gluon amplitudes of Dixon and Esterlis (DE), arXiv:1602.02107 (MHV, through five loops) and
the erratum JHEP 08 (2016) 131 (NMHV, through four loops).  Paper 2, Sec. 5.1, Tables 5 and 6.

Structure tested (Paper 2, eqs. (5.3)-(5.6), conformal limit, no running coupling):
    K~(b) = a K0(b) [1 + r0 a + r0b a^2 + r0c a^3 + r0d a^4] exp[ -(Gamma(a)/8) X^2 + B(a) X ],
    X = ln(Q_*^2 e^{-i pi} b^2/b_0^2),  Gamma(a) = gamma_K = 4a - 4 z2 a^2 + 22 z4 a^3 + G4 a^4,
    B(a) = B2 a^2 + B3 a^3 + B4 a^4 (no one-loop non-cusp term in N = 4; B1 = beta_0 + D_gamma = 0),
(r0, r0b, ... are kappa_1, kappa_2, ... of the paper) with the transforms of Paper 2, eq. (5.5):
    nu = 0 (MHV, |Delta J_z| = 0):  X^n -> -C_1 F_n(G),  F_n = n! [alpha^n] { Phi_0(alpha) - 1/(2 alpha) },
        Phi_0(alpha) = (1/(2 alpha)) exp(-2 G alpha - sum_{k odd >= 3} (2 zeta_k/k) alpha^k),
        so that the one-loop DPS term is C_1 G;
    nu = 2 (NMHV, |Delta J_z| = 2):  X^n -> K R_n(G),  R_n = n! [alpha^n] { (1/(1 - alpha)) exp(same) },
        so that the one-loop term is K.
Dictionary (Secs. 2.3 and 3 of DE): delta = p_T^2/Q_*^2 with Q_*^2 = s_36 x y (1-x)(1-y) (DE eq. 2.22
= Paper 2's m_1^2 m_2^2/s); for 3 -> 3 kinematics G = ln|delta|/2 and C_1 = 2 pi i; for 2 -> 4,
G = ln(delta)/2 + i pi/2 and C_1 = -2 pi i.  The one-loop terms define C_1 and K (normalisation);
the unknown constants are solved from the lowest loop order at which each appears in the MHV
data; everything else, MHV and NMHV (constants k = 0 included), is a prediction.  Also computed:
what the three-loop ln^3 tests, the exact relation between the two components, the five-loop
prediction for the |Delta J_z| = 2 component, and the 2 -> 4 two-loop coefficients.

    python3 stage_d.py [de_data.json] [results_stage_d.json]       (about 10 s)
"""
import json, pathlib, sys
import sympy as sp

HERE = pathlib.Path(__file__).resolve().parent
DATA = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else HERE/'de_data.json'
OUT = pathlib.Path(sys.argv[2]) if len(sys.argv) > 2 else HERE/'results_stage_d.json'

a, X, G, al, lnd, v = sp.symbols('a X G alpha lnd v')
# odd zeta values are independent symbols; even ones are the exact powers of pi, so that products
# such as zeta_2 zeta_4 = 7 zeta_6/4 are reduced; for printing, powers of pi are turned back into
# the single even zeta value of the same weight (to_zeta)
z = {k: (sp.zeta(k) if k % 2 == 0 else sp.Symbol(f'zeta_{k}')) for k in range(2, 10)}
PIPOW = {2*n: sp.pi**(2*n)/sp.zeta(2*n) for n in range(1, 5)}      # 6, 90, 945, 9450

def to_zeta(expr):
    """pi^(2n) -> (pi^(2n)/zeta_(2n)) zeta_(2n) with zeta_(2n) a symbol, for display"""
    expr = sp.expand(expr)
    out = 0
    for term in sp.Add.make_args(expr):
        c, rest = term.as_independent(sp.pi)
        deg = sp.degree(rest, sp.pi) if rest.has(sp.pi) else 0
        if deg == 0:
            out += term
        else:
            assert deg % 2 == 0 and rest == sp.pi**deg, term
            out += c*PIPOW[deg]*sp.Symbol(f'zeta_{deg}')
    return out
I, pi = sp.I, sp.pi
LMAX = 5                                     # loops compared (DE give E^(1) ... E^(5))
NA = LMAX - 1                                # order in a relative to the tree-level kernel

# ---------------------------------------------------------------- the b-space structure
r0, r0b, r0c, r0d, B2, B3, B4, G4 = sp.symbols('r0 r0b r0c r0d B2 B3 B4 G4')
Gamma = 4*a - 4*z[2]*a**2 + 22*z[4]*a**3 + G4*a**4      # G4: the four-loop cusp, solved from the data
# B_n absorbs the boundary terms r_1^(n) X of the same order, which the data cannot separate from the non-cusp anomalous dimension
expo = -(Gamma/8)*X**2 + (B2*a**2 + B3*a**3 + B4*a**4)*X
const = 1 + r0*a + r0b*a**2 + r0c*a**3 + r0d*a**4
bspace = sp.series(sp.exp(expo)*const, a, 0, NA + 1).removeO()
bspace = sp.expand(bspace)

# ---------------------------------------------------------------- the transform X^n -> -C_1 F_n(G)
NMAX = 2*NA
gen = sp.exp(-2*G*al - sum(2*z[k]/k*al**k for k in (3, 5, 7)))/(2*al) - 1/(2*al)
gen = sp.series(gen, al, 0, NMAX + 1).removeO()
F = [sp.expand(sp.factorial(n)*gen.coeff(al, n)) for n in range(NMAX + 1)]
assert F[0] == -G and sp.expand(F[1] - G**2) == 0 and sp.expand(F[2] + sp.Rational(4, 3)*G**3 + sp.Rational(2, 3)*z[3]) == 0

def transform(poly_in_X):
    """polynomial in X (coefficients may contain a) -> amplitude in units of -C_1"""
    p = sp.Poly(sp.expand(poly_in_X), X)
    return sp.expand(sum(c*F[m[0]] for m, c in zip(p.monoms(), p.coeffs())))

amp = transform(bspace)                      # sum_L a^{L-1} A^(L) in units of -C_1
A = {L: sp.expand(amp.coeff(a, L - 1)) for L in range(1, LMAX + 1)}   # A^(L)/(-C_1)
assert A[1] == -G                            # one-loop DPS term: -C_1 (-G) = C_1 G

# ---------------------------------------------------------------- DE, 3 -> 3: E^(L) = 2 pi i [ sum_k c_k ln^k|delta| + ... ]
data = json.loads(DATA.read_text())
zz = {f'z{k}': z[k] for k in range(2, 10)}
pred33 = {L: sp.expand(sp.expand(A[L]*(-2*pi*I)/(2*pi*I)).subs(G, lnd/2)) for L in A}   # C_1 = 2 pi i, G = ln|delta|/2
raw = {}                                     # the coefficients as printed, for the table
unknowns = [r0, r0b, r0c, r0d, B2, B3, B4, G4]
eqs, solved, report = [], {}, []
for L in range(1, LMAX + 1):
    de = data['three_to_three'][str(L)]['coefficients']
    for k in range(2*L - 1, 0, -1):
        p = pred33[L].coeff(lnd, k)
        p = p.subs(solved)
        if str(k) not in de:
            continue
        if de[str(k)] == 'unknown':
            report.append((L, k, 'not transcribed', str(sp.simplify(p))))
            continue
        target = sp.sympify(de[str(k)], locals=zz)
        raw[(L, k)] = de[str(k)]
        diff = sp.expand(p - target)
        free = [u for u in unknowns if diff.has(u)]
        if (L, k) == (1, 1):                     # the one-loop term defines C_1 and G: not a prediction
            assert diff == 0
            report.append((L, k, 'normalisation', str(target)))
        elif diff == 0:
            report.append((L, k, 'PASS', str(target)))
        elif free:
            sol = sp.solve(diff, free[0])
            assert len(sol) == 1, (L, k, free, sol)
            solved[free[0]] = sp.expand(sol[0])
            report.append((L, k, f'fixes {free[0]} = {to_zeta(solved[free[0]])}', str(target)))
        else:
            report.append((L, k, f'FAIL: predicted {sp.simplify(p)}, DE {target}', str(target)))

print('3 -> 3 kinematics, coefficients of 2 pi i ln^k|delta| in E^(L):')
for L, k, status, target in report:
    print(f'  L={L} k={k}: {status}')
print('constants solved from the data:', {str(k): str(to_zeta(val)) for k, val in solved.items()})
bes = -sp.Rational(219, 2)*z[6] - 4*z[3]**2
g4_ok = G4 in solved and sp.expand(solved[G4] - bes) == 0
if G4 in solved:
    print('four-loop cusp from the data:', to_zeta(solved[G4]), '; Beisert-Eden-Staudacher: -(219/2) z6 - 4 z3^2 ->', g4_ok)

# ---------------------------------------------------------------- the |Delta J_z| = 2 component: NMHV (erratum, App. D)
# transform nu = 2: X^n -> K R_n(G), R_n = n-th derivative at 0 of Phi_2(alpha)/Phi_2(0) =
# e^{2(gamma_E - G) alpha} Gamma(1+alpha)/Gamma(2-alpha); the one-loop term is K, which is
# pi i (1+w)/(1+w*) in rho of the erratum; the same constants as for nu = 0.  With
# ln Gamma(1+alpha) = -gamma_E alpha + sum_{k>=2} (-1)^k zeta_k alpha^k/k the ratio is
# (1/(1-alpha)) exp(-2 G alpha - sum_{k odd >= 3} (2 zeta_k/k) alpha^k): the same exponential
# as for nu = 0, with 1/(1-alpha) in place of 1/(2 alpha); the closed form is checked against
# the Gamma functions through alpha^4
Rgen = sp.series(sp.exp(-2*G*al - sum(2*z[k]/k*al**k for k in (3, 5, 7)))/(1 - al), al, 0, NMAX + 1).removeO()
Rn = [sp.expand(sp.factorial(n)*Rgen.coeff(al, n)) for n in range(NMAX + 1)]
Rgam = sp.series(sp.exp(2*(sp.EulerGamma - G)*al)*sp.gamma(1+al)/sp.gamma(2-al), al, 0, 5).removeO()
Rgam = Rgam.subs({sp.zeta(3): z[3], sp.zeta(5): z[5]})
assert all(sp.expand(sp.factorial(n)*Rgam.coeff(al, n) - Rn[n]) == 0 for n in range(5))
assert Rn[0] == 1 and sp.expand(Rn[1] - (1 - 2*G)) == 0 and sp.expand(Rn[2] - (4*G**2 - 4*G + 2)) == 0
p2 = sp.Poly(bspace, X)
amp2 = sp.expand(sum(c*Rn[m[0]] for m, c in zip(p2.monoms(), p2.coeffs())))
nm = data['nmhv_three_to_three']
if {k for k in nm if not k.startswith('_')}!={str(L) for L in range(1,5)} or any(
        set(nm[str(L)]['coefficients'])!={str(k) for k in range(2*L-1)} for L in range(1,5)):
    raise ValueError('Incomplete NMHV comparison data')
report2, raw2 = [], {}
for L in range(1, 5):
    pred = sp.expand(sp.expand(amp2.coeff(a, L - 1)).subs(G, lnd/2).subs(solved))
    for k in range(2*L - 2, -1, -1):             # including the constant, k = 0
        if str(k) not in nm[str(L)]['coefficients']:
            continue
        target = sp.sympify(nm[str(L)]['coefficients'][str(k)], locals=zz)
        raw2[(L, k)] = nm[str(L)]['coefficients'][str(k)]
        diff = sp.expand(pred.coeff(lnd, k) - target)
        if (L, k) == (1, 0):                     # the one-loop term defines K: not a prediction
            assert diff == 0 and target == 1
            report2.append((L, k, 'normalisation'))
            continue
        report2.append((L, k, 'PASS' if diff == 0 else f'FAIL: predicted {to_zeta(pred.coeff(lnd, k))}, erratum {to_zeta(target)}'))
print('NMHV, |Delta J_z| = 2 component, coefficients of pi i (1+w)/(1+w*) ln^k|delta| in rho^(L), no new constant:')
for L, k, s in report2:
    print(f'  L={L} k={k}: {s}')
n_pass2 = sum(1 for _, _, s in report2 if s == 'PASS')
n_checked2 = sum(1 for _, _, s in report2 if s != 'normalisation')

# ---------------------------------------------------------------- LL tower in closed form (Gaussian integral)
# -G int_0^1 dt exp(-2 a G^2 t^2) = -G sum_n (-2 a G^2)^n / (n! (2n+1))
LLclosed = sp.expand(-G*sum((-2*a*G**2)**n/(sp.factorial(n)*(2*n+1)) for n in range(NA + 1)))
ll_ok = all(sp.expand(LLclosed.coeff(a, L-1).coeff(G, 2*L-1) - A[L].coeff(G, 2*L-1)) == 0 for L in range(1, LMAX + 1))
print('LL tower equals -G int_0^1 dt exp(-2 a G^2 t^2) (Gaussian integral) through', LMAX, 'loops:', ll_ok)

# ---------------------------------------------------------------- DE, 2 -> 4 at two loops from R_6^(2) and Y (eq. 2.13)
d24 = data['two_to_four']
Li2 = lambda t: sp.polylog(2, t)
H2v = Li2(1 - v)         # H_2^v = H_{0,1}(1 - v) = Li_2(1 - v): DE write H^u_w for H_w(1 - u) (their eq. (3.9) is the
                         # value at v = 1, where every H^v vanishes); until 20 October 2026 this was read as Li_2(v)
Y = sp.sympify(d24['Y']['expr'], locals={'lnd': lnd, 'v': v, 'Li2': Li2, 'I': I, 'pi': pi})
R62 = sum(sp.sympify(c, locals={'I': I, 'pi': pi, 'H2v': H2v, 'v': v, **zz})*lnd**int(k) for k, c in d24['R6_2']['lnd_terms'].items())
gK = sp.sympify(d24['gamma_K']['expr'], locals={'a': a, 'G4': G4, **zz})
E24 = sp.expand(sp.series(sp.exp(a**2*R62 - gK/8*Y), a, 0, 3).removeO())      # R_6 starts at two loops
E24_1, E24_2 = sp.expand(E24.coeff(a, 1)), sp.expand(E24.coeff(a, 2))
# prediction: C_1 = -2 pi i, G = ln(delta)/2 + i pi/2; the one-loop constant r0 = z2 taken from 3 -> 3
pred24 = {L: sp.expand(sp.expand(A[L]*(2*pi*I)).subs(G, lnd/2 + I*pi/2).subs(solved)) for L in (1, 2)}
rows24 = []
for k in (1,):
    rows24.append(('1', k, sp.simplify(E24_1.coeff(lnd, k)), sp.simplify(pred24[1].coeff(lnd, k))))
for k in (3, 2, 1):
    rows24.append(('2', k, sp.simplify(E24_2.coeff(lnd, k)), sp.simplify(pred24[2].coeff(lnd, k))))
print('2 -> 4 kinematics (coefficient of ln^k delta in E^(L)): DE from eqs. (2.13), (2.32), (3.10) | prediction')
diffs24 = []
for L, k, de_c, pr in rows24:
    diff = sp.simplify(sp.expand(de_c - pr))
    diffs24.append(diff)
    print(f'  L={L} k={k}: DE = {de_c}  |  prediction = {pr}  |  difference = {diff}')
# every coefficient agrees, the two-loop single logarithm included: -2 i pi^3/3, independent of v (with H_2^v read
# as Li_2(v) the single logarithm would exceed the prediction by i pi [Li_2(1-v) - Li_2(v)]); the 2 -> 4 single
# logarithm fixes the finite part of the Glauber region of the kernel to zero (r0/theory/glauber.py)
single_log_24 = all(d == 0 for d in diffs24) and not rows24[-1][2].has(v)
print('2 -> 4 kinematics: all coefficients agree, the two-loop single logarithm is independent of v:', single_log_24)

n_pass = sum(1 for _, _, s, _ in report if s == 'PASS')
n_norm = sum(1 for _, _, s, _ in report if s == 'normalisation')
assert not any(s.startswith('FAIL') for _, _, s, _ in report)

# ---------------------------------------------------------------- what the vanishing three-loop ln^3 tests
# with the two-loop cusp coefficient g2 left free, the three-loop coefficient of G^3 (units of -C_1)
g2 = sp.Symbol('g2')
b3 = sp.expand(sp.series(sp.exp(-(4*a + g2*a**2)/8*X**2)*(1 + r0*a), a, 0, 3).removeO())
ln3 = sp.expand(transform(b3).coeff(a, 2).coeff(G, 3))
assert sp.expand(ln3 - sp.Rational(2, 3)*(r0 + g2/4)) == 0         # M_1 = -C_1 amp: C_1 (-2/3)(kappa_1 + gamma_K^(2)/4) G^3
assert sp.expand(ln3.subs({r0: solved[r0], g2: -4*z[2]})) == 0

# ---------------------------------------------------------------- the two components are related exactly
# d/dG of the nu = 0 generating function removes the factor 1/(2 alpha); 1/(1 - alpha) = int_0^oo ds e^{-s(1-alpha)}
# and e^{s alpha} shifts G by -s/2: M_2(G)/K = int_0^oo ds e^{-s} (d M_1/dG)(G - s/2)/C_1
sv = sp.Symbol('s')
dM1 = -sp.diff(amp, G)                            # d(M_1/C_1)/dG, amp being M_1/(-C_1)
shifted = sp.Poly(sp.expand(dM1.subs(G, G - sv/2)), sv)
rel = sum(c*sp.factorial(m[0]) for m, c in zip(shifted.monoms(), shifted.coeffs()))   # int_0^oo ds e^{-s} s^m = m!
components_related = sp.expand(rel - amp2) == 0
print('M_2(G)/K = int_0^oo ds e^{-s} (dM_1/dG)(G - s/2)/C_1 through', LMAX, 'loops:', components_related)

# ---------------------------------------------------------------- five-loop prediction for the |Delta J_z| = 2 component
pred5 = sp.expand(sp.expand(amp2.coeff(a, 4)).subs(G, lnd/2).subs(solved))
five = {k: str(to_zeta(pred5.coeff(lnd, k))) for k in range(8, -1, -1)}
print('five-loop |Delta J_z| = 2 prediction, coefficients of K ln^k|delta|:', five)
# the nu = 2 exponent is a plain exponential at leading logarithm: K exp(-2 a G^2)
LL2 = sp.expand(sp.series(sp.exp(-2*a*G**2), a, 0, NA + 1).removeO())
ll2_ok = all(sp.expand(LL2.coeff(a, L-1).coeff(G, 2*L-2) - amp2.coeff(a, L-1).coeff(G, 2*L-2)) == 0 for L in range(1, LMAX + 1))
print('nu = 2 LL tower equals exp(-2 a G^2) through', LMAX, 'loops:', ll2_ok)
out = {'loops_compared': LMAX, 'source': data['_source'],
       'three_to_three': [{'L': L, 'k': k, 'status': s, 'DE': raw.get((L, k), ''), 'DE_eq': data['three_to_three'][str(L)]['eq']}
                          for L, k, s, t in report],
       'n_checked': len(report), 'n_pass': n_pass, 'n_fixed': len(solved), 'n_norm': n_norm,
       'three_loop_ln3_coefficient': 'C_1 (-2/3) (kappa_1 + gamma_K^(2)/4) G^3, in the notation of the paper',
       'components_related': bool(components_related), 'nmhv_five_loop_prediction': five,
       'solved_constants': {str(k): str(to_zeta(val)) for k, val in solved.items()},
       'solved_order': [str(k) for k in solved],
       'four_loop_cusp_agrees': bool(g4_ok), 'LL_Gaussian_integral': bool(ll_ok),
       'nmhv': [{'L': L, 'k': k, 'status': s, 'DE': raw2[(L, k)], 'DE_eq': nm[str(L)]['eq']} for L, k, s in report2],
       'nmhv_n_checked': n_checked2, 'nmhv_n_pass': n_pass2, 'nmhv_LL_exponential': bool(ll2_ok),
       'two_to_four_two_loops': [{'L': L, 'k': k, 'DE': str(d), 'prediction': str(p)} for L, k, d, p in rows24],
       'two_to_four_single_log_difference': '0' if single_log_24 else 'other'}
if (any(s.startswith('FAIL') for _,_,s in report2) or
        any(s.startswith('FAIL') for _,_,s,_ in report) or
        not all((g4_ok,ll_ok,ll2_ok,components_related,single_log_24)) or
        n_pass+n_norm+len(solved)!=len(report) or n_pass2!=15 or n_checked2!=15):
    raise ValueError('Self-crossing coefficient comparison failed; results were not written')
OUT.write_text(json.dumps(out, indent=1))
print('wrote', OUT.name)
