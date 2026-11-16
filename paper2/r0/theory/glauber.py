"""The imaginary part of the one-loop DPS kernel: the Coulomb phase of the soft function and the Glauber region
(Paper 2, task 5; workbook step VII.15).  Writes glauber.json.  About a minute; no amplitude evaluations.

    python3 glauber.py

Result.  The kernel depends on the hard scales only through its rapidity logarithm (appendix E), whose argument is
Q_*^2 e^{-i pi}: the phase of -s - i0 of the two incoming parents.  The regions of appendix F are real and give the
real part of that logarithm; the imaginary part of the kernel is -i pi times its coefficient,

    Im K^(1)/K_0 = -2 pi i Dc R(eps),   R(eps) = e^{-eps gamma_E} Gamma(-eps) (mu^2 b^2/b0^2)^eps,
                 = 2 pi i Dc [1/eps + X_mu + eps (X_mu^2/2 + zeta_2/2) + O(eps^2)],

the pole 2 pi i Dc/eps of 2 I_K (check 11), the X_mu term the imaginary part of the bracket of eq. (5.3), and at
mu = mu_b no finite part, for the same reason as r_1 = 0 (step VI.5): R has no eps^0 term for b0 = 2 e^{-gamma_E}.

Where the phase comes from (parts 1-3 below):
 1  The eikonal (soft-function) exchange between the cross pair, 1/(l^- + i0) for the incoming n line that absorbs
    the gluon and 1/(-l^+ + i0) for the incoming nbar line that emits it: the l^+ residue leaves l^- > 0 only, with the
    eikonal pole at the end point l^- = 0 of the range.  With the eta regulator |l^+ - l^-|^-eta the integrand there
    is (l^-)^(eta-1), the i0 drops out and the soft function is real, eq. (F.3): the Glauber corner of the eikonal
    integrand gives the real soft function of eq. (F.3), with the endpoint treated before the
    regulator is removed.  (With the
    delta regulator the two i pi/2 of the two lines cancel, kernel/soft_function.py: also real.)
 2  In the two-loop DPS graph the daughters are off shell.  Use the exchange routing of kernel/model_kres.py:
    kbar = p1-k, r = p3-k-l and rbar = p2-p3+k+l. The first hard vertex receives k+l and r;
    the second receives kbar and rbar-l. At the DPS residues (k^- at the pole of kbar^2, k^+ at
    the pole of r^2) the two lines that carry the gluon momentum are, at leading power in the Glauber scaling
    l^+- ~ lambda^2, l_perp ~ lambda,
        (k + l)^2 -> x (l^- - sigma_1) + i0,  sigma_1 = k_perp^2/(1-x) + (k+l)_perp^2/x > 0,
        (rbar - l)^2 -> -x (l^+ + sigma_2) + i0,  sigma_2 = (k+l-q)_perp^2/(1-x) + (k-q)_perp^2/x > 0,
    eikonal lines with the light-cone energy denominators of the splittings, both positive for two incoming parents
    (the signs are fixed by the pinch: at the residue of kbar^2 the line k is below threshold).  The exact integral
        J = int dl^+ dl^- / [(l^+ l^- - A + i0)(l^- - sigma_1 + i0)(-l^+ - sigma_2 + i0)],   A = l_perp^2,
          = 2 pi i [ln(A/(sigma_1 sigma_2)) - i pi]/(A + sigma_1 sigma_2):
    the pole at l^- = sigma_1 now lies inside the range and its i0 gives the -i pi.  The logarithm is the rapidity
    logarithm (sigma_1 sigma_2 ~ mu_b^4/Q_*^2), and the -i pi accompanies it: ln(Q_*^2 e^{-i pi}/l_perp^2).  For an
    incoming-outgoing pair (sigma_2 < 0, the 3 -> 3 hexagon) the second pole enters the range too and the two i pi's
    cancel: J/(2 pi i) is real.  Zero bin: the eta-regulated soft and collinear functions carry the real part of the logarithm
    (their rapidity divergences cancel and leave ln(Q_*^2/mu_b^2)) and no phase (part 1), so the Glauber region adds
    -i pi times the coefficient of the rapidity logarithm with nothing to subtract.
 2b The same by regions, with the regulator (glauber_regions.py): J is expanded at fixed l_perp into its n-collinear
    (l^+ ~ 1), nbar-collinear (l^- ~ 1), soft (l^+- ~ lambda) and Glauber (l^+- ~ lambda^2) regions, each with the
    rapidity regulator of appendix F at its own scaling (|l^+|^-eta, |l^-|^-eta, |l^+ - l^-|^-eta, |l^+ - l^-|^-eta),
    all four in closed form, the zero bins scaleless: the 1/eta poles cancel (collinear -2, soft +2 in units of
    2 pi i/(A eta)), the sum tends to J as eta -> 0 up to the power correction sigma_1 sigma_2/A, and the Glauber region
    alone is the -i pi term, 2 pi^2/A: with 1/(x + i0) = PV - i pi delta it is the delta-delta piece, pi^2/A (all that a
    product regulator |l^+|^-eta |l^-|^-eta gives), plus the cross term of the two principal values, pi^2 (sigma_1 +
    sigma_2)^-eta/A, which the symmetric regulator couples.  The collinear and soft regions are purely imaginary (real
    in the normalisation of the kernel, where 2 pi i is the lowest-order factor).  For sigma_2 < 0 the Glauber region
    keeps its 2 pi^2/A and the collinear function of the outgoing parent, whose daughter pole lies inside its range,
    cancels it.  In b space, with the transverse integral in 2 - 2 eps dimensions and eta -> 0 before eps -> 0, the
    Glauber region is -i pi x 2 Dc e^{-eps gamma_E} Gamma(-eps) (mu^2 b^2/b0^2)^eps: part 0, no constant at mu = mu_b.
 3  The explicit two-loop graph of kernel/model_fast.py (scalar lines, the DPS numerator, the gluon between rbar and
    k, the l^+- integrals exact with their i0's): at fixed l_perp = lambda lhat, lambda^4 F = f1 ln lambda + f0 +
    O(lambda^2) with f1 = -2 (2 pi i)^3 N/(x^2 k^2 a l_perp^2) imaginary (model_numerics.py).  The Glauber region of
    this graph from its leading-power integrand, the DPS residues at l^+- = 0 (R^- -> 2 pi i/(x k^2), R^+ -> 2 pi i/(a x)),
    the static gluon -1/l_perp^2 and the two eikonal lines with sigma_1, sigma_2, is N (2 pi i)^2/(x^2 k^2 a) times
    2 pi^2/l_perp^2 (part 2b), and it equals the real part of the exact graph,
        Re f0 = G_graph = -(pi/2) Im f1,   i.e.  lambda^4 F = -(f1/2) ln(e^{-i pi}/lambda^2) + i (real),
    for eighteen configurations (three (x, khat), six lhat), extrapolated in lambda^2 from lambda = 0.004 and 0.002.
Tests (parts 4, 5):
 4  QCD: Im r0 = 0 within the fit uncertainties, from the stored fits of gg -> HH (scan and blind velocities),
    gg -> ZZ (scan and blind, all eight amplitudes) and gg -> AA (common/dps/results_aa.json, if present).
 5  Planar N = 4 in 2 -> 4 kinematics (s_36 > 0): the assembly of n4_hexagon.py continued with l_12, l_45 -> l - i pi:
    the imaginary single poles of h_1 + h_2 (-2 i pi/eps in units of a) and of -M6hat/2 (+ i pi/eps) cancel only with
    the kernel's + i pi/eps (2 pi i Dc/eps in a_s units), and the finite part is the 3 -> 3 form with X = X_r - i pi:
    no v, no scattering angle.  Dixon-Esterlis' two-loop single logarithm there, from E = exp[R_6 - gamma_K Y/8] with
    Y of eq. (2.32) and R_6^(2) of eq. (3.10), in which H^v_w = H_w(1-v) (their eq. (3.9) at v = 1 and
    "H^u_{3,1} = H_{0,0,1,1}(1-u)"), is -2 i pi^3/3, independent of v, and equals C_1[-2 G^3/3 - zeta_3/3 + kappa_1 G]
    with G = ln(delta)/2 + i pi/2, C_1 = -2 pi i and kappa_1 = zeta_2 + c_G, c_G the finite part of the Glauber
    region in units of a: the coefficient fixes c_G = 0.  (Reading H_2^v as Li_2(v), as de_data.json did until
    20 October 2026, produces the spurious v-dependent difference i pi [Li_2(1-v) - Li_2(v)].)  The ln^2 delta
    coefficient, -pi^2/2, tests Im X = -pi.
"""
import importlib.util, json, os, sys
import numpy as np
import sympy as sp
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
KERNEL = os.path.normpath(os.path.join(HERE, '..', '..', 'kernel'))
N4DIR = os.path.normpath(os.path.join(HERE, '..', '..', 'n4'))
DATA = os.path.normpath(os.path.join(HERE, '..', 'data'))
DPS = os.path.normpath(os.path.join(HERE, '..', '..', '..', 'common', 'dps'))

eps = sp.Symbol('epsilon')
z2, z3 = sp.zeta(2), sp.zeta(3)
I, pi = sp.I, sp.pi
ok = lambda e: sp.simplify(sp.expand(e)) == 0
out = {}


def report(n, text):
    print(f'{n}  {text}')


# ---------------------------------------------------------------------------------------- 0. the b-space result
def part0():
    Dc, Lb = sp.symbols('Delta_C X_mu', real=True)
    R = sp.exp(-eps*sp.EulerGamma)*sp.gamma(-eps)*sp.exp(eps*Lb)            # (mu^2 b^2/b0^2)^eps = e^{eps X_mu}
    imK = sp.expand(sp.series(-2*pi*I*Dc*R, eps, 0, 2).removeO())
    c = {j: sp.simplify(imK.coeff(eps, j)) for j in (-1, 0, 1)}
    assert ok(c[-1] - 2*pi*I*Dc) and ok(c[0] - 2*pi*I*Dc*Lb) and ok(c[1] - 2*pi*I*Dc*(Lb**2/2 + z2/2))
    at_mub = {j: sp.simplify(v.subs(Lb, 0)) for j, v in c.items()}
    assert at_mub[0] == 0
    # the pole of 2 I_K (check 11 of check_formulas.py): the imaginary single pole is 2 pi i Dc/eps
    CA, b0s = sp.symbols('C_A beta0', real=True)
    m2 = sp.Symbol('m2', positive=True)
    catani = lambda C, g, Q: -sp.exp(eps*sp.EulerGamma)/sp.gamma(1 - eps)*(C/eps**2 + g/eps)*sp.exp(I*pi*eps)*Q**(-eps)
    twoIK = 2*catani(CA, b0s/2, 1) - 2*2*catani(CA, b0s/2, m2)
    p1 = sp.series(sp.expand(twoIK), eps, 0, 0).removeO().coeff(eps, -1)
    assert ok(sp.im(p1) - 2*pi*CA)
    # the imaginary part of the bracket of eq. (5.3), Dc X_mu^2 - 2 Dc X X_mu + B X_mu, with Im X = -pi
    X, B = sp.symbols('X B', real=True)
    br = Dc*Lb**2 - 2*Dc*(X - I*pi)*Lb + B*Lb
    assert ok(sp.im(sp.expand(br)) - 2*pi*Dc*Lb)
    report(0, 'Im K^(1)/K_0 = -2 pi i Dc e^{-eps gamma} Gamma(-eps) (mu^2 b^2/b0^2)^eps = 2 pi i Dc [1/eps + X_mu + '
              'eps (X_mu^2/2 + zeta_2/2)]: pole of 2 I_K, X_mu term of the bracket of eq. (5.3), no eps^0 at mu = mu_b')
    out['b_space'] = {'Im K/K0 in units of 2 pi i Dc': {'eps^-1': str(c[-1]/(2*pi*I*Dc)), 'eps^0': str(c[0]/(2*pi*I*Dc)),
                                                         'eps^1': str(c[1]/(2*pi*I*Dc))},
                      'at mu = mu_b': {str(j): str(v) for j, v in at_mub.items()},
                      'finite part at mu = mu_b': 0, 'pole equals that of 2 I_K': True,
                      'X_mu term equals Im of the RG bracket': True}


# ---------------------------------------------------------------------------------------- 1. the eikonal end point
def part1():
    """the l^- integral of the on-shell eikonal exchange with the eta regulator is real; the delta regulator too"""
    mp.mp.dps = 30
    A = mp.mpf('0.7')
    # eta regulator at the gluon pole, l^+ = A/l^-: int_0^Lambda dl^- (l^-)^-1 |A/l^- - l^-|^-eta is real for any eta > 0
    eta = mp.mpf('0.3')
    f = lambda lm: (1/lm)*abs(A/lm - lm)**(-eta)
    val = mp.quad(f, [0, mp.sqrt(A), 10, mp.inf])
    assert mp.im(val) == 0
    # with a finite i0 instead of the regulator: int_0^Lambda dl^-/(l^- + i e0) = ln(Lambda/e0) - i pi/2: the phase
    # sits at the end point and is a property of the regulator, not of the exchange (step VI.5, soft_function.py)
    e0, Lam = mp.mpf('1e-6'), mp.mpf('1e3')
    g = mp.quad(lambda lm: 1/(lm + 1j*e0), [0, 10*e0, 1, Lam])
    assert abs(g - (mp.log(Lam/e0) - 1j*mp.pi/2)) < 1e-6
    report(1, 'on-shell eikonal lines: the eikonal pole sits at the end point l^- = 0 of the range; with the eta '
              'regulator the integrand is (l^-)^(eta-1) and real; with an i0 cutoff the phase -i pi/2 is the regulator\'s')
    out['eikonal_end_point'] = {'eta-regulated l^- integral imaginary part': float(mp.im(val)),
                                'i0-cutoff integral minus [ln(Lambda/e0) - i pi/2]': float(abs(g - (mp.log(Lam/e0) - 1j*mp.pi/2)))}


# ---------------------------------------------------------------------------------------- 2. the off-shell eikonal pair
def J_closed(A, s1, s2):
    """int dl^+ dl^- / [(l^+ l^- - A + i0)(l^- - s1 + i0)(-l^+ - s2 + i0)], the l^+ integral by residues (l^- > 0),
    the l^- integral in closed form: 2 pi i [ln(A/(s1 s2)) - i pi]/(A + s1 s2) for s1, s2 > 0, real for s1 s2 < 0"""
    A,s1,s2=map(mp.mpf,(A,s1,s2))
    if not all(mp.isfinite(v) for v in (A,s1,s2)) or A<=0 or s1*s2==0:
        raise ValueError('A must be positive and offsets finite and nonzero')
    if A+s1*s2==0:
        return 2j*mp.pi/A     # removable limit for opposite-sign offsets
    phase=mp.pi*(mp.sign(s1)+mp.sign(s2))/2
    return 2j*mp.pi*(mp.log(A/abs(s1*s2))-1j*phase)/(A+s1*s2)


def J_numeric(A, s1, s2, e0):
    """the same with the i0 as a finite e0 in the eikonal lines, l^+ by residues (exact: for l^- > 0 the gluon pole lies
    below and the eikonal pole -l^+ - s2 + i0 = 0 above the axis; for l^- < 0 both lie above), l^- by quadrature on
    (0, oo) with break points around the Lorentzian peaks of width e0"""
    f = lambda lm: 2j*mp.pi/((A + s2*lm - 1j*e0)*(lm - s1 + 1j*e0))
    pts = {mp.mpf(0), mp.mpf(1), mp.mpf(10), mp.mpf(100), mp.mpf(1000)}
    for p0 in ([s1] if s1 > 0 else []) + ([A/abs(s2)] if s2 < 0 else []):
        pts.add(p0)
        for k in range(0, 8):
            for sgn in (1, -1):
                v = p0 + sgn*e0*mp.mpf(10)**k
                if v > 0:
                    pts.add(v)
    return mp.quad(f, sorted(pts) + [mp.inf])


def part2():
    mp.mp.dps = 25
    worst = 0
    rows = []
    for A, s1, s2 in ((mp.mpf('0.7'), mp.mpf('0.2'), mp.mpf('0.5')), (mp.mpf('1.3'), mp.mpf('0.05'), mp.mpf('2.0')),
                      (mp.mpf('0.7'), mp.mpf('0.2'), mp.mpf('-0.5')), (mp.mpf('2.0'), mp.mpf('-0.3'), mp.mpf('0.9'))):
        # the l^- integral with the i0's as finite e0, at two values of e0, extrapolated linearly to e0 = 0
        v1, v2 = J_numeric(A, s1, s2, mp.mpf('1e-6')), J_numeric(A, s1, s2, mp.mpf('5e-7'))
        num = 2*v2 - v1
        cl = J_closed(A, s1, s2)
        d = abs(num - cl)/abs(cl)
        worst = max(worst, d)
        rows.append({'A': float(A), 'sigma1': float(s1), 'sigma2': float(s2), 'closed form': [float(mp.re(cl)), float(mp.im(cl))],
                     'numerical': [float(mp.re(num)), float(mp.im(num))], 'rel diff': float(d),
                     'phase relative to the log': 'ln(A/sigma1 sigma2) - i pi' if s1 > 0 and s2 > 0 else 'none (real)'})
    assert worst < 1e-9, worst
    report(2, f'off-shell eikonal pair: J = 2 pi i [ln(A/sigma1 sigma2) - i pi]/(A + sigma1 sigma2) for sigma1, sigma2 > 0 '
              f'(two incoming parents), real for sigma1 sigma2 < 0 (3 -> 3); closed form against quadrature to {float(worst):.0e}')
    out['off_shell_eikonal_pair'] = {'rows': rows, 'worst rel diff': float(worst)}


# ---------------------------------------------------------------------------------------- 2b. the regions with the regulator
def part2b():
    spec = importlib.util.spec_from_file_location('glauber_regions', os.path.join(HERE, 'glauber_regions.py'))
    GR = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(GR)
    res = GR.results()
    two_pi2 = 2*float(mp.pi)**2
    worst_sum = max(c['rel diff'] for c in res['cases'] if c['A/|s1 s2|'] >= 1e5)
    worst_G = max(abs(c['Glauber region x A, eta -> 0'][0]/two_pi2 - 1) for c in res['cases'])
    worst_other = max(abs(c['real part of C_n + C_nbar + S, eta -> 0']) for c in res['cases'] if c['s2'] > 0)
    worst_cancel = max(abs(c['real part of C_n + C_nbar + S, eta -> 0']/two_pi2 + 1) for c in res['cases'] if c['s2'] < 0)
    assert worst_sum < 1e-5 and worst_G < 1e-8 and worst_other < 1e-20 and worst_cancel < 1e-5, (worst_sum, worst_G, worst_other, worst_cancel)
    assert all(abs(c['pole coefficients x eta/(2 pi i/A) at the smallest eta (collinear, soft)'][0] + 2) < 0.1 and
               abs(c['pole coefficients x eta/(2 pi i/A) at the smallest eta (collinear, soft)'][1] - 2) < 0.1 for c in res['cases'])
    report('2b', f'the four eta-regulated regions of J (|l^+|^-eta, |l^-|^-eta, |l^+ - l^-|^-eta, |l^+ - l^-|^-eta; zero bins '
                 f'scaleless) sum to J as eta -> 0 to {worst_sum:.0e} at A/|sigma1 sigma2| = 5e5 (the power correction), the '
                 f'1/eta poles cancelling (collinear -2, soft +2); the Glauber region alone is the -i pi term 2 pi^2/A (to '
                 f'{worst_G:.0e}): pi^2/A from the delta functions of the two eikonal lines and pi^2/A from the cross term of '
                 f'their principal values; collinear + soft real parts {worst_other:.0e}; for sigma2 < 0 the collinear '
                 f'function of the outgoing parent cancels it (to {worst_cancel:.0e})')
    res |= {'worst rel diff of the sum from J at A/|s1 s2| = 5e5': worst_sum, 'worst rel deviation of G from 2 pi^2/A': worst_G,
            'largest real part of C_n + C_nbar + S (s2 > 0)': worst_other, 'worst rel deviation of the 3 -> 3 cancellation': worst_cancel,
            'in b space': 'G = -i pi x 2 Dc e^{-eps gamma_E} Gamma(-eps) (mu^2 b^2/b0^2)^eps (eta -> 0 before eps -> 0): part 0'}
    out['eta_regulated_regions'] = res


# ---------------------------------------------------------------------------------------- 3. the explicit two-loop graph
def part3():
    spec = importlib.util.spec_from_file_location('model_fast', os.path.join(KERNEL, 'model_fast.py'))
    MF = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(MF)

    def F(lam, lt, x, KH, QH):
        kx, ky = lam*KH[0], lam*KH[1]; qx, qy = lam*QH[0], lam*QH[1]
        L = max(1000.0, 100.0*(lt[0]**2 + lt[1]**2))
        return MF.lplus_integral(lt[0], lt[1], kx, ky, qx, qy, x, 1e-4*lam*lam, L=L)

    rows, worst, worst_raw, worst_slope = [], 0.0, 0.0, 0.0
    lams = (0.004, 0.002)
    for x, KH, QH in ((0.75, (0.6, 0.3), (1.0, 0.0)), (0.4, (-0.5, 0.8), (1.0, 0.0)), (0.9, (0.3, -0.2), (1.0, 0.0))):
        for lh in ((0.2, -0.5), (-3.0, 2.0), (0.15, 0.1), (0.9, 0.4), (-0.3, -1.1), (2.0, 0.3)):
            k2 = KH[0]**2 + KH[1]**2
            a = (KH[0] + lh[0] - QH[0])**2 + (KH[1] + lh[1] - QH[1])**2
            N = KH[0]*(KH[0] + lh[0] - QH[0]) + KH[1]*(KH[1] + lh[1] - QH[1])
            l2 = lh[0]**2 + lh[1]**2
            f1 = -2*(2j*np.pi)**3*N/(x**2*k2*a*l2)                  # the eikonal coefficient of ln lambda, imaginary
            # the Glauber region of the graph from its leading-power integrand: N R^- R^+/D7 with R^- -> 2 pi i/(x k^2),
            # R^+ -> 2 pi i/(a x), D7 -> -l_perp^2 and the l^+- integral of the two eikonal lines, 2 pi^2/l_perp^2 x (-1)
            # (part 2b: G = (1/A)[PP + pi^2 DD] with the static gluon -1/A); lambda^4 cancels (N, k^2, a, l^2 ~ lambda^2)
            G_graph = (N*(2j*np.pi)**2/(x**2*k2*a)*2*np.pi**2/l2).real
            assert abs(G_graph/(-np.pi/2*f1.imag) - 1) < 1e-14          # the same number as -(pi/2) Im f1, by the algebra of part 2b
            target = G_graph                                        # Re f0 = G_graph = -(pi/2) Im f1
            v = [lam**4*F(lam, (lam*lh[0], lam*lh[1]), x, KH, QH) for lam in lams]
            rich = (4*v[1] - v[0])/3                                # the real part converges like lambda^2
            slope = ((v[1] - v[0])/np.log(lams[1]/lams[0])).imag   # the coefficient of ln lambda (converges like lambda)
            dev = rich.real/target - 1
            worst = max(worst, abs(dev)); worst_raw = max(worst_raw, abs(v[1].real/target - 1))
            worst_slope = max(worst_slope, abs(slope/f1.imag - 1))
            rows.append({'x': x, 'khat': KH, 'qhat': QH, 'lhat': lh, 'Im f1 (eikonal)': f1.imag,
                         'Im f1 (slope, lambda = 0.004, 0.002)': slope,
                         'Re lambda^4 F at lambda = 0.004, 0.002': [v[0].real, v[1].real],
                         'extrapolated Re f0': rich.real, 'Glauber region from the leading-power integrand': G_graph,
                         '-(pi/2) Im f1': -np.pi/2*f1.imag, 'rel deviation': dev})
    assert worst < 1e-3 and worst_slope < 0.1, (worst, worst_slope)
    report(3, f'explicit two-loop graph: Re f0 = Glauber region of the graph from its leading-power integrand (= -(pi/2) Im f1) '
              f'at {len(rows)} configurations, to {worst:.0e} (extrapolated in lambda^2 from lambda = 0.004, 0.002; raw '
              f'{worst_raw:.0e} at 0.002); the slope d(lambda^4 F)/d ln lambda is the eikonal coefficient to {worst_slope:.0e} '
              f'(converges like lambda)')
    out['explicit_two_loop_graph'] = {'statement': 'lambda^4 F = f1 [ln lambda + i pi/2] + i (real) = -(f1/2) ln(e^{-i pi}/lambda^2) + i (real); '
                                                   'Re f0 = N (2 pi i)^2/(x^2 k^2 a) x 2 pi^2/l_perp^2, the Glauber region of the graph',
                                      'rows': rows, 'worst rel deviation of Re f0 from -(pi/2) Im f1': worst,
                                      'worst at lambda = 0.002 without extrapolation': worst_raw,
                                      'worst rel deviation of the slope from the eikonal coefficient': worst_slope}


# ---------------------------------------------------------------------------------------- 4. QCD: Im r0 from the data
def part4():
    rows = []

    def add(process, tag, comp, label, value, unc):
        rows.append({'process': process, 'where': tag, 'component': comp, 'colour': label,
                     'Im r0': float(value), 'fit uncertainty': float(unc), 'Im r0 over uncertainty': float(abs(value)/unc) if unc > 0 else None})

    hh = json.load(open(os.path.join(DATA, 'r0_from_data.json')))
    for b, rec in hh.items():
        for nu in ('nu0', 'nu2'):
            for key, lab in (('r0_N', 'C_A'), ('r0_nf', 'T_F n_f')):
                v = rec[nu][key]
                add('gg -> HH', f'beta = {b}', nu, lab, float(v[1]), float(v[2]))
    for r in json.load(open(os.path.join(DATA, 'blind_test_hh.json'))):
        add('gg -> HH', f'beta = {r["beta"]} (blind)', r['component'], 'C_A', r['C_A'][1], r['shift_C_A'])
        add('gg -> HH', f'beta = {r["beta"]} (blind)', r['component'], 'T_F n_f', r['n_f'][1], r['shift_n_f'])
    zz = json.load(open(os.path.join(DATA, 'r0_from_zz.json')))
    for b, rec in zz.items():
        for nu in ('nu0', 'nu2'):
            for hel, hr in rec[nu].items():
                for key in ('2l_CF', '2l_CA'):
                    v = hr[key]
                    add('gg -> ZZ', f'beta = {b}, {hel}', nu, key[3:], float(v[1]), float(v[2]))
    for r in json.load(open(os.path.join(DATA, 'blind_test_zz.json'))):
        for key in ('2l_CF', '2l_CA'):
            add('gg -> ZZ', f'beta = {r["beta"]}, {r["hel"]} (blind)', r['component'], key[3:], r[key]['r0'][1], r[key]['shift'])
    aa_file = os.path.join(DPS, 'results_aa.json')
    aa_present = os.path.exists(aa_file)
    if aa_present:
        aa = json.load(open(aa_file))
        for b, rec in aa['betas'].items():
            for nu in ('nu0', 'nu2'):
                for c in ('N', 'TFnf'):
                    e = rec['r0'][nu][c]
                    add('gg -> AA', f'beta = {b}', nu, c, e['imag'], float(e['shift']))
    worst = max(r['Im r0'] and abs(r['Im r0']) for r in rows)
    ratio = max(r['Im r0 over uncertainty'] for r in rows if r['Im r0 over uncertainty'] is not None)
    assert worst < 1e-8 and ratio < 1.0, (worst, ratio)
    procs = sorted({r['process'] for r in rows})
    report(4, f'QCD: Im r0 = 0 in {", ".join(procs)}: largest |Im r0| = {worst:.0e} over {len(rows)} fitted values, '
              f'at most {ratio:.2f} of the fit uncertainty' + ('' if aa_present else ' (results_aa.json not present)'))
    out['qcd_imaginary_parts'] = {'rows': rows, 'largest |Im r0|': worst, 'largest |Im r0|/uncertainty': ratio,
                                  'processes': procs, 'gg -> AA included': aa_present}


# ---------------------------------------------------------------------------------------- 5. planar N = 4, 2 -> 4
def part5():
    # (a) the hexagon assembly of n4_hexagon.symbolic continued to s_36 > 0: l_12, l_45 -> l - i pi, X = X_r - i pi
    Xr, lx, l1x, ly, l1y, lT1, lT2, cG = sp.symbols('X_r lnx ln1mx lny ln1my lamT1 lamT2 c_G', real=True)
    ipi = I*pi
    X = Xr - ipi
    lS = {1: Xr - lx - ly - ipi, 2: Xr + lT1 - l1x, 3: Xr + lT2 - lx, 4: Xr - l1x - l1y - ipi, 5: Xr + lT2 - ly, 6: Xr + lT1 - l1y}
    lT = {1: Xr + lT1, 2: Xr + lT2}
    pole = lambda l: -sp.exp(-eps*l)/eps**2
    h = {k: pole(lS[s]) + pole(lT[k]) + sp.Rational(1, 2)*(lS[s] - lT[k])**2 + 4*z2 for k, s in ((1, 1), (2, 4))}
    Q = (sp.Rational(1, 2)*sum(lS[i]**2 for i in lS) - sum(lS[i]*lS[i % 6 + 1] for i in lS)
         + sp.Rational(1, 2)*sum(lS[i]*lS[(i + 2) % 6 + 1] for i in lS))
    M6hat = sum(pole(lS[i]) for i in lS) + Q + 9*z2
    hard = sp.expand(sp.series(h[1] + h[2] - M6hat/2, eps, 0, 1).removeO())
    im_pole_hard = sp.simplify(sp.im(hard.coeff(eps, -1)))
    assert ok(im_pole_hard + pi)                                     # -i pi/eps from the hard functions and M6hat
    kernel = 1/eps**2 + (-X)/eps + cG                                 # L_* = -X, the Glauber finite part c_G (units of a)
    total = sp.expand(sp.series(sp.expand(h[1] + h[2] - M6hat/2 + kernel), eps, 0, 1).removeO())
    assert sp.simplify(total.coeff(eps, -2)) == 0 and sp.simplify(total.coeff(eps, -1)) == 0
    Bfin = sp.expand(total.coeff(eps, 0))
    B33 = -X**2/2 + lx*l1x + ly*l1y + sp.Rational(7, 2)*z2 + cG
    assert ok(Bfin - B33) and not Bfin.has(lT1) and not Bfin.has(lT2)
    # without the kernel's i pi/eps the imaginary pole survives
    total_re = sp.expand(sp.series(sp.expand(h[1] + h[2] - M6hat/2 + 1/eps**2 - Xr/eps), eps, 0, 1).removeO())
    assert ok(sp.im(total_re.coeff(eps, -1)) + pi)
    # (b) Dixon-Esterlis, 2 -> 4: E = exp[R_6 - gamma_K Y/8], Y and R_6^(2) from de_data.json with H^v_w = H_w(1-v)
    de = json.load(open(os.path.join(N4DIR, 'de_data.json')))['two_to_four']
    lnd, v, a = sp.symbols('lnd v a')
    Li2 = lambda t: sp.polylog(2, t)
    H2v = Li2(1 - v)                                                  # H_2^v = H_{0,1}(1 - v) = Li_2(1 - v)
    Y = sp.sympify(de['Y']['expr'], locals={'lnd': lnd, 'v': v, 'Li2': Li2, 'I': I, 'pi': pi})
    R62 = sum(sp.sympify(c, locals={'I': I, 'pi': pi, 'H2v': H2v, 'v': v, 'z2': z2, 'z3': z3})*lnd**int(k)
              for k, c in de['R6_2']['lnd_terms'].items())
    gK = 4*a - 4*z2*a**2
    E2 = sp.expand(sp.series(sp.exp(a**2*R62 - gK/8*Y), a, 0, 3).removeO().coeff(a, 2))
    de_coeff = {k: sp.simplify(E2.coeff(lnd, k)) for k in (3, 2, 1)}
    # the structure: C_1 [-2 G^3/3 - zeta_3/3 + kappa_1 G], G = ln(delta)/2 + i pi/2, C_1 = -2 pi i, kappa_1 = zeta_2 + c_G
    G = lnd/2 + ipi/2
    C1 = -2*pi*I
    pred = sp.expand(C1*(-sp.Rational(2, 3)*G**3 - z3/3 + (z2 + cG)*G))
    pred_coeff = {k: sp.simplify(pred.coeff(lnd, k)) for k in (3, 2, 1)}
    assert ok(de_coeff[3] - pred_coeff[3]) and ok(de_coeff[2] - pred_coeff[2])
    diff1 = sp.simplify(sp.expand(de_coeff[1] - pred_coeff[1]))
    sol = sp.solve(diff1, cG)
    assert len(sol) == 1 and sol[0] == 0 and not de_coeff[1].has(v), (diff1, sol)
    # the misreading H_2^v = Li_2(v) would give the v-dependent difference i pi [Li_2(1-v) - Li_2(v)]
    R62_wrong = sum(sp.sympify(c, locals={'I': I, 'pi': pi, 'H2v': Li2(v), 'v': v, 'z2': z2, 'z3': z3})*lnd**int(k)
                    for k, c in de['R6_2']['lnd_terms'].items())
    E2w = sp.expand(sp.series(sp.exp(a**2*R62_wrong - gK/8*Y), a, 0, 3).removeO().coeff(a, 2))
    spurious = sp.simplify(sp.expand(E2w.coeff(lnd, 1) - pred_coeff[1].subs(cG, 0)))
    assert ok(spurious - I*pi*(Li2(1 - v) - Li2(v)))
    report(5, 'planar N = 4, 2 -> 4: poles cancel only with the kernel\'s + i pi/eps; finite part -X^2/2 + ln x ln(1-x) + '
              'ln y ln(1-y) + 7 zeta_2/2 + c_G with X = X_r - i pi, no v, no angles; Dixon-Esterlis single logarithm '
              f'{de_coeff[1]} = structure with kappa_1 = zeta_2 and c_G = 0 (ln^2: {de_coeff[2]}, ln^3: {de_coeff[3]})')
    out['n4_two_to_four'] = {'imaginary single pole of h1 + h2 - M6hat/2 (units of a)': str(im_pole_hard) + ' i',
                             'poles cancel with the kernel L_*/eps, L_* = -X_r + i pi': True,
                             'finite part': str(B33), 'depends on t1, t2 or v': False,
                             'DE coefficients of ln^k delta in E^(2) (H^v_w = H_w(1-v))': {str(k): str(de_coeff[k]) for k in (3, 2, 1)},
                             'structure (kappa_1 = zeta_2 + c_G)': {str(k): str(pred_coeff[k]) for k in (3, 2, 1)},
                             'c_G from the single logarithm': str(sol[0]),
                             'spurious difference with H_2^v read as Li_2(v)': str(spurious)}


def main():
    part0(); part1(); part2(); part2b(); part3(); part4(); part5()
    out['summary'] = ('Coulomb phase + Glauber region = -i pi x (coefficient of the rapidity logarithm) = '
                      '2 pi i Dc/eps + 0 eps^0 at mu = mu_b; Im r0 = 0 in QCD; kappa_1 = zeta_2 with c_G = 0 in 2 -> 4')
    json.dump(out, open(os.path.join(HERE, 'glauber.json'), 'w'), indent=1)
    print('wrote glauber.json')
    print('passed')


if __name__ == '__main__':
    main()
