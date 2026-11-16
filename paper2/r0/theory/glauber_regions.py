#!/usr/bin/env python3
"""The eta-regulated regions of the cross-pair exchange (step 1 of the Glauber-region check).

The scalar integral of the two light-cone components, with A = l_perp^2 and the light-cone energy
denominators s1, s2 of the two daughters (both positive for two incoming parents),

    I = int dl+ dl- / [(l+ l- - A + i0)(l- - s1 + i0)(-l+ - s2 + i0)]
      = 2 pi i [ln(A/(s1 s2)) - i pi] / (A + s1 s2)                          (exact; eq. (8) of the Letter)

is expanded by regions with the rapidity regulator of the paper, nu = 1, at each region's scaling:
  n-collinear    (l+ ~ 1, l- ~ lambda^2): gluon exact, k exact (s1 kept), rbar -> 1/(-l+ + i0), |l+|^-eta
  nbar-collinear (l- ~ 1, l+ ~ lambda^2): gluon exact, rbar exact (s2 kept), k -> 1/(l- + i0), |l-|^-eta
  soft           (l+, l- ~ lambda):       gluon exact, both eikonal, |l+ - l-|^-eta
  Glauber        (l+, l- ~ lambda^2):     gluon static, -1/A, both offsets kept, |l+ - l-|^-eta
with the zero bins (a region's integrand at a smaller region's scaling) scaleless. All four are closed forms:
  C_n   = -(2 pi i/A) (s1/A)^eta pi/sin(pi eta)                    (residue in l-, then int_0^inf du u^(-1-eta)/(A + s1 u))
  C_nb  = (2 pi i) int_0^inf dl (l)^(-1-eta)/(A + s2 l - i0)       (= C_n with s1 -> s2 for s2 > 0; for s2 < 0 the pole
                                                                      at l = A/|s2| lies inside the range)
  S     = -4 [Y_s R_s + Y_c R_c], in light-cone polar coordinates l+ = +-rho e^y, l- = +-rho e^-y: the regulator
          is (2 rho |sinh y|)^-eta in the quadrants with l+ l- > 0 and (2 rho cosh y)^-eta in the others, and the
          integrand factorises; Y_s = Gamma(eta/2)Gamma(1-eta)/Gamma(1-eta/2), Y_c = 2^-eta sqrt(pi) Gamma(eta/2)/Gamma((1+eta)/2),
          R_s = A^(-1-eta/2) [pi cot(pi eta/2) - i pi]/2, R_c = -A^(-1-eta/2) pi/(2 sin(pi eta/2))
  G     = (1/A) [PP + pi^2 DD + i pi (P1 D2 - D1 P2)] from 1/(x + i0) = PV - i pi delta:
          DD = |s1 + s2|^-eta, PP = pi^2 |s1 + s2|^-eta (the cross term of the two principal values, which the
          symmetric regulator couples; a product regulator |l+|^-eta |l-|^-eta gives PP = O(eta^2)),
          P1 D2 - D1 P2 = -2 pi tan(pi eta/2) sign(s1 + s2) |s1 + s2|^-eta.
The checks: the 1/eta poles cancel (C_n + C_nb -> -2, S -> +2 in units of 2 pi i/(A eta)); the sum of the four
regions tends to I as eta -> 0 up to the power correction s1 s2/A; the real part of the sum, 2 pi^2/A, is the
Glauber region alone, the collinear and soft regions being purely imaginary (real in the normalisation of the
kernel, where the overall 2 pi i is the lowest-order factor); for s2 < 0 (an incoming and an outgoing parent)
the Glauber region keeps its 2 pi^2/A and the collinear function of the outgoing parent cancels it, so that
the sum is purely imaginary, as the exact integral.

    python3 glauber_regions.py        (seconds; mpmath only)
"""
import mpmath as mp

mp.mp.dps = 30
PI = mp.pi


def C_n(A, s1, eta):
    return -(2j*PI/A)*(s1/A)**eta*PI/mp.sin(PI*eta)


def C_nb(A, s2, eta):
    if s2 > 0:
        return -(2j*PI/A)*(s2/A)**eta*PI/mp.sin(PI*eta)
    # s2 < 0: 2 pi i int_0^inf dl l^(-1-eta) / (A - |s2| l - i0), with 1/(A - |s2| l - i0) = -(1/|s2|)[PV 1/(l - A/|s2|) - i pi delta(l - A/|s2|)]:
    # the PV part by analytic continuation (t = |s2| l/A: PV int_0^inf t^(-1-eta)/(1 - t) dt = pi cot(-pi eta)), the delta part at l = A/|s2|
    a = abs(s2)
    pv = -(1/A)*(a/A)**eta*PI/mp.tan(PI*eta)
    dl = 1j*PI*(1/a)*(A/a)**(-1 - eta)
    return 2j*PI*(pv + dl)


def S(A, eta):
    Ys = mp.gamma(eta/2)*mp.gamma(1 - eta)/mp.gamma(1 - eta/2)
    Yc = 2**(-eta)*mp.sqrt(PI)*mp.gamma(eta/2)/mp.gamma((1 + eta)/2)
    Rs = mp.mpf(1)/2*A**(-1 - eta/2)*(PI/mp.tan(PI*eta/2) - 1j*PI)
    Rc = -mp.mpf(1)/2*A**(-1 - eta/2)*PI/mp.sin(PI*eta/2)
    return -4*(Ys*Rs + Yc*Rc)


def G(A, s1, s2, eta):
    s = s1 + s2
    PP = PI**2*abs(s)**(-eta)
    DD = abs(s)**(-eta)
    cross = -2*PI*mp.tan(PI*eta/2)*mp.sign(s)*abs(s)**(-eta)
    return (PP + PI**2*DD + 1j*PI*cross)/A


def exact(A, s1, s2):
    # for s1 s2 < 0 the logarithm's argument is negative: ln(A/(s1 s2)) - i pi = ln(A/|s1 s2|) with the i pi cancelled
    A,s1,s2=map(mp.mpf,(A,s1,s2))
    if not all(mp.isfinite(v) for v in (A,s1,s2)) or A<=0 or s1*s2==0:
        raise ValueError('A must be positive and offsets finite and nonzero')
    if A+s1*s2==0:
        return 2j*PI/A
    phase=PI*(mp.sign(s1)+mp.sign(s2))/2
    return 2j*PI*(mp.log(A/abs(s1*s2))-1j*phase)/(A+s1*s2)


def exact_numeric(A, s1, s2):
    # the l- integral after the l+ residue, 2 pi i int_0^inf dl / [(l - s1 + i0)(A + s2 l - i0)], by quadrature with
    # the poles treated as PV + delta terms (a check of the closed form, both signs of s2)
    f = lambda l: 1/((l - s1)*(A + s2*l))
    pts = sorted({s1} | ({A/abs(s2)} if s2 < 0 else set()))
    d = mp.mpf('1e-12')
    segs, prev = [], mp.mpf(0)
    for p in pts:
        segs.append([prev, (prev + p)/2, p - d]); prev = p + d
    segs.append([prev, prev + 1, mp.inf])
    val = sum(mp.quad(f, s) for s in segs)
    val += -1j*PI/(A + s2*s1)                               # delta at l = s1
    if s2 < 0:
        val += 1j*PI/(abs(s2)*(A/abs(s2) - s1))             # delta at l = A/|s2| (from A + s2 l - i0)
    return 2j*PI*val


def extrapolate(etas, vals):
    """polynomial extrapolation to eta = 0 (the regions are Laurent series in eta without logarithms)"""
    n = len(etas)
    return sum(v*mp.fprod((0 - etas[j])/(etas[i] - etas[j]) for j in range(n) if j != i) for i, v in enumerate(vals))


def results():
    """the checks as a dictionary (glauber.py, part 2b): for s2 > 0 and s2 < 0, at two values of A/|s1 s2|, the
    extrapolated sum of the four regions against the exact integral, and the Glauber region alone"""
    s1 = mp.mpf(1)
    out = {'regulator': 'nu^eta |l+|^-eta (n-collinear), |l-|^-eta (nbar-collinear), |l+ - l-|^-eta (soft, Glauber); zero bins scaleless',
           'closed forms': {'C_n': '-(2 pi i/A) (s1/A)^eta pi/sin(pi eta)',
                            'C_nbar': 'the same with s2 (s2 > 0); for s2 < 0 the daughter pole at l = A/|s2| lies inside the range and adds -2 pi^2/A',
                            'S': '-4 [Y_s R_s + Y_c R_c] in light-cone polar coordinates (see the docstring)',
                            'G': '(1/A) [PP + pi^2 DD + i pi (P1 D2 - D1 P2)], PP = pi^2 |s1+s2|^-eta, DD = |s1+s2|^-eta'},
           'cases': []}
    etas = [mp.mpf(1)/2**k for k in range(6, 10)]
    for s2 in (mp.mpf(2), mp.mpf(-2)):
        for A in (mp.mpf(10)**4, mp.mpf(10)**6):
            ex = exact(A, s1, s2)
            sums, Gs, poles = [], [], []
            for eta in etas:
                parts = (C_n(A, s1, eta), C_nb(A, s2, eta), S(A, eta), G(A, s1, s2, eta))
                sums.append(sum(parts)*A); Gs.append(parts[3]*A)
                poles.append((mp.im((parts[0] + parts[1])*A*eta)/(2*PI), mp.im(parts[2]*A*eta)/(2*PI)))
            tot, g0 = extrapolate(etas, sums), extrapolate(etas, Gs)
            out['cases'].append({'s1': float(s1), 's2': float(s2), 'A': float(A), 'A/|s1 s2|': float(A/abs(s1*s2)),
                                 'exact x A': [float(mp.re(ex*A)), float(mp.im(ex*A))],
                                 'sum of the four regions x A, eta -> 0': [float(mp.re(tot)), float(mp.im(tot))],
                                 'rel diff': float(abs(tot - ex*A)/abs(ex*A)),
                                 'Glauber region x A, eta -> 0': [float(mp.re(g0)), float(mp.im(g0))],
                                 'real part of C_n + C_nbar + S, eta -> 0': float(mp.re(tot - g0)),
                                 'pole coefficients x eta/(2 pi i/A) at the smallest eta (collinear, soft)': [float(poles[-1][0]), float(poles[-1][1])]})
    return out


def main():
    s1 = mp.mpf(1)
    for s2 in (mp.mpf(2), mp.mpf(-2)):
        print(f'\n=== s1 = {s1}, s2 = {s2}: {"two incoming parents" if s2 > 0 else "an incoming and an outgoing parent (3 -> 3)"}')
        for A in (mp.mpf(10)**4, mp.mpf(10)**6):
            ex = exact(A, s1, s2)
            exn = exact_numeric(A, s1, s2)
            print(f'A/|s1 s2| = {float(A/abs(s1*s2)):.0e}:  exact x A = {mp.nstr(ex*A, 12)}  (quadrature: {mp.nstr(exn*A, 12)})')
            etas = [mp.mpf(1)/2**k for k in range(6, 10)]
            sums, Gs = [], []
            for eta in etas:
                parts = (C_n(A, s1, eta), C_nb(A, s2, eta), S(A, eta), G(A, s1, s2, eta))
                sums.append(sum(parts)*A); Gs.append(parts[3]*A)
                print(f'   eta = {float(eta):8.5f}: poles x eta/(2 pi i/A): C_n + C_nb = {mp.nstr(mp.im((parts[0] + parts[1])*A*eta)/(2*PI), 7)}, '
                      f'S = {mp.nstr(mp.im(parts[2]*A*eta)/(2*PI), 7)};  sum x A = {mp.nstr(sums[-1], 12)}')
            tot, g0 = extrapolate(etas, sums), extrapolate(etas, Gs)
            print(f'   eta -> 0: sum x A = {mp.nstr(tot, 12)}   exact x A = {mp.nstr(ex*A, 12)}   '
                  f'rel diff {mp.nstr(abs(tot - ex*A)/abs(ex*A), 2)}  (power correction |s1 s2|/A = {float(abs(s1*s2)/A):.0e})')
            print(f'             Glauber region x A = {mp.nstr(g0, 10)} (2 pi^2 = {mp.nstr(2*PI**2, 10)});  '
                  f'real part of the sum = {mp.nstr(mp.re(tot), 10)};  real part of C_n + C_nb + S = {mp.nstr(mp.re(tot - g0), 4)}')


if __name__ == '__main__':
    main()
