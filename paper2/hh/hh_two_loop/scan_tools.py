# Shared helpers for the colour/scale scan of the class-B gg -> HH remainders
# (scan_colour_scale.py). Used by all analysis scripts in this directory.
#
# Every point file holds A = i Gauge/(4 M2), M2 = 1, for the four colour terms
#     Nf2:       coefficient of n_f^2
#     Nf_over_N: coefficient of n_f/N
#     N_Nf:      coefficient of N n_f
#     N2:        coefficient of N^2
# at L = ln(m_H^2/mu_R^2) = -1, 0, +1. Each is a quadratic in L, f = f0 + f1 L + f2 L^2.
#
# Fits: at fixed beta the forward remainder is a polynomial in G plus power
# corrections of relative order p_T^2 (times powers of G). Two models:
#   fits():      the leading-power polynomial alone, interpolated through deg+1
#                consecutive deltas;
#   fit_model(): the polynomial plus p_T^2 (d_0 + d_1 G + ... + d_pdeg G^pdeg),
#                interpolated through deg+pdeg+2 consecutive deltas.
# Window 0 is the one ending at the smallest delta. Interpolation leaves no
# residual by construction, so it is not an error estimate; the change between
# window 0 and window 1 (shifted one step up in delta) is quoted instead.
import json, pathlib
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 50
COLOURS = ('N2', 'N_Nf', 'Nf_over_N', 'Nf2')
COLOUR_WEIGHT = {'N2': lambda N, nf: N**2, 'N_Nf': lambda N, nf: N*nf,
                 'Nf_over_N': lambda N, nf: mp.mpf(nf)/N, 'Nf2': lambda N, nf: nf**2}


def load(points_dir):
    """Return {(gauge, beta): [records sorted by decreasing delta]} using the
    regular-precision points only (the extra-digit control points are returned separately)."""
    pts = [json.loads(p.read_text()) for p in sorted(pathlib.Path(points_dir).glob('*.json'))]
    regular, control = {}, []
    for r in pts:
        k = r['delta_exp']
        if r['dps'] != 30+3*k:
            control.append(r)
            continue
        regular.setdefault((r['gauge'], r['beta']), []).append(r)
    for v in regular.values():
        v.sort(key=lambda r: r['delta_exp'])
    return regular, control


def cplx(pair):
    return mp.mpc(mp.mpf(pair[0]), mp.mpf(pair[1]))


def L_parts(rec, colour):
    """(f0, f1, f2) of the quadratic in L for one colour term at one point."""
    f = {int(L): cplx(v) for L, v in rec['A'][colour].items()}
    f0 = f[0]
    f1 = (f[1]-f[-1])/2
    f2 = (f[1]+f[-1])/2-f0
    return f0, f1, f2


def G_of(rec):
    return cplx(rec['G'])


def window_fit(Gs, vals, deg):
    """Coefficients c_0..c_deg of sum_j c_j G^j through deg+1 points (exact interpolation)."""
    n = deg+1
    assert len(Gs) == n
    M = mp.matrix(n, n)
    for i, g in enumerate(Gs):
        for j in range(n):
            M[i, j] = g**j
    sol = mp.lu_solve(M, mp.matrix(vals))
    return [sol[j] for j in range(n)]


def fits(records, colour, part, deg):
    """Window fits for one colour term and one power of L, smallest-delta window first.
    Returns (coefficients of the smallest window, coefficients of the next window up)."""
    Gs = [G_of(r) for r in records]
    vals = [L_parts(r, colour)[part] for r in records]
    n = deg+1
    best = window_fit(Gs[-n:], vals[-n:], deg)
    nxt = window_fit(Gs[-n-1:-1], vals[-n-1:-1], deg) if len(Gs) > n else None
    return best, nxt


def fit_values(records, vals, deg, pdeg, window=0):
    """Leading coefficients c_0..c_deg of sum_j c_j G^j + p_T^2 sum_k d_k G^k, interpolated
    through deg+pdeg+2 consecutive points of `records` (with values `vals`); window counts
    up from the smallest delta. Returns None if there are not enough points."""
    n = deg+pdeg+2
    end = len(records)-window
    if end-n < 0:
        return None
    M = mp.matrix(n, n)
    v = mp.matrix(n, 1)
    for i in range(n):
        r = records[end-n+i]
        G = G_of(r)
        p = mp.mpf(r['pT2'])
        for j in range(deg+1):
            M[i, j] = G**j
        for k in range(pdeg+1):
            M[i, deg+1+k] = p*G**k
        v[i] = vals[end-n+i]
    sol = mp.lu_solve(M, v)
    return [sol[j] for j in range(deg+1)]


def fit_model(records, colour, part, deg, pdeg, window=0):
    """fit_values for one colour term and one power of L."""
    return fit_values(records, [L_parts(r, colour)[part] for r in records], deg, pdeg, window)


def beta_number(beta):
    b = Fraction(beta)
    return mp.mpf(b.numerator)/b.denominator


def weights(beta):
    """w = (1+beta)/2, F_S(w) and Ktilde(w) = (w^4+(1-w)^4)/(w(1-w))."""
    w = (1+beta_number(beta))/2
    v = w*(1-w)
    return w, 1/v+2*v, (w**4+(1-w)**4)/v


def rel(a, b):
    return abs(a-b)/abs(b) if b != 0 else abs(a)


def fmt(z, n=10):
    return mp.nstr(z, n)
