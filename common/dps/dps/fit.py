"""Interpolation fits on the forward line.

At fixed velocity a forward amplitude is a polynomial in G plus power corrections of relative
order p_T^2, themselves polynomials in G:

    sum_{j<=deg} c_j G^j + p_T^2 sum_{l<=pdeg} d_l G^l .

`leading(records, vals, deg, pdeg, window)` fixes the deg + pdeg + 2 coefficients by exact
interpolation through that many consecutive points and returns c_0 ... c_deg.  Window 0 ends at
the smallest delta of the scan, window 1 is shifted one point up in delta.  Interpolation leaves no
residual, so the change between the two windows is quoted as the uncertainty of a coefficient.

The matrix is built and solved at the current mpmath precision (the analysis scripts of Paper 2
use 50 digits), in exactly the arithmetic of hh/hh_two_loop/scan_tools.py and
zz/scripts/analyse_zz.py, so that the stored results are reproduced to the last digit.
"""
import mpmath as mp


def leading(records, vals, deg, pdeg, window=0):
    """c_0 ... c_deg from the records (sorted by decreasing delta, each with 'G' and 'pT2' as
    mpmath numbers) and the values `vals`; None if the window does not fit in the scan."""
    n = deg + pdeg + 2
    if len(records) != len(vals):
        raise ValueError('Different numbers of fit records and values')
    if any(not mp.isfinite(v) for v in vals):
        raise ValueError('Nonfinite amplitude in fit')
    if any(not mp.isfinite(r[k]) for r in records for k in ('G', 'pT2')):
        raise ValueError('Nonfinite fit kinematics')
    end = len(records) - window
    if end - n < 0:
        return None
    M = mp.matrix(n, n)
    v = mp.matrix(n, 1)
    for i in range(n):
        r = records[end - n + i]
        G, p = r['G'], r['pT2']
        for j in range(deg + 1):
            M[i, j] = G ** j
        for l in range(pdeg + 1):
            M[i, deg + 1 + l] = p * G ** l
        v[i] = vals[end - n + i]
    sol = mp.lu_solve(M, v)
    if any(not mp.isfinite(c) for c in sol):
        raise ValueError('Nonfinite fitted coefficient')
    return [sol[j] for j in range(deg + 1)]


def plain(Gs, vals, deg):
    """Coefficients c_0 ... c_deg of the leading-power polynomial alone through deg + 1 points."""
    n = deg + 1
    if len(Gs) != n or len(vals) != n:
        raise ValueError('Wrong number of polynomial fit points')
    if any(not mp.isfinite(v) for v in [*Gs, *vals]):
        raise ValueError('Nonfinite polynomial fit input')
    M = mp.matrix(n, n)
    for i, g in enumerate(Gs):
        for j in range(n):
            M[i, j] = g ** j
    sol = mp.lu_solve(M, mp.matrix(vals))
    if any(not mp.isfinite(c) for c in sol):
        raise ValueError('Nonfinite fitted coefficient')
    return [sol[j] for j in range(n)]


def with_shift(records, vals, deg, pdeg):
    """(coefficients of window 0, |window 1 - window 0| per coefficient); the shifts are None
    when the scan has no second window."""
    c0 = leading(records, vals, deg, pdeg, 0)
    c1 = leading(records, vals, deg, pdeg, 1)
    if c0 is None:
        return None, None
    shifts = [abs(a - b) for a, b in zip(c1, c0)] if c1 is not None else None
    return c0, shifts
