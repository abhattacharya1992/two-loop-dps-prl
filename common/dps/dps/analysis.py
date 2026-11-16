"""The two analyses of Paper 2 for any process: the leading logarithms against the
renormalisation-group prediction (Secs. 3 and 6), and the constant r_0 of the boundary term
(Sec. 7).

All fits are interpolations through consecutive points of the scan (fit.leading); the quoted
uncertainty is the change when the window moves one point up in delta.  Gluon splittings
(gg -> HH, gg -> AA): the two-loop remainder is split into the colour monomials N^2, N n_f, n_f^2,
n_f/N, and the N^2 and N n_f terms carry the DPS logarithms; the one-loop terms C_1 = 2 pi N F sigma
and K = pi N W are analytic (for gg -> AA they are also fitted from the one-loop files).  Quark
splittings (gg -> ZZ): C_F and C_A parts of the two-loop helicity amplitudes, one-loop C_1 and K
fitted, VVamp units (alpha_s/(2 pi) per loop, so two-loop/one-loop ratios are doubled).
"""
from fractions import Fraction

import mpmath as mp

from . import fit, theory
from .kinematics import momentum_fraction
from .points import series

mp.mp.dps = 50

# fit models (degree in G, degree of the p_T^2 correction) ---------------------------------------
GLUON_MODELS = {1: (3, 2), 2: (2, 2)}        # by component: M_1 cubic, M_2 quadratic, + p_T^2 (1, G, G^2)
QUARK_MODEL = (3, 3)                         # two-loop gg -> ZZ: cubic + p_T^2 cubic
QUARK_ONE_LOOP = (1, 2)
HEL = {'+': 1, '-': -1, '0': 0}


def jz_change(h):
    l1, l2, l3, l4 = (HEL[c] for c in h)
    return (l1 - l2) - (l3 - l4)


def rel(a, b):
    return abs(a - b) / abs(b) if b != 0 else abs(a)


# gluon splittings ---------------------------------------------------------------------------
def leading_logs_gluon(process, regular, beta, L=0, colours=theory.GLUON_COLOURS):
    """Fitted coefficients of G^3, G^2 (component 1) and G^2, G^1 (component 2) at mu = m for each
    colour monomial, with the window shift and the prediction of eq. (5.5); the key 'L' picks the
    power of ln(m^2/mu^2) of the gg -> HH points (only L = 0 exists for gg -> AA)."""
    x = momentum_fraction(beta)
    pred = theory.predictions_gluon(x, 1, 1, process.sigma, process.F0)
    rows = []
    for comp, name in ((1, 'M1'), (2, 'M2')):
        recs = regular.get((comp, beta))
        if not recs:
            continue
        deg, pdeg = GLUON_MODELS[comp]
        for power in ((3, 2) if comp == 1 else (2, 1)):
            key = f'G{power}'
            for col in colours:
                vals = series(recs, f'{col}:L{L}' if process.name == 'hh' else col)
                c0, shifts = fit.with_shift(recs, vals, deg, pdeg)
                p = pred[name][key].get(col, mp.mpf(0))
                rows.append({'component': comp, 'power': power, 'colour': col, 'value': c0[power],
                             'shift': shifts[power] if shifts else None, 'prediction': p,
                             'deviation': rel(c0[power], p) if p != 0 else abs(c0[power])})
    return rows


def r0_gluon(process, regular, beta):
    """r_0 for the two components from the N^2 and N n_f terms (eq. 7.1 with the known parts of
    theory.known_parts_gluon), and the n_f/N and n_f^2 terms of the same coefficients, which must
    vanish; values in units of the one-loop terms."""
    x = momentum_fraction(beta)
    known = theory.known_parts_gluon(beta, process)
    out = {'x': x}
    for comp, nu, power in ((1, 'nu0', 1), (2, 'nu2', 0)):
        recs = regular.get((comp, beta))
        if not recs:
            continue
        deg, pdeg = GLUON_MODELS[comp]
        norm = process.c1_over_N(x) if comp == 1 else process.K_over_N(x)
        res = {}
        for col in theory.GLUON_COLOURS:
            vals = series(recs, f'{col}:L0' if process.name == 'hh' else col)
            c0, shifts = fit.with_shift(recs, vals, deg, pdeg)
            res[col] = {'over_oneloop': c0[power] / norm, 'shift': shifts[power] / abs(norm) if shifts else None}
        res['r0_N'] = res['N2']['over_oneloop'] - known[nu]['N']
        res['r0_nf'] = res['N_Nf']['over_oneloop'] - known[nu]['nf']
        res['known_N'], res['known_nf'] = known[nu]['N'], known[nu]['nf']
        out[nu] = res
    return out


def r0_closed_gluon(process, beta):
    """The closed form (7.2) for the N and n_f parts of both components at this velocity."""
    x = momentum_fraction(beta)
    return {nu_name: {'N': theory.r0_closed(process, nu, x, 1, 0), 'nf': theory.r0_closed(process, nu, x, 0, 1)}
            for nu, nu_name in ((0, 'nu0'), (2, 'nu2'))}


# quark splittings (gg -> ZZ) -----------------------------------------------------------------
def one_loop_quark(recs):
    """(c1, K) of the one-loop helicity amplitude: coefficient of G and constant, window 0 and 1."""
    vals = series(recs, '1l')
    f0 = fit.leading(recs, vals, *QUARK_ONE_LOOP, 0)
    f1 = fit.leading(recs, vals, *QUARK_ONE_LOOP, 1)
    return f0, f1


def leading_logs_quark(regular, beta, pdeg=QUARK_MODEL[1]):
    """Two-loop coefficients of the gg -> ZZ helicity amplitudes against the quark-splitting
    prediction (analyse_zz.py): for dJz = 0 amplitudes with a one-loop logarithm, f3/c1 and
    f2/c1; for |dJz| = 2, f2/K and f1/K; C_F and C_A parts separately."""
    x = momentum_fraction(beta)
    Lx = mp.log(x * (1 - x))
    pred = {'2l_CF': {'dJz0': {3: -mp.mpf(4) / 3, 2: -2 * Lx - 3}, 'dJz2': {2: -4, 1: -4 * Lx - 2}},
            '2l_CA': {'dJz0': {3: mp.mpf(2) / 3, 2: 2 * Lx}, 'dJz2': {2: 2, 1: 4 * Lx - 2}}}
    hels = sorted(h for (h, b) in regular if b == beta)
    cmax = max(abs(one_loop_quark(regular[(h, beta)])[0][1]) for h in hels)
    rows = []
    for h in hels:
        recs = regular[(h, beta)]
        one, one_w = one_loop_quark(recs)
        c1, K = one[1], one[0]
        dj = jz_change(h)
        if dj == 0 and abs(c1) > 1e-3 * cmax:
            kind, norm, powers = 'dJz0', c1, (3, 2)
        elif abs(dj) == 2 and abs(K) > 1e-6:
            kind, norm, powers = 'dJz2', K, (2, 1)
        else:
            continue
        for part in ('2l_CF', '2l_CA'):
            vals = series(recs, part)
            f0, shifts = fit.with_shift(recs, vals, 3, pdeg)
            if f0 is None:
                continue
            for power in powers:
                ratio = f0[power] / norm
                rows.append({'helicity': h, 'kind': kind, 'part': part, 'power': power, 'c1': c1, 'K': K,
                             'ratio': ratio, 'shift': shifts[power] / abs(norm) if shifts else None,
                             'prediction': pred[part][kind][power],
                             'deviation': rel(ratio, pred[part][kind][power])})
    return rows


def r0_quark(regular, beta):
    """r_0 for quark splittings from the gg -> ZZ fits (r0/data/r0_from_zz.py), in a_s units."""
    x = momentum_fraction(beta)
    hels = sorted(h for (h, b) in regular if b == beta)
    cmax = max(abs(one_loop_quark(regular[(h, beta)])[0][1]) for h in hels)
    out = {'x': x, 'nu0': {}, 'nu2': {}}
    for h in hels:
        recs = regular[(h, beta)]
        one, _ = one_loop_quark(recs)
        c1, K = one[1], one[0]
        dj = jz_change(h)
        if dj == 0 and abs(c1) > 1e-3 * cmax:
            comp, norm, j = 'nu0', c1, 1
        elif abs(dj) == 2 and abs(K) > 1e-6:
            comp, norm, j = 'nu2', K, 0
        else:
            continue
        row = {}
        for part, (CF, CA) in (('2l_CF', (1, 0)), ('2l_CA', (0, 1))):
            f0, shifts = fit.with_shift(recs, series(recs, part), *QUARK_MODEL)
            kn = theory.known_parts_quark(x, CF, CA)[0 if comp == 'nu0' else 1]
            row[part] = {'r0': 2 * f0[j] / norm - kn, 'shift': 2 * shifts[j] / abs(norm) if shifts else None}
        out[comp][h] = row
    return out


def betas_of(regular):
    return sorted({b for (_, b) in regular}, key=Fraction)
