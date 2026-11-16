"""Laurent expansion of the one-loop collinear function J1 (g -> q qbar) at a fixed x, diagram by diagram, in
b-space at L_b = 0 (mu = mu_b, units p1^+ = mu = nu = 1), from the parametric integrands of run_gqq.py.

    python3 assemble_gqq.py 3/10        (reads gqq_integrands_x3o10.pkl)

For every diagram d and helicity structure st:  j_d = a_s FTF(eps) * pref_d * sum_pieces int (weight) I,
FTF(eps) = e^{-gamma eps} Gamma(1-2eps)/(Gamma(1+eps)Gamma(1-eps)) (Fourier factor of (k2)^(-eps)), the
weight z^(-1-eta) for the Wilson-line diagrams.  Output: coefficients of eta^i eps^j (i = -1, 0; j = -2..0)
for the C_F and C_A parts, written to gqq_laurent_x...json.
"""
import sys, json, pickle, time
import numpy as np
import sympy as sp
import mpmath as mp
from lc_engine import eps, D
from factors import eta, normal_form
from secdec import decompose, Sector, laurent
from gqq_diagrams import DIAGRAMS, z
from cache_utils import (arguments, atomic_write, checkpoint, laurent_provenance,
                         load_integrands, new_laurent, quadrature, sector_cache, selected_names)

s_, t1, t2 = sp.symbols('s t1 t2')


def setup(xval, *, allow_legacy=False):
    global tag, INPUT_PATH, data, CACHE_FILE, CACHE_PROVENANCE, _CACHE
    tag = str(xval).replace('/', 'o')
    INPUT_PATH = f'gqq_integrands_x{tag}.pkl'
    data = load_integrands(INPUT_PATH, 'gqq', xval, allow_legacy=allow_legacy)
    CACHE_FILE = f'secdec_cache_x{tag}.pkl'
    CACHE_PROVENANCE, _CACHE = sector_cache(CACHE_FILE)


def piece_sectors(expr, iv, prm, eik):
    key = (sp.srepr(expr), str(iv), str(prm), eik)
    if key in _CACHE:
        const, cube, raw = _CACHE[key]
        return [Sector(j, m, r, cube, const) for j, m, r in raw]
    return _piece_sectors(expr, iv, prm, eik)


def _piece_sectors(expr, iv, prm, eik):
    za, zb = sp.Rational(iv[0]), sp.Rational(iv[1])
    dz = zb - za
    e = expr.subs(D, 2 - 2*eps)
    jac = dz
    if len(prm) == 2:
        a1, a2 = prm
        e = e.subs({a1: t1, a2: (1 - t1)*t2}, simultaneous=True)
        jac *= (1 - t1)
        cube = [s_, t1, t2]
    elif len(prm) == 1:
        e = e.subs(prm[0], t1)
        cube = [s_, t1]
    else:
        cube = [s_]
    e = e.subs(z, za + dz*s_)*jac
    if eik:
        assert za >= 0
        e = e*(za + dz*s_)**(-1 - eta)
    test = {v: sp.Rational(2, 7) for v in cube}
    const, facs = normal_form(e, cube, test)
    secs = []
    raw = []
    for j, m, r, vm in decompose(facs, cube):
        secs.append(Sector(j, m, r, cube, const))
        raw.append((j, m, r))
    _CACHE[(sp.srepr(expr), str(iv), str(prm), eik)] = (const, cube, raw)
    atomic_write(CACHE_FILE, {'provenance': CACHE_PROVENANCE, 'entries': _CACHE}, binary=True)
    return secs


def ftf(ep):
    return np.exp(-np.euler_gamma*ep)*complex(mp.gamma(1 - 2*ep)/(mp.gamma(1 + ep)*mp.gamma(1 - ep)))


def main():
    args = arguments(__doc__, assembly=True)
    setup(sp.Rational(args.x), allow_legacy=args.allow_legacy_input)
    # points on the eps circle |eps| = 0.2: FTF(eps) has a pole at eps = 1/2, so the aliasing error of the
    # Laurent coefficients is ~ (2 x 0.2)^NEPS: 4e-7 for 16 points (seen), 2e-13 for 32 (default)
    nq, neps, suffix = quadrature()
    fn = f'gqq_laurent_x{tag}{suffix}.json'
    metadata = laurent_provenance(data, INPUT_PATH, 'assemble_gqq.py', nq, neps)
    out = new_laurent(fn, data, metadata, recompute=args.recompute)
    only = args.only
    selected_names(only, data['results'])
    for name, pieces in data['results'].items():
        if only is not None and name not in only:
            continue
        if only is None and name in out['diagrams']:
            continue
        diag = DIAGRAMS[name]
        t0 = time.time()
        res = {}
        for st in ('ek', 'ke'):
            secs = []
            for iv, pcs in pieces:
                for prm, o in pcs:
                    if o[st] == 0:
                        continue
                    secs.extend(piece_sectors(o[st], iv, prm, diag['eik']))
            f = lambda ep, et: ftf(ep)*sum(S.value(ep, et, n=nq[len(S.vs)]) for S in secs)
            c = laurent(f, has_eta=diag['eik'], n_eps=neps)
            res[st] = {f'{i},{j}': [v.real, v.imag] for (i, j), v in c.items()}
        out['diagrams'][name] = {'pref': {k: str(v) for k, v in diag['pref'].items()}, 'laurent': res}
        print(name, f'{time.time()-t0:.1f}s', {st: {k: round(v[0], 10) for k, v in r.items() if abs(v[0]) > 1e-9}
                                               for st, r in res.items()}, flush=True)
        checkpoint(fn, out, 'diagrams')
    checkpoint(fn, out, 'diagrams')


if __name__ == '__main__':
    main()
