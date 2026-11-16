"""Laurent expansion of the one-loop collinear function J1 for g -> gg at a fixed x, graph by graph, in b space at
L_b = 0 (mu = mu_b, units p1^+ = mu = nu = 1), from the parametric integrands of run_gg.py; as assemble_gqq.py.

    python3 assemble_gg.py 3/10 [names]        (reads gg_integrands_x3o10.pkl, writes gg_laurent_x3o10.json)
    NQ=96,72,40 python3 assemble_gg.py 1/2     (higher quadrature orders: gg_laurent_x1o2_nq96-72-40.json)
    NEPS=16 python3 assemble_gg.py 1/2         (16 points on the eps circle: gg_laurent_x1o2_neps16.json)

For every graph and structure s in (S1, S2, S3, T4): j = a_s FT_s(eps) * int (weight) I_s, in units of the tree
coefficient of S1; FT_s = FTF(eps) for the linear structures S1-S3 and for T4 as well (the decomposition of T4
into its linear and its spin-3 parts is done in the combination, which multiplies the spin-3 part by
FT3/FTF); the weight z^(-1-eta) for the Wilson-line graphs.  Colour weights {CA, nf} from gg_diagrams.colour.
"""
import sys, json, pickle, time, os
import numpy as np
import sympy as sp
import mpmath as mp
from lc_engine import eps, D
from factors import eta, normal_form
from secdec import decompose, Sector, laurent
import gg_diagrams as G
from gg_driver import z, STRUCT
from cache_utils import (arguments, atomic_write, checkpoint, laurent_provenance,
                         load_integrands, new_laurent, quadrature, sector_cache, selected_names)

s_, t1, t2 = sp.symbols('s t1 t2')


def setup(xval, *, allow_legacy=False):
    global tag, INPUT_PATH, data, CACHE_FILE, CACHE_PROVENANCE, _CACHE
    tag = str(xval).replace('/', 'o')
    INPUT_PATH = f'gg_integrands_x{tag}.pkl'
    data = load_integrands(INPUT_PATH, 'gg', xval, allow_legacy=allow_legacy)
    CACHE_FILE = f'secdec_cache_gg_x{tag}.pkl'
    CACHE_PROVENANCE, _CACHE = sector_cache(CACHE_FILE)


def piece_sectors(expr, iv, prm, eik):
    key = (sp.srepr(expr), str(iv), str(prm), eik)
    if key in _CACHE:
        const, cube, raw = _CACHE[key]
        return [Sector(j, m, r, cube, const) for j, m, r in raw]
    za, zb = sp.Rational(iv[0]), sp.Rational(iv[1])
    dz = zb - za
    e = sp.sympify(expr).subs(D, 2 - 2*eps)
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
        assert za >= 0, ('Wilson-line graph with support at z < 0', iv)
        e = e*(za + dz*s_)**(-1 - eta)
    test = {v: sp.Rational(2, 7) for v in cube}
    const, facs = normal_form(e, cube, test)
    secs, raw = [], []
    for j, m, r, vm in decompose(facs, cube):
        secs.append(Sector(j, m, r, cube, const))
        raw.append((j, m, r))
    _CACHE[key] = (const, cube, raw)
    atomic_write(CACHE_FILE, {'provenance': CACHE_PROVENANCE, 'entries': _CACHE}, binary=True)
    return secs


def ftf(ep):
    return np.exp(-np.euler_gamma*ep)*complex(mp.gamma(1 - 2*ep)/(mp.gamma(1 + ep)*mp.gamma(1 - ep)))


def main():
    args = arguments(__doc__, assembly=True)
    setup(sp.Rational(args.x), allow_legacy=args.allow_legacy_input)
    nq, neps, suffix = quadrature()
    fn = f'gg_laurent_x{tag}{suffix}.json'
    metadata = laurent_provenance(data, INPUT_PATH, 'assemble_gg.py', nq, neps)
    out = new_laurent(fn, data, metadata, recompute=args.recompute)
    out['tree'] = {k: [complex(sp.N(v)).real, complex(sp.N(v)).imag] for k, v in data['tree'].items()}
    only = args.only
    selected_names(only, data['results'])
    for name, pieces in data['results'].items():
        if only is not None and name not in only:
            continue
        if only is None and name in out['diagrams']:
            continue
        diag = G.DIAGRAMS[name]
        col = G.colour(diag['colour'])
        t0 = time.time()
        res = {}
        for st in STRUCT:
            secs = []
            for iv, pcs in pieces:
                for prm, o in pcs:
                    if o[st] == 0:
                        continue
                    secs.extend(piece_sectors(o[st], iv, prm, diag['eik']))
            if not secs:
                res[st] = {}
                continue
            f = lambda ep, et: ftf(ep)*sum(S.value(ep, et, n=nq[len(S.vs)]) for S in secs)
            c = laurent(f, has_eta=diag['eik'], n_eps=neps)
            res[st] = {f'{i},{j}': [v.real, v.imag] for (i, j), v in c.items()}
        out['diagrams'][name] = {'colour': {k: [complex(v).real, complex(v).imag] for k, v in col.items()},
                                 'laurent': res}
        print(name, f'{time.time()-t0:.1f}s',
              {st: {k: (round(v[0], 8), round(v[1], 8)) for k, v in r.items() if abs(complex(*v)) > 1e-8}
               for st, r in res.items()}, flush=True)
        checkpoint(fn, out, 'diagrams')
    checkpoint(fn, out, 'diagrams')


if __name__ == '__main__':
    main()
