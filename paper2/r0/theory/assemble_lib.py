"""The assembly of assemble_gg.py as functions, for the scalar graphs (assemble_n4.py) and the scheme check
(check_fdh.py): parametric integrands (run_gg.diagram_integrands) -> sector decomposition -> Laurent coefficients in
eta and eps, graph by graph, b space at L_b = 0.  Same quadrature orders and the same 32-point eps circle as
assemble_gg.py; the sector cache is per file name."""
import json, pickle, time, os
import numpy as np
import sympy as sp
import mpmath as mp
from lc_engine import eps, D
from factors import eta, normal_form
from secdec import decompose, Sector, laurent
from gg_driver import z, STRUCT
from cache_utils import (atomic_write, checkpoint, laurent_provenance, new_laurent,
                         sector_cache, selected_names)

s_, t1, t2 = sp.symbols('s t1 t2')


class Assembler:
    def __init__(self, cache_file):
        self.cache_file = cache_file
        self.provenance, self.cache = sector_cache(cache_file)

    def piece_sectors(self, expr, iv, prm, eik):
        key = (sp.srepr(expr), str(iv), str(prm), eik)
        if key in self.cache:
            const, cube, raw = self.cache[key]
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
        self.cache[key] = (const, cube, raw)
        atomic_write(self.cache_file, {'provenance': self.provenance, 'entries': self.cache}, binary=True)
        return secs


def ftf(ep):
    return np.exp(-np.euler_gamma*ep)*complex(mp.gamma(1 - 2*ep)/(mp.gamma(1 + ep)*mp.gamma(1 - ep)))


def laurent_of_graph(asm, pieces, eik, nq=None, neps=32):
    """{structure: {'i,j': [re, im]}} for one graph from its integrand pieces [(interval, [(params, {s: I})])]"""
    nq = nq or {1: 64, 2: 48, 3: 28}
    res = {}
    for st in STRUCT:
        secs = []
        for iv, pcs in pieces:
            for prm, o in pcs:
                if o[st] == 0:
                    continue
                secs.extend(asm.piece_sectors(o[st], iv, prm, eik))
        if not secs:
            res[st] = {}
            continue
        f = lambda ep, et: ftf(ep)*sum(S.value(ep, et, n=nq[len(S.vs)]) for S in secs)
        c = laurent(f, has_eta=eik, n_eps=neps)
        res[st] = {f'{i},{j}': [v.real, v.imag] for (i, j), v in c.items()}
    return res


def assemble(data, diagrams, colour_fn, out_fn, cache_file, only=None, nq=None, neps=32, verbose=True,
             *, input_path, script='assemble_n4.py', recompute=False):
    """data = pickle of run_*.py; diagrams = registry dict; writes/extends out_fn (json)"""
    nq = nq or {1: 64, 2: 48, 3: 28}
    selected_names(only, data['results'])
    metadata = laurent_provenance(data, input_path, script, nq, neps)
    out = new_laurent(out_fn, data, metadata, recompute=recompute)
    asm = Assembler(cache_file)
    out['tree'] = {k: [complex(sp.N(v)).real, complex(sp.N(v)).imag] for k, v in data['tree'].items()}
    for name, pieces in data['results'].items():
        if only is not None and name not in only:
            continue
        if only is None and name in out['diagrams']:
            continue
        diag = diagrams[name]
        col = colour_fn(diag['colour'])
        t0 = time.time()
        res = laurent_of_graph(asm, pieces, diag['eik'], nq=nq, neps=neps)
        out['diagrams'][name] = {'colour': {k: [complex(v).real, complex(v).imag] for k, v in col.items()},
                                 'laurent': res}
        if verbose:
            print(name, f'{time.time()-t0:.1f}s',
                  {st: {k: (round(v[0], 8), round(v[1], 8)) for k, v in r.items() if abs(complex(*v)) > 1e-8}
                   for st, r in res.items()}, flush=True)
        checkpoint(out_fn, out, 'diagrams')
    checkpoint(out_fn, out, 'diagrams')
    return out
