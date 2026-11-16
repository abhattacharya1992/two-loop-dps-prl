"""Parametric integrands of the scalar-loop graphs of J1 for g -> gg (n4_diagrams.py) at a fixed x, as run_gg.py.

    python3 run_n4.py 3/10 [names]       (writes n4_integrands_x3o10.pkl)
    python3 run_n4.py 3/10 --fdh         (the gluon-loop graphs of gg_diagrams.py with g^mu_mu = 4: fdh_integrands_x3o10.pkl)
"""
import sys, pickle, time
import sympy as sp
import vec_engine
import gg_diagrams
import n4_diagrams
from gg_driver import diagram_integrands, tree_coefficients, u, w
from cache_utils import arguments, checkpoint, new_integrands, selected_names

args = arguments(__doc__, fdh=True)
fdh = args.fdh
xval = sp.Rational(args.x)
only = args.only
tag = str(xval).replace('/', 'o')
if fdh:
    vec_engine.METRIC_TRACE = 4
    DIAGRAMS = {k: v for k, v in gg_diagrams.DIAGRAMS.items() if k.startswith('V')}   # the graphs with closed loops
    fn = f'fdh_integrands_x{tag}.pkl'
else:
    DIAGRAMS = n4_diagrams.DIAGRAMS
    fn = f'n4_integrands_x{tag}.pkl'
selected = selected_names(only, DIAGRAMS)
store = new_integrands(fn, 'fdh' if fdh else 'n4', xval, recompute=args.recompute)
store['failed_diagrams'] = [name for name in store['failed_diagrams'] if name not in selected]
tree = tree_coefficients(gg_diagrams.tree(), xval)
store['tree'] = tree
norm = tree['S1']
b = sorted({sp.Rational(-1), -xval, xval - 1, sp.Integer(0), xval, 1 - xval, sp.Integer(1)})
intervals = [(sp.Rational(-3, 2), -1)] + list(zip(b[:-1], b[1:])) + [(1, sp.Rational(3, 2))]
for name, diag in DIAGRAMS.items():
    if only is not None and name not in only:
        continue
    if only is None and name in store['results']:
        continue
    store['results'].pop(name, None)
    t0 = time.time()
    try:
        res = diagram_integrands(diag, xval, intervals, norm, order=(u, w))
    except ValueError as e1:
        try:
            res = diagram_integrands(diag, xval, intervals, norm, order=(w, u))
        except ValueError as e2:
            print(name, 'FAILED', e1, '|', e2, flush=True)
            store['failed_diagrams'].append(name)
            checkpoint(fn, store, 'results')
            continue
    nz = [(iv, pcs) for iv, pcs in res if any(v != 0 for prm, o in pcs for v in o.values())]
    store['results'][name] = nz
    checkpoint(fn, store, 'results')
    print(name, f'{time.time()-t0:.1f}s', 'support:', [iv for iv, _ in nz], 'params:',
          [[len(p) for p, _ in pcs] for _, pcs in nz], flush=True)
print('tree (S1 normalised to its own value):', {k: sp.nsimplify(v/norm) for k, v in tree.items()})
checkpoint(fn, store, 'results')
if store['failed_diagrams']:
    sys.exit('integrand generation failed for: ' + ', '.join(store['failed_diagrams']))
