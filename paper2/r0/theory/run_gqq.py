import sys, pickle, time
import sympy as sp
import mpmath as mp
import lc_engine as E
from gqq_diagrams import DIAGRAMS, x, z, u, w
from driver import diagram_integrands
from cache_utils import arguments, checkpoint, new_integrands, selected_names

args = arguments('Generate g -> q qbar integrands with verified resumable checkpoints.')
xval = sp.Rational(args.x)
only = args.only
selected = selected_names(only, DIAGRAMS)
fn = f'gqq_integrands_x{str(xval).replace("/", "o")}.pkl'
store = new_integrands(fn, 'gqq', xval, recompute=args.recompute)
store['failed_diagrams'] = [name for name in store['failed_diagrams'] if name not in selected]
b = sorted({sp.Rational(-1), -xval, xval - 1, sp.Integer(0), xval, 1 - xval, sp.Integer(1)})
intervals = [(sp.Rational(-3, 2), -1)] + list(zip(b[:-1], b[1:])) + [(1, sp.Rational(3, 2))]
results = store['results']
for name, diag in DIAGRAMS.items():
    if name not in selected or (only is None and name in results):
        continue
    results.pop(name, None)
    t0 = time.time()
    try:
        res = diagram_integrands(diag, xval, intervals, order=(u, w))
    except ValueError as e1:
        try:
            res = diagram_integrands(diag, xval, intervals, order=(w, u))
        except ValueError as e2:
            print(name, 'FAILED', e1, '|', e2)
            store['failed_diagrams'].append(name)
            checkpoint(fn, store, 'results')
            continue
    nz = [(iv, pcs) for iv, pcs in res if any(v != 0 for prm, o in pcs for v in o.values())]
    results[name] = nz
    checkpoint(fn, store, 'results')
    print(name, f'{time.time()-t0:.1f}s', 'support:', [iv for iv, _ in nz], 'params:', [[len(p) for p, _ in pcs] for _, pcs in nz], flush=True)
checkpoint(fn, store, 'results')
if store['failed_diagrams']:
    sys.exit('integrand generation failed for: ' + ', '.join(store['failed_diagrams']))
