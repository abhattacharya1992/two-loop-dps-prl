"""Scheme check of the scalar graphs (n4_diagrams.py): the gluon-loop graphs of J1 in the four-dimensional-helicity
scheme (metric traces g^mu_mu = 4, loop momenta in d dimensions; fdh_laurent_x*.json from run_n4.py/assemble_n4.py
--fdh) minus the same graphs in conventional dimensional regularisation (gg_laurent_x*.json) must equal eps times the
graphs of one complex adjoint scalar (n4_laurent_x*.json): the 2 eps epsilon-scalars of dimensional reduction are eps
complex ones.  Checked colour weight by colour weight for the Laurent coefficients in eps (the 1/eps^2 and 1/eps of
FDH - CDR against nothing and the 1/eps^2 of the scalar graphs, the eps^0 of FDH - CDR against the 1/eps of the scalar
graphs), for every structure.  The fermion loops and the ghosts have no metric trace and must not change.

    python3 check_fdh.py 1/2        (writes fdh_check_x1o2.json)
"""
import sys, json, math
from fractions import Fraction
import numpy as np
from cache_utils import atomic_write, diagram_names, provenance, sources
from laurent_inputs import load_laurent

STRUCT = ('S1', 'S2', 'S3', 'T4')
xval = str(Fraction(sys.argv[1] if len(sys.argv) > 1 else '1/2'))
tag = xval.replace('/', 'o')
paths = {kind: f'{kind}_laurent_x{tag}.json' for kind in ('gg', 'fdh', 'n4')}
records = {kind: load_laurent(path, kind, expected_x=xval) for kind, path in paths.items()}
cdr, fdh, ns = (records[kind]['diagrams'] for kind in ('gg', 'fdh', 'n4'))
for kind, record in records.items():
    if str(record.get('x')) != xval:
        sys.exit(f'{paths[kind]}: recorded x does not match requested x={xval}')
    if set(record['diagrams']) != set(diagram_names(kind)):
        sys.exit(f'{paths[kind]}: incomplete diagram set; scheme check cannot certify an incomplete calculation')
    if 'complete' in record and not record['complete']:
        sys.exit(f'{paths[kind]}: incomplete calculation')


def weighted(dd, col):
    """{structure: {eps order j: value}} of the eta^0 terms, colour-weighted"""
    out = {s: {} for s in STRUCT}
    for name, d in dd.items():
        w = complex(*d['colour'].get(col, [0, 0]))
        if w == 0:
            continue
        for s in STRUCT:
            for key, (re, im) in d['laurent'][s].items():
                i, j = map(int, key.split(','))
                if i != 0:
                    continue
                out[s][j] = out[s].get(j, 0) + w*complex(re, im)
    return out


worst = 0
finite = True
for col in ('CA', 'nf'):
    a = weighted({k: v for k, v in cdr.items() if k in fdh}, col)
    b = weighted(fdh, col)
    sc = weighted(ns, 'ns') if col == 'CA' else {s: {} for s in STRUCT}
    for s in STRUCT:
        diff = {j: b[s].get(j, 0) - a[s].get(j, 0) for j in (-2, -1, 0)}
        # eps * scalar: coefficient of eps^j is scalar's eps^(j-1)
        pred = {j: sc[s].get(j - 1, 0) for j in (-2, -1, 0)}
        for j in (-2, -1, 0):
            dev = abs(diff[j] - pred[j])
            finite = finite and math.isfinite(dev)
            worst = max(worst, dev)
            print(f'{col} {s} eps^{j:+d}: FDH - CDR = {diff[j].real:+.9f}  eps x scalar = {pred[j].real:+.9f}  dev {dev:.1e}')
print('largest deviation', f'{worst:.1e}')
passed = finite and worst < 1e-9
print('FDH = CDR + eps x (one complex adjoint scalar):', 'OK' if passed else 'FAILED')
metadata = provenance.fingerprint('fdh-scheme-check', inputs={path: path for path in paths.values()},
    sources=sources(('check_fdh.py', 'laurent_inputs.py', 'gg_diagrams.py', 'n4_diagrams.py',
                     'gg_driver.py', 'vec_engine.py', 'lc_engine.py')),
    parameters={'x': tag.replace('o', '/'), 'absolute_tolerance': 1e-9,
                'input_provenance_status': {kind: ('tagged' if 'provenance' in record else 'legacy-unverified')
                                            for kind, record in records.items()}},
    libraries=('numpy', 'sympy', 'mpmath'))
atomic_write(f'fdh_check_x{tag}.json',
    {'x': tag.replace('o', '/'), 'largest deviation': worst if finite else None,
     'graphs': sorted(fdh), 'passed': passed, 'absolute_tolerance': 1e-9, 'provenance': metadata,
     'statement': 'gluon-loop graphs of J1 in the FDH scheme minus CDR = eps x the graphs of one complex adjoint scalar'})
if not passed:
    sys.exit(1)
