#!/usr/bin/env python3
"""Recomputes results_stage_b.json (the numbers of appendix E of Paper 2) and checks that the stored file is
reproduced and that the statements of the appendix hold:
  - the ultrasoft cancellation: both colour components within 1e-4 at lambda = 1e-4, falling with lambda,
    while one attachment alone is of order one; no cancellation for Glauber momenta;
  - the explicit two-loop graph: lambda^6 F converges at l_perp ~ lambda^2; the coefficient of the single
    logarithm at l_perp ~ lambda is the eikonal one within 1e-2; no logarithm in the real part;
  - the bookkeeping: r2 = c_u/2, r1 = c_u kappa_u + d_u + g1; the rapidity coefficient has no eps^0 term.
    python3 test_kernel.py          (about 10 s)
"""
import json, math
import stage_b_results as SB

def close(a, b, rel=1e-8, path=''):
    if isinstance(a, dict):
        assert set(a) == set(b), path
        for k in a:
            close(a[k], b[k], rel, path + '/' + k)
    elif isinstance(a, (list, tuple)):
        assert len(a) == len(b), path
        for i, (u, v) in enumerate(zip(a, b)):
            close(u, v, rel, f'{path}[{i}]')
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=rel, abs_tol=1e-300), f'{path}: {a} != {b}'
    else:
        assert a == b, f'{path}: {a!r} != {b!r}'

stored = json.load(open('results_stage_b.json'))
new = json.loads(json.dumps({'telescoping': SB.telescoping(), 'model': SB.model(), 'logs': SB.logs()}))
close(new, stored)
t, m, l = new['telescoping'], new['model'], new['logs']
us = t['ultrasoft']
assert max(us[-1]['dev_Tk'], us[-1]['dev_Tkbar']) < 1e-4
assert all(us[i+1]['dev_Tk'] < us[i]['dev_Tk'] for i in range(len(us) - 1))
assert all(r['single_attachment'] > 0.1 for r in us)
assert min(min(g['dev_Tk'], g['dev_Tkbar']) for g in t['glauber']) > 0.05
assert t['residues_vs_quadrature_worst'] < 1e-5
for c, r in m['ultrasoft_lambda6F'].items():
    ch = r['rel_changes']
    assert all(ch[i+1] < 0.6*ch[i] for i in range(len(ch) - 1)) and ch[-1] < 5e-2
assert max(s['rel_diff'] for s in m['log_slopes']) < 1e-2
assert max(s['real_part_change'] for s in m['log_slopes']) < 1e-3
assert l['r2'] == 'c_u/2' and l['r1'] == 'c_u*kappa_u + d_u + g1' and l['rapidity_coefficient_eps0'] == '0'
print('results_stage_b.json reproduced; all statements of appendix E hold')
