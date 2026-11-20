#!/usr/bin/env python3
"""Reruns the comparison with planar N = 4 (stage_d.py) into a temporary file and checks that it
reproduces results_stage_d.json, the file the manuscript reads, and the statements of Sec. 5.1.

    python3 test_n4.py        (or `make test`)
"""
import json, pathlib, subprocess, sys, tempfile

HERE = pathlib.Path(__file__).resolve().parent
stored = json.loads((HERE/'results_stage_d.json').read_text())
with tempfile.TemporaryDirectory() as tmp:
    out = pathlib.Path(tmp)/'results.json'
    log = subprocess.run([sys.executable, str(HERE/'stage_d.py'), str(HERE/'de_data.json'), str(out)],
                         capture_output=True, text=True, check=True).stdout
    new = json.loads(out.read_text())
assert new == stored, 'stage_d.py no longer reproduces results_stage_d.json'
mhv, nmhv = stored['three_to_three'], stored['nmhv']
assert stored['loops_compared'] == 5 and len(mhv) == 25
assert stored['n_pass'] == 18 and stored['n_fixed'] == 6 and stored['n_external'] == 2 and stored['n_norm'] == 1
assert not any(r['status'].startswith('FAIL') for r in mhv + nmhv)
assert stored['four_loop_cusp_agrees'] and stored['LL_Gaussian_integral'] and stored['nmhv_LL_exponential']
assert stored['components_related'] and stored['two_to_four_single_log_difference'] == '0'
assert all(r['DE'] == r['prediction'] for r in stored['two_to_four_two_loops'])   # 2 -> 4: every coefficient, no v
assert len(nmhv) == 16 and stored['nmhv_n_checked'] == stored['nmhv_n_pass'] == 15
assert [(r['L'], r['k']) for r in mhv if r['status'] == 'normalisation'] == [(1, 1)]
assert [(r['L'], r['k']) for r in nmhv if r['status'] == 'normalisation'] == [(1, 0)]
assert stored['solved_constants']['r0'] == 'zeta_2' and stored['solved_constants']['B2'] == '3*zeta_3'
assert stored['solved_constants']['B3'] == '-2*zeta_2*zeta_3 - 10*zeta_5'
assert stored['external_inputs']['values'] == {'B2': '3*zeta_3', 'B3': '-2*zeta_2*zeta_3 - 10*zeta_5'}
assert {(r['L'],r['k']) for r in mhv if r['status'] == 'PASS' and (r['L'],r['k']) in ((3,2),(4,2))} == {(3,2),(4,2)}
assert len(stored['nmhv_five_loop_prediction']) == 9
assert 'wrote' in log
print('n4: stage_d.py reproduces results_stage_d.json; MHV 18 of 18 tests agree, 6 constants fitted and B2/B3 externally fixed; '
      'NMHV 15 of 15, no new constant; four-loop cusp agrees; the two components are related exactly')
