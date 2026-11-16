"""Numerical regression gates: disagreement, incomplete data and shell failures."""
import json
import copy
from pathlib import Path
import runpy
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import mpmath as mp
import secdec
from numeric_checks import close,finite_array,interior_values
from compare_stored import compare,validate_laurent

HERE=Path(__file__).resolve().parent


class NumericalTests(unittest.TestCase):
    def test_laurent_and_assembled_mismatches_fail(self):
        original=json.loads((HERE/'gqq_laurent_x1o2.json').read_text())
        def record(value):
            out=copy.deepcopy(original)
            next(iter(out['diagrams'].values()))['laurent']['ek']['0,0']=[value,0]
            return out
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'stored').mkdir()
            old=p/'stored/gqq_laurent_x1o2.json';new=p/'gqq_laurent_x1o2.json'
            old.write_text(json.dumps(record(1)));new.write_text(json.dumps(record(1e100)))
            run=lambda:subprocess.run([sys.executable,str(HERE/'compare_stored.py'),'gqq'],cwd=p,capture_output=True,text=True)
            self.assertNotEqual(run().returncode,0)
            new.write_bytes(old.read_bytes());self.assertEqual(run().returncode,0)
            (p/'stored/r0_theory_quarks.json').write_text('{"r0":1}')
            (p/'r0_theory_quarks.json').write_text('{"r0":2}')
            self.assertNotEqual(run().returncode,0)
            with self.assertRaises(ValueError):compare({'x':[1]},{'x':[]},1e-7,1e-6)
            for bad in ({},{'ek':[]},{'ek':{'0,0':['nan',0]}},{'ek':{'0,0':[float('nan'),0]}}):
                with self.assertRaises(ValueError):validate_laurent(bad,'test')

    def test_zero_scalar_structures_pass_but_missing_structures_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'stored').mkdir()
            old=p/'stored/n4_laurent_x3o10.json';new=p/'n4_laurent_x3o10.json'
            old.write_bytes((HERE/old.name).read_bytes());new.write_bytes(old.read_bytes())
            run=lambda:subprocess.run([sys.executable,str(HERE/'compare_stored.py'),'n4'],cwd=p,capture_output=True,text=True)
            passed=run();self.assertEqual(passed.returncode,0,passed.stdout+passed.stderr)
            incomplete=json.loads(old.read_text())
            incomplete['diagrams']['Vsbuba']['laurent'].pop('T4')
            # Identical omissions in old and fresh data must still be rejected.
            old.write_text(json.dumps(incomplete));new.write_bytes(old.read_bytes())
            self.assertNotEqual(run().returncode,0)

    def test_fake_laurent_values_fail_exact_examples(self):
        with patch.object(secdec,'laurent',return_value={(-1,-2):1e20,(0,-1):1e20,(0,0):1e20}),self.assertRaises(ValueError):
            runpy.run_path(str(HERE/'test_secdec.py'),run_name='__main__')

    def test_failed_python_stops_logging_targets(self):
        with tempfile.TemporaryDirectory() as tmp:
            for target in ('glauber','threeloop'):
                run=subprocess.run(['make','-f',str(HERE/'Makefile'),target,'PY=false'],cwd=tmp,capture_output=True,text=True)
                self.assertNotEqual(run.returncode,0)

    def test_nonfinite_and_disagreeing_quadratures_fail(self):
        for value in (float('nan'),float('inf')):
            with self.assertRaises(ValueError):close(value,1)
            with self.assertRaises(ValueError):finite_array([1,value])
        with self.assertRaises(ValueError):close(2,1)
        close(1+1e-7,1)

    def test_endpoint_cancellation_is_evaluated_without_discarding_nodes(self):
        import sympy as sp
        x,y=sp.symbols('x y')
        node=np.array([1-1e-14])
        value=interior_values(1/(1-x-y+x*y),[x,y],[node,node])
        expected=1/((1-node)*(1-node))
        np.testing.assert_allclose(value,expected,rtol=1e-14)

    def test_eikonal_offset_signs_against_quadrature(self):
        import glauber,glauber_regions
        with mp.workdps(35):
            A=mp.mpf('.7')
            for first,second in [('.2','.5'),('-.2','-.5'),('.2','-.5')]:
                s1,s2=mp.mpf(first),mp.mpf(second)
                coarse=glauber.J_numeric(A,s1,s2,mp.mpf('1e-6'))
                fine=glauber.J_numeric(A,s1,s2,mp.mpf('5e-7'))
                value=2*fine-coarse
                expected=glauber.J_closed(A,s1,s2)
                self.assertLess(abs(value-expected),mp.mpf('1e-8'))
                self.assertLess(abs(glauber_regions.exact(A,s1,s2)-expected),mp.mpf('1e-30'))


if __name__=='__main__':unittest.main()
