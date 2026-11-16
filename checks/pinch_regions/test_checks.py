"""Positive, negative and replay controls for the curated scalar checks."""
from __future__ import annotations
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import sympy as s
import param_space as p
import exact_checks as e
from graphs import Graph, all_graphs

class Algebra(unittest.TestCase):
    def test_constructive_sign_witnesses(self):
        f={(1,1,0,0):p.Q(2),(0,0,1,1):p.Q(-17)}
        points=p.witnesses(f)
        self.assertGreater(p.ev(f,points[1]),0)
        self.assertLess(p.ev(f,points[-1]),0)
        root,_=p.root_point(f)
        self.assertEqual(p.ev(f,root),0)
        self.assertTrue(all(q>0 for q in root))

    def test_definite_sign_and_rejected_domains(self):
        self.assertEqual(set(p.witnesses({(1,0):1,(0,1):2})),{1})
        for invalid in ({},{(2,0):1,(1,1):-2,(0,2):1},{(1,0):1,(1,1):-1}):
            with self.assertRaises(ArithmeticError):p.witnesses(invalid)

    def test_normal_coordinate_recentering(self):
        # F=z1*z2+lambda*(z1+z2)+lambda^2+lambda^3: critical order 3, not 2.
        self.assertEqual(p.critical_order({0:p.Q(1)},{1:p.Q(1)},{1:p.Q(1)},{2:p.Q(1),3:p.Q(1)}),3)
        with self.assertRaises(ArithmeticError):p.critical_order({0:1},{1:1},{1:1},{2:1})

    def test_jet_split_negative_controls(self):
        self.assertTrue(p.jet_split([{'1','S'},{'2'}]))
        self.assertTrue(p.jet_split([{'2'},{'1'}]))
        for wrong in ([{'H'},{'2'}],[{'H'},{'H'}],[{'1'},{'1'}],[{'1','2'},{'2'}],[{'1'}]):
            self.assertFalse(p.jet_split(wrong))

    def test_factorisation_without_labels(self):
        f={(1,0,1,0):p.Q(1),(1,0,0,1):p.Q(-1),(0,1,1,0):p.Q(-1),(0,1,0,1):p.Q(1)}
        factors=p.factors(f)
        self.assertEqual(len(factors),2)
        self.assertTrue(all(len(p.witnesses(q))==2 for q in factors))
        self.assertFalse(set(p.support(factors[0])) & set(p.support(factors[1])))

    def test_graph_routing_and_box(self):
        self.assertEqual(len(all_graphs()),54)
        self.assertEqual(e.box_check()['scalar_power_d4'],-2)
        with self.assertRaises(ValueError):Graph('a1|a1')

    def test_production_symanzik_matches_independent_symbolic_polynomials(self):
        # Validate the *production* scalar-polynomial implementation, not just
        # forest-versus-matrix identities within the independent checker.
        point={e.x:s.Rational(3,7), e.y:s.Rational(2,5)}
        def frac(z): return s.Rational(z.numerator,z.denominator)
        for g in all_graphs():
            aa,ue,fe=e.expressions(g)
            u,f,_=p.symanzik(g.name)
            def mon(m): return s.prod(a**k for a,k in zip(aa,m))
            up=s.expand(sum(frac(c)*mon(m) for m,c in u.items()))
            fp=s.expand(sum((frac(c0)+e.t*frac(c2))*mon(m) for m,(c0,c2) in f.items()))
            self.assertEqual(s.expand(ue-up),0,g.name+' U')
            self.assertEqual(s.expand(fe.subs(point)-fp),0,g.name+' F')

    def test_stationary_positive_control_not_excluded(self):
        class Toy: name='positive interior control'
        aa=s.symbols('a:4')
        F=(aa[0]-aa[1])*(aa[2]-aa[3])
        with self.assertRaises(ArithmeticError):e.no_stationary(Toy(),F,aa)

    def test_wrong_polynomial_fails_identity(self):
        g=Graph('k|rb');aa,U,F=e.expressions(g);Um,Fm=e.matrix_polynomials(g,aa)
        self.assertEqual(s.expand(F-Fm),0)
        with self.assertRaises(ArithmeticError):p.require(s.expand(F+aa[0]*aa[1]*aa[2]-Fm)==0,'mutated polynomial')

    def test_corrupt_input_hash_fails(self):
        with tempfile.TemporaryDirectory() as d:
            Path(d,'vectors.json.xz').write_bytes(b'not the input data')
            with patch.object(p,'HERE',Path(d)):
                with self.assertRaisesRegex(ArithmeticError,'hash mismatch'):p.load_rows()

    def test_changed_power_and_classifier_fail(self):
        pack,rows=p.load_rows()
        bad=copy.deepcopy(rows)
        index=next(i for i,r in enumerate(bad) if r['old_omega'] is not None and r['kind']=='D')
        bad[index]['old_omega']+=1
        with patch.object(p,'load_rows',return_value=(pack,bad)):
            with self.assertRaisesRegex(ArithmeticError,'local saddle power changed'):p.run()
        with patch.object(p,'jet_split',return_value=True):
            with self.assertRaisesRegex(ArithmeticError,'SPS_rows_passing_jet_split'):p.run()

    def test_replay_across_hash_seed_and_optimisation(self):
        out=[]
        for seed,flags in [('1',[]),('987',['-O'])]:
            env=dict(os.environ,PYTHONHASHSEED=seed)
            run=subprocess.run([sys.executable,*flags,str(p.HERE/'param_space.py')],capture_output=True,text=True,env=env,timeout=180)
            self.assertEqual(run.returncode,0,run.stderr)
            out.append(run.stdout)
        self.assertEqual(out[0],out[1])
        self.assertEqual(json.loads(out[0])['selected_scalar_power_checks'],127)
        self.assertEqual(json.loads(out[0])['scalar_power_checks'],142)

if __name__=='__main__':unittest.main(verbosity=2)
