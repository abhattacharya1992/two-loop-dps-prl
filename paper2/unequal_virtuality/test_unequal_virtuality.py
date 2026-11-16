#!/usr/bin/env python3
"""Offline checks: kinematics, polarizations, crossings, legacy limit, fits, cache."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import mpmath as mp
import sympy as sp
from uv_core import *
from run_scan import run_point, validate_record
from analyse_scan import fit, analyse
from zz_points import all_coefficients as equal_complete, helicity_amplitudes as equal_contract


class Tests(unittest.TestCase):
    def setUp(self):
        self.context=mp.workdps(100)
        self.context.__enter__()
        self.config=configuration("15/2","1","2")
        self.kin=kinematics(self.config,6)

    def tearDown(self):
        self.context.__exit__(None,None,None)

    def near(self,a,b,tol="1e-85"):
        self.assertLess(abs(a-b),mp.mpf(tol)*max(mp.mpf(1),abs(a),abs(b)))

    def test_invariant_relations_and_sudakov_fractions(self):
        for c in (self.config,configuration("17","2","3"),configuration("17","3","2")):
            q=kinematics(c,8)
            self.near(q['s']+q['t']+q['u'],q['ma2']+q['mb2'])
            self.near(q['s']*q['pT2'],q['t']*q['u']-q['ma2']*q['mb2'])
            self.near(q['ma2'],q['s']*q['x1']*q['x2'])
            self.near(q['mb2'],q['s']*(1-q['x1'])*(1-q['x2']))
            self.near(q['st']**2+q['ct']**2,1)

    def test_four_momenta_and_all_polarizations(self):
        for h in ALL_HELICITIES:
            ps,es=vectors(self.kin,h)
            self.near(dot(ps[0],ps[0]),0)
            self.near(dot(ps[2],ps[2]),1)
            self.near(dot(ps[3],ps[3]),2)
            for j in range(4):
                self.near(ps[0][j]+ps[1][j],ps[2][j]+ps[3][j])
            for p,e in zip(ps,es):
                self.near(dot(p,e),0)
                self.near(dot(e,[mp.conj(z) for z in e]),-1)
            self.near(dot(es[0],ps[1]),0)
            self.near(dot(es[1],ps[0]),0)
            dp=[ps[0][j]-ps[2][j] for j in range(4)]
            self.near(dot(dp,dp),self.kin['t'])

    def test_scale_covariance(self):
        q=kinematics(configuration("750","100","200"),6)
        self.near(q['G'],self.kin['G'])
        for name in ('s','t','u','ma2','mb2','pT2','root'):
            self.near(q[name],100*self.kin[name])
        p=predictions(self.config)
        pp=predictions(configuration("750","100","200"))
        self.near(p['L'],pp['L'])

    def test_mass_exchange_kinematics(self):
        q=kinematics(configuration("15/2","2","1"),6)
        for n in ('root','t','u','pT2','G'):
            self.near(q[n],self.kin[n])
        self.near(q['E3'],self.kin['E4'])
        self.near(q['x1'],1-self.kin['x2'])

    def test_reject_bad_inputs(self):
        for c in (("5","1","2"),("4","1","1"),("10","0","1"),("10","-1","2")):
            with self.assertRaises(ValueError):
                configuration(*c)
        with self.assertRaises(ValueError):
            canonical(0.1)
        with self.assertRaises(ValueError):
            kinematics(self.config,0)
        with self.assertRaises(ValueError):
            vectors(self.kin,"0000")

    def test_equal_mass_crossing_reduction(self):
        a={j:mp.mpf(j) for j in DIRECT}
        b={j:mp.mpf(j*j) for j in DIRECT}
        self.assertEqual(complete_coefficients(a,b,b,a),equal_complete(a,b))

    def test_unequal_mass_crossing_slots(self):
        c=[{j:1000*n+j for j in DIRECT} for n in range(4)]
        a=complete_coefficients(*c)
        self.assertEqual(a[10],-2008)
        self.assertEqual(a[14],-3009)
        self.assertEqual(a[19],1018)
        # Negative control: deliberately applying the equal-mass shortcut changes these.
        wrong=complete_coefficients(c[0],c[1],c[1],c[0])
        self.assertNotEqual(wrong[10],a[10])
        self.assertNotEqual(wrong[14],a[14])

    def test_equal_mass_stored_amplitudes_all_36(self):
        f=HERE.parent/'zz/points/zz_b1o2_d4_p50.json'
        r=json.loads(f.read_text())
        q=kinematics(configuration("16/3","1","1"),4)
        A={p:{int(j):complex_pair(v) for j,v in vals.items()} for p,vals in r['A'].items()}
        out=contract(A,q,ALL_HELICITIES)
        old=equal_contract(A,q)
        for h in ALL_HELICITIES:
            for p in PARTS:
                self.near(out[h][p],old[h][p])
                self.near(out[h][p],complex_pair(r['M'][h][p]),'1e-35')

    def test_symbolic_scale_shift_and_colour_split(self):
        cf,ca,la,lb,ls=sp.symbols('CF CA loga logb logs')
        dc=2*cf-ca
        L=(la+lb)/2-ls
        shift=(la+lb-ls-(cf*(la+lb)-ca*ls)/dc)/2
        self.assertEqual(sp.simplify(shift+(ca-cf)*L/dc),0)
        f2=-2*dc*shift-3*cf
        self.assertEqual(sp.simplify(f2-(cf*(-2*L-3)+ca*2*L)),0)

    def test_prediction_equal_limit_and_mass_symmetry(self):
        p=predictions(configuration("16/3","1","1"))
        L=mp.log(mp.mpf(3)/16)
        self.near(p['2l_CF'][2],-2*L-3)
        a=predictions(self.config)
        b=predictions(configuration("15/2","2","1"))
        self.near(a['L'],b['L'])
        self.near(a['2l_CF'][3],-mp.mpf(4)/3)

    def synthetic(self):
        records=[]
        c1=mp.mpc(2,-7)
        for k in range(4,13):
            q=kinematics(self.config,k)
            G,rho=q['G'],q['pT2']/q['q2']
            # Deliberately NOT the theoretical prediction: verifies unconstrained fitting.
            v=c1*(mp.mpf('1.2345')*G**3-mp.mpf('2.3456')*G**2+G+2)+rho*(1+2*G+3*G**2+4*G**3)
            kin={n:cpair(z) if n=='G' else mp.nstr(z,100) for n,z in q.items()}
            records.append(dict(configuration=self.config,k=k,dps=90,kinematics=kin,
                         M={'++++':{'1l':cpair(c1*G+5+rho*(1+G+G**2)),
                                    '2l_CF':cpair(v),'2l_CA':cpair(2*v)}}))
        return records,c1

    def test_unconstrained_fit_and_wrong_prediction_negative_control(self):
        r,c1=self.synthetic()
        out=analyse(r)
        f=next(x for x in out['fits'] if x['part']=='2l_CF' and x['power_correction_degree']==3)
        self.near(complex_pair(f['comparison']['3']['extracted']),mp.mpf('1.2345'),'1e-32')
        self.assertGreater(mp.mpf(f['comparison']['3']['absolute_deviation']),1)
        v=[complex_pair(x['M']['++++']['2l_CF']) for x in r]
        ls=fit(r,v,3,3,extra=1)
        self.near(ls['coeff'][3]/c1,mp.mpf('1.2345'),'1e-80')

    def test_inadequate_or_corrupt_fit_rejected(self):
        r,_=self.synthetic()
        with self.assertRaises(ValueError):
            analyse(r[:4])
        r[0]['kinematics']['t']='-99'
        with self.assertRaises(ValueError):
            analyse(r)

    def test_cache_resume_invalidation_and_failure_preserves_output(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); build=root/'build'; build.mkdir()
            for name in ('p0002.py','algebra.marshal','rationals.json','coefficients.json'):
                (build/name).write_text('synthetic')
            class Fake: pass
            g=Fake();g.build=build;g.input_hashes=provenance.file_hashes(translated_inputs(build))
            def fake_eval(g,q,hs):
                H={h:{p:mp.mpc(1,2) for p in PARTS} for h in hs}
                A={p:{j:mp.mpc(j) for j in range(1,21)} for p in PARTS}
                return H,A,{'li22_bad_order':0}
            with patch('run_scan.evaluate',side_effect=fake_eval) as mocked:
                path,_=run_point(g,self.config,4,50,root)
                before=Path(path).read_bytes()
                self.assertEqual(run_point(g,self.config,4,50,root)[1],'cached')
                self.assertEqual(mocked.call_count,1)
            with patch('run_scan.evaluate',side_effect=RuntimeError('intentional')):
                with self.assertRaises(RuntimeError):
                    run_point(g,self.config,4,50,root,recompute=True)
            self.assertEqual(Path(path).read_bytes(),before)
            (build/'rationals.json').write_text('changed')
            with self.assertRaises(provenance.CacheMismatch):
                run_point(g,self.config,4,50,root)


if __name__=='__main__':
    unittest.main(verbosity=2)
