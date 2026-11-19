"""Kinematic and input guards for optional fresh HH and ZZ evaluations."""
import tempfile
from pathlib import Path
import unittest
import mpmath as mp
import exact_evaluators as exact

class ExactEvaluationTests(unittest.TestCase):
    def test_fixed_angle_kinematics(self):
        row=exact.fixed_angle_kinematics("1/10","1/3",digits=60)
        s,t,u,pt=[mp.mpf(row[k]) for k in ("s","t","u","pT2")]
        self.assertAlmostEqual(float(s+t+u),2.0)
        self.assertGreater(pt,0)
        self.assertAlmostEqual(float(s),400/99)
        self.assertAlmostEqual(float(mp.mpf(row["x"])),9/11)
        self.assertAlmostEqual(float(t-u),float(s)/30)

    def test_bad_kinematics_fail(self):
        for beta,cos in (("0","0"),("1","0"),("1/2","-1/3"),
                         ("1/2","1"),("2","0")):
            with self.subTest(beta=beta,cos=cos),self.assertRaises(ValueError):
                exact.fixed_angle_kinematics(beta,cos,digits=50)

    def test_missing_external_inputs_are_not_silently_replaced(self):
        kin=exact.fixed_angle_kinematics("1/10","0",digits=40)
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(FileNotFoundError):
                exact.exact_hh(kin,40,1,Path(d))
            with self.assertRaises(FileNotFoundError):
                exact.exact_zz(kin,40,1,Path(d),Path(d))
if __name__=="__main__":
    unittest.main()
