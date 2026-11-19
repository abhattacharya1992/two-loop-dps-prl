"""Precision sweep invariants and archived-reference consistency."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import precision_sweep as sweep
import full_performance as perf


class PrecisionTests(unittest.TestCase):
    def test_increasing_levels(self):
        self.assertEqual(sweep.parse_precisions("40,60,80"), (40,60,80))
        for s in ("", "40", "40,40", "60,40", "10,30",
                  "40,130", "40,x", "20,40,60,80,100,120"):
            with self.subTest(s=s), self.assertRaises(ValueError):
                sweep.parse_precisions(s)

    def test_low_precision_must_agree_with_highest_to_qualify(self):
        points = [{"A1": mp.mpc("1.00001"), "A2": mp.mpc("2")},
                  {"A1": mp.mpc("1.00000000001"), "A2": mp.mpc("2")},
                  {"A1": mp.mpc("1"), "A2": mp.mpc("2")}]
        summary = sweep.numeric_summary(points, [2.,4.,7.],
                                        (40,60,80), 1e-9)
        self.assertEqual(summary["component_count"],2)
        self.assertEqual(summary["accuracy_preserving_digits"],60)
        self.assertEqual(summary["accuracy_preserving_seconds"],4.)
        self.assertEqual(summary["raw_baseline_seconds"],2.)
        self.assertAlmostEqual(
            summary["raw_over_accuracy_preserving_runtime_ratio"], .5)

    def test_last_two_dps_must_converge(self):
        points = [{"a": mp.mpc(1.05)}, {"a": mp.mpc(1.01)},
                  {"a": mp.mpc(1)}]
        with self.assertRaisesRegex(ValueError,"do not converge"):
            sweep.numeric_summary(points,[1.,2.,3.],(40,60,80),1e-9)

    def test_independent_archive_check_fails_for_bad_component(self):
        good = {"++++:2l_CF": mp.mpc(1)}
        close = {"++++:2l_CF": mp.mpc("1.00000000001")}
        self.assertTrue(sweep._close_against_archive(close,good,1e-8)["pass"])
        with self.assertRaisesRegex(ValueError,"FAILED"):
            sweep._close_against_archive(
                {"++++:2l_CF": mp.mpc(2)},good,1e-8)
        with self.assertRaisesRegex(ValueError,"different component"):
            sweep._close_against_archive({"wrong":mp.mpc(1)},good,1e-8)

    def test_original_stage1_point_must_match(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"stage1.json"
            payload={"status":"evaluated_published_full_amplitude",
                     "process":"zz",
                     "configuration":{"beta":"1/10","cos_theta":"1/3","digits":40},
                     "verified_external_inputs":{"test":{"sha256":"1"*64}}}
            perf.write(path,perf.seal(payload))
            orig,digest=sweep.verify_stage1(
                path,"zz",{"beta":"1/10","cos_theta":"1/3"})
            self.assertEqual(orig["process"],"zz")
            self.assertEqual(digest,hashlib.sha256(path.read_bytes()).hexdigest())
            with self.assertRaisesRegex(ValueError,"disagree"):
                sweep.verify_stage1(path,"zz",
                                    {"beta":"1/10","cos_theta":"2/3"})

    def test_stage2_rejects_changed_precision_certificate(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"precision.json"
            payload={"status":"precision_convergence_verified",
                     "process":"hh", "original_stage1_sha256":"0"*64,
                     "precision":{"precision_stability_passed":True},
                     "archived_reference":{"pass":True}}
            perf.write(path,perf.seal(payload))
            self.assertEqual(
                perf.verified_precision_report(path,"hh","0"*64)["process"],"hh")
            with self.assertRaisesRegex(ValueError,"Unverified"):
                perf.verified_precision_report(path,"hh","1"*64)
            raw=json.loads(path.read_text())
            raw["payload"]["archived_reference"]["pass"]=False
            path.write_text(json.dumps(raw))
            with self.assertRaisesRegex(ValueError,"integrity"):
                perf.verified_precision_report(path,"hh","0"*64)

if __name__=="__main__":
    unittest.main()
