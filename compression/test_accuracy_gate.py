"""Fast accuracy and integrity tests for optional full-amplitude benchmarks."""
import json
import tempfile
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import full_performance as perf


class AccuracyGateTests(unittest.TestCase):
    def test_unreliable_aa_result_fails_and_leaves_diagnostics(self):
        measurement = {"relative_errors": {"M1": 1.0, "M2": 1.0},
                       "evaluation_median_seconds": {"exact": 150, "beta_series": .02}}
        archive = {"status": "compared",
                   "exact_relative_errors": {"M1": 1.0, "M2": 1.0}}
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "AA numerical accuracy failed"):
                perf.enforce_aa_accuracy(measurement, archive, Path(tmp))
            saved = json.loads((Path(tmp)/"accuracy_failure.json").read_text())
            self.assertFalse(saved["gate"]["passed"])
            self.assertEqual(saved["status"], "failed_numerical_accuracy")

    def test_independently_validated_aa_errors_pass(self):
        measurement = {"relative_errors": {"M1": 6.3816e-9, "M2": 8.3391e-8}}
        archive = {"status": "compared",
                   "exact_relative_errors": {"M1": 1e-15, "M2": 2e-14}}
        self.assertTrue(perf.check_aa_accuracy(measurement, archive)["passed"])

    def test_nan_inf_and_missing_components_rejected(self):
        independent = {"status": "not_tabulated", "exact_relative_errors": {}}
        for bad in ({"M1": float("nan"), "M2": 1e-9},
                    {"M1": float("inf"), "M2": 1e-9}):
            self.assertFalse(perf.check_aa_accuracy(
                {"relative_errors": bad}, independent)["passed"])
        with self.assertRaisesRegex(ValueError, "Missing AA"):
            perf.check_aa_accuracy({"relative_errors": {"M1": 1e-8}}, independent)

    def test_aa_singular_angle_disallowed(self):
        with self.assertRaisesRegex(ValueError, "removable 0/0"):
            perf.validated_args("aa", "1/10", "0", 120, 1)

    def test_bad_saved_aa_stage2_cannot_report_valid_speedup(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"stage1.json"
            payload = {
                "status": "evaluated_published_full_amplitude",
                "process": "aa",
                "verified_external_inputs": {"source": {"sha256": "1"*64}},
                "measurement": {
                    "relative_errors": {"M1": 1, "M2": 1},
                    "archived_exact_comparison": {
                        "status": "not_tabulated", "exact_relative_errors": {}},
                    "accuracy_gate": {"passed": False},
                },
            }
            perf.write(path, perf.seal(payload))
            with self.assertRaisesRegex(ValueError, "numerical accuracy gate"):
                perf.stage2(path, "aa", Path(tmp))


if __name__ == "__main__":
    unittest.main()
