"""Unit tests for the shared two-stage HH/AA/ZZ pipeline."""
from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import unified_pipeline as pipeline


def valid_payload():
    rows = {
        p: {"process": p, "storage": {
            "verified_byte_identical": True,
            "records": 1, "raw_bytes": 1000, "compressed_bytes": 500,
            "percent_smaller": 50.0,
            "raw_json_parse_median_ms": 0.04,
            "compressed_verify_and_parse_median_ms": 0.07,
        }, "analytic": {"status": "not_derived"},
            "original_exact": {"status": "not_run"}}
        for p in pipeline.PROCESSES
    }
    rows["AA"]["analytic"] = {
        "status": "evaluated_validated_analytic_series",
        "median_seconds": {"beta4": 0.002, "beta6": 0.003},
        "preparation_seconds": 0.01,
        "series_input_bytes": 100,
    }
    rows["AA"]["historical_exact_comparison"] = {
        "expression_size_ratio_exact_over_beta6": 70,
        "timing_rows": [{"reported_warm_speedup": 19.1},
                        {"reported_warm_speedup": 18.2}],
        "validation": {"largest_relative_component_error_beta6": 0.00015},
    }
    rows["AA"]["original_exact"] = {"status": "historical_comparison_only"}
    return {
        "repository_commit": "a"*40,
        "aa_repository_commit": "b"*40,
        "configuration": {"beta": "1/10", "run_exact": False},
        "peak_rss_kib_process": 10000,
        "memory_scope": "whole process",
        "processes": rows,
    }


class IntegrityTests(unittest.TestCase):
    def test_payload_roundtrip_and_report(self):
        payload = valid_payload()
        sealed = pipeline.seal(payload)
        self.assertEqual(pipeline.verify_record(sealed), payload)
        report = pipeline.summary(payload)
        self.assertEqual(report["processes"]["HH"]["analytical_status"], "not_derived")
        self.assertEqual(report["processes"]["AA"]["historical_warm_speedup_range"], [18.2, 19.1])
        self.assertIn("not_derived", pipeline.markdown(report))

    def test_tamper_rejected_even_if_schema_unchanged(self):
        sealed = pipeline.seal(valid_payload())
        sealed["payload"]["processes"]["ZZ"]["storage"]["compressed_bytes"] = 1
        with self.assertRaisesRegex(ValueError, "changed"):
            pipeline.verify_record(sealed)

    def test_missing_process_rejected(self):
        sealed = pipeline.seal(valid_payload())
        del sealed["payload"]["processes"]["HH"]
        with self.assertRaisesRegex(ValueError, "missing"):
            pipeline.verify_record(sealed)

    def test_unverified_data_rejected(self):
        payload = valid_payload()
        payload["processes"]["AA"]["storage"]["verified_byte_identical"] = False
        with self.assertRaisesRegex(ValueError, "Unverified"):
            pipeline.verify_record(pipeline.seal(payload))

    def test_stage2_is_independent_and_does_not_modify_stage1(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            input_file = root/"stage1.json"
            input_file.write_text(json.dumps(pipeline.seal(valid_payload())))
            before = input_file.read_bytes()
            with patch.object(pipeline, "stage1", side_effect=AssertionError("must not rerun")):
                report = pipeline.stage2(input_file, root/"out")
            self.assertEqual(input_file.read_bytes(), before)
            self.assertEqual(report["processes"]["HH"]["point_records"], 1)
            self.assertTrue((root/"out"/"report.md").is_file())
            self.assertTrue((root/"out"/"stage2.json").is_file())

    def test_stage2_rejects_modified_stage1(self):
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp)/"stage1.json"
            record = pipeline.seal(valid_payload())
            record["payload"]["processes"]["AA"]["process"] = "ZZ"
            source.write_text(json.dumps(record))
            with self.assertRaises(ValueError):
                pipeline.stage2(source, Path(temp)/"out")


if __name__ == "__main__":
    unittest.main()
