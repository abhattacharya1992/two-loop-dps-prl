"""Unit tests for the two optional historical-benchmark report generators."""
import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark_size as sizes
import benchmark_pipeline as pipeline


class TreeSizeTests(unittest.TestCase):
    def test_nul_terminated_tree_and_space_in_name(self):
        h = b"a" * 40
        raw = (b"100644 blob " + h + b" 10\tfoo bar.py\0"
               + b"100644 blob " + h + b" 7\tresults.json\0"
               + b"120000 blob " + h + b" 4\tsymlink\0")
        self.assertEqual(sizes.parse_tree(raw), {
            "foo bar.py": 10, "results.json": 7, "symlink": 4})

    def test_disjoint_categories_and_code_not_same_as_total_size(self):
        files = {
            "src/analysis.py": 40,
            "Makefile": 10,
            "data/results.json": 100,
            "paper2/r0/theory/secdec_cache_x.pkl": 500,
            "paper2/manuscript_prl/letter.tex": 200,
            "logs/output.log": 30,
            "LICENSE": 20,
        }
        info = sizes.describe(files)
        self.assertEqual(info["tracked_files"], 7)
        self.assertEqual(info["tracked_bytes"], 900)
        self.assertEqual(info["categories"]["source_code"]["bytes"], 50)
        self.assertEqual(info["categories"]["numerical_data"]["bytes"], 100)
        self.assertEqual(info["categories"]["optional_sector_caches"]["bytes"], 500)
        self.assertEqual(sum(v["bytes"] for v in info["categories"].values()), 900)

    def test_reduction_signs_and_zero_reference(self):
        a = sizes.describe({"script.py": 100, "data.json": 200})
        b = sizes.describe({"script.py": 110, "data.json": 100})
        result = sizes.reduction(a, b)
        self.assertAlmostEqual(result["percent_smaller"], 30.0)
        self.assertEqual(result["source_code_bytes_difference"], -10)
        self.assertIsNone(sizes.reduction(sizes.describe({}), b)["percent_smaller"])


class PipelineTests(unittest.TestCase):
    def test_baseline_to_candidate_speedup(self):
        raw = {
            "shared_dps": {
                "baseline_median_seconds": 10.0,
                "candidate_median_seconds": 5.0,
                "peak_rss_kb": {"baseline": 200, "candidate": 210},
                "failed": False, "possible_regression": False,
            },
            "matching": {
                "baseline_median_seconds": 4.0,
                "candidate_median_seconds": 8.0,
                "peak_rss_kb": {"baseline": 200, "candidate": 205},
                "failed": False, "possible_regression": True,
            },
        }
        result = pipeline.make_summary(raw, "original_main")
        self.assertEqual(result["shared_dps"]["speedup_baseline_over_candidate"], 2.0)
        self.assertEqual(result["shared_dps"]["candidate_percent_faster"], 50.0)
        self.assertEqual(result["matching"]["speedup_baseline_over_candidate"], 0.5)
        self.assertEqual(result["matching"]["candidate_percent_faster"], -100.0)
        self.assertTrue(result["matching"]["possible_regression"])

    def test_scope_warning_is_in_report(self):
        report = {
            "baseline_label": "original_main",
            "baseline": "a" * 40,
            "candidate": "b" * 40,
            "python": "python",
            "platform": "Linux",
            "runs": 5,
            "summary": {},
            "errors": [],
        }
        doc = pipeline.markdown(report)
        self.assertIn("NOT time complete HH, AA or ZZ amplitude", doc)
        self.assertIn("All benchmark commands returned zero", doc)


if __name__ == "__main__":
    unittest.main()
