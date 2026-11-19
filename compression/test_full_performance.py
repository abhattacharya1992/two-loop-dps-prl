"""Fast, deterministic guards for optional original full-amplitude benchmarks."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import full_performance as full


class FullPerformanceTests(unittest.TestCase):
    def test_sha256_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            inp = root/"inputs"
            inp.mkdir()
            good = hashlib.sha256(b"original").hexdigest()
            (inp/"SHA256SUMS").write_text(good + "  Finite_2Loop_P1\n")
            folder = inp/"arXiv-1909.08993/anc"
            folder.mkdir(parents=True)
            (folder/"Finite_2Loop_P1").write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                full.verified_sources("aa",root,root)

    def test_missing_inputs_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaises(FileNotFoundError):
                full.verified_sources("hh",Path(temp),Path(temp))

    def test_zz_pinned_hash_verification(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            inp=root/"paper2/zz/inputs"
            inp.mkdir(parents=True)
            digest={}
            for name in full.ZZ_REQUIRED:
                raw=("fixture-"+name).encode("utf-8")
                (inp/name).write_bytes(raw)
                digest[name]=hashlib.sha256(raw).hexdigest()
            (inp/"SHA256SUMS").write_text("".join(
                f"{hash}  {name}\n" for name,hash in digest.items()))
            verified=full.verified_sources("zz",root,root)
            self.assertEqual(set(verified),set(full.ZZ_REQUIRED))
            self.assertEqual(verified[full.ZZ_REQUIRED[0]]["sha256"],
                             digest[full.ZZ_REQUIRED[0]])

    def test_hh_stage2_reports_no_beta_expansion(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp)
            payload={
                "status":"evaluated_published_full_amplitude",
                "process":"hh","repo_commit":"a"*40,"aa_repo_commit":None,
                "configuration":{"beta":"1/10","cos_theta":"0","digits":40,"repeats":1},
                "verified_external_inputs":{"Gauge1Simp.in":{
                    "sha256":"b"*64,"bytes":100}},
                "measurement":{"evaluation_median_seconds":2.0,
                               "evaluation_seconds":[2.0],
                               "preparation_seconds":0.5,"input_total_bytes":100},
                "whole_process_peak_rss_kib":50000,
                "total_wall_seconds":2.9,"cpu_user_seconds":2.0,
                "cpu_system_seconds":0.4,
                "interpretation":"published HH class-B exact only",
            }
            file=folder/"stage1.json"
            full.write(file, full.seal(payload))
            out=full.stage2(file,"hh",folder/"out")
            self.assertEqual(out["exact_median_seconds"],2)
            self.assertEqual(out["analytic_beta_series_status"],"not_derived")
            self.assertNotIn("fresh_exact_over_beta6_speedup",out)
            self.assertTrue((folder/"out/report.md").is_file())

    def test_stage1_tamper_rejected(self):
        record=full.seal({"status":"evaluated_published_full_amplitude",
                          "process":"aa","verified_external_inputs":{
                              "source":{"sha256":"a"*64}}})
        with self.assertRaisesRegex(ValueError, "Wrong benchmark process"):
            full.unseal(record,"zz")
        record["payload"]["process"]="zz"
        with self.assertRaisesRegex(ValueError, "modified"):
            full.unseal(record,"zz")

    def test_fails_closed_for_missing_full_amplitudes(self):
        args=argparse.Namespace(beta="1/10",costheta="0",digits=40,
               repeats=1,aa_repo=Path("/nonexistent-aa"),
               repo=Path("/nonexistent-prl"),output=Path("/nonexistent"))
        with self.assertRaises(FileNotFoundError):
            full.stage1("hh",args)

    def test_common_validity_domain(self):
        bad=[
            ("invalid","1/10","0",40,1),
            ("aa","1/4","0",40,1),
            ("hh","0","0",40,1),
            ("zz","1/10","-1/10",40,1),
            ("zz","1/10","1",40,1),
            ("aa","1/10","0",121,1),
            ("hh","1/10","0",40,6),
        ]
        for p,b,c,d,r in bad:
            with self.subTest(process=p,beta=b):
                with self.assertRaises(ValueError):
                    full.validated_args(p,b,c,d,r)

if __name__=="__main__":
    unittest.main()
