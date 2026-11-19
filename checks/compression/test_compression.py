"""Fail-closed checks for HH/AA/ZZ compression and AA historical results."""
import csv
import json
from pathlib import Path
import tempfile
import unittest
import aa_threshold as aa
import point_benchmark as point

class PointTests(unittest.TestCase):
    def test_round_trip_byte_identical(self):
        raw = b'{"beta":"1/2","A":{"x":["1.000000","-0.0"]}}'
        packed = point.compress_point(raw)
        self.assertEqual(point.decompress_point(packed), raw)

    def test_corruption_rejected(self):
        raw = b'{"beta":"1/2","A":{"x":["1.000000","-0.0"]}}'
        packed = point.compress_point(raw)
        variants = [packed[:-1], packed+b"junk",
                    packed[:len(point.MAGIC)+8]+b"\0"*32+packed[point.HEADER_LENGTH:],
                    b"wrong"+packed]
        for item in variants:
            with self.subTest(size=len(item)), self.assertRaises(ValueError):
                point.decompress_point(item)

    def test_all_process_record_types(self):
        with tempfile.TemporaryDirectory() as d:
            for process, key in point.AMPLITUDE_KEY.items():
                path = Path(d)/f"{process}.json"
                path.write_text(json.dumps({"beta":"1/10",key:{"x":["0.00000001","-0.0"]}}))
                out = point.measure(process,[path],2,9)
                self.assertEqual(out["records"],1)
                self.assertTrue(out["verified_byte_identical"])

class AATests(unittest.TestCase):
    def test_complex_relative_error(self):
        self.assertAlmostEqual(aa.relative_error(1+1j,1+1.01j),0.01/(2**0.5))
        self.assertEqual(aa.relative_error(0,0.3),0.3)

    def test_historical_fixture(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            (root/"data/timing").mkdir(parents=True)
            (root/"data/validation").mkdir(parents=True)
            (root/"data/timing/timing_run.json").write_text(json.dumps({
                "representation_bytes":{"unexpanded":1000,"beta6":20},
                "environment":{"cpu":"fixture"}}))
            cols=["Beta","CosTheta","WorkingPrecision","WarmRepeats",
                  "ExactFirstSeconds","SeriesFirstSeconds","ExactWarmMedianSeconds",
                  "SeriesWarmMedianSeconds","FirstCallSpeedup","WarmSpeedup",
                  "P1RelativeError","P2RelativeError"]
            with (root/"data/timing/timing_benchmark.csv").open("w",newline="") as f:
                out=csv.DictWriter(f,fieldnames=cols);out.writeheader()
                out.writerow(dict(zip(cols,["1/10","0",50,3,20,2,2,0.1,10,20,0.0001,0.001])))
            cols=["quantity","beta","costheta","re_exact","im_exact",
                  "re_series4","im_series4","re_series6","im_series6"]
            with (root/"data/validation/validation_1080.csv").open("w",newline="") as f:
                out=csv.DictWriter(f,fieldnames=cols);out.writeheader()
                out.writerow(dict(zip(cols,["M1_nf0","1/10","0","2","0","1.8","0","1.99","0"])))
            report=aa.historical(root)
            self.assertEqual(report["expression_size_ratio_exact_over_beta6"],50)
            self.assertEqual(report["validation"]["rows_beta_le_0_2"],1)
            self.assertAlmostEqual(report["validation"]["largest_relative_component_error_beta6"],0.005)
            self.assertEqual(len(report["input_sha256"]),3)

    def test_missing_exact_inputs_fail(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(FileNotFoundError):
                aa.fresh(Path(d),beta="1/10",costheta="0",digits=50,
                         repeats=1,mA="125",muR="125",Nc=3,nf=5,order=6)

if __name__=="__main__":
    unittest.main()
