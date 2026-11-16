#!/usr/bin/env python3
"""Exercise cache decisions with temporary synthetic inputs, without VVamp data."""
import json
import marshal
import pathlib
import tempfile
import unittest
from unittest.mock import Mock, patch

import mpmath as mp

from ggvv import GGVV
import zz_cache
import zz_points
import zz_routeB
import zz_wide

provenance = zz_cache.cache_provenance


class CacheTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # GGVV's real loader receives only synthetic, no-op translated functions.
        definitions = "\n".join(f"def p{k:04d}(r, Q): pass" for k in range(3, 2464))
        cls.algebra = marshal.dumps([compile(definitions, "synthetic-algebra", "exec")])

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = pathlib.Path(self.temporary.name)
        self.previous_dps = mp.mp.dps
        self.build = self.root / "build"
        self.build.mkdir()
        (self.build / "p0002.py").write_text("def p0002(r, Q, F): pass\n")
        (self.build / "algebra.marshal").write_bytes(self.algebra)
        (self.build / "rationals.json").write_text("{}")
        (self.build / "coefficients.json").write_text("{}")
        self.g = self.load_evaluator()
        self.inputs = self.root / "inputs"
        self.inputs.mkdir()
        (self.inputs / "gg_Aj_A_qT-1.0.inc").write_text("\n".join(
            f"g [A{j}_A_{loop}l_qT] = {j};"
            for loop in (1, 2) for j in range(1, 21)))
        (self.inputs / "vvpmasters_opt-1.0.m").write_text('M["P12", 1, 0] -> 1\n,\n')

    def tearDown(self):
        mp.mp.dps = self.previous_dps
        self.temporary.cleanup()

    def load_evaluator(self):
        g = GGVV(self.build)
        g.compute = Mock(side_effect=lambda *args: {
            j: {part: mp.mpf(j) for part in zz_points.PARTS} for j in zz_points.DIRECT})
        return g

    def point(self, kind, out, *, recompute=False, dps=40):
        if kind == "forward":
            return zz_points.run_point(self.g, "1/2", 4, dps, out, recompute=recompute)
        return zz_wide.run_point(self.g, "1/2", 1, out, recompute=recompute)

    def route(self, out, parts=("1l",), *, recompute=False):
        return zz_routeB.evaluate(self.inputs, out, "1/2", 4, 40, list(parts), recompute=recompute)

    def test_matching_forward_and_wide_results_skip_evaluation(self):
        for kind in ("forward", "wide"):
            with self.subTest(kind=kind):
                out = self.root / kind
                path, _ = self.point(kind, out)
                record = json.loads(path.read_text())
                self.assertEqual(record["provenance"]["inputs"], self.g.input_hashes)
                self.g.compute.reset_mock()
                before = path.read_bytes()
                self.assertEqual(self.point(kind, out), (path, "cached"))
                self.g.compute.assert_not_called()
                self.assertEqual(path.read_bytes(), before)

    def test_legacy_points_require_explicit_recomputation(self):
        for kind in ("forward", "wide"):
            with self.subTest(kind=kind):
                out = self.root / kind
                path, _ = self.point(kind, out)
                record = json.loads(path.read_text())
                del record["provenance"]
                path.write_text(json.dumps(record))
                legacy = path.read_bytes()
                self.g.compute.reset_mock()
                with self.assertRaisesRegex(provenance.CacheMismatch, "legacy"):
                    self.point(kind, out)
                self.g.compute.assert_not_called()
                self.assertEqual(path.read_bytes(), legacy)
                self.point(kind, out, recompute=True)
                self.assertIn("provenance", json.loads(path.read_text()))

    def test_tagged_but_incomplete_or_nonfinite_payload_is_rejected(self):
        for kind in ("forward", "wide", "route_b"):
            with self.subTest(kind=kind):
                out = self.root / kind
                call = (lambda: self.route(out)) if kind == "route_b" else (lambda: self.point(kind, out))
                path, _ = call()
                original = json.loads(path.read_text())
                malformed = dict(original)
                malformed.pop("M")
                path.write_text(json.dumps(malformed))
                with self.assertRaisesRegex(provenance.CacheMismatch, "incomplete scientific result"):
                    call()
                malformed = json.loads(json.dumps(original))
                amplitude = next(iter(malformed["M"].values()))
                amplitude[next(iter(amplitude))][0] = "nan"
                path.write_text(json.dumps(malformed))
                with self.assertRaisesRegex(provenance.CacheMismatch, "finite real number"):
                    call()

    def test_failed_recomputation_preserves_old_point(self):
        for kind in ("forward", "wide"):
            with self.subTest(kind=kind):
                self.g = self.load_evaluator()
                out = self.root / kind
                path, _ = self.point(kind, out)
                before = path.read_bytes()
                self.g.compute.side_effect = RuntimeError("synthetic evaluation failure")
                with self.assertRaisesRegex(RuntimeError, "synthetic evaluation failure"):
                    self.point(kind, out, recompute=True)
                self.assertEqual(path.read_bytes(), before)

    def test_changed_translated_bytes_reject_loaded_state_and_old_output(self):
        out = self.root / "points"
        path, _ = self.point("forward", out)
        before = path.read_bytes()
        (self.build / "rationals.json").write_text('{"1/1": 0}')
        with self.assertRaisesRegex(provenance.CacheMismatch, "after GGVV was loaded"):
            self.point("forward", out)
        self.g = self.load_evaluator()
        with self.assertRaisesRegex(provenance.CacheMismatch, "changed inputs"):
            self.point("forward", out)
        self.g.compute.assert_not_called()
        self.assertEqual(path.read_bytes(), before)
        # A fresh directory accepts the changed inputs without relabelling the old result.
        fresh, _ = self.point("forward", self.root / "fresh")
        self.assertNotEqual(json.loads(fresh.read_text())["provenance"]["inputs"],
                            json.loads(before)["provenance"]["inputs"])

    def test_precision_is_checked_even_if_a_file_is_renamed(self):
        out = self.root / "points"
        path, _ = self.point("forward", out)
        renamed = out / path.name.replace("_p40.json", "_p45.json")
        renamed.write_bytes(path.read_bytes())
        with self.assertRaises(provenance.CacheMismatch):
            self.point("forward", out, dps=45)
        self.point("forward", out, dps=45, recompute=True)
        self.assertEqual(json.loads(renamed.read_text())["provenance"]["parameters"]["dps"], 45)

    def test_evaluator_source_and_runtime_changes_invalidate_results(self):
        synthetic_source = self.root / "evaluator.py"
        synthetic_source.write_text("version = 1\n")
        original_sources = zz_cache.source_files
        def sources(*names):
            return dict(original_sources(*names), synthetic_evaluator=synthetic_source)
        with patch.object(zz_cache, "source_files", side_effect=sources):
            out = self.root / "points"
            path, _ = self.point("forward", out)
            before = path.read_bytes()
            synthetic_source.write_text("version = 2\n")
            with self.assertRaises(provenance.CacheMismatch):
                self.point("forward", out)
            synthetic_source.write_text("version = 1\n")
            with patch.object(provenance.platform, "python_version", return_value="changed-runtime"):
                with self.assertRaises(provenance.CacheMismatch):
                    self.point("forward", out)
            self.assertEqual(path.read_bytes(), before)

    def test_split_cache_checks_inputs_and_every_generated_file(self):
        cache = self.root / "_split"
        zz_routeB.split_inputs(self.inputs, cache)
        manifest = (cache / "manifest.json").read_bytes()
        zz_routeB.split_inputs(self.inputs, cache)
        self.assertEqual((cache / "manifest.json").read_bytes(), manifest)
        coefficient = cache / "A1_A_1l_qT.txt"
        coefficient.write_text("altered")
        with self.assertRaisesRegex(provenance.CacheMismatch, "changed split contents"):
            zz_routeB.split_inputs(self.inputs, cache)
        zz_routeB.split_inputs(self.inputs, cache, recompute=True)
        coefficient.unlink()
        with self.assertRaisesRegex(provenance.CacheMismatch, "changed split files"):
            zz_routeB.split_inputs(self.inputs, cache)
        zz_routeB.split_inputs(self.inputs, cache, recompute=True)
        (self.inputs / "vvpmasters_opt-1.0.m").write_text('M["P12", 1, 0] -> 2\n,\n')
        with self.assertRaisesRegex(provenance.CacheMismatch, "changed inputs"):
            zz_routeB.split_inputs(self.inputs, cache)
        zz_routeB.split_inputs(self.inputs, cache, recompute=True)
        self.assertEqual(json.loads((cache / "rules.json").read_text())["P12,1,0"], "2")

    def test_master_mapping_change_invalidates_route_b_point_and_split_cache(self):
        mapping=self.root/'mapping.py'
        mapping.write_text('crossing = "current"\n')
        original=zz_cache.source_files
        def sources(*names):
            files=original(*names)
            if 'zz_master_families.py' in files:
                files['zz_master_families.py']=mapping
            return files
        with patch.object(zz_cache,'source_files',side_effect=sources):
            out=self.root/'route_b'
            path,_=self.route(out)
            before=path.read_bytes()
            mapping.write_text('crossing = "previous"\n')
            with self.assertRaises(provenance.CacheMismatch):self.route(out)
            with self.assertRaises(provenance.CacheMismatch):
                zz_routeB.split_inputs(self.inputs,out/'_split')
            self.assertEqual(path.read_bytes(),before)

    def test_legacy_split_marker_requires_a_clean_rebuild(self):
        cache = self.root / "_split"
        cache.mkdir()
        (cache / "done").write_text("ok")
        (cache / "rules.json").write_text('{"old": "record"}')
        with self.assertRaisesRegex(provenance.CacheMismatch, "legacy split cache"):
            zz_routeB.split_inputs(self.inputs, cache)
        self.assertEqual((cache / "rules.json").read_text(), '{"old": "record"}')
        zz_routeB.split_inputs(self.inputs, cache, recompute=True)
        self.assertFalse((cache / "done").exists())
        self.assertTrue((cache / "manifest.json").is_file())

    def test_failed_split_rebuild_preserves_previous_cache(self):
        cache = self.root / "_split"
        zz_routeB.split_inputs(self.inputs, cache)
        before = {p.name: p.read_bytes() for p in cache.iterdir()}
        (self.inputs / "gg_Aj_A_qT-1.0.inc").write_text("invalid FORM input")
        with self.assertRaises(ValueError):
            zz_routeB.split_inputs(self.inputs, cache, recompute=True)
        self.assertEqual({p.name: p.read_bytes() for p in cache.iterdir()}, before)

    def test_route_b_point_reuse_checks_parts_and_legacy_status(self):
        out = self.root / "route"
        path, _ = self.route(out)
        self.assertEqual(self.route(out), (path, "cached"))
        before = path.read_bytes()
        with self.assertRaises(provenance.CacheMismatch):
            self.route(out, ("2l_CF",))
        self.assertEqual(path.read_bytes(), before)
        self.route(out, ("2l_CF",), recompute=True)
        record = json.loads(path.read_text())
        self.assertEqual(record["provenance"]["parameters"]["parts"], ["2l_CF"])
        del record["provenance"]
        path.write_text(json.dumps(record))
        legacy = path.read_bytes()
        with self.assertRaisesRegex(provenance.CacheMismatch, "legacy"):
            self.route(out, ("2l_CF",))
        self.assertEqual(path.read_bytes(), legacy)
        self.route(out, ("2l_CF",), recompute=True)
        self.assertIn("provenance", json.loads(path.read_text()))


if __name__ == "__main__":
    unittest.main()
