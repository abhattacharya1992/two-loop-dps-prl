"""Cache regressions using tiny inputs and mocked expensive amplitude evaluators."""
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

PKG = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PKG))
from dps import provenance
from dps.backends import aa
from dps import layout


class AmplitudeCacheTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.anc = self.root / 'inputs/arXiv-1909.08993/anc'
        self.anc.mkdir(parents=True)
        self.sources = self.root / 'src/python/ggaa'
        self.sources.mkdir(parents=True)
        for name in ('twoloop.py', 'li22.py', 'cache.py'):
            (self.sources / name).write_text('# fixture evaluator\n')
        self.files = {(1, 1): 'one', (1, 2): 'two'}
        (self.anc / 'one').write_text('1')
        (self.anc / 'two').write_text('2')
        self.manifest()
        owner = self
        self.calls = 0

        class Evaluator:
            def __init__(self, symbols):
                self._cache = {}

            def evaluate(self, expression):
                owner.calls += 1
                value = int(expression)
                return (aa.N + aa.NF if value == 1 else
                        value * (aa.N**2 + aa.N*aa.NF + aa.NF**2 + aa.NF/aa.N))

        evaluator = SimpleNamespace(FILES=self.files, Evaluator=Evaluator,
                                    read_expression=lambda path: path.read_text())
        patches = patch.multiple(aa, P1_ROOT=self.root, ANC=self.anc, p1_package=lambda: evaluator)
        patches.start()
        self.addCleanup(patches.stop)
        aa._TEXT.clear()
        self.out = self.root / 'points'

    def manifest(self):
        (self.root / 'inputs/SHA256SUMS').write_text(''.join(
            f'{provenance.sha256_file(self.anc / name)}  {name}\n' for name in self.files.values()))

    def evaluate(self, **options):
        return aa.evaluate(1, '1/2', 3, 40, out_dir=self.out, **options)

    def test_aa_checks_actual_input_before_cache_reuse(self):
        path, _ = self.evaluate()
        path = Path(path)
        original = path.read_bytes()
        first = json.loads(original)
        self.assertEqual(first['inputs']['two'], provenance.sha256_file(self.anc / 'two'))
        self.assertEqual(self.evaluate()[1], 'cached')
        self.assertEqual(self.calls, 2)
        (self.anc / 'two').write_text('3')
        with self.assertRaisesRegex(ValueError, 'actual SHA-256'):
            self.evaluate()
        self.assertEqual(path.read_bytes(), original)
        self.manifest()
        with self.assertRaises(provenance.CacheMismatch):
            self.evaluate()
        self.assertEqual(path.read_bytes(), original)
        self.evaluate(recompute=True)
        second = json.loads(path.read_text())
        self.assertNotEqual(second['A'], first['A'])  # in-memory parsed expression was invalidated too
        self.assertEqual(second['inputs']['two'], provenance.sha256_file(self.anc / 'two'))
        self.assertEqual(self.calls, 4)

    def test_aa_source_precision_and_job_parameters_invalidate(self):
        path, _ = self.evaluate()
        original = Path(path).read_bytes()
        (self.out / aa.file_name(1, '1/2', 3, 41)).write_bytes(original)
        with self.assertRaises(provenance.CacheMismatch):
            aa.evaluate(1, '1/2', 3, 41, out_dir=self.out)
        with self.assertRaises(provenance.CacheMismatch):
            self.evaluate(regular=False)
        (self.sources / 'li22.py').write_text('# changed evaluator algorithm\n')
        with self.assertRaises(provenance.CacheMismatch):
            self.evaluate()
        self.assertEqual(self.calls, 2)
        self.assertEqual(Path(path).read_bytes(), original)

    def test_aa_legacy_result_is_preserved_until_explicit_recompute(self):
        self.out.mkdir()
        path = self.out / aa.file_name(1, '1/2', 3, 40)
        path.write_text('{"A": {"N": ["9", "0"]}}')
        original = path.read_bytes()
        with self.assertRaisesRegex(provenance.CacheMismatch, 'legacy'):
            self.evaluate()
        self.assertEqual(self.calls, 0)
        self.assertEqual(path.read_bytes(), original)
        self.evaluate(recompute=True)
        self.assertIn('provenance', json.loads(path.read_text()))

    def test_aa_matching_metadata_does_not_certify_missing_or_nonfinite_results(self):
        path, _ = self.evaluate()
        path = Path(path)
        valid = json.loads(path.read_text())
        for invalid in ({'provenance': valid['provenance']},
                        {**valid, 'A': {**valid['A'], 'N2': ['nan', '0']}}):
            path.write_text(json.dumps(invalid))
            before = path.read_bytes()
            with self.assertRaisesRegex(provenance.CacheMismatch, 'malformed'):
                self.evaluate()
            self.assertEqual(path.read_bytes(), before)
        self.assertEqual(self.calls, 2)

    def test_hh_input_change_invalidates_result_and_compiled_expression(self):
        sys.path.insert(0, str(layout.path('hh_points').parent))
        import scan_colour_scale as scan
        input_dir = self.root / 'hh'
        (input_dir / 'anc').mkdir(parents=True)
        source = input_dir / 'anc/Gauge1Simp.in'
        source.write_text('2')
        scan._CODES.clear()
        def compile_terms(gauge):
            return [compile(source.read_text(), source.name, 'eval')] * 4
        with patch.object(scan.ev, 'DIR', input_dir), patch.object(scan, 'OUT', self.out), \
             patch.object(scan, 'split_terms', side_effect=compile_terms) as compiler, \
             patch.object(scan, 'identify', return_value=list(scan.MONOMIALS)):
            path, _ = scan.run_job((1, '1/2', 3, 40))
            before = Path(path).read_bytes()
            self.assertEqual(scan.run_job((1, '1/2', 3, 40))[1], 'cached')
            source.write_text('3')
            with self.assertRaises(provenance.CacheMismatch):
                scan.run_job((1, '1/2', 3, 40))
            self.assertEqual(Path(path).read_bytes(), before)
            scan.run_job((1, '1/2', 3, 40), recompute=True)
            self.assertEqual(compiler.call_count, 2)
            self.assertNotEqual(json.loads(before)['A'], json.loads(Path(path).read_text())['A'])
            # A worker receives the explicit output directory; it need not inherit OUT.
            fresh = self.root / 'worker-output'
            scan.OUT = self.root / 'wrong-worker-default'
            worker_path, _ = scan.run_job((1, '1/2', 3, 40), out_dir=fresh)
            self.assertEqual(Path(worker_path).parent, fresh)
            self.assertFalse(scan.OUT.exists())

    def test_failed_atomic_write_preserves_previous_result(self):
        path = self.root / 'result.json'
        path.write_text('{"previous": true}\n')
        before = path.read_bytes()
        with self.assertRaises(ValueError):
            provenance.atomic_write(path, {'result': float('nan')})
        self.assertEqual(path.read_bytes(), before)


if __name__ == '__main__':
    unittest.main()
