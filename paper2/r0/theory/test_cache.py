"""Small regressions for region cache provenance, partial runs and failed checks.

Expensive symbolic diagrams/quadrature are mocked; all outputs use temporary
directories. Bundled numerical datasets are never rewritten by these tests.
"""
from contextlib import ExitStack, contextmanager
import json
import os
from pathlib import Path
import pickle
import runpy
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import warnings

import mpmath as mp
import sympy as sp
import cache_utils as cache
import assemble_lib
import archive_run
import gg_diagrams
import gg_driver

HERE = Path(__file__).resolve().parent


@contextmanager
def working_directory(path):
    before = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(before)


class RegionCacheTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.input_path = self.root / 'fdh_integrands_x1o2.pkl'
        self.output_path = self.root / 'fdh_laurent_x1o2.json'
        self.sector_path = self.root / 'sector.pkl'
        self.x = sp.Rational(1, 2)

    def integrands(self):
        data = cache.new_integrands(self.input_path, 'fdh', self.x)
        data['results'] = {name: [] for name in cache.diagram_names('fdh')}
        data['tree'] = {name: sp.Integer(1) for name in assemble_lib.STRUCT}
        cache.checkpoint(self.input_path, data, 'results')
        return cache.load_integrands(self.input_path, 'fdh', self.x)

    def assemble(self, data, **kwargs):
        return assemble_lib.assemble(data, gg_diagrams.DIAGRAMS, gg_diagrams.colour,
            self.output_path, self.sector_path, input_path=self.input_path, verbose=False, **kwargs)

    def source_copy(self):
        copy = self.root / 'sources'
        copy.mkdir(exist_ok=True)
        for name in set(cache.ASSEMBLY_SOURCES).union(*map(set, cache.INTEGRAND_SOURCES.values())):
            shutil.copy(HERE / name, copy / name)
        return copy

    def test_legacy_inputs_require_explicit_opt_in_and_are_not_relabelled(self):
        legacy = {'x': self.x, 'results': {name: [] for name in cache.diagram_names('fdh')},
                  'tree': {name: sp.Integer(1) for name in assemble_lib.STRUCT}}
        cache.atomic_write(self.input_path, legacy, binary=True)
        before = self.input_path.read_bytes()
        with self.assertRaisesRegex(cache.CacheMismatch, 'legacy'):
            cache.new_integrands(self.input_path, 'fdh', self.x)
        with self.assertRaisesRegex(cache.CacheMismatch, 'legacy'):
            cache.load_integrands(self.input_path, 'fdh', self.x)
        with warnings.catch_warnings(record=True):
            data = cache.load_integrands(self.input_path, 'fdh', self.x, allow_legacy=True)
        self.assertEqual(data['_input_provenance_status'], 'legacy-unverified')
        self.assertNotIn('provenance', data)
        self.assertEqual(self.input_path.read_bytes(), before)
        new = cache.new_integrands(self.input_path, 'fdh', self.x, recompute=True)
        self.assertFalse(new['complete'])
        self.assertEqual(new['results'], {})

    def test_partial_or_malformed_integrands_cannot_be_assembled(self):
        data = cache.new_integrands(self.input_path, 'fdh', self.x)
        data['tree'] = {name: sp.Integer(1) for name in assemble_lib.STRUCT}
        data['results'][data['expected_diagrams'][0]] = []
        cache.checkpoint(self.input_path, data, 'results')
        with self.assertRaisesRegex(cache.CacheMismatch, 'incomplete'):
            cache.load_integrands(self.input_path, 'fdh', self.x)
        data['results'][data['expected_diagrams'][0]] = None
        cache.atomic_write(self.input_path, data, binary=True)
        with self.assertRaisesRegex(cache.CacheMismatch, 'malformed'):
            cache.new_integrands(self.input_path, 'fdh', self.x)

    def test_integrand_source_change_rejects_old_checkpoint(self):
        self.integrands()
        before = self.input_path.read_bytes()
        copy = self.source_copy()
        with patch.object(cache, 'HERE', copy):
            cache.load_integrands(self.input_path, 'fdh', self.x)
            with (copy / 'vec_engine.py').open('a') as stream:
                stream.write('\n# source change for invalidation regression\n')
            with self.assertRaises(cache.CacheMismatch):
                cache.load_integrands(self.input_path, 'fdh', self.x)
        self.assertEqual(self.input_path.read_bytes(), before)

    def test_laurent_resume_binds_input_source_and_precision(self):
        data = self.integrands()
        fake = {name: {} for name in assemble_lib.STRUCT}
        with patch.object(assemble_lib, 'laurent_of_graph', return_value=fake) as numerical:
            self.assemble(data)
            count = numerical.call_count
            self.assertGreater(count, 0)
            original = self.output_path.read_bytes()
            self.assemble(data)
            self.assertEqual(numerical.call_count, count)
            for options in ({'neps': 16}, {'nq': {1: 96, 2: 72, 3: 40}}):
                with self.assertRaises(cache.CacheMismatch):
                    self.assemble(data, **options)
            with mp.workdps(mp.mp.dps + 5):
                with self.assertRaises(cache.CacheMismatch):
                    self.assemble(data)
            copy = self.source_copy()
            with (copy / 'secdec.py').open('a') as stream:
                stream.write('\n# changed quadrature implementation\n')
            with patch.object(cache, 'HERE', copy), self.assertRaises(cache.CacheMismatch):
                self.assemble(data)
            loaded = pickle.loads(self.input_path.read_bytes())
            loaded['note'] = 'same filename, different actual input bytes'
            cache.atomic_write(self.input_path, loaded, binary=True)
            with self.assertRaises(cache.CacheMismatch):
                self.assemble(data)
            self.assertEqual(self.output_path.read_bytes(), original)
            self.assemble(data, recompute=True)
            self.assertEqual(numerical.call_count, count * 2)
            runtime = json.loads(self.output_path.read_text())['provenance']['runtime']
            self.assertIn('scipy', runtime)

    def test_matching_metadata_does_not_certify_malformed_laurent(self):
        data = self.integrands()
        self.assemble(data)  # empty integrand pieces are exact zero; no quadrature
        record = json.loads(self.output_path.read_text())
        name = next(iter(record['diagrams']))
        del record['diagrams'][name]['laurent']['S1']
        self.output_path.write_text(json.dumps(record))
        before = self.output_path.read_bytes()
        with self.assertRaisesRegex(cache.CacheMismatch, 'malformed'):
            self.assemble(data)
        self.assertEqual(self.output_path.read_bytes(), before)

    def test_missing_or_nonfinite_colour_weights_cannot_be_resumed(self):
        data = self.integrands()
        self.assemble(data)
        original = json.loads(self.output_path.read_text())
        name = next(iter(original['diagrams']))
        for weights in ({}, {'CA': [1, 0]}, {'CA': [float('nan'), 0], 'nf': [0, 0]}):
            record = json.loads(json.dumps(original))
            record['diagrams'][name]['colour'] = weights
            self.output_path.write_text(json.dumps(record))
            before = self.output_path.read_bytes()
            with self.assertRaisesRegex(cache.CacheMismatch, 'malformed'):
                self.assemble(data)
            self.assertEqual(self.output_path.read_bytes(), before)

    def test_sector_cache_discards_legacy_and_changed_source(self):
        cache.atomic_write(self.sector_path, {'old-expression': 'old-sectors'}, binary=True)
        with warnings.catch_warnings(record=True) as caught:
            metadata, entries = cache.sector_cache(self.sector_path)
        self.assertTrue(caught)
        self.assertEqual(entries, {})
        cache.atomic_write(self.sector_path, {'provenance': metadata, 'entries': {'sentinel': 1}}, binary=True)
        self.assertEqual(cache.sector_cache(self.sector_path)[1], {'sentinel': 1})
        copy = self.source_copy()
        with (copy / 'factors.py').open('a') as stream:
            stream.write('\n# changed sector normal form\n')
        with patch.object(cache, 'HERE', copy), warnings.catch_warnings(record=True):
            self.assertEqual(cache.sector_cache(self.sector_path)[1], {})

    def test_failed_generation_exits_nonzero_and_checkpoint_can_resume(self):
        import driver, gqq_diagrams, n4_diagrams, vec_engine
        for script, module, generator, kind, flags in (
                ('run_gg.py', gg_diagrams, gg_driver, 'gg', []),
                ('run_gqq.py', gqq_diagrams, driver, 'gqq', []),
                ('run_n4.py', gg_diagrams, gg_driver, 'fdh', ['--fdh'])):
            with self.subTest(script=script), ExitStack() as stack:
                directory = self.root / kind
                directory.mkdir()
                stack.enter_context(working_directory(directory))
                stack.enter_context(patch.object(module, 'DIAGRAMS', {'Vfixture': {}}))
                stack.enter_context(patch.object(gg_driver, 'tree_coefficients',
                                             return_value={name: sp.Integer(1) for name in assemble_lib.STRUCT}))
                stack.enter_context(patch.object(vec_engine, 'METRIC_TRACE', vec_engine.METRIC_TRACE))
                stack.enter_context(patch.object(sys, 'argv', [script, '1/2'] + flags))
                with patch.object(generator, 'diagram_integrands', side_effect=ValueError('fixture failure')):
                    with self.assertRaises(SystemExit) as failed:
                        runpy.run_path(str(HERE / script), run_name='__main__')
                self.assertTrue(failed.exception.code)
                path = directory / f'{kind}_integrands_x1o2.pkl'
                checkpoint = pickle.loads(path.read_bytes())
                self.assertFalse(checkpoint['complete'])
                self.assertEqual(checkpoint['failed_diagrams'], ['Vfixture'])
                with patch.object(generator, 'diagram_integrands', return_value=[]):
                    runpy.run_path(str(HERE / script), run_name='__main__')
                checkpoint = pickle.loads(path.read_bytes())
                self.assertTrue(checkpoint['complete'])
                self.assertEqual(checkpoint['failed_diagrams'], [])

    def test_archiving_preserves_both_current_and_previous_integrands(self):
        current = self.root / 'gg_integrands_x1o2.pkl'
        stored = self.root / 'stored'
        stored.mkdir()
        current_bytes = pickle.dumps({'complete': True, 'generation': 'current'})
        previous_bytes = pickle.dumps({'complete': True, 'generation': 'previous'})
        current.write_bytes(current_bytes)
        (stored / current.name).write_bytes(previous_bytes)
        archive_run.archive('gg', self.root)
        self.assertFalse(current.exists())
        self.assertEqual((stored / current.name).read_bytes(), current_bytes)
        history = list((stored / 'history').glob(f'gg-*/{current.name}'))
        self.assertEqual(len(history), 1)
        self.assertEqual(history[0].read_bytes(), previous_bytes)

    def test_interrupted_or_empty_rerun_retains_full_comparison_baseline(self):
        stored = self.root / 'stored'
        stored.mkdir()
        name = 'gg_laurent_x1o2.json'
        baseline = b'{"complete": true, "diagrams": {"one": {}, "two": {}}}'
        (stored / name).write_bytes(baseline)
        (self.root / name).write_text('{"complete": false, "diagrams": {"one": {}}}')
        archive_run.archive('gg', self.root)
        self.assertEqual((stored / name).read_bytes(), baseline)
        self.assertEqual(len(list((stored / 'history').glob(f'gg-*/incomplete-{name}'))), 1)
        archive_run.archive('gg', self.root)
        self.assertEqual((stored / name).read_bytes(), baseline)
        run = subprocess.run([sys.executable, '-B', str(HERE / 'compare_stored.py'), 'gqq'],
                             cwd=self.root, capture_output=True, text=True)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn('comparison is unavailable', run.stderr)

    def test_fdh_checker_hashes_actual_inputs_and_fails_wrong_x_or_incomplete_set(self):
        for kind in ('gg', 'fdh', 'n4'):
            shutil.copy(HERE / f'{kind}_laurent_x1o2.json', self.root)
        command = [sys.executable, '-B', str(HERE / 'check_fdh.py'), '1/2']
        run = subprocess.run(command, cwd=self.root, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        checked = json.loads((self.root / 'fdh_check_x1o2.json').read_text())
        for name, digest in checked['provenance']['inputs'].items():
            self.assertEqual(digest, cache.provenance.sha256_file(self.root / name))
        path = self.root / 'fdh_laurent_x1o2.json'
        record = json.loads(path.read_text())
        original = json.dumps(record)
        record['x'] = '3/10'
        path.write_text(json.dumps(record))
        run = subprocess.run(command, cwd=self.root, capture_output=True, text=True)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn('recorded x', run.stderr)
        record['x'] = '1/2'
        record['diagrams'].pop(next(iter(record['diagrams'])))
        path.write_text(json.dumps(record))
        run = subprocess.run(command, cwd=self.root, capture_output=True, text=True)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn('incomplete', run.stderr)
        record = json.loads(original)
        name = next(name for name, diagram in record['diagrams'].items()
                    if complex(*diagram['colour']['CA']) != 0 and '0,0' in diagram['laurent']['S1'])
        record['diagrams'][name]['laurent']['S1']['0,0'][0] += 1
        path.write_text(json.dumps(record))
        run = subprocess.run(command, cwd=self.root, capture_output=True, text=True)
        self.assertNotEqual(run.returncode, 0)
        failed = json.loads((self.root / 'fdh_check_x1o2.json').read_text())
        self.assertFalse(failed['passed'])
        self.assertGreater(failed['largest deviation'], failed['absolute_tolerance'])


if __name__ == '__main__':
    unittest.main()
