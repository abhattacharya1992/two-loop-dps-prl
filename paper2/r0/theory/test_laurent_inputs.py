"""Downstream kernels must reject partial or nonfinite Laurent records."""
from contextlib import contextmanager
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

from laurent_inputs import load_laurent
import combine_gg
import combine_n4
import r0_theory_quarks
with patch.object(sys, 'argv', ['test_laurent_inputs.py']):
    import combine_gqq

HERE = Path(__file__).resolve().parent


@contextmanager
def working_directory(path):
    before = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(before)


class LaurentConsumerTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.original = {}
        for kind in ('gg', 'gqq', 'n4', 'fdh'):
            name = f'{kind}_laurent_x1o2.json'
            record = json.loads((HERE / name).read_text())
            for key in ('provenance', 'complete', 'expected_diagrams', 'failed_diagrams'):
                record.pop(key, None)
            self.original[kind] = json.dumps(record)
            (self.root / name).write_text(self.original[kind])
        self.consumers = (
            ('gg', lambda: combine_gg.load('1o2', 'CA')),
            ('gqq', lambda: combine_gqq.load('1o2')),
            ('n4', lambda: combine_n4.load_n4('1o2')),
            ('gqq', lambda: r0_theory_quarks.load('1o2')),
        )

    def write(self, kind, record):
        path = self.root / f'{kind}_laurent_x1o2.json'
        path.write_text(json.dumps(record))
        return path

    def test_complete_legacy_data_remain_usable_and_unmodified(self):
        with working_directory(self.root):
            for kind, consumer in self.consumers:
                with self.subTest(kind=kind, consumer=consumer):
                    self.assertTrue(consumer())
            for kind in self.original:
                path = self.root / f'{kind}_laurent_x1o2.json'
                result = load_laurent(path, kind, expected_x='1/2')
                self.assertNotIn('provenance', result)
                self.assertEqual(path.read_text(), self.original[kind])

    def test_complete_tagged_records_preserve_their_existing_metadata(self):
        for kind in self.original:
            record = json.loads(self.original[kind])
            record.update(provenance={'stage': 'test fixture'}, complete=True,
                          expected_diagrams=sorted(record['diagrams']))
            path = self.write(kind, record)
            self.assertEqual(load_laurent(path, kind, expected_x='1/2'), record)

    def test_all_consumers_reject_a_missing_diagram_or_incomplete_flag(self):
        with working_directory(self.root):
            for kind, consumer in self.consumers:
                for missing in (True, False):
                    with self.subTest(kind=kind, missing=missing, consumer=consumer):
                        record = json.loads(self.original[kind])
                        if missing:
                            record['diagrams'].pop(next(iter(record['diagrams'])))
                        else:
                            record['complete'] = False
                        path = self.write(kind, record)
                        before = path.read_bytes()
                        with self.assertRaisesRegex(ValueError, 'invalid or incomplete'):
                            consumer()
                        self.assertEqual(path.read_bytes(), before)
                        path.write_text(self.original[kind])

    def test_all_consumers_reject_nonfinite_coefficients(self):
        with working_directory(self.root):
            for kind, consumer in self.consumers:
                with self.subTest(kind=kind, consumer=consumer):
                    record = json.loads(self.original[kind])
                    coefficients = next(coefficients for graph in record['diagrams'].values()
                                        for coefficients in graph['laurent'].values() if coefficients)
                    coefficients[next(iter(coefficients))][0] = float('nan')
                    path = self.write(kind, record)
                    with self.assertRaisesRegex(ValueError, 'finite numbers'):
                        consumer()
                    path.write_text(self.original[kind])

    def test_missing_colour_weights_wrong_x_and_missing_orders_are_rejected(self):
        for kind in self.original:
            for defect in ('colour', 'x', 'orders'):
                with self.subTest(kind=kind, defect=defect):
                    record = json.loads(self.original[kind])
                    first = next(iter(record['diagrams'].values()))
                    if defect == 'colour':
                        first['pref' if kind == 'gqq' else 'colour'] = {}
                    elif defect == 'x':
                        record['x'] = '3/10'
                    else:
                        coefficients = next(coefficients for graph in record['diagrams'].values()
                                            for coefficients in graph['laurent'].values() if coefficients)
                        coefficients.pop(next(iter(coefficients)))
                    path = self.write(kind, record)
                    with self.assertRaisesRegex(ValueError, 'invalid or incomplete'):
                        load_laurent(path, kind, expected_x='1/2')


if __name__ == '__main__':
    unittest.main()
