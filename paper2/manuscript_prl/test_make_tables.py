"""Regression checks for nested manuscript inputs and their citations."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import make_tables as mt
import make_arxiv as ma

class NestedInputs(unittest.TestCase):
    def test_nested_labels_citations_and_cycles(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            (root/'supplement.tex').write_text(r'\input{section}\ref{sec:child}')
            (root/'section.tex').write_text(r'\label{sec:child}\cite{Nested}\input{supplement}')
            with patch.object(mt,'HERE',root):
                text=mt.source('supplement')
                self.assertIn('sec:child',mt.labels(text))
                self.assertEqual(mt.cited(root/'supplement.tex'),{'Nested'})
                self.assertEqual(text.count(r'\label{sec:child}'),1)

    def test_input_cannot_escape_manuscript_folder(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'letter.tex').write_text(r'\input{../outside}')
            with patch.object(mt,'HERE',root):
                with self.assertRaises(ValueError):mt.source('letter')

    def test_internal_review_files_not_in_ancillary_selection(self):
        source = ['revision/README.md', 'revision/apply_letter.py',
                  '.github/workflows/prl-revision.yml', 'README.md',
                  'checks/pinch_regions/param_space.py',
                  'paper2/manuscript_prl/letter.tex']
        selected, dirty = ma.select_ancillary(source, source,
                                               'paper2/manuscript_prl', {'letter.tex'})
        expected = ['README.md', 'checks/pinch_regions/param_space.py']
        self.assertEqual(selected, expected)
        self.assertEqual(dirty, expected)

    def test_no_internal_rewrite_markers_in_public_tex(self):
        self.assertNotIn('% Achievement-led PRL revision', (mt.HERE/'letter.tex').read_text())
        self.assertNotIn('% Audited pinch cross-checks and scope revision',
                         (mt.HERE/'supplement.tex').read_text())

    def test_new_pinch_source_is_archived(self):
        self.assertIn('pinch_checks.tex',{p.name for p in ma.source_list()})
        self.assertIn('Gardi:2024dissecting',mt.cited(mt.HERE/'supplement.tex'))

    def test_planar_parallel_panels_preserve_41_coefficients(self):
        # Exercise the final table renderer directly from independent structured
        # rows. The former test only exercised a now-retired TeX reparser.
        from types import SimpleNamespace
        renderer = mt.load(mt.FULL/'make_tables_letter.py', 'letter_direct_regression')
        rows = {'MHV': [], 'NMHV': []}
        expressions = []
        for mode, loops in (('MHV', range(1, 6)), ('NMHV', range(1, 5))):
            for L in loops:
                powers = (range(2 * L - 1, 0, -1) if mode == 'MHV'
                          else range(2 * L - 2, -1, -1))
                for k in powers:
                    value = rf'\mathrm{{{mode}}}_{{{L}{k}}}'
                    rows[mode].append({'L': L, 'k': k, 'DE': value, 'status': 'PASS'})
                    expressions.append(value)
        fake_mt = SimpleNamespace(
            load=lambda _: {'three_to_three': rows['MHV'],
                            'nmhv': rows['NMHV'], 'loops_compared': 5},
            sp=None, n4_tex=lambda expr: expr, N4_NAMES={})
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'tables').mkdir()
            with patch.object(renderer, 'HERE', root):
                count = renderer.letter_table(fake_mt)
            result = (root/'tables'/'table_n4_letter.tex').read_text()
        self.assertEqual(count, 41)
        self.assertEqual(len(expressions), 41)
        for expression in expressions:
            self.assertEqual(result.count('$' + expression + '$'), 1, expression)
        self.assertEqual(result.count(r'\begin{tabular}{@{}r@{'), 1)
        self.assertEqual(result.count(r'\end{tabular}'), 1)
        self.assertEqual(result.count(r'\addlinespace[3pt]'), 4)
        self.assertIn(
            '4 & 7 & $' + r'\mathrm{MHV}_{47}' +
            '$ & 4 & 6 & $' + r'\mathrm{NMHV}_{46}' + '$', result)
        last = result.split(
            '5 & 9 & $' + r'\mathrm{MHV}_{59}' + '$', 1)[1].split(
            r'\bottomrule', 1)[0]
        five_rows = [
            '5 & 9 & $' + r'\mathrm{MHV}_{59}' + '$' +
            last.splitlines()[0], *last.splitlines()[1:]]
        self.assertEqual(len(five_rows), 9)
        self.assertTrue(all(row.count('&') == 5 for row in five_rows))
        for k in range(9, 0, -1):
            self.assertIn('& ' + str(k) + ' & $' +
                          rf'\mathrm{{MHV}}_{{5{k}}}' + '$ & & & ' + r'\\', result)

if __name__=='__main__':unittest.main()
