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
        rows = [
            r'\begin{table*}[t]', r'\centering', r'\small',
            r'\begin{tabular}{cp{0.9\textwidth}}',
            r'\toprule', r'$L$ & MHV: coefficients\\', r'\midrule',
        ]
        expressions = []
        for mode, loops in (('MHV', range(1, 6)), ('NMHV', range(1, 5))):
            if mode == 'NMHV':
                rows += [r'\midrule', r'$L$ & NMHV: coefficients\\', r'\midrule']
            for L in loops:
                powers = (range(2 * L - 1, 0, -1) if mode == 'MHV'
                          else range(2 * L - 2, -1, -1))
                entries = [rf'\mathrm{{{mode}}}_{{{L}{k}}}' for k in powers]
                expressions.extend(entries)
                rows.append(f'{L} & $' + r',\ '.join(entries) + '$' + r'\\')
        rows += [r'\bottomrule', r'\end{tabular}',
                 r'\caption{ONE CAPTION}', r'\label{tab:n4letter}',
                 r'\end{table*}']
        result = mt.format_planar_rows('\n'.join(rows) + '\n')
        self.assertEqual(len(expressions), 41)
        for expression in expressions:
            self.assertEqual(result.count('$' + expression + '$'), 1, expression)
        self.assertEqual(result.count(r'\begin{tabular}{@{}r@{'), 1)
        self.assertEqual(result.count(r'\end{tabular}'), 1)
        self.assertEqual(result.count(r'\caption{ONE CAPTION}'), 1)
        self.assertEqual(result.count(r'\addlinespace[3pt]'), 3)
        self.assertEqual(result.count(r'\multicolumn{4}{@{}l@{}}'), 0)
        self.assertIn(r'\multicolumn{6}{@{}l}{\textbf{MHV}, $L=5$', result)
        self.assertNotIn(r'\begin{tabular*}{\textwidth}',result)
        self.assertIn(
            '4 & 7 & $' + r'\mathrm{MHV}_{47}' +
            '$ & 4 & 6 & $' + r'\mathrm{NMHV}_{46}' + '$',
            result)
        self.assertIn(
            '5 & 9 & $' + r'\mathrm{MHV}_{59}' +
            '$ & 5 & 4 & $' + r'\mathrm{MHV}_{54}' + '$', result)
        self.assertIn(
            ' & 6 & $' + r'\mathrm{MHV}_{56}' +
            '$ &  & 1 & $' + r'\mathrm{MHV}_{51}' + '$', result)
        self.assertIn(
            ' & 5 & $' + r'\mathrm{MHV}_{55}' +
            '$ & & & ' + r'\\', result)
        five_loop = result.split(r'\cmidrule(l){1-6}', 1)[1].split(r'\bottomrule', 1)[0]
        five_rows = [line for line in five_loop.splitlines() if '&' in line]
        self.assertEqual(len(five_rows), 5)
        self.assertTrue(all(row.count('&') == 5 for row in five_rows))
        self.assertNotIn('$ &  & &',result)

if __name__=='__main__':unittest.main()
