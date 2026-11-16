"""Regression cases from the APS length guide; run with python3 -m unittest discover -s . -v."""
import unittest
import zlib
from unittest.mock import patch

from count_words import count, page_size


def document(body):
    return r'\begin{document}' + body + r'\end{document}'


def core(source, **kwargs):
    """the measures of the core (count() returns the core and the End Matter)"""
    return count(source, **kwargs)[0]


def total(result):
    return sum(words for _, words in result.values())


class LengthEstimate(unittest.TestCase):
    def test_excludes_front_matter_acknowledgments_and_references(self):
        result = core(document(r'''
            \title{Title with \emph{nested formatting}}
            \author{Author Name}\affiliation{Institute Name}
            \begin{abstract}An abstract.\end{abstract}
            Body text.\footnote{A body footnote.}
            \begin{acknowledgments}
                Thanks for help. \begin{equation}x = y\end{equation}
            \end{acknowledgments}
            Data availability.
            \bibliographystyle{apsrev4-2}
            These reference words must not count.
        '''))
        self.assertEqual(result['title (excluded)'], (4, 0))
        self.assertEqual(result['abstract (excluded)'], (2, 0))
        self.assertEqual(result['equation rows (1 col)'], (0, 0))
        self.assertEqual(total(result), 7)  # body (2), footnote (3), data availability (2)

    def test_equation_width_and_rows(self):
        result = core(document(r'''
            \begin{align}a &= b\\ c &= d\end{align}
            \begin{widetext}
                \begin{equation}e = f\end{equation}
            \end{widetext}
        '''))
        self.assertEqual(result['equation rows (1 col)'], (2, 32))
        self.assertEqual(result['equation rows (2 col)'], (1, 32))
        self.assertEqual(total(result), 64)

    def test_table_base_allowance_width_and_wrapped_lines(self):
        source = document(r'''
            \begin{table}
                \begin{tabular}{c}A\\B\\\end{tabular}
                \caption{One column.}
            \end{table}
            \begin{table*}
                \begin{tabular}{p{0.9\textwidth}}A wrapped cell\\\end{tabular}
                \caption{Two columns.}
            \end{table*}
        ''')
        estimated = core(source)
        rendered = core(source, table_lines=[2, 3])
        self.assertEqual(estimated['table source rows'], (3, 65))
        self.assertEqual(rendered['table rendered lines'], (5, 91))
        self.assertEqual(total(rendered), 95)  # 26 + 65 table allowance, plus four caption words
        with self.assertRaisesRegex(ValueError, 'one positive integer per table'):
            count(source, table_lines=[3])
        with self.assertRaisesRegex(ValueError, 'one positive integer per table'):
            count(source, table_lines=[0, 3])

    def test_tabular_star_with_nested_columns_and_booktabs_rules(self):
        source = document(r'''
            \begin{table*}
                \begin{tabular*}{\textwidth}{@{}r@{\hspace{5pt}}l@{\extracolsep{\fill}}r@{}}
                \toprule
                \multicolumn{3}{c}{MHV and NMHV}\\
                \cmidrule(lr){1-2}\cmidrule(l){3-3}
                $L$ & $k$ & coefficient\\
                \midrule
                1 & 1 & $a$\\
                \addlinespace[3pt]
                2 & 3 & $b$\\
                \bottomrule
                \end{tabular*}
                \caption{Logarithmic coefficients.}
            \end{table*}
        ''')
        estimated = core(source)
        self.assertEqual(estimated['table source rows'], (4, 26 + 13 * 4))
        self.assertEqual(estimated['table captions'], (2, 2))
        self.assertEqual(total(estimated), 80)
        rendered = core(source, table_lines=[6])
        self.assertEqual(rendered['table rendered lines'], (6, 26 + 13 * 6))
        self.assertEqual(total(rendered), 106)

    def test_rule_only_rows_are_not_counted(self):
        source = document(r'''
            \begin{table}
              \begin{tabular}{cc}
                A & B\\
                \cmidrule(lr){1-2}
                \addlinespace[2pt]
                \bottomrule
              \end{tabular}
              \caption{Small table.}
            \end{table}
        ''')
        self.assertEqual(core(source)['table source rows'], (1, 19.5))

    def test_end_matter_measured_apart_from_the_core(self):
        source = document(r'''
            Body text here.
            \begin{table}\begin{tabular}{c}A\\\end{tabular}\caption{Core table.}\end{table}
            \bibliographystyle{apsrev4-2}
            \bibliography{refs}
            \section*{End Matter}
            \textit{Appendix: Something.} Four more words here.
            \begin{table*}\begin{tabular}{c}A\\B\\\end{tabular}\caption{End table.}\end{table*}
        ''')
        body, end = count(source)
        self.assertEqual(total(body), 3 + 2 + 19.5)    # body text, caption, 13 + 6.5 per row
        self.assertEqual(total(end), 6 + 2 + 52)       # run-in heading and text, caption, 26 + 13 per row
        body, end = count(source, table_lines=[1, 3])  # in source order, End Matter included
        self.assertEqual(end['table rendered lines'], (3, 65))
        self.assertIsNone(count(document('Body text.'))[1])

    def test_figure_allowance_and_caption(self):
        source = document(r'''
            \begin{figure}\includegraphics[width=\columnwidth]{a.pdf}\caption{One column figure.}\end{figure}
            \begin{figure*}\includegraphics{b.pdf}\caption{Two columns.}\end{figure*}
        ''')
        result = core(source, figure_aspects=[1.5, 3.0])
        self.assertEqual(result['figures'], (2, 120 + 240))      # 150/1.5 + 20 and 300/(0.5*3) + 40
        self.assertEqual(result['figure captions'], (5, 5))
        self.assertEqual(total(result), 365)
        with self.assertRaisesRegex(ValueError, 'one positive number per figure'):
            count(source, figure_aspects=[1.5])

    def test_two_panel_figure_counts_both_displayed_images(self):
        source = document(r'''
            \begin{figure*}
            \includegraphics[width=0.63\textwidth]{left.pdf}
            \includegraphics[width=0.35\textwidth]{right.pdf}
            \caption{Two panel figure.}
            \end{figure*}
        ''')
        with patch('count_words.aspect_ratio', side_effect=[3.0, 1.0]):
            result = core(source)
        expected_aspect = 0.98 / max(0.63 / 3.0, 0.35 / 1.0)
        expected_words = round(300/(0.5*expected_aspect) + 40, 1)
        self.assertEqual(result['figures'], (1, expected_words))
        self.assertEqual(result['figure captions'], (3, 3))
        # Explicit measured aspect: one number per composite float.
        self.assertEqual(core(source, figure_aspects=[expected_aspect])['figures'],
                         (1, expected_words))

    def test_invalid_multi_panel_figure_fails(self):
        for options in ('width=0.75\\textwidth', 'scale=0.5'):
            source = document(r'\begin{figure*}\includegraphics[width=0.63\textwidth]{a.pdf}'
                              + '\\includegraphics[' + options + r']{b.pdf}\end{figure*}')
            with self.subTest(options=options), patch('count_words.aspect_ratio', return_value=2.0), self.assertRaises(ValueError):
                core(source)

    def test_page_size_read_from_plain_and_compressed_page_objects(self):
        plain = b'%PDF-1.4\n1 0 obj << /Type /Page /MediaBox [0 0 474.216 158.49] >> endobj'
        self.assertEqual(page_size(plain), (474.216, 158.49))
        packed = (b'%PDF-1.5\n5 0 obj << /Type /ObjStm >>\nstream\n'
                  + zlib.compress(b'<< /Type /Page /MediaBox [ 0 0 504 396 ] >>') + b'\nendstream')
        self.assertEqual(page_size(packed), (504.0, 396.0))
        with self.assertRaises(ValueError):
            page_size(b'%PDF-1.4 no page box')

    def test_unsupported_layouts_fail_instead_of_silently_undercounting(self):
        for source in (r'\begin{figure}\caption{A figure without a file.}\end{figure}', r'\[x = y\]'):
            with self.subTest(source=source), self.assertRaises(ValueError):
                count(document(source))


if __name__ == '__main__':
    unittest.main()
