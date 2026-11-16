#!/usr/bin/env python3
"""Estimate the Letter's length using https://journals.aps.org/authors/length-guide.

Count body text (including data availability and footnotes) and the captions of tables and
figures; exclude the title, author block, abstract, acknowledgments and references. Displayed
equations contribute 16 words per row, or 32 within widetext. Tables contribute 13 + 6.5 per
line in one column, or 26 + 13 per line across two columns. Figures contribute 150/aspect + 20
in one column, or 300/(0.5 aspect) + 40 across two columns, with aspect the width of the
figure divided by its height.

End Matter (a section headed \\section*{End Matter} after the references) does not count
against the core limit (https://journals.aps.org/prl/authors); it is measured with the same
rules and reported separately, against its limit of about two pages, which at the core's
3750 words for about four pages is about 1875 word equivalents.

This is a source-based estimate, not the journal's official word count. Inline formulae
are treated as one word; LaTeX text tokenization and nested mathematical layouts are
approximate. Source table rows can wrap: inspect the final PDF and pass the rendered
line counts with --table-lines (one count per table, in source order, End Matter
included). The aspect ratio of a figure is read from the page size of its PDF file;
--figure-aspects gives the ratios by hand (one per figure, in source order).

    python3 count_words.py                  estimate using source table rows
    python3 count_words.py --table-lines 12  use the visually counted table lines

Both forms exit with an error if the estimate exceeds the 3750-word PRL limit.
"""
import argparse
import pathlib
import re
import sys
import zlib

HERE = pathlib.Path(__file__).resolve().parent
LIMIT = 3750
END_MATTER_LIMIT = 1875     # about two pages at the core's 3750 words for about four pages
EQ_WORDS, WIDE_EQ_WORDS = 16, 32
BRACED = r'((?:[^{}]|\{(?:[^{}]|\{[^{}]*\})*\})*)'


def expand_inputs(tex):
    def repl(m):
        return (HERE/(m.group(1) if m.group(1).endswith('.tex') else m.group(1) + '.tex')).read_text()
    return re.sub(r'\\input\{([^}]*)\}', repl, tex)


def strip_comments(tex):
    return re.sub(r'(?<!\\)%.*', '', tex)


def apply_macros(tex):
    """expand the \\newcommand macros without arguments (numbers.tex and the notation), so that the
    words they hold are counted as what appears on the page"""
    macros = dict(re.findall(r'\\newcommand\{\\([A-Za-z]+)\}\{((?:[^{}]|\{(?:[^{}]|\{[^{}]*\})*\})*)\}', tex))
    tex = re.sub(r'\\newcommand\{\\[A-Za-z]+\}(\[\d\])?\{((?:[^{}]|\{(?:[^{}]|\{[^{}]*\})*\})*)\}', '', tex)
    for _ in range(4):
        for k, v in macros.items():
            tex = re.sub(r'\\' + k + r'(?![A-Za-z])', lambda m: v, tex)     # a function replacement is used literally
    return tex


def words_in(text):
    """words of running text: inline formulae count as one word each, commands as nothing"""
    text = re.sub(r'\$\$(.*?)\$\$', ' W ', text, flags=re.S)
    text = re.sub(r'\\\((.*?)\\\)', ' W ', text, flags=re.S)
    text = re.sub(r'\$(.*?)\$', ' W ', text, flags=re.S)
    text = re.sub(r'\\(cite|ref|label|eqref|email|affiliation|bibliographystyle|bibliography|includegraphics|input)\*?\{[^}]*\}', ' ', text)
    text = re.sub(r'\\(begin|end)\{[^}]*\}', ' ', text)
    text = re.sub(r'\\[A-Za-z]+\*?', ' ', text)
    text = re.sub(r'[{}~]', ' ', text)
    return len([w for w in re.split(r'\s+', text) if re.search(r'[A-Za-z0-9]', w)])


def page_size(data):
    """(width, height) of the first page box found in the bytes of a PDF file; the page objects of
    a file written with object streams are compressed, so the streams are searched as well"""
    chunks = [data]
    for m in re.finditer(rb'stream\r?\n(.*?)\r?\nendstream', data, flags=re.S):
        try:
            chunks.append(zlib.decompress(m.group(1)))
        except zlib.error:
            pass
    number = rb'(-?\d+(?:\.\d+)?)'
    for chunk in chunks:
        box = re.search(rb'/MediaBox\s*\[\s*' + rb'\s+'.join([number]*4) + rb'\s*\]', chunk)
        if box:
            x0, y0, x1, y1 = (float(v) for v in box.groups())
            return x1 - x0, y1 - y0
    raise ValueError('no page size found in the PDF file')


def aspect_ratio(path):
    """width divided by height of a figure, from its PDF file"""
    path = pathlib.Path(path)
    if path.suffix.lower() != '.pdf':
        raise ValueError(f'{path.name}: only PDF figures are measured; use --figure-aspects')
    width, height = page_size(path.read_bytes())
    return width/height


def split(tex):
    """(core, End Matter): the body up to the references, and what follows \\section*{End Matter}"""
    body = tex.split('\\begin{document}', 1)[1]
    parts = re.split(r'\\bibliographystyle|\\bibliography|\\begin\{thebibliography\}|\\end\{document\}', body, maxsplit=1)
    rest = parts[1] if len(parts) > 1 else ''
    end = re.search(r'\\section\*\{End Matter\}(.*?)(?:\\end\{document\}|$)', rest, flags=re.S)
    return parts[0], (end.group(1) if end else None)


def count(tex, table_lines=None, figure_aspects=None):
    """Return ((item count, word equivalent) per category of the core, the same for the End Matter
    or None); excluded items have zero weight. table_lines and figure_aspects replace the source
    rows of the tables and the page sizes of the figure files, one value per table or figure in
    source order."""
    tex = apply_macros(strip_comments(expand_inputs(strip_comments(tex))))
    body, end_matter = split(tex)
    n_tables = len(re.findall(r'\\begin\{table\*?\}', body)) + (len(re.findall(r'\\begin\{table\*?\}', end_matter)) if end_matter else 0)
    if table_lines is not None:
        if len(table_lines) != n_tables or any(type(n) is not int or n <= 0 for n in table_lines):
            raise ValueError('--table-lines needs one positive integer per table, in source order')
    lines = list(table_lines) if table_lines is not None else None
    n_figures = len(re.findall(r'\\begin\{figure\*?\}', body)) + (len(re.findall(r'\\begin\{figure\*?\}', end_matter)) if end_matter else 0)
    if figure_aspects is not None:
        if len(figure_aspects) != n_figures or any(not a > 0 for a in figure_aspects):
            raise ValueError('--figure-aspects needs one positive number per figure, in source order')
    aspects = list(figure_aspects) if figure_aspects is not None else None
    body = re.sub(r'\\(?:author|email|affiliation|date|keywords|pacs|thanks)\{' + BRACED + r'\}', '', body)
    out = {}
    title = re.search(r'\\title(?:\[[^]]*\])?\{' + BRACED + r'\}', body)
    out['title (excluded)'] = (words_in(title.group(1)), 0) if title else (0, 0)
    body = re.sub(r'\\title(?:\[[^]]*\])?\{' + BRACED + r'\}', ' ', body)
    for environment in ('abstract', 'acknowledgments'):
        pattern = r'\\begin\{' + environment + r'\}(.*?)\\end\{' + environment + r'\}'
        matches = re.findall(pattern, body, flags=re.S)
        out[environment + ' (excluded)'] = (sum(words_in(text) for text in matches), 0)
        body = re.sub(pattern, ' ', body, flags=re.S)
    out.update(measure(body, lines, aspects))
    return out, (measure(end_matter, lines, aspects) if end_matter is not None else None)


def measure(body, lines, aspects=None):
    """displayed equations, tables (rows and captions), figures (area and captions) and running
    text of one part of the Letter; lines and aspects are the lists of rendered table line counts
    and of figure aspect ratios still to be used, or None"""
    out = {}
    # figures: the APS allowance for the area, and the caption as text
    figure_pattern = r'\\begin\{(figure\*?)\}(.*?)\\end\{\1\}'
    n_figures, figure_words, figure_captions = 0, 0, 0
    for m in re.finditer(figure_pattern, body, flags=re.S):
        graphics = re.findall(r'\\includegraphics(?:\[([^]]*)\])?\{([^}]*)\}', m.group(2))
        if not graphics:
            raise ValueError('each figure must include at least one graphics file')
        if aspects is not None:
            aspect = aspects.pop(0)
        elif len(graphics) == 1:
            aspect = aspect_ratio(HERE/graphics[0][1])
        else:
            # A two-panel figure occupies one float. Measure its displayed
            # bounding box, including both images, rather than counting only one.
            widths = []
            heights = []
            for options, name in graphics:
                match = re.fullmatch(r'\s*width\s*=\s*(0?\.\d+|1(?:\.0+)?)\\textwidth\s*', options)
                if match is None:
                    raise ValueError('multi-panel figures require a numeric width in units of textwidth')
                width = float(match.group(1))
                if width <= 0:
                    raise ValueError('figure widths must be positive')
                widths.append(width)
                heights.append(width / aspect_ratio(HERE/name))
            if sum(widths) > 1.000001:
                raise ValueError('multi-panel figure widths exceed textwidth')
            aspect = sum(widths) / max(heights)
        n_figures += 1
        figure_words += 300/(0.5*aspect) + 40 if m.group(1).endswith('*') else 150/aspect + 20
        cap = re.search(r'\\caption\{' + BRACED + r'\}', m.group(2), flags=re.S)
        figure_captions += words_in(cap.group(1)) if cap else 0
    body = re.sub(figure_pattern, ' ', body, flags=re.S)
    out['figures'] = (n_figures, round(figure_words, 1))
    out['figure captions'] = (figure_captions, figure_captions)
    # displayed equations: lines of equation/align/multline/gather bodies, with aligned/split inside
    eq_lines, wide_eq_lines = 0, 0
    equation_pattern = r'\\begin\{((?:equation|align|multline|gather)\*?)\}(.*?)\\end\{\1\}'
    for m in re.finditer(equation_pattern, body, flags=re.S):
        inner = re.sub(r'\\(begin|end)\{(aligned|split|array|cases)\}', '', m.group(2))
        inner = re.sub(r'\\label\{[^}]*\}', '', inner).strip()
        rows = sum(bool(row.strip()) for row in re.split(r'\\\\(?:\[[^]]*\])?', inner))
        before = body[:m.start()]
        wide = before.rfind('\\begin{widetext}') > before.rfind('\\end{widetext}')
        if wide:
            wide_eq_lines += rows
        else:
            eq_lines += rows
    body = re.sub(equation_pattern, ' ', body, flags=re.S)
    if re.search(r'\\\[|\$\$', body):
        raise ValueError('use equation/align/multline/gather environments for displayed math')
    out['equation rows (1 col)'] = (eq_lines, eq_lines*EQ_WORDS)
    out['equation rows (2 col)'] = (wide_eq_lines, wide_eq_lines*WIDE_EQ_WORDS)
    # tables: rows of tabular or tabular*, and the caption as text
    table_rows, table_words, captions = 0, 0, 0
    tables = list(re.finditer(r'\\begin\{(table\*?)\}(.*?)\\end\{\1\}', body, flags=re.S))
    for m in tables:
        star = m.group(1).endswith('*')
        # tabular* has a width argument before its column specification.
        # Keep nested @{} column declarations intact.
        tabular = re.search(
            r'\\begin\{(tabular\*?)\}'
            r'(?:\{' + BRACED + r'\})?'
            r'\{' + BRACED + r'\}'
            r'(.*?)\\end\{\1\}', m.group(2), flags=re.S)
        if tabular is None:
            raise ValueError('each table must contain a tabular or tabular* environment')
        tab = tabular.group(4)
        # Booktabs rules and spacing commands do not constitute data rows.
        rule_pattern = (r'\\(?:toprule|midrule|bottomrule|addlinespace)(?:\[[^]]*\])?'
                        r'|\\cmidrule(?:\([^)]*\))?\{[^}]*\}'
                        r'|\\(?:hline|cline)(?:\{[^}]*\})?')
        rows = [r for r in tab.split('\\\\')
                if re.search(r'[A-Za-z0-9]', re.sub(rule_pattern, '', r))]
        n = len(rows) if lines is None else lines.pop(0)
        table_rows += n
        table_words += 26 + 13*n if star else 13 + 6.5*n
        cap = re.search(r'\\caption\{' + BRACED + r'\}', m.group(2), flags=re.S)
        captions += words_in(cap.group(1)) if cap else 0
    body = re.sub(r'\\begin\{(table\*?)\}.*?\\end\{\1\}', ' ', body, flags=re.S)
    out['table source rows' if lines is None else 'table rendered lines'] = (table_rows, table_words)
    out['table captions'] = (captions, captions)
    out['text'] = (words_in(body),)*2
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--table-lines', metavar='N[,N,...]',
                        help='rendered line counts read from the final PDF, one per table in source order')
    parser.add_argument('--figure-aspects', metavar='A[,A,...]',
                        help='width divided by height, one per figure in source order, instead of the PDF page sizes')
    args = parser.parse_args()
    try:
        table_lines = [int(n) for n in args.table_lines.split(',')] if args.table_lines is not None else None
        aspects = [float(a) for a in args.figure_aspects.split(',')] if args.figure_aspects is not None else None
        out, end_matter = count((HERE/'letter.tex').read_text(), table_lines, aspects)
    except ValueError as error:
        parser.error(str(error))
    total = sum(v[1] for v in out.values())
    print('APS length estimate (not an official journal word count)')
    for k, (n, w) in out.items():
        print(f'  {k:28s} {n:5d}' + (f'  -> {w:7g} words' if w != n else ''))
    print(f'  {"estimated total":28s} {total:7g} words (limit {LIMIT})')
    if end_matter is not None:
        em = sum(v[1] for v in end_matter.values())
        print('End Matter (after the references; not counted against the core limit)')
        for k, (n, w) in end_matter.items():
            print(f'  {k:28s} {n:5d}' + (f'  -> {w:7g} words' if w != n else ''))
        print(f'  {"End Matter total":28s} {em:7g} word equivalents (about two pages, ~{END_MATTER_LIMIT})')
        if em > END_MATTER_LIMIT:
            print(f'the End Matter estimate is {em - END_MATTER_LIMIT:g} word equivalents over about two pages')
            sys.exit(1)
    if table_lines is None:
        print('Table source rows may wrap; inspect the final PDF and use --table-lines for rendered lines.')
    print('Text tokenization and equation rows remain estimates; check the final PDF against the APS guide.')
    if total > LIMIT:
        print(f'the estimate is {total - LIMIT:g} words over the PRL limit')
        sys.exit(1)


if __name__ == '__main__':
    main()
