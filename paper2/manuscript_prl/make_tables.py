#!/usr/bin/env python3
# Writes numbers.tex, tables/ and sm.bib for the single-paper layout: the Letter (letter.tex) and its
# Supplemental Material (supplement.tex), which holds the derivations and the QCD amplitudes.
#
# The numbers and the tables are those of the full manuscript: this script runs the generator
# ../_build/make_tables.py with its output directed here, so that every macro and every
# table comes from the same result files, through the same assertions. It then
#   - writes the two tables of the End Matter: tables/table_n4_letter.tex, with the function of
#     ../_build/make_tables_letter.py, and tables/table_threeloop_letter.tex, the rows of the
#     three-loop table of the full manuscript with a caption for the Letter;
#   - keeps the tables that the two documents input and removes the others; their captions are
#     brought to the conventions of the two documents ("Sec.", "Eq.", "Ref."), and in the tables of
#     the Supplement a reference to an equation of the Letter becomes \Lref (see letter_labels.py);
#   - checks that no label is defined in both documents and that every \Lref points to a label
#     of the Letter;
#   - writes sm.bib and sm_nocite.tex, retaining Supplement-only citations
#     in the Letter bibliography while keeping the SM note concise.
#
#   python3 make_tables.py          (or `make tables`)
import importlib.util
from fractions import Fraction
import pathlib
import re

HERE = pathlib.Path(__file__).resolve().parent
FULL = HERE.parent/'_build'
LETTER_I = HERE.parent/'_build'
SM_TEXT = (
    'See Supplemental Material at [URL will be inserted by publisher] for the two-loop DPS derivation, '
    'matching calculations, numerical validation, and '
    'higher-order consistency checks.'
)



def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def source(doc):
    """Document plus non-table TeX inputs, without comments.

    Tables stay explicit because keep_used_tables applies their house style.
    Keeping each input command also preserves dependency discovery.
    """
    seen = set()
    def read(path):
        path = path.resolve()
        if not path.is_relative_to(HERE):
            raise ValueError(f'TeX input escapes manuscript directory: {path}')
        if path in seen:
            return ''
        seen.add(path)
        text = re.sub(r'(?<!\\)%.*', '', path.read_text())
        children = []
        for name in re.findall(r'\\input\{([^}]+)\}', text):
            if name.startswith('tables/'):
                continue
            name += '' if name.endswith('.tex') else '.tex'
            child = HERE / name
            # The initial table pass precedes these generated manuscript inputs.
            if name in {'numbers.tex', 'sm_nocite.tex', 'letter_labels.tex'} and not child.exists():
                continue
            children.append(read(child))
        return text + '\n' + '\n'.join(children)
    return read(HERE / f'{doc}.tex')

def labels(text):
    return set(re.findall(r'\\label\{([^}]*)\}', text))


def house_style(s):
    """the conventions of the two documents: Sec., Eq., Ref. and Table with capitals"""
    s = s.replace('(appendix~\\ref{', '(Sec.~\\ref{').replace('appendix~\\ref{', 'Sec.~\\ref{')
    s = s.replace('section~\\ref{', 'Sec.~\\ref{').replace('sections~\\ref{', 'Secs.~\\ref{')
    s = re.sub(r'\bsec\.~', 'Sec.~', s)
    s = re.sub(r'\beq(s?)\.~', r'Eq\1.~', s)
    s = re.sub(r'\bref(s?)\.~', r'Ref\1.~', s)
    return re.sub(r'\btable~\\ref', r'Table~\\ref', s)



def letter_threeloop_table(mt):
    """Render the final six checked coefficients, without parsing an old TeX table."""
    numbers = mt.load('three_loop')['numbers gg -> HH at mu = m']
    xs = sorted(numbers, key=Fraction)
    assert xs == ['5/8', '3/4', '5/6'], xs
    cols = ('N^2', 'N n_f', 'n_f^2')
    fmt = lambda v: f'{v[0]:.3f}' + (f'{v[1]:+.3f}\\,i' if abs(v[1]) > 5e-4 else '')
    t = [r'\begin{table}[t]', r'\centering', r'\footnotesize',
         r'\setlength{\tabcolsep}{7pt}', r'\renewcommand{\arraystretch}{1.35}',
         r'\begin{tabular}{lccc}', r'\toprule',
         r'$x$ & $N_c^2$ & $N_cn_f$ & $n_f^2$\\', r'\midrule']
    for name, heading in (('nu=0 G^3', r'$[\G^3]\,\M_1^{(3)}/C_1$'),
                          ('nu=2 G^2', r'$[\G^2]\,\M_2^{(3)}/K$')):
        if name != 'nu=0 G^3':
            t.append(r'\midrule')
        t.append(r'\multicolumn{4}{@{}l}{' + heading + r'}\\')
        for x in xs:
            t.append('$' + mt.frac_tex(Fraction(x)) + '$ & ' +
                     ' & '.join('$' + fmt(numbers[x][name][c]) + '$' for c in cols) + r'\\')
    t += [r'\bottomrule', r'\end{tabular}',
          r'\caption{Conditional three-loop DPS coefficients for $gg\to HH$ at $\mu=m_H$ in minimal subtraction; '
          r'the boundary assumptions are stated in the text.}',
          r'\label{tab:threeloop}', r'\end{table}']
    (HERE/'tables'/'table_threeloop_letter.tex').write_text('\n'.join(t) + '\n')


def keep_used_tables():
    docs = {d: source(d) for d in ('letter', 'supplement')}
    used = {d: set(re.findall(r'\\input\{tables/(table_[a-z0-9_]+)\.tex\}', s)) for d, s in docs.items()}
    assert not used['letter'] & used['supplement'], 'a table is input by both documents'
    present = {p.stem for p in (HERE/'tables').glob('table_*.tex')}
    assert used['letter'] | used['supplement'] <= present, (used['letter'] | used['supplement']) - present
    for name in sorted(present - used['letter'] - used['supplement']):
        (HERE/'tables'/f'{name}.tex').unlink()
    if 'tables/eq_threeloop.tex' not in docs['letter'] + docs['supplement']:
        (HERE/'tables'/'eq_threeloop.tex').unlink(missing_ok=True)
    text = {d: docs[d] + ''.join((HERE/'tables'/f'{n}.tex').read_text() for n in sorted(used[d])) for d in docs}
    own = {d: labels(text[d]) for d in docs}
    assert not own['letter'] & own['supplement'], sorted(own['letter'] & own['supplement'])
    for doc in docs:
        for name in sorted(used[doc]):
            p = HERE/'tables'/f'{name}.tex'
            s = house_style(p.read_text())
            s = re.sub(r'(\\begin\{table\*?\})(?:\[[^\]]*\])?', r'\1[!htb]', s)
            if doc == 'supplement':
                # the relations of the full manuscript's eq. (relations) are Eq. (4) of the Letter
                s = s.replace('the relations of Eq.~(\\ref{eq:relations})', 'the coefficients of Eq.~(\\Lref{eq:rg})')
                s = re.sub(r'\\ref\{([^}]*)\}', lambda m: ('\\Lref{' if m.group(1) in own['letter'] else '\\ref{') + m.group(1) + '}', s)
            p.write_text(s)
    # every reference resolves: within its own document, or to the Letter through \Lref
    supp = docs['supplement'] + ''.join((HERE/'tables'/f'{n}.tex').read_text() for n in sorted(used['supplement']))
    lett = docs['letter'] + ''.join((HERE/'tables'/f'{n}.tex').read_text() for n in sorted(used['letter']))
    bad = {k for k in re.findall(r'\\Lref\{([^}]*)\}', supp) if k not in own['letter']}
    bad |= {k for k in re.findall(r'\\ref\{([^}]*)\}', supp) if k not in own['supplement']}
    bad |= {k for k in re.findall(r'\\ref\{([^}]*)\}', lett) if k not in own['letter']}
    assert not bad, f'references without a label: {sorted(bad)}'
    return sorted(used['letter'] | used['supplement'])


def cited(path):
    keys = set()
    for m in re.findall(r'\\cite\{([^}]*)\}', source(path.stem)):
        keys |= {k.strip() for k in m.split(',')}
    return keys



def write_sm_bib():
    sm_only = sorted(
        cited(HERE/'supplement.tex')
        - cited(HERE/'letter.tex') - {'SM'}
    )

    (HERE/'sm.bib').write_text(
        '% generated by make_tables.py; do not edit\n'
        '@misc{SM,\n'
        f'    note = "{SM_TEXT}"\n'
        '}\n'
    )

    # Preserve Supplement-only references in the Letter bibliography.
    (HERE/'sm_nocite.tex').write_text(
        '% generated by make_tables.py; needed for clean arXiv compilation\n'
        + ('\\nocite{' + ','.join(sm_only) + '}\n'
           if sm_only else '')
    )
    return sm_only
    


if __name__ == '__main__':
    mt = load(FULL/'make_tables.py', 'make_tables_full')
    mt.HERE = HERE                   # numbers.tex and tables/ are written here; SOURCES stay as registered there
    mt.main()
    numbers = (HERE/'numbers.tex').read_text()
    for macro in ('rzFitAA', 'rzBlindAA', 'rzAAlead'):      # Table I of the Letter has a column for gg -> AA
        assert f'\\newcommand{{\\{macro}}}' in numbers, f'{macro} is not defined: the gg -> AA results are missing'
    one = load(LETTER_I/'make_tables_letter.py', 'make_tables_letter')
    one.HERE = HERE
    n = one.letter_table(mt)
    letter_threeloop_table(mt)
    tables = keep_used_tables()
    sm_only = write_sm_bib()
    print(f'numbers.tex, {len(tables)} tables ({", ".join(t[6:] for t in tables)}; the planar N = 4 table of the '
          f'End Matter holds {n} coefficients) and sm.bib (the Supplement adds {len(sm_only)} references) written')
