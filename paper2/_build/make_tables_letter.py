"""Planar comparison table used by the Letter.

The numerical coefficients are loaded and checked by the shared table generator.
This module only lays them out for the Letter.
"""
import pathlib

HERE = pathlib.Path(__file__).resolve().parent


def letter_table(mt):
    """tables/table_n4_letter.tex: the planar N = 4 comparison, one row per loop order."""
    res = mt.load('n4')
    sp = mt.sp

    def coeff(r):
        tex = mt.n4_tex(r['DE']).replace('\\frac', '\\tfrac')
        if r['status'].startswith('fixes'):
            tex += '\\,[' + mt.N4_NAMES[r['status'].split()[1]] + ']'
        else:
            assert r['status'] in ('PASS', 'normalisation'), r
        return tex

    def row(rows, L, kmax, kmin):
        rs = {r['k']: r for r in rows if r['L'] == L}
        assert set(rs) == set(range(kmin, kmax + 1)), (L, sorted(rs))
        return ',\\ '.join(coeff(rs[k]) for k in range(kmax, kmin - 1, -1))

    mhv, nmhv = res['three_to_three'], res['nmhv']
    loops = res['loops_compared']
    assert loops == 5 and {r['L'] for r in nmhv} == {1, 2, 3, 4}
    # Construct the final two-panel PRL v2 layout directly from checked rows,
    # instead of parsing and rewriting an intermediate TeX table.
    t = [r'\begin{table*}[t]', r'\centering', r'\small',
         r'\setlength{\tabcolsep}{5pt}',
         r'\begin{tabular}{@{}r@{\hspace{4pt}}r@{\hspace{6pt}}l@{\hspace{12pt}}r@{\hspace{4pt}}r@{\hspace{6pt}}l@{}}',
         r'\toprule',
         r'\multicolumn{3}{c}{\textbf{MHV}: $2\pi i\,\ell^k$} & \multicolumn{3}{c}{\textbf{NMHV}: $\pi i\,\frac{1+w}{1+w^*}\,\ell^k$}\\',
         r'\cmidrule(lr){1-3}\cmidrule(l){4-6}',
         r'$L$ & $k$ & coefficient & $L$ & $k$ & coefficient\\', r'\midrule']
    for L in range(1, 5):
        left, right = row(mhv, L, 2*L - 1, 1).split(r',\ '), row(nmhv, L, 2*L - 2, 0).split(r',\ ')
        assert len(left) == len(right) == 2*L - 1
        for i, (cl, cr) in enumerate(zip(left, right)):
            lh = f'{L if i == 0 else ""} & {2*L-1-i} & $' + cl + '$'
            rh = f'{L if i == 0 else ""} & {2*L-2-i} & $' + cr + '$'
            t.append(lh + ' & ' + rh + r'\\')
        if L < 4:
            t.append(r'\addlinespace[3pt]')
    t.append(r'\addlinespace[3pt]')
    last = row(mhv, 5, 9, 1).split(r',\ ')
    assert len(last) == 9
    for i, coeff_ in enumerate(last):
        t.append(f'{5 if i == 0 else ""} & {9-i} & $' + coeff_ + '$ & & & ' + r'\\')
    t += [r'\bottomrule', r'\end{tabular}',
          r'\caption{Planar $\mathcal N=4$: the coefficients of the powers of $\ell=\ln|\delta_*|$ in the $L$-loop '
          r'MHV amplitude $\mathcal E^{(L)}$, Eqs.~(3.17)--(3.21) of Ref.~\cite{Dixon:2016epj}, and in the next-to-MHV '
          r'amplitude $\rho^{(L)}$, Eqs.~(D.14)--(D.17) of its erratum, in $3\to3$ kinematics, against '
          r'Eqs.~(\ref{eq:kernel}), (\ref{eq:Fm}) and~(\ref{eq:n4exp}) with $\nu=0$ and $\nu=2$. The one-loop '
          r'coefficients define $C_1$ and $K$. The inputs $B^{(2)}$ and $B^{(3)}$ are fixed '
          r'from the published planar virtual anomalous dimension~\cite{CaronHuot:2019vjl}; '
          r'brackets mark the other constants fitted from MHV coefficients. '
          r'The remaining 18 MHV coefficients test the representation and agree.}',
          r'\label{tab:n4letter}', r'\end{table*}']
    (HERE/'tables'/'table_n4_letter.tex').write_text('\n'.join(t) + '\n')
    return len(mhv) + len(nmhv)


if __name__ == '__main__':
    import subprocess
    import sys

    manuscript = HERE.parent / 'manuscript_prl'
    raise SystemExit(subprocess.call(
        [sys.executable, str(manuscript / 'make_tables.py')], cwd=manuscript
    ))
