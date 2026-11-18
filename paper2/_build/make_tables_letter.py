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
    t = ['\\begin{table*}[t]', '\\centering', '\\small', '\\setlength{\\tabcolsep}{5pt}',
         '\\begin{tabular}{cp{0.9\\textwidth}}', '\\toprule',
         '$L$ & MHV: coefficients of $2\\pi i\\,\\ell^k$, $k=2L-1,\\dots,1$\\\\', '\\midrule']
    t += [f'{L} & ${row(mhv, L, 2*L - 1, 1)}$\\\\' for L in range(1, loops + 1)]
    t += ['\\midrule', '$L$ & NMHV: coefficients of $\\pi i\\,\\frac{1+w}{1+w^*}\\,\\ell^k$, $k=2L-2,\\dots,0$\\\\', '\\midrule']
    t += [f'{L} & ${row(nmhv, L, 2*L - 2, 0)}$\\\\' for L in range(1, 5)]
    t += ['\\bottomrule', '\\end{tabular}',
          '\\caption{Planar $\\mathcal N=4$: the coefficients of the powers of $\\ell=\\ln|\\delta_*|$ in the $L$-loop '
          'MHV amplitude $\\mathcal E^{(L)}$, eqs.~(3.17)--(3.21) of ref.~\\cite{Dixon:2016epj}, and in the next-to-MHV '
          'amplitude $\\rho^{(L)}$, eqs.~(D.14)--(D.17) of its erratum, in $3\\to3$ kinematics, against '
          'eqs.~(\\ref{eq:kernel}), (\\ref{eq:Fm}) and~(\\ref{eq:n4exp}) with $\\nu=0$ and $\\nu=2$. The one-loop '
          'coefficients define $C_1$ and $K$; a coefficient followed by a constant in brackets fixes that constant, '
          'and every other coefficient is a prediction. All agree.}',
          '\\label{tab:n4letter}', '\\end{table*}']
    (HERE/'tables'/'table_n4_letter.tex').write_text('\n'.join(t) + '\n')
    return len(mhv) + len(nmhv)


if __name__ == '__main__':
    import subprocess
    import sys

    manuscript = HERE.parent / 'manuscript_prl'
    raise SystemExit(subprocess.call(
        [sys.executable, str(manuscript / 'make_tables.py')], cwd=manuscript
    ))
