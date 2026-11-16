"""Gate A1 for gluon splittings: the tree-level collinear function of g -> gg and the two components of the kernel.

The tree is t1 [S1 - S2/(1-x) - S3/x] with no T4, and the contraction of the two sectors (combine_gg) has the
x-dependence F_0 = 1/(x(1-x)) + 2x(1-x) for nu = 0 and F_2 = (x^4 + (1-x)^4)/(x(1-x)) for nu = 2, in the ratio
F_0/(2 F_2) (the nu = 0 part of (e1.b)(e2.b) is half the trace).

    python3 test_gg_tree.py
"""
import sympy as sp
import gg_diagrams as G
from gg_driver import tree_coefficients
import combine_gg as C

for xv in (sp.Rational(1, 10), sp.Rational(3, 10), sp.Rational(1, 2)):
    t = tree_coefficients(G.tree(), xv)
    tau = {k: sp.nsimplify(v/t['S1']) for k, v in t.items()}
    assert tau == {'S1': 1, 'S2': -1/(1 - xv), 'S3': -1/xv, 'T4': 0}, tau
    ty = tree_coefficients(G.tree(), 1 - xv)
    tauy = {k: complex(v/ty['S1']) for k, v in ty.items()}
    taux = {k: complex(v) for k, v in tau.items()}
    den = {nu: sum(taux[s]*tauy[u]*C.CT[nu][(s, u)] for s in C.STRUCT for u in C.STRUCT) for nu in (0, 2)}
    xf = float(xv)
    F0 = 1/(xf*(1 - xf)) + 2*xf*(1 - xf)
    F2 = (xf**4 + (1 - xf)**4)/(xf*(1 - xf))
    r = (den[0]/den[2]).real
    print(f'x = {xv}: tree {tau}, nu=0/nu=2 contraction {r:.12f}, F_0/(2 F_2) {F0/(2*F2):.12f}')
    assert abs(r - F0/(2*F2)) < 1e-12
print('passed')
