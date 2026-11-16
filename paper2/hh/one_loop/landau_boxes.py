#!/usr/bin/env python3
# Leading Landau singularities of the one-loop topologies of gg -> YY (Y = A or H,
# two effective vertices), massless internal lines, external masses (0, 0, m^2, m^2).
#
# For a one-loop graph with massless propagators the leading Landau singularity sits
# where the modified Cayley determinant det Y vanishes, Y_ij = -(q_i - q_j)^2 with q_i
# the dual (region) momenta, i.e. Y_ij = -(sum of the external momenta between
# propagators i and j)^2. We list det Y for every ordering of the external legs.
#
# Expectation used in paper Sec. 2.5: only the box with the two massive legs at opposite
# corners (channels t and u, the crossed box of the DPS pinch) is singular at
# tu = m^4, i.e. at p_T = 0 (tu - m^4 = s p_T^2). The boxes with adjacent massive legs
# are singular only at s t = 0 or s u = 0; the triangles only where one of their
# external invariants vanishes, which for a massless gluon leg is the collinear
# configuration present at every p_T (it gives the infrared poles).
#
#   python3 landau_boxes.py               (sympy; one second)
import sympy as S

s, t, u, m2 = S.symbols('s t u m2')
# external momenta, all incoming: p1, p2 gluons; p3, p4 the massive bosons (incoming -p)
inv = {frozenset([1]): 0, frozenset([2]): 0, frozenset([3]): m2, frozenset([4]): m2,
       frozenset([1, 2]): s, frozenset([3, 4]): s,
       frozenset([1, 3]): t, frozenset([2, 4]): t,
       frozenset([1, 4]): u, frozenset([2, 3]): u}


def cayley(order):
    """det Y for the one-loop graph with external legs attached in cyclic `order`."""
    n = len(order)
    Y = S.zeros(n, n)
    for i in range(n):
        for j in range(i+1, n):
            legs = frozenset(order[i:j])             # legs between propagators i and j
            Y[i, j] = Y[j, i] = -inv[legs] if legs in inv else -inv[frozenset({1, 2, 3, 4})-legs]
    return S.factor(Y.det().subs(u, 2*m2-s-t))


print('boxes (cyclic order of the external legs; the three orderings up to reflection):')
for order in ((1, 2, 3, 4), (1, 2, 4, 3), (1, 3, 2, 4)):
    note = '   <- massive legs at opposite corners: zero at tu = m^4' if order == (1, 3, 2, 4) else ''
    print(f'  {order}: det Y = {cayley(order)}{note}')

# check the statements, with u = 2 m^2 - s - t
tu_m4 = S.expand((t*(2*m2-s-t)-m2**2))
assert S.simplify(cayley((1, 3, 2, 4))-tu_m4**2) == 0
assert S.simplify(cayley((1, 2, 3, 4))-(s*(2*m2-s-t))**2) == 0
assert S.simplify(cayley((1, 2, 4, 3))-(s*t)**2) == 0
print('  ok: (1,3,2,4) gives (tu - m^4)^2; (1,2,3,4) gives (s u)^2; (1,2,4,3) gives (s t)^2')

print('triangles (two adjacent legs joined at one vertex):')
for order in [((1,), (2,), (3, 4)), ((1,), (3,), (2, 4)), ((1,), (4,), (2, 3)), ((2,), (3,), (1, 4)),
              ((2,), (4,), (1, 3)), ((3,), (4,), (1, 2))]:
    vals = [inv[frozenset(o)] for o in order]
    Y = S.Matrix([[0, -vals[0], -vals[2]], [-vals[0], 0, -vals[1]], [-vals[2], -vals[1], 0]])
    d = S.factor(Y.det())
    why = ('vanishes for all kinematics (massless leg: collinear singularity, gives the 1/eps poles)'
           if d == 0 else 'zero only at s = 0')
    print(f'  external invariants {vals}: det Y = {d}, {why}')
print('Each triangle and bubble depends on one of s, t, u only. On the forward line at fixed')
print('beta, s stays fixed above 4 m^2 while t -> -m^2 (1-beta)/(1+beta) and')
print('u -> -m^2 (1+beta)/(1-beta) stay away from zero, so none of them is singular at p_T = 0;')
print('nor are the two boxes with adjacent massive legs. Only the box (1,3,2,4) depends on t and')
print('u together, and it is singular exactly at tu = m^4.')
