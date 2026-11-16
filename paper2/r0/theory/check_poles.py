"""Diagnostic: does the residue structure (which poles are picked, the i0 of the remaining denominators) and the
definiteness of the transverse denominators stay the same across each z-interval?  run_gqq.py decides both at
the midpoint only.

    python3 check_poles.py 3/10 Vqq
"""
import sys
import numpy as np
import sympy as sp
from gqq_diagrams import DIAGRAMS, x, z, u, w, k2, l2, lk

xval = sp.Rational(sys.argv[1]) if len(sys.argv) > 1 else sp.Rational(3, 10)
names = sys.argv[2].split(',') if len(sys.argv) > 2 else list(DIAGRAMS)
b = sorted({sp.Rational(-1), -xval, xval - 1, sp.Integer(0), xval, 1 - xval, sp.Integer(1)})
intervals = list(zip(b[:-1], b[1:]))


def step(dens, var, zt):
    """one residue step at z = zt: returns the list of (picked pole index, new dens) and the side pattern"""
    coef = [sp.cancel(sp.diff(e, var)) for e, s in dens]
    vd = [(e, s, c) for (e, s), c in zip(dens, coef) if c != 0]
    other = [(e, s) for (e, s), c in zip(dens, coef) if c == 0]
    sides = []
    for e, s, c in vd:
        loc = sp.N((s/c).subs(z, zt))
        sides.append('L' if loc > 0 else 'U')
    nlow = sides.count('L')
    side = 'L' if nlow <= len(sides) - nlow else 'U'
    out = []
    for i, (e, s, c) in enumerate(vd):
        if sides[i] != side:
            continue
        v0 = sp.cancel(var - e/c)
        newd = [(sp.cancel(sp.together(ej.subs(var, v0))), sp.cancel(sj - cj*s/c))
                for j, (ej, sj, cj) in enumerate(vd) if j != i]
        out.append((i, other + newd))
    return out, ''.join(sides)


def definite(d, zt):
    dd = sp.expand(sp.cancel(sp.together(d)).subs(z, zt))
    A, B, C = dd.coeff(l2), dd.coeff(lk), dd.coeff(k2)
    A, B, C = float(A), float(B), float(C)
    if A != 0:
        disc = B*B - 4*A*C
        return 'def' if disc < -1e-12 else ('semi' if abs(disc) <= 1e-12 else 'INDEF'), np.sign(A)
    return ('def' if B == 0 else 'LIN'), np.sign(C)


if __name__ != "__main__":
    names = []
for name in names:
    diag = DIAGRAMS[name]
    dens0 = [(sp.expand(d.subs(x, xval)), 1) for d in diag['dens']]
    print(name)
    for za, zb in intervals:
        pats = []
        for f in np.linspace(0.02, 0.98, 25):
            zt = za + (zb - za)*sp.Rational(int(f*1000), 1000)
            o1, s1 = step(dens0, u, zt)
            pat = [s1]
            defs = []
            for i, d1 in o1:
                o2, s2 = step(d1, w, zt)
                pat.append(f'{i}:{s2}')
                for j, d2 in o2:
                    for e, s in d2:
                        if sp.diff(e, l2) != 0 or sp.diff(e, lk) != 0:
                            defs.append(definite(e, zt))
            pats.append((' '.join(pat), tuple(defs)))
        uniq = sorted(set(p for p, d in pats))
        udef = sorted(set(d for p, d in pats))
        flag = '' if len(uniq) == 1 and len(udef) == 1 else '   <-- CHANGES'
        print(f'  ({float(za):+.2f},{float(zb):+.2f})  poles {uniq}  transverse {[[(a, int(s)) for a, s in d] for d in udef]}{flag}')
