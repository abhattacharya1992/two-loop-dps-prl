"""Turn a diagram (Dirac strings, light-cone denominators) into parametric integrands over the loop momentum
fraction z and Feynman parameters, at a fixed numerical momentum fraction x.

For each diagram and each z-interval with support, returns the integrands of the two helicity structures
'ek' and 'ke' (coefficients of (eps.Gamma)(Gamma.k) and (Gamma.k)(eps.Gamma)), normalised as the ratio to the
tree, divided by a_s (mu^2 e^gamma/k2)^eps, and without the eikonal factor nu^eta |z|^-eta/(z + i0):
    ratio = a_s (mu^2 e^gamma/k2)^eps  int dz [eik]  int dF  I(z, F),
where int dF runs over the Feynman-parameter simplex.  The (4 pi)^(-1+eps) of the transverse integral and
the 1/(2 pi) of the z measure are absorbed in this normalisation (see gqq_diagrams.py for the tree).
"""
import sympy as sp
import lc_engine as E
from gqq_diagrams import x, z, u, w, k2, l2, lk


def ncoeff_terms(strings):
    tot = {}
    for c, s in strings:
        words = E.dirac_ncoeff(s)
        for word, cw in words.items():
            if isinstance(c, str):
                sign = -1 if c.startswith('-') else 1
                for pre in (('e', 'l'), ('l', 'e')):
                    key = pre + word
                    tot[key] = sp.expand(tot.get(key, 0) + sign*cw)
            else:
                tot[word] = sp.expand(tot.get(word, 0) + c*cw)
    return {k: v for k, v in tot.items() if v != 0}


def residues(terms, dens, at, order=(u, w)):
    """double residue integral; returns {tuple(dens): {word: numerator}}"""
    out = {}
    d0 = [(d, 1) for d in dens]
    for n1, d1 in E.lc_integrate_multi(terms, d0, order[0], at):
        for n2, d2 in E.lc_integrate_multi(n1, d1, order[1], at):
            key = tuple(sp.factor(e) for e, s in d2)
            out.setdefault(key, {})
            for word, c in n2.items():
                out[key][word] = sp.cancel(out[key].get(word, 0) + c)
    return out


def diagram_integrands(diag, xval, intervals, order=(u, w)):
    """returns list of (z-interval, params, {'ek': I, 'ke': I}) with x = xval"""
    terms = ncoeff_terms(diag['strings'])
    terms = {k: sp.expand(v.subs(x, xval)) for k, v in terms.items()}
    dens = [sp.expand(d.subs(x, xval)) for d in diag['dens']]
    res = []
    for (za, zb) in intervals:
        zt = (sp.Rational(za) + sp.Rational(zb))/2
        at = {z: zt, k2: 1, l2: sp.Rational(7, 3), lk: sp.Rational(1, 5)}
        rr = residues(terms, dens, at, order)
        n = 2 - 2*E.eps
        strip = (4*sp.pi)**(n/2)
        tree = {'ek': -xval/2, 'ke': (1 - xval)/2}
        pieces = []
        for dl, tw in rr.items():
            # separate l-independent denominators
            dep = [d for d in dl if sp.diff(d, l2) != 0 or sp.diff(d, lk) != 0]
            ind = [d for d in dl if d not in dep]
            fac = 1/sp.Mul(*ind) if ind else 1
            tw2 = {kw: sp.expand(v*fac) for kw, v in tw.items()}
            if not dep:
                raise ValueError('no transverse denominator')
            r, prm = E.transverse_integrate(tw2, dep, l2, lk, k2, test=at)
            # normalisation: ratio/(a_s (mu^2 e^g/k2)^eps) = p * int dz eik * T/tree, T without (4 pi)^(-n/2)
            out = {st: (r[st]*strip/tree[st] if r[st] != 0 else 0) for st in ('ek', 'ke')}
            pieces.append((prm, out))
        res.append(((za, zb), pieces))
    return res
