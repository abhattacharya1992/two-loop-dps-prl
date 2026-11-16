"""Factor a parametric integrand into C(eps) * prod_j f_j(v)^(n_j + m_j eps + l_j eta), with f_j irreducible
polynomials over Q, positive on the integration domain.

Bases raised to non-integer (eps-dependent) powers come from positive quantities (the combined transverse
denominators), so their factors are replaced by absolute values; for integer powers the signs are kept."""
import sympy as sp
from lc_engine import eps, D

eta = sp.Symbol('eta')


def split_factors(expr, vars_, test):
    """returns (const, {f: exponent}), expr = const * prod f^exponent on the domain containing `test`"""
    expr = expr.subs(D, 2 - 2*eps)
    state = {'const': sp.Integer(1), 'facs': {}}

    def has_vars(e):
        return any(e.has(v) for v in vars_)

    def poly_part(p, ex):
        """p a polynomial-like expression (no powers of variable expressions) raised to ex"""
        integer = sp.sympify(ex).is_integer
        num, den = sp.fraction(sp.together(p))
        for part, sgn in ((num, 1), (den, -1)):
            c, lst = sp.factor_list(sp.expand(part), *vars_)
            state['const'] *= (c if integer else sp.Abs(c))**(sgn*ex)
            for f, mult in lst:
                if not has_vars(f):
                    state['const'] *= (f if integer else sp.Abs(f))**(sgn*mult*ex)
                    continue
                val = sp.N(f.subs(test).subs({eps: sp.Rational(1, 100), eta: sp.Rational(1, 1000)}))
                if val < 0:
                    f = sp.expand(-f)
                    if integer:
                        state['const'] *= (-1)**(sgn*mult*ex)
                state['facs'][f] = sp.expand(state['facs'].get(f, 0) + sgn*mult*ex)

    def walk(e, ex):
        if not has_vars(e):
            state['const'] *= e**ex
            return
        if isinstance(e, sp.Mul):
            for t in e.args:
                walk(t, ex)
            return
        if isinstance(e, sp.Pow):
            b, p = e.args
            walk(b, sp.expand(p*ex))
            return
        if isinstance(e, sp.Add):
            # a sum: either a polynomial in the variables, or a sum containing powers -> combine first
            if any(isinstance(a, sp.Pow) and has_vars(a.args[0]) and not sp.sympify(a.args[1]).is_integer
                   for a in sp.preorder_traversal(e)):
                raise ValueError('sum with non-integer powers: ' + str(e)[:200])
            poly_part(e, ex)
            return
        poly_part(e, ex)

    walk(expr, sp.Integer(1))
    facs = {f: e for f, e in state['facs'].items() if e != 0}
    return sp.simplify(state['const']), facs


def _canon_factors(b, vars_, test):
    """b (rational in vars_, independent of eps and eta) -> (c, {f: k}) with b = c * prod f^k, every f an
    irreducible primitive polynomial over Z with f(test) > 0, and c free of vars_"""
    num, den = sp.fraction(sp.together(b))
    c = sp.Integer(1)
    out = {}
    for part, sgn in ((num, 1), (den, -1)):
        cc, lst = sp.factor_list(sp.expand(part), *vars_)
        c *= cc**sgn
        for f, mult in lst:
            if not any(f.has(v) for v in vars_):
                c *= f**(sgn*mult)
                continue
            if sp.N(f.subs(test)) < 0:
                f = sp.expand(-f)
                c *= sp.Integer(-1)**(sgn*mult)          # (v2 had the Python int -1: (-1)**-1 is a float)
            out[f] = out.get(f, 0) + sgn*mult
    return c, out


def _numeric_powers(e):
    """c^(m eps + l eta) for numbers c -> exp((m eps + l eta) sum_p k_p ln p), so that equal powers written with
    different bases (4^eps, 2^(2 eps)) are recognised as equal"""
    def is_np(a):
        return (isinstance(a, sp.Pow) and a.base.is_number and a.base.is_positive
                and (a.exp.has(eps) or a.exp.has(eta)))
    def plog(c):
        """ln c with every rational factored into primes"""
        out = 0
        for t in sp.Add.make_args(sp.expand_log(sp.log(c), force=True)):
            co, rest = t.as_coeff_Mul()
            if isinstance(rest, sp.log) and rest.args[0].is_Rational:
                r = sp.Rational(rest.args[0])
                for p, k in sp.factorint(r.p).items():
                    out += co*k*sp.log(p)
                for p, k in sp.factorint(r.q).items():
                    out -= co*k*sp.log(p)
            else:
                out += t
        return out
    e = e.replace(is_np, lambda a: sp.exp(sp.expand(a.exp*plog(a.base))))
    e = sp.powsimp(e, combine='exp', deep=True)
    return e.replace(lambda a: isinstance(a, sp.exp), lambda a: sp.exp(sp.expand(a.args[0])))


def normal_form(expr, vars_, test):
    """expr -> (const, {f: exponent}), expr = const * prod f^exponent on the domain containing `test`.

    Every base^(n + m eps + l eta) with a base that depends on the variables is first written as
    base^n * c^(m eps + l eta) * prod_f Ye_f^(k_f m) Yh_f^(k_f l), base = c prod_f f^k_f in irreducible factors
    (one symbol per irreducible factor, so that equal factors written differently in different terms get the
    same symbol); the Y-monomial must then be common to all terms and is taken out by factor().
    v1 assigned one symbol per base expression: two forms of one base (e.g. pi*P and P) got different symbols,
    the sum could not be factored, and the eps-powers of some terms were dropped."""
    # no expand_power_base: a positive base may be a product of negative factors, whose phases must not appear
    expr = sp.powdenest(expr.subs(D, 2 - 2*eps), force=True)
    ymap = {}
    Ys = []

    def ysym(f, which):
        key = (f, which)
        if key not in ymap:
            ymap[key] = sp.Symbol(f'Y{which}{len(Ys)}', positive=True)
            Ys.append((ymap[key], f, eps if which == 'e' else eta))
        return ymap[key]

    def repl(e):
        if isinstance(e, sp.Pow):
            b, p = e.args
            if any(b.has(v) for v in vars_) and (p.has(eps) or p.has(eta)):
                pe = sp.expand(p)
                n = pe.subs({eps: 0, eta: 0})
                m = pe.coeff(eps)
                l = pe.coeff(eta)
                assert sp.expand(pe - n - m*eps - l*eta) == 0, pe
                assert sp.sympify(n).is_integer, ('non-integer power', pe)
                assert not (b.has(eps) or b.has(eta)), ('eps-dependent base', b)
                assert not any(isinstance(a, sp.Pow) and not sp.sympify(a.args[1]).is_integer
                               for a in sp.preorder_traversal(b)), ('nested non-integer power', b)
                c, fl = _canon_factors(b, vars_, test)
                assert sp.N(c) > 0, ('base negative at the test point', b)
                out = b**n*c**(m*eps + l*eta)
                for f, k in fl.items():
                    if m != 0:
                        out *= ysym(f, 'e')**(k*m)
                    if l != 0:
                        out *= ysym(f, 'h')**(k*l)
                return out
            return repl(b)**p
        if e.args:
            return e.func(*[repl(a) for a in e.args])
        return e

    e2 = _numeric_powers(repl(expr))
    e2 = sp.factor(sp.together(e2))
    # e2 must now be a product: (rational in vars, with eps-dependent constants) * (monomial in the Y's)
    for a in sp.preorder_traversal(e2):
        if isinstance(a, sp.Add) and any(a.has(Y) for Y, _, _ in Ys):
            raise ValueError('Y-monomials not common to all terms: ' + str(a)[:300])
    const, facs = sp.Integer(1), {}
    c0, f0 = split_factors(e2.subs({Y: 1 for Y, _, _ in Ys}), vars_, test)
    const *= c0
    for f, ex in f0.items():
        facs[f] = sp.expand(facs.get(f, 0) + ex)
    powdict = sp.powsimp(e2, force=True).as_powers_dict()
    for Y, f, var in Ys:
        m = sum(ex for base, ex in powdict.items() if base == Y)
        if m == 0:
            continue
        facs[f] = sp.expand(facs.get(f, 0) + m*var)
    facs = {f: e for f, e in facs.items() if e != 0}
    return sp.simplify(const), facs
