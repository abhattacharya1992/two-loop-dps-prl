"""Light-cone tools for the one-loop collinear functions of the DPS kernel (gate C, theory).

Conventions: v^+ = v^0 + v^3, v^- = v^0 - v^3, v.w = (v^+ w^- + v^- w^+)/2 - vperp.wperp (Euclidean
transverse vectors), d^4v = (1/2) dv^+ dv^- d^2vperp.  n = (1,0,0,1), nbar = (1,0,0,-1), so that
nslash = gamma.n carries v^- ... more precisely vslash = (v^+/2) nslash + (v^-/2) nbarslash + vslash_perp,
with vslash_perp = -(v_perp . Gamma) where Gamma^i are the Euclidean transverse gammas,
{Gamma^i, Gamma^j} = -2 delta^{ij}, sum_i Gamma^i Gamma^i = -D, D = d - 2.

1. Dirac strings.  A string is a list of factors; dirac_ncoeff(string) returns the coefficient of nslash
   (nslash to the left) as a dict {perp word: coefficient}, where a perp word is a tuple of tokens:
   'e' (eps.Gamma), a vector name ('k', 'l', ...) meaning (vector . Gamma), or ('i', label) for one of a pair
   of contracted Euclidean indices.  This is all that the projections onto the good components keep:
   (nslash nbarslash/4) X (nbarslash nslash/4) = nslash X_1, with X_1 the nslash-coefficient of X.
2. Perp words.  perp_contract removes contracted index pairs; perp_canonical reduces a word with one 'e'
   and vectors 'k' to (-k2)^m (Gamma.k)(eps.Gamma) or (eps.Gamma)(Gamma.k).
3. Light-cone energy integrals.  lc_integrate does int dv/(2 pi) over the minus components by residues, with
   the i0 of every denominator tracked through the substitutions.
"""
import itertools
import sympy as sp

# ----------------------------------------------------------------------------- 1. Dirac strings
# A factor is one of:
#   ('slash', plus, minus, {vec: coeff})   the slash of a momentum with these light-cone components
#   ('eps',)                               epsslash for a transverse physical polarisation (= -(eps.Gamma))
#   ('nbar',)                              nbarslash (eikonal vertex of an n-collinear Wilson line)
#   ('n',)                                 nslash
#   ('gam', label)                         gamma^alpha with a Lorentz index contracted (Feynman gauge
#                                          -g_{alpha beta}) with the other factor of the same label
def _expand_factor(f):
    """basis expansion of one factor: list of (coeff, element), element in 'N', 'Nb', ('p', token)"""
    if f[0] == 'slash':
        _, plus, minus, perp = f
        out = []
        if plus != 0:
            out.append((sp.Rational(1, 2)*plus, 'N'))
        if minus != 0:
            out.append((sp.Rational(1, 2)*minus, 'Nb'))
        for vec, c in perp.items():
            if c != 0:
                out.append((-c, ('p', vec)))
        return out
    if f[0] == 'eps':
        return [(-1, ('p', 'e'))]
    if f[0] == 'nbar':
        return [(1, 'Nb')]
    if f[0] == 'n':
        return [(1, 'N')]
    raise ValueError(f)


def _split_contractions(string):
    """gamma^a X gamma_a = (1/2)(n X nbar + nbar X n) + gperp: returns list of (coeff, string without 'gam')"""
    labels = sorted({f[1] for f in string if f[0] == 'gam'})
    results = [(sp.Integer(1), list(string))]
    for lab in labels:
        new = []
        for c, s in results:
            pos = [i for i, f in enumerate(s) if f[0] == 'gam' and f[1] == lab]
            assert len(pos) == 2
            p, q = pos
            for coeff, a, b in ((sp.Rational(1, 2), ('n',), ('nbar',)), (sp.Rational(1, 2), ('nbar',), ('n',))):
                t = list(s); t[p] = a; t[q] = b
                new.append((c*coeff, t))
            # perp part: gamma_perp^a X gamma_perp a = -sum_i Gamma^i X Gamma^i
            t = list(s); t[p] = ('pidx', lab); t[q] = ('pidx', lab)
            new.append((-c, t))
        results = new
    return results


_LC_REDUCE = {}


def _reduce_lc(word):
    """reduce a word in N, Nb (N N = Nb Nb = 0, N Nb N = 4N, Nb N Nb = 4Nb): returns (coeff, word) or (0, None)"""
    w = list(word)
    c = 1
    changed = True
    while changed:
        changed = False
        for i in range(len(w) - 1):
            if w[i] == w[i + 1]:
                return 0, None
        for i in range(len(w) - 2):
            if w[i] == w[i + 2] and w[i] != w[i + 1]:
                # X Y X = 4 X  (X Y X = X (4 - X Y) = 4X - XXY)
                w = w[:i] + [w[i]] + w[i + 3:]
                c *= 4
                changed = True
                break
    return c, tuple(w)


def dirac_ncoeff(string):
    """coefficient of nslash (to the left) of a Dirac string, as {perp word: coeff}"""
    out = {}
    for c0, s in _split_contractions(string):
        expanded = []
        for f in s:
            if f[0] == 'pidx':
                expanded.append([(1, ('p', ('i', f[1])))])
            else:
                expanded.append(_expand_factor(f))
        for combo in itertools.product(*expanded):
            coeff = c0
            lc, perp, sign = [], [], 1
            for c, el in combo:
                coeff = coeff*c
                if el in ('N', 'Nb'):
                    # move to the left past all perp elements collected so far
                    if len(perp) % 2:
                        sign = -sign
                    lc.append(el)
                else:
                    perp.append(el[1])
            if coeff == 0:
                continue
            k, word = _reduce_lc(lc)
            if word != ('N',):
                continue
            key = tuple(perp)
            out[key] = sp.expand(out.get(key, 0) + sign*k*coeff)
    return {k: v for k, v in out.items() if v != 0}


# ----------------------------------------------------------------------------- 2. perp words
D = sp.Symbol('D')          # D = d - 2


def _dot(a, b, dots):
    """Euclidean dot product of two tokens (vectors or index tokens)"""
    if isinstance(a, tuple) and isinstance(b, tuple):
        return 1 if a == b else 0          # delta^{ij}: only used when equal labels are contracted later
    if isinstance(a, tuple) or isinstance(b, tuple):
        raise ValueError('dot of an index with a vector should not occur here')
    key = tuple(sorted((a, b)))
    return dots[key]


def perp_contract(word, dots):
    """remove all contracted index pairs: returns {word without indices: coeff}.
    sum_i Gamma^i (Y X') Gamma^i = -Y (sum_i Gamma^i X' Gamma^i) - 2 X' Y,  sum_i Gamma^i Gamma^i = -D."""
    idx = [t for t in word if isinstance(t, tuple)]
    if not idx:
        return {tuple(word): sp.Integer(1)}
    lab = idx[0]
    p = word.index(lab)
    q = p + 1 + word[p + 1:].index(lab)
    head, mid, tail = list(word[:p]), list(word[p + 1:q]), list(word[q + 1:])

    def inner(X):
        """sum_i Gamma^i X Gamma^i as {word: coeff}"""
        if not X:
            return {(): -D}
        Y, Xp = X[0], X[1:]
        out = {}
        for w2, c in inner(Xp).items():
            key = (Y,) + w2
            out[key] = out.get(key, 0) - c
        key = tuple(Xp) + (Y,)
        out[key] = out.get(key, 0) - 2
        return out

    result = {}
    for w2, c in inner(mid).items():
        for w3, c3 in perp_contract(tuple(head) + w2 + tuple(tail), dots).items():
            result[w3] = sp.expand(result.get(w3, 0) + c*c3)
    return {k: v for k, v in result.items() if v != 0}


def perp_canonical(word, k2):
    """a word with one 'e' and only 'k' otherwise -> (coeff, 'ek' or 'ke'): (Gamma.k)^2 = -k2"""
    assert word.count('e') == 1 and all(t in ('e', 'k') for t in word), word
    i = word.index('e')
    a, b = i, len(word) - i - 1
    c = (-k2)**(a//2 + b//2)
    if a % 2 == 1 and b % 2 == 0:
        return c, 'ke'
    if a % 2 == 0 and b % 2 == 1:
        return c, 'ek'
    raise ValueError('even number of k: ' + str(word))


# ----------------------------------------------------------------------------- 3. light-cone integrals
def lc_integrate(num, dens, var, signs_at):
    """int dvar/(2 pi) num / prod (den + i0 s) by residues.
    dens: list of (expr, s) with expr linear in var, s the coefficient of i0 (a sympy number/expression);
    signs_at: dict of numerical values used to decide the half-plane of each pole.
    Returns a list of (num', dens') terms (one per residue), with the i0 coefficients of the remaining
    denominators updated. Requires num/prod to fall off at least like 1/var^2."""
    dens = [(sp.cancel(sp.together(e)), s) for e, s in dens]
    vd = [(e, s) for e, s in dens if sp.simplify(sp.diff(e, var)) != 0]
    other = [(e, s) for e, s in dens if sp.simplify(sp.diff(e, var)) == 0]
    deg = sp.Poly(sp.expand(num), var).degree() if sp.expand(num) != 0 else -10
    if deg - len(vd) > -2:
        raise ValueError(f'integrand falls off too slowly in {var}: degree {deg}, {len(vd)} denominators')
    if not vd:
        return []
    poles = []
    for i, (e, s) in enumerate(vd):
        b = sp.simplify(sp.diff(e, var))
        loc = sp.simplify((s/b).subs(signs_at))
        poles.append((i, b, 'lower' if loc > 0 else 'upper'))
    nlow = sum(1 for p in poles if p[2] == 'lower')
    side = 'lower' if nlow <= len(poles) - nlow else 'upper'
    pref = -sp.I if side == 'lower' else sp.I         # (1/2pi) * (-+ 2 pi i)
    out = []
    for i, b, where in poles:
        if where != side:
            continue
        e, s = vd[i]
        v0 = sp.solve(sp.Eq(e, 0), var)[0]
        newd = []
        for j, (ej, sj) in enumerate(vd):
            if j == i:
                continue
            bj = sp.simplify(sp.diff(ej, var))
            newd.append((sp.cancel(sp.together(ej.subs(var, v0))), sp.simplify(sj - bj*s/b)))
        out.append((sp.cancel(sp.together(pref*num.subs(var, v0)/b)), other + newd))
    return out


# ----------------------------------------------------------------------------- 4. transverse integrals
eps = sp.Symbol('eps')


def _pairings(items):
    """all perfect matchings of a list of items"""
    if not items:
        yield []
        return
    a = items[0]
    for i in range(1, len(items)):
        rest = items[1:i] + items[i + 1:]
        for m in _pairings(rest):
            yield [(a, items[i])] + m


def transverse_integrate(terms, dens, l2, lk, k2, lvec='l', kvec='k', test=None):
    """int d^n l/(2 pi)^n  sum_w c_w(l2, lk, k2) w / prod_j den_j,  n = 2 - 2 eps.
    terms: {word: coeff}, words contain 'e', kvec and lvec tokens; dens: list of expressions linear in
    (l2, lk, k2) (all of one sign at the test point). Feynman parameters a1..a_{m-1} (a_m = 1 - sum).
    Returns ({'ek': integrand, 'ke': integrand}, params), the integrand to be integrated over the simplex
    sum a_i <= 1 with measure da_1...da_{m-1}; the result carries the factor (k2)^(-1-eps)."""
    m = len(dens)
    params = sp.symbols(f'a1:{m}')
    alphas = list(params) + [1 - sum(params)]
    # make every denominator positive definite (flip its sign), and check definiteness at the test point
    signs = []
    for d in dens:
        dd = sp.expand(sp.cancel(sp.together(d)))
        A_ = dd.coeff(l2); B_ = dd.coeff(lk); C_ = dd.coeff(k2)
        tv = {kk: vv for kk, vv in (test or {}).items() if kk not in (l2, lk, k2)}
        An, Bn, Cn = (sp.N(c.subs(tv)) for c in (A_, B_, C_))
        if An != 0:
            assert Bn**2 - 4*An*Cn < 1e-12*abs(An*Cn) + 1e-30 and An*Cn >= 0, ('not definite', dd)
            signs.append(1 if An > 0 else -1)
        else:
            assert Bn == 0, ('linear in lk', dd)
            signs.append(1 if Cn > 0 else -1)
    sgn = 1
    dens = [sg*sp.expand(sp.cancel(sp.together(d))) for sg, d in zip(signs, dens)]
    sgn_total = 1
    for sg in signs:
        sgn_total *= sg
    comb = sp.expand(sum(a*d for a, d in zip(alphas, dens)))
    A = sp.expand(comb.coeff(l2))
    B = sp.expand(comb.coeff(lk))
    C = sp.expand(comb.coeff(k2))
    assert sp.expand(comb - A*l2 - B*lk - C*k2) == 0
    beta = sp.cancel(-B/(2*A))             # l = lp + beta k
    M = sp.cancel(C - B**2/(4*A))          # combined = A lp^2 + M k2
    n = 2 - 2*eps
    lp2, lpk = sp.symbols('lp2 lpk')
    out = {'ek': 0, 'ke': 0}
    for word, c in terms.items():
        # substitute l -> lp + beta k in the scalar coefficient
        cs = sp.expand(c.subs({l2: lp2 + 2*beta*lpk + beta**2*k2, lk: lpk + beta*k2}, simultaneous=True))
        # expand the word: each lvec token -> 'lp' + beta kvec
        lpos = [i for i, t in enumerate(word) if t == lvec]
        for choice in itertools.product((0, 1), repeat=len(lpos)):
            w = list(word)
            fac = 1
            for pos, ch in zip(lpos, choice):
                if ch:
                    w[pos] = kvec; fac *= beta
                else:
                    w[pos] = 'lp'
            # scalar polynomial in lp2, lpk
            poly = sp.Poly(cs, lp2, lpk)
            for (p2, pk), cc in poly.terms():
                nl = pk + w.count('lp')            # number of free lp vectors
                if nl % 2:
                    continue
                P = p2 + nl//2                     # total power of lp^2 after the angular average
                # angular average: (lp^i1 ... lp^iN) -> lp^(N) sym / (n (n+2) ... (n+N-2))
                denom = sp.Integer(1)
                for j in range(nl//2):
                    denom *= (n + 2*j)
                items = [('kdot', i) for i in range(pk)] + [('w', i) for i, t in enumerate(w) if t == 'lp']
                acc = {}
                for match in _pairings(items):
                    ww = list(w)
                    sc = sp.Integer(1)
                    lab = 0
                    for a_, b_ in match:
                        if a_[0] == 'kdot' and b_[0] == 'kdot':
                            sc *= k2
                        elif a_[0] == 'kdot' or b_[0] == 'kdot':
                            pos = b_[1] if a_[0] == 'kdot' else a_[1]
                            ww[pos] = kvec
                        else:
                            lab += 1
                            ww[a_[1]] = ('i', 1000 + lab); ww[b_[1]] = ('i', 1000 + lab)
                    for w2, c2 in perp_contract(tuple(ww), {}).items():
                        c3, st = perp_canonical(w2, k2)
                        acc[st] = acc.get(st, 0) + sc*c2*c3
                # radial integral: int d^n lp/(2pi)^n (lp^2)^P / (A lp^2 + M k2)^m * Gamma(m)
                rad = (4*sp.pi)**(-n/2)*sp.gamma(P + n/2)*sp.gamma(m - P - n/2)/sp.gamma(n/2) \
                    * A**(-P - n/2)*(M)**(P + n/2 - m)
                for st, v in acc.items():
                    out[st] += sgn_total*fac*cc*v/denom*rad
    # every term scales as k2^(-1-eps) (the words carry one power of k): return the value at k2 = 1
    res = {st: expr.subs(k2, 1) for st, expr in out.items()}
    return res, params


def lc_integrate_multi(nums, dens, var, signs_at):
    """as lc_integrate, for a dict of numerators sharing the same denominators:
    returns list of ({key: num'}, dens')"""
    dens = [(sp.cancel(sp.together(e)), s) for e, s in dens]
    coef = [sp.cancel(sp.diff(e, var)) for e, s in dens]
    vd = [(e, s, b) for (e, s), b in zip(dens, coef) if b != 0]
    other = [(e, s) for (e, s), b in zip(dens, coef) if b == 0]
    if not vd:
        return []
    for key, num in nums.items():
        nn = sp.expand(num)
        if nn == 0:
            continue
        deg = sp.Poly(nn, var).degree()
        if deg - len(vd) > -2:
            raise ValueError(f'integrand falls off too slowly in {var}: degree {deg}, {len(vd)} denominators')
    poles = []
    for i, (e, s, b) in enumerate(vd):
        loc = sp.N((s/b).subs(signs_at))
        poles.append((i, 'lower' if loc > 0 else 'upper'))
    nlow = sum(1 for p in poles if p[1] == 'lower')
    side = 'lower' if nlow <= len(poles) - nlow else 'upper'
    pref = -sp.I if side == 'lower' else sp.I
    out = []
    for i, where in poles:
        if where != side:
            continue
        e, s, b = vd[i]
        v0 = sp.cancel(var - e/b)              # e = b var + rest -> var = -rest/b
        newd = []
        for j, (ej, sj, bj) in enumerate(vd):
            if j == i:
                continue
            newd.append((sp.cancel(sp.together(ej.subs(var, v0))), sp.cancel(sj - bj*s/b)))
        newn = {key: sp.cancel(sp.together(pref*num.subs(var, v0)/b)) for key, num in nums.items()}
        out.append((newn, other + newd))
    return out
