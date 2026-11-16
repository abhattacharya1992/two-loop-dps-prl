#!/usr/bin/env python3
# Numerical evaluation of the two-loop class-B gg->HH finite remainders of
# Banerjee, Borowka, Dhani, Gehrmann, Ravindran (arXiv:1809.05388), read
# directly from the authors' ancillary files anc/Gauge{1,2}Simp.in.
# Returns the raw Gauge_i expression; the dimensionless remainder used everywhere
# (paper Sec. 3.1) is A_i = i Gauge_i/(4 M2), the authors' own plotting convention
# (eval_gauge*.m in their ancillary files).
#
# Kinematics: s/M2 = (1+x)^2/x, -t/M2 = y, -u/M2 = z = x + 1/x - y.
# Forward limit at fixed beta: x = (1-beta)/(1+beta), y = x + delta.
# The only non-standard function in the files is Li22new[a,b], the weight-4
# MPL Li_{2,2}(a,b) = sum_{m>n>=1} a^m b^n/(m^2 n^2) of the lievaluate
# package (Frellesvig, Tommasini, Wever). It is computed here from the
# one-dimensional integral in li22() below (with integration points that accumulate at
# the end of the path, where a branch point approaches it in the forward limit; see
# there and check_li22_values.py). Which argument carries the
# outer sum is a convention; the ordering adopted here is the one that
# reproduces Bose symmetry and the Li_{2,2} series and stuffle identities, all
# checked in an earlier evaluation (tests/test_release.py --full repeats one of those
# points). It has not been checked against lievaluate itself, which is not available
# without Mathematica.
import ast, pathlib, sys, time, functools
import mpmath as mp

DIR = pathlib.Path(__file__).parent
_CACHE_DPS = None   # the lru caches below hold mpf values, so they are flushed whenever mp.dps changes

# The .in files are Mathematica. After replacing ^ -> **, [ ] -> ( ) they are
# valid Python expressions in these symbols only. The AST walk below is a
# check that the text substitution produced nothing but arithmetic and the
# four expected function heads before the expression is handed to exec().
VALID_FUNCTIONS = {'Log', 'PolyLog', 'Li22new', 'Zeta'}
VALID_PARAMS = {'x', 'y', 'M2', 'muR2', 'N', 'Nf', 'I', 'Pi'}
VALID_OPS = (ast.Expression, ast.Module, ast.Expr, ast.Assign, ast.Load, ast.Store, ast.Name,
             ast.Constant, ast.Call, ast.BinOp, ast.UnaryOp, ast.Add, ast.Sub, ast.Mult, ast.Div,
             ast.Pow, ast.UAdd, ast.USub)


def compile_amp(n):
    raw = (DIR/'anc'/f'Gauge{n}Simp.in').read_text(encoding='utf8').split('=', 1)[1]
    raw = raw.replace('^', '**').replace('[', '(').replace(']', ')')
    root = ast.parse('result='+raw)
    assert all(isinstance(z, VALID_OPS) for z in ast.walk(root)), 'unexpected Mathematica syntax'
    names = {z.id for z in ast.walk(root) if isinstance(z, ast.Name)}
    allowed = VALID_FUNCTIONS | VALID_PARAMS | {'result'}
    assert names <= allowed, names-allowed
    for node in ast.walk(root):
        if isinstance(node, ast.Call):
            assert isinstance(node.func, ast.Name) and node.func.id in VALID_FUNCTIONS
    return compile(root, f'Gauge{n}Simp.in', 'exec')


LI22_GUARD = 15     # digits added to the working precision inside li22()
LI22_RATIO = 4      # the integration points 1 - d of li22_points() shrink d by this factor


def li22_points(a, w):
    """0, 1/2 and points 1 - d, d = 1/2, 1/8, 1/32, ..., accumulating at t = 1 until d is no larger than
    the distance from t = 1 of the nearest singular point of the integrand beyond it: the pole
    t = 1/a of 1/(1 - a t) and the branch point t = 1/w of Li_2(w t), for 0 < a < 1, 0 < w < 1. Every
    interval then ends at least 2/3 of its length before the singular point. The same at t = 0 for
    a singular point in (-1, 0), that is a < -1 or w < -1 (not met on the forward lines used here)."""
    near1, near0 = mp.mpf(1), mp.mpf(1)
    for v in (a, w):
        if 0 < v < 1:
            near1 = min(near1, 1/v - 1)
        elif v < -1:
            near0 = min(near0, -1/v)
    points, d = [mp.mpf(1)/2], mp.mpf(1)/2
    while d > near0:
        d /= LI22_RATIO
        points.insert(0, d)
    points.insert(0, mp.mpf(0))
    d = mp.mpf(1)/2
    while d > near1:
        d /= LI22_RATIO
        points.append(1 - d)
    points.append(mp.mpf(1))
    return points


@functools.lru_cache(maxsize=None)
def li22(a, b):
    # Li_{2,2}(a,b) = -a int_0^1 log(t) Li_2(abt)/(1-at) dt (expand 1/(1-at) and Li_2 and integrate).
    # The integrand has a pole at t = 1/a and Li_2 has its cut at abt >= 1, so
    # the real path [0,1] is usable only when neither a nor ab reaches 1.
    # The guards below stop the evaluation if a call ever lands on the cut; none does on
    # the forward lines used here, where all arguments are real (complex ones are refused, since
    # the integration points below are chosen for real arguments). For a < 0 or ab < 0 the
    # integral continues smoothly.
    # Towards the forward limit the branch point t = 1/(ab) approaches t = 1 from beyond: ab = 1 - eps
    # with b > 1 and eps of the order of delta occurs. mpmath's tanh-sinh quadrature then reports
    # convergence before it resolves the branch point; split only at t = 1/2, the integral lost up
    # to 21 digits there (check_li22_values.py). It is therefore split at the points of li22_points()
    # and evaluated with LI22_GUARD more digits, because mpmath's error estimate remains optimistic
    # by a few digits. Written as a w int_0^1 (-log t) [Li_2(w t)/w]/(1 - a t) dt, w = ab, the
    # integrand is positive and of order one, so that mpmath's absolute error target is a relative one.
    if mp.im(a) != 0 or mp.im(b) != 0:
        raise ValueError('Li22 with a complex argument: the integration points are chosen for real ones')
    a, b = mp.re(a), mp.re(b)
    with mp.workdps(mp.mp.dps + LI22_GUARD):
        w = a*b
        if a >= 1: raise ValueError('Li22 first argument on integration cut')
        if w >= 1: raise ValueError('Li22 product on integration cut')
        if a == 0 or w == 0:
            return mp.mpf(0)

        def integrand(t):
            if not t: return mp.mpf(0)
            return -mp.log(t)*(mp.polylog(2, w*t)/w)/(1 - a*t)
        value = a*w*mp.quad(integrand, li22_points(a, w))
    return +value


@functools.lru_cache(maxsize=None)
def LiClassical(n, z): return mp.polylog(int(n), z)
@functools.lru_cache(maxsize=None)
def Log(x): return mp.log(x)   # principal branch; log(-r) = log r + i pi
@functools.lru_cache(maxsize=None)
def Zeta(n): return mp.zeta(n)


def eval_amp(code, x, y, M2=1, muR2=mp.mpf('0.5'), N=3, Nf=5):
    global _CACHE_DPS
    if _CACHE_DPS != mp.mp.dps:
        for fun in (li22, LiClassical, Log, Zeta): fun.cache_clear()
        _CACHE_DPS = mp.mp.dps
    env = dict(x=mp.mpf(x), y=mp.mpf(y), M2=mp.mpf(M2), muR2=mp.mpf(muR2), N=mp.mpf(N), Nf=mp.mpf(Nf),
               I=mp.j, Pi=mp.pi, Log=Log, PolyLog=LiClassical, Li22new=li22, Zeta=Zeta)
    exec(code, {'__builtins__': {}}, env)
    return env['result']


def direct_forward(n, beta, delta, dps=40, M2=1, muR2=mp.mpf('0.5'), N=3, Nf=5):
    mp.mp.dps = dps
    b = mp.mpf(str(beta)); x = (1-b)/(1+b); y = x+mp.mpf(str(delta))
    code = compile_amp(n)
    t = time.time(); a = eval_amp(code, x, y, M2, muR2, N, Nf)
    return x, y, a, time.time()-t


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='one forward-limit point of Gauge1/2 at fixed beta')
    parser.add_argument('--gauge', type=int, choices=(1, 2), default=1)
    parser.add_argument('--beta', type=str, default='0.4')
    parser.add_argument('--delta', type=str, default='0.05')
    parser.add_argument('--dps', type=int, default=40)
    parser.add_argument('--muR2', default='0.5', help='mu_R^2/M2; the authors use 1/2')
    parser.add_argument('--nf', type=int, default=5)
    args = parser.parse_args()
    x, y, a, t = direct_forward(args.gauge, args.beta, args.delta, args.dps, muR2=args.muR2, Nf=args.nf)
    print('x =', mp.nstr(x, 20), 'y =', mp.nstr(y, 20), 'dps =', mp.mp.dps)
    print('Gauge expression / M2 =', mp.nstr(a, 30))
    print('A = i/4 * Gauge / M2 =', mp.nstr(mp.j*a/4, 30))
    print('evaluation time (s)', round(t, 2))
    print('cache hits/misses', li22.cache_info(), LiClassical.cache_info(), Log.cache_info())
