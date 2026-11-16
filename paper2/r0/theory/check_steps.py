"""Locate the error in the parametric pipeline, piece by piece, at one real eps where everything converges:
(1) original integrand on the cube vs normal form const * prod f^e at random points;
(2) direct tanh-sinh integral of the normal form over the cube;
(3) sum over sectors of the tanh-sinh integral of each sector integrand (no subtraction needed);
(4) Sector.value (with subtractions, Gauss-Legendre).

    python3 check_steps.py 3/10 Vqq ek -3/10
"""
import sys, warnings, itertools
warnings.filterwarnings('ignore')
import numpy as np
import sympy as sp
import scipy.special as ss
from lc_engine import eps, D
from factors import eta, normal_form
from secdec import decompose, Sector
from gqq_diagrams import z
import assemble_gqq as A
from numeric_checks import close,interior_values

args = [a for a in sys.argv[1:] if a != '--allow-legacy-input']
name, st = args[1], args[2]
epsv = sp.Rational(args[3])
A.setup(sp.Rational(args[0]), allow_legacy='--allow-legacy-input' in sys.argv)
rng = np.random.default_rng(1)


def ts_nodes(h=1/16, tmax=3.6):
    t = np.arange(-tmax, tmax + h/2, h)
    u = 0.5*np.pi*np.sinh(t)
    x = 1/(1 + np.exp(-2*u))
    wgt = h*0.5*np.pi*np.cosh(t)/(2*np.cosh(u)**2)
    interior=(x>0)&(x<1)
    return x[interior], wgt[interior]


xs, ws = ts_nodes()


def cube_integral(expr, vs):
    grids = np.meshgrid(*([xs]*len(vs)), indexing='ij')
    W = np.ones_like(grids[0])
    for k in range(len(vs)):
        sh = [1]*len(vs); sh[k] = len(xs)
        W = W*ws.reshape(sh)
    v=interior_values(expr,vs,[g+0j for g in grids])
    return np.sum(W*v)


s_, t1, t2 = A.s_, A.t1, A.t2
checked=0
for iv, pcs in A.data['results'][name]:
    for prm, o in pcs:
        if o[st] == 0:
            continue
        expr = o[st]
        za, zb = sp.Rational(iv[0]), sp.Rational(iv[1])
        dz = zb - za
        e = expr.subs(D, 2 - 2*eps)
        jac = dz
        if len(prm) == 2:
            a1, a2 = prm
            e = e.subs({a1: t1, a2: (1 - t1)*t2}, simultaneous=True); jac *= (1 - t1); cube = [s_, t1, t2]
        elif len(prm) == 1:
            e = e.subs(prm[0], t1); cube = [s_, t1]
        else:
            cube = [s_]
        e = e.subs(z, za + dz*s_)*jac
        test = {v: sp.Rational(2, 7) for v in cube}
        const, facs = normal_form(e, cube, test)
        nf = const*sp.Mul(*[f**ex for f, ex in facs.items()])
        # (1) pointwise
        e0 = e.subs(eps, epsv); nf0 = nf.subs(eps, epsv)
        fe = sp.lambdify(cube, e0, [{'gamma': ss.gamma}, 'numpy'])
        fn = sp.lambdify(cube, nf0, [{'gamma': ss.gamma}, 'numpy'])
        pts = rng.random((5, len(cube)))
        rat = [complex(fn(*p))/complex(fe(*p)) for p in pts]
        close(np.asarray(rat),np.ones(len(rat)),atol=1e-10,rtol=0,label='normal form pointwise')
        # (2) direct integral of the original and of the normal form
        I_orig = cube_integral(e0, cube)
        I_nf = cube_integral(nf0, cube)
        # (3) sectors, integrated directly
        I_sec = 0
        secs = []
        for j, m, r, vm in decompose(facs, cube):
            integrand = j*sp.Mul(*[v**ex for v, ex in m.items()])*sp.Mul(*[Q**ex for Q, ex in r.items()])*const
            I_sec += cube_integral(integrand.subs({eps: epsv, eta: 0}), cube)
            secs.append(Sector(j, m, r, cube, const))
        # (4) Sector.value
        I_val = sum(S.value(complex(epsv), 0.0, n={1: 64, 2: 48, 3: 28}[len(S.vs)]) for S in secs)
        for label,value in [('normal form',I_nf),('sectors',I_sec),('Sector.value',I_val)]:
            close(value,I_orig,label=label)
        checked+=1
        print(f'{iv} {len(prm)}p: nf/orig at points {[f"{r.real:.6f}" for r in rat]}', flush=True)
        print(f'     orig {I_orig.real:+.8f}  normal form {I_nf.real:+.8f}  sectors (tanh-sinh) {I_sec.real:+.8f}'
              f'  Sector.value {I_val.real:+.8f}', flush=True)
if not checked:raise ValueError('No nonzero integrands were compared')
