"""The parametric pipeline for g -> gg, piece by piece: normal form against the integrand at random points
(eik weight z^(-1-eta) included), and the sector integrals against direct tanh-sinh quadrature at real
(eps, eta) where they converge.

    python3 check_steps_gg.py 3/10 [names]
"""
import sys, warnings, pickle
warnings.filterwarnings('ignore')
import numpy as np
import sympy as sp
import scipy.special as ss
from lc_engine import eps, D
from factors import eta, normal_form
from secdec import decompose, Sector
import gg_diagrams as G
from gg_driver import z, STRUCT
from numeric_checks import close,interior_values

xval = sp.Rational(sys.argv[1])
tag = str(xval).replace('/', 'o')
only = sys.argv[2].split(',') if len(sys.argv) > 2 else None
data = pickle.load(open(f'gg_integrands_x{tag}.pkl', 'rb'))
EPS, ETA = sp.Rational(-3, 10), sp.Rational(-3, 10)
rng = np.random.default_rng(2)


def ts_nodes(h=1/16, tmax=3.6):
    t = np.arange(-tmax, tmax + h/2, h)
    u = 0.5*np.pi*np.sinh(t)
    xx = 1/(1 + np.exp(-2*u))
    wgt = h*0.5*np.pi*np.cosh(t)/(2*np.cosh(u)**2)
    interior=(xx>0)&(xx<1)
    return xx[interior], wgt[interior]


xs, ws = ts_nodes()


def cube_integral(expr, vs):
    grids = np.meshgrid(*([xs]*len(vs)), indexing='ij')
    W = np.ones_like(grids[0])
    for k in range(len(vs)):
        sh = [1]*len(vs); sh[k] = len(xs)
        W = W*ws.reshape(sh)
    v=interior_values(expr,vs,[g+0j for g in grids])
    return np.sum(W*v)


s_, t1, t2 = sp.symbols('s t1 t2')
worst_pt, worst_int = 0, 0
checked=0
for name, pieces in data['results'].items():
    if only and name not in only:
        continue
    eik = G.DIAGRAMS[name]['eik']
    for iv, pcs in pieces:
        for prm, o in pcs:
            for st in STRUCT:
                if o[st] == 0:
                    continue
                za, zb = sp.Rational(iv[0]), sp.Rational(iv[1]); dz = zb - za
                e = sp.sympify(o[st]).subs(D, 2 - 2*eps)
                jac = dz
                if len(prm) == 2:
                    a1, a2 = prm
                    e = e.subs({a1: t1, a2: (1 - t1)*t2}, simultaneous=True); jac *= (1 - t1); cube = [s_, t1, t2]
                elif len(prm) == 1:
                    e = e.subs(prm[0], t1); cube = [s_, t1]
                else:
                    cube = [s_]
                e = e.subs(z, za + dz*s_)*jac
                if eik:
                    e = e*(za + dz*s_)**(-1 - eta)
                const, facs = normal_form(e, cube, {v: sp.Rational(2, 7) for v in cube})
                nf = const*sp.Mul(*[f**ex for f, ex in facs.items()])
                e0 = e.subs({eps: EPS, eta: ETA}); nf0 = nf.subs({eps: EPS, eta: ETA})
                fe = sp.lambdify(cube, e0, [{'gamma': ss.gamma}, 'numpy'])
                fn = sp.lambdify(cube, nf0, [{'gamma': ss.gamma}, 'numpy'])
                pts = rng.random((4, len(cube)))
                dev = max(abs(complex(fn(*p))/complex(fe(*p)) - 1) for p in pts)
                close(dev,0,atol=1e-10,rtol=0,label=f'{name} {st} normal form')
                checked+=1
                worst_pt = max(worst_pt, dev)
                line = f'{name:6s} {st} {str(iv):14s} {len(prm)}p: |nf/orig - 1| {dev:.1e}'
                if len(cube) <= 2:
                    I_orig = cube_integral(e0, cube)
                    secs = [Sector(j, m, r, cube, const) for j, m, r, vm in decompose(facs, cube)]
                    I_val = sum(S.value(complex(EPS), complex(ETA), n=64) for S in secs)
                    rel = abs(I_val - I_orig)/max(abs(I_orig), 1e-12)
                    close(I_val,I_orig,label=f'{name} {st} sector integral')
                    worst_int = max(worst_int, rel)
                    line += f'   integral: direct {I_orig.real:+.8f}{I_orig.imag:+.1e}i  sectors {I_val.real:+.8f}  rel {rel:.1e}'
                print(line, flush=True)
print(f'worst pointwise {worst_pt:.1e}, worst integral {worst_int:.1e}')
if not checked:raise ValueError('No nonzero integrands were compared')
