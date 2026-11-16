"""Abelian vertex correction to J1(g -> q qbar): diagram Vqq, a gluon exchanged between the quark and the antiquark.

Conventions of paper2/r0/theory/gqq_diagrams.py: p1 = (1, 0, 0), K = (x, K^-, k), l = (z, l^-, l_perp), q = K - l,
v^2 = v+ v- - v_perp^2. Units k^2 = mu^2 = 1, a_s = alpha_s/(4 pi), colour factor C_F - C_A/2.
The two helicity structures eslash kslash and kslash eslash are written as
    Vqq = A T + (B/2) (e.k),   T = (x/2) eslash kslash - ((1-x)/2) kslash eslash   (the tree splitting),
so that their coefficients, normalised to the tree splitting, are A - B/(2x) and A + B/(2(1-x)).
Result:   A = 1/eps - ln(x(1-x)),   B = -4 x ln x + 4 (1-x) ln(1-x).

Steps checked here:
  1  Dirac algebra in d = 4, 6, 8, 10 dimensions (the result is linear in d, so this fixes it for all d)
  2  both minus components by residues, against quadrature in two variables
  3  charge conjugation, which gives the region -(1-x) < z < 0 from the region 0 < z < x
  4  the parameter integrals at finite eps against their Laurent expansion
  5  the result against the stored Laurent coefficients
Usage: python3 check_qq_Vqq.py <repository root>"""
import json, sys, pathlib
import numpy as np
from scipy.special import gamma as Gam
from scipy.integrate import quad
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from dirac import Space, charge_conjugation, dot, mx
from report import check, require, finish
repo = pathlib.Path(sys.argv[1])
EG = 0.5772156649015329
rng = np.random.default_rng(7)


def string(sp, x, u, k, z, w, l, e, contracted=True):
    """P Kslash gamma^a (K-l)slash eslash (K-l-p1)slash gamma_a (K-p1)slash Pbar, or the same without gamma^a ... gamma_a"""
    q = k - l; up = u - w
    inner = sp.slash(x - z, up, q) @ sp.perp(e) @ sp.slash(x - z - 1, up, q)
    if contracted: inner = sp.contract(inner)
    return sp.P @ sp.slash(x, u, k) @ inner @ sp.slash(x - 1, u, k) @ sp.Pb


def N12(sp, x, k, z, q, e, beta):
    """coefficients of nslash: N1 from (q-p1)slash eslash qslash and N2 from qslash eslash (q-p1)slash; beta = q^-/2"""
    xb = 1 - x; al = (x - z)/2; alp = al - 0.5
    K_, Q_, E_ = sp.perp(k), sp.perp(q), sp.perp(e); QEQ = Q_ @ E_ @ Q_; eq = dot(e, q)*sp.one
    N1 = 2*x*xb*beta*eq + x/2*(QEQ - 4*al*beta*E_) @ K_ - xb/2*K_ @ (QEQ - 4*alp*beta*E_) - K_ @ (alp*E_ @ Q_ + al*Q_ @ E_) @ K_
    N2 = 2*x*xb*beta*eq + x/2*(QEQ - 4*alp*beta*E_) @ K_ - xb/2*K_ @ (QEQ - 4*al*beta*E_) - K_ @ (al*E_ @ Q_ + alp*Q_ @ E_) @ K_
    return N1, N2


print('1  Dirac algebra in d dimensions: gamma^a X gamma_a = -2 (reversed X) + (4 - d) X for three gamma matrices')
x, z = 0.37, 0.11; xb = 1 - x; ze = z/x; u, w = 0.83, -0.41
dev = {k_: 0.0 for k_ in ('string', 'K-', 'q-', 'pole', 'tensor')}
for d in (4, 6, 8, 10):
    sp = Space(d); k, l, e = (rng.normal(size=d - 2) for _ in range(3)); q = k - l; q2, k2 = dot(q, q), dot(k, k)
    N1, N2 = N12(sp, x, k, z, q, e, (u - w)/2)
    dev['string'] = max(dev['string'], mx(string(sp, x, u, k, z, w, l, e) - sp.nsl @ (-2*N1 + (4 - d)*N2)),
                        mx(string(sp, x, u, k, z, w, l, e, contracted=False) - sp.nsl @ N2))
    dev['K-'] = max(dev['K-'], mx(string(sp, x, u, k, z, w, l, e) - string(sp, x, u + 1.7, k, z, w + 1.7, l, e)))
    fa, fb, fc = (string(sp, x, u, k, z, u - up, l, e) for up in (0.2, 0.9, 1.6))
    dev['q-'] = max(dev['q-'], mx(fa - 2*fb + fc))
    # on the pole q^- = q_perp^2/(x - z)
    al = x*(1 - ze)/2; alp = al - 0.5; K_, Q_, E_ = sp.perp(k), sp.perp(q), sp.perp(e); eq = dot(e, q)
    N1p, N2p = N12(sp, x, k, z, q, e, q2/(2*(x - z)))
    N1h = xb*q2/(1 - ze)*(eq*sp.one - K_ @ E_/(2*x)) - eq*(x*Q_ @ K_ - xb*K_ @ Q_) - K_ @ (alp*E_ @ Q_ + al*Q_ @ E_) @ K_
    N2h = q2/(1 - ze)*(xb*eq*sp.one + E_ @ K_/2) - eq*(x*Q_ @ K_ - xb*K_ @ Q_) - K_ @ (al*E_ @ Q_ + alp*Q_ @ E_) @ K_
    dev['pole'] = max(dev['pole'], mx(N1p - N1h), mx(N2p - N2h))
    # shift q = a k - l' and average over the directions of l' in d - 2 dimensions
    a, r = 0.43, 0.77; pts = sp.basis(r)
    avg = lambda f: sum(f(a*k - lp) for lp in pts)/len(pts)
    t1 = avg(lambda qq: dot(e, qq)*sp.perp(qq) @ K_); t2 = avg(lambda qq: dot(e, qq)*K_ @ sp.perp(qq))
    t3 = avg(lambda qq: K_ @ E_ @ sp.perp(qq) @ K_); t4 = avg(lambda qq: K_ @ sp.perp(qq) @ E_ @ K_)
    ek = dot(e, k)*sp.one; n = d - 2
    dev['tensor'] = max(dev['tensor'], mx(t1 - (-a*a*k2*ek + r*r/n*E_ @ K_)), mx(t2 - (-a*a*k2*ek + r*r/n*K_ @ E_)),
                        mx(t3 + a*k2*K_ @ E_), mx(t4 + a*k2*E_ @ K_))
check('Dirac string = nslash [-2 N1 + (4 - d) N2], d = 4, 6, 8, 10', dev['string'], 1e-10)
check('numerator independent of K^- at fixed q^-', dev['K-'], 1e-10)
check('numerator linear in q^-', dev['q-'], 1e-10)
check('N1 and N2 on the pole q^- = q_perp^2/(x - z)', dev['pole'], 1e-10)
check('tensor reduction to the two helicity structures in d - 2 transverse dimensions', dev['tensor'], 1e-10)

print('2  both minus components by residues, against quadrature in two variables with finite imaginary parts')


def lc_numeric(x, z, m, p):
    xb = 1 - x
    def inner(up):
        f = lambda uu: 1/((x*uu - m[0])*(-xb*uu - m[1])*(z*(uu - up) - m[4]))
        return quad(f, -np.inf, np.inf, complex_func=True, epsabs=1e-12, epsrel=1e-11, limit=400)[0]
    g = lambda up: up**p*inner(up)/(((x - z)*up - m[2])*(-(xb + z)*up - m[3]))
    return quad(g, -np.inf, np.inf, complex_func=True, epsabs=1e-10, epsrel=1e-9, limit=400)[0]


def lc_residue(x, z, m, p):
    xb = 1 - x
    if 0 < z < x:            # l^- on the pole of (K-p1)^2, q^- on the pole of q^2
        u3 = m[2]/(x - z)
        return (2j*np.pi/(x*m[1] + xb*m[0]))*(-2j*np.pi/(x - z))*u3**p/((-(xb + z)*u3 - m[3])*(-z*u3 - z*m[1]/xb - m[4]))
    if -xb < z < 0:          # l^- on the pole of K^2, q^- on the pole of (q-p1)^2
        u4 = -m[3]/(xb + z)
        return (2j*np.pi/(xb*m[0] + x*m[1]))*(-2j*np.pi/(xb + z))*u4**p/(((x - z)*u4 - m[2])*(-z*u4 + z*m[0]/x - m[4]))
    return 0.0


m = np.array([1.0 - 0.40j, 1.0 - 0.50j, 0.7 - 0.30j, 0.7 - 0.60j, 0.45 - 0.35j])
for zz, where in ((0.12, '0 < z < x'), (-0.25, '-(1-x) < z < 0'), (0.45, 'z > x'), (-0.85, 'z < -(1-x)')):
    d_ = max(abs(lc_numeric(0.3, zz, m, p) - lc_residue(0.3, zz, m, p)) for p in (0, 1))
    check(f'x = 0.3, z = {zz:+.2f} ({where}), numerators 1 and q^-', d_, 1e-9)
m0 = np.array([1.0, 1.0, 0.7, 0.7, 0.45]); zz = 0.12; zeta = zz/0.3
closed = (1 - zeta)/(1.0*0.7*((1 - zeta)*0.45 + zeta*0.7 - zeta*(1 - zeta)*1.0 + zeta*(1 - zeta)*1.0/0.7))
check('residues/(2 pi)^2 = (1-zeta)/(k^2 q_perp^2 [(l_perp - zeta k)^2 + zeta(1-zeta) k^2/(1-x)])',
      (lc_residue(0.3, zz, m0, 0)/(2*np.pi)**2).real - closed, 1e-12)

print('3  charge conjugation: C String(K, l) C^-1 = -String(p1 - K, -l)^T')
sp = Space(4); C = charge_conjugation(sp); k, l, e = (rng.normal(size=2) for _ in range(3))
def raw(x, u, k, z, w, l):
    q = k - l; up = u - w
    return sp.slash(x, u, k) @ sp.contract(sp.slash(x - z, up, q) @ sp.perp(e) @ sp.slash(x - z - 1, up, q)) @ sp.slash(x - 1, u, k)
check('identity for random momenta', mx(C @ raw(x, u, k, z, w, l) @ np.linalg.inv(C) + raw(1 - x, -u, -k, -z, -w, -l).T), 1e-12)

print('4  region 0 < z < x: parameter integrals at finite eps against the Laurent expansion')


def region(x, ep):
    """coefficients of (eslash kslash, kslash eslash) of e^(eps gamma_E) Int dz Int d^(2-2eps)l/pi^(1-eps) N/(q^2 z E)"""
    xb = 1 - x; Tv = np.array([x/2, -xb/2]); ek = np.array([-0.5, -0.5])
    D = lambda t, ze: t*(1 - ze) + ze/xb
    tint = lambda ze, p1, p2: quad(lambda t: (1 - t)**p1*D(t, ze)**p2, 0, 1, epsabs=1e-13, epsrel=1e-12, limit=200)[0]
    def integrand(ze, comp, which):
        al = x*(1 - ze)/2; alp = al - 0.5
        J0 = (1 - ze)**(-ep)*tint(ze, -ep, -ep)
        J1 = (1 - ze)**(1 - ep)*tint(ze, 1 - ep, -1 - ep)
        J2 = (1 - ze)**(-ep)*tint(ze, -ep, -1 - ep)
        D0 = (ze*(1 - ze)/xb)**(-ep)
        if which == 1: Na = xb*ek + np.array([0, -xb/(2*x*(1 - ze))]); lin = np.array([al, alp])
        else: Na = xb*ek + np.array([1/(2*(1 - ze)), 0]); lin = np.array([alp, al])
        val = Gam(ep)*D0*Na - Tv*Gam(ep)*J0 + Gam(1 + ep)*((x - xb)*ek*J1 + lin*J2)
        return x*(1 - ze)*val[comp]
    I = lambda which: np.array([quad(integrand, 0, 1, args=(c_, which), epsabs=1e-10, epsrel=1e-10, limit=200)[0] for c_ in (0, 1)])
    return np.exp(EG*ep)*(-2*I(1) + 2*ep*I(2))


def region_series(x):
    """pole and finite part of the region"""
    xb = 1 - x
    return np.array([x/2, xb]), np.array([2*x + (1 + xb)/2*np.log(xb), 2*xb + xb/2*np.log(xb)])


for x in (0.3, 0.7):
    pole, fin = region_series(x)
    d_ = [region(x, ep) - pole/ep - fin for ep in (0.04, 0.02, 0.01)]
    for c_, name in ((0, 'eslash kslash'), (1, 'kslash eslash')):
        check(f'x = {x}, {name}: difference at eps = 0.01', d_[2][c_], 0.05)
        require(f'x = {x}, {name}: difference halves with eps (ratios {d_[0][c_]/d_[1][c_]:.3f}, {d_[1][c_]/d_[2][c_]:.3f})',
                1.9 < d_[0][c_]/d_[1][c_] < 2.1 and 1.9 < d_[1][c_]/d_[2][c_] < 2.1)
    # the two regions together: the second is the first at 1 - x with the structures exchanged and a sign
    xb = 1 - x; p2, f2 = region_series(xb)
    Vp = -(pole + np.array([-p2[1], -p2[0]])); Vf = -(fin + np.array([-f2[1], -f2[0]]))
    check(f'x = {x}: pole of A equals 1', 2*(Vp[0] - Vp[1]) - 1, 1e-13)
    check(f'x = {x}: pole of B vanishes', -4*xb*Vp[0] - 4*x*Vp[1], 1e-13)
    check(f'x = {x}: finite part of A equals -ln(x(1-x))', 2*(Vf[0] - Vf[1]) + np.log(x*xb), 1e-13)
    check(f'x = {x}: B equals -4 x ln x + 4 (1-x) ln(1-x)', -4*xb*Vf[0] - 4*x*Vf[1] - (-4*x*np.log(x) + 4*xb*np.log(xb)), 1e-13)

print('5  result against the stored Laurent coefficients (five x, two helicity structures, six orders)')
worst = 0
for tag, x in {'1o10': 0.1, '3o10': 0.3, '1o2': 0.5, '7o10': 0.7, '9o10': 0.9}.items():
    st = json.load(open(repo/f'paper2/r0/theory/gqq_laurent_x{tag}.json'))['diagrams']['Vqq']['laurent']
    y = 1 - x; A = {'0,-1': 1.0, '0,0': -np.log(x*y)}; B = {'0,0': -4*x*np.log(x) + 4*y*np.log(y)}
    for key in ('-1,-2', '-1,-1', '-1,0', '0,-2', '0,-1', '0,0'):
        a, b = A.get(key, 0), B.get(key, 0)
        worst = max(worst, abs(st['ek'][key][0] - (a - b/(2*x))), abs(st['ke'][key][0] - (a + b/(2*y))))
check('largest deviation (the stored values of this diagram are accurate to about 1e-7)', worst, 5e-7)
finish()
