"""Non-Abelian vertex correction and Wilson line couplings to the parent gluon in J1(g -> q qbar):
diagrams Vggq, W1g and W2g, all with colour factor C_A/2.

Conventions of paper2/r0/theory/gqq_diagrams.py: p1 = (1, 0, 0), K = (x, K^-, k), l = (z, l^-, l_perp), q = K - l.
Units k^2 = mu^2 = nu_r = 1, a_s = alpha_s/(4 pi). Each diagram is (C_A/2) [A T + (B/2) (e.k)] with the tree
splitting T = (x/2) eslash kslash - ((1-x)/2) kslash eslash. Results at b = b0, with L = ln(x(1-x)):
    Vggq:  A = (2 - 2L)/eps + 4 - ln^2 x - ln^2(1-x),
           B = (-2x ln x + 2(1-x) ln(1-x))/eps - x ln^2 x + (1-x) ln^2(1-x)
    W1g:   A = -2/eps^2 - (2 ln(1-x) + x)/eps - 3 zeta2 - ln(1-x) - ln^2(1-x) + 2 Li2(x) - 2x,
           B = (2x(1-x) - 2(1-x) ln(1-x))/eps + 4x(1-x) - (1-x) ln^2(1-x)
    W2g:   W1g with x <-> 1-x and B -> -B (charge conjugation)

Steps checked here:
  1  Dirac algebra in d = 4, 6, 8, 10 dimensions and charge conjugation
  2  both minus components by residues, against quadrature in two variables
  3  the parameter integrals at finite eps against their Laurent expansions
  4  the results against the stored Laurent coefficients
Usage: python3 check_qq_Vggq_W1g_W2g.py <repository root>"""
import json, sys, pathlib, warnings
import numpy as np
from scipy.special import gamma as Gam, spence, beta as Beta
from scipy.integrate import quad
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from dirac import Space, charge_conjugation, dot, mx
from report import check, require, finish
warnings.filterwarnings('ignore')
repo = pathlib.Path(sys.argv[1])
EG, Z2 = 0.5772156649015329, np.pi**2/6
Li2 = lambda t: spence(1 - t)
rng = np.random.default_rng(11)


def raw_V(sp, x, u, k, z, w, l, e):
    """Kslash [e q r + 2 (e.l) gamma^a q gamma_a + s q e] (K-p1)slash with q = K-l, r = p1+l, s = l-2p1"""
    E = sp.perp(e); qs = sp.slash(x - z, u - w, k - l)
    c = -(E @ sp.perp(l) + sp.perp(l) @ E)
    inner = E @ qs @ sp.slash(1 + z, w, l) + c @ sp.contract(qs) + sp.slash(z - 2, w, l) @ qs @ E
    return sp.slash(x, u, k) @ inner @ sp.slash(x - 1, u, k)


def raw_W1(sp, x, u, k, z, w, l, e):
    E = sp.perp(e); c = E @ sp.perp(l) + sp.perp(l) @ E; qs = sp.slash(x - z, u - w, k - l); Qs = sp.slash(x - 1, u, k)
    return (2 - z)*qs @ E @ Qs + c @ qs @ sp.nbsl @ Qs


def raw_W2(sp, x, u, k, z, w, l, e):
    E = sp.perp(e); c = E @ sp.perp(l) + sp.perp(l) @ E; Ks = sp.slash(x, u, k); Ql = sp.slash(x - 1 + z, u + w, k + l)
    return (2 - z)*Ks @ E @ Ql + c @ Ks @ sp.nbsl @ Ql


print('1  Dirac algebra in d = 4, 6, 8, 10 dimensions (gamma^a qslash gamma_a = -(d - 2) qslash)')
x, z = 0.37, 0.11; xb = 1 - x; ze = z/x; u, w = 0.83, -0.41
al, rho, sig = (x - z)/2, (1 + z)/2, (z - 2)/2
dev = {k_: 0.0 for k_ in ('V', 'Vpole', 'Vshift', 'W', 'Wshift')}
for d in (4, 6, 8, 10):
    sp = Space(d); k, l, e = (rng.normal(size=d - 2) for _ in range(3)); q = k - l; k2, l2 = dot(k, k), dot(l, l)
    K_, L_, Q_, E_ = sp.perp(k), sp.perp(l), sp.perp(q), sp.perp(e); one = sp.one
    T = x/2*E_ @ K_ - xb/2*K_ @ E_

    def NV(beta, omega):
        """coefficient of nslash of the Vggq string; beta = q^-/2, omega = l^-/2"""
        Ni = x*xb*E_ @ (beta*L_ - omega*Q_) + x/2*(4*beta*rho*E_ + E_ @ Q_ @ L_) @ K_ - xb/2*K_ @ (E_ @ Q_ @ L_ + 4*al*omega*E_) \
            + K_ @ E_ @ (al*L_ - rho*Q_) @ K_
        Niii = -x*xb*(omega*Q_ - beta*L_) @ E_ + x/2*(4*omega*al*E_ + L_ @ Q_ @ E_) @ K_ - xb/2*K_ @ (L_ @ Q_ @ E_ + 4*sig*beta*E_) \
            - K_ @ (sig*Q_ - al*L_) @ E_ @ K_
        Nii = -2*(d - 2)*dot(e, l)*((-x*xb*beta + al*k2)*one + x/2*Q_ @ K_ - xb/2*K_ @ Q_)
        return Ni + Nii + Niii

    def Nb(lv):
        """part of the numerator on the poles that keeps both transverse denominators"""
        Lv = sp.perp(lv); el = dot(e, lv)
        return (x*k2*el*one - x*(1 + z)*k2/(2*xb)*E_ @ K_ + (z - 2)*k2/2*K_ @ E_ + x/2*(E_ @ K_ @ Lv + Lv @ K_ @ E_) @ K_
                - xb/2*K_ @ (E_ @ K_ @ Lv + Lv @ K_ @ E_) + (1 + x)/2*K_ @ E_ @ Lv @ K_ + rho*k2*K_ @ E_ + sig*k2*E_ @ K_
                - (1 + xb)/2*K_ @ Lv @ E_ @ K_ - (d - 2)*(2*x - z)*k2*el*one - 2*(d - 2)*el*((xb - x)*k2/2*one - x/2*Lv @ K_ + xb/2*K_ @ Lv))

    dev['V'] = max(dev['V'], mx(sp.P @ raw_V(sp, x, u, k, z, w, l, e) @ sp.Pb - sp.nsl @ NV((u - w)/2, w/2)))
    om = l2/(2*z); be = -k2/(2*xb) - om
    pole = om*(-2*x*(xb + z)*E_ @ K_ - 2*xb*(1 + x - z)*K_ @ E_ - 2*(d - 2)*x*xb*dot(e, l)*one) + 2*l2*T + Nb(l)
    dev['Vpole'] = max(dev['Vpole'], mx(NV(be, om) - pole))
    s_, r_ = 0.23, 0.77; pts = sp.basis(r_)
    avg = sum(Nb(s_*k + lp) for lp in pts)/len(pts)
    e1 = (-(1 + xb) + z*(xb - x))/(2*xb) + s_*(2 - 3*x)/2
    e2 = (2*z - 1)/2 + s_*(1 - 3*x)/2
    e3 = s_*x - (d - 2)*s_*(1 - z + (x - xb)*s_)
    dev['Vshift'] = max(dev['Vshift'], mx(avg - k2*(e1*E_ @ K_ + e2*K_ @ E_ + e3*dot(e, k)*one) - 2*r_*r_*T))
    NW = (2 - z)*((x - z)/2*E_ @ K_ - xb/2*Q_ @ E_) - xb*(x - z)*(E_ @ L_ + L_ @ E_)
    dev['W'] = max(dev['W'], mx(sp.P @ raw_W1(sp, x, u, k, z, w, l, e) @ sp.Pb - sp.nsl @ NW))
    ks = s_*k
    NWs = (2 - z)*((x - z)/2*E_ @ K_ - xb/2*sp.perp(k - ks) @ E_) - xb*(x - z)*(E_ @ sp.perp(ks) + sp.perp(ks) @ E_)
    An = (2 - x*ze)*(1 - x*ze - xb*s_); Bn = 2*x*xb*ze*((1 - x*ze - xb*s_) + (1 - s_))
    dev['Wshift'] = max(dev['Wshift'], mx(NWs - An*T - Bn/2*dot(e, k)*one))
check('Vggq: Dirac string = nslash N(q^-, l^-)', dev['V'], 1e-10)
check('Vggq: numerator on the poles, split by powers of l_perp^2', dev['Vpole'], 1e-10)
check("Vggq: shift l = l' + s k and tensor reduction in d - 2 transverse dimensions", dev['Vshift'], 1e-10)
check('W1g: Dirac string = nslash N, free of the minus components', dev['W'], 1e-10)
check('W1g: at l = s k, N = A_n T + (B_n/2) e.k', dev['Wshift'], 1e-10)
sp = Space(4); C = charge_conjugation(sp); Ci = np.linalg.inv(C); k, l, e = (rng.normal(size=2) for _ in range(3))
check('charge conjugation maps Vggq(K, l) to Vggq(p1 - K, p1 - l)',
      mx(C @ raw_V(sp, x, u, k, z, w, l, e) @ Ci + raw_V(sp, 1 - x, -u, -k, 1 - z, -w, -l, e).T), 1e-12)
check('charge conjugation maps W1g(K, l) to W2g(p1 - K, l)',
      mx(C @ raw_W1(sp, x, u, k, z, w, l, e) @ Ci + raw_W2(sp, 1 - x, -u, -k, z, w, l, e).T), 1e-12)

print('2  both minus components by residues, against quadrature in two variables with finite imaginary parts')
QK = dict(complex_func=True, epsabs=1e-11, epsrel=1e-10, limit=400)


def V_numeric(x, z, m, pu, pw):
    xb = 1 - x
    inner = lambda u: quad(lambda w: (u - w)**pu*w**pw/(((x - z)*(u - w) - m[2])*(z*w - m[3])*(-(1 - z)*w - m[4])), -np.inf, np.inf, **QK)[0]
    return quad(lambda u: inner(u)/((x*u - m[0])*(-xb*u - m[1])), -np.inf, np.inf, **QK)[0]


def V_residue(x, z, m, pu, pw):
    xb = 1 - x
    if 0 < z < x:          # K^- on the pole of (K-p1)^2, l^- on the pole of l^2
        ub, wd = -m[1]/xb, m[3]/z
        return (2j*np.pi/(x*m[1] + xb*m[0]))*(-2j*np.pi/z)*(ub - wd)**pu*wd**pw/(((x - z)*(ub - wd) - m[2])*(-(1 - z)*wd - m[4]))
    if x < z < 1:          # K^- on the pole of K^2, l^- on the pole of (p1-l)^2
        u1, w5 = m[0]/x, -m[4]/(1 - z)
        return (2j*np.pi/(xb*m[0] + x*m[1]))*(-2j*np.pi/(1 - z))*(u1 - w5)**pu*w5**pw/(((x - z)*(u1 - w5) - m[2])*(z*w5 - m[3]))
    return 0.0


m = np.array([1.0 - 0.40j, 1.0 - 0.50j, 0.7 - 0.30j, 0.45 - 0.60j, 0.45 - 0.35j])
for zz, where in ((0.12, '0 < z < x'), (0.60, 'x < z < 1'), (-0.30, 'z < 0'), (1.40, 'z > 1')):
    d_ = max(abs(V_numeric(0.3, zz, m, pu, pw) - V_residue(0.3, zz, m, pu, pw)) for pu, pw in ((0, 0), (1, 0), (0, 1)))
    check(f'Vggq, x = 0.3, z = {zz:+.2f} ({where}), numerators 1, q^-, l^-', d_, 2e-5)


def W_numeric(x, z, m):
    xb = 1 - x
    inner = lambda w: quad(lambda u: 1/(((x - z)*(u - w) - m[0])*(-xb*u - m[1])), -np.inf, np.inf, **QK)[0]
    return quad(lambda w: inner(w)/((-(1 - z)*w - m[2])*(z*w - m[3])), -np.inf, np.inf, **QK)[0]


def W_residue(x, z, m):
    xb = 1 - x
    if 0 < z < x:
        ub, wd = -m[1]/xb, m[3]/z
        return (-2j*np.pi/xb)*(-2j*np.pi/z)/(((x - z)*(ub - wd) - m[0])*(-(1 - z)*wd - m[2]))
    return 0.0


mw = np.array([0.7 - 0.30j, 1.0 - 0.50j, 0.45 - 0.60j, 0.45 - 0.35j])
for zz, where in ((0.12, '0 < z < x'), (0.60, 'z > x'), (-0.30, 'z < 0')):
    check(f'W1g, x = 0.3, z = {zz:+.2f} ({where})', abs(W_numeric(0.3, zz, mw) - W_residue(0.3, zz, mw)), 1e-6)
kk, qq, ll, zz = 1.0, 0.7, 0.45, 0.12; zeta = zz/0.3
D = (1 - zeta)*ll + zeta*qq - zeta*(1 - zeta)*kk + zeta*(1 - zeta)*kk/0.7          # (l_perp - zeta k)^2 + zeta(1-zeta) k^2/(1-x)
check('Vggq: residues/(2 pi)^2 = zeta/(k^2 l_perp^2 D)', (V_residue(0.3, zz, np.array([kk, kk, qq, ll, ll]), 0, 0)/(2*np.pi)**2).real - zeta/(kk*ll*D), 1e-12)
check('W1g: residues/(2 pi)^2 = -zeta/((1-x) l_perp^2 D)', (W_residue(0.3, zz, np.array([qq, kk, ll, ll]))/(2*np.pi)**2).real + zeta/(0.7*ll*D), 1e-12)

print('3  parameter integrals at finite eps against the Laurent expansions')
QR = dict(epsabs=1e-12, epsrel=1e-11, limit=300)


def W1g_A(x, ep):
    """-e^(eps gamma_E) Gamma(1+eps) (1-x)^eps Int_0^1 ds s^(-1-eps) Int_s^1 dzeta/zeta (2 - x zeta) V^(-eps), with s = e^(-y)"""
    xb = 1 - x
    def g(y):
        s = np.exp(-y); c0 = (1 - xb*s)**(-ep)
        a = quad(lambda ze: 2/ze*((1 - x*ze - xb*s)**(-ep) - c0), s, 1, **QR)[0]
        b = quad(lambda ze: -x*(1 - x*ze - xb*s)**(-ep), s, 1, **QR)[0]
        return a + b + 2*y*c0
    return -np.exp(EG*ep)*Gam(1 + ep)*xb**ep*quad(lambda y: np.exp(ep*y)*g(y), 0, np.inf, epsabs=1e-10, epsrel=1e-10, limit=400)[0]


def W1g_A_series(x, ep):                 # momentum space
    xb = 1 - x; l_ = np.log(xb)
    return -2/ep**2 - (2*l_ + x)/ep - Z2 - l_ - l_**2 + 2*Li2(x) - 2*x


def W1g_B(x, ep):                        # all orders in eps
    xb = 1 - x
    return -2*xb*np.exp(EG*ep)*Gam(1 + ep)*Beta(-ep, 2 - ep)*((xb**ep - xb)/(1 - ep) + (1 - xb**ep)/ep)


def W1g_B_direct(x, ep):
    xb = 1 - x
    inner = lambda ze: quad(lambda s: s**(-1 - ep)*((1 - x*ze - xb*s)**(-ep) + (1 - s)*(1 - x*ze - xb*s)**(-1 - ep)), 0, ze, **QR)[0]
    return -np.exp(EG*ep)*Gam(1 + ep)*xb**ep*2*x*xb*quad(inner, 0, 1, epsabs=1e-10, epsrel=1e-10, limit=300)[0]


def W1g_B_series(x, ep):
    xb = 1 - x; l_ = np.log(xb)
    return (2*x*xb - 2*xb*l_)/ep + 4*x*xb - xb*l_**2


def Vggq_region(x, ep):
    """region 0 < z < x of Vggq: coefficients of (eslash kslash, kslash eslash); eps < 0"""
    xb = 1 - x; Tv = np.array([x/2, -xb/2]); pre = np.exp(EG*ep)*xb**ep
    b0, b1 = Beta(1 - ep, 1 - ep), Beta(2 - ep, 1 - ep)
    Ra = Gam(ep)*x*np.array([-xb*(b0 - b1) - ep*xb*b1, -(xb/x)*(1 + x)*b0 + xb*b1 - ep*xb*b1])
    Vf = lambda ze, s: 1 - x*ze - xb*s
    i2 = quad(lambda ze: quad(lambda s: (s*Vf(ze, s))**(-ep), 0, ze, **QR)[0], 0, 1, **QR)[0]
    Rb2 = 2*(1 - ep)*Gam(ep)*x*Tv*i2
    def ecoef(ze, s, c_):
        e3 = s*x - 2*(1 - ep)*s*(1 - x*ze + (x - xb)*s)
        if c_ == 0: return (-(1 + xb) + x*ze*(xb - x))/(2*xb) + s*(2 - 3*x)/2 - e3/2
        return (2*x*ze - 1)/2 + s*(1 - 3*x)/2 - e3/2
    def inner(ze, c_):
        F = lambda s: Vf(ze, s)**(-1 - ep)*ecoef(ze, s, c_)
        F0 = F(0.0)
        return -ze**(-ep)*F0/ep + quad(lambda s: s**(-1 - ep)*(F(s) - F0), 0, ze, **QR)[0]
    Rb1 = Gam(1 + ep)*x*xb*np.array([quad(inner, 0, 1, args=(c_,), epsabs=1e-10, epsrel=1e-10, limit=300)[0] for c_ in (0, 1)])
    return pre*(Ra + Rb2 + Rb1)


def Vggq_region_series(x, ep):
    xb = 1 - x; l_ = np.log(xb)
    return np.array([-(1 + x)*l_/(2*ep) - (1 + x)*l_**2/4, -xb/ep - 2*xb + xb*l_/(2*ep) + xb*l_**2/4])


for x in (0.3, 0.7):
    d_ = [W1g_A(x, ep) - W1g_A_series(x, ep) for ep in (-0.04, -0.02, -0.01)]
    r_ = (2*d_[1] - d_[0])/(2*d_[2] - d_[1])       # one Richardson step removes the term linear in eps
    require(f'W1g, x = {x}: A approaches its expansion (differences {d_[0]:+.4f}, {d_[1]:+.4f}, {d_[2]:+.4f}; second order ratio {r_:.2f})',
            abs(d_[2]) < 0.05 and 3.5 < r_ < 4.5)
    d_ = [W1g_B(x, ep) - W1g_B_series(x, ep) for ep in (-0.04, -0.02, -0.01)]
    require(f'W1g, x = {x}: B to all orders approaches its expansion (ratios {d_[0]/d_[1]:.3f}, {d_[1]/d_[2]:.3f})',
            1.9 < d_[0]/d_[1] < 2.1 and 1.9 < d_[1]/d_[2] < 2.1)
    check(f'W1g, x = {x}: B to all orders against quadrature at eps = -0.3', W1g_B(x, -0.3) - W1g_B_direct(x, -0.3), 1e-8)
    d_ = [Vggq_region(x, ep) - Vggq_region_series(x, ep) for ep in (-0.04, -0.02, -0.01)]
    for c_, name in ((0, 'eslash kslash'), (1, 'kslash eslash')):
        require(f'Vggq, x = {x}, {name}: region 0 < z < x approaches its expansion (ratios {d_[0][c_]/d_[1][c_]:.3f}, {d_[1][c_]/d_[2][c_]:.3f})',
                abs(d_[2][c_]) < 0.05 and 1.8 < d_[0][c_]/d_[1][c_] < 2.1 and 1.8 < d_[1][c_]/d_[2][c_] < 2.1)
wA = wB = 0
for x in rng.uniform(0.05, 0.95, 5):
    xb = 1 - x; lx, lb = np.log(x), np.log(xb)
    r1 = np.array([-(1 + x)*lb, -2*xb + xb*lb])/2; f1 = np.array([-(1 + x)*lb**2/4, -2*xb + xb*lb**2/4])        # region at x: pole, finite part
    r2 = np.array([-(1 + xb)*lx, -2*x + x*lx])/2; f2 = np.array([-(1 + xb)*lx**2/4, -2*x + x*lx**2/4])         # region at 1 - x
    pole = r1 + np.array([-r2[1], -r2[0]]); fin = f1 + np.array([-f2[1], -f2[0]])                              # second region by charge conjugation
    wA = max(wA, abs(2*(pole[0] - pole[1]) - (2 - 2*(lx + lb))), abs(2*(fin[0] - fin[1]) - (4 - lx**2 - lb**2)))
    wB = max(wB, abs(-4*xb*pole[0] - 4*x*pole[1] - (-2*x*lx + 2*xb*lb)), abs(-4*xb*fin[0] - 4*x*fin[1] - (-x*lx**2 + xb*lb**2)))
check('Vggq: the two regions together give A', wA, 1e-12)
check('Vggq: the two regions together give B', wB, 1e-12)


def results(x):
    """(A, B) by Laurent order 'eta order, eps order', at b = b0"""
    y = 1 - x; lx, ly = np.log(x), np.log(y); L = lx + ly
    return {'Vggq': {'0,-1': (2 - 2*L, -2*x*lx + 2*y*ly), '0,0': (4 - lx**2 - ly**2, -x*lx**2 + y*ly**2)},
            'W1g': {'0,-2': (-2, 0), '0,-1': (-(2*ly + x), 2*x*y - 2*y*ly), '0,0': (-3*Z2 - ly - ly**2 + 2*Li2(x) - 2*x, 4*x*y - y*ly**2)},
            'W2g': {'0,-2': (-2, 0), '0,-1': (-(2*lx + y), -2*x*y + 2*x*lx), '0,0': (-3*Z2 - lx - lx**2 + 2*Li2(y) - 2*y, -4*x*y + x*lx**2)}}


print('4  results against the stored Laurent coefficients (five x, two helicity structures, six orders)')
worst = {}
for tag, x in {'1o10': 0.1, '3o10': 0.3, '1o2': 0.5, '7o10': 0.7, '9o10': 0.9}.items():
    st = json.load(open(repo/f'paper2/r0/theory/gqq_laurent_x{tag}.json'))['diagrams']
    for name, F in results(x).items():
        for key in ('-1,-2', '-1,-1', '-1,0', '0,-2', '0,-1', '0,0'):
            a, b = F.get(key, (0, 0))
            for h, pred in (('ek', a - b/(2*x)), ('ke', a + b/(2*(1 - x)))):
                worst[name] = max(worst.get(name, 0), abs(-st[name]['laurent'][h][key][0] - pred))     # the stored prefactor is -C_A/2
for name, v in worst.items(): check(name, v, 1e-8)

finish()
