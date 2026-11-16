"""Wilson line couplings in J1(g -> q qbar) with colour factor C_F: diagrams W1q and W2qb.
The gluon from the Wilson line of an operator attaches to the daughter created by the same operator.

Units k^2 = mu^2 = 1, a_s = alpha_s/(4 pi). There is no rapidity divergence. Normalised to the tree splitting,
in impact parameter space at b = b0,
    W1q  = C_F G(eps) (1-x)^eps,      W2qb = C_F G(eps) x^eps,
    G(eps) = -2 F_eps Gamma(eps) Gamma(-eps) Gamma(2-eps)/Gamma(2-2eps)
           = (2/eps^2) e^(-eps gamma_E) Gamma(1-eps) (1-eps)/(1-2eps) = 2/eps^2 + 2/eps + 4 + zeta2 + O(eps),
with F_eps = e^(-eps gamma_E) Gamma(1-2eps)/(Gamma(1+eps) Gamma(1-eps)) the factor of Eq. (eq:Feps) of the Supplement.
The powers (1-x)^eps and x^eps come from the squared momenta -k^2/(1-x) of the quark and -k^2/x of the antiquark
at the DPS residue.

Steps checked here:
  1  the two forms of G(eps), the Feynman parameter integral behind it, and its Laurent expansion
  2  the result against the stored Laurent coefficients
Usage: python3 check_qq_W1q_W2qb.py <repository root>"""
import json, sys, pathlib
import numpy as np
from scipy.special import gamma, beta
from scipy.integrate import quad
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from supplement_formulas import F_eps, laurent, EG, Z2
from report import check, finish
repo = pathlib.Path(sys.argv[1])
G = lambda e: (2/e**2)*np.exp(-EG*e)*gamma(1 - e)*(1 - e)/(1 - 2*e)
G_momentum = lambda e: -2*gamma(e)*gamma(-e)*gamma(2 - e)/gamma(2 - 2*e)       # before the Fourier transform


print('1  the function G(eps)')
e0 = 0.3 + 0.2j
check('G = F_eps times the momentum space factor, at eps = 0.3 + 0.2i', abs(G(e0) - F_eps(e0)*G_momentum(e0)), 1e-12)
e1 = -0.37
check('Int_0^1 dt t^(-1-eps) (1-t)^(1-eps) = B(-eps, 2-eps) at eps = -0.37',
      quad(lambda t: t**(-1 - e1)*(1 - t)**(1 - e1), 0, 1, epsabs=1e-13, epsrel=1e-13)[0] - beta(-e1, 2 - e1), 1e-8)
c = laurent(G)
check('Laurent expansion 2/eps^2 + 2/eps + 4 + zeta2', max(abs(c[-2] - 2), abs(c[-1] - 2), abs(c[0] - 4 - Z2)), 1e-11)

print('2  result against the stored Laurent coefficients (five x, two helicity structures, six orders)')
worst = {'W1q': 0.0, 'W2qb': 0.0}
for tag, x in {'1o10': 0.1, '3o10': 0.3, '1o2': 0.5, '7o10': 0.7, '9o10': 0.9}.items():
    st = json.load(open(repo/f'paper2/r0/theory/gqq_laurent_x{tag}.json'))['diagrams']
    for name, arg, sign in (('W2qb', x, 1), ('W1q', 1 - x, -1)):          # the stored prefactor of W1q is -C_F
        c = laurent(lambda e: G(e)*arg**e)
        for h in ('ek', 'ke'):
            for j in (-2, -1, 0):
                worst[name] = max(worst[name], abs(sign*st[name]['laurent'][h][f'0,{j}'][0] - c[j]),
                                  abs(st[name]['laurent'][h][f'-1,{j}'][0]))                # no pole in eta
for name, v in worst.items(): check(name, v, 1e-9)
finish()
