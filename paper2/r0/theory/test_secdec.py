import numpy as np, sympy as sp, mpmath as mp
from secdec import decompose, Sector, laurent
from lc_engine import eps
from factors import eta
from numeric_checks import close
zz, tt = sp.symbols('s t')
def total(facs, vs, const=1):
    secs = [Sector(j, m, r, vs, const) for j, m, r, vm in decompose(facs, vs)]
    return (lambda ep, et: sum(S.value(ep, et) for S in secs)), len(secs)
# test 1: int z^(-1-eta-eps)(1-z)^(1-eps) = B(-eta-eps, 2-eps)
f, ns = total({zz: -1 - eta - eps, 1 - zz: 1 - eps}, [zz])
c = laurent(f)
def check_coefficients(observed, expected):
    if not observed:
        raise ValueError('No Laurent coefficients were computed')
    for key in observed.keys() | expected.keys():
        close(complex(observed.get(key,0)),complex(expected.get(key,0)),
              atol=1e-8,rtol=0,label=f'Laurent coefficient {key}')
check_coefficients(c,{(0,-1):-1,(0,0):-1})
print('test1 sectors', ns, {k: complex(v) for k, v in c.items() if abs(v) > 1e-12})
print('   exact: -1/eps - 1')
# test 2: pure eta pole
f, ns = total({zz: -1 - eta}, [zz])
c=laurent(f)
check_coefficients(c,{(-1,0):-1})
print('test2', {k: complex(v) for k, v in c.items() if abs(v) > 1e-12})
# test 3: corner  int dz dt z^(-1-eta) t^(-eps) (z+t)^(-1+2eps)
f, ns = total({zz: -1 - eta, tt: -eps, zz + tt: -1 + 2*eps}, [zz, tt])
c = laurent(f)
check_coefficients(c,{(-1,-1):-1,(0,-2):-1,(0,0):mp.pi**2/4})
print('test3 sectors', ns, {k: complex(v) for k, v in c.items() if abs(v) > 1e-12})
