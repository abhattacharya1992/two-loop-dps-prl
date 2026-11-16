"""Formulas of the Supplement that enter r0 besides the collinear function: the factor F_eps of Eq. (eq:Feps), the
soft function of Eq. (eq:softeta) at nu_r = 1 and b = b0, and the daughter self-energy of Eq. (eq:sigma)."""
import numpy as np
from scipy.special import gamma, digamma

EG, Z2 = np.euler_gamma, np.pi**2/6


def F_eps(e):
    return np.exp(-EG*e)*gamma(1 - 2*e)/(gamma(1 + e)*gamma(1 - e))


def soft(e):
    """S/Delta_C = 2 e^(-eps gamma_E) Gamma(-eps) [2/eta - gamma_E - psi(-eps)]: (coefficient of 1/eta, coefficient of eta^0)"""
    pre = 2*np.exp(-EG*e)*gamma(-e)
    return 2*pre, pre*(-EG - digamma(-e))


def sigma(e, x, daughter):
    """self-energy of a daughter with momentum fraction x: the coefficient of C_F for a quark,
    the coefficients of (C_A, T_F n_f) for a gluon"""
    common = F_eps(e)*gamma(e)*gamma(1 - e)**2/gamma(2 - 2*e)*(1 - x)**e
    if daughter == 'quark': return -(1 - e)*common
    return (5 - 3*e)/(3 - 2*e)*common, -4*(1 - e)/(3 - 2*e)*common


def laurent(f, orders=(-2, -1, 0), radius=0.05, n=64):
    """Laurent coefficients of f(eps) from a discrete Fourier sum on a circle around eps = 0"""
    e = radius*np.exp(2j*np.pi*(np.arange(n) + 0.5)/n); v = f(e)
    return {j: float(np.mean(v*e**(-j)).real) for j in orders}
