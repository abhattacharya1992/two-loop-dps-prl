"""The forward line at fixed velocity, in units of the mass m of the two produced particles.

With the velocity beta of the produced particles in the partonic centre-of-mass frame,

    xk = (1 - beta)/(1 + beta),   s = (1 + xk)^2/xk,   -t = y = xk + delta,   -u = z = xk + 1/xk - y,

so that delta -> 0 is the forward limit, t -> -xk m^2.  The variables of the forward fits are

    p_T^2 = (t u - m^4)/s = (y z - 1) xk/(1 + xk)^2,
    G     = (1/2) ln[(m^4 - t u)/m^4] = (1/2) ln(1 - y z)      (principal branch: 1 - y z < 0, Im G = pi/2),

and the momentum fraction of the DPS configuration is x = (1 + beta)/2 = 1/(1 + xk).  The same
kinematics is written in the three native forms: (xk, y) for the gg -> HH files, (s, t) for VVamp,
and (beta, cos theta) for the gg -> AA files, with cos theta = [1 - 2(1 + y)/s]/beta.

Everything is computed at the current mpmath precision; `point(beta, k)` takes beta as a string
fraction and delta = 10^-k.
"""
from fractions import Fraction

import mpmath as mp


def beta_number(beta):
    """mpf of a velocity given as a fraction string such as '1/2'."""
    b = Fraction(beta)
    return mp.mpf(b.numerator) / b.denominator


def momentum_fraction(beta):
    """x = (1 + beta)/2 of the DPS configuration."""
    return (1 + beta_number(beta)) / 2


def point(beta, k):
    """All kinematic variables of the point delta = 10^-k on the forward line at velocity beta."""
    b = beta_number(beta)
    xk = (1 - b) / (1 + b)
    delta = mp.mpf(10) ** (-k)
    y = xk + delta
    z = xk + 1 / xk - y
    s = (1 + xk) ** 2 / xk
    return dict(beta=b, xk=xk, delta=delta, y=y, z=z, s=s, t=-y, u=-z,
                pT2=(y * z - 1) * xk / (1 + xk) ** 2,
                G=mp.log(1 - y * z) / 2,
                costh=(1 - 2 * (1 + y) / s) / b)


def from_costh(beta, costh):
    """The same variables from (beta, cos theta), the variables of the gg -> AA files of P1."""
    b = beta_number(beta) if isinstance(beta, str) else mp.mpmathify(beta)
    c = mp.mpmathify(costh)
    xk = (1 - b) / (1 + b)
    y = (1 + b ** 2 - 2 * b * c) / (1 - b ** 2)
    z = (1 + b ** 2 + 2 * b * c) / (1 - b ** 2)
    s = (1 + xk) ** 2 / xk
    return dict(beta=b, xk=xk, delta=y - xk, y=y, z=z, s=s, t=-y, u=-z,
                pT2=(y * z - 1) * xk / (1 + xk) ** 2,
                G=mp.log(1 - y * z) / 2, costh=c)
