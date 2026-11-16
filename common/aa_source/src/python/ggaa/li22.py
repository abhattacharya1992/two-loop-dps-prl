"""The weight-four multiple polylogarithm

    Li_{2,2}(x, y) = sum_{i > j >= 1} x^i y^j / (i^2 j^2),

in the convention of Frellesvig, Tommasini and Wever, JHEP 03 (2016) 189, which is the
one used by the ancillary files of arXiv:1909.08993.

For max(|x|, |x y|) < 0.99 the defining series is summed directly.  Closer to the edge
of convergence, which happens near threshold, the integral representation

    Li_{2,2}(x, y) = -x * Integral_0^1 dt ln(t) Li_2(x y t) / (1 - x t)

is used; it follows from the definition by writing 1/i^2 = -Integral_0^1 t^(i-1) ln t dt.
"""

import mpmath as mp

SERIES_LIMIT = 0.99


def li22_series(x, y):
    x, y = mp.mpmathify(x), mp.mpmathify(y)
    if not (mp.isfinite(x) and mp.isfinite(y)):
        raise ValueError('Li22 arguments must be finite')
    if max(abs(x), abs(x * y)) >= 1:
        raise ValueError(f"Li22 series does not converge at ({x}, {y})")
    if x == 0 or y == 0:
        return mp.mpf(0)
    with mp.extradps(10):
        ax, aw = abs(x), abs(x * y)
        scaled_inner = mp.mpf(0)   # x^i sum_{j<i} y^j/j^2
        w_power = x * y            # (xy)^i
        total = mp.mpf(0)
        i = 1
        while True:
            total += scaled_inner / i**2
            scaled_inner = x * (scaled_inner + w_power / i**2)
            w_power *= x * y
            # After term i, keep S_i fixed in the rest of the outer sum and
            # separately bound all j>i additions. Both bounds use absolute
            # geometric tails; cancellation of S_i cannot hide the second one.
            tail = (abs(scaled_inner) / ((i + 1)**2 * (1 - ax)) +
                    ax * abs(w_power) / ((1 - ax) * (1 - aw) *
                                          (i + 1)**2 * (i + 2)**2))
            if tail <= mp.eps * abs(total) / 4:
                break
            i += 1
    return +total


def li22_integral(x, y):
    w = x * y
    integrand = lambda t: -mp.log(t) * mp.polylog(2, w * t) / (1 - x * t)
    return x * mp.quad(integrand, [0, 0.5, 0.9, 0.99, 1])


def li22(x, y):
    x, y = mp.mpmathify(x), mp.mpmathify(y)
    if not (mp.isfinite(x) and mp.isfinite(y)):
        raise ValueError('Li22 arguments must be finite')
    if max(abs(x), abs(x * y)) < SERIES_LIMIT:
        return li22_series(x, y)
    return li22_integral(x, y)
