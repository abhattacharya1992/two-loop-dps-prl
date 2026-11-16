"""Evaluate the gg -> V1 V2 amplitude coefficients A_j^[A] of von Manteuffel and Tancredi
(arXiv:1503.08835) at any precision, from the translated ggvvamp code (translate_ggvvamp.py).

Conventions of the source (qT scheme, mu^2 = s, a factor alpha_s/(2 pi) per loop, class A):
    s = (p1+p2)^2, t = (p1-p3)^2, u = (p2-p3)^2, ma2 = p3^2, mb2 = p4^2.
The two-loop coefficients are returned split into their C_F and C_A parts.
"""
import json
import marshal
import pathlib
import time

import mpmath as mp
from zz_cache import cache_provenance, translated_inputs

RSIZE = 493625


class Functions:
    """The special functions of p0002, real parts as in ggvvamp (GiNaC numerics, real part).

    li22(a, b) = sum_{i>j>=1} a^i b^j / (i^2 j^2): first argument with the outer sum, as in
    eq. (5.14) of arXiv:1503.04812. `order` records which ordering satisfies the convergence
    conditions |a| < 1, |a b| < 1 of that equation, as a check of the convention."""

    def __init__(self, swap_li22=False):
        self.swap = swap_li22
        self.cache = {}
        self.order_ok = 0
        self.order_bad = 0
        self.max_rate = 0.0

    def li22(self, a, b):
        if self.swap:
            a, b = b, a
        key = ("li22", mp.mp.prec, a, b)
        if key not in self.cache:
            rate = max(abs(a), abs(a * b))
            if rate < 1:
                self.order_ok += 1
            else:
                self.order_bad += 1
            self.max_rate = max(self.max_rate, float(rate))
            self.cache[key] = mp.re(li22(a, b))
        return self.cache[key]

    def lin(self, n, x):
        key = ("li", mp.mp.prec, n, x)
        if key not in self.cache:
            self.cache[key] = mp.re(mp.polylog(n, x))
        return self.cache[key]

    def log(self, x):
        return mp.log(abs(x))          # real part of the principal logarithm

    def zeta(self, n):
        return mp.zeta(n)


def li22_series(x, y):
    """sum_{i>j>=1} x^i y^j/(i^2 j^2) for max(|x|, |x y|) < 1."""
    rate = max(abs(x), abs(x * y))
    tol = mp.eps * (1 - rate) / 16
    inner = mp.mpf(0)
    xp, yp = mp.mpf(1), mp.mpf(1)
    total = mp.mpf(0)
    i = 1
    while True:
        xp *= x
        term = xp * inner / (i * i)
        total += term
        yp *= y
        inner += yp / (i * i)
        if i > 3 and abs(term) <= tol * abs(total) and abs(xp * inner) <= tol * abs(total) * i * i:
            return total
        i += 1


def li22_integral(x, y):
    """-x Integral_0^1 ln(t) Li2(x y t)/(1 - x t) dt, valid for x, x y not on [1, inf)."""
    w = x * y
    f = lambda t: -mp.log(t) * mp.polylog(2, w * t) / (1 - x * t)
    return x * mp.quad(f, [0, mp.mpf(1) / 2, mp.mpf(9) / 10, mp.mpf(99) / 100, 1])


def li22(x, y):
    x, y = mp.mpmathify(x), mp.mpmathify(y)
    if max(abs(x), abs(x * y)) < mp.mpf("0.95"):
        return li22_series(x, y)
    return li22_integral(x, y)


class GGVV:
    def __init__(self, build):
        build = pathlib.Path(build).resolve()
        self.build = build
        self.input_hashes = cache_provenance.file_hashes(translated_inputs(build))
        ns = {}
        exec(compile((build / "p0002.py").read_text(), "p0002", "exec"), ns)
        self.p0002 = ns["p0002"]
        ns = {}
        for code in marshal.loads((build / "algebra.marshal").read_bytes()):
            exec(code, ns)
        self.algebra = [ns[f"p{k:04d}"] for k in range(3, 2464)]
        rat = json.loads((build / "rationals.json").read_text())
        self.rationals = sorted(((v, k) for k, v in rat.items()))
        self.table = json.loads((build / "coefficients.json").read_text())
        if self.input_hashes != cache_provenance.file_hashes(translated_inputs(build)):
            raise cache_provenance.CacheMismatch("Translated ZZ inputs changed while GGVV was loading")

    def variables(self, s, t, ma2, mb2):
        u = ma2 + mb2 - s - t
        root = mp.sqrt((s - ma2 - mb2) ** 2 - 4 * ma2 * mb2)
        t0 = (ma2 + mb2 - s - root) / 2
        t1 = (ma2 + mb2 - s + root) / 2
        if not (t0 <= t <= t1):
            raise ValueError("t outside the physical region")
        x = -t0 / ma2
        y = t1 * (u + t) / (ma2 * mb2) - 1
        z = t1 * t / (ma2 * mb2)
        return u, x, y, z

    def compute(self, s, t, ma2, mb2, functions=None):
        """Return {j: {"1l": A, "2l_CF": A, "2l_CA": A}} (complex) for the nine computed j."""
        s, t, ma2, mb2 = (mp.mpmathify(v) for v in (s, t, ma2, mb2))
        u, x, y, z = self.variables(s, t, ma2, mb2)
        F = functions or Functions()
        Q = [mp.mpf(int(k.split("/")[0])) / int(k.split("/")[1]) for _, k in self.rationals]
        r = [None] * (RSIZE + 1)
        r[1], r[2], r[3], r[4], r[5], r[6], r[7] = ma2, t / ma2, u / ma2, mb2 / ma2, x, y, z
        t0 = time.time()
        self.p0002(r, Q, F)
        t1 = time.time()
        for f in self.algebra:
            f(r, Q)
        t2 = time.time()
        self.timing = (t1 - t0, t2 - t1)
        out = {}
        for key, entry in self.table.items():
            j, loop, reim = map(int, key.split(","))
            for part, idx in entry["terms"].items():
                v = r[idx] * (mp.pi if entry["pi"] else 1)
                label = "1l" if loop == 0 else ("2l_" + part)
                d = out.setdefault(j, {})
                d[label] = d.get(label, mp.mpc(0)) + (v if reim == 0 else mp.mpc(0, 1) * v)
        self.functions = F
        return out
