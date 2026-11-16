"""Unequal-timelike-virtuality gg -> V1* V2* amplitudes, VVamp class A.

The published massless-quark evaluator uses alpha_s/(2*pi) per loop,
qT finite remainders and mu^2=s. No electroweak couplings are included.
This module does not modify the equal-mass implementation or its datasets.
"""
from fractions import Fraction
from itertools import product
from pathlib import Path
import sys

import mpmath as mp

HERE = Path(__file__).resolve().parent
ZZ = HERE.parent / "zz" / "scripts"
sys.path.insert(0, str(ZZ))
from ggvv import GGVV, Functions  # noqa: E402
from zz_points import pol, dot, PARTS, DIRECT  # noqa: E402
from zz_cache import cache_provenance as provenance, translated_inputs  # noqa: E402

SCHEMA = 1
LOG_HELICITIES = ("++++", "++--", "----", "--++")
ALL_HELICITIES = tuple("".join(h) for h in product("+-", "+-", "+0-", "+0-"))


def canonical(value):
    """Exact decimal/rational CLI inputs; avoid silently using binary floats."""
    if isinstance(value, float):
        raise ValueError("Pass decimal/rational strings, not binary float inputs")
    return str(Fraction(str(value)))


def real(value):
    f = Fraction(canonical(value))
    return mp.mpf(f.numerator) / f.denominator


def configuration(s, ma2, mb2):
    config = {"s": canonical(s), "ma2": canonical(ma2), "mb2": canonical(mb2)}
    ss, aa, bb = (real(config[k]) for k in ("s", "ma2", "mb2"))
    if aa <= 0 or bb <= 0 or ss <= (mp.sqrt(aa) + mp.sqrt(bb))**2:
        raise ValueError("Require positive timelike ma2, mb2 and s > (sqrt(ma2)+sqrt(mb2))^2")
    return config


def kinematics(config, k):
    """Forward t=t_plus-d, d=sqrt(ma2*mb2)*10^-k (mass-squared units).

    Stable product forms avoid subtracting tu-ma2*mb2 near the pinch.
    Only 0 < d < sqrt(lambda)/2 is scanned; exact endpoints are excluded.
    """
    if not isinstance(k, int) or isinstance(k, bool) or k < 1:
        raise ValueError("k must be a positive integer")
    c = configuration(**config)
    s, a, b = (real(c[name]) for name in ("s", "ma2", "mb2"))
    q2 = mp.sqrt(a*b)
    root = mp.sqrt((s-a-b)**2 - 4*a*b)
    tminus = -(s-a-b+root)/2
    tplus = a*b/tminus
    d = q2*mp.power(10, -k)
    if not 0 < d < root/2:
        raise ValueError("point is not on the forward half; increase k")
    t, u = tplus-d, tminus+d
    pT2 = d*(root-d)/s
    G = mp.log(s*pT2/(a*b))/2 + mp.j*mp.pi/2
    E = mp.sqrt(s)/2
    E3, E4 = (s+a-b)/(4*E), (s+b-a)/(4*E)
    km = root/(4*E)
    ct = 1 - 2*d/root
    st = mp.sqrt(pT2)/km
    return dict(s=s, ma2=a, mb2=b, root=root, tplus=tplus, tminus=tminus,
                d=d, t=t, u=u, pT2=pT2, q2=q2, G=G, E=E, E3=E3, E4=E4,
                kmod=km, ct=ct, st=st,
                x1=(s+a-b+root)/(2*s), x2=(s+a-b-root)/(2*s))


def predictions(config):
    """Conditional top two log coefficients / one-loop log in VVamp units.

    Assumes the same r1=r2=0 two-loop boundary as the Letter. No new
    two-variable finite matching function is assumed or evaluated.
    The predictions are independent of amplitude data and fit choices.
    """
    c = configuration(**config)
    s, a, b = (real(c[name]) for name in ("s", "ma2", "mb2"))
    L = mp.log(mp.sqrt(a*b)/s)
    return {"L": L, "2l_CF": {3: -mp.mpf(4)/3, 2: -2*L-3},
            "2l_CA": {3: mp.mpf(2)/3, 2: 2*L}}


def complete_coefficients(A, A12, A34, A1234):
    """Exact crossing map from published ggvvamp-1.0/ggvvampprec.cpp.

    A12=(s,u,a,b), A34=(s,u,b,a), A1234=(s,t,b,a).
    Distinguishing all four evaluations is essential for unequal masses.
    """
    out = dict(A)
    out[3], out[6], out[7] = A12[2], A12[5], A12[4]
    out[10], out[11] = -A34[8], -A34[9]
    out[12], out[13] = A12[9], A12[8]
    out[14], out[15] = -A1234[9], -A1234[8]
    out[19], out[20] = A12[18], A12[17]
    return out


def vectors(kin, helicity):
    if helicity not in ALL_HELICITIES:
        raise ValueError(f"invalid helicity: {helicity}")
    E, km, ct, st = (kin[n] for n in ("E", "kmod", "ct", "st"))
    theta = mp.atan2(st, ct)
    p1, p2 = [E, 0, 0, E], [E, 0, 0, -E]
    p3 = [kin["E3"], km*st, 0, km*ct]
    p4 = [kin["E4"], -km*st, 0, -km*ct]
    l1, l2, l3, l4 = ({"+": 1, "-": -1, "0": 0}[h] for h in helicity)
    e1, e2 = pol(0, 0, l1, E, E), pol(mp.pi, 0, l2, E, E)
    e3 = [mp.conj(z) for z in pol(theta, 0, l3, kin["E3"], km, mp.sqrt(kin["ma2"]))]
    e4 = [mp.conj(z) for z in pol(mp.pi-theta, mp.pi, l4, kin["E4"], km, mp.sqrt(kin["mb2"]))]
    return (p1, p2, p3, p4), (e1, e2, e3, e4)


def tensor_weights(kin, helicity):
    (p1, p2, p3, _), (e1, e2, e3, e4) = vectors(kin, helicity)
    e12, e13, e14, e23, e24, e34 = (dot(e1,e2),dot(e1,e3),dot(e1,e4),dot(e2,e3),dot(e2,e4),dot(e3,e4))
    a, b = dot(e1,p3), dot(e2,p3)
    c, d, e, f = dot(e3,p1), dot(e3,p2), dot(e4,p1), dot(e4,p2)
    return {1:e12*e34, 2:e13*e24, 3:e14*e23,
            4:e12*c*e, 5:e12*c*f, 6:e12*d*e, 7:e12*d*f,
            8:e13*b*e, 9:e13*b*f, 10:e14*b*c, 11:e14*b*d,
            12:e23*a*e, 13:e23*a*f, 14:e24*a*c, 15:e24*a*d,
            16:e34*a*b, 17:a*b*c*e, 18:a*b*c*f, 19:a*b*d*e, 20:a*b*d*f}


def contract(coefficients, kin, helicities):
    out = {}
    for h in helicities:
        T = tensor_weights(kin, h)
        out[h] = {part: mp.fsum(coefficients[part][j]*T[j] for j in range(1,21)) for part in PARTS}
    return out


def evaluate(g, kin, helicities=LOG_HELICITIES):
    """Four independent crossing configurations, sharing special-function cache."""
    s, t, u, a, b = (kin[n] for n in ("s", "t", "u", "ma2", "mb2"))
    F = Functions()
    calls = [(s,t,a,b), (s,u,a,b), (s,u,b,a), (s,t,b,a)]
    computed, timing = {}, []
    for args in calls:
        if args not in computed:
            computed[args] = g.compute(*args, functions=F)
            timing.append(list(g.timing))
    if F.order_bad:
        raise ArithmeticError("Li22 argument outside validated convergence region")
    raw = [computed[args] for args in calls]
    A = {part: complete_coefficients(*[{j:C[j][part] for j in DIRECT} for C in raw]) for part in PARTS}
    H = contract(A, kin, helicities)
    if any(not mp.isfinite(z) for pieces in H.values() for z in pieces.values()):
        raise ArithmeticError("nonfinite helicity amplitude")
    return H, A, {"distinct_crossings":len(computed), "li22_bad_order":F.order_bad,
                  "li22_max_rate":F.max_rate, "timing_function_algebra_seconds":timing}


def cpair(z, digits=None):
    n = digits or mp.mp.dps
    return [mp.nstr(mp.re(z), n), mp.nstr(mp.im(z), n)]


def complex_pair(pair):
    provenance.require_pair(pair)
    return mp.mpc(*[mp.mpf(v) for v in pair])


def sources():
    paths = {p.name:p for p in HERE.glob("*.py") if p.name in ("uv_core.py", "run_scan.py")}
    paths.update({"zz/"+n:ZZ/n for n in ("ggvv.py", "zz_points.py", "zz_cache.py", "translate_ggvvamp.py")})
    paths["dps/provenance.py"] = Path(provenance.__file__)
    return paths
