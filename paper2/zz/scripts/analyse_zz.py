#!/usr/bin/env python3
"""Leading logarithms of the forward gg -> ZZ amplitudes and the quark-splitting prediction
of Paper 2 (Secs. 5 and 6).

Usage: analyse_zz.py POINTS_DIR [OUT_JSON]

At fixed beta the helicity amplitudes on the forward line are fitted to
    sum_{j<=J} f_j G^j + p_T^2 sum_{l<=P} h_l G^l
through J + P + 2 consecutive points; window 0 ends at the smallest delta, window 1 is shifted
one point up, and the change between the two is quoted as the uncertainty.

Normalisation of the source: a factor alpha_s/(2 pi) per loop, i.e. M^(n) = 2^n x (VVamp), with
the Paper 2 expansion in a_s = alpha_s/(4 pi). With DC = 2 C_F - C_A, B = 6 C_F,
G_q = G + d_q, d_q = -((C_A - C_F)/DC) L, L = ln[x(1 - x)], x = (1 + beta)/2, eq. (5.x) of
Paper 2 gives, in VVamp units,
  one-loop c1 G:      M^(2) = (c1/2) [ -(4/3) DC G_q^3 - B G_q^2 ] + O(G):
      C_F part:  f3/c1 = -4/3,  f2/c1 = -2 L - 3
      C_A part:  f3/c1 = +2/3,  f2/c1 = 2 L
  one-loop constant K (|J_z| changes by 2):  M^(2) = (K/2) [ -4 DC G_q^2 + (4 DC - 2 B) G_q ] + O(1):
      C_F part:  f2/K = -4,  f1/K = -4 L - 2
      C_A part:  f2/K = +2,  f1/K =  4 L - 2
"""
import json
import pathlib
import sys
from fractions import Fraction

import mpmath as mp

mp.mp.dps = 50
PARTS = ("1l", "2l_CF", "2l_CA")
PDEG = 3          # degree in G of the p_T^2 correction in the two-loop fits
HEL = {"+": 1, "-": -1, "0": 0}


def cplx(p):
    return mp.mpc(mp.mpf(p[0]), mp.mpf(p[1]))


def jz_change(h):
    l1, l2, l3, l4 = (HEL[c] for c in h)
    return (l1 - l2) - (l3 - l4)


def load(d):
    pts = [json.loads(p.read_text()) for p in sorted(pathlib.Path(d).glob("zz_*.json"))]
    regular, control = {}, []
    for r in pts:
        if r["dps"] != 30 + 5 * r["delta_exp"]:
            control.append(r)
        else:
            regular.setdefault(r["beta"], []).append(r)
    for v in regular.values():
        v.sort(key=lambda r: r["delta_exp"])
    return regular, control


def fit(records, vals, deg, pdeg, window=0):
    n = deg + pdeg + 2
    end = len(records) - window
    if end - n < 0:
        return None
    M = mp.matrix(n, n)
    v = mp.matrix(n, 1)
    for i in range(n):
        r = records[end - n + i]
        G, p = cplx(r["G"]), mp.mpf(r["pT2"])
        for j in range(deg + 1):
            M[i, j] = G ** j
        for k in range(pdeg + 1):
            M[i, deg + 1 + k] = p * G ** k
        v[i] = vals[end - n + i]
    sol = mp.lu_solve(M, v)
    return [sol[j] for j in range(deg + 1)]


def bnum(beta):
    b = Fraction(beta)
    return mp.mpf(b.numerator) / b.denominator


def c2(z):
    """Complex number as two strings with 30 significant digits (no float rounding)."""
    return [mp.nstr(mp.re(z), 30), mp.nstr(mp.im(z), 30)]


def analyse(regular):
    out = {}
    for beta, recs in sorted(regular.items(), key=lambda kv: float(Fraction(kv[0]))):
        x = (1 + bnum(beta)) / 2
        L = mp.log(x * (1 - x))
        res = {"n_points": len(recs), "deltas": [r["delta_exp"] for r in recs], "x": float(x), "L": float(L),
               "hel": {}}
        for h in recs[0]["M"]:
            vals = {p: [cplx(r["M"][h][p]) for r in recs] for p in PARTS}
            one = fit(recs, vals["1l"], 1, 2)
            if one is None:
                continue
            one_w = fit(recs, vals["1l"], 1, 2, window=1)
            e = {"dJz": jz_change(h), "c1": c2(one[1]), "K": c2(one[0]),
                 "c1_w1": c2(one_w[1]) if one_w else None, "K_w1": c2(one_w[0]) if one_w else None}
            for part in ("2l_CF", "2l_CA"):
                # main model: cubic in G plus p_T^2 times a cubic in G (8 points);
                # window 1 (shifted one point up) and the model with a quadratic as checks
                f = fit(recs, vals[part], 3, PDEG)
                if f is None:
                    continue
                fw = fit(recs, vals[part], 3, PDEG, window=1)
                fp = fit(recs, vals[part], 3, PDEG - 1)
                e[part] = {"f": [c2(c) for c in f], "f_w1": [c2(c) for c in fw] if fw else None,
                           "f_lower_model": [c2(c) for c in fp] if fp else None}
            res["hel"][h] = e
        out[beta] = res
    return out


def controls(regular, control):
    rows = []
    for r in control:
        base = [q for q in regular.get(r["beta"], []) if q["delta_exp"] == r["delta_exp"]]
        if not base:
            continue
        q = base[0]
        worst = {}
        for h in r["M"]:
            for p in PARTS:
                a, b = cplx(r["M"][h][p]), cplx(q["M"][h][p])
                if abs(a) > 0:
                    worst[p] = max(worst.get(p, 0.0), float(abs(a - b) / abs(a)))
        rows.append({"beta": r["beta"], "delta_exp": r["delta_exp"], "dps": r["dps"], "base_dps": q["dps"],
                     "max_rel_diff": worst})
    return rows


def C(v):
    return mp.mpc(mp.mpf(v[0]), mp.mpf(v[1])) if v is not None else None


def report(res, ctrl):
    print("precision controls (same point with 15 more digits):", ctrl)
    for beta, r in res.items():
        L, x = mp.mpf(r["L"]), mp.mpf(r["x"])
        print(f"\n=== beta = {beta}  (x = {r['x']:.6f}, L = ln[x(1-x)] = {r['L']:.8f}; points k = {r['deltas']})")
        hel = r["hel"]
        if not hel:
            print("  not enough points yet")
            continue
        cmax = max(abs(C(e["c1"])) for e in hel.values())
        print("one loop, J_z conserved (dJz = 0) and |dJz| = 2:")
        for h, e in sorted(hel.items(), key=lambda kv: -abs(C(kv[1]["c1"]))):
            if abs(e["dJz"]) in (0, 2):
                c1, K = C(e["c1"]), C(e["K"])
                sh = abs(C(e["c1"]) - C(e["c1_w1"])) if e["c1_w1"] else None
                print(f"  {h} dJz={e['dJz']:+d}  c1 = {mp.nstr(c1, 12)}  |c1|/(8 pi x^2) = {mp.nstr(abs(c1) / (8 * mp.pi * x**2), 10)}"
                      f"  |c1|/(8 pi (1-x)^2) = {mp.nstr(abs(c1) / (8 * mp.pi * (1 - x)**2), 10)}  K = {mp.nstr(K, 12)}"
                      f"  |K|/(pi x(1-x)) = {mp.nstr(abs(K) / (mp.pi * x * (1 - x)), 10)}  (c1 window shift {mp.nstr(sh, 2) if sh is not None else '-'})")
        pred_log = {"2l_CF": (-mp.mpf(4) / 3, -2 * L - 3), "2l_CA": (mp.mpf(2) / 3, 2 * L)}
        pred_K = {"2l_CF": (mp.mpf(0), mp.mpf(-4), -4 * L - 2), "2l_CA": (mp.mpf(0), mp.mpf(2), 4 * L - 2)}
        print("two loops, amplitudes with a one-loop logarithm: f3/c1, f2/c1 [prediction]")
        for h, e in hel.items():
            c1 = C(e["c1"])
            if e["dJz"] != 0 or abs(c1) < 1e-3 * cmax or "2l_CF" not in e:
                continue
            for part in ("2l_CF", "2l_CA"):
                f = [C(v) for v in e[part]["f"]]
                fw = [C(v) for v in e[part]["f_w1"]] if e[part]["f_w1"] else None
                p3, p2 = pred_log[part]
                d3 = abs(f[3] - fw[3]) / abs(c1) if fw else None
                d2 = abs(f[2] - fw[2]) / abs(c1) if fw else None
                print(f"  {h} {part}: f3/c1 = {mp.nstr(f[3] / c1, 14)} [{mp.nstr(p3, 14)}]  "
                      f"f2/c1 = {mp.nstr(f[2] / c1, 14)} [{mp.nstr(p2, 14)}]  window shifts {mp.nstr(d3, 2)}, {mp.nstr(d2, 2)}")
        print("two loops, |dJz| = 2 with a one-loop constant K: f3/K, f2/K, f1/K [prediction]")
        for h, e in hel.items():
            K = C(e["K"])
            if abs(e["dJz"]) != 2 or abs(K) < 1e-6 or "2l_CF" not in e:
                continue
            for part in ("2l_CF", "2l_CA"):
                f = [C(v) for v in e[part]["f"]]
                fw = [C(v) for v in e[part]["f_w1"]] if e[part]["f_w1"] else None
                p3, p2, p1 = pred_K[part]
                d = [abs(f[j] - fw[j]) / abs(K) for j in range(4)] if fw else None
                print(f"  {h} {part}: f3/K = {mp.nstr(f[3] / K, 6)} [0]  f2/K = {mp.nstr(f[2] / K, 14)} [{mp.nstr(p2, 14)}]  "
                      f"f1/K = {mp.nstr(f[1] / K, 14)} [{mp.nstr(p1, 14)}]  window shifts {[mp.nstr(v, 2) for v in d] if d else '-'}")
        print("two loops, dJz = 0 without a one-loop logarithm: f3/K, f2/K (no prediction of a DPS term)")
        for h, e in hel.items():
            c1, K = C(e["c1"]), C(e["K"])
            if e["dJz"] != 0 or abs(c1) >= 1e-3 * cmax or abs(K) < 1e-6 or "2l_CF" not in e:
                continue
            for part in ("2l_CF", "2l_CA"):
                f = [C(v) for v in e[part]["f"]]
                print(f"  {h} {part}: f3/K = {mp.nstr(f[3] / K, 6)}  f2/K = {mp.nstr(f[2] / K, 10)}  f1/K = {mp.nstr(f[1] / K, 10)}")


def main():
    regular, control = load(sys.argv[1])
    res = analyse(regular)
    ctrl = controls(regular, control)
    if len(sys.argv) > 2:
        pathlib.Path(sys.argv[2]).write_text(json.dumps({"betas": res, "control": ctrl}, indent=1))
    report(res, ctrl)


if __name__ == "__main__":
    main()
