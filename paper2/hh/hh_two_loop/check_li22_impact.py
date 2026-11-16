#!/usr/bin/env python3
# Paper App. A: where check_li22_values.py finds Li_{2,2} differences far above the working precision,
# which side is off, and whether the stored amplitude values are affected.
#
# For each point BETA:K (delta = 10^-K, 30 + 3K digits, as in the scan) and each expression:
#   1. the Li22new arguments are recorded as in check_li22_values.py;
#   2. for each argument, evaluate_hh.li22 at the working precision (the value the scan used), li22 at
#      twice the digits (a more precise value of the same representation) and the independent
#      reference of check_li22_values.py are compared; the three worst arguments are printed with the
#      error of each side relative to li22 at twice the digits;
#   3. the four colour terms at L = -1, 0, 1 are evaluated as the scan evaluates them, once with the
#      working-precision values (this must reproduce the stored point, which is checked) and once with
#      the values at twice the digits; the largest relative change of a colour term is the effect of the
#      Li_{2,2} evaluation on the stored amplitude values.
# Needs the ancillary files (make fetch in ..).
#
# On the points computed with the earlier evaluate_hh.li22 (split at t = 1/2 only) it found that the
# working-precision value was the one off, by up to 10^-29, and that the colour terms changed by up
# to 3 x 10^-16 (n_f/N of Gauge2Simp at beta = 1/2, delta = 10^-7); all points were
# recomputed with the present li22 on 3 November 2026.
#
#   python3 check_li22_impact.py 1/2:6 1/2:7 1/2:8
import sys, json, pathlib, time
from fractions import Fraction
import mpmath as mp
import evaluate_hh as ev
from scan_colour_scale import L_VALUES, working_digits, split_terms, identify, env
from check_li22_values import reference, arguments

HERE = pathlib.Path(__file__).resolve().parent


def stored(gauge, beta, k):
    """the stored point of the scan at this gauge, beta and delta, at the regular precision"""
    p = HERE/'points'/f'g{gauge}_b{beta.replace("/", "o")}_d{k}_p{working_digits(k)}.json'
    return json.loads(p.read_text())


def colour_terms(codes, names, x, y, li22):
    """A = i Gauge/4 for each colour term at L = -1, 0, 1, with Li22new = li22, as scan_colour_scale.run_job"""
    out = {}
    for L in L_VALUES:
        muR2 = mp.e**(-L)
        for name, c in zip(names, codes):
            e = env(x, y, muR2, 1, 1)
            e['Li22new'] = li22
            out[(name, L)] = mp.j*eval(c, {'__builtins__': {}}, e)/4
    return out


def main(specs):
    prepared = {}
    for n in (1, 2):
        codes = split_terms(n)
        prepared[n] = (codes, identify(codes))        # identify() works at 20 digits; the precision is reset below
    for spec in specs:
        beta, k = spec.split(':')
        k = int(k)
        for n, (codes, names) in prepared.items():
            t = time.time()
            dps = working_digits(k)
            mp.mp.dps = dps
            for f in (ev.li22, ev.LiClassical, ev.Log, ev.Zeta):
                f.cache_clear()
            fb = Fraction(beta)
            bv = mp.mpf(fb.numerator)/fb.denominator
            x = (1 - bv)/(1 + bv)
            y = x + mp.mpf(10)**(-k)
            pairs = arguments(codes, x, y)
            work, fine, rows = {}, {}, []
            for a, b in pairs:
                a, b = mp.re(a), mp.re(b)
                work[(a, b)] = ev.li22(a, b)                          # as in the scan
                with mp.workdps(2*dps):
                    ev.li22.cache_clear()
                    v2 = ev.li22(mp.mpf(a), mp.mpf(b))
                ev.li22.cache_clear()
                fine[(a, b)] = v2
                ref = reference(a, b)
                scale = max(abs(v2), mp.mpf(10)**(-dps))
                rows.append((max(abs(work[(a, b)] - v2), abs(ref - v2))/scale, a, b, v2,
                             abs(work[(a, b)] - v2)/scale, abs(ref - v2)/scale))
            rows.sort(key=lambda r: -r[0])
            print(f'Gauge{n}Simp, beta = {beta}, delta = 10^-{k}, {dps} digits, {len(pairs)} arguments; worst three '
                  f'(error of the working-precision li22 and of the reference, relative to li22 at {2*dps} digits):')
            for _, a, b, v2, e_work, e_ref in rows[:3]:
                print(f'    a = {mp.nstr(a, 12)}, b = {mp.nstr(b, 12)}, Li22 = {mp.nstr(v2, 6)}: '
                      f'li22 {mp.nstr(e_work, 3)}, reference {mp.nstr(e_ref, 3)}')
            # effect on the colour terms: the scan's values against those with the more precise Li22 values
            mp.mp.dps = dps
            A_work = colour_terms(codes, names, x, y, lambda a, b: work[(mp.re(a), mp.re(b))])
            A_fine = colour_terms(codes, names, x, y, lambda a, b: +fine[(mp.re(a), mp.re(b))])
            rec = stored(n, beta, k)
            # each colour term relative to its largest modulus over the three values of L (a term can vanish at one L)
            scale = {c: max(max(abs(A_fine[(c, L)]) for L in L_VALUES), mp.mpf(10)**(-dps)) for c in names}
            repro = max(abs(A_work[(c, L)] - mp.mpc(*map(mp.mpf, rec['A'][c][str(L)])))/scale[c] for (c, L) in A_work)
            change = {c: max(abs(A_work[(c, L)] - A_fine[(c, L)]) for L in L_VALUES)/scale[c] for c in names}
            print(f'    stored point reproduced to {mp.nstr(repro, 3)}; largest relative change of a colour term with '
                  f'the more precise Li22 values: ' + ', '.join(f'{c} {mp.nstr(v, 3)}' for c, v in change.items())
                  + f' ({time.time() - t:.0f} s)', flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
