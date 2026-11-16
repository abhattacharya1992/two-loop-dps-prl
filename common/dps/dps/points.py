"""Point records: the stored forward points of the three processes in one form.

A record is a dict with

    process      'hh', 'aa' or 'zz'
    component    the helicity projection: 1 or 2 (tensor T_1, T_2) for hh and aa, the helicity
                 string ('++++', '+-++', ...) for zz
    beta         velocity as a fraction string
    delta_exp    k, with delta = 10^-k
    dps          working precision of the evaluation
    G, pT2       mpmath numbers
    values       {name: mpc}: the amplitude parts at this point (see each loader)
    regular      True for the points of the scan, False for the precision-control points

Native files are read as they are: hh/hh_two_loop/points (the colour terms at three
renormalisation scales, turned into the coefficients of the powers of ln(m^2/mu^2)), zz/points
(36 helicity amplitudes for the parts 1l, 2l_CF, 2l_CA, in VVamp units) and dps/points_aa (this
package's gg -> AA points, colour monomials at mu = m).  Loaders return
({component_or_(component, beta): [records sorted by decreasing delta]}, [control records]).
"""
import json
import pathlib

import mpmath as mp


def cplx(pair):
    if not isinstance(pair,(list,tuple)) or len(pair)!=2:
        raise ValueError('Complex point value must be a real/imaginary pair')
    value = mp.mpc(mp.mpf(pair[0]), mp.mpf(pair[1]))
    if not mp.isfinite(value):
        raise ValueError('Nonfinite complex point value')
    return value


def positive_pT2(value):
    value = mp.mpf(value)
    if not mp.isfinite(value) or value <= 0:
        raise ValueError('pT2 must be positive and finite')
    return value


def _finish(regular):
    for v in regular.values():
        v.sort(key=lambda r: r['delta_exp'])
    return regular


def load_hh(points_dir):
    """gg -> HH: values 'COLOUR:L0', 'COLOUR:L1', 'COLOUR:L2' for the four colour monomials, the
    coefficients of ln(m_H^2/mu_R^2)^n of A = i Gauge/(4 M2); regular points have dps = 30 + 3k."""
    regular, control = {}, []
    for p in sorted(pathlib.Path(points_dir).glob('g*_b*_d*_p*.json')):
        r = json.loads(p.read_text())
        rec = {'process': 'hh', 'component': r['gauge'], 'beta': r['beta'], 'delta_exp': r['delta_exp'],
               'dps': r['dps'], 'G': cplx(r['G']), 'pT2': positive_pT2(r['pT2']), 'values': {},
               'regular': r['dps'] == 30 + 3 * r['delta_exp'], 'file': p.name}
        for colour, byL in r['A'].items():
            f = {int(L): cplx(v) for L, v in byL.items()}
            rec['values'][f'{colour}:L0'] = f[0]
            rec['values'][f'{colour}:L1'] = (f[1] - f[-1]) / 2
            rec['values'][f'{colour}:L2'] = (f[1] + f[-1]) / 2 - f[0]
        if rec['regular']:
            regular.setdefault((r['gauge'], r['beta']), []).append(rec)
        else:
            control.append(rec)
    return _finish(regular), control


def load_zz(points_dir):
    """gg -> ZZ: one record per helicity amplitude with values '1l', '2l_CF', '2l_CA' (VVamp units:
    alpha_s/(2 pi) per loop); regular points have dps = 30 + 5k."""
    regular, control = {}, []
    for p in sorted(pathlib.Path(points_dir).glob('zz_b*_d*_p*.json')):
        r = json.loads(p.read_text())
        for hel, parts in r['M'].items():
            rec = {'process': 'zz', 'component': hel, 'beta': r['beta'], 'delta_exp': r['delta_exp'],
                   'dps': r['dps'], 'G': cplx(r['G']), 'pT2': positive_pT2(r['pT2']),
                   'values': {part: cplx(v) for part, v in parts.items()},
                   'regular': r['dps'] == 30 + 5 * r['delta_exp'], 'file': p.name}
            if rec['regular']:
                regular.setdefault((hel, r['beta']), []).append(rec)
            else:
                control.append(rec)
    return _finish(regular), control


def load_aa(points_dir):
    """gg -> AA: values are the colour monomials ('N2', 'N_Nf', 'Nf2', 'Nf_over_N', ...) of the
    two-loop remainder in the normalisation of P1 (file/64, mu_R = m_A = 1) and, under '1l:',
    those of the one-loop remainder ('1l:N', '1l:Nf'); the record says itself whether it belongs
    to the scan ('regular')."""
    regular, control = {}, []
    for p in sorted(pathlib.Path(points_dir).glob('aa_g*_b*_d*_p*.json')):
        r = json.loads(p.read_text())
        rec = {'process': 'aa', 'component': r['gauge'], 'beta': r['beta'], 'delta_exp': r['delta_exp'],
               'dps': r['dps'], 'G': cplx(r['G']), 'pT2': positive_pT2(r['pT2']),
               'values': {name: cplx(v) for name, v in r['A'].items()},
               'regular': r.get('regular', True), 'file': p.name}
        for name, v in r.get('A1', {}).items():
            rec['values'][f'1l:{name}'] = cplx(v)
        if rec['regular']:
            regular.setdefault((r['gauge'], r['beta']), []).append(rec)
        else:
            control.append(rec)
    return _finish(regular), control


LOADERS = {'hh': load_hh, 'zz': load_zz, 'aa': load_aa}


def compare_stored(points_dir, name_a, name_b):
    """Largest relative difference between the amplitude parts of two point files of the same
    point, at the full precision of the stored digits (the loaders above parse the values at the
    working precision of the analysis, 50 digits, at which a control point evaluated with more
    digits coincides with the point it repeats)."""
    a, b = (json.loads((pathlib.Path(points_dir)/n).read_text()) for n in (name_a, name_b))

    def leaves(d, prefix=()):          # the [re, im] pairs of the amplitude parts, however they are nested
        if isinstance(d, dict):
            for k, v in d.items():
                yield from leaves(v, prefix + (k,))
        else:
            yield prefix, d

    with mp.workdps(max(int(a['dps']), int(b['dps'])) + 20):
        worst = mp.mpf(0)
        if ('A1' in a) != ('A1' in b):
            raise ValueError('Different one-loop component sets in precision control')
        va_all = {**{('A',) + k: v for k, v in leaves(a['A'])},
                  **{('A1',) + k: v for k, v in leaves(a.get('A1', {}))}}
        vb = {**{('A',) + k: v for k, v in leaves(b['A'])},
              **{('A1',) + k: v for k, v in leaves(b.get('A1', {}))}}
        if not va_all or va_all.keys() != vb.keys():
            raise ValueError('Different or empty amplitude component sets in precision control')
        for k, va in va_all.items():
            za, zb = cplx(va), cplx(vb[k])
            if abs(zb) > 0:
                worst = max(worst, abs(za - zb) / abs(zb))
            elif za != 0:
                return mp.inf   # nonzero against exact zero has no finite relative error
        return worst


def series(records, name):
    """The values of one named part along a list of records."""
    return [r['values'][name] for r in records]
