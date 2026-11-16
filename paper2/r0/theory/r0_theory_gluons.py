"""r0 for gluon splittings from the regions, against the closed form found in the gg -> HH data, for the C_A and the
T_F n_f parts and the two components nu = 0, 2 (combine_gg.kernel), and the precision checks at x = 1/2 (higher
quadrature orders, 16 points on the eps circle).  Reads gg_laurent_x*.json, writes r0_theory_gluons.json.

    python3 r0_theory_gluons.py
"""
import json, os
import combine_gg as C

PAIRS = [(1, 10), (3, 10), (1, 2)]


def main():
    out = {'x': {}}
    print(f'{"x":>5} {"col":>3} {"nu":>2} {"eps^-2":>12} {"exp":>7} {"eps^-1":>13} {"exp":>13} {"r0 (regions)":>14} '
          f'{"r0 (data)":>14} {"diff":>9} {"1/eta":>8}')
    for num, den in PAIRS:
        if not (os.path.exists(f'gg_laurent_x{num}o{den}.json') and os.path.exists(f'gg_laurent_x{den-num}o{den}.json')):
            print(f'x = {num}/{den}: not computed')
            continue
        xx = num/den
        rec = {}
        for col in ('CA', 'nf'):
            res = C.kernel(num, den, col)
            ep = C.expected_poles(xx, col)
            rec[col] = {}
            for nu, r in res.items():
                data = C.r0_closed(xx, col, nu)
                rec[col][f'nu={nu}'] = dict(r, **{'expected poles': {'eps^-2': ep[0], 'eps^-1': ep[1]},
                                                  'r0 closed form (data)': data,
                                                  'difference': r['eps^0'] - data})
                print(f'{xx:5.2f} {col:>3} {nu:2d} {r["eps^-2"]:12.8f} {ep[0]:7.3f} {r["eps^-1"]:13.8f} {ep[1]:13.8f}'
                      f' {r["eps^0"]:14.9f} {data:14.9f} {r["eps^0"] - data:9.1e} {r["max |1/eta coefficient|"]:8.1e}')
        out['x'][f'{num}/{den}'] = rec
    # precision at x = 1/2: (a) higher quadrature orders, (b) 16 instead of 32 points on the eps circle (aliasing)
    checks = {}
    for label, suffix in (('quadrature orders 96, 72, 40 instead of 64, 48, 28', '_nq96-72-40'),
                          ('16 instead of 32 points on the eps circle', '_neps16')):
        if os.path.exists('gg_laurent_x1o2.json') and os.path.exists(f'gg_laurent_x1o2{suffix}.json'):
            ch = {}
            for col in ('CA', 'nf'):
                base, new = C.kernel(1, 2, col), C.kernel(1, 2, col, suffix)
                for nu in (0, 2):
                    ch[f'{col}, nu={nu}'] = new[nu]['eps^0'] - base[nu]['eps^0']
            checks[label] = ch
            print(f'x = 1/2, {label}: change of r0', {k: f'{v:+.1e}' for k, v in ch.items()})
    out['precision checks at x = 1/2 (change of r0)'] = checks
    json.dump(out, open('r0_theory_gluons.json', 'w'), indent=1)
    print('wrote r0_theory_gluons.json')


if __name__ == '__main__':
    main()
