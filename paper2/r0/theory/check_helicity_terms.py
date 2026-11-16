"""Where the helicity-dependent term (2/3)(C_A - 2 T_F n_f)/F_nu of r0 for gluon splittings comes from.

The finite part of J1 (structures S1, S2, S3, T4 of gg_driver.py, in units of the tree coefficient of S1) is split
into its spin-1 part, the three linear structures plus the trace part (S1 + S2 + S3)/4 of T4, and its spin-3 part,
T4 minus that trace part: the helicity flip, in which both daughters have the helicity opposite to the parent's.
The spin-1 part is the tree t = S1 - S2/(1-x) - S3/x times a number plus

    R = (c/6) (S1 - S2 - S3),     c = 1 for the C_A part, -2 for the T_F n_f part,

which is x(1-x) c/3 times the tree in the helicity configuration where both daughters keep the helicity of the
parent and vanishes in the other three: the rational term of the one-loop splitting amplitude, Bern, Del Duca,
Kilgore, Schmidt, Phys. Rev. D 60 (1999) 116001, eq. (4.33).  In b space the spin-1 part keeps its form and the
spin-3 part changes sign (combine_gg.to_b).  The script checks the decomposition at x = 1/10, 3/10, 1/2 and that
R gives (2/3) c/F_0 to nu = 0 and nothing to nu = 2, the spin-3 part nothing to nu = 0 and (2/3) c/F_2 to nu = 2.

    python3 check_helicity_terms.py        (writes helicity_terms.json)
"""
import json
import combine_gg as C

STRUCT = C.STRUCT


def F(xx, nu):
    return 1/(xx*(1 - xx)) + 2*xx*(1 - xx) if nu == 0 else (xx**4 + (1 - xx)**4)/(xx*(1 - xx))


def parts(J, tau, c):
    """k-space coefficients of the tree-like part, of R and of the spin-3 part; deviation of the rest from the tree"""
    a = {s: J[s].get('0,0', 0) for s in STRUCT}
    spin1 = {'S1': a['S1'] + a['T4']/4, 'S2': a['S2'] + a['T4']/4, 'S3': a['S3'] + a['T4']/4, 'T4': 0}
    spin3 = {'S1': -a['T4']/4, 'S2': -a['T4']/4, 'S3': -a['T4']/4, 'T4': a['T4']}
    R = {'S1': c/6, 'S2': -c/6, 'S3': -c/6, 'T4': 0}
    rest = {s: spin1[s] - R[s] for s in STRUCT}
    lam = rest['S1']/tau['S1']
    tree = {s: lam*tau[s] for s in STRUCT}
    dev = max(abs(rest[s] - tree[s]) for s in STRUCT)
    return tree, R, spin3, dev


def helicity_check(xx):
    """R/t in the four helicity configurations (parent +; daughters ++, +-, -+, --), from the structures"""
    out = []
    for S1, S2, S3, T4 in ((0, 1, 1, 0.5), (1, 1, 0, 0.5), (1, 0, 1, 0.5), (0, 0, 0, 1)):
        t = S1 - S2/(1 - xx) - S3/xx
        r = (S1 - S2 - S3)/6
        out.append(r/t if t != 0 else (None if r == 0 else float('inf')))
    return out


worst = 0.0
record = {}
for num, den in ((1, 10), (3, 10), (1, 2)):
    xx = num/den
    h = helicity_check(xx)
    assert abs(h[0] - xx*(1 - xx)/3) < 1e-15 and h[1] == 0 and h[2] == 0 and h[3] is None, h
    for col, c in (('CA', 1.0), ('nf', -2.0)):
        J1, tau1 = C.load(f'{num}o{den}', col)
        J2, tau2 = C.load(f'{den - num}o{den}', col)
        p1, p2 = parts(J1, tau1, c), parts(J2, tau2, c)
        dev = max(p1[3], p2[3])
        line = [f'x = {num}/{den} {col:2s}  spin-1 part minus R, deviation from the tree {dev:.1e}']
        worst = max(worst, dev)
        for nu in (0, 2):
            contrib = {name: C.ratio(nu, C.to_b(p1[k]), tau1, C.to_b(p2[k]), tau2).real
                       for k, name in ((1, 'R'), (2, 'spin3'))}
            target = (2/3)*c/F(xx, nu)
            want = {'R': target, 'spin3': 0.0} if nu == 0 else {'R': 0.0, 'spin3': target}
            worst = max(worst, *(abs(contrib[k] - want[k]) for k in want))
            record.setdefault(f'{num}/{den}', {})[f'{col}, nu={nu}'] = {
                'R': contrib['R'], 'spin 3': contrib['spin3'], '(2/3) c/F_nu': target}
            line.append(f'nu={nu}: R {contrib["R"]:+.9f}, spin 3 {contrib["spin3"]:+.9f}, (2/3)c/F {target:+.9f}')
        print('   '.join(line))
print(f'largest deviation {worst:.1e}')
assert worst < 1e-6
json.dump({'x': record, 'largest deviation': worst}, open('helicity_terms.json', 'w'), indent=1)
print('passed: the nu = 0 term comes from R, the nu = 2 term from the helicity flip')
