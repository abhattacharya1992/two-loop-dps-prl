#!/usr/bin/env python3
# Writes numbers.tex (every number the text quotes, as a macro) and tables/table_*.tex
# (Tables 1-11) from the result files, and checks what the text says about the figures.
# Nothing in the tables or in the quoted numbers is typed by hand: results come from the
# analyses in ../hh, ../zz, ../n4, ../kernel, ../r0 and the shared code dps/ (../dps or ../../common/dps), the scan settings from the scan points, and the few numbers
# taken from P1 or from an earlier evaluation from small data files that say where each one
# comes from. make_figure.py reads the result files listed in SOURCES below.
#
#   python3 make_tables.py          (or `make tables`; `make figures` in ../hh runs it too)
import json, math, pathlib
from fractions import Fraction
import mpmath as mp
import sympy as sp

HERE = pathlib.Path(__file__).resolve().parent
R = HERE.parent
DPS = R/'dps' if (R/'dps').is_dir() else R.parent/'common'/'dps'     # the shared code: in paper2/ or in <master>/common/
SOURCES = {
    'colour': R/'hh'/'hh_two_loop'/'results_colour.json',      # Sec. 3.3
    'one_loop': R/'hh'/'hh_two_loop'/'results_one_loop.json',  # Sec. 3.4
    'scale': R/'hh'/'hh_two_loop'/'results_scale.json',        # Sec. 3.5
    'forward': R/'hh'/'hh_two_loop'/'results_forward.json',    # Sec. 2.5
    'stage2a': R/'hh'/'stage2a'/'results_stage2a.json',        # Sec. 4
    'scan': R/'hh'/'hh_two_loop'/'points',                     # the scan settings, read from the points
    'p1': HERE/'p1_numbers.json',                              # numbers quoted from P1, with their sources
    'bose': R/'hh'/'hh_two_loop'/'bose_symmetry.json',         # Bose-symmetry check of the evaluator
    'hh_li22': R/'hh'/'hh_two_loop'/'logs'/'check_li22_values.log',  # App. A: Li_{2,2} against eq. (li22sum)
    'hh_eval': R/'hh'/'hh_two_loop'/'evaluate_hh.py',          # App. A: the guard digits of its Li_{2,2}
    'zz_li22_scan': R/'zz'/'logs'/'check_li22_scan.log',       # App. A: gg -> ZZ Li_{2,2} at the precision of the scan
    'zz': R/'zz'/'results_zz.json',                            # Sec. 6.2, the fits of the gg -> ZZ amplitudes
    'zz_summary': R/'zz'/'summary.json',                       # Sec. 6.2, App. A: scan settings and cross-checks
    'windows': R/'hh'/'hh_two_loop'/'results_windows.json',    # Fig. 5: gg -> HH fits through every window and model
    'zz_windows': R/'zz'/'results_windows.json',               # Fig. 5: the same for gg -> ZZ
    'zz_forward': R/'zz'/'results_forward.json',               # Sec. 6.2: one-loop gg -> ZZ approach to the forward limit
    'aa_curve': R/'hh'/'one_loop'/'results_aa_curve.json',     # Fig. 2: exact one-loop gg -> AA from 90 degrees
    'zz_wide': R/'zz'/'results_wide.json',                     # Fig. 2: one-loop gg -> ZZ from 90 degrees
    'n4': R/'n4'/'results_stage_d.json',                       # Sec. 5.1, Tables 5, 6: planar N = 4 (Dixon-Esterlis)
    'kernel': R/'kernel'/'results_stage_b.json',               # App. E: the boundary term at one loop
    'r0_fits': R/'r0'/'data'/'r0_closed_forms.json',           # Sec. 7.1: the closed forms of r0 against the fits
    'r0_scan_hh': R/'r0'/'data'/'r0_from_data.json',           # Sec. 7.1: r0 of gg -> HH at the three scan velocities, with bounds
    'r0_scan_zz': R/'r0'/'data'/'r0_from_zz.json',             # Sec. 7.1, Table r0_zz: r0 of the eight gg -> ZZ amplitudes
    'r0_blind_hh': R/'r0'/'data'/'blind_test_hh.json',         # Sec. 7.1: blind test of the HH closed form
    'r0_blind_zz': R/'r0'/'data'/'blind_test_zz.json',         # Sec. 7.1: blind test of the ZZ closed form
    'r0_quarks': R/'r0'/'theory'/'r0_theory_quarks.json',      # Sec. 7.2, App. F, Table 9: g -> q qbar from the regions
    'r0_gluons': R/'r0'/'theory'/'r0_theory_gluons.json',      # Sec. 7.2, App. F, Table 9: g -> gg from the regions
    'n4reg': R/'r0'/'theory'/'kappa1_n4.json',                 # Sec. 7.3, Table 11: kappa_1 of planar N = 4 from the regions
    'n4fdh': R/'r0'/'theory'/'fdh_check_x1o2.json',             # Sec. 7.3: the scheme check of the scalar graphs
    'glauber': R/'r0'/'theory'/'glauber.json',                 # Sec. 7.2, 7.3, App. F: the imaginary part of the kernel
    'three_loop': R/'r0'/'theory'/'three_loop.json',           # Paper II Sec. 7.2: the three-loop prediction (equation and table)
    'dps_aa': DPS/'results_aa.json',                           # Sec. 7.1, Table 10: two loop gg -> AA on the forward line
    'dps_hh': DPS/'results_hh.json',                           # Table 10: the gg -> HH fits of the shared code
    'dps_aa_check': DPS/'aa_digit_check.json',                 # App. A: digit loss of the gg -> AA files
    'n4_primary': R/'n4'/'results_primary_nmhv.json',          # Sec. 5.1: five-loop next-to-MHV coefficients from the multi-Regge amplitude
    'ultrasoft': R.parent/'checks'/'ultrasoft_numerators.json',  # App. E: ultrasoft attachments with the complete numerators
    'ultrasoft_d': R.parent/'checks'/'ultrasoft_dimensions_20261112',   # App. E: the same in d dimensions (a folder of records)
    'uv_matching': R/'unequal_virtuality'/'results'/'matching_checks.json',  # separate unequal-virtuality matching tests
    'collins_soper': R.parent/'checks'/'collins_soper_kernel.json',     # App. E, Sec. 5.1: the Collins-Soper equation of the kernel
}


NUMBER_WORDS = ['no', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten',
                'eleven', 'twelve', 'thirteen', 'fourteen', 'fifteen', 'sixteen', 'seventeen', 'eighteen']


def load(key):
    return json.loads(SOURCES[key].read_text())


def sci(x, digits=1, floor=None):
    """1.2e-11 -> 1.2\\times10^{-11}; values below `floor` are written as the floor."""
    x = float(x)
    if floor is not None and x < floor:
        x = floor
    if x == 0:
        return '0'
    e = math.floor(math.log10(abs(x)))
    m = round(x/10**e, digits)
    if m >= 10:
        m /= 10; e += 1
    ms = f'{m:.{digits}f}'                 # keep trailing zeros: 2.0e-12 is quoted as 2.0, like 2.3
    return f'{ms}\\times10^{{{e}}}' if ms != '1' else f'10^{{{e}}}'


def hh_li22_check():
    """App. A: the log of hh/hh_two_loop/check_li22_values.py, which compares the Li_{2,2} of evaluate_hh.py
    with eq. (li22sum) at twice the digits for every argument of the scan, its precision control and the
    blind test: the largest relative difference, and the largest in units of 10^-digits of the point"""
    import re
    text = SOURCES['hh_li22'].read_text()
    rows = re.findall(r'^Gauge(\d)Simp, beta = (\S+), delta = 10\^-(\d+), (\d+) digits: (\d+) distinct arguments, '
                      r'nearest singular point of the integrand at t = 1 \+ (\S+), largest relative difference (\S+) ',
                      text, flags=re.M)
    points = {(b, int(k), int(d)) for _, b, k, d, *_ in rows}
    # every point of the scan (three velocities, delta = 10^-4 ... 10^-11), its control and the blind test
    scan = {(b, k, 30 + 3*k) for b in ('1/2', '1/4', '2/3', '1/10', '4/5') for k in range(4, 12)} | {('1/2', 10, 75)}
    assert points == scan and len(rows) == 2*len(scan), sorted(points ^ scan)
    total, = re.findall(r'^largest relative difference over all points: (\S+)$', text, flags=re.M)
    worst = max(float(r[6]) for r in rows)
    assert abs(worst - float(total)) <= 1e-3*worst
    ratio = max(float(r[6])*10**int(r[3]) for r in rows)
    guard, = re.findall(r'^LI22_GUARD = (\d+)', SOURCES['hh_eval'].read_text(), flags=re.M)
    return {'max_rel_diff': worst, 'max_in_working_precision': ratio, 'guard': int(guard),
            'nearest': min(float(r[5]) for r in rows)}


def zz_li22_scan_check():
    """App. A: the log of zz/scripts/check_li22_scan.py, which compares at the working precision of every point of
    the gg -> ZZ scan, its control and the blind test each real argument on the integral branch of ggvv.li22 with
    the integral split as in hh/hh_two_loop/evaluate_hh.li22: the number of arguments and the largest absolute
    difference in units of the last digit (the values are at most of order one)"""
    import re
    text = SOURCES['zz_li22_scan'].read_text()
    rows = re.findall(r'^beta = (\S+), delta = 10\^-(\d+), (\d+) digits: \d+ distinct arguments, (\d+) real ones on the '
                      r'integral branch(, \(\d+ complex)?.*; largest absolute difference (\S+) x 10\^-(\d+) ', text, flags=re.M)
    points = {(b, int(k), int(d)) for b, k, d, *_ in rows}
    expected = {(b, k, 30 + 5*k) for b in ('1/4', '1/2', '2/3', '4/5') for k in range(4, 13)} | {('1/2', 10, 95)}
    assert points == expected and len(rows) == len(expected), sorted(points ^ expected)
    assert not any(r[4] for r in rows), 'complex arguments on the integral branch were not compared'
    assert all(r[2] == r[6] for r in rows)
    units = max(float(r[5]) for r in rows)
    total, = re.findall(r'largest absolute difference in units of 10\^-digits: (\S+)$', text, flags=re.M)
    assert abs(units - float(total)) <= 1e-2*units
    return {'arguments': sum(int(r[3]) for r in rows), 'units': units}


def up(x):
    """a deviation rounded up to one significant digit, so the text quotes an upper bound"""
    x = float(x)
    e = math.floor(math.log10(x))
    m = math.ceil(round(x/10**e, 9))       # the inner round stops 2.0000000001 becoming 3
    if m == 10:
        m, e = 1, e+1
    return f'{m}\\times10^{{{e}}}' if m != 1 else f'10^{{{e}}}'


def frac_tex(q):
    """Fraction(1, 2) -> \\frac12, Fraction(3, 10) -> \\frac{3}{10}, Fraction(2) -> 2"""
    q = Fraction(q)
    if q.denominator == 1:
        return str(q.numerator)
    if q.numerator < 10 and q.denominator < 10:
        return f'\\frac{q.numerator}{q.denominator}'
    return f'\\frac{{{q.numerator}}}{{{q.denominator}}}'


def coeff_tex(expr):
    """'-4*C_A/3' -> -\\frac43C_A, '2' -> +2: an entry of Table 4, sign always written"""
    CA = sp.Symbol('C_A')
    e = sp.sympify(expr, locals={'C_A': CA})
    c, rest = e.as_coeff_Mul()
    assert rest in (1, CA), expr
    c = Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q))
    return ('-' if c < 0 else '+') + frac_tex(abs(c)) + ('C_A' if rest == CA else '')


def summarise_points(points_dir):
    """Scan settings from the point files: velocities, delta = 10^-k, and the digits rule."""
    pts = [json.loads(p.read_text()) for p in sorted(pathlib.Path(points_dir).glob('*.json'))]
    betas = sorted({p['beta'] for p in pts}, key=Fraction)
    ks = sorted({p['delta_exp'] for p in pts})
    regular = {k: min(p['dps'] for p in pts if p['delta_exp'] == k) for k in ks}
    per_decade = (regular[ks[-1]]-regular[ks[0]])//(ks[-1]-ks[0])
    base = regular[ks[0]]-per_decade*ks[0]
    assert all(regular[k] == base+per_decade*k for k in ks), 'working digits are not linear in k'
    controls = [{'gauge': p['gauge'], 'beta': p['beta'], 'delta_exp': p['delta_exp'],
                 'extra_digits': p['dps']-regular[p['delta_exp']]}
                for p in pts if p['dps'] != regular[p['delta_exp']]]
    L_values = sorted(int(L) for L in pts[0]['A']['N2'])
    pT2 = [float(p['pT2']) for p in pts]                 # p_T^2/m_H^2 of each point
    return {'points': len(pts), 'betas': betas, 'delta_exponents': ks, 'L_values': L_values,
            'digits': {'base': base, 'per_decade': per_decade}, 'controls': controls,
            'pT2_range': [min(pT2), max(pT2)]}


def preferred_beta(betas):
    """the velocity of Tables 2 and 7 and of Fig. 5: 1/2 if it was scanned"""
    return '1/2' if '1/2' in betas else sorted(betas, key=Fraction)[0]


def scan_summary():
    src = SOURCES['scan']
    return summarise_points(src) if src.is_dir() else json.loads(src.read_text())


# ---------------------------------------------------------------- gg -> ZZ (section 6.2, Table 7)
ZZ_LOG = ('++++', '++--', '----', '--++')      # amplitudes with a one-loop logarithm C_1 G
ZZ_K = ('+-++', '+---', '-+++', '-+--')        # |Delta J_z| = 2, one-loop forward limit K
ZZ_NOLOG = ('++00', '--00', '+-+-', '-+-+')    # Delta J_z = 0, forward limit without logarithm
LAM = sp.Symbol('L_x')                         # L_x = ln[x(1-x)]


def ratio_spread(curve, from_exp):
    """(M - K)/(K p_T^2) at the points with delta <= 10^-from_exp of an approach curve (Secs. 2.5, 6.2):
    the largest difference between two of its values relative to their mean, and the number of
    points"""
    r = [mp.mpc(mp.mpf(c['rel_re']), mp.mpf(c['rel_im']))/mp.mpf(c['pT2']) for c in curve
         if c['delta_exp'] >= from_exp]
    return max(abs(a-b) for a in r for b in r)/abs(sum(r)/len(r)), len(r)


def digits(spread):
    """'constant to n digits': the number of leading digits the spread leaves unchanged"""
    return int(mp.floor(-mp.log10(spread)))


def zz_predictions():
    """Eq. (5.6) with Delta_C = 2C_F - C_A, B = 6C_F and G -> G_q = G - (C_A - C_F)/Delta_C L_x
    (eq. 6.1), expanded in G: {(kind, power, part): coefficient}, linear in L_x, in units of
    C_1 (kind 'log') or K (kind 'K'), for the coefficient of C_F or C_A."""
    CF, CA, G = sp.symbols('C_F C_A G')
    DC, B = 2*CF-CA, 6*CF
    Gq = G-(CA-CF)/DC*LAM
    forms = {'log': (-sp.Rational(4, 3)*DC*Gq**3-B*Gq**2, (3, 2)),
             'K': (-4*DC*Gq**2+(4*DC-2*B)*Gq, (3, 2, 1))}
    out = {}
    for kind, (expr, powers) in forms.items():
        poly = sp.Poly(sp.expand(expr), G)
        for j in powers:
            c = sp.expand(sp.cancel(poly.coeff_monomial(G**j)))
            a, b = sp.expand(c.coeff(CF)), sp.expand(c.coeff(CA))
            assert sp.simplify(c-a*CF-b*CA) == 0, 'not linear in the Casimirs'
            out[(kind, j, 'CF')], out[(kind, j, 'CA')] = a, b
    return out


def lin_coeffs(expr):
    """a + b*L_x -> (Fraction(a), Fraction(b))"""
    p = sp.Poly(expr, LAM)
    return tuple(Fraction(str(p.coeff_monomial(m))) for m in (1, LAM))


def mpq(q):
    return mp.mpf(q.numerator)/q.denominator


def lin_tex(expr):
    """-4*L_x - 6 -> -4L_x-6, -8/3 -> -\\frac83"""
    a, b = lin_coeffs(expr)
    s = ''
    if b:
        s = ('-' if b < 0 else '')+(frac_tex(abs(b)) if abs(b) != 1 else '')+'L_x'
    if a or not s:
        s += ('-' if a < 0 else ('+' if s else ''))+frac_tex(abs(a))
    return s


def zz_part(zz, zs, zf):
    """Macros of section 6.2 and appendices A and B, and tables/table_zz.tex. The source
    expands in alpha_s/(2 pi) per loop, the paper in a_s = alpha_s/(4 pi): the n-loop
    coefficients are multiplied by 2^n, so the two-loop/one-loop ratios by 2."""
    mp.mp.dps = 30
    R = lambda v: mp.mpc(mp.mpf(v[0]), mp.mpf(v[1]))
    pred = zz_predictions()
    betas = sorted(zz['betas'], key=Fraction)
    scan = zs['scan']
    assert betas == scan['betas']
    rows = {}            # (kind, j, part) -> {'dev', 'shift', 'value' at the table velocity}
    tb = preferred_beta(betas)
    kcubic = nolog = 0
    for b in betas:
        r = zz['betas'][b]
        x = (1+mpq(Fraction(b)))/2
        lam = mp.log(x*(1-x))
        assert abs(lam-mp.mpf(r['L'])) < 1e-12
        for kind, hels, key, powers in (('log', ZZ_LOG, 'c1', (3, 2)), ('K', ZZ_K, 'K', (3, 2, 1))):
            for h in hels:
                one = R(r['hel'][h][key])
                for part in ('CF', 'CA'):
                    fit = r['hel'][h]['2l_'+part]
                    f, fw = [R(v) for v in fit['f']], [R(v) for v in fit['f_w1']]
                    for j in powers:
                        ratio = 2*f[j]/one
                        a0, a1 = lin_coeffs(pred[(kind, j, part)])
                        p = mpq(a0)+mpq(a1)*lam
                        row = rows.setdefault((kind, j, part), {'dev': 0, 'shift': 0})
                        row['dev'] = max(row['dev'], abs(ratio-p))
                        row['shift'] = max(row['shift'], 2*abs(f[j]-fw[j])/abs(one))
                        if b == tb and h == hels[0]:
                            row['value'], row['pred'] = ratio, p
                        if kind == 'K' and j == 3:
                            kcubic = max(kcubic, abs(ratio))
        for h in ZZ_NOLOG:
            K = R(r['hel'][h]['K'])
            for part in ('CF', 'CA'):
                f = [R(v) for v in r['hel'][h]['2l_'+part]['f']]
                nolog = max(nolog, *(abs(2*f[j]/K) for j in (1, 2, 3)))
    for key, row in rows.items():
        assert row['dev'] <= row['shift'], f'{key}: deviation larger than the fit-window shift'
        assert abs(mp.im(row['value'])) < 1e-9
    dev = lambda kind: max(row['dev'] for (k, j, p), row in rows.items() if k == kind and pred[(k, j, p)] != 0)
    shift = max(row['shift'] for row in rows.values())
    ctrl, = scan['controls']
    chk, fwd = zs['checks'], zs['route_b_forward']['CF']
    cf, = fwd
    ca_rows = sorted(zs['route_b_forward']['CA_plus4'], key=lambda r: r['delta_exp'])
    ca_matches=zs['route_b_forward']['CA']
    assert len(ca_rows)>=3 and len(ca_matches)>=2
    assert all(row['max_rel_diff']<1e-15 for row in ca_matches)
    ca_cf,=[row for row in ca_matches if (row['beta'],row['delta_exp'])==(cf['beta'],cf['delta_exp'])]
    assert ca_cf['dps']==cf['dps'] and ca_cf['amplitudes']==cf['amplitudes']
    odd = zs['one_loop']['odd_dJz_effective_power']
    assert 0.6 < odd[0] and odd[1] < 1.01, 'the text says that the odd-Delta J_z amplitudes vanish like p_T'
    assert zs['one_loop']['no_log_c1_max_rel'] < 1e-15, 'the text says that these amplitudes have no logarithm'
    cxx, liv, p2 = chk['cxx'], chk['li22_values'], chk['p0002']
    assert p2['differ'] == 0 and p2['statements_py'] == p2['statements']
    bench, = [p for p in cxx['points'] if p.get('point') == 'benchmark']
    fwd, = [p for p in cxx['points'] if p.get('point') != 'benchmark']
    lfwd = [p for p in liv['points'] if p.get('point') != 'benchmark']
    assert len(lfwd) == len(liv['points']) - 1 == len(betas)
    assert sorted((p['beta'] for p in lfwd), key=Fraction) == betas
    assert {p['delta_exp'] for p in lfwd} == {max(scan['delta_exponents'])}
    li_where = (f'the benchmark point and at $\\delta=10^{{-{lfwd[0]["delta_exp"]}}}$, the smallest value of the '
                f'scan, for all {NUMBER_WORDS[len(betas)]} velocities')
    alist = lambda names: ', '.join(f'$A_{{{a[1:]}}}$' for a in names[:-1]) + \
        (' and ' if len(names) > 1 else '') + f'$A_{{{names[-1][1:]}}}$'
    bench_cf = chk['benchmark_two_loop_CF']
    li = chk['li22']
    # Sec. 6.2: (M - K)/(K p_T^2) of the one-loop +-++ amplitude is constant to n digits over the
    # points with delta <= 10^-6, as in Sec. 2.5 for gg -> HH
    fwd_spreads = [ratio_spread(v['curve'], v['spread_from_delta_exp']) for v in zf.values()]
    assert len({n for _, n in fwd_spreads}) == 1
    fwd_digits, fwd_points = digits(max(s for s, _ in fwd_spreads)), fwd_spreads[0][1]
    assert li['convergent_as_given'] == li['calls']        # benchmark point (logs of make checks)
    zls = zz_li22_scan_check()
    assert 1 < zls['units'] < 1000, 'the text says "all but at most the last few digits"'
    assert scan['li22']['bad_order'] == 0 and scan['li22']['max_rate'] < 1   # every scan point
    macros = {
        'zzBetas': ','.join(frac_tex(Fraction(b)) for b in betas),
        'zzDeltaFirst': f'10^{{-{min(scan["delta_exponents"])}}}',
        'zzDeltaLast': f'10^{{-{max(scan["delta_exponents"])}}}',
        'zzDigitsBase': str(scan['digits']['base']), 'zzDigitsPerDecade': str(scan['digits']['per_decade']),
        'zzCtrlDelta': f'10^{{-{ctrl["delta_exp"]}}}', 'zzCtrlBeta': frac_tex(Fraction(ctrl['beta'])),
        'zzPrecisionCtrl': up(ctrl['max_rel_diff']),
        'zzOneLoopDev': up(zs['one_loop']['closed_form_max_rel_dev']),
        'zzRouteBCAbeta': frac_tex(Fraction(ca_rows[0]['beta'])),
        'zzVanishMax': up(zs['one_loop']['vanishing_max_rel']),
        'zzLogDevMax': up(dev('log')), 'zzKDevMax': up(dev('K')), 'zzDevMax': up(max(dev('log'), dev('K'))),
        'zzShiftMax': up(shift), 'zzKcubicMax': up(kcubic), 'zzNoLogMax': up(nolog),
        'zzLowerModelDev': sci(2*zs['worst_log']['f2_lower_model'], 0),          # the text says "about"
        'zzTableBeta': frac_tex(Fraction(tb)),
        'zzTensorDev': up(chk['tensor_basis_max_rel_diff']),
        'zzBenchOneLoop': up(chk['benchmark_one_loop']['max_rel_diff']),
        'zzBenchCF': up(chk['benchmark_two_loop_CF']['max_rel_diff']),
        'zzBenchCA': up(chk['benchmark_two_loop_CA']['max_rel_diff']),
        'zzForwardCA': up(ca_cf['max_rel_diff']),
        'zzForwardCF': up(cf['max_rel_diff']), 'zzForwardCFBeta': frac_tex(Fraction(cf['beta'])),
        'zzForwardCFDelta': f'10^{{-{cf["delta_exp"]}}}',
        'zzCxxBenchOneLoop': up(bench['max_rel_diff']['1l']),
        'zzCxxBenchTwo': up(max(bench['max_rel_diff']['2l_CF'], bench['max_rel_diff']['2l_CA'])),
        'zzCxxFwdOneLoop': up(fwd['max_rel_diff']['1l']), 'zzCxxFwdCF': up(fwd['max_rel_diff']['2l_CF']),
        'zzCxxFwdCA': up(fwd['max_rel_diff']['2l_CA']),
        'zzCxxFwdBeta': frac_tex(Fraction(fwd['beta'])), 'zzCxxFwdDelta': f'10^{{-{fwd["delta_exp"]}}}',
        'zzLiValueDev': up(liv['max_rel_diff']),
        'zzpTsqMin': sci(min(float(c['pT2']) for v in zf.values() for c in v['curve']), 1),
        'zzpTsqMax': sci(max(float(c['pT2']) for v in zf.values() for c in v['curve']), 1),
    }
    words = {
        'zzVelocityCount': NUMBER_WORDS[len(betas)],
        'zzDeltaCount': NUMBER_WORDS[len(scan['delta_exponents'])],
        'zzCtrlExtraDigits': str(ctrl['extra_digits']),
        'zzLiCalls': str(li['calls']), 'zzLiSwapped': str(li['convergent_if_swapped']),
        'zzBenchOneLoopCount': str(chk['benchmark_one_loop']['coefficients']),
        'zzForwardCFDigits': str(cf['dps']), 'zzForwardCFCount': str(cf['amplitudes']),
        'zzBenchCFList': alist(bench_cf['benchmark']), 'zzBenchCFListX': alist(bench_cf['t_u_exchanged']),
        'zzLiValuePoints': li_where, 'zzPzeroStatements': str(p2['statements']),
        'zzLiScanArgs': str(zls['arguments']),
        'zzLiScanDigits': NUMBER_WORDS[math.ceil(math.log10(zls['units']))] if zls['units'] > 1 else 'none',
        'zzFwdDigits': NUMBER_WORDS[fwd_digits], 'zzFwdDecades': NUMBER_WORDS[fwd_points-1],
    }

    # Table 7
    names = {'CF': '$C_F$', 'CA': '$C_A$'}
    hel = lambda hs: ', '.join('$\\hel'+h+'$' for h in hs)
    heads = {'log': 'amplitudes '+hel(ZZ_LOG)+', in units of $C_1$', 'K': 'amplitudes '+hel(ZZ_K)+', in units of $K$'}
    t = ['\\begin{table}[t]', '\\centering', '\\small',
         '\\begin{tabular}{lllllll}', '\\toprule',
         f'&&&\\multicolumn{{2}}{{c}}{{$\\beta={frac_tex(Fraction(tb))}$}}&\\multicolumn{{2}}{{c}}{{largest}}\\\\',
         '\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}',
         'coefficient & part & prediction & predicted & extracted & deviation & shift\\\\']
    for kind, powers in (('log', (3, 2)), ('K', (3, 2, 1))):
        t += ['\\midrule', f'\\multicolumn{{7}}{{l}}{{{heads[kind]}}}\\\\']
        for j in powers:
            for i, part in enumerate(('CF', 'CA')):
                row, pr = rows[(kind, j, part)], pred[(kind, j, part)]
                coeff = ('$\\G^3$', '$\\G^2$', '$\\G$')[3-j] if i == 0 else ''
                if pr == 0:
                    t.append(f'{coeff} & {names[part]} & $0$ & $0$ & ${sci(abs(row["value"]), 1)}$ & '
                             f'${sci(row["dev"], 1)}$ & ${sci(row["shift"], 1)}$\\\\')
                else:
                    t.append(f'{coeff} & {names[part]} & ${lin_tex(pr)}$ & ${float(row["pred"]):.10f}$ & '
                             f'${float(mp.re(row["value"])):.10f}$ & ${sci(row["dev"], 1)}$ & ${sci(row["shift"], 1)}$\\\\')
    t += ['\\bottomrule', '\\end{tabular}',
          '\\caption{Coefficients of the leading logarithms of the two-loop $gg\\to ZZ$ amplitudes in the forward '
          'limit, by colour part, in units of the one-loop coefficients, both expanded in $a_s$, '
          'compared with eq.~(\\ref{eq:zzpred}); $L_x=\\ln[x(1-x)]$. The extracted values are those of the first '
          f'amplitude of each group at $\\beta={frac_tex(Fraction(tb))}$: the real part where the prediction is '
          'nonzero, the modulus where it vanishes. The last two columns give the largest '
          'modulus of the difference from the prediction and the largest change when the fit window is moved one '
          f'step up in $\\delta$, over the four amplitudes of each group and the {NUMBER_WORDS[len(betas)]} velocities.}}',
          '\\label{tab:zz}', '\\end{table}']
    (HERE/'tables'/'table_zz.tex').write_text('\n'.join(t)+'\n')
    return macros, words



# ---------------------------------------------------------------- Fig. 2: from 90 degrees to the forward limit
def forward_part(aa, zw, zz, tol=mp.mpf(1)/100, below=6):
    """Checks what Secs. 2.3, 2.5 and 6.2 and the caption of Fig. 2 say about the curves, and
    returns their macros:
      dpsOnset: the largest p_T^2/m^2 = 10^-k below which, at all velocities and for both
        processes, the derivative of the amplitude with Delta J_z = 0 with respect to G, in units
        of C_1, and the amplitude with |Delta J_z| = 2, in units of K, differ from 1 by less than
        one per cent (the derivative taken as in the figure: the difference quotient between
        neighbouring points, at the geometric mean of their p_T^2);
      aaRealLogDev: along each line, Re M_1/C_1 of gg -> AA equals (1/2) ln(p_T^2/m_T^2) plus a
        constant, m_T^2 = m^2 + p_T^2, to this accuracy (the largest spread of the difference),
        which is why the real parts in Fig. 2(a) coincide;
      aaApproach...: (M_2 - K)/(K p_T^2) of gg -> AA at p_T^2 <= 10^-below m^2 is constant to so
        many digits over so many decades, as for gg -> HH in Sec. 2.5."""
    mp.mp.dps = 30
    R = lambda v: mp.mpc(mp.mpf(v[0]), mp.mpf(v[1])) if isinstance(v, list) else mp.mpf(v)
    betas = sorted(zz['betas'], key=Fraction)
    assert sorted(aa, key=Fraction) == sorted(zw, key=Fraction) == betas
    outside = []            # p_T^2 of the values that differ from 1 by one per cent or more
    real_log, approach, decades = 0, 0, []
    for proc, data in (('AA', aa), ('ZZ', zw)):
        for b in betas:
            v = data[b]
            beta = mpq(Fraction(b))
            x, s = (1+beta)/2, 4/(1-beta**2)              # s in units of m^2
            rows = v['rows']
            pT2 = [mp.mpf(r['pT2']) for r in rows]
            assert all(p > q for p, q in zip(pT2, pT2[1:])), 'ordered from large to small p_T^2'
            assert abs(pT2[0]*(1-beta**2)/beta**2-1) < 1e-20, 'each curve begins at 90 degrees'
            G = [R(r['G']) for r in rows]
            assert all(abs(g-mp.log(s*p)/2-mp.pi*1j/2) < 1e-20 for g, p in zip(G, pT2)), 'G of eq. (2.8)'
            m1 = [R(r['M1_over_C1_minus_G'])+g for r, g in zip(rows, G)]          # M/C_1
            d1 = [((m1[k+1]-m1[k])/(G[k+1]-G[k]), mp.sqrt(pT2[k]*pT2[k+1])) for k in range(len(rows)-1)]
            r2 = [(R(r['M2_over_K']), p) for r, p in zip(rows, pT2)]
            outside += [p for val, p in d1+r2 if abs(val-1) >= tol]
            assert abs(d1[-1][0]-1) < 1e-6 and abs(r2[-1][0]-1) < 1e-6, 'both tend to 1'
            if proc == 'AA':        # C_1 = 2 pi N_c m^2 F_A(x) and K = pi N_c m^2 W(x)
                FA, W = 1/(x*(1-x))-2*x*(1-x), (x**4+(1-x)**4)/(x*(1-x))
                assert abs(R(v['C1'])/(6*mp.pi*FA)-1) < 1e-25 and abs(R(v['K'])/(3*mp.pi*W)-1) < 1e-25
                c = [mp.re(m)-mp.log(p/(1+p))/2 for m, p in zip(m1, pT2)]
                real_log = max(real_log, max(c)-min(c))
                eta = [(val-1)/p for val, p in r2 if p <= mp.mpf(10)**-below]
                mean = sum(eta)/len(eta)
                approach = max(approach, max(abs(e1-e2) for e1 in eta for e2 in eta)/abs(mean))
                sel = [p for p in pT2 if p <= mp.mpf(10)**-below]
                decades.append(mp.log10(sel[0]/sel[-1]))
            else:                   # the closed forms of Sec. 6.2 in the normalisation of the source,
                C1, K = -8j*mp.pi*x**2, -4j*mp.pi*x*(1-x)          # and the fitted coefficients
                hz = zz['betas'][b]['hel']
                assert abs(R(v['C1'])-C1) < 1e-25 and abs(R(v['K'])-K) < 1e-25
                assert abs(R(hz['++++']['c1'])/C1-1) < 1e-15 and abs(R(hz['+-++']['K'])/K-1) < 1e-15
    assert real_log < 1e-25, 'Fig. 2(a): Re M_1/C_1 - (1/2) ln(p_T^2/m_T^2) is constant'
    k = int(mp.ceil(-mp.log10(min(outside))))
    return {'dpsOnset': f'10^{{-{k}}}', 'aaRealLogDev': up(real_log), 'aaApproachBelow': f'10^{{-{below}}}'},         {'aaApproachDigits': NUMBER_WORDS[digits(approach)],
         'aaApproachDecades': NUMBER_WORDS[int(mp.floor(min(decades)+mp.mpf(10)**-20))]}


# ---------------------------------------------------------------- appendix B: Fig. 5, the fit windows
def windows_part(hw, zw, colour, zz):
    """Checks the statements of appendix B about Fig. 5 on all coefficients and velocities, not
    only those shown, and returns the velocity of the figure. The rows of each model are the
    windows in the order of the scan, so that the last row is the window of smallest delta."""
    last = lambda rows, m: float([r for r in rows if r['model'] == m][-1]['rel_dev'])
    # the last window of the model of the paper is the fit of Tables 2 and 7 (caption of Fig. 5)
    for b in hw:
        for k, v in hw[b].items():
            row, = [r for r in colour[b]['rows'] if r['coefficient'] == k and r['colour'] == 'N2']
            assert math.isclose(float(row['deviation']), last(v['rows'], 'pT2 G^2'), rel_tol=1e-6), (b, k)
    mp.mp.dps = 30
    R = lambda v: mp.mpc(mp.mpf(v[0]), mp.mpf(v[1]))
    pred = zz_predictions()
    for b in zw:
        x = (1+mpq(Fraction(b)))/2
        for k, v in zw[b].items():
            h, part = k.split()                       # e.g. '++++', '2l_CA'
            kind, key = ('log', 'c1') if h in ZZ_LOG else ('K', 'K')
            hz = zz['betas'][b]['hel'][h]
            for power in {r['power'] for r in v['rows']}:
                a0, a1 = lin_coeffs(pred[(kind, power, part[3:])])
                pr = mpq(a0)+mpq(a1)*mp.log(x*(1-x))
                dev = abs(2*R(hz[part]['f'][power])/R(hz[key])-pr)/abs(pr)
                assert math.isclose(float(dev), last([r for r in v['rows'] if r['power'] == power], 'pT2 G^3'),
                                    rel_tol=1e-6), (b, k, power)
    sets = [(b, k, v['rows'], 'pT2 G^2') for b in hw for k, v in hw[b].items()]
    sets += [(b, f'{k} G^{p}', [r for r in v['rows'] if r['power'] == p], 'pT2 G^3')
             for b in zw for k, v in zw[b].items() for p in {r['power'] for r in v['rows']}]
    for b, k, rows, paper in sets:
        models = list(dict.fromkeys(r['model'] for r in rows))
        for m in models:                    # the deviation decreases as the window moves to smaller delta
            dev = [float(r['rel_dev']) for r in rows if r['model'] == m]
            assert all(y < x for x, y in zip(dev, dev[1:])), (b, k, m)
        # the paper model is better than the fit without p_T^2 terms by orders of magnitude
        assert last(rows, paper) < 1e-3*last(rows, 'plain'), (b, k)
        if paper == 'pT2 G^2' or k.startswith('++++'):    # gg -> HH, and the first gg -> ZZ group:
            assert all(last(rows, paper) <= last(rows, m) for m in models), (b, k)   # the paper model is best
        else:                               # the second gg -> ZZ group: P = 2 does better than P = 3
            assert last(rows, 'pT2 G^2') < last(rows, 'pT2 G^3'), (b, k)
    assert set(hw) == set(zw)
    # the statements about the two gg -> ZZ groups hold for all their amplitudes: results_zz.json
    # has the fit with P = 2 through the window of smallest delta for each of them
    for b in zz['betas']:
        x = (1+mpq(Fraction(b)))/2
        for kind, hels, key, powers in (('log', ZZ_LOG, 'c1', (3, 2)), ('K', ZZ_K, 'K', (2, 1))):
            for h in hels:
                hz = zz['betas'][b]['hel'][h]
                for part in ('CF', 'CA'):
                    for power in powers:
                        a0, a1 = lin_coeffs(pred[(kind, power, part)])
                        pr = mpq(a0)+mpq(a1)*mp.log(x*(1-x))
                        dev = lambda f: abs(2*R(f[power])/R(hz[key])-pr)/abs(pr)
                        p3, p2 = dev(hz['2l_'+part]['f']), dev(hz['2l_'+part]['f_lower_model'])
                        if kind == 'log':           # first group: P = 3 better by orders of magnitude
                            assert p3 < 1e-3*p2, (b, h, part, power)
                        else:                       # second group: P = 2 better
                            assert p2 < p3, (b, h, part, power)
    return {'windowBeta': frac_tex(Fraction(preferred_beta(list(hw))))}


# ---------------------------------------------------------------- Sec. 5.1, Tables 5 and 6: planar N = 4
N4_NAMES = {'r0': '\\kappa_1', 'r0b': '\\kappa_2', 'r0c': '\\kappa_3', 'r0d': '\\kappa_4', 'B2': 'B^{(2)}', 'B3': 'B^{(3)}',
            'B4': 'B^{(4)}', 'G4': '\\gamma_K^{(4)}'}


def n4_tex(expr):
    """a coefficient of the Dixon-Esterlis data or a solved constant, in the notation z2, zeta_2, ..., as LaTeX"""
    loc = {f'z{k}': sp.Symbol(f'zeta_{k}') for k in range(2, 10)} | {f'zeta_{k}': sp.Symbol(f'zeta_{k}') for k in range(2, 10)}
    e = sp.sympify(expr, locals=loc)
    return sp.latex(sp.nsimplify(e)) if e.is_number else sp.latex(e)


def n4_lines(expr, per_line=6, first=''):
    """LaTeX of a sum, broken into lines of at most per_line terms (for the aligned environments)"""
    loc = {f'zeta_{k}': sp.Symbol(f'zeta_{k}') for k in range(2, 10)}
    e = sp.expand(sp.sympify(expr, locals=loc))
    weight = lambda t: sum(k*sp.degree(t, loc[f'zeta_{k}']) for k in range(2, 10) if t.has(loc[f'zeta_{k}']))
    terms = sorted(sp.Add.make_args(e), key=lambda t: (-weight(t), str(t)))       # by transcendental weight
    pieces = []
    for t in terms:
        tex = sp.latex(t)
        pieces.append(('-' + tex[1:].strip() if tex.startswith('-') else '+' + tex.strip()))
    lines = [''.join(pieces[i:i + per_line]) for i in range(0, len(pieces), per_line)]
    lines[0] = lines[0][1:] if lines[0].startswith('+') else lines[0]
    return [first + lines[0]] + ['&\\quad ' + l for l in lines[1:]]


def n4_part(res):
    """Checks the statements of Sec. 5.1 about the comparison with planar N = 4 and writes Tables 5 and 6."""
    rows, nm = res['three_to_three'], res['nmhv']
    assert res['n_pass'] + res['n_fixed'] + res['n_norm'] == res['n_checked'] == len(rows)
    assert not any(r['status'].startswith('FAIL') for r in rows)
    assert res['four_loop_cusp_agrees'] and res['LL_Gaussian_integral'] and res['nmhv_LL_exponential']
    assert res['components_related'] and res['two_to_four_single_log_difference'] == '0'
    assert all(r['DE'] == r['prediction'] for r in res['two_to_four_two_loops'])   # 2 -> 4: ln^3, ln^2, ln, no v
    two4 = {(r['L'], r['k']): r['DE'] for r in res['two_to_four_two_loops']}
    assert two4 == {('1', 1): '-I*pi', ('2', 3): 'I*pi/6', ('2', 2): '-pi**2/2', ('2', 1): '-2*I*pi**3/3'}, two4
    assert res['three_loop_ln3_coefficient'].startswith('C_1 (-2/3) (kappa_1 + gamma_K^(2)/4) G^3')
    norm = [r for r in rows if r['status'] == 'normalisation']
    assert len(norm) == 1 and (norm[0]['L'], norm[0]['k']) == (1, 1) and norm[0]['DE'] == '1/2'
    assert not any(r['status'].startswith('FAIL') for r in nm)
    norm = [r for r in nm if r['status'] == 'normalisation']
    assert len(norm) == 1 and (norm[0]['L'], norm[0]['k']) == (1, 0) and norm[0]['DE'] == '1'
    assert res['nmhv_n_pass'] == res['nmhv_n_checked'] == len(nm) - 1
    assert {(r['L'], r['k']) for r in nm} == {(L, k) for L in range(1, 5) for k in range(2*L - 1)}   # all coefficients, k = 0 included
    five = res['nmhv_five_loop_prediction']
    assert sorted(five, key=int) == [str(k) for k in range(9)] and five['8'] == '1/384'
    # the NLL tower vanishes: every k = 2L-2 row is a PASS with coefficient 0
    for r in rows:
        if r['k'] == 2*r['L'] - 2:
            assert r['status'] == 'PASS' and r['DE'] == '0', r
    first = {}                                   # the order at which each constant is fixed
    for r in rows:
        if r['status'].startswith('fixes'):
            first[r['status'].split()[1]] = (r['L'], r['k'])

    def status(r, norm):
        if r['status'] == 'PASS':
            return 'agrees'
        if r['status'].startswith('fixes'):
            return f'fixes ${N4_NAMES[r["status"].split()[1]]}$'
        assert r['status'] == 'normalisation', r
        return f'defines ${norm}$'

    def block(rs, norm):
        return [f'{r["L"]} & {r["k"]} & ${n4_tex(r["DE"])}$ & {status(r, norm)}' for r in rs]

    def table(rs, label, caption, name, norm):
        half = (len(rs) + 1)//2
        left, right = block(rs[:half], norm), block(rs[half:], norm) + [' & & & ']*(half - len(rs[half:]))
        t = ['\\begin{table}[t]', '\\centering', '\\small', '\\setlength{\\tabcolsep}{4pt}',
             '\\begin{tabular}{ccll@{\\hspace{14pt}}|@{\\hspace{14pt}}ccll}', '\\toprule',
             '$L$ & $k$ & coefficient & result & $L$ & $k$ & coefficient & result\\\\', '\\midrule']
        t += [f'{l} & {r}\\\\' for l, r in zip(left, right)]
        t += ['\\bottomrule', '\\end{tabular}', f'\\caption{{{caption}}}', f'\\label{{{label}}}', '\\end{table}']
        (HERE/'tables'/f'{name}.tex').write_text('\n'.join(t)+'\n')

    table(rows, 'tab:n4',
          'The coefficients of $2\\pi i\\ln^k|\\delta_*|$ in the $L$-loop MHV amplitude $\\mathcal E^{(L)}$ of ref.~\\cite{Dixon:2016epj}, '
          'eqs.~(3.17)--(3.21), in $3\\to3$ kinematics, and the comparison with the structure of section~\\ref{sec:rg} '
          'in the conformal theory: each fitted constant is fixed by the first coefficient in which it appears '
          '(the constants are given in eq.~(\\ref{eq:n4const})); every other coefficient is a prediction.',
          'table_n4', 'C_1')
    table(nm, 'tab:n4nmhv',
          'The same for the next-to-MHV amplitude: coefficients of $\\pi i\\,\\frac{1+w}{1+w^*}\\ln^k|\\delta_*|$ in $\\rho^{(L)}$, '
          'including the constants $k=0$, eqs.~(D.14)--(D.17) of the erratum to ref.~\\cite{Dixon:2016epj}, compared with the '
          '$|\\Delta J_z|=2$ component of the kernel, transformed with $\\nu=2$ in eq.~(\\ref{eq:Fm}), with the constants of '
          'table~\\ref{tab:n4} and no new one; the one-loop term defines the normalisation $K$.',
          'table_n4nmhv', 'K')
    g4 = res['solved_constants']['G4']
    # the constants have uniform transcendental weight: 2n for kappa_n, 2n - 1 for B^(n) (Sec. 5.1)
    loc = {f'zeta_{k}': sp.Symbol(f'zeta_{k}') for k in range(2, 10)}
    weight = lambda t: sum(k*sp.degree(t, loc[f'zeta_{k}']) for k in range(2, 10) if t.has(loc[f'zeta_{k}']))
    for name, w in (('r0', 2), ('r0b', 4), ('r0c', 6), ('r0d', 8), ('B2', 3), ('B3', 5), ('B4', 7), ('G4', 6)):
        e = sp.expand(sp.sympify(res['solved_constants'][name], locals=loc))
        assert all(weight(t) == w for t in sp.Add.make_args(e)), (name, e)
    rows = (('r0', None), ('r0b', 'B2'), ('r0c', 'B3'), ('r0d', 'B4'))       # one line of eq. (5.8) per loop order
    consts = []
    for i, (k, b) in enumerate(rows):
        line = f'{N4_NAMES[k]}&={n4_tex(res["solved_constants"][k])},'
        if b:
            line += f' &\\qquad {N4_NAMES[b]}&={n4_tex(res["solved_constants"][b])}' + (',' if i < len(rows) - 1 else '.')
        consts.append(line)
    consts = '\\\\\n'.join(consts)
    # the five-loop prediction for the |Delta J_z| = 2 component, one line per power, long coefficients broken
    lines = []
    for k in range(8, -1, -1):
        power = {0: '', 1: '\\ln|\\delta_*|'}.get(k, f'\\ln^{k}|\\delta_*|')
        coef = [c.replace('\\frac', '\\tfrac') for c in n4_lines(five[str(k)], per_line=8)]
        if len(coef) == 1 and '+' not in coef[0][1:] and '-' not in coef[0][1:]:     # a single term
            body = [coef[0] + ('\\,' + power if power else '')]
        else:
            body = ['\\big(' + coef[0]] + coef[1:]
            body[-1] += '\\big)' + ('\\,' + power if power else '')
        body[0] = ('&' if k == 8 else '&+') + body[0] if not body[0].startswith('-') else '&' + body[0]
        lines += body
    five_tex = '\\notag\\\\\n'.join(lines)       # for an align environment that may break across pages
    return {'nFourCuspFour': n4_tex(g4), 'nFourConstants': consts}, \
        {'nFourChecked': str(res['n_checked']), 'nFourPass': str(res['n_pass']), 'nFourFixed': NUMBER_WORDS[res['n_fixed']],
         'nFourLoops': NUMBER_WORDS[res['loops_compared']], 'nFourNmhvChecked': str(res['nmhv_n_checked']),
         'nFourNmhvPass': str(res['nmhv_n_pass']), 'nFourNmhvFive': five_tex}      # the last one is math, for an aligned environment


# ---------------------------------------------------------------- Sec. 5.1: the five-loop next-to-MHV coefficients
def nmhv_primary_part(prim, res):
    """The comparison of n4/check_primary_nmhv.py with the multi-Regge amplitude: the counts the text quotes,
    and the check that the five-loop coefficients compared there are those printed in eq. (n4five)."""
    comps = prim['comparisons']
    assert prim['passed'] is True and all(str(c['difference']) == '0' for c in comps)
    low = [c for c in comps if c['loop'] <= 4]
    five = {c['power']: c for c in comps if c['loop'] == 5}
    assert len(low) + len(five) == len(comps)
    assert len(low) == prim['lower_loop_coefficients'] == res['nmhv_n_checked'] + 1, \
        'the text counts the one-loop coefficient, which defines K, among those through four loops'
    assert len(five) == prim['five_loop_coefficients'] and set(five) == {int(k) for k in res['nmhv_five_loop_prediction']}
    pi2 = sp.pi**2
    even = {2: pi2/6, 4: pi2**2/90, 6: pi2**3/945, 8: pi2**4/9450}
    odd = {k: sp.Symbol(f'z{k}') for k in (3, 5, 7, 9)}
    printed_names = {f'zeta_{k}': (even[k] if k in even else odd[k]) for k in range(2, 10)}
    primary_names = {'p2': pi2, **{f'z{k}': v for k, v in odd.items()}}
    for k, c in five.items():
        printed = sp.sympify(res['nmhv_five_loop_prediction'][str(k)], locals=printed_names)
        assert sp.expand(printed - sp.sympify(c['target'], locals=primary_names)) == 0, \
            f'the coefficient of ln^{k} compared with the multi-Regge amplitude is not that of eq. (n4five)'
        assert sp.expand(sp.sympify(c['primary'], locals=primary_names) - printed) == 0, k
    return {'nFourNmhvFiveCount': NUMBER_WORDS[len(five)], 'nFourNmhvLower': str(len(low))}


# ---------------------------------------------------------------- appendix E: ultrasoft attachments with the complete numerators
def ultrasoft_part(us):
    """Checks what the text says about checks/ultrasoft_numerators.py (exact arithmetic, stored in
    checks/ultrasoft_numerators.json) and writes the counts it quotes."""
    assert us['passed'] is True and us['arithmetic'].startswith('exact')
    one_q, one_g = us['C g->qqbar, one gluon'], us['D g->gg, one gluon']
    two_q, two_g = us['E g->qqbar, two gluons'], us['F g->gg, two gluons']
    assert us['F generic g->gg builder reproduces section D'] is True
    # one gluon: the three attachments and the Wilson line are leading, the contact terms down by lambda^2
    rel = one_g['relative power (at least)']
    contact = [v for k, v in rel.items() if k[0] in 'de']
    assert all(rel[k] == 0 for k in ('a', 'b', 'c', 'W')) and len(contact) == 5 and min(contact) == 2, \
        'the text says that the contact terms are suppressed by lambda^2'
    # the other polarisations of the gluon: G and the transverse ones down by lambda, p1 by lambda^2
    pols = [one_q['relative power of the other polarisations (at least)'],
            {k.replace('polarisation ', ''): v for k, v in rel.items() if k.startswith('polarisation ')}]
    for pol in pols:
        assert pol == {'G': 1, 'transverse x': 1, 'transverse y': 1, 'p1': 2}, pol
    # two gluons: twelve leading configurations for either splitting; every other graph down by lambda^2 or zero
    orders = two_q['numerator orders (None: vanishes identically)']
    lead_q = [k for k, v in orders.items() if v is not None and not k.startswith('W')]
    assert [k for k, v in orders.items() if v is None] == ['s1', 's2', 's3'] and len(set(orders[k] for k in lead_q)) == 1
    lead_g, other = two_g['leading-power graphs'], two_g['relative power of the other graphs (at least)']
    zero = two_g['vanish identically (all colour structures)']
    assert len(lead_q) == len(lead_g) == 12, 'the text says twelve'
    assert two_g['graphs'] == len(lead_g) + len(other) + len(zero) and min(other.values()) == 2 and len(zero) == 1, \
        'the text says that the other graphs are suppressed by lambda^2 or vanish'
    return {'usTwoGluonGraphs': str(two_g['graphs']), 'usTwoGluonLeading': NUMBER_WORDS[len(lead_g)]}


def ultrasoft_dimension_part(folder, four):
    """The runs of checks/ultrasoft_numerators.py with the algebra in d dimensions (records in
    checks/ultrasoft_dimensions_20261112) against what the text says, and the dimensions it quotes."""
    rec = {p.stem: json.loads(p.read_text()) for p in sorted(folder.glob('*.json'))}
    lo, hi, dep, lin = rec['n6_d6'], rec['n6_d7'], rec['dependence_n6'], rec['dependence_n5']
    assert all(r['passed'] is True for r in rec.values())
    assert lo['dimension'] == lo['vector dimension'] == hi['vector dimension'] == dep['vector dimension'], \
        'the text says that the momenta and polarisations have the lower of the two dimensions'
    assert dep['dimensions'] == [lo['dimension'], hi['dimension']] and hi['dimension'] == lo['dimension'] + 1
    assert lo['seed'] == hi['seed'] == four['seed'] and lo['generators'] == hi['generators'], 'the same points'
    assert set(lo['sections']) == set('ABCDEF') and {'D', 'F'} <= set(hi['sections']), \
        'every section in the lower dimension, the g -> gg sections in the higher'
    # the structure of the four-dimensional run: contact terms, other polarisations, the twelve leading graphs
    one, two = 'D g->gg, one gluon', 'F g->gg, two gluons'
    for r in (lo, hi):
        rel = r[one]['relative power (at least)']
        ref = four[one]['relative power (at least)']
        assert all(rel[k] == v for k, v in ref.items()), 'the powers of the four-dimensional run'
        extra = [v for k, v in rel.items() if k not in ref]
        assert extra and min(extra) == 1, 'the polarisations along the additional directions are down by lambda'
        for key in ('graphs', 'leading-power graphs', 'relative power of the other graphs (at least)',
                    'vanish identically (all colour structures)'):
            assert r[two][key] == four[two][key], key
    quark = lo['C g->qqbar, one gluon']['relative power of the other polarisations (at least)']
    assert quark['G'] == 1 and quark['p1'] == 2 and all(v == 1 for k, v in quark.items() if k.startswith('transverse'))
    assert lo['E g->qqbar, two gluons']['numerator orders (None: vanishes identically)'] \
        == four['E g->qqbar, two gluons']['numerator orders (None: vanishes identically)']
    # the g -> gg numerators depend on the dimension, linearly
    assert all(v['of which depend on the dimension'] > 0 for v in dep['sections'].values())
    assert dep['sections'][one]['of which depend on the dimension'] == dep['sections'][one]['Wilson-line terms compared']
    assert len(lin['dimensions']) >= 3 and all(v.get('linear in the dimension') is True for v in lin['sections'].values()), \
        'the text says that a gluon numerator is linear in d'
    return {'usDimLow': NUMBER_WORDS[lo['dimension']], 'usDimHigh': NUMBER_WORDS[hi['dimension']]}


# ---------------------------------------------------------------- appendix E, Sec. 5.1: the Collins-Soper equation of the kernel
def cs_poly_tex(rows, key, loops):
    """sum_n c_n a^n from the records of checks/collins_soper_kernel.py: {monomial: coefficient} per loop order"""
    zeta = {f'zeta_{k}': sp.Symbol(f'zeta_{k}') for k in range(2, 10)}
    pieces = []
    for n in loops:
        terms = rows[f'B^({n})'][key]
        e = sum(sp.Rational(c)*sp.Mul(*[zeta[z] for z in mono.split()]) for mono, c in terms.items())
        negative = all(Fraction(c) < 0 for c in terms.values())
        tex = sp.latex(-e if negative else e).replace('\\frac', '\\tfrac')
        if len(terms) > 1:
            tex = '\\big(' + tex + '\\big)'
        pieces.append(('-' if negative else '+') + tex + f'a^{n}')
    out = ''.join(pieces)
    return out[1:] if out.startswith('+') else out


def collins_soper_part(cs, res):
    """checks/collins_soper_kernel.py (stored in checks/collins_soper_kernel.json) against what the text says about
    the colour structure of the collinear functions, the soft function, the rapidity anomalous dimension and the
    two parts of B in planar N = 4; writes the expressions and counts the text quotes."""
    assert cs['passed'] is True and cs['arithmetic'].startswith('exact')
    inv = cs['A1 invariant tensors']
    assert inv['(8; 3, 3bar)'] == 1 and inv['(8; 8, 8)'] == 2, 'the text says: one tensor for g -> q qbar, f and d for g -> gg'
    assert cs['A3 charge conjugation'] == {'f^ABC': 'even', 'd^ABC': 'odd'}, 'the text says that charge conjugation excludes d'
    p = cs['A4 the collinear tensor']['p = (C_a - C_k - C_kbar)/2']
    assert p['g -> q qbar (t^A)'] == '1/6' and p['g -> g g (f^ABC)'] == '-3/2', 'T_k.T_kbar = -Delta_C/2 on the tensor'
    for kind in ('g -> q qbar', 'g -> g g'):
        four = cs[f'A5 four lines, {kind}']
        assert four['Psi is an eigenvector of D'] is False and not four['<D^2> - <D>^2'].startswith('0 '), \
            'the text says that the exchange between the two positions is not diagonal and <D^2> differs from <D>^2'
    assert inv['singlets of q qbar qbar q'] == 2 and 'matrix of D on the two singlets (columns: images)' in cs['A5 four lines, g -> q qbar']
    red = cs['A6 reduction with the vertex conditions']
    n_words = red['colour factors of two exchanges']
    assert len(red['two exchanges']) == n_words and len(red['three-gluon vertex on three lines']) == 4
    # planar N = 4: the two parts of B add up to the constants printed in eq. (n4const)
    B = cs['B anomalous dimensions']
    rows, loops = B['planar N = 4, coupling a = 2 g^2'], B['loops through which B = D_gamma + gamma_rap is checked']
    zeta = {f'zeta_{k}': sp.Symbol(f'zeta_{k}') for k in range(2, 10)}
    for n in range(2, loops + 1):
        row = rows[f'B^({n})']
        assert row['difference'] == '0' and row['sum'] == row['paper2/n4/results_stage_d.json']
        total = sum(sp.Rational(c)*sp.Mul(*[zeta[z] for z in mono.split()]) for mono, c in row['sum'].items())
        printed = sp.sympify(res['solved_constants'][f'B{n}'], locals=zeta)
        assert sp.expand(total - printed) == 0, f'B^({n}) of eq. (n4const)'
        parts = [sum(sp.Rational(c)*sp.Mul(*[zeta[z] for z in mono.split()]) for mono, c in row[k].items())
                 for k in ('D_gamma', 'gamma_rap')]
        assert sp.expand(parts[0] + parts[1] - printed) == 0
    assert loops == 4 and set(rows) == {'B^(2)', 'B^(3)', 'B^(4)'}, 'the text says B^(2), B^(3) and B^(4)'
    # QCD: the two-loop coefficient the text prints
    d2 = B['QCD']['d_2 = -gamma_rap^(2)/Delta_C']
    assert set(d2) == {'C_A', 'C_A zeta_3', 'T_F n_f'} and Fraction(d2['C_A zeta_3']).denominator == 1
    assert Fraction(d2['C_A']) > 0 > Fraction(d2['C_A zeta_3']) and Fraction(d2['T_F n_f']) < 0
    d2_tex = (f'\\big({frac_tex(d2["C_A"])}-{abs(Fraction(d2["C_A zeta_3"]))}\\zeta_3\\big)C_A'
              f'-{frac_tex(abs(Fraction(d2["T_F n_f"])))}T_Fn_f')
    model = cs['D three-loop model']['fourth logarithm']
    assert model['nu = 0']['fourth power of G'] == 2 and model['nu = 2']['fourth power of G'] == 1, \
        'the text says: the coefficients of G^2 in M_1 and of G in M_2'
    assert model['nu = 0']['coefficient'].startswith('-1*(DG1 + beta1 + r1_3)')
    assert model['nu = 2']['coefficient'].startswith('-2*(DG1 + beta1 + r1_3)')
    return {'csDgammaNFour': cs_poly_tex(rows, 'D_gamma', range(2, loops + 1)),
            'csGammaRapNFour': cs_poly_tex(rows, 'gamma_rap', range(2, loops + 1)),
            'csDtwo': d2_tex}, \
        {'csLoops': NUMBER_WORDS[loops], 'csColourFactors': str(n_words)}


# ---------------------------------------------------------------- appendix E: the boundary term at one loop
def kernel_part(res):
    """Checks the statements of appendix E (task K, Stage B) and writes the numbers it quotes."""
    t, mo, lg = res['telescoping'], res['model'], res['logs']
    # bookkeeping of the region scales: only ultrasoft regions give r2, and r1 is the eps^0 term of the
    # rapidity coefficient plus ultrasoft terms; the rapidity coefficient has no eps^0 term
    assert lg['r2'] == 'c_u/2' and lg['r1'] == 'c_u*kappa_u + d_u + g1'
    assert lg['rapidity_coefficient_eps0'] == '0' and lg['triangle_cu'] == -1 and lg['triangle_X2'] == '-1/2'
    # the ultrasoft cancellation: both colour components, corrections falling with lambda, one attachment O(1)
    us = t['ultrasoft']
    dev = [max(r['dev_Tk'], r['dev_Tkbar']) for r in us]
    assert all(dev[i+1] < dev[i]/5 for i in range(len(dev) - 1)) and all(r['single_attachment'] > 0.1 for r in us)
    assert min(min(g['dev_Tk'], g['dev_Tkbar']) for g in t['glauber']) > 0.05       # fails for Glauber momenta
    assert all(s['parent_over_daughter'] < 2*s['lambda'] for s in t['soft'])         # parent suppressed for soft
    assert t['residues_vs_quadrature_worst'] < 1e-5
    # the explicit two-loop graph
    lams = sorted({l for s in mo['log_slopes'] for l in s['lambda']})                  # the lambdas of the logarithm test
    for r in mo['ultrasoft_lambda6F'].values():
        ch = r['rel_changes']
        assert all(ch[i+1] < 0.6*ch[i] for i in range(len(ch) - 1))                  # lambda^6 F converges
    slope = max(s['rel_diff'] for s in mo['log_slopes'])
    real = max(s['real_part_change'] for s in mo['log_slopes'])
    assert slope < 1e-2 and real < 1e-3
    def pow10(v):
        e = round(math.log10(v)); assert abs(v - 10**e) < 1e-12*10**e; return f'10^{{{e}}}'
    return {'kernX': f'{mo["x"]:g}', 'kernLamMax': pow10(max(lams)), 'kernLamMin': pow10(min(lams)),
            'kernSlopeDev': up(slope), 'kernRealDev': up(real), 'kernUsDev': up(dev[-1]),
            'kernUsLam': pow10(us[-1]['lambda'])}, \
        {'kernSlopePoints': NUMBER_WORDS[len(mo['log_slopes'])]}


def r0_closed(x, which, nu=0):
    """the closed forms of eqs. (7.2)-(7.3) at momentum fraction x: 'CA', 'nf' (coefficient of 2 T_F n_f) for gluon
    splittings, 'CF', 'CAq' for quark splittings (as r0/data/r0_closed_forms.py)"""
    x = mp.mpf(x)
    z2 = mp.zeta(2)
    lx, ly = mp.log(x), mp.log(1 - x)
    F = 1/(x*(1 - x)) + 2*x*(1 - x) if nu == 0 else (x**4 + (1 - x)**4)/(x*(1 - x))
    return {'CA': mp.mpf(134)/9 - 5*z2 - 4*lx*ly + mp.mpf(2)/3/F, 'nf': -mp.mpf(20)/9 - mp.mpf(2)/3/F,
            'CF': 12 - 6*z2 + 2*(lx - ly)**2, 'CAq': 2 + z2 - 2*(lx**2 + ly**2)}[which]


def r0_fits_part(shh, szz, bhh, bzz):
    """Sec. 7.1: the fitted r0 of gg -> HH and of the eight gg -> ZZ amplitudes at the three scan velocities against
    the closed forms, with the window-shift bounds; the identification of the closed forms re-done from the scan
    values (exact solves for the rational coefficients, PSLQ for the constants); Table r0_zz. Returns the macros."""
    mp.mp.dps = 30
    # ---- gg -> HH: value, imaginary part, window shift per velocity, component and colour part
    dev_hh = bound_hh = ratio_hh = 0.0
    for b, r in shh.items():
        for comp, nu in (('nu0', 0), ('nu2', 2)):
            for key, which in (('r0_N', 'CA'), ('r0_nf', 'nf')):
                v = [mp.mpf(s) for s in r[comp][key]]
                dev = abs(mp.mpc(v[0] - r0_closed(r['x'], which, nu), v[1]))
                assert dev <= v[2], (b, comp, key, dev, v[2])                      # every deviation below its bound
                dev_hh, bound_hh, ratio_hh = max(dev_hh, dev), max(bound_hh, v[2]), max(ratio_hh, dev/v[2])
    # ---- gg -> ZZ: the same for the eight helicity amplitudes; the table shows the first amplitude of each group
    dev_zz = bound_zz = ratio_zz = spread_zz = 0.0
    groups = {'nu0': None, 'nu2': None}
    cols = [('nu0', '2l_CF'), ('nu0', '2l_CA'), ('nu2', '2l_CF'), ('nu2', '2l_CA')]
    col_dev = {c: 0.0 for c in cols}
    col_bound = {c: 0.0 for c in cols}
    rows = []
    for b in sorted(szz, key=Fraction):
        r = szz[b]
        x = Fraction(r['x']).limit_denominator(100)
        fitted, closed = [], []
        for comp, key in cols:
            hels = list(r[comp])
            assert len(hels) == 4, hels
            groups[comp] = groups[comp] or hels
            assert groups[comp] == hels
            which = 'CF' if key == '2l_CF' else 'CAq'
            cf = r0_closed(r['x'], which)
            vals = [[mp.mpf(s) for s in r[comp][h][key]] for h in hels]
            for v in vals:
                dev = abs(mp.mpc(v[0] - cf, v[1]))
                assert dev <= v[2], (b, comp, key, dev, v[2])
                dev_zz, bound_zz, ratio_zz = max(dev_zz, dev), max(bound_zz, v[2]), max(ratio_zz, dev/v[2])
                col_dev[comp, key] = max(col_dev[comp, key], dev)
                col_bound[comp, key] = max(col_bound[comp, key], v[2])
            spread_zz = max(spread_zz, max(abs(v[0] - vals[0][0]) for v in vals))
            fitted.append(f'{float(vals[0][0]):.8f}')
            closed.append(f'{float(cf):.8f}')
        rows.append((f'$x={frac_tex(x)}$', fitted, closed))
    assert dev_hh < 1e-8 and dev_zz < 1e-8 and spread_zz < 1e-9 and ratio_hh < 0.1 and ratio_zz < 0.1
    # ---- the identification of the closed forms, re-done from the scan values
    xs = [mp.mpf(r['x']) for r in shh.values()]
    assert len(xs) == 3
    z2 = mp.zeta(2)
    ident_res = ident_c = 0.0

    def pslq_check(c, expect):
        """PSLQ on the basis {1, zeta_2} for a constant known to ~1e-9: the relation must be the expected (p, q) with
        c = p + q zeta_2 (coefficients below 200; tolerance 1e-8, above the fit precision), and |c - p - q zeta_2| is
        returned"""
        rel = mp.pslq([c, 1, z2], tol=mp.mpf(10)**-8, maxcoeff=200, maxsteps=10**5)
        assert rel is not None and (Fraction(-rel[1], rel[0]), Fraction(-rel[2], rel[0])) == expect, (rel, expect)
        return abs(c - expect[0] - expect[1]*z2)
    for comp, nu in (('nu0', 0), ('nu2', 2)):
        F = [1/(x*(1 - x)) + 2*x*(1 - x) if nu == 0 else (x**4 + (1 - x)**4)/(x*(1 - x)) for x in xs]
        # C_A part: c + a ln x ln(1-x) + b/F_nu, three unknowns from the three velocities (exact solve)
        A = mp.matrix([[1, mp.log(x)*mp.log(1 - x), 1/f] for x, f in zip(xs, F)])
        y = mp.matrix([mp.mpf(r[comp]['r0_N'][0]) for r in shh.values()])
        c, a, bb = mp.lu_solve(A, y)
        ident_res = max(ident_res, abs(a + 4), abs(bb - mp.mpf(2)/3))
        ident_c = max(ident_c, pslq_check(c, (Fraction(134, 9), Fraction(-5))))
        # 2 T_F n_f part: c' + b'/F_nu, two unknowns; the third velocity is a check
        A = mp.matrix([[1, 1/f] for f in F[:2]])
        y = mp.matrix([mp.mpf(r[comp]['r0_nf'][0]) for r in list(shh.values())[:2]])
        c2, b2 = mp.lu_solve(A, y)
        third = mp.mpf(list(shh.values())[2][comp]['r0_nf'][0])
        ident_res = max(ident_res, abs(b2 + mp.mpf(2)/3), abs(c2 + b2/F[2] - third))
        ident_c = max(ident_c, pslq_check(c2, (Fraction(-20, 9), Fraction(0))))
    # quark splittings: C_F part c + a ln^2[x/(1-x)], C_A part c + a [ln^2 x + ln^2(1-x)]; the first helicity amplitude
    xz = [mp.mpf(r['x']) for r in szz.values()]
    for key, fn, a_exp, c_exp in (('2l_CF', lambda x: (mp.log(x) - mp.log(1 - x))**2, 2, (12, -6)),
                                  ('2l_CA', lambda x: mp.log(x)**2 + mp.log(1 - x)**2, -2, (2, 1))):
        vals = [mp.mpf(r['nu0'][groups['nu0'][0]][key][0]) for r in szz.values()]
        A = mp.matrix([[1, fn(x)] for x in xz[:2]])
        c, a = mp.lu_solve(A, mp.matrix(vals[:2]))
        ident_res = max(ident_res, abs(a - a_exp), abs(c + a*fn(xz[2]) - vals[2]))
        ident_c = max(ident_c, pslq_check(c, (Fraction(c_exp[0]), Fraction(c_exp[1]))))
    assert ident_res < 1e-7 and ident_c < 1e-8, (ident_res, ident_c)
    # ---- Table r0_zz
    t = ['\\begin{table}[t]', '\\centering', '\\footnotesize', '\\setlength{\\tabcolsep}{5pt}',
         '\\begin{tabular}{llcccc}', '\\toprule',
         '& & \\multicolumn{2}{c}{$\\nu=0$: $\\M_{\\hel++++}$} & \\multicolumn{2}{c}{$\\nu=2$: $\\M_{\\hel+-++}$}\\\\',
         '\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}',
         '& & $C_F$ & $C_A$ & $C_F$ & $C_A$\\\\', '\\midrule']
    for label, fitted, closed in rows:
        t.append(label + ' & amplitudes & ' + ' & '.join(f'${v}$' for v in fitted) + '\\\\')
        t.append(' & closed form & ' + ' & '.join(f'${v}$' for v in closed) + '\\\\')
    t += ['\\midrule',
          '\\multicolumn{2}{l}{largest difference} & ' + ' & '.join(f'${sci(col_dev[c], 1)}$' for c in cols) + '\\\\',
          '\\multicolumn{2}{l}{window shift} & ' + ' & '.join(f'${sci(col_bound[c], 1)}$' for c in cols) + '\\\\',
          '\\bottomrule', '\\end{tabular}',
          '\\caption{The constant $r_0$ of the one-loop kernel read off the two-loop $gg\\to ZZ$ amplitudes with '
          'eq.~(\\ref{eq:r0read}), by colour factor, at the three velocities of the scans, against eq.~(\\ref{eq:r0q}). '
          'The values shown are those of the first amplitude of each group; the other three amplitudes of each group '
          f'agree with them to ${sci(spread_zz, 1)}$. The row marked largest difference gives, for each column, the '
          'largest difference between the amplitudes and the closed form over the four amplitudes of the group and '
          'the three velocities, and the row marked window shift the largest change of the fitted value when the fit window '
          'is moved one step up in $\\delta$ (appendix~\\ref{app:fits}).}',
          '\\label{tab:r0zz}', '\\end{table}']
    (HERE/'tables'/'table_r0_zz.tex').write_text('\n'.join(t)+'\n')
    assert groups['nu0'][0] == '++++' and groups['nu2'][0] == '+-++', groups
    macros = {'rzBoundHH': sci(bound_hh, 1), 'rzBoundZZ': sci(bound_zz, 1),
              'rzRatioHH': f'{float(ratio_hh):.2f}', 'rzRatioZZ': f'{float(ratio_zz):.2f}',
              'rzIdentRes': sci(ident_res, 1), 'rzIdentC': sci(ident_c, 1), 'rzZZspread': sci(spread_zz, 1)}
    return macros, {'rzScanX': '$x=' + '$, $'.join(frac_tex(Fraction(float(x)).limit_denominator(100)) for x in xs[:-1])
                    + f'$ and ${frac_tex(Fraction(float(xs[-1])).limit_denominator(100))}$'}


def r0_part(fits, bhh, bzz, rq, rg):
    """Sec. 7 and App. F: the constant r0, from the fits (closed forms, blind tests) and from the regions
    (Table 9); checks what the text says and returns the macros."""
    blind_hh = max(max(r['dev_C_A'], r['dev_n_f']) for r in bhh)
    blind_hh_unc = max(max(r['shift_C_A'], r['shift_n_f']) for r in bhh)
    blind_zz = max(max(r['2l_CF']['dev'], r['2l_CA']['dev']) for r in bzz)
    blind_zz_unc = max(max(r['2l_CF']['shift'], r['2l_CA']['shift']) for r in bzz)
    assert all(abs(r['C_A'][1]) < 1e-8 and abs(r['n_f'][1]) < 1e-8 for r in bhh)          # imaginary parts cancel
    assert all(abs(r[k]['r0'][1]) < 1e-8 for r in bzz for k in ('2l_CF', '2l_CA'))
    assert blind_hh < 1e-8 and blind_zz < 1e-8 and fits['HH'] < 1e-8 and fits['ZZ'] < 1e-8
    xs = sorted(set(rq['x']) & set(rg['x']), key=Fraction)
    assert len(xs) == 3, xs
    eta, pole, imag, dq, dg, spread = 0, 0, 0, 0, 0, 0
    rows = []
    for xk in xs:
        q, g = rq['x'][xk], rg['x'][xk]
        row = [f'{float(Fraction(xk)):.1f}']
        for c in ('CF', 'CA'):
            t = q[c]['theory, J1 structure ek, J2 structure ek']
            eta = max(eta, t['max |1/eta coefficient|']); imag = max(imag, t['max |imag|'])
            ex = q[c]['expected poles']
            pole = max(pole, abs(t['eps^-2'] - ex['eps^-2']), abs(t['eps^-1'] - ex['eps^-1']))
            spread = max(spread, q[c]['max spread of eps^0 over the four helicity configurations'])
            dq = max(dq, abs(q[c]['difference']))
            row.append(t['eps^0'])
        for nu in ('nu=0', 'nu=2'):
            for c in ('CA', 'nf'):
                t = g[c][nu]
                eta = max(eta, t['max |1/eta coefficient|']); imag = max(imag, t['max |imag|'])
                ex = t['expected poles']
                pole = max(pole, abs(t['eps^-2'] - ex['eps^-2']), abs(t['eps^-1'] - ex['eps^-1']))
                dg = max(dg, abs(t['difference']))
                row.append(t['eps^0'])
        rows.append(row)
    assert eta < 1e-9 and pole < 1e-8 and imag < 1e-9 and dq < 1e-6 and dg < 1e-6 and spread < 1e-6
    cols = [('quarks', 'CF'), ('quarks', 'CA'), ('nu=0', 'CA'), ('nu=0', 'nf'), ('nu=2', 'CA'), ('nu=2', 'nf')]
    diffs = []
    for grp, c in cols:
        if grp == 'quarks':
            diffs.append(max(abs(rq['x'][xk][c]['difference']) for xk in xs))
        else:
            diffs.append(max(abs(rg['x'][xk][c][grp]['difference']) for xk in xs))
    t = ['\\begin{table}[t]', '\\centering', '\\small', '\\setlength{\\tabcolsep}{4pt}',
         '\\begin{tabular}{lcccccc}', '\\toprule',
         '& \\multicolumn{2}{c}{$g\\to q\\bar q$} & \\multicolumn{2}{c}{$g\\to gg$, $\\nu=0$} '
         '& \\multicolumn{2}{c}{$g\\to gg$, $\\nu=2$}\\\\',
         '\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}',
         '$x$ & $C_F$ & $C_A$ & $C_A$ & $T_Fn_f$ & $C_A$ & $T_Fn_f$\\\\', '\\midrule']
    for r in rows:
        t.append(r[0] + ' & ' + ' & '.join(f'${v:.8f}$' for v in r[1:]) + '\\\\')
    t += ['\\midrule', 'difference & ' + ' & '.join(f'${sci(v, 1)}$' for v in diffs) + '\\\\',
          '\\bottomrule', '\\end{tabular}',
          '\\caption{The constant $r_0$ of the one-loop kernel from the regions (appendix~\\ref{app:regions}), '
          'coefficient by colour factor, for quark splittings ($gg\\to ZZ$) and for the two components of '
          'gluon splittings ($gg\\to HH$); by the symmetry $x\\to1-x$, $x=0.7$ and $0.9$ give the same. The '
          'last row gives the largest difference from eqs.~(\\ref{eq:r0g}) and~(\\ref{eq:r0q}), which were '
          'identified in the two-loop amplitudes, over the three values of $x$.}',
          '\\label{tab:r0}', '\\end{table}']
    (HERE/'tables'/'table_r0.tex').write_text('\n'.join(t)+'\n')
    # the finite part of J1 is the same for the two helicity structures of g -> q qbar, at x and at 1 - x
    jspread = max(abs(rq['x'][xk][c][f'J1{y}finite part, ek'] - rq['x'][xk][c][f'J1{y}finite part, ke'])
                  for xk in xs for c in ('CF', 'CA') for y in (' ', ' at 1 - x, '))
    assert jspread < 1e-6
    # precision checks at x = 1/2, for both splittings: higher quadrature orders, 16 points on the eps circle
    prec = {}
    for name, res in (('quarks', rq), ('gluons', rg)):
        ch = res['precision checks at x = 1/2 (change of r0)']
        assert set(ch) == {'quadrature orders 96, 72, 40 instead of 64, 48, 28',
                           '16 instead of 32 points on the eps circle'}, (name, set(ch))
        for lab, v in ch.items():
            key = 'quad' if lab.startswith('quadrature') else 'alias'
            prec[key] = max(prec.get(key, 0), max(abs(w) for w in v.values()))
    assert prec['quad'] < 1e-6 and prec['alias'] < 1e-5
    # at x = 1/2 the higher quadrature orders move every colour part towards the closed forms
    lab = 'quadrature orders 96, 72, 40 instead of 64, 48, 28'
    moved = []
    for c in ('CF', 'CA'):
        d = rq['x']['1/2'][c]['difference']
        moved.append((d, d + rq['precision checks at x = 1/2 (change of r0)'][lab][c]))
    for c in ('CA', 'nf'):
        for nu in (0, 2):
            d = rg['x']['1/2'][c][f'nu={nu}']['difference']
            moved.append((d, d + rg['precision checks at x = 1/2 (change of r0)'][lab][f'{c}, nu={nu}']))
    assert all(abs(new) < abs(old) for old, new in moved), moved
    quad_dev = max(abs(new) for _, new in moved)
    macros = {'rzFitHH': sci(fits['HH'], 1), 'rzFitZZ': sci(fits['ZZ'], 1),
              'rzBlindHH': sci(blind_hh, 1), 'rzBlindZZ': sci(blind_zz, 1),
              'rzBlindHHunc': sci(blind_hh_unc, 1), 'rzBlindZZunc': sci(blind_zz_unc, 1),
              'rzRegQ': sci(dq, 1), 'rzRegG': sci(dg, 1), 'rzEta': sci(eta, 1), 'rzPole': sci(pole, 1),
              'rzSpread': sci(spread, 1), 'rzJspread': sci(jspread, 1), 'rzImag': sci(imag, 1),
              'rzQuad': sci(prec['quad'], 1), 'rzAlias': sci(prec['alias'], 1), 'rzQuadDev': sci(quad_dev, 1)}
    vals = [f'{float(Fraction(xk)):.1f}' for xk in xs]
    words = {'rzXvalues': '$x=' + '$, $'.join(vals[:-1]) + f'$ and ${vals[-1]}$'}
    return macros, words


# ---------------------------------------------------------------- Sec. 7.3, Table 11: kappa_1 of planar N = 4 from the regions
def n4reg_part(res, fdh):
    """kappa_1 = zeta_2 from the regions (r0/theory/n4_hexagon.py): checks what the text says, writes Table 11,
    returns the macros."""
    assert res['symbolic']['poles'] == 'cancel' and res['symbolic']['X^2'] == '-1/2' and res['symbolic']['X^1'] == '0'
    assert not res['symbolic']['depends on t1, t2']
    z2 = res['zeta_2']
    pairs = res['pairs']
    fracs = sorted({p['x1'] for p in pairs} | {p['x2'] for p in pairs}, key=Fraction)
    assert len(fracs) == 5 and len(pairs) == 15, (fracs, len(pairs))
    eta = imag = pole = closed = dev = spread_nu = 0.0
    rows = []
    for p in pairs:
        ex = p['expected poles (a_s C_A units)']
        row = [f'{float(Fraction(p["x1"])):.1f}', f'{float(Fraction(p["x2"])):.1f}']
        for nu in ('nu=0', 'nu=2'):
            t = p[nu]
            eta = max(eta, t['max |1/eta coefficient|']); imag = max(imag, t['max |imag|'])
            pole = max(pole, abs(t['eps^-2'] - ex['eps^-2']), abs(t['eps^-1'] - ex['eps^-1']))
            closed = max(closed, abs(t['r0 - [-5 z2 - 2 ln x ln(1-x) terms]']))
            dev = max(dev, abs(t['kappa_1 - zeta_2']))
        spread_nu = max(spread_nu, abs(p['nu=0']['eps^0'] - p['nu=2']['eps^0']))
        row += [p['nu=0']['eps^0'], p['nu=2']['eps^0'], p['nu=0']['kappa_1'], p['nu=2']['kappa_1']]
        rows.append(row)
    assert eta < 1e-9 and imag < 1e-8 and pole < 1e-8 and closed < 1e-6 and dev < 1e-6 and spread_nu < 1e-6
    assert fdh['largest deviation'] < 1e-9 and fdh['x'] == '1/2'
    prec = {}
    for k, v in res.items():
        if k.startswith('change of r0/C_A at x = 1/2'):
            key = 'quad' if 'quadrature' in k else 'alias'
            prec[key] = max(abs(w) for w in v.values())
    t = ['\\begin{table}[t]', '\\centering', '\\small', '\\setlength{\\tabcolsep}{5pt}',
         '\\begin{tabular}{cccccc}', '\\toprule',
         '& & \\multicolumn{2}{c}{$r_0^{\\mathcal N=4}/C_A$} & \\multicolumn{2}{c}{$\\kappa_1$}\\\\',
         '\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}',
         '$x_1$ & $x_2$ & $\\nu=0$ & $\\nu=2$ & $\\nu=0$ & $\\nu=2$\\\\', '\\midrule']
    for r in rows:
        t.append(f'{r[0]} & {r[1]} & ' + ' & '.join(f'${v:.8f}$' for v in r[2:]) + '\\\\')
    t += ['\\midrule', f'eq.~(\\ref{{eq:r0n4}}), $\\zeta_2$ & & \\multicolumn{{2}}{{c}}{{to ${sci(closed, 1)}$}} '
          f'& \\multicolumn{{2}}{{c}}{{${z2:.8f}$, to ${sci(dev, 1)}$}}\\\\',
          '\\bottomrule', '\\end{tabular}',
          '\\caption{Planar $\\mathcal N=4$ super Yang--Mills theory from the regions: the finite part of the '
          'one-loop kernel at $\\mu=\\mu_b$, in units of $a_sC_A$, for the two components and all pairs of the '
          'momentum fractions of the two splittings, and the constant $\\kappa_1$ of eq.~(\\ref{eq:kappaone}) '
          'that follows with the hard functions and the normalisation of ref.~\\cite{Dixon:2016epj}. The last '
          'row gives the largest difference from eq.~(\\ref{eq:r0n4}) and from $\\zeta_2$.}',
          '\\label{tab:n4reg}', '\\end{table}']
    (HERE/'tables'/'table_n4reg.tex').write_text('\n'.join(t)+'\n')
    macros = {'kappaOneReg': sci(dev, 1), 'nFourFDH': sci(fdh['largest deviation'], 1), 'rzNfourPole': sci(pole, 1),
              'rzNfourEta': sci(eta, 1), 'rzNfourImag': sci(imag, 1), 'rzNfourClosed': sci(closed, 1),
              'rzNfourSpread': sci(spread_nu, 1)}
    if prec:
        macros |= {'rzNfourQuad': sci(prec['quad'], 1), 'rzNfourAlias': sci(prec['alias'], 1)}
    words = {'nFourRegPairs': NUMBER_WORDS[len(pairs)] if len(pairs) < len(NUMBER_WORDS) else str(len(pairs)),
             'nFourRegFracs': '$x=' + '$, $'.join(frac_tex(Fraction(f)) for f in fracs[:-1]) + f'$ and ${frac_tex(Fraction(fracs[-1]))}$'}
    return macros, words


# ---------------------------------------------------------------- Paper II Sec. 7.2: the three-loop prediction
def three_loop_part(tl):
    """The three-loop coefficients of r0/theory/three_loop.py: the display equation (tables/eq_threeloop.tex) in the
    general form and the table of the x-dependent coefficients for gg -> HH at mu = m (tables/table_threeloop.tex);
    checks the structure and returns the macros."""
    names = {'Delta_C': sp.Symbol(r'\DD'), 'B': sp.Symbol('B'), 'beta0': sp.Symbol(r'\beta_0'), 'g1': sp.Symbol(r'\Gamma_1'),
             'h': sp.Symbol('h'), 'r0': sp.Symbol('r_0'), 'L_*': sp.Symbol('L_*'), 'zeta3': sp.Symbol(r'\zeta_3')}
    loc = {k: v for k, v in names.items()}
    DD, Bs, b0, G1, hs, r0s, Ls = (names[k] for k in ('Delta_C', 'B', 'beta0', 'g1', 'h', 'r0', 'L_*'))

    def coeff(nu, power):
        rec = tl[nu][f'G^{power}']
        assert not rec['open'], (nu, power, rec['open'])
        return sp.sympify(rec['coefficient'].replace('L_*', 'Lstar'), locals=loc | {'Lstar': Ls})

    c5, c4, c3 = coeff('nu=0', 5), coeff('nu=0', 4), coeff('nu=0', 3)
    d4, d3, d2 = coeff('nu=2', 4), coeff('nu=2', 3), coeff('nu=2', 2)
    # the structure the text states
    assert sp.expand(c5 - sp.Rational(8, 5)*DD**2) == 0 and sp.expand(d4 - 8*DD**2) == 0
    assert sp.expand(c4 - (2*Bs*DD + sp.Rational(4, 3)*DD*b0)) == 0
    assert sp.expand(d3 - (8*Bs*DD - 16*DD**2 + sp.Rational(16, 3)*DD*b0)) == 0
    c4_form = 2*DD*(Bs + sp.Rational(2, 3)*b0)
    c3_form = (sp.Rational(2, 3)*Bs*(Bs + b0) - sp.Rational(1, 3)*DD*G1 - sp.Rational(8, 3)*DD*hs - sp.Rational(4, 3)*DD*r0s
               - sp.Rational(4, 3)*DD*Ls*(Bs + b0 + DD*Ls))
    assert sp.expand(c4 - c4_form) == 0 and sp.expand(c3 - c3_form) == 0
    # the nu = 2 coefficients in terms of the nu = 0 ones (the two transforms of eq. (5.5))
    assert sp.expand(d3 - (4*c4_form - 16*DD**2)) == 0
    assert sp.expand(d2 - (3*c3_form - 4*DD*(3*Bs + 2*b0 - 6*DD))) == 0
    # the next coefficients are open: they contain DG1, beta1, r1_3
    assert set(tl['nu=0']['G^2']['open']) >= {'DG1', 'beta1', 'r1_3'} and set(tl['nu=2']['G^1']['open']) >= {'DG1', 'beta1', 'r1_3'}
    # the component label: three_loop.py evaluates r_0 of each component (r0[0], r0[2]); the display must say so
    eq = ['\\begin{equation}', '\\begin{aligned}',
          '\\M_1^{(3)}&\\supset C_1\\Big[\\tfrac85\\DD^2\\,\\G^5+c_4\\G^4+c_3^{(0)}\\G^3\\Big],\\\\',
          '\\M_2^{(3)}&\\supset K\\Big[8\\DD^2\\,\\G^4+\\big(4c_4-16\\DD^2\\big)\\G^3'
          '+\\big(3c_3^{(2)}-4\\DD(3B+2\\beta_0-6\\DD)\\big)\\G^2\\Big],\\\\',
          'c_4&=2\\DD\\big(B+\\tfrac23\\beta_0\\big),\\\\',
          'c_3^{(\\nu)}&=\\tfrac23B(B+\\beta_0)-\\tfrac13\\DD\\Gamma_1-\\tfrac83\\DD h-\\tfrac43\\DD r_0^{(\\nu)}'
          '-\\tfrac43\\DD L_*\\big(B+\\beta_0+\\DD L_*\\big),',
          '\\end{aligned}', '\\label{eq:threeloop}', '\\end{equation}']
    (HERE/'tables'/'eq_threeloop.tex').write_text('\n'.join(eq) + '\n')
    # the gg -> HH numbers at mu = m
    num = tl['numbers gg -> HH at mu = m']
    xs = sorted(num, key=Fraction)
    assert xs == ['5/8', '3/4', '5/6'], xs
    cols = ['N^2', 'N n_f', 'n_f^2']
    fmt = lambda v: (f'{v[0]:.3f}' + (f'{v[1]:+.3f}\\,i' if abs(v[1]) > 5e-4 else ''))
    t = ['\\begin{table}[t]', '\\centering', '\\small', '\\setlength{\\tabcolsep}{5pt}',
         '\\begin{tabular}{lccc}', '\\toprule', '& $N^2$ & $Nn_f$ & $n_f^2$\\\\', '\\midrule',
         '\\multicolumn{4}{l}{coefficient of $\\G^3$ in $\\M_1^{(3)}/C_1$}\\\\']
    for xk in xs:
        t.append(f'$x={frac_tex(Fraction(xk))}$ & ' + ' & '.join(f'${fmt(num[xk]["nu=0 G^3"][c])}$' for c in cols) + '\\\\')
    t += ['\\midrule', '\\multicolumn{4}{l}{coefficient of $\\G^2$ in $\\M_2^{(3)}/K$}\\\\']
    for xk in xs:
        t.append(f'$x={frac_tex(Fraction(xk))}$ & ' + ' & '.join(f'${fmt(num[xk]["nu=2 G^2"][c])}$' for c in cols) + '\\\\')
    t += ['\\bottomrule', '\\end{tabular}',
          '\\caption{Conditional predictions under the boundary ansatz of sec.~\\ref{sec:threeloop}: the '
          '$x$-dependent coefficients of eq.~(\\ref{eq:threeloop}) for the three-loop $gg\\to HH$ amplitudes with '
          'effective vertices, in minimal subtraction at $\\mu=m_H$, by colour structure, at the three velocities of the '
          'scans; the real part first, the imaginary part marked $i$. The $x$-independent coefficients are those of '
          'eq.~(\\ref{eq:threeloop}) with $\\DD=C_A$ and $B=2\\beta_0$: $\\tfrac85C_A^2$ and '
          '$\\tfrac{176}9C_A^2-\\tfrac{32}9C_An_f$ for $\\G^5$ and $\\G^4$ in $\\M_1^{(3)}/C_1$, '
          '$8C_A^2$ and $\\tfrac{560}9C_A^2-\\tfrac{128}9C_An_f$ for $\\G^4$ and $\\G^3$ in $\\M_2^{(3)}/K$.}',
          '\\label{tab:threeloop}', '\\end{table}']
    (HERE/'tables'/'table_threeloop.tex').write_text('\n'.join(t) + '\n')
    # the x-independent gluon coefficients quoted in the caption, checked against the script's values
    CA, nf = sp.symbols('C_A n_f')
    gl = {k: sp.sympify(v, locals={'C_A': CA, 'n_f': nf}) for k, v in tl['cases']['gluon splittings (gg -> HH)'].items()
          if k in ('nu=0 G^5', 'nu=0 G^4', 'nu=2 G^4', 'nu=2 G^3')}
    assert sp.expand(gl['nu=0 G^5'] - sp.Rational(8, 5)*CA**2) == 0 and sp.expand(gl['nu=2 G^4'] - 8*CA**2) == 0
    assert sp.expand(gl['nu=0 G^4'] - (sp.Rational(176, 9)*CA**2 - sp.Rational(32, 9)*CA*nf)) == 0
    assert sp.expand(gl['nu=2 G^3'] - (sp.Rational(560, 9)*CA**2 - sp.Rational(128, 9)*CA*nf)) == 0
    return {}, {}


# ---------------------------------------------------------------- Sec. 7.1, Table 10: r0 from gg -> HH and gg -> AA
def aa_part(aa, hh, check):
    """The two-loop gg -> AA amplitudes on the forward line (dps/): r0 of both components against
    eq. (7.2) with F_A (nu = 0) and with the opposite sign of the helicity-flip term (nu = 2), next to
    gg -> HH; checks what the text says and writes Table 10."""
    betas = [b for b in sorted(aa['betas'], key=Fraction) if b in ('1/4', '1/2', '2/3')]
    blind = [b for b in sorted(aa['betas'], key=Fraction) if b not in betas]
    assert betas == sorted(hh['betas'], key=Fraction) == ['1/4', '1/2', '2/3'], betas
    cols = [(nu, c) for nu in ('nu0', 'nu2') for c in ('N', 'TFnf')]
    rows, diffs = [], {('HH', k): 0.0 for k in cols} | {('AA', k): 0.0 for k in cols}
    bounds = {('HH', k): 0.0 for k in cols} | {('AA', k): 0.0 for k in cols}
    fit_aa = fit_nu2 = fit_nu0 = imag = vanish = lead = lead0 = oneloop = oneloop_nf = parity = overlap = 0.0
    lead_window_shift = 0.0
    for b in sorted(aa['betas'], key=Fraction):           # the difference from gg -> HH, at every velocity
        a = aa['betas'][b]['r0']['against_hh']
        parity = max(parity, a['nu2']['N'], a['nu2']['TFnf'])
        overlap = max(overlap, a['nu0']['N'], a['nu0']['TFnf'])
    assert parity < 1e-8 and overlap < 1e-8, (parity, overlap)
    for b in betas:
        x = Fraction(1 + Fraction(b), 2)
        r = {'HH': hh['betas'][b]['r0'], 'AA': aa['betas'][b]['r0']}
        fitted = [f'{float(r[pr][nu][c]["value"][0]):.8f}' for nu, c in cols for pr in ('HH', 'AA')]
        closed = [f'{float(r[pr][nu][c]["closed"]):.8f}' for nu, c in cols for pr in ('HH', 'AA')]
        rows.append((f'$x={frac_tex(x)}$', fitted, closed))
        for nu, c in cols:
            for pr in ('HH', 'AA'):
                e = r[pr][nu][c]
                diffs[pr, (nu, c)] = max(diffs[pr, (nu, c)], e['deviation'])
                sh = float(e['shift'][0] if isinstance(e['shift'], list) else e['shift'])
                assert e['deviation'] <= sh < 1e-7, (pr, b, nu, c)
                bounds[pr, (nu, c)] = max(bounds[pr, (nu, c)], sh)
            e = r['AA'][nu][c]
            fit_aa = max(fit_aa, e['deviation']); imag = max(imag, e['imag'])
            if nu == 'nu2':
                fit_nu2 = max(fit_nu2, e['deviation'])
            else:
                fit_nu0 = max(fit_nu0, e['deviation'])
            vanish = max(vanish, r['AA'][nu]['Nf2_over_oneloop'], r['AA'][nu]['Nf_over_N_over_oneloop'])
        for row in aa['betas'][b]['leading_logs']:
            if float(row['prediction']) != 0:
                lead = max(lead, row['deviation'])
                lead_window_shift = max(lead_window_shift, float(row['shift']) / abs(float(row['prediction'])))
            else:
                lead0 = max(lead0, row['deviation'] / abs(float(
                    [q for q in aa['betas'][b]['leading_logs'] if q['component'] == row['component']
                     and q['power'] == row['power'] and q['colour'] == 'N2'][0]['prediction'])))
        for key in ('C1/N', 'K/N'):
            o = aa['betas'][b]['one_loop'][key]
            oneloop = max(oneloop, o['deviation']); oneloop_nf = max(oneloop_nf, o['nf_term_over_fit'])
    assert fit_aa < 1e-8 and imag < 1e-8 and vanish < 1e-8 and lead < 1e-8 and lead0 < 1e-8, (fit_aa, imag, vanish, lead, lead0)
    assert oneloop < 1e-8 and oneloop_nf < 1e-8, (oneloop, oneloop_nf)
    ctrl = max((c['max_rel_diff'] for c in aa['control']), default=0.0)
    assert aa['control'] and ctrl < 1e-20, ctrl
    # Table 10: the two components one above the other, four numbers per row
    t = ['\\begin{table}[t]', '\\centering', '\\footnotesize', '\\setlength{\\tabcolsep}{5pt}',
         '\\begin{tabular}{llcccc}', '\\toprule',
         '& & \\multicolumn{2}{c}{$C_A$} & \\multicolumn{2}{c}{$T_Fn_f$}\\\\',
         '\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}',
         '& & $HH$ & $AA$ & $HH$ & $AA$\\\\', '\\midrule']
    for k, nu in enumerate(('nu0', 'nu2')):
        if k:
            t.append('\\midrule')
        t.append(f'\\multicolumn{{6}}{{l}}{{$\\nu={nu[-1]}$}}\\\\')
        sl = slice(4*k, 4*k + 4)
        for label, fitted, closed in rows:
            t.append(label + ' & amplitudes & ' + ' & '.join(f'${v}$' for v in fitted[sl]) + '\\\\')
            t.append(' & closed form & ' + ' & '.join(f'${v}$' for v in closed[sl]) + '\\\\')
        t.append('\\multicolumn{2}{l}{largest difference} & ' +
                 ' & '.join(f'${sci(diffs[pr, c], 1)}$' for c in cols[2*k:2*k + 2] for pr in ('HH', 'AA')) + '\\\\')
        t.append('\\multicolumn{2}{l}{window shift} & ' +
                 ' & '.join(f'${sci(bounds[pr, c], 1)}$' for c in cols[2*k:2*k + 2] for pr in ('HH', 'AA')) + '\\\\')
    t += ['\\bottomrule', '\\end{tabular}',
          '\\caption{The constant $r_0$ of the one-loop kernel read off the two-loop $gg\\to HH$ and $gg\\to AA$ '
          'amplitudes with eq.~(\\ref{eq:r0read}), for the two components and the two colour factors of '
          'gluon splittings, at the three velocities of the scans, against the closed forms: eq.~(\\ref{eq:r0g}) '
          'with $F_0=F_S$ for $HH$, and eq.~(\\ref{eq:r0aa}) for $AA$, which has $F_0=F_A$ for $\\nu=0$ and the '
          'helicity-flip term with the opposite sign for $\\nu=2$. The two processes differ by the helicity term '
          'alone in both components. The rows marked largest difference give, for each component, the largest '
          'difference between the amplitudes and the closed forms over the three velocities, and the rows marked '
          'window shift the largest change of the fitted value when the fit window is moved one step up in $\\delta$ '
          '(appendix~\\ref{app:fits}).}',
          '\\label{tab:r0amp}', '\\end{table}']
    (HERE/'tables'/'table_r0_amp.tex').write_text('\n'.join(t)+'\n')
    pts = aa['betas'][betas[0]]
    macros = {'rzFitAA': sci(fit_aa, 1), 'rzAAnuTwo': sci(fit_nu2, 1), 'rzAAnuZero': sci(fit_nu0, 1),
              'rzAAimag': sci(imag, 1), 'rzAAvanish': sci(vanish, 1), 'rzAAlead': sci(lead, 1),
              'rzAALogShift': sci(lead_window_shift, 1),
              'rzAAleadVanish': sci(lead0, 1), 'rzAAoneLoop': sci(oneloop, 1), 'rzAAoneLoopNf': sci(oneloop_nf, 1),
              'rzAActrl': sci(ctrl, 1), 'rzBoundAA': sci(max(bounds['AA', c] for c in cols), 1),
              'aaDeltaFirst': f'10^{{-{min(pts["deltas"]["1"])}}}',
              'aaDeltaLast': f'10^{{-{max(pts["deltas"]["1"])}}}', 'rzAAparity': sci(parity, 1),
              'rzAAoverlap': sci(overlap, 1)}
    if blind:
        bdev = max(aa['betas'][b]['r0'][nu][c]['deviation'] for b in blind for nu, c in cols)
        bunc = max(float(aa['betas'][b]['r0'][nu][c]['shift'][0] if isinstance(aa['betas'][b]['r0'][nu][c]['shift'], list)
                         else aa['betas'][b]['r0'][nu][c]['shift']) for b in blind for nu, c in cols)
        assert bdev < 1e-8 and set(blind) == {'1/10', '4/5'}, (bdev, blind)
        macros |= {'rzBlindAA': sci(bdev, 1), 'rzBlindAAunc': sci(bunc, 1),
                   'aaBlindBetas': ','.join(frac_tex(Fraction(b)) for b in blind)}
    c4, c8 = check['g1_b1o2_k4'], check['g1_b1o2_k8']
    assert c4['beta'] == c8['beta'] == '1/2' and c4['delta_exp'] == 4 and c8['delta_exp'] == 8
    words_loss = {'aaLossFour': str(round(c4['loss']['N2'])), 'aaLossEight': str(round(c8['loss']['N2']))}
    words = {'aaPointCount': str(sum(len(v['deltas']['1']) + len(v['deltas']['2']) for v in aa['betas'].values()))}
    words |= words_loss
    return macros, words


# ---------------------------------------------------------------- Table 8: the three processes compared
def rg_coefficients(DC, B):
    """The four coefficients of eq. (5.6): G^3 and G^2 of the Delta J_z = 0 amplitude in units of
    C_1, G^2 and G of the |Delta J_z| = 2 amplitude in units of K, in the variable of Q_*."""
    return [-Fraction(4, 3)*DC, -B, -4*DC, 4*DC-2*B]


def compare_table(colour, zz):
    N, nf = 3, 5
    CA, CF = Fraction(N), Fraction(N*N-1, 2*N)
    beta0 = Fraction(11, 3)*CA-Fraction(2, 3)*nf
    gluon = rg_coefficients(CA, 2*beta0)
    quark = rg_coefficients(2*CF-CA, 6*CF)
    assert set(colour) == set(zz['betas'])
    vel = ','.join(frac_tex(Fraction(b)) for b in sorted(colour, key=Fraction))
    sgn = lambda q: ('-' if q < 0 else '')+frac_tex(abs(q))
    both = lambda a: f'\\multicolumn{{2}}{{c}}{{{a}}}'          # the same entry for AA and HH
    rows = [
        ('final state couples to', both('gluons, effective vertex'), 'quarks, massless loop'),
        ('splitting', both('$g\\to gg$'), '$g\\to q\\bar q$'),
        ('$\\DD$', both('$C_A$'), '$2C_F-C_A$'),
        ('$B$', both('$2\\beta_0$'), '$6C_F$'),
        ('$Q_*^2$', both('$m^4/s$'), '$s^{N_c^2}/m^{2(N_c^2-1)}$'),
        ('variable $\\G_*$', both('$\\G$'), '$\\G_q$'),
    ]
    t = ['\\begin{table}[t]', '\\centering', '\\small', '\\begin{tabular}{lccc}', '\\toprule',
         '& $gg\\to AA$ & $gg\\to HH$ & $gg\\to ZZ$\\\\', '\\midrule']
    t += [' & '.join(r)+'\\\\' for r in rows]
    t += ['\\midrule', '\\multicolumn{4}{l}{one loop}\\\\',
          '$C_1$ & $2\\pi N_cm^2F_A(x)$ & $-2\\pi N_cm^2F_S(x)$ & $|C_1|=16\\pi x^2,\\ 16\\pi(1-x)^2$\\\\',
          f'$K$ & {both("$\\pi N_cm^2W(x)$")} & $|K|=8\\pi x(1-x)$\\\\']
    for head, items in (('two loops, $\\Delta J_z=0$, units of $C_1$', (('$\\G_*^3$', 0), ('$\\G_*^2$', 1))),
                        ('two loops, $|\\Delta J_z|=2$, units of $K$', (('$\\G_*^2$', 2), ('$\\G_*$', 3)))):
        t += ['\\midrule', f'\\multicolumn{{4}}{{l}}{{{head}}}\\\\']
        t += [f'{lab} & {both("$"+sgn(gluon[i])+"$")} & ${sgn(quark[i])}$\\\\' for lab, i in items]
    t += ['\\midrule',
          f'tested in & ref.~\\cite{{P1}}, section~\\ref{{sec:r0}} & section~\\ref{{sec:hh}} & section~\\ref{{sec:zz}}\\\\',
          '\\bottomrule', '\\end{tabular}',
          '\\caption{The DPS terms in the three processes. The two-loop coefficients follow from eq.~(\\ref{eq:rg}) '
          'with the colour factors above, in the variable $\\G_*=\\ln(p_T/Q_*)+i\\pi/2$ of the scale $Q_*$, '
          'which is $\\G$ of eq.~(\\ref{eq:Gdef}) for gluon splittings and $\\G_q$ of eq.~(\\ref{eq:Gq}) for '
          f'quark splittings, for $N_c={N}$ and $n_f={nf}$. '
          'For $gg\\to ZZ$ the one-loop coefficients are moduli in the normalisation of section~\\ref{sec:zz}, '
          'and the amplitudes with $\\Delta J_z=0$ and $|\\Delta J_z|=2$ are those of the first and second '
          'group there; for $gg\\to AA$ and $gg\\to HH$ they are $\\M_1$ and $\\M_2$. Ref.~\\cite{P1} tested the '
          '$gg\\to AA$ relations at threshold and, for the $\\G^3$ and $\\G^2$ coefficients, in the forward limit; '
          'section~\\ref{sec:r0} tests the constant of the boundary term there. The $gg\\to HH$ and '
          f'$gg\\to ZZ$ amplitudes are tested at $\\beta={vel}$.}}',
          '\\label{tab:compare}', '\\end{table}']
    (HERE/'tables'/'table_compare.tex').write_text('\n'.join(t)+'\n')

def glauber_part(g):
    """The imaginary part of the kernel (r0/theory/glauber.py, task 5): checks what Secs. 7.2, 7.3 and App. F say
    and returns the macros."""
    b = g['b_space']
    assert b['Im K/K0 in units of 2 pi i Dc']['eps^-1'] == '1' and b['finite part at mu = mu_b'] == 0
    assert b['pole equals that of 2 I_K'] and b['X_mu term equals Im of the RG bracket']
    pair = g['off_shell_eikonal_pair']
    assert pair['worst rel diff'] < 1e-9
    assert all(r['phase relative to the log'] == ('ln(A/sigma1 sigma2) - i pi' if r['sigma1'] > 0 and r['sigma2'] > 0
                                                  else 'none (real)') for r in pair['rows'])
    r = g['eta_regulated_regions']                   # part 2b: the four regions of the off-shell exchange with the regulator
    assert r['worst rel diff of the sum from J at A/|s1 s2| = 5e5'] < 1e-5 and r['worst rel deviation of G from 2 pi^2/A'] < 1e-8
    assert r['largest real part of C_n + C_nbar + S (s2 > 0)'] < 1e-20 and r['worst rel deviation of the 3 -> 3 cancellation'] < 1e-5
    assert len(r['cases']) == 4 and {c['s2'] > 0 for c in r['cases']} == {True, False}
    ratio = max(c['A/|s1 s2|'] for c in r['cases'])
    assert ratio == 5e5
    m = g['explicit_two_loop_graph']
    graph = m['worst rel deviation of Re f0 from -(pi/2) Im f1']
    assert graph < 1e-3 and len(m['rows']) == 18
    assert all(abs(row['Glauber region from the leading-power integrand']/row['-(pi/2) Im f1'] - 1) < 1e-12 for row in m['rows'])
    assert len({(r['x'], tuple(r['khat'])) for r in m['rows']}) == 3 and len({tuple(r['lhat']) for r in m['rows']}) == 6
    q = g['qcd_imaginary_parts']
    assert q['largest |Im r0|'] < 1e-8 and q['largest |Im r0|/uncertainty'] < 1
    procs = q['processes']
    assert {'gg -> HH', 'gg -> ZZ'} <= set(procs)
    h = g['n4_two_to_four']
    assert h['poles cancel with the kernel L_*/eps, L_* = -X_r + i pi'] and not h['depends on t1, t2 or v']
    assert h['c_G from the single logarithm'] == '0'
    de = h['DE coefficients of ln^k delta in E^(2) (H^v_w = H_w(1-v))']
    assert sp.simplify(sp.sympify(de['1']) + 2*sp.I*sp.pi**3/3) == 0 and sp.simplify(sp.sympify(de['2']) + sp.pi**2/2) == 0
    names = {'gg -> HH': '$gg\\to HH$', 'gg -> ZZ': '$gg\\to ZZ$', 'gg -> AA': '$gg\\to AA$'}
    order = [p for p in ('gg -> HH', 'gg -> AA', 'gg -> ZZ') if p in procs]
    plist = ', '.join(names[p] for p in order[:-1]) + ' and ' + names[order[-1]]
    macros = {'glauberGraph': sci(graph, 1), 'glauberImag': sci(q['largest |Im r0|'], 1),
              'glauberImagRatio': f"{q['largest |Im r0|/uncertainty']:.2f}", 'glauberPair': sci(pair['worst rel diff'], 1),
              'glauberSum': sci(r['worst rel diff of the sum from J at A/|s1 s2| = 5e5'], 1),
              'glauberRatio': sci(ratio, 0), 'glauberG': sci(r['worst rel deviation of G from 2 pi^2/A'], 1),
              'glauberCancel': sci(r['worst rel deviation of the 3 -> 3 cancellation'], 1)}
    words = {'glauberConfigs': NUMBER_WORDS[len(m['rows'])], 'glauberProcs': plist,
             'glauberValues': str(len(q['rows']))}
    return macros, words




def unequal_matching_part(data):
    """New fits of retained unequal-mass amplitudes, including both helicities."""
    results=data['results']
    assert len(results)==2 and all(set(r['helicities'])=={'++++','++--'} for r in results)
    macros={}
    for i, r in enumerate(results,1):
        leading=float(r['leading_log_max_deviation'])
        p2=float(r['leading_log_P2_change'])
        window=float(r['leading_log_window_shift'])
        matching=float(r['finite_matching_max_deviation'])
        c1=max(float(h['c1_absolute_deviation']) for h in r['helicities'].values())
        assert leading<1e-8 and p2<1e-6 and window<1e-7 and matching<1e-8 and c1<1e-12
        suffix=('One','Two')[i-1]
        macros.update({f'uvLead{suffix}':sci(leading,1),f'uvPtwo{suffix}':sci(p2,1),
                       f'uvWindow{suffix}':sci(window,1),f'uvMatching{suffix}':sci(matching,1),
                       f'uvCOne{suffix}':sci(c1,1)})
    return macros

def main():
    colour, oneloop, scale, forward = load('colour'), load('one_loop'), load('scale'), load('forward')
    stage2a, p1, bose, scan = load('stage2a'), load('p1'), load('bose'), scan_summary()
    BETAS = sorted(colour, key=Fraction)
    rows_all = [r for b in BETAS for r in colour[b]['rows']]
    nonzero = [r for r in rows_all if float(r['prediction']) != 0]

    # ------------------------------------------------------------ numbers quoted in the text
    lead_dev = max(float(r['deviation']) for r in nonzero)
    vanish = max(float(r['deviation']) for r in rows_all if float(r['prediction']) == 0)
    shift = max(float(r['window_shift']) for r in rows_all if r['window_shift'] != 'nan')
    plain_max = max(float(r['plain_fit_deviation']) for r in nonzero)
    ol_slope = max(float(oneloop[b]['one_loop']['spin0_slope_rel_dev']) for b in BETAS)
    ol_K = max(float(oneloop[b]['one_loop']['spin2_forward_rel_dev']) for b in BETAS)
    ol_shift = max(float(oneloop[b]['one_loop']['spin0_slope_shift']) for b in BETAS)   # fit bound (window shift)
    tree = max(float(oneloop[b]['tree_max_rel_diff']) for b in BETAS)
    absent = max(float(v) for d in scale['absent_terms'].values() for c in d.values() for v in c.values())
    # the spin-0 single log moves with mu_R by design; its slope is checked separately
    spread = max(float(v['spread']) for b in scale['three_scales'] for k, v in scale['three_scales'][b].items()
                 if not k.startswith('spin0 G^1'))
    slope_dev = max(float(scale['three_scales'][b]['spin0 G^1 slope in L']['rel_dev']) for b in scale['three_scales'])
    stored = float(scale['stored_points']['max_rel_diff'])
    prec = max(float(c['rel_diff']) for c in scale['precision_control'])
    ctrl = scale['precision_control'][0]
    assert all(c['delta_exp'] == ctrl['delta_exp'] and c['beta'] == ctrl['beta'] and c['dps'] == ctrl['dps']
               for c in scale['precision_control']), 'the text describes one control setting'
    kdig = min(int(v['constant_digits']) for v in forward.values())
    # constant_digits refers to the real part of (M - K)/p_T^2; the text does not say so, so check
    # that the complex ratio is constant to as many digits, over the same points
    for v in forward.values():
        sp, n = ratio_spread(v['curve'], 6)
        assert digits(sp) >= kdig and n == len(v['r_values']), 'forward_approach.py stores r for delta <= 10^-6'
    kslope = max(float(v['dr_dpT2_spread']) for v in forward.values())
    table2_beta = preferred_beta(BETAS)
    # Appendix B: imaginary parts of the fitted leading coefficients (real in the prediction),
    # and one power of G more than expected, both relative to the leading N^2 coefficient
    imag = max(abs(float(r['im']))/abs(float(r['prediction'])) for r in nonzero)

    def lead(b, coeff):
        return abs(float([r for r in colour[b]['rows'] if r['coefficient'] == coeff and r['colour'] == 'N2'][0]['prediction']))
    extra = max(float(v)/lead(b, 'spin0 G^3' if k.startswith('spin0') else 'spin2 G^2')
                for b in BETAS for k, v in colour[b]['extra_powers'].items())
    pT2_lo, pT2_hi = scan['pT2_range']
    hhli = hh_li22_check()
    assert hhli['max_in_working_precision'] < 1, 'the text says that Li_{2,2} agrees within the working precision'

    macros = {
        # task B
        'devLeadMax': up(lead_dev), 'vanishMax': up(vanish), 'shiftMax': up(shift),
        'plainDevMax': sci(plain_max, 0),                      # the text says "about"
        'colourTableBeta': frac_tex(Fraction(table2_beta)),
        'imagMax': up(imag), 'extraPowerMax': up(extra),
        # task C
        'devOneLoopMax': up(max(ol_slope, ol_K)), 'devKMax': up(ol_K), 'treeMax': up(tree),
        'oneLoopShiftMax': up(ol_shift),
        # task D
        'absentMax': up(absent), 'scaleSpread': up(spread), 'scaleSlopeDev': up(slope_dev),
        'storedMaxDiff': up(stored), 'precisionCtrl': up(prec),
        'ctrlDelta': f'10^{{-{ctrl["delta_exp"]}}}', 'ctrlBeta': frac_tex(Fraction(ctrl['beta'])),
        'ctrlExtraDigits': str(ctrl['dps'][1]-ctrl['dps'][0]),
        # task F
        'kOneSlopeSpread': up(kslope),
        # scan settings, from the points themselves
        'scanBetas': ','.join(frac_tex(Fraction(b)) for b in scan['betas']),
        'scanDeltaFirst': f'10^{{-{min(scan["delta_exponents"])}}}',
        'scanDeltaLast': f'10^{{-{max(scan["delta_exponents"])}}}',
        'digitsBase': str(scan['digits']['base']), 'digitsPerDecade': str(scan['digits']['per_decade']),
        'pTsqMin': sci(pT2_lo, 1), 'pTsqMax': sci(pT2_hi, 1),
        # quoted from P1 and from the earlier package (data/p1_numbers.json, data/bose_symmetry.json)
        'pOneKoneDev': sci(p1['AA_K1']['value'], 0),
        'pOneAALogDev': sci(p1['AA_log']['value'], 0),
        'pOneHHLogDev': sci(p1['HH_log']['value'], 0),
        'boseDevMone': up(bose['1']['rel']), 'boseDevMtwo': up(bose['2']['rel']),
        # App. A: Li_{2,2} against eq. (li22sum); the text says "within the working precision"
        'hhLiValueDev': up(hhli['max_rel_diff']),
    }
    # counts the text spells out as words (not math)
    words = {
        'hhLiGuard': str(hhli['guard']),
        'velocityCount': NUMBER_WORDS[len(BETAS)],
        'lValueCount': NUMBER_WORDS[len(scan['L_values'])],
        'kOneApproachDecades': NUMBER_WORDS[min(len(v['r_values']) for v in forward.values())-1],
        'storedCount': str(scale['stored_points']['compared']),
        'kOneApproachDigits': NUMBER_WORDS[kdig],
    }
    # gg -> ZZ (task J): section 6.2, Table 7 and the gg -> ZZ parts of appendices A and B
    (HERE/'tables').mkdir(exist_ok=True)
    zz_macros, zz_words = zz_part(load('zz'), load('zz_summary'), load('zz_forward'))
    assert len({v['spread_from_delta_exp'] for v in load('zz_forward').values()}) == 1
    compare_table(colour, load('zz'))
    assert not set(zz_macros) & set(macros) and not set(zz_words) & set(words)
    macros |= zz_macros
    macros |= windows_part(load('windows'), load('zz_windows'), colour, load('zz'))
    words |= zz_words
    fwd_macros, fwd_words = forward_part(load('aa_curve'), load('zz_wide'), load('zz'))
    macros |= fwd_macros
    words |= fwd_words
    n4_macros, n4_words = n4_part(load('n4'))
    macros |= n4_macros
    words |= n4_words
    k_macros, k_words = kernel_part(load('kernel'))
    assert not set(k_macros) & set(macros) and not set(k_words) & set(words)
    macros |= k_macros
    words |= k_words
    (HERE/'tables').mkdir(exist_ok=True)
    r_macros, r_words = r0_part(load('r0_fits'), load('r0_blind_hh'), load('r0_blind_zz'), load('r0_quarks'),
                                load('r0_gluons'))
    assert not set(r_macros) & set(macros) and not set(r_words) & set(words)
    macros |= r_macros
    words |= r_words
    uv_macros=unequal_matching_part(load('uv_matching'))
    assert not set(uv_macros) & set(macros)
    macros |= uv_macros
    f_macros, f_words = r0_fits_part(load('r0_scan_hh'), load('r0_scan_zz'), load('r0_blind_hh'), load('r0_blind_zz'))
    assert not set(f_macros) & set(macros) and not set(f_words) & set(words)
    macros |= f_macros
    words |= f_words
    nr_macros, nr_words = n4reg_part(load('n4reg'), load('n4fdh'))
    assert not set(nr_macros) & set(macros) and not set(nr_words) & set(words)
    macros |= nr_macros
    words |= nr_words
    g_macros, g_words = glauber_part(load('glauber'))
    assert not set(g_macros) & set(macros) and not set(g_words) & set(words)
    macros |= g_macros
    words |= g_words
    three_loop_part(load('three_loop'))
    for more in (nmhv_primary_part(load('n4_primary'), load('n4')), ultrasoft_part(load('ultrasoft')),
                 ultrasoft_dimension_part(SOURCES['ultrasoft_d'], load('ultrasoft'))):
        assert not set(more) & set(macros) and not set(more) & set(words)
        words |= more
    cs_macros, cs_words = collins_soper_part(load('collins_soper'), load('n4'))
    assert not set(cs_macros) & set(macros) and not set(cs_words) & set(words)
    macros |= cs_macros
    words |= cs_words
    if all(SOURCES[k].exists() for k in ('dps_aa', 'dps_hh', 'dps_aa_check')):
        a_macros, a_words = aa_part(load('dps_aa'), load('dps_hh'), load('dps_aa_check'))
        assert not set(a_macros) & set(macros) and not set(a_words) & set(words)
        macros |= a_macros
        words |= a_words
    else:
        print('gg -> AA results of dps/ not present: Table 10 and the gg -> AA passages are left out '
              '(the LaTeX source tests \\ifdefined\\rzFitAA)')
    with open(HERE/'numbers.tex', 'w') as f:
        f.write('% generated by make_tables.py from the result files; do not edit\n')
        for k, v in macros.items():
            f.write(f'\\newcommand{{\\{k}}}{{\\ensuremath{{{v}}}}}\n')
        for k, v in words.items():
            f.write(f'\\newcommand{{\\{k}}}{{{v}}}\n')
    (HERE/'tables').mkdir(exist_ok=True)

    # ------------------------------------------------------------ Table 1: one loop
    rng = lambda e: f'${p1[e]["beta_range"][0]}\\le\\beta\\le{p1[e]["beta_range"][1]}$'
    rows = [
        ('$gg\\to AA$, $\\ln p_T^2$ in $\\M_1^{(1)}$', '$\\pi N_cm^2F_A(x)$',
         'exact one-loop remainders \\cite{P1}', f'${sci(p1["AA_log"]["value"], 0)}$'),
        ('$gg\\to AA$, forward limit of $\\M_2^{(1)}$', '$m^2K_1(x)$',
         'exact one-loop remainders \\cite{P1}', f'${sci(p1["AA_K1"]["value"], 0)}$'),
        ('$gg\\to HH$, $\\ln p_T^2$ in $\\M_1^{(1)}$', '$-\\pi N_cm^2F_S(x)$',
         'one-loop amplitude of ref.~\\cite{P1}', f'${sci(p1["HH_log"]["value"], 0)}$'),
        ('', '', 'scale dependence of the two-loop amplitudes, section~\\ref{sec:hhoneloop}', f'${sci(ol_slope, 1)}$'),
        ('$gg\\to HH$, forward limit of $\\M_2^{(1)}$', '$m^2K_1(x)$',
         'scale dependence of the two-loop amplitudes, section~\\ref{sec:hhoneloop}', f'${sci(ol_K, 1)}$'),
    ]
    t = ['\\begin{table}[t]', '\\centering', '\\small',
         '\\begin{tabular}{p{4.2cm}p{2.4cm}p{4.6cm}c}', '\\toprule',
         'quantity & derived & compared with & rel.\\ difference\\\\', '\\midrule']
    t += [' & '.join(r)+'\\\\' for r in rows]
    t += ['\\bottomrule', '\\end{tabular}',
          '\\caption{One-loop DPS terms, eqs.~(\\ref{eq:FA})--(\\ref{eq:K1}), compared with exact amplitudes. '
          f'The $AA$ comparisons are those of ref.~\\cite{{P1}}, for {rng("AA_log")} ($F_A$) and '
          f'{rng("AA_K1")} ($K_1$). The comparison with the one-loop $gg\\to HH$ amplitude of ref.~\\cite{{P1}} '
          f'covers {rng("HH_log")}, and those based on the two-loop amplitudes are at '
          f'$\\beta={",".join(BETAS)}$.}}', '\\label{tab:oneloop}', '\\end{table}']
    (HERE/'tables'/'table_oneloop.tex').write_text('\n'.join(t)+'\n')

    # ------------------------------------------------------------ Table 2: colour terms at one velocity
    names = {'N2': '$N^2$', 'N_Nf': '$Nn_f$', 'Nf_over_N': '$n_f/N$', 'Nf2': '$n_f^2$'}
    coeffs = {'spin0 G^3': '$c_{13}$, $\\G^3$ in $A_1$', 'spin0 G^2': '$c_{12}$, $\\G^2$ in $A_1$',
              'spin2 G^2': '$c_{22}$, $\\G^2$ in $A_2$', 'spin2 G^1': '$c_{21}$, $\\G$ in $A_2$'}
    b = table2_beta
    t = ['\\begin{table}[t]', '\\centering', '\\small',
         '\\begin{tabular}{llrrc}', '\\toprule',
         'coefficient & colour factor & prediction & extracted & deviation\\\\']
    for ck, cname in coeffs.items():
        t.append('\\midrule')
        for j, col in enumerate(('N2', 'N_Nf', 'Nf_over_N', 'Nf2')):
            r = [x for x in colour[b]['rows'] if x['coefficient'] == ck and x['colour'] == col][0]
            pred = float(r['prediction'])
            mod = math.hypot(float(r['re']), float(r['im']))
            dev = float(r['deviation'])
            ps = f'${pred:.10f}$' if pred != 0 else '$0$'          # math mode, so that minus signs are minus signs
            es = f'${float(r["re"]):.10f}$' if pred != 0 else f'${sci(mod, 1)}$'
            ds = f'${sci(dev, 1)}$' if pred != 0 else f'(${sci(dev, 1)}$)'
            t.append(f'{cname if j == 0 else ""} & {names[col]} & {ps} & {es} & {ds}\\\\')
    t += ['\\bottomrule', '\\end{tabular}',
          f'\\caption{{Coefficients of the leading logarithms of the two-loop $gg\\to HH$ remainders in the '
          f'forward limit, by colour factor, at $\\beta={b}$ and $\\mu_R=m_H$, compared with eq.~(\\ref{{eq:colourpred}}). The deviations are '
          'relative. For terms predicted to vanish, the modulus of the extracted coefficient is given, '
          'with its size relative to the $N^2$ coefficient in brackets. Moving the fit window one step '
          f'up in $\\delta$ changes the coefficients by at most ${up(shift)}$ of the leading coefficient '
          f'at all {NUMBER_WORDS[len(BETAS)]} velocities.}}',
          '\\label{tab:colour}', '\\end{table}']
    (HERE/'tables'/'table_colour.tex').write_text('\n'.join(t)+'\n')

    # ------------------------------------------------------------ Table 3: one loop from the scale dependence
    t = ['\\begin{table}[t]', '\\centering', '\\small',
         '\\begin{tabular}{ll' + 'c'*len(BETAS) + '}', '\\toprule',
         'quantity & expected & ' + ' & '.join(f'$\\beta={x}$' for x in BETAS) + '\\\\', '\\midrule']
    ol = lambda key: [float(oneloop[x]['one_loop'][key]) for x in BETAS]
    t.append('$\\bar c_1$, coefficient of $N\\G$ in $A_1^{(1)}$ & $-2\\pi F_S$ & '
             + ' & '.join(f'${sci(v, 1)}$' for v in ol('spin0_slope_rel_dev')) + '\\\\')
    t.append('$\\bar K_1$, coefficient of $N$ in $A_2^{(1)}$ at $p_T=0$ & $\\pi W$ & '
             + ' & '.join(f'${sci(v, 1)}$' for v in ol('spin2_forward_rel_dev')) + '\\\\')
    t.append('coefficient of $n_f\\G$ in $A_1^{(1)}$ & 0 & '
             + ' & '.join(f'${sci(v, 1, 1e-30)}$' for v in ol('spin0_slope_mF_abs')) + '\\\\')
    labels = {'c13(N^2)/slope': ('$c_{13}(N^2)/\\bar c_1$', '$-4/3$'), 'c12(N^2)/slope': ('$c_{12}(N^2)/\\bar c_1$', '$-22/3$'),
              'c12(N n_f)/slope': ('$c_{12}(Nn_f)/\\bar c_1$', '$4/3$'), 'c22(N^2)/K': ('$c_{22}(N^2)/\\bar K_1$', '$-4$'),
              'c21(N^2)/K': ('$c_{21}(N^2)/\\bar K_1$', '$-32/3$'), 'c21(N n_f)/K': ('$c_{21}(Nn_f)/\\bar K_1$', '$8/3$')}
    for k, (lab, ex) in labels.items():
        t.append(f'{lab} & {ex} & '
                 + ' & '.join(f'${sci(float(oneloop[x]["ratios"][k]["rel_dev"]), 1)}$' for x in BETAS) + '\\\\')
    t += ['\\bottomrule', '\\end{tabular}',
          '\\caption{One-loop $gg\\to HH$ coefficients obtained from the coefficient of $L$ in the two-loop '
          'remainders, and the relations of eq.~(\\ref{eq:relations}) tested with numbers taken from the same '
          'expressions. Here $\\bar c_1=c_1/N$ and $\\bar K_1=K_1/N$ are the one-loop coefficients per unit $N$, '
          'and $c_{ij}(N^2)$ and $c_{ij}(Nn_f)$ are the coefficients of $N^2$ and $Nn_f$ in $c_{ij}$. The '
          'entries are relative deviations from the expected values, except in the third row, which gives '
          'the absolute value.}',
          '\\label{tab:hhoneloop}', '\\end{table}']
    (HERE/'tables'/'table_hhoneloop.tex').write_text('\n'.join(t)+'\n')

    # ------------------------------------------------------------ Table 4: what the subtraction adds (task A)
    row_names = {'spin 0, C_1 G^3': '$\\G^3$ in $\\M_1$, units of $C_1$',
                 'spin 0, C_1 beta0 G^2': '$\\beta_0\\G^2$ in $\\M_1$, units of $C_1$',
                 'spin 2, K G^2': '$\\G^2$ in $\\M_2$, units of $K$',
                 'spin 2, K beta0 G': '$\\beta_0\\G$ in $\\M_2$, units of $K$'}
    assert stage2a['same_in_both_eps_conventions']
    t = ['\\begin{table}[t]', '\\centering', '\\small',
         '\\begin{tabular}{lccc}', '\\toprule',
         'coefficient & finite & added by renormalisation & two-loop amplitude\\\\',
         '& remainder & and $-2I^{(1)}\\M^{(1)}$ & at $O(\\epsilon^0)$\\\\', '\\midrule']
    for r in stage2a['rows']:
        t.append(f'{row_names[r["coefficient"]]} & ${coeff_tex(r["finite_remainder"])}$ & '
                 f'${coeff_tex(r["added_by_uv_and_catani"])}$ & ${coeff_tex(r["two_loop_amplitude"])}$\\\\')
    t += ['\\bottomrule', '\\end{tabular}',
          '\\caption{Leading coefficients of the two-loop finite remainder, the contribution of '
          'ultraviolet renormalisation and infrared subtraction acting on the one-loop DPS term, and the '
          'resulting coefficients of the unrenormalised, unsubtracted two-loop amplitude at '
          '$O(\\epsilon^0)$.}', '\\label{tab:stage2a}', '\\end{table}']
    (HERE/'tables'/'table_stage2a.tex').write_text('\n'.join(t)+'\n')

    print('numbers.tex:')
    for k, v in {**macros, **words}.items():
        print(f'  {k:20s} {v}')
    print('tables: ' + ', '.join(f'table_{n}.tex' for n in ('oneloop', 'colour', 'hhoneloop', 'stage2a', 'zz', 'compare', 'n4', 'n4nmhv', 'r0', 'n4reg')))


if __name__ == '__main__':
    main()
