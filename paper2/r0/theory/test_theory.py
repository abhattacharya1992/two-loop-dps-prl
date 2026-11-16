"""Recombine the stored Laurent coefficients of the collinear functions into r0 (quark and gluon splittings) and
check that the stored results r0_theory_quarks.json and r0_theory_gluons.json are reproduced, and that the
statements of the paper hold: 1/eta poles cancel, poles as expected, finite parts equal to the closed forms
found in the data, helicity independence for quarks.  Then the same for planar N = 4 (kappa1_n4.json from the
stored n4_laurent_x*.json and gg_laurent_x*.json): poles 2/eps^2 + 2L_*/eps, kappa_1 = zeta_2 for every pair of
momentum fractions and both components, and the scheme check fdh_check_x1o2.json.  Then the imaginary part of
the kernel (glauber.py -> glauber.json, workbook VII.15): the pole 2 pi i Dc/eps with no finite part at mu = mu_b,
the four eta-regulated regions of the off-shell exchange (glauber_regions.py: their sum is the exact exchange and the
Glauber region alone its -i pi term), the explicit two-loop graph, Im r0 = 0 in the data, and planar N = 4 in 2 -> 4
kinematics (c_G = 0, no v).

    python3 test_theory.py        (about fifteen seconds; no recomputation of the collinear functions)
"""
import json, os, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TOL = 1e-12


def walk(a, b, path=''):
    """largest numerical difference between two nested JSON structures with the same keys"""
    if isinstance(a, dict):
        assert set(a) == set(b), (path, set(a) ^ set(b))
        return max([walk(a[k], b[k], path + '/' + k) for k in a] or [0])
    if isinstance(a, list):
        assert len(a) == len(b), path
        return max([walk(x, y, path) for x, y in zip(a, b)] or [0])
    if isinstance(a, (int, float)):
        return abs(a - b)
    assert a == b, (path, a, b)
    return 0


def walk_rel(a, b, path=''):
    """largest relative difference between two nested JSON structures with the same keys (strings must agree)"""
    if isinstance(a, dict):
        assert set(a) == set(b), (path, set(a) ^ set(b))
        return max([walk_rel(a[k], b[k], path + '/' + k) for k in a] or [0])
    if isinstance(a, (list, tuple)):
        assert len(a) == len(b), path
        return max([walk_rel(x, y, path) for x, y in zip(a, b)] or [0])
    if isinstance(a, bool) or a is None:
        assert a == b, (path, a, b)
        return 0
    if isinstance(a, (int, float)):
        return abs(a - b)/max(abs(a), abs(b), 1e-300) if (a != 0 or b != 0) else 0
    if a != b and isinstance(a, str):
        assert sp_equal(a, b), (path, a, b)              # sympy prints expressions differently across versions
    else:
        assert a == b, (path, a, b)
    return 0


def sp_equal(x, y):
    import sympy as sp
    try:
        return sp.simplify(sp.sympify(x) - sp.sympify(y)) == 0
    except Exception:
        return False


def main():
    tmp = tempfile.mkdtemp()
    try:
        for name in ('r0_theory_quarks', 'r0_theory_gluons'):
            stored = json.load(open(os.path.join(HERE, name + '.json')))
            shutil.copy(os.path.join(HERE, name + '.json'), os.path.join(tmp, name + '.json'))
            subprocess.run([sys.executable, name + '.py'], cwd=HERE, check=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
            new = json.load(open(os.path.join(HERE, name + '.json')))
            d = walk(stored, new)
            print(f'{name}.json reproduced, largest difference {d:.1e}')
            assert d < TOL
        q = json.load(open(os.path.join(HERE, 'r0_theory_quarks.json')))
        g = json.load(open(os.path.join(HERE, 'r0_theory_gluons.json')))
        worst = 0
        for xk, rec in q['x'].items():
            for c in ('CF', 'CA'):
                t = rec[c]['theory, J1 structure ek, J2 structure ek']
                assert t['max |1/eta coefficient|'] < 1e-9 and t['max |imag|'] < 1e-9
                ex = rec[c]['expected poles']
                assert abs(t['eps^-2'] - ex['eps^-2']) < 1e-8 and abs(t['eps^-1'] - ex['eps^-1']) < 1e-8
                assert rec[c]['max spread of eps^0 over the four helicity configurations'] < 1e-6
                worst = max(worst, abs(rec[c]['difference']))
        for xk, rec in g['x'].items():
            for c in ('CA', 'nf'):
                for nu, t in rec[c].items():
                    assert t['max |1/eta coefficient|'] < 1e-9 and t['max |imag|'] < 1e-9
                    ex = t['expected poles']
                    assert abs(t['eps^-2'] - ex['eps^-2']) < 1e-8 and abs(t['eps^-1'] - ex['eps^-1']) < 1e-8
                    worst = max(worst, abs(t['difference']))
        assert worst < 1e-6
        print(f'poles, 1/eta cancellation and helicity independence as stated; r0 equals the closed forms to {worst:.1e}')
        # planar N = 4: kappa_1 from the regions
        name = 'kappa1_n4'
        if os.path.exists(os.path.join(HERE, name + '.json')):
            stored = json.load(open(os.path.join(HERE, name + '.json')))
            shutil.copy(os.path.join(HERE, name + '.json'), os.path.join(tmp, name + '.json'))
            subprocess.run([sys.executable, 'n4_hexagon.py'], cwd=HERE, check=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
            new = json.load(open(os.path.join(HERE, name + '.json')))
            d = walk(stored, new)
            print(f'{name}.json reproduced, largest difference {d:.1e}')
            assert d < TOL
            assert new['symbolic']['poles'] == 'cancel' and new['symbolic']['X^2'] == '-1/2' and new['symbolic']['X^1'] == '0'
            z2 = new['zeta_2']
            assert len(new['pairs']) >= 15, len(new['pairs'])
            for rec in new['pairs']:
                ex = rec['expected poles (a_s C_A units)']
                for nu in ('nu=0', 'nu=2'):
                    t = rec[nu]
                    assert t['max |1/eta coefficient|'] < 1e-9 and t['max |imag|'] < 1e-8
                    assert abs(t['eps^-2'] - ex['eps^-2']) < 1e-8 and abs(t['eps^-1'] - ex['eps^-1']) < 1e-8
                    assert abs(t['kappa_1'] - z2) < 1e-6, (rec['x1'], rec['x2'], nu, t['kappa_1'] - z2)
            fdh = json.load(open(os.path.join(HERE, 'fdh_check_x1o2.json')))
            assert fdh['largest deviation'] < 1e-9
            print(f'planar N = 4: poles as stated, kappa_1 = zeta_2 to {new["largest |kappa_1 - zeta_2|"]:.1e} for '
                  f'{len(new["pairs"])} pairs of fractions and both components; FDH = CDR + eps x scalar to '
                  f'{fdh["largest deviation"]:.1e}')
        else:
            print('kappa1_n4.json not present (make n4)')
        # the imaginary part of the kernel: the Coulomb phase and the Glauber region (glauber.py, workbook VII.15)
        name = 'glauber'
        stored = json.load(open(os.path.join(HERE, name + '.json')))
        shutil.copy(os.path.join(HERE, name + '.json'), os.path.join(tmp, name + '.json'))
        subprocess.run([sys.executable, 'glauber.py'], cwd=HERE, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        new = json.load(open(os.path.join(HERE, name + '.json')))
        # the data part depends on which result files exist (gg -> AA once common/dps/results_aa.json is there)
        d = walk_rel({k: v for k, v in stored.items() if k != 'qcd_imaginary_parts'},
                     {k: v for k, v in new.items() if k != 'qcd_imaginary_parts'})
        print(f'{name}.json reproduced, largest relative difference {d:.1e}')
        assert d < 1e-6                 # the model graph is double precision: the last digits may differ between machines
        g = new['b_space']
        assert g['Im K/K0 in units of 2 pi i Dc']['eps^-1'] == '1' and g['finite part at mu = mu_b'] == 0
        assert g['pole equals that of 2 I_K'] and g['X_mu term equals Im of the RG bracket']
        assert new['off_shell_eikonal_pair']['worst rel diff'] < 1e-9
        assert all(r['phase relative to the log'] == ('ln(A/sigma1 sigma2) - i pi' if r['sigma1'] > 0 and r['sigma2'] > 0 else 'none (real)')
                   for r in new['off_shell_eikonal_pair']['rows'])
        r = new['eta_regulated_regions']
        assert r['worst rel diff of the sum from J at A/|s1 s2| = 5e5'] < 1e-5 and r['worst rel deviation of G from 2 pi^2/A'] < 1e-8
        assert r['largest real part of C_n + C_nbar + S (s2 > 0)'] < 1e-20 and r['worst rel deviation of the 3 -> 3 cancellation'] < 1e-5
        assert len(r['cases']) == 4 and all(abs(c['pole coefficients x eta/(2 pi i/A) at the smallest eta (collinear, soft)'][0] + 2) < 0.1 for c in r['cases'])
        m = new['explicit_two_loop_graph']
        assert m['worst rel deviation of Re f0 from -(pi/2) Im f1'] < 1e-3 and len(m['rows']) == 18
        assert all(abs(row['Glauber region from the leading-power integrand']/row['-(pi/2) Im f1'] - 1) < 1e-12 for row in m['rows'])
        q = new['qcd_imaginary_parts']
        assert q['largest |Im r0|'] < 1e-8 and q['largest |Im r0|/uncertainty'] < 1 and {'gg -> HH', 'gg -> ZZ'} <= set(q['processes'])
        h = new['n4_two_to_four']
        assert h['poles cancel with the kernel L_*/eps, L_* = -X_r + i pi'] and not h['depends on t1, t2 or v']
        assert h['c_G from the single logarithm'] == '0'
        de, st = h['DE coefficients of ln^k delta in E^(2) (H^v_w = H_w(1-v))'], h['structure (kappa_1 = zeta_2 + c_G)']
        assert all(sp_equal(de[k], st[k].replace('c_G', '0')) for k in ('3', '2', '1'))
        assert sp_equal(de['1'], '-2*I*pi**3/3') and sp_equal(de['2'], '-pi**2/2') and sp_equal(de['3'], 'I*pi/6')
        print(f'imaginary part of the kernel: 2 pi i Dc/eps and no finite part at mu = mu_b; the four eta-regulated regions sum to the '
              f'exact exchange to {r["worst rel diff of the sum from J at A/|s1 s2| = 5e5"]:.0e} with the Glauber region alone the -i pi term '
              f'(to {r["worst rel deviation of G from 2 pi^2/A"]:.0e}); the explicit two-loop graph has '
              f'Re f0 = Glauber region = -(pi/2) Im f1 to {m["worst rel deviation of Re f0 from -(pi/2) Im f1"]:.0e}; Im r0 = 0 in '
              f'{", ".join(q["processes"])} (largest {q["largest |Im r0|"]:.0e}, at most {q["largest |Im r0|/uncertainty"]:.2f} of '
              f'the fit uncertainty); planar N = 4 in 2 -> 4: c_G = {h["c_G from the single logarithm"]}, no v dependence')
        # the conditional three-loop prediction (three_loop.py -> three_loop.json, Paper II Sec. 7.2): rerun and compare;
        # within the boundary ansatz, the coefficients contain no open quantity and the two-loop checks must pass
        name = 'three_loop'
        stored = json.load(open(os.path.join(HERE, name + '.json')))
        shutil.copy(os.path.join(HERE, name + '.json'), os.path.join(tmp, name + '.json'))
        subprocess.run([sys.executable, 'three_loop.py'], cwd=HERE, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        new = json.load(open(os.path.join(HERE, name + '.json')))
        assert new == stored, 'three_loop.json differs from the stored file'
        assert all(not new[nu][k]['open'] for nu, ks in (('nu=0', ('G^5', 'G^4', 'G^3')), ('nu=2', ('G^4', 'G^3', 'G^2'))) for k in ks)
        assert new['nu=0']['G^5']['coefficient'] == '8*Delta_C**2/5' and new['nu=2']['G^4']['coefficient'] == '8*Delta_C**2'
        print('three_loop.json reproduced: the G^5, G^4, G^3 (nu = 0) and G^4, G^3, G^2 (nu = 2) coefficients of the three-loop '
              'amplitude follow from Delta_C, beta_0, B, the two-loop cusp, r_0 and h within the boundary ansatz '
              'with no independent quadratic or higher-power boundary logarithm; the next coefficients remain open')
        print('passed')
    finally:
        for name in ('r0_theory_quarks', 'r0_theory_gluons', 'kappa1_n4', 'glauber', 'three_loop'):
            if os.path.exists(os.path.join(tmp, name + '.json')):
                shutil.copy(os.path.join(tmp, name + '.json'), os.path.join(HERE, name + '.json'))
        shutil.rmtree(tmp)


if __name__ == '__main__':
    main()
