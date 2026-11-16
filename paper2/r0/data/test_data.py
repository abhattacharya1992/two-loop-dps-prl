"""Rerun the extraction of r0 from the two-loop amplitudes (the stored scan points and fits) and the blind tests,
and check that the stored results are reproduced and that the closed forms hold as the paper states.

    python3 test_data.py        (a few seconds)
"""
import json, os, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
FILES = ['r0_from_data', 'r0_from_zz', 'r0_closed_forms', 'blind_test_hh', 'blind_test_zz']


def walk(a, b, path=''):
    if isinstance(a, dict):
        assert set(a) == set(b), (path, set(a) ^ set(b))
        return max([walk(a[k], b[k], path + '/' + k) for k in a] or [0])
    if isinstance(a, list):
        assert len(a) == len(b), path
        return max([walk(x, y, path) for x, y in zip(a, b)] or [0])
    if isinstance(a, (int, float)):
        return abs(a - b)
    if isinstance(a, str):
        try:
            return abs(float(a) - float(b))
        except ValueError:
            assert a == b, (path, a, b)
            return 0
    assert a == b, (path, a, b)
    return 0


def main():
    tmp = tempfile.mkdtemp()
    try:
        for name in FILES:
            shutil.copy(os.path.join(HERE, name + '.json'), os.path.join(tmp, name + '.json'))
        for name in FILES:
            stored = json.load(open(os.path.join(tmp, name + '.json')))
            subprocess.run([sys.executable, name + '.py'], cwd=HERE, check=True, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
            new = json.load(open(os.path.join(HERE, name + '.json')))
            d = walk(stored, new)
            print(f'{name}.json reproduced, largest difference {d:.1e}')
            assert d < 1e-12
        cf = json.load(open(os.path.join(HERE, 'r0_closed_forms.json')))
        assert cf['HH'] < 1e-8 and cf['ZZ'] < 1e-8
        hh = json.load(open(os.path.join(HERE, 'blind_test_hh.json')))
        zz = json.load(open(os.path.join(HERE, 'blind_test_zz.json')))
        assert max(max(r['dev_C_A'], r['dev_n_f']) for r in hh) < 1e-8
        assert max(max(r['2l_CF']['dev'], r['2l_CA']['dev']) for r in zz) < 1e-8
        print(f"closed forms: HH {cf['HH']:.1e}, ZZ {cf['ZZ']:.1e}; blind tests below 1e-8")
        print('passed')
    finally:
        for name in FILES:
            if os.path.exists(os.path.join(tmp, name + '.json')):
                shutil.copy(os.path.join(tmp, name + '.json'), os.path.join(HERE, name + '.json'))
        shutil.rmtree(tmp)


if __name__ == '__main__':
    main()
