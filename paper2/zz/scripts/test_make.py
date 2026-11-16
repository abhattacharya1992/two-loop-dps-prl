#!/usr/bin/env python3
"""A failed point evaluation must fail the scan instead of being hidden by later points."""
from pathlib import Path
import subprocess
import sys
import tempfile


def test_scan_stops_on_failure():
    makefile = Path(__file__).resolve().parents[1] / 'Makefile'
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        evaluator = root / 'fail_first.py'
        evaluator.write_text(
            "from pathlib import Path\n"
            "import sys\n"
            "counter = Path('calls')\n"
            "n = int(counter.read_text()) if counter.exists() else 0\n"
            "counter.write_text(str(n + 1))\n"
            "sys.exit(1 if n == 0 else 0)\n"
        )
        result = subprocess.run(
            ['make', '--no-print-directory', '-f', str(makefile), 'scan',
             f'PY={sys.executable} {evaluator}'],
            cwd=root, capture_output=True, text=True,
        )
        assert result.returncode != 0, 'make scan hid a failed point evaluation'
        assert (root / 'calls').read_text() == '1', 'scan continued after an evaluator failed'


def test_checks_do_not_bypass_ca_disagreement():
    makefile=Path(__file__).resolve().parents[1]/'Makefile'
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp);(root/'logs').mkdir()
        evaluator=root/'ca_failure.py'
        evaluator.write_text(
            'import sys\n'
            'args=sys.argv[1:]\n'
            'if args and args[0].endswith("crosscheck_form.py") and "1l" not in args:\n'
            '    print("2l_CA fixture", flush=True)\n'
            '    sys.exit(0 if "--allow-unresolved-ca" in args else 1)\n'
            'print("fixture passed", flush=True)\n'
        )
        result=subprocess.run(['make','--no-print-directory','-f',str(makefile),'checks',
                               f'PY={sys.executable} {evaluator}'],cwd=root,capture_output=True,text=True)
        assert result.returncode!=0,'make checks bypassed a CA disagreement'


if __name__ == '__main__':
    test_scan_stops_on_failure()
    test_checks_do_not_bypass_ca_disagreement()
    print('OK: failed points stop the scan and CA disagreements fail make checks')
