"""The ggvvamp transcription check must reject empty or incomplete comparisons."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


class TranscriptionCheck(unittest.TestCase):
    def check(self, cpp, generated, *, optimise=False):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'p0002.cpp').write_text(cpp)
            (root/'p0002.py').write_text(generated)
            command = [sys.executable] + (['-O'] if optimise else [])
            command += [str(Path(__file__).with_name('check_p0002.py')), str(root), str(root)]
            return subprocess.run(command, capture_output=True, text=True)

    def statements(self):
        cpp = ''.join(f'r[{i}] = n(1);\n' for i in range(8, 1621))
        generated = ''.join(f'    r[{i}] = (1)\n' for i in range(8, 1621))
        return cpp, generated

    def test_complete_matching_transcription_passes(self):
        result = self.check(*self.statements())
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('1613 statements', result.stdout)

    def test_empty_and_truncated_comparisons_fail_even_with_optimisation(self):
        cpp, generated = self.statements()
        for optimise in (False, True):
            for a, b in (('', ''), (cpp, ''), ('', generated),
                         ('r[8] = n(1);\n', '    r[8] = (1)\n')):
                with self.subTest(optimise=optimise, length=len(a)):
                    self.assertNotEqual(self.check(a, b, optimise=optimise).returncode, 0)

    def test_duplicate_malformed_and_different_assignments_fail(self):
        cpp, generated = self.statements()
        cases = [(cpp+'r[8] = n(1);\n', generated),
                 (cpp, generated+'    r[8] = (1)\n'),
                 (cpp+'r[1621] += n(1);\n', generated),
                 (cpp, generated+'    r[1621] += (1)\n'),
                 (cpp, generated.replace('r[8] = (1)', 'r[8] = (2)'))]
        for a, b in cases:
            with self.subTest(lengths=(len(a), len(b))):
                self.assertNotEqual(self.check(a, b).returncode, 0)


if __name__ == '__main__':
    unittest.main()
