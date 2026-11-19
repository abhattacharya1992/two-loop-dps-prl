"""Check that multiple reported coefficients reuse the same interpolation fits."""
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from dps import analysis, layout, points, theory  # noqa: E402


class FitReuseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.hh, _ = points.load_hh(layout.path('hh_points'))
        cls.zz, _ = points.load_zz(layout.path('zz_points'))

    def test_gluon_fit_supplies_both_logarithmic_powers(self):
        original = analysis.fit.with_shift
        with patch.object(analysis.fit, 'with_shift', wraps=original) as fitted:
            rows = analysis.leading_logs_gluon(theory.HH, self.hh, '1/2')
        self.assertEqual(len(rows), 4 * len(theory.GLUON_COLOURS))
        self.assertEqual(fitted.call_count, 2 * len(theory.GLUON_COLOURS))

    def test_quark_one_loop_fits_are_reused_for_leading_logs(self):
        expected = sum(beta == '1/2' for _, beta in self.zz)
        original = analysis.one_loop_quark
        with patch.object(analysis, 'one_loop_quark', wraps=original) as fitted:
            self.assertTrue(analysis.leading_logs_quark(self.zz, '1/2'))
        self.assertEqual(fitted.call_count, expected)

    def test_quark_one_loop_fits_are_reused_for_matching(self):
        expected = sum(beta == '1/2' for _, beta in self.zz)
        original = analysis.one_loop_quark
        with patch.object(analysis, 'one_loop_quark', wraps=original) as fitted:
            self.assertTrue(analysis.r0_quark(self.zz, '1/2'))
        self.assertEqual(fitted.call_count, expected)


if __name__ == '__main__':
    unittest.main()
