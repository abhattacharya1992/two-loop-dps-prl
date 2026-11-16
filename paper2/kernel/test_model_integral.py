"""Regression tests for the standalone longitudinal integral and its comparisons."""
import unittest
from unittest.mock import patch
import mpmath as mp
import model_exact2 as exact
import model_fast as fast


class LongitudinalIntegral(unittest.TestCase):
    def setUp(self):
        self.precision = mp.workdps(30)
        self.precision.__enter__()
        self.addCleanup(self.precision.__exit__, None, None, None)
        self.x = mp.mpf('.75')
        self.k = (mp.mpf('.018'), mp.mpf('.009'))
        self.q = (mp.mpf('.03'), mp.mpf(0))
        self.lt = (mp.mpf('.006'), mp.mpf('-.015'))
        self.eps0 = mp.mpf('9e-8')

    def test_both_negative_pinches_satisfy_the_pole_equation(self):
        roots = exact.pinch_points(self.lt, self.k, self.q, self.x)
        self.assertEqual(len(roots), 2)
        k2 = sum(v*v for v in self.k)
        kl2 = sum((a+b)**2 for a, b in zip(self.k, self.lt))
        a = sum((k+l-q)**2 for k, l, q in zip(self.k, self.lt, self.q))
        a6 = sum((k-q)**2 for k, q in zip(self.k, self.q))
        for lp in roots:
            self.assertTrue(-(1-self.x) < lp < -a6/self.x)
            z6 = 1-self.x+self.x*a/(self.x*lp+a6)
            z2 = kl2/self.x+k2/(1-self.x+lp)
            self.assertLess(z2, 1-self.x)
            self.assertLess(abs(z6-z2), mp.mpf('1e-25'))

    def test_partition_stays_inside_the_requested_interval(self):
        cap = mp.mpf('.2')
        points = exact.segments(self.lt, self.k, self.q, self.x, self.eps0, cap)
        self.assertEqual(points[0], -(1-self.x))
        self.assertEqual(points[-1], cap)
        self.assertTrue(all(a < b for a, b in zip(points, points[1:])))
        for lp in exact.pinch_points(self.lt, self.k, self.q, self.x):
            self.assertIn(lp, points)
        for eps in (0, -1, mp.inf, mp.nan):
            with self.assertRaises(ValueError):
                exact.segments(self.lt, self.k, self.q, self.x, eps, cap)

    def test_standalone_reference_and_production_integral(self):
        # No production partition or pinch finder is used by the reference.
        with patch.object(fast, 'segments', side_effect=RuntimeError('shared partition')):
            reference = exact.lplus_integral(self.lt, self.k, self.q, self.x, self.eps0)
        expected = mp.mpc('29315454560.6104453555181653514',
                          '31675394086.5394042936001713342')
        self.assertLess(abs(reference-expected)/abs(expected), mp.mpf('1e-20'))
        value = fast.lplus_integral(*map(float, (*self.lt, *self.k, *self.q,
                                               self.x, self.eps0)))
        self.assertLess(fast.check_reference(value, complex(reference)), 1e-8)

    def test_inaccurate_and_nonfinite_comparisons_fail(self):
        for value, reference in ((2, 1), (complex('nan'), 1), (1, complex('inf'))):
            with self.assertRaises(ValueError):
                fast.check_reference(value, reference)
        self.assertEqual(fast.check_reference(0, 0), 0)
        with self.assertRaises(ValueError):
            fast.check_reference(1, 1, rtol=-1)


if __name__ == '__main__':
    unittest.main()
