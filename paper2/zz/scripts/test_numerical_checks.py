"""Numerical disagreement, missing coverage and nonfinite outputs must fail."""
import unittest
import mpmath as mp
from check_cxx import compare_output
from numeric_checks import require_close
from zz_routeB import MasterFunctions
from ggvv import Functions, li22


class ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.py={j:{p:mp.mpf(1) for p in ('1l','2l_CF','2l_CA')} for j in range(1,10)}
        self.rows=[f'{j} {p} 1 0' for j,v in self.py.items() for p in v]

    def test_complete_matching_output(self):
        self.assertEqual(compare_output(self.py,'\n'.join(self.rows),rtol='1e-14'),
                         {'1l':0,'2l_CF':0,'2l_CA':0})

    def test_full_100_percent_error_fails(self):
        rows=[r.rsplit('1 0',1)[0]+'2 0' for r in self.rows]
        with self.assertRaises(ValueError):compare_output(self.py,'\n'.join(rows),rtol='1e-14')

    def test_missing_duplicate_and_nonfinite_results(self):
        cases=['','\n'.join(self.rows[:-1]),'\n'.join(self.rows+[self.rows[0]])]
        for value in ('nan','inf','-inf'):
            cases.append('\n'.join([f'1 1l {value} 0']+self.rows[1:]))
        for output in cases:
            with self.subTest(output=output[:30]),self.assertRaises(ValueError):
                compare_output(self.py,output,rtol='1e-14')

    def test_zero_and_invalid_tolerances(self):
        self.assertEqual(require_close(0,0,rtol=0),0)
        with self.assertRaises(ValueError):require_close(1,0,rtol='1e-14')
        for value in ('nan','inf','-1'):
            with self.assertRaises(ValueError):require_close(1,1,rtol=value)

    def test_master_cache_keeps_distinct_high_precision_arguments(self):
        with mp.workdps(40):
            a=mp.mpf(1)+mp.mpf('1e-36')
            b=mp.mpf(1)+mp.mpf('2e-36')
            functions=MasterFunctions()
            first=functions('Log',[a]);second=functions('Log',[b])
            self.assertNotEqual(first,second)
            self.assertEqual(second,mp.log(b))
        with mp.workdps(70):
            value=functions('PolyLog',[2,mp.mpf('0.3')])
            self.assertEqual(value,mp.polylog(2,mp.mpf('0.3')))

    def test_special_function_cache_respects_working_precision(self):
        functions = Functions()
        with mp.workdps(20):
            low_li = functions.lin(2, mp.mpf('0.5'))
            low_li22 = functions.li22(mp.mpf('0.5'), mp.mpf('0.25'))
        with mp.workdps(70):
            high_li = functions.lin(2, mp.mpf('0.5'))
            high_li22 = functions.li22(mp.mpf('0.5'), mp.mpf('0.25'))
            self.assertNotEqual(low_li, high_li)
            self.assertNotEqual(low_li22, high_li22)
            self.assertEqual(high_li, mp.polylog(2, mp.mpf('0.5')))
            self.assertEqual(high_li22, li22(mp.mpf('0.5'), mp.mpf('0.25')))


if __name__=='__main__':unittest.main()
