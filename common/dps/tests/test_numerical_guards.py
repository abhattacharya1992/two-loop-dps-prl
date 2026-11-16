"""Reject bad point records and incomplete or nonfinite numerical comparisons."""
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import mpmath as mp
from dps import points, fit, layout
from dps.colour import N


class NumericalGuards(unittest.TestCase):
    def test_zero_reference_and_component_sets(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            def write(name, values, **extra):
                (root/name).write_text(json.dumps({'dps':80,'A':values,**extra}))
            write('a', {'N2':['1','0']});write('b',{'N2':['0','0']})
            self.assertEqual(points.compare_stored(root,'a','b'), mp.inf)
            write('a',{'N2':['0','0']})
            self.assertEqual(points.compare_stored(root,'a','b'), 0)
            write('a',{'N2':['0','0'],'Nf2':['0','0']})
            with self.assertRaises(ValueError):points.compare_stored(root,'a','b')
            write('a',{'N2':['0','0']}, A1={'N':['1','0']})
            write('b',{'N2':['0','0']}, A1={'N':['0','0']})
            self.assertEqual(points.compare_stored(root,'a','b'), mp.inf)

    def test_actual_aa_regression_rejects_nan(self):
        spec=importlib.util.spec_from_file_location('scientific_tests',Path(__file__).with_name('test_dps.py'))
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for p in layout.AA_POINTS.glob('*.json'):shutil.copy2(p,root/p.name)
            target=root/'aa_g1_b4o5_d11_p149.json'
            record=json.loads(target.read_text());record['A']['N2'][0]='nan';target.write_text(json.dumps(record))
            with patch.object(layout, 'AA_POINTS',root), self.assertRaises(ValueError):module.test_aa()

    def test_nonfinite_data_and_fit_coefficients(self):
        for value in ('nan','inf','-inf'):
            with self.assertRaises(ValueError):points.cplx([value,'0'])
            with self.assertRaises(ValueError):points.positive_pT2(value)
        records=[{'G':mp.mpf(i),'pT2':mp.mpf(i+1)} for i in range(3)]
        with self.assertRaises(ValueError):fit.leading(records,[1,mp.nan,3],0,1)
        with self.assertRaises(ValueError):fit.plain([0,1],[1,mp.inf],1)
        with patch.object(mp, 'lu_solve', return_value=mp.matrix([mp.nan,1,1])), self.assertRaises(ValueError):
            fit.leading(records,[1,2,3],0,1)

    def test_colour_power_must_be_integral(self):
        with self.assertRaises(ValueError):N**mp.mpf('1.5')
        with self.assertRaises(ValueError):N**mp.mpf('-1.5')
        self.assertEqual((N**2).at(4,0),16)
        self.assertEqual((N**-1).at(4,0),mp.mpf(1)/4)


if __name__ == '__main__':unittest.main()
