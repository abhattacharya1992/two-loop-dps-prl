"""Colour proportionality must include components where the tree vanishes."""
import json
from pathlib import Path
import unittest
import numpy as np
import gg_diagrams as G
import n4_diagrams as N4

HERE=Path(__file__).resolve().parent

class ColourValidation(unittest.TestCase):
    def ratio_tensor(self,value):
        out=np.zeros_like(G._F,dtype=complex)
        out[np.abs(G._F)>1e-9]=value
        return out

    def test_zero_tree_component_is_checked(self):
        value=self.ratio_tensor(.5)
        self.assertEqual(G._ratio(value),.5)
        value[0,0,0]=1
        with self.assertRaises(ValueError):G._ratio(value)

    def test_nonfinite_and_malformed_tensors_fail(self):
        for number in (float('nan'),float('inf')):
            for index in ((0,0,0),(0,1,2)):
                value=self.ratio_tensor(1);value[index]=number
                with self.assertRaises(ValueError):G._ratio(value)
        with self.assertRaises(ValueError):G._ratio(np.zeros((8,8)))

    def test_all_production_factors_match_preserved_records(self):
        for kind,module in (('gg',G),('n4',N4)):
            records=json.loads((HERE/f'{kind}_laurent_x3o10.json').read_text())['diagrams']
            self.assertEqual(set(records),set(module.DIAGRAMS))
            for name,diagram in module.DIAGRAMS.items():
                with self.subTest(channel=kind,diagram=name):
                    actual=module.colour(diagram['colour']);expected=records[name]['colour']
                    self.assertEqual(set(actual),set(expected))
                    for key,value in actual.items():
                        self.assertLessEqual(abs(value-complex(*expected[key])),1e-12)

if __name__=='__main__':unittest.main()
