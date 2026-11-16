"""Regression against coefficients derived independently from the bare export."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'checks'))
from check_zz_ca_mapping import libraries, verify_mapping
from zz_master_families import FAMILY

class NonplanarCrossings(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.proof=json.loads((Path(__file__).parent/'fixtures/nonplanar_weight4_residues.json').read_text())

    def test_finite_coefficients_match_bare_weight_four_residues(self):
        verify_mapping(FAMILY,self.proof)

    def test_original_interchange_fails_the_primary_residue_comparison(self):
        wrong=dict(FAMILY,MN13m='N13x34',MN13cm='N13x34x12')
        with self.assertRaisesRegex(ValueError,'bare residue disagrees'):
            verify_mapping(wrong,self.proof)

    def test_missing_or_nonfinite_primary_coefficients_fail(self):
        missing=copy.deepcopy(self.proof);missing['finite_coefficients'].pop()
        with self.assertRaisesRegex(ValueError,'Incomplete'):
            verify_mapping(FAMILY,missing)
        invalid=copy.deepcopy(self.proof);invalid['finite_coefficients'][0]['value'][0]='nan'
        with self.assertRaises(ValueError):verify_mapping(FAMILY,invalid)

    def test_native_library_provenance_requires_both_dependencies(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary)
            ginac=root/'libginac.so.11';cln=root/'libcln.so.6'
            ginac.write_bytes(b'GiNaC fixture');cln.write_bytes(b'CLN fixture')
            output=f'libginac.so.11 => {ginac} (0x1)\nlibcln.so.6 => {cln} (0x2)\n'
            with patch('check_zz_ca_mapping.subprocess.run') as run:
                run.return_value.stdout=output
                self.assertEqual(libraries(root/'reader'),
                                 {'shared/libginac.so.11':ginac,'shared/libcln.so.6':cln})
                cln.unlink()
                with self.assertRaisesRegex(ValueError,'shared libraries were not identified'):
                    libraries(root/'reader')

if __name__=='__main__':unittest.main()
