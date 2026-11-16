"""A failed or incomplete scientific comparison must not write a result."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

HERE=Path(__file__).resolve().parent


class FailureTests(unittest.TestCase):
    def test_nmhv_disagreement_and_missing_coefficient_fail(self):
        for missing in (False,True):
            with self.subTest(missing=missing),tempfile.TemporaryDirectory() as tmp:
                p=Path(tmp);data=json.loads((HERE/'de_data.json').read_text())
                coefficients=data['nmhv_three_to_three']['2']['coefficients']
                if missing:del coefficients['1']
                else:coefficients['1']='2'
                (p/'bad.json').write_text(json.dumps(data))
                (p/'output.json').write_text('preserved output')
                run=subprocess.run([sys.executable,str(HERE/'stage_d.py'),str(p/'bad.json'),str(p/'output.json')],capture_output=True,text=True)
                self.assertNotEqual(run.returncode,0,run.stdout)
                self.assertEqual((p/'output.json').read_text(),'preserved output')


if __name__=='__main__':unittest.main()
