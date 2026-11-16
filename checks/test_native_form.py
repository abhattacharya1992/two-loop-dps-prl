"""Optional integrity checks of a compiled native_form reader."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import mpmath as mp

EXE=None
if __name__=='__main__':
    if len(sys.argv)!=2:raise SystemExit('usage: test_native_form.py /path/to/native_form')
    EXE=Path(sys.argv.pop()).resolve()

@unittest.skipIf(EXE is None,'Pass the compiled native_form executable to run these optional checks')
class NativeForm(unittest.TestCase):
    def evaluate(self,text):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);expr=root/'expression';syms=root/'symbols'
            expr.write_text(text)
            syms.write_text('Pi '+mp.nstr(mp.pi,90)+' 0\ni_ 0 1\n')
            return subprocess.run([str(EXE),str(expr),str(syms),'90'],cwd=root,capture_output=True,text=True)

    def test_precedence_exact_division_functions_and_colour(self):
        with mp.workdps(90):
            cases=[('-2^2+5',{(0,0):1}),('1/3',{(0,0):mp.mpf(1)/3}),
                   ('Num(2*(3+4))*Den(7)',{(0,0):2}),
                   ('Log(1)+i_*Pi',{(0,0):mp.j*mp.pi}),
                   ('CA*CF+2*N^-1',{(1,1):1,(-1,0):2})]
            for text,expected in cases:
                with self.subTest(expression=text):
                    r=self.evaluate(text);self.assertEqual(r.returncode,0,r.stderr)
                    actual={}
                    for line in r.stdout.splitlines():
                        a,b,x,y=line.split();key=int(a),int(b)
                        self.assertNotIn(key,actual);actual[key]=mp.mpc(x,y)
                    self.assertEqual(set(actual),set(expected))
                    for key,value in expected.items():
                        self.assertLess(abs(actual[key]-value),mp.mpf('1e-85'))

    def test_unsupported_malformed_and_singular_input_fails(self):
        for text in ['Den(0)','Num(1,2)','Log(1,2)','PolyLog(2,1,0)',
                     'PolyLog(1/2,1/3)','Missing(1)','unknown','1.2','1+','(1','CA+@']:
            with self.subTest(expression=text):
                self.assertNotEqual(self.evaluate(text).returncode,0)

if __name__=='__main__':unittest.main()
