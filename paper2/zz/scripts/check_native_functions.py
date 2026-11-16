"""Compare original GiNaC p0002 with its Python transcription at fresh points.

Usage: check_native_functions.py NATIVE_EXECUTABLE BUILD_GGV OUTPUT POINT...
POINT is bench or BETA:K. Build the executable from cxx/native_functions.cpp,
the original ggvvamp.cpp and p0002.cpp; see README.md for modern GiNaC APIs.
This script never feeds Python special-function values to the native program.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import mpmath as mp

sys.path.insert(0,str(Path(__file__).resolve().parent))
from ggvv import Functions,RSIZE
from zz_points import kinematics
from numeric_checks import require_close

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable',type=Path)
    parser.add_argument('build',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('points',nargs='+')
    parser.add_argument('--digits',type=int,default=75)
    parser.add_argument('--rtol',default='1e-55')
    args=parser.parse_args()
    if args.digits<60:raise ValueError('At least 60 digits required')
    mp.mp.dps=args.digits
    ns={};exec(compile((args.build/'p0002.py').read_text(),'p0002','exec'),ns)
    ratios=json.loads((args.build/'rationals.json').read_text())
    Q=[mp.mpf(int(k.split('/')[0]))/int(k.split('/')[1])
       for _,k in sorted((v,k) for k,v in ratios.items())]
    records=[]
    for spec in args.points:
        if spec=='bench':
            values=list(map(mp.mpf,('15625','-2384.89727896833','6467.215561','625')))
        else:
            beta,k=spec.split(':');p=kinematics(beta,int(k))
            values=[p['s'],p['t'],mp.mpf(1),mp.mpf(1)]
        # Serialize and reparse the same decimal inputs before either evaluation.
        decimal=[mp.nstr(v,args.digits) for v in values]
        result=subprocess.run([str(args.executable.resolve()),str(args.digits),*decimal],
                              capture_output=True,text=True,check=True)
        native={}
        for line in result.stdout.splitlines():
            if line.startswith('R '):
                _,key,value=line.split();key=int(key)
                if key in native:raise ValueError('Duplicate native intermediate')
                native[key]=mp.mpf(value)
        if set(native)!=set(range(1,1621)):
            raise ValueError('Incomplete native special-function stage')
        s,t,ma,mb=map(mp.mpf,decimal);u=ma+mb-s-t
        root=mp.sqrt((s-ma-mb)**2-4*ma*mb)
        t0=(ma+mb-s-root)/2;t1=(ma+mb-s+root)/2
        r=[None]*(RSIZE+1)
        r[1:8]=[ma,t/ma,u/ma,mb/ma,-t0/ma,t1*(u+t)/(ma*mb)-1,t1*t/(ma*mb)]
        ns['p0002'](r,Q,Functions())
        worst=mp.mpf(0);worst_index=None
        for key in range(1,1621):
            # Mixed absolute/relative normalization is explicit at zero.
            require_close(r[key],native[key],rtol=args.rtol,atol=args.rtol,
                          label='native p0002 r['+str(key)+']')
            dev=abs(r[key]-native[key])/max(1,abs(native[key]))
            if dev>worst:worst,worst_index=dev,key
        records.append({'point':spec,'decimal_inputs':decimal,'intermediates':1620,
                        'worst_index':worst_index,'normalized_difference':mp.nstr(worst,25)})
        print(spec,1620,'intermediates; worst',mp.nstr(worst,6),flush=True)
    output={'scope':'Native GiNaC special-function stage, independent of Python values',
            'digits':args.digits,'rtol':args.rtol,'atol':args.rtol,
            'executable_sha256':hashlib.sha256(args.executable.read_bytes()).hexdigest(),
            'transcription_sha256':hashlib.sha256((args.build/'p0002.py').read_bytes()).hexdigest(),
            'points':records,'passed':True}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2)+'\n')

if __name__=='__main__':main()
