"""Determine the crossed nonplanar labels from bare weight-four residues.

No finite C++ amplitude is used. The bare export defines INT=(Ceps/2)^2
times M at epsilon=0, hence the coefficient of a nonplanar M[...,4] in
the finite export equals one quarter of the bare epsilon^-4 residue.
The one-loop amplitude and its counterterms contain no nonplanar masters.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'paper2/zz/scripts'))
from formeval import letters
from ggvv import GGVV
from zz_master_families import FAMILY

REF=re.compile(r'M\(\s*(\w+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)')
BARE_REF=re.compile(r'INT\(\s*(\w+)\s*,\s*(\d+)\s*\)')
KINEMATICS=['15625','-2384.89727896833','6467.215561','625']

def libraries(executable):
    """Identify the native shared libraries recorded in the result provenance.

    Retained from DPS_papers/checks/regenerate_native_pairs.py at commit
    65efa9a70c73d56aa3f55cc612e3d412603d63f5, without its unrelated amplitude driver.
    """
    result=subprocess.run(['ldd',str(executable)],check=True,capture_output=True,text=True)
    paths=[Path(x) for x in re.findall(r'(/[^\s]+)',result.stdout)]
    paths=[p.resolve() for p in paths if p.is_file()]
    if not all(any(name in p.name for p in paths) for name in ('ginac','cln')):
        raise ValueError('Native GiNaC/CLN shared libraries were not identified')
    return {'shared/'+p.name:p for p in paths}

def verify_mapping(mapping,proof):
    expected={(family,index) for family in ('MN13m','MN13cm') for index in (40,42,43)}
    rows=proof['finite_coefficients']
    if len(rows)!=len(expected) or {(r['family'],r['index']) for r in rows}!=expected:
        raise ValueError('Incomplete nonplanar residue comparison')
    bare={family+','+str(index) for family in ('N13x34','N13x34x12') for index in (40,42,43)}
    if set(proof['bare_residues_over_four'])!=bare:
        raise ValueError('Incomplete bare residue set')
    with mp.workdps(80):
        for row in proof['finite_coefficients']:
            family=mapping[row['family']]
            target=proof['bare_residues_over_four'][family+','+str(row['index'])]
            a,b=mp.mpc(*row['value']),mp.mpc(*target)
            deviation=abs(a-b)/max(1,abs(a))
            if not mp.isfinite(deviation) or deviation>mp.mpf(proof['tolerance']):
                raise ValueError(f"{row['family']},{row['index']}: bare residue disagrees with {family}: {mp.nstr(deviation,8)}")

def extract(path,name):
    found=False;lines=[]
    with gzip.open(path,'rt') as source:
        for line in source:
            head=re.match(r'^g \[(A\d+_A_\dl(?:_qT)?)\] =',line)
            if head:
                if found:break
                found=head.group(1)==name
            if found:lines.append(line)
    if not found:raise ValueError('Missing coefficient '+name)
    text=''.join(lines).split('=',1)[1]
    return text[:text.rindex(';')]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--finite',type=Path,required=True,help='original gg_Aj_A_qT-1.0.inc.gz')
    parser.add_argument('--bare',type=Path,required=True,help='original gg_Aj_A-1.0.inc.gz')
    parser.add_argument('--prefix',type=Path,default=Path('/usr'))
    parser.add_argument('--work',type=Path,required=True,help='new build/output directory')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();work=args.work.resolve()
    if work.exists():raise ValueError('Use a new work directory')
    hashes={n:h for h,n in (line.split() for line in (ROOT/'paper2/zz/inputs/SHA256SUMS').read_text().splitlines())}
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    for path in (args.finite,args.bare):
        if hashes.get(path.name)!=sha(path):raise ValueError('Upstream input checksum mismatch: '+str(path))
    work.mkdir(parents=True)
    cpp=ROOT/'checks/native_form.cpp';exe=work/'native_form'
    prefix=args.prefix.resolve();lib=prefix/'lib/x86_64-linux-gnu'
    command=['g++','-O2','-std=c++17','-I'+str(prefix/'include'),str(cpp),
             '-L'+str(lib),'-Wl,-rpath,'+str(lib),'-lginac','-lcln','-o',str(exe)]
    subprocess.run(command,check=True)
    finite=extract(args.finite,'A1_A_2l_qT');bare=extract(args.bare,'A1_A_2l')
    finite_refs=sorted(set(REF.findall(finite)));bare_refs=sorted(set(BARE_REF.findall(bare)))
    fids={key:'M'+str(i) for i,key in enumerate(finite_refs)}
    bids={key:'J'+str(i) for i,key in enumerate(bare_refs)}
    fp=work/'finite.expression';bp=work/'bare.expression'
    fp.write_text(REF.sub(lambda m:fids[m.groups()],finite))
    bp.write_text(BARE_REF.sub(lambda m:bids[m.groups()],bare).replace('[','(').replace(']',')'))
    mp.mp.dps=160;s,t,a,b=map(mp.mpf,KINEMATICS)
    u,x,y,z=GGVV.variables(None,s,t,a,b)
    symbols=dict(letters(x,y,z),t=t,u=u,ma2=a,mb2=b,root=mp.sqrt((s-a-b)**2-4*a*b),Pi=mp.pi,i_=mp.j)
    pair=lambda v:[mp.nstr(mp.re(v),80),mp.nstr(mp.im(v),80)]
    kinetic=[[name,*pair(v)] for name,v in symbols.items()]
    counter=0
    def evaluate(expression,rows,digits):
        nonlocal counter
        path=work/f'symbols_{counter}.tsv';counter+=1
        path.write_text('\n'.join(' '.join(row) for row in rows)+'\n')
        result=subprocess.run([str(exe),str(expression),str(path),str(digits)],check=True,capture_output=True,text=True)
        terms={}
        for line in result.stdout.splitlines():
            ca,cf,r,i=line.split();key=int(ca),int(cf)
            if key in terms:raise ValueError('Duplicate native colour monomial')
            terms[key]=mp.mpc(r,i)
            if not mp.isfinite(terms[key]):raise ValueError('Nonfinite native coefficient')
        if not terms:raise ValueError('Empty native coefficient')
        return terms
    baseline=evaluate(fp,kinetic+[[fids[key],'0','0'] for key in finite_refs],80)
    rows=[];bare_values={}
    for index in (40,42,43):
        for family in ('MN13m','MN13cm'):
            active=(family,str(index),'4')
            unit=kinetic+[[fids[key],'1' if key==active else '0','0'] for key in finite_refs]
            value=evaluate(fp,unit,80).get((1,0),0)-baseline.get((1,0),0)
            rows.append({'family':family,'index':index,'value':pair(value)})
        for family in ('N13x34','N13x34x12'):
            values=[]
            for sign in (-1,1):
                eps=sign*mp.mpf('1e-25')
                unit=kinetic+[['d',mp.nstr(4-2*eps,160),'0']]
                unit += [[bids[key],'1' if key==(family,str(index)) else '0','0'] for key in bare_refs]
                terms=evaluate(bp,unit,160)
                values.append(eps**4*(terms.get((1,0),0)+terms.get((-1,0),0))/4)
            bare_values[family+','+str(index)]=pair(sum(values)/2)
            print('Derived bare residue',family,index,flush=True)
    proof={'decimal_inputs':KINEMATICS,'finite_coefficients':rows,'bare_residues_over_four':bare_values,
           'tolerance':'1e-45','bare_working_dps':160,'finite_working_dps':80,
           'symmetric_epsilon':'1e-25','method':'One quarter of the symmetric epsilon^-4 residue; the symmetric-limit error begins at epsilon^2',
           'scope':'Six nonplanar weight-four coefficients; no finite C++ amplitude values are used'}
    verify_mapping(FAMILY,proof)
    old=dict(FAMILY,MN13m='N13x34',MN13cm='N13x34x12')
    try:verify_mapping(old,proof)
    except ValueError:proof['old_mapping_rejected']=True
    else:raise ValueError('The old mapping unexpectedly passes the primary residue check')
    proof.update(passed=True,mapping=dict(FAMILY),
        hashes={str(p):sha(p) for p in (args.finite,args.bare,cpp,Path(__file__),ROOT/'paper2/zz/scripts/zz_master_families.py',exe)},
        library_hashes={name:sha(p) for name,p in libraries(exe).items()},
        compiler=subprocess.run(['g++','--version'],capture_output=True,text=True,check=True).stdout,
        compile_command=command)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(proof,indent=2)+'\n')
    print('PASS: six primary residue comparisons; old mapping rejected',flush=True)

if __name__=='__main__':main()
