#!/usr/bin/env python3
"""Replay Stage-2 checks and require deliberate scientific errors to fail."""
from __future__ import annotations
import argparse, hashlib, json, shutil, subprocess, sys, tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():p.error('output must be new')
    source=(HERE/'verify_stage2.py').read_text(); rows=[]
    def run(script,output,opt=False):
        return subprocess.run([sys.executable,*(['-O'] if opt else []),str(script),'--output',str(output)],capture_output=True,text=True,timeout=120)
    with tempfile.TemporaryDirectory(prefix='dps_stage2_guard_') as td:
        root=Path(td); shutil.copy2(HERE/'inputs.json',root/'inputs.json')
        script=root/'verify_stage2.py';script.write_text(source)
        expected=json.loads((HERE/'results.json').read_text())
        expected.pop('runtime',None)
        for optimized in (False,True):
            out=root/f'replay_{optimized}.json';r=run(script,out,optimized)
            if r.returncode:raise RuntimeError('fresh replay failed: '+r.stderr)
            actual=json.loads(out.read_text());actual.pop('runtime',None)
            if actual!=expected:raise RuntimeError('replay differs from retained results')
            rows.append({'case':'optimized replay' if optimized else 'normal replay','passed':True})
        saved=(root/'replay_False.json').read_bytes()
        r=run(script,root/'replay_False.json')
        if r.returncode==0 or (root/'replay_False.json').read_bytes()!=saved:raise RuntimeError('overwrite protection failed')
        rows.append({'case':'existing output is rejected and preserved','passed':True})
        mutants=[
          ('wrong MHV-to-NMHV normalization','predicted=s.expand(2*sum','predicted=s.expand(1*sum'),
          ('wrong helicity principal symbol','ratio=(pp/2-2*h/pp)/(pm/2-2*h/pm)','ratio=(pm/2-2*h/pm)/(pp/2-2*h/pp)'),
          ('wrong spin-three Fourier sign','spin3_b=[c/6,c/6,c/6,-2*c/3]','spin3_b=[c/6,c/6,c/6,2*c/3]'),
          ('wrong pseudoscalar parity',"('A',FA,-1)","('A',FA,1)"),
          ('wrong known Sudakov normalization','Z[2]*gamma[j]/4+2*gamma3.get(j,0)','Z[2]*gamma[j]/4+gamma3.get(j,0)')]
        for i,(name,old,new) in enumerate(mutants):
            if source.count(old)!=1:raise RuntimeError('mutation anchor is not unique: '+name)
            script.write_text(source.replace(old,new));out=root/f'mutant_{i}.json';r=run(script,out)
            if r.returncode==0 or out.exists():raise RuntimeError('incorrect calculation accepted: '+name)
            if 'ArithmeticError' not in r.stderr:raise RuntimeError('mutant failed for a nonscientific reason: '+r.stderr)
            rows.append({'case':name,'rejected':True,'error':r.stderr.strip().splitlines()[-1]})
    record={'checker_sha256':hashlib.sha256(source.encode()).hexdigest(),
            'guard_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'checks':rows,'passed':True}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('x') as f:json.dump(record,f,indent=2);f.write('\n')
    print('PASS: two replays, overwrite protection, five scientific failure injections')
if __name__=='__main__':main()
