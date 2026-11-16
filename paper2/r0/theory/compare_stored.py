"""Compare fresh region integrals with the preserved results.

Defaults allow 1e-7 absolute plus 1e-6 relative differences, consistent with the
precision of the production quadratures. Override explicitly for tighter runs:
    python compare_stored.py gqq --atol 1e-8 --rtol 1e-7
"""
import argparse
import json
import math
from pathlib import Path
from laurent_inputs import load_laurent


def validate_laurent(value, label):
    """Check payload shape; {} denotes a zero structure, not a missing key.

    load_laurent checks the complete file's required diagrams and structures
    before this helper is called.
    """
    if not isinstance(value,dict) or not value:
        raise ValueError(f'{label}: empty Laurent structures')
    for structure,coefficients in value.items():
        if not isinstance(coefficients,dict):
            raise ValueError(f'{label}/{structure}: missing Laurent coefficient mapping')
        for key,pair in coefficients.items():
            if (not isinstance(pair,list) or len(pair)!=2 or
                    not all(isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v) for v in pair)):
                raise ValueError(f'{label}/{structure}/{key}: malformed or nonfinite coefficient')


def compare(a, b, atol, rtol, label='result'):
    """Compare complete numerical trees; missing entries and nonfinite values fail."""
    if isinstance(a, dict):
        if not isinstance(b, dict) or a.keys() != b.keys():
            raise ValueError(f'{label}: different fields')
        return max((compare(a[k], b[k], atol, rtol, f'{label}/{k}') for k in a), default=0)
    if isinstance(a, list):
        if not isinstance(b, list) or len(a) != len(b):
            raise ValueError(f'{label}: different list lengths')
        return max((compare(x,y,atol,rtol,f'{label}/{i}') for i,(x,y) in enumerate(zip(a,b))), default=0)
    if isinstance(a, bool) or isinstance(b, bool):
        if a is not b:
            raise ValueError(f'{label}: different flags')
        return 0
    if isinstance(a, (int,float)) or isinstance(b, (int,float)):
        if not isinstance(a,(int,float)) or not isinstance(b,(int,float)):
            raise ValueError(f'{label}: incompatible numerical values')
        if not (math.isfinite(a) and math.isfinite(b)):
            raise ValueError(f'{label}: nonfinite value')
        delta=abs(a-b)
        if delta > atol + rtol*max(abs(a),abs(b)):
            raise ValueError(f'{label}: difference {delta:.3g} exceeds atol={atol:g}, rtol={rtol:g}')
        return delta
    if a != b:
        raise ValueError(f'{label}: different values')
    return 0


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kind',choices=('gqq','gg','n4'),nargs='?',default='gqq')
    parser.add_argument('--atol',type=float,default=1e-7)
    parser.add_argument('--rtol',type=float,default=1e-6)
    args=parser.parse_args()
    if not all(math.isfinite(t) and t>=0 for t in (args.atol,args.rtol)):
        parser.error('tolerances must be finite and nonnegative')
    failures=[];compared=0;worst=0
    for old in sorted(Path('stored').glob(f'{args.kind}_laurent_x*.json')):
        new=Path(old.name)
        try:
            x=old.stem.split('_x',1)[1].split('_',1)[0].replace('o','/')
            a=load_laurent(old,args.kind,expected_x=x)
            b=load_laurent(new,args.kind,expected_x=x)
            if a.get('complete') is False or b.get('complete') is False:
                raise ValueError(f'{new}: incomplete result')
            da,db=a.get('diagrams',{}),b.get('diagrams',{})
            if not da or da.keys()!=db.keys():
                raise ValueError(f'{new}: incomplete or incompatible diagram sets')
            compare(a.get('x'),b.get('x'),args.atol,args.rtol,f'{new}/x')
            delta=0
            for name in da:
                validate_laurent(da[name]['laurent'],f'{old}/{name}')
                validate_laurent(db[name]['laurent'],f'{new}/{name}')
                delta=max(delta,compare(da[name]['laurent'],db[name]['laurent'],args.atol,args.rtol,f'{new}/{name}'))
                colour='pref' if args.kind=='gqq' else 'colour'
                compare(da[name][colour],db[name][colour],args.atol,args.rtol,f'{new}/{name}/{colour}')
            print(f'{new}: largest change of a Laurent coefficient {delta:.1e}')
            compared+=1;worst=max(worst,delta)
        except (OSError,ValueError,KeyError,TypeError) as exc:
            failures.append(str(exc))
    result={'gqq':'r0_theory_quarks.json','gg':'r0_theory_gluons.json','n4':'kappa1_n4.json'}[args.kind]
    if (Path('stored')/result).exists():
        try:
            delta=compare(json.loads((Path('stored')/result).read_text()),json.loads(Path(result).read_text()),
                          args.atol,args.rtol,result)
            print(f'{result}: agrees within atol={args.atol:g}, rtol={args.rtol:g}; largest change {delta:.1e}')
        except (OSError,ValueError,TypeError) as exc:
            failures.append(str(exc))
    if not compared:failures.append('no compatible Laurent files were compared; comparison is unavailable')
    print(f'compared {compared} Laurent files; largest change {worst:.1e}')
    if failures:raise SystemExit('\n'.join(failures))


if __name__=='__main__':main()
