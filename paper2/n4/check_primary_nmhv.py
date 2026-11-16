"""Brown SVHPL endpoint coefficients from shuffle-regularized associators.

No fitted kernel constants enter the endpoint construction. The primary
amplitude is read from the authors' hexMRKL1-7.m (arXiv:1903.10890).
The associator construction follows arXiv:1207.0186, section 3.3.
The first holomorphic endpoint term at weight nine needs regularized
associator coefficients only through weight eight.
"""
from collections import Counter, defaultdict
from functools import cache
from itertools import product
from pathlib import Path
import argparse, hashlib, json, re
import sympy as s
from sympy.parsing.mathematica import parse_mathematica

P=s.Symbol('p2'); z3,z5,z7,z53=s.symbols('z3 z5 z7 z53')
u,a=s.symbols('u a')
def clean(x):return s.expand(x)
def zeta(n):
    return clean(s.zeta(n).subs(s.pi**2,P)) if n%2==0 else {3:z3,5:z5,7:z7}[n]

@cache
def shuffles(x,y):
    if not x:return Counter({y:1})
    if not y:return Counter({x:1})
    r=Counter()
    for w,c in shuffles(x[1:],y).items():r[(x[0],)+w]+=c
    for w,c in shuffles(x,y[1:]).items():r[(y[0],)+w]+=c
    return r

@cache
def stuffles(x,y):
    if not x:return Counter({y:1})
    if not y:return Counter({x:1})
    r=Counter()
    for w,c in stuffles(x[1:],y).items():r[(x[0],)+w]+=c
    for w,c in stuffles(x,y[1:]).items():r[(y[0],)+w]+=c
    for w,c in stuffles(x[1:],y[1:]).items():r[(x[0]+y[0],)+w]+=c
    return r

def indices(w):
    out=[];n=1
    for b in w:
        if b==0:n+=1
        else:out.append(n);n=1
    assert n==1
    return tuple(out)

def word(ns):return sum(((0,)*(n-1)+(1,) for n in ns),())
def convergent(n):return [(0,)+w+(1,) for w in product((0,1),repeat=n-2)]
VALUES={():s.Integer(1)}

@cache
def regularize(w):
    """Linear combination of convergent words; H_0(1)=H_1(1)_reg=0."""
    if not w:return {():s.Integer(1)}
    if len(set(w))==1:return {}
    if w[-1]==0:terms=shuffles(w[:-1],(0,))
    elif w[0]==1:terms=shuffles((1,),w[1:])
    else:return {w:s.Integer(1)}
    mult=terms[w];out=defaultdict(lambda:s.Integer(0))
    for v,c in terms.items():
        if v!=w:
            for q,r in regularize(v).items():out[q]-=s.Rational(c,mult)*r
    return {q:clean(v) for q,v in out.items() if v!=0}

def H(w, vals=VALUES):return clean(sum(c*vals[q] for q,c in regularize(w).items()))

def mzv_table(n):
    words=convergent(n);vs=s.symbols('v0:'+str(len(words)))
    vals={**VALUES,**dict(zip(words,vs))}
    eq=[vals[(0,)*(n-1)+(1,)]-zeta(n)]
    if n==8:eq.append(vals[word((5,3))]-z53)
    for w in words:
        dual=tuple(1-b for b in w[::-1]);eq.append(vals[w]-vals[dual])
    for m in range(2,n-1):
        for x in convergent(m):
            for y in convergent(n-m):
                target=H(x)*H(y)
                eq.append(sum(c*H(w,vals) for w,c in shuffles(x,y).items())-target)
                eq.append(sum(c*H(word(ns),vals) for ns,c in stuffles(indices(x),indices(y)).items())-target)
    # Extended double shuffle with a single regularized zeta(1).
    if n>=3:
        for y in convergent(n-1):
            eq.append(sum(c*H(word(ns),vals) for ns,c in stuffles((1,),indices(y)).items()))
    eq=[clean(v) for v in eq if v!=0]
    sol=s.linsolve(eq,vs)
    if sol is s.EmptySet:raise ArithmeticError('inconsistent MZV equations '+str(n))
    row=next(iter(sol))
    if any(v.has(*vs) for v in row):raise ArithmeticError('undetermined MZV weight '+str(n))
    VALUES.update({w:clean(v) for w,v in zip(words,row)})
    print('MZV weight',n,'words',len(words),'relations',len(eq),flush=True)

def add(x,y,scale=1):
    r=dict(x)
    for w,v in y.items():r[w]=r.get(w,0)+scale*v
    return {w:clean(v) for w,v in r.items() if clean(v)!=0}
def mul(x,y,N):
    r=defaultdict(lambda:s.Integer(0))
    for wx,vx in x.items():
        for wy,vy in y.items():
            if len(wx)+len(wy)<=N:r[wx+wy]+=vx*vy
    return {w:clean(v) for w,v in r.items() if clean(v)!=0}
def inverse(x,N):
    assert x.get((),0)==1
    q=add(x,{():-1});out={():s.Integer(1)};p=out
    for k in range(1,N+1):
        p=mul(p,q,N)
        if not p:break
        out=add(out,p,(-1)**k)
    return out
def substitute(x,y1,N):
    y0={(0,):s.Integer(1)}
    @cache
    def sub(w):
        if not w:return {():s.Integer(1)}
        return mul(sub(w[:-1]),y0 if w[-1]==0 else y1,N)
    out={}
    for w,v in x.items():out=add(out,sub(w),v)
    return out

def associator_endpoint(N=8):
    Z={w:H(w) for n in range(N+1) for w in product((0,1),repeat=n)}
    Z={w:v for w,v in Z.items() if v!=0}
    R={w[::-1]:v for w,v in Z.items()}
    xi={(1,):s.Integer(1)};y=xi
    target=mul(mul(inverse(Z,N),xi,N),Z,N)
    for i in range(N):
        Ry=substitute(R,y,N)
        new=mul(mul(inverse(Ry,N),target,N),Ry,N)
        delta=add(new,y,-1)
        print('associator iteration',i,'changed',len(delta),flush=True)
        y=new
        if not delta:break
    else:raise ArithmeticError('associator did not stabilize')
    A=mul(Z,substitute(R,y,N),N)
    if add(mul(A,y,N),mul(xi,A,N),-1):raise ArithmeticError('monodromy equation failed')
    # First nontrivial y_1 terms, eq. (3.14) of arXiv:1207.0186.
    expected={tuple(map(int,w)):-z3*c for w,c in
        [('0011',2),('0101',-4),('0111',2),('1010',4),('1011',-6),('1100',-2),('1101',6),('1110',-2)]}
    assert {w:v for w,v in y.items() if len(w)==4}==expected
    return A,y

def endpoint_first(A,N=9):
    @cache
    def constant(w):
        ans=0
        for k in range(len(w)+1):
            if all(v==1 for v in w[:k]):ans+=(-u)**k/s.factorial(k)*A.get(w[k:],0)
        return clean(ans)
    @cache
    def first(w):
        if not w:return s.Integer(0)
        q= -constant(w[1:]) if w[0]==0 else -first(w[1:])
        out=0
        for k in range(s.degree(q,u)+1 if q else 0):out+=(-1)**k*s.diff(q,u,k)
        return clean(out)
    return constant,first

def main():
    here=Path(__file__).resolve().parent
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--primary',type=Path,default=here/'inputs/hexMRKL1-7.m')
    parser.add_argument('--data',type=Path,default=here/'de_data.json')
    parser.add_argument('--prediction',type=Path,default=here/'results_stage_d.json')
    parser.add_argument('--output',type=Path,default=here/'results_primary_nmhv.json')
    args=parser.parse_args()
    primary_bytes=args.primary.read_bytes()
    digest=hashlib.sha256(primary_bytes).hexdigest()
    if digest!='56d7ae8ef863b218613ba549dfb86155d03ceb3dc08f8fe8caeb1ce323c65fcc':
        raise ValueError('Unexpected primary ancillary data; expected the recorded 2019 source')
    for n in range(2,9):mzv_table(n)
    A,y=associator_endpoint()
    constant,first=endpoint_first(A)
    src=primary_bytes.decode()
    pred={}
    for L in range(1,6):
        m=re.search(r'FMtres\['+str(L)+r', 0\] = (.*?)(?=\nFMtres\[|\Z)',src,re.S)
        expr=parse_mathematica(m.group(1).strip())/s.pi
        atoms=[v for v in expr.atoms(s.Function) if v.func.__name__=='cL']
        def unpack(v):return word(tuple(map(int,v.args))) if tuple(v.args)!=(0,) else (0,)
        # cL[0] is the sole zero-tail compressed index in these expressions.
        repl={v:s.Symbol('cl_'+('_'.join(str(x) for x in v.args))) for v in atoms}
        pol=expr.xreplace(repl).replace(
            lambda q:q.is_Function and q.func.__name__ in ('Zeta','zeta'),
            lambda q:zeta(int(q.args[0])))
        zero={v:constant(unpack(cl)) for cl,v in repl.items()}
        At=sum(s.diff(pol,v).subs(zero)*first(unpack(cl)) for cl,v in repl.items())
        pred[L]=-clean(At.subs({s.pi**2:P,s.zeta(3):z3,s.zeta(5):z5,s.zeta(7):z7}))
        print('raw phase L',L,pred[L],flush=True)
    gam=4*a-4*zeta(2)*a*a+22*zeta(4)*a**3+(-s.Rational(219,2)*zeta(6)-4*z3**2)*a**4
    norm=s.series(s.exp(gam*zeta(2)/2),a,0,5).removeO()
    full=clean(sum(a**L*pred[L] for L in pred)*norm)
    data=json.loads(args.data.read_text())
    repo=json.loads(args.prediction.read_text())
    loc={f'z{k}':zeta(k) for k in range(2,9)}
    rows=[]
    for L in range(1,6):
        expected=data['nmhv_three_to_three'][str(L)]['coefficients'] if L<5 else repo['nmhv_five_loop_prediction']
        if set(expected)!={str(k) for k in range(2*L-1)}:
            raise ValueError('Incomplete target coefficient inventory at loop '+str(L))
        if s.degree(full.coeff(a,L),u)>2*L-2:
            raise ArithmeticError('Unexpected higher endpoint logarithm')
        for k,v in expected.items():
            if L==5:v=re.sub(r'zeta_(\d+)',r'z\1',v)
            target=s.sympify(v,locals=loc);got=full.coeff(a,L).coeff(u,int(k));diff=clean(got-target)
            rows.append({'loop':L,'power':int(k),'primary':str(got),'target':str(target),'difference':str(diff)})
            print('COMPARE',L,k,diff,flush=True)
    if any(r['difference']!='0' for r in rows):raise ArithmeticError('comparison failed')
    result={'scope':'Selected NMHV self-crossing component in the multi-Regge overlap',
            'primary_sha256':digest,'primary_reference':'arXiv:1903.10890v2, hexMRKL1-7.m',
            'method_reference':'arXiv:1207.0186, Brown associator construction',
            'sympy_version':s.__version__,'comparisons':rows,'passed':True,
            'lower_loop_coefficients':16,'five_loop_coefficients':9}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: 16 complete lower-loop and 9 complete five-loop coefficients',flush=True)

if __name__=='__main__':main()
