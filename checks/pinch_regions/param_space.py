#!/usr/bin/env python3
"""Deterministic local scalar-pinch checks, not a region-completeness theorem.
Input vectors and momentum labels are supplied by the pinned momentum-space
analysis. Polynomial construction, exact sign witnesses and local saddle powers
are computed here. No four-dimensional divergent scalar integral is integrated.
"""
from __future__ import annotations
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import lzma, hashlib, json
from graphs import Graph, all_graphs, forest_terms

HERE=Path(__file__).resolve().parent
X,Y=Q(3,7),Q(2,5)
VECTORS_SHA256='6c8df4034f253442f2fbd2cdbc39cc4a8464a8dcd991b659f8f9c93381b9c33f'

def require(condition,message):
    if not condition:raise ArithmeticError(message)

def prod(values):
    out=Q(1)
    for v in values:out*=v
    return out

def ev(poly,point):return sum(c*prod(point[i]**k for i,k in enumerate(m) if k) for m,c in poly.items())
def support(poly):return sorted({i for m in poly for i,k in enumerate(m) if k})

def witnesses(poly):
    """Exact positive-orthant sign witnesses for homogeneous multiaffine polynomials."""
    poly={m:Q(c) for m,c in poly.items() if c}
    require(bool(poly),'zero polynomial has no strict-sign witnesses')
    require(len({sum(m) for m in poly})==1 and all(k in (0,1) for m in poly for k in m),'requires homogeneous multiaffine polynomial')
    result={}
    for sign in (-1,1):
        targets=[(m,c) for m,c in sorted(poly.items()) if sign*c>0]
        if not targets:continue
        m,c=targets[0];M=2+sum(abs(t) for mm,t in poly.items() if mm!=m)/abs(c)
        p=[M if k else Q(1) for k in m]
        require(sign*ev(poly,p)>0,'invalid sign witness');result[sign]=p
    return result

def factors(poly):
    """Exact disjoint-support rank-one factorisation; no momentum labels."""
    vars_=support(poly)
    if len(vars_)<=1:return [poly]
    for n in range(len(vars_)-1):
        for tail in combinations(vars_[1:],n):
            left=set((vars_[0],)+tail);right=set(vars_)-left
            table={(tuple(k if i in left else 0 for i,k in enumerate(m)),tuple(k if i in right else 0 for i,k in enumerate(m))):c for m,c in poly.items()}
            rows=sorted({a for a,b in table});cols=sorted({b for a,b in table})
            (a0,b0),c0=next(iter(sorted(table.items())))
            if all(table.get((a,b),0)*c0==table.get((a,b0),0)*table.get((a0,b),0) for a in rows for b in cols):
                A={a:table[(a,b0)] for a in rows if (a,b0) in table}
                B={b:table[(a0,b)]/c0 for b in cols if (a0,b) in table}
                return factors(A)+factors(B)
    return [poly]

@lru_cache(None)
def symanzik(name):
    g=Graph(name);U=defaultdict(Q);F=defaultdict(lambda:[Q(0),Q(0)])
    for kind,m,p in forest_terms(g):
        if kind=='U':U[m]+=1
        else:
            a,b,c=p;F[m][0]-=4*(a*b+a*c*Y+b*c*X+c*c*X*Y);F[m][1]+=4*c*c
    return dict(U),{m:tuple(v) for m,v in F.items() if any(v)},len(g.edges)-len(g.vertices)+1

def leading(U,F,w):
    ub=defaultdict(dict);fb=defaultdict(dict)
    for m,c in U.items():ub[-sum(a*b for a,b in zip(w,m))][m]=Q(c)
    for m,cs in F.items():
        p=-sum(a*b for a,b in zip(w,m))
        for off,c in zip((0,2),cs):
            if c:fb[p+off][m]=fb[p+off].get(m,0)+c
    fb={p:{m:c for m,c in terms.items() if c} for p,terms in fb.items()};fb={p:v for p,v in fb.items() if v}
    return min(ub),dict(sorted(fb.items()))

def root_point(f):
    ws=witnesses(f);require(len(ws)==2,'factor does not change sign')
    p=list(ws[1]);target=ws[-1];prev=ev(f,p)
    for i in support(f):
        nxt=p.copy();nxt[i]=target[i];now=ev(f,nxt)
        if now<=0:
            slope=(now-prev)/(nxt[i]-p[i]);require(slope!=0,'zero slope')
            p[i]-=prev/slope
            require(p[i]>0 and ev(f,p)==0,'invalid positive root')
            return p,i
        p=nxt;prev=now
    raise ArithmeticError('root not found')

def critical_order(A,B,C,D):
    """Order of (A D - B C)/A, with exact rational coefficients."""
    def mul(p,q):
        out=defaultdict(Q)
        for i,c in p.items():
            for j,d in q.items():out[i+j]+=c*d
        return out
    det=mul(A,D)
    for p,c in mul(B,C).items():det[p]-=c
    nonzero=[p for p,c in det.items() if c]
    aden=[p for p,c in A.items() if c]
    require(bool(nonzero) and bool(aden),'degenerate chart; choose another positive tangent point')
    return min(nonzero)-min(aden)

def local_power(U,F,L,w,van):
    muU,fb=leading(U,F,w);muF=min(fb);E=len(w);b=E-2*L;a=E-2*(L+1)
    require(len(van)==2,'local chart requires two factors')
    point=[Q(1)]*E;norm=[]
    for f in van:
        p,i=root_point(f);norm.append(i)
        for j in support(f):point[j]=p[j]
    i,j=norm
    D={p-muF:ev(f,point) for p,f in fb.items()}
    pi=point.copy();pi[i]+=1;pj=point.copy();pj[j]+=1;pij=pi.copy();pij[j]+=1
    B={p-muF:ev(f,pi)-D[p-muF] for p,f in fb.items()}
    C={p-muF:ev(f,pj)-D[p-muF] for p,f in fb.items()}
    A={p-muF:ev(f,pij)-D[p-muF]-B[p-muF]-C[p-muF] for p,f in fb.items()}
    delta=critical_order(A,B,C,D)
    facet=-sum(w)+a*muU-b*muF
    return facet+delta*(1-b),delta,b

def jet_split(factor_labels):
    """Check supplied momentum labels, not an independent DPS classifier."""
    if len(factor_labels)!=2:return False
    return {frozenset(x for x in f if x!='S') for f in factor_labels}=={frozenset({'1'}),frozenset({'2'})}

def load_rows():
    raw=(HERE/'vectors.json.xz').read_bytes()
    require(hashlib.sha256(raw).hexdigest()==VECTORS_SHA256,'supplied vector archive hash mismatch')
    pack=json.loads(lzma.decompress(raw))
    rows=[]
    for line in pack['records']:
        gid,ws,kind,km,lm,pden,omega,minimum,labels=line.split(';')
        rows.append(dict(graph=pack['graphs'][int(gid)],w=[Q(int(n),2) for n in ws.split(',')],kind=kind,K=pack['K'][int(km)],L=pack['L'][int(lm)],pden=Q(pden),old_omega=Q(omega) if omega!='-' else None,minimum=Q(minimum) if minimum!='-' else None,labels=list(labels)))
    return pack,rows

def run():
    pack,rows=load_rows();hidden=[];selected=[];all_scalar=[];checks=0;factor_count=0
    for r in rows:
        U,F,L=symanzik(r['graph']);_,fb=leading(U,F,r['w']);van=[]
        for f in factors(next(iter(fb.values()))):
            ws=witnesses(f);factor_count+=1
            if len(ws)==2:van.append(f)
        if len(van)>1:
            require(len(van)==2,'more than two vanishing factors')
            g=Graph(r['graph']);external=[];labels=[]
            for f in van:
                parents={g.edges[i][3] for i in support(f)}
                jets={1 for p in parents if p in ('k','kb','leg1')}|{2 for p in parents if p in ('r','rb','leg2')}
                require(len(jets)==1,'factor mixes incoming families');external.append(next(iter(jets)))
                labels.append({r['labels'][i] for i in support(f)})
            require(set(external)=={1,2},'missing incoming family')
            power,delta,b=local_power(U,F,L,r['w'],van)
            require(power==r['old_omega'],'local saddle power changed: '+r['graph'])
            r=dict(r,power=power,b=b,jet_split=jet_split(labels));hidden.append(r)
        if r['K']=='dps' and r['kind']=='D' and r['L'] in ('hard','c1+','c1-','c2+','c2-','usoft'):
            require(len(van)==2,'selected region has no pair product')
            require(power==r['pden'],'selected denominator power mismatch');all_scalar.append(r)
            if r['pden']<=0:selected.append(r)
        checks+=1
    key=lambda r:(r['graph'],tuple(v-min(r['w']) for v in r['w']))
    result={'source_commit':pack['source_commit'],'input_rows':checks,'distinct_vectors':len({key(r) for r in rows}),'hidden_rows':len(hidden),'distinct_hidden_vectors':len({key(r) for r in hidden}),'two_factor_rows':len(hidden),'hidden_DPS_labels':sum(r['kind']=='D' for r in hidden),'hidden_SPS_labels':sum(r['kind']=='S' for r in hidden),'SPS_rows_passing_jet_split':sum(r['kind']=='S' and r['jet_split'] for r in hidden),'selected_scalar_power_checks':len(selected),'scalar_power_checks':len(all_scalar),'b1_logarithmic_rows':sum(r['b']==1 for r in hidden),'local_power_less_than_grid':sum(r['power']<r['minimum'] for r in hidden),'local_power_greater_than_grid':sum(r['power']>r['minimum'] for r in hidden),'factors_checked':factor_count,'scope':'Supplied momentum-space vectors; exact polynomial and local scalar checks. No numerator powers, global dissection, or regulated overlap cancellation.'}
    expected={'input_rows':2823,'distinct_vectors':2363,'hidden_rows':1709,'distinct_hidden_vectors':1305,'hidden_DPS_labels':1031,'hidden_SPS_labels':678,'SPS_rows_passing_jet_split':0,'selected_scalar_power_checks':127,'scalar_power_checks':142,'b1_logarithmic_rows':140,'local_power_less_than_grid':422,'local_power_greater_than_grid':85}
    for k,v in expected.items():require(result[k]==v,f'{k}: got {result[k]}, expected {v}')
    return result

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);args=ap.parse_args()
    text=json.dumps(run(),indent=2,sort_keys=True)+'\n'
    if args.output:args.output.write_text(text)
    print(text,end='')
