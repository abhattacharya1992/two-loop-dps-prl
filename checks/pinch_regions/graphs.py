"""DPS-box graph routings, adapted from the pinned extension_v2 graph generator.
Source: DPS_papers@2542267940bd82143bf71cdc01867619234b10af,
checks/pinch_completeness_20261113/graphs.py (BSD-3-Clause).
Only the canonical routing is needed: Symanzik polynomials do not depend on
which two independent internal momenta were named K,L in the input census.
"""
from collections import defaultdict
from itertools import combinations

LINES = {'k': ('a1','H1'), 'kb': ('a1','H2'), 'r': ('a2','H1'), 'rb': ('a2','H2')}
SITES = ['k','kb','r','rb','leg1','leg2','a1','a2','H1','H2']
TIPS = [('T1','H1'),('T2','H2'),('T1','T2'),('T1','H2'),('H1','T2')]

class Graph:
    def __init__(self, name):
        self.name=name
        A,B=name.split('@')[0].split('|')
        if A==B and A not in LINES: raise ValueError(name)
        cuts=defaultdict(list); legcut={}; ends=[]
        tips={s for s in (A,B) if s in ('T1','T2')}
        corner={'H1':'C1' if 'T1' in tips else 'H1','H2':'C2' if 'T2' in tips else 'H2'}
        for idx,site in enumerate((A,B)):
            if site in ('T1','T2'): ends.append('H'+site[1])
            elif site in ('H1','H2'): ends.append(corner[site])
            elif site in LINES:
                v=f'u{idx}';cuts[site].append(v);ends.append(v)
            elif site in ('leg1','leg2'):
                v=f'u{idx}';legcut[site]=v;ends.append(v)
            else: ends.append(site)
        self.edges=[]
        for name,(t,h) in LINES.items():
            chain=[t]+cuts[name]+[corner[h]]
            for i in range(len(chain)-1):
                self.edges.append((name if len(chain)==2 else f'{name}_{i+1}',chain[i],chain[i+1],name))
        self.ext=defaultdict(list)
        for leg,vertex,sym in (('leg1','a1','p1'),('leg2','a2','p2')):
            if leg in legcut:
                self.edges.append((leg+'_in',legcut[leg],vertex,leg));self.ext[legcut[leg]].append((sym,1))
            else:self.ext[vertex].append((sym,1))
        for tip in sorted(tips):self.edges.append(('tip'+tip[1],'C'+tip[1],'H'+tip[1],'tip'+tip[1]))
        self.ext['H1'].append(('Q1',-1));self.ext['H2'].append(('Q2',-1))
        self.edges.append(('ell',ends[0],ends[1],'ell'))
        self.vertices=sorted({v for e in self.edges for v in e[1:3]})
        chord=next(e[0] for e in self.edges if e[3]=='k' and e[2] in ('H1','C1'))
        self.mom={chord:{'K':1},'ell':{'L':1}}
        pending={e[0] for e in self.edges}-set(self.mom)
        while pending:
            before=len(pending)
            for v in self.vertices:
                inc=[e for e in self.edges if v in e[1:3] and e[1]!=e[2]]
                unknown=[e for e in inc if e[0] in pending]
                if len(unknown)!=1:continue
                total=self.external(v)
                for e in inc:
                    if e[0] in pending:continue
                    for q,c in self.mom[e[0]].items():total[q]+=(1 if e[2]==v else -1)*c
                e=unknown[0];sign=-1 if e[2]==v else 1
                self.mom[e[0]]={q:sign*c for q,c in total.items() if c};pending.remove(e[0])
            if len(pending)==before:raise ArithmeticError('Cannot route '+self.name)
        for v in self.vertices:
            total=self.external(v)
            for e in self.edges:
                for q,c in self.mom[e[0]].items():total[q]+=((e[2]==v)-(e[1]==v))*c
            if any(total.values()):raise ArithmeticError('Momentum conservation: '+self.name)

    def external(self,v):
        out=defaultdict(int)
        for sym,sign in self.ext[v]:
            for q,c in ((('p1',1),('p2',1),('Q1',-1)) if sym=='Q2' else ((sym,1),)):out[q]+=sign*c
        return out

def all_graphs():
    return [Graph(A+'|'+B) for i,A in enumerate(SITES) for B in SITES[i:] if A!=B or A in LINES]+[Graph(A+'|'+B) for A,B in TIPS]

def forest_terms(g):
    """(complement monomial, component momentum) for trees and two-forests."""
    V=len(g.vertices);E=len(g.edges)
    for size in (V-1,V-2):
        for selected in combinations(range(E),size):
            parent={v:v for v in g.vertices}
            def find(v):
                while parent[v]!=v:v=parent[v]
                return v
            cyclic=False
            for i in selected:
                _,t,h,_=g.edges[i];rt,rh=find(t),find(h)
                if rt==rh:cyclic=True;break
                parent[rt]=rh
            if cyclic:continue
            comps=defaultdict(list)
            for v in g.vertices:comps[find(v)].append(v)
            mono=tuple(0 if i in selected else 1 for i in range(E))
            if size==V-1 and len(comps)==1:yield 'U',mono,None
            if size==V-2 and len(comps)==2:
                p=defaultdict(int)
                for v in next(iter(comps.values())):
                    for q,c in g.external(v).items():p[q]+=c
                yield 'F',mono,tuple(p[q] for q in ('p1','p2','Q1'))
