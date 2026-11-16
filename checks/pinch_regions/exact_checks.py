#!/usr/bin/env python3
"""Exact Symanzik identities and positive-interior stationary exclusions.
These are statements about scalar graph polynomials, not regulated amplitude
integrals. All incoming fractions remain symbolic with 0<x,y<1.
"""
from __future__ import annotations
import json
from pathlib import Path
import sympy as s
from graphs import Graph, all_graphs, forest_terms
from param_space import require
x,y,t,u,v=s.symbols('x y t u v')

def expressions(g):
    aa=s.symbols('a:'+str(len(g.edges)));U=0;F=0
    for kind,m,p in forest_terms(g):
        term=s.prod(a**k for a,k in zip(aa,m))
        if kind=='U':U+=term
        else:
            a,b,c=p;F-=4*(a*b+a*c*y+b*c*x+c*c*x*y-c*c*t)*term
    return aa,s.expand(U),s.expand(F)

def matrix_polynomials(g,aa):
    gram=s.Matrix([[0,2,2*y],[2,0,2*x],[2*y,2*x,4*x*y-4*t]])
    A=s.zeros(2);B=s.zeros(2,3);C=0
    for alpha,edge in zip(aa,g.edges):
        q=g.mom[edge[0]];cl=s.Matrix([q.get('K',0),q.get('L',0)]);ce=s.Matrix([q.get(k,0) for k in ('p1','p2','Q1')])
        A+=alpha*cl*cl.T;B+=alpha*cl*ce.T;C+=alpha*(ce.T*gram*ce)[0]
    U=A.det();F=-U*C+sum(A.adjugate()[i,j]*(B.row(i)*gram*B.row(j).T)[0] for i in range(2) for j in range(2))
    return s.expand(U),s.expand(F)

def strict_sign(expr):
    expr=s.cancel(expr.subs({x:u/(1+u),y:v/(1+v)}));num,den=s.fraction(expr)
    def sign(poly):
        symbols=sorted(poly.free_symbols,key=str)
        cs=s.Poly(poly,*symbols).coeffs() if symbols else [poly]
        if cs and all(c>0 for c in cs):return 1
        if cs and all(c<0 for c in cs):return -1
        return 0
    return sign(num)*sign(den)

def no_stationary(g,F,aa):
    F0=F.subs(t,0);grad=[s.diff(F0,a) for a in aa]
    for a,d in zip(aa,grad):
        sign=strict_sign(d)
        if sign:return {'graph':g.name,'method':'strict derivative','derivative':str(a),'sign':sign}
    for a,d in zip(aa,grad):
        coeff,fs=s.factor_list(d,*aa)
        if not strict_sign(coeff):continue
        mixed=[f for f,n in fs if not strict_sign(f)]
        if not mixed:continue
        branches=[]
        for f in mixed:
            found=None
            for b in aa:
                if s.degree(f,b)!=1:continue
                sol=s.cancel(-f.subs(b,0)/s.diff(f,b))
                if strict_sign(sol)!=1:continue
                for c,dd in zip(aa,grad):
                    reduced=s.cancel(dd.subs(b,sol));sign=strict_sign(reduced)
                    if sign:
                        found={'factor':str(f),'solve_for':str(b),'solution':str(sol),'derivative':str(c),'sign':sign};break
                if found:break
            if found is None:break
            branches.append(found)
        if len(branches)==len(mixed):return {'graph':g.name,'method':'stationary branches excluded','necessary_derivative':str(a),'branches':branches}
    raise ArithmeticError('No exact stationary-exclusion certificate for '+g.name)

def box_check():
    class Box:
        vertices=['a1','a2','H1','H2']
        edges=[('k','a1','H1','k'),('kb','a1','H2','kb'),('r','a2','H1','r'),('rb','a2','H2','rb')]
        def external(self,v):
            return {'a1':{'p1':1},'a2':{'p2':1},'H1':{'Q1':-1},'H2':{'p1':-1,'p2':-1,'Q1':1}}[v]
    aa,U,F=expressions(Box());a,b,c,d=aa
    f1=x*a-(1-x)*b;f2=y*c-(1-y)*d
    require(s.expand(F/4+f1*f2-t*(a+b)*(c+d))==0,'box identity')
    pinch={a:(1-x)*b/x,c:(1-y)*d/y}
    require(F.subs(t,0).subs(pinch).simplify()==0,'box zero')
    for q in aa:require(s.diff(F.subs(t,0),q).subs(pinch).simplify()==0,'box gradient')
    return {'identity':'F/s=-f1*f2+t*(alpha_k+alpha_kb)*(alpha_r+alpha_rb)','domain':'0<x,y<1; all alpha>0','scalar_power_d4':-2,'numerator_included':False}

def run():
    interjet={'k','kb','leg1','a1'},{'r','rb','leg2','a2'}
    certificates=[];count=0
    for g in all_graphs():
        aa,U,F=expressions(g);Um,Fm=matrix_polynomials(g,aa)
        require(s.expand(U-Um)==0 and s.expand(F-Fm)==0,'forest/matrix disagreement '+g.name);count+=1
        A,B=g.name.split('|')
        if (A in interjet[0] and B in interjet[1]) or (B in interjet[0] and A in interjet[1]):
            certificates.append(no_stationary(g,F,aa))
    require(count==54 and len(certificates)==16,'graph census changed')
    return {'one_loop':box_check(),'canonical_graph_polynomial_identities':count,'interjet_no_stationary_certificates':certificates,'scope':'All symbolic interior fractions, strictly positive Feynman parameters, equal virtualities of the 16 interjet-exchange graphs. Does not exclude boundary pinches or regulator-defined Glauber phases.'}

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);args=ap.parse_args();text=json.dumps(run(),indent=2,sort_keys=True)+'\n'
    if args.output:args.output.write_text(text)
    print(text,end='')
