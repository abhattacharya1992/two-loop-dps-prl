#!/usr/bin/env python3
"""Independent algebra for the N=4 endpoint comparison and finite DPS boundary.

No production DPS modules are imported. The source inputs include stored-data
transcriptions, not freshly integrated amplitudes. All scope limits are recorded.
Requirements: Python >=3.10, sympy, mpmath. Existing outputs are never overwritten.
"""
from __future__ import annotations
import argparse, hashlib, json, platform
from pathlib import Path
import sympy as s
import mpmath as mp

HERE = Path(__file__).resolve().parent
CHECKS: list[str] = []
P2=s.Symbol('p2', positive=True)
Z={n: s.Symbol(f'z{n}') for n in (2,3,4,5,6,7,8,9)}
EVEN={Z[2]:P2/6,Z[4]:P2**2/90,Z[6]:P2**3/945,Z[8]:P2**4/9450}

def reduce(e): return s.factor(s.expand(e).subs(EVEN))
def equal(name, a, b=0):
    entries=list(a-b) if isinstance(a,s.MatrixBase) else [a-b]
    if any(reduce(e)!=0 for e in entries): raise ArithmeticError(name+': '+str(a-b))
    CHECKS.append(name)
def parse(text): return s.sympify(text,locals={'p2':P2,**{f'z{k}':v for k,v in Z.items()}})
def exp_coeffs(log_coeff,n):
    """Coefficients of exp(sum log_coeff[k] t**k), exact recursive arithmetic."""
    result=[s.Integer(1)]
    for m in range(1,n+1):
        result.append(s.expand(sum(k*log_coeff.get(k,0)*result[m-k] for k in range(1,m+1))/m))
    return result

def n4_checks(data):
    ell,G,X,alpha=s.symbols('ell G X alpha')
    f={int(L):sum(parse(c)*ell**int(k) for k,c in coeff.items()) for L,coeff in data['mhv'].items()}
    g={int(L):sum(parse(c)*ell**int(k) for k,c in coeff.items()) for L,coeff in data['nmhv'].items()}
    comparison=[]
    for L in range(1,6):
        predicted=s.expand(2*sum((-1)**j*s.diff(f[L],ell,j+1) for j in range(2*L-1)))
        equal(f'L{L} endpoint differential equation',g[L]+s.diff(g[L],ell),2*s.diff(f[L],ell))
        for k in range(2*L-1):
            equal(f'L{L} NMHV coefficient ell^{k}',predicted.coeff(ell,k),g[L].coeff(ell,k))
            comparison.append({'loop':L,'power':k,'derived':str(reduce(predicted.coeff(ell,k))), 'difference':'0'})
    z,zb=s.symbols('z zbar',positive=True)
    # Independent Wirtinger variables; do not set zbar=conjugate(z) before differentiating.
    for L in range(1,6):
        A=2*f[L].subs(ell,s.log(z*zb)); H=z/zb*g[L].subs(ell,s.log(z*zb))
        equal(f'L{L} transverse helicity differential identity',s.diff(H,z),s.diff(A,zb))
    coeff=s.symbols('c0:10')
    poly=sum(coeff[j]*ell**j for j in range(10))
    mapped=2*sum((-1)**j*s.diff(poly,ell,j+1) for j in range(9))
    equal('generic degree-nine polynomial inverse',mapped+s.diff(mapped,ell),2*s.diff(poly,ell))
    equal('homogeneous solution is excluded inverse-power term',s.diff(s.exp(-ell),ell)+s.exp(-ell))
    C=s.Symbol('common_normalization')
    equal('endpoint-independent normalization commutes',C*(mapped+s.diff(mapped,ell)),2*s.diff(C*poly,ell))
    defect=s.expand((ell*mapped+s.diff(ell*mapped,ell))-2*s.diff(ell*poly,ell))
    if reduce(defect)==0: raise ArithmeticError('dependent normalization negative control')
    CHECKS.append('endpoint-dependent normalization is not silently allowed')

    # Known 2019 anomalous dimensions, with a=2*g_literature**2 and cosmic rho undone.
    gamma={1:4,2:-4*Z[2],3:22*Z[4],4:-4*Z[3]**2-s.Rational(219,2)*Z[6]}
    gamma3={2:s.Rational(7,8)*Z[4],3:-s.Rational(115,24)*Z[6]-s.Rational(7,3)*Z[3]**2,
            4:s.Rational(1701,64)*Z[8]+s.Rational(65,4)*Z[3]*Z[5]+s.Rational(7,4)*Z[2]*Z[3]**2}
    B={2:3*Z[3],3:-10*Z[5]-2*Z[2]*Z[3],4:s.Rational(175,4)*Z[7]+5*Z[2]*Z[5]+s.Rational(21,2)*Z[3]*Z[4]}
    lk={j:Z[2]*gamma[j]/4+2*gamma3.get(j,0) for j in range(1,5)}
    kappa=exp_coeffs(lk,4)
    target={1:Z[2],2:Z[4]/2,3:-s.Rational(14,3)*Z[3]**2-s.Rational(13,24)*Z[6],
            4:-s.Rational(13,6)*Z[2]*Z[3]**2+s.Rational(65,2)*Z[3]*Z[5]+s.Rational(257,36)*Z[8]}
    for j in range(1,5): equal(f'known representation kappa{j}',kappa[j],target[j])
    kernel=exp_coeffs({j:lk[j]-gamma[j]*X**2/8+B.get(j,0)*X for j in range(1,5)},4)
    mellin_log={1:-2*G,3:-2*Z[3]/3,5:-2*Z[5]/5,7:-2*Z[7]/7,9:-2*Z[9]/9}
    ec=exp_coeffs(mellin_log,9)
    t0=[-s.factorial(j)*ec[j+1]/2 for j in range(9)]
    t2=[s.factorial(j)*sum(ec[:j+1]) for j in range(9)]
    for j in range(9): equal(f'Bessel component identity X^{j}',t2[j]+s.diff(t2[j],G)/2,s.diff(t0[j],G))
    for L in range(1,6):
        w=s.Poly(kernel[L-1],X)
        transform0=s.expand(sum(w.nth(j)*t0[j] for j in range(2*L-1))).subs(G,ell/2).expand()
        transform2=s.expand(sum(w.nth(j)*t2[j] for j in range(2*L-1))).subs(G,ell/2).expand()
        for k in range(1,2*L): equal(f'known Sudakov MHV L{L} ell^{k}',transform0.coeff(ell,k),f[L].coeff(ell,k))
        equal(f'known Sudakov plus endpoint NMHV L{L}',transform2,g[L])

    # Principal symbol of the known helicity factor. The full BFKL integral is not rerun.
    pp,pm,h=s.symbols('p_plus p_minus h',nonzero=True)
    # x(Rp/2)=Rp/2-2g^2/(Rp)+O(R^-3); h=R^-2.
    ratio=(pp/2-2*h/pp)/(pm/2-2*h/pm)
    equal('Zhukovsky principal helicity symbol',s.limit(ratio,h,0),pp/pm)
    equal('Fourier derivative intertwiner',pm*(pp/pm),pp)
    mp.mp.dps=50
    records=[]
    def zh(v):
        root=mp.sqrt(v*v-4)
        if abs(root-v)>abs(-root-v): root=-root
        return (v+root)/2
    for angle in ('0.31','1.13','2.2','-0.61'):
        angle=mp.mpf(angle); err=[]
        for R in (100,1000):
            v=mp.mpf(R)/2*mp.exp(1j*angle)
            exact=zh(v)/zh(mp.conj(v)); leading=mp.exp(2j*angle)
            err.append(abs(exact-leading))
        scaling=err[0]/err[1]
        if not (mp.mpf('99')<scaling<mp.mpf('101')): raise ArithmeticError('principal symbol scaling')
        records.append({'angle':str(angle),'R100_error':mp.nstr(err[0],12),'R1000_error':mp.nstr(err[1],12),'ratio':mp.nstr(scaling,12)})
    CHECKS.append('four numerical Zhukovsky subleading-power checks')
    return {'coefficient_comparisons':comparison,'known_kappa':{str(k):str(reduce(kappa[k])) for k in range(1,5)},
            'B_from_known_virtual_anomalous_dimension':{str(k):str(v) for k,v in B.items()},
            'principal_symbol_checks':records,
            'scope':'leading nonanalytic selected helicity component in the perturbative self-crossing/MRK overlap; not full NMHV or nonperturbative resummation',
            'primary_ancillary_freshly_expanded':False}

def r0_checks(data):
    x,c=s.symbols('x c',real=True); ex,ey,fx,fy=s.symbols('ex ey fx fy')
    def structures(a,b):return [s.eye(2)*a,s.Matrix([[a,0],[b,0]]),s.Matrix([[a,b],[0,0]]),s.Matrix([[a,0],[0,0]])]
    def tensor(st,weights):return sum((w*m for w,m in zip(weights,st)),s.zeros(2))
    st1,st2=structures(ex,ey),structures(fx,fy)
    tree=[1,-1/(1-x),-1/x,0]; tree2=[1,-1/x,-1/(1-x),0]
    rational=[c/6,-c/6,-c/6,0]; spin3_b=[c/6,c/6,c/6,-2*c/3]
    equal('finite rational b-space tensor simplifies',s.Matrix(rational)+s.Matrix(spin3_b),s.Matrix([c/3,0,0,-2*c/3]))
    A,B=tensor(st1,tree),tensor(st2,tree2)
    R1,R2=tensor(st1,rational),tensor(st2,rational)
    T1,T2=tensor(st1,spin3_b),tensor(st2,spin3_b)
    def con(U,V,parity): return s.trace(U*V.T) if parity=='S' else s.trace(U)*s.trace(V)-s.trace(U*V)
    def proj(v):
        v=s.expand(v); u=v.coeff(ex).coeff(fx); w=v.coeff(ey).coeff(fy)
        return {0:s.factor((u+w)/2),2:s.factor(u-w)}
    sigma=x*(1-x); FS=1/sigma+2*sigma; FA=1/sigma-2*sigma; W=(x**4+(1-x)**4)/sigma
    projections={}
    for parity,F,sign in [('S',FS,1),('A',FA,-1)]:
        den=proj(con(A,B,parity))
        for name,U,V in [('R',R1,R2),('spin3',T1,T2),('total',R1+T1,R2+T2)]:
            num=proj(con(U,B,parity)+con(A,V,parity))
            for nu in (0,2):
                actual=s.factor(num[nu]/den[nu])
                target=2*c/(3*F) if nu==0 and name in ('R','total') else (sign*2*c/(3*W) if nu==2 and name in ('spin3','total') else 0)
                equal(f'arbitrary-x {parity} {name} harmonic{nu}',actual,target)
                projections[f'{parity}_{name}_{nu}']=str(actual)
    equal('supersymmetric rational cancellation',c.subs(c,1-4+3))
    equal('CP-even versus CP-odd spin-zero difference',2*c/3*(1/FA-1/FS),8*c*sigma**3/(3*(1-4*sigma**4)))
    equal('CP-even versus CP-odd spin-two difference',-2*c/(3*W)-2*c/(3*W),-4*c/(3*W))
    equal('threshold rational scalar spin-zero', (2*c/(3*FS)).subs(x,s.Rational(1,2)),4*c/27)
    equal('threshold rational pseudoscalar spin-zero',(2*c/(3*FA)).subs(x,s.Rational(1,2)),4*c/21)
    equal('threshold rational scalar spin-two',(2*c/(3*W)).subs(x,s.Rational(1,2)),4*c/3)
    equal('threshold rational pseudoscalar spin-two',(-2*c/(3*W)).subs(x,s.Rational(1,2)),-4*c/3)
    dilog=s.polylog(2,x)+s.polylog(2,1-x)+s.log(x)*s.log(1-x)-s.zeta(2)
    dd=s.diff(dilog,x).subs({s.polylog(1,x):-s.log(1-x),s.polylog(1,1-x):-s.log(x)})
    equal('symmetric weight-two dilogarithm reduction derivative',dd)
    # The integration constant is zero at x->0. This concerns the stated alphabet
    # and constant coefficients; it does not establish the loop integral's alphabet.
    e=s.Symbol('eps'); lx,ly=s.symbols('lx ly',real=True)
    # Exact Gamma cancellation; use recurrence rather than fragile float series.
    Fe=s.exp(-s.EulerGamma*e)*s.gamma(1-2*e)/(s.gamma(1+e)*s.gamma(1-e))
    bubble=s.gamma(e)*s.gamma(1-e)**2/s.gamma(2-2*e)
    simple=s.exp(-s.EulerGamma*e)*s.gamma(1-e)/(e*(1-2*e))
    equal('Gamma recurrence for bubble times Fourier factor',s.combsimp(Fe*bubble/simple),1)
    # Numerator exp(-gamma eps)Gamma(1-eps) has no linear term.
    Q=(1+Z[2]*e**2/2)/(e*(1-2*e))
    def finite(f):return s.series(f,e,0,1).removeO().expand().coeff(e,0)
    qse=finite(-(1-e)*Q*((1+e*lx)+(1+e*ly)))*2
    gse=finite((5-3*e)/(3-2*e)*Q*((1+e*lx)+(1+e*ly)))*2
    nfse=finite(-4*(1-e)/(3-2*e)*Q*((1+e*lx)+(1+e*ly)))*2
    equal('four quark daughter finite self-energies',qse,-4-2*(lx+ly))
    equal('four gluon CA finite self-energies',gse,s.Rational(124,9)+s.Rational(10,3)*(lx+ly))
    equal('four gluon nf finite self-energies',nfse,-s.Rational(80,9)-s.Rational(8,3)*(lx+ly))
    # gamma(-eps)e^-gamma eps = -1/eps-zeta2 eps/2+...;
    # -gamma-psi(-eps)=-1/eps+zeta2 eps+... . Their product times 2 Delta.
    soft=finite(2*(-1/e-Z[2]*e/2)*(-1/e+Z[2]*e))
    equal('finite soft coefficient per Delta',soft,-Z[2])
    jqcf=8-2*Z[2]+(lx-ly)**2+lx+ly
    jqca=1-lx**2-ly**2
    jgca=s.Rational(5,9)-2*Z[2]-2*lx*ly-s.Rational(5,3)*(lx+ly)
    jgnf=s.Rational(20,9)+s.Rational(4,3)*(lx+ly)
    targets={'quark_CF':jqcf,'quark_CA':jqca,'gluon_CA':jgca,'gluon_TFnf':jgnf}
    equal('quark CF assembly from remaining target',2*jqcf+qse+2*soft,12-6*Z[2]+2*(lx-ly)**2)
    equal('quark CA assembly from remaining target',2*jqca-soft,2+Z[2]-2*lx**2-2*ly**2)
    equal('gluon CA assembly from remaining target',2*jgca+gse+soft,s.Rational(134,9)-5*Z[2]-4*lx*ly)
    equal('gluon nf assembly from remaining target',2*jgnf+nfse,-s.Rational(40,9))
    # These j functions are inferred integration targets, NOT independent arbitrary-x loop results.
    mp.mp.dps=50
    records=[]
    worst=mp.mpf(0)
    for xs,cols in data['quark_collinear_finite'].items():
        rat=s.Rational(xs); xx=mp.mpf(int(rat.p))/int(rat.q)
        pred={'CF':8-2*mp.zeta(2)+mp.log(xx/(1-xx))**2+mp.log(xx*(1-xx)),
              'CA':1-mp.log(xx)**2-mp.log(1-xx)**2}
        for col,vals in cols.items():
            err=max(abs(mp.mpf(v)-pred[col]) for v in vals); worst=max(worst,err)
            if err>mp.mpf('2e-7'):raise ArithmeticError('stored collinear coefficients')
            records.append({'x':xs,'colour':col,'max_abs_error':mp.nstr(err,12)})
    CHECKS.append('24 stored quark collinear finite parts support inferred targets')
    # Even a low-weight symmetric function can vanish at all sampled fractions.
    samples=[s.Rational(1,10),s.Rational(3,10),s.Rational(1,2),s.Rational(11,20),s.Rational(5,8),s.Rational(3,4),s.Rational(5,6),s.Rational(9,10)]
    sigmas=sorted(set(t*(1-t) for t in samples))
    ambiguity=s.prod(sigma-t for t in sigmas)
    equal('symmetric reconstruction ambiguity',ambiguity.subs(x,1-x),ambiguity)
    for xv in samples:equal(f'finite-data ambiguity vanishes at {xv}',ambiguity.subs(x,xv))
    if ambiguity.subs(x,s.Rational(2,5))==0:raise ArithmeticError('ambiguity negative control')
    CHECKS.append('finite-data ambiguity survives at unmeasured fraction')
    equal('cusp constant relation is known input',2*((s.Rational(67,9)-2*Z[2]))-Z[2],s.Rational(134,9)-5*Z[2])
    return {'symbolic_rational_projections':projections,
            'self_energies':{'q_CF':str(qse),'g_CA':str(gse),'g_TFnf':str(nfse)},'soft_per_Delta':str(soft),
            'remaining_collinear_targets':{k:str(v) for k,v in targets.items()},
            'target_status':'Inferred by exact subtraction from reconstructed r0; not an independent integration of the common collinear functions.',
            'stored_quark_support':records,'max_abs_stored_quark_error':mp.nstr(worst,12),
            'symmetric_ambiguity':str(ambiguity),
            'full_common_collinear_integrals_analytically_evaluated':False}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():parser.error('output must be new')
    inp=HERE/'inputs.json'; data=json.loads(inp.read_text())
    record={'baseline':data['baseline'],'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'inputs_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),
            'runtime':{'python':platform.python_version(),'sympy':s.__version__,'mpmath':mp.__version__}}
    record['n4']=n4_checks(data);record['r0']=r0_checks(data)
    record['checks']=CHECKS;record['number_of_named_checks']=len(CHECKS);record['passed']=True
    record['not_established']=['complete DPS factorisation','full arbitrary-x evaluation of common collinear loop integrals','new general N=4 all-order law','exhaustive literature priority','fresh raw-amplitude or primary-ancillary regeneration']
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(record,f,indent=2);f.write('\n')
    print(f'PASS: {len(CHECKS)} bounded checks; scope and exclusions recorded in {args.output}')
if __name__=='__main__':main()
