#!/usr/bin/env python3
"""Unconstrained complex fits; predictions are compared only AFTER extraction.

Fit y=sum(f_j G^j,j=0..J)+rho*sum(h_j G^j,j=0..P),
rho=pT2/sqrt(ma2*mb2). P=-1 disables the power correction.
Window shifts and model variations are diagnostics, not statistical errors.
"""
import argparse
import gzip
import json
from pathlib import Path
import mpmath as mp
from uv_core import (HERE, LOG_HELICITIES, PARTS, complex_pair, cpair, predictions,
                     provenance, real, kinematics)
from run_scan import validate_record


def load_records(path):
    path=Path(path)
    if path.suffix == ".csv" or path.name.endswith(".csv.gz"):
        from stored_points import read_csv
        records=read_csv(path)
    else:
        files=sorted(path.glob("uv_*.json")) if path.is_dir() else [path]
        if not files:
            raise ValueError("no uv_*.json point files")
        records=[]
        for f in files:
            text=gzip.open(f,"rt").read() if f.suffix==".gz" else f.read_text()
            data=json.loads(text)
            values=data if isinstance(data,list) else [data]
            for r in values:
                validate_record(r)
                records.append(r)
    configs={json.dumps(r["configuration"],sort_keys=True) for r in records}
    if len(configs)!=1:
        raise ValueError("Do not mix physical configurations in one fit")
    # At a repeated k, retain highest precision and report the discrepancy.
    groups={}
    for r in records:
        groups.setdefault(r["k"],[]).append(r)
    chosen=[]
    controls=[]
    for k,group in sorted(groups.items()):
        group.sort(key=lambda r:r["dps"])
        best=group[-1]
        for low in group[:-1]:
            hs=set(low["M"]) & set(best["M"])
            deviation={}
            for part in PARTS:
                deviation[part]=max(abs(complex_pair(best["M"][h][part])-complex_pair(low["M"][h][part])) /
                                    max(mp.mpf(1),abs(complex_pair(best["M"][h][part]))) for h in hs)
            controls.append(dict(k=k,dps_low=low["dps"],dps_high=best["dps"],
                                 scaled_difference={p:mp.nstr(v,15) for p,v in deviation.items()}))
        chosen.append(best)
    return chosen,controls


def fit(records,values,degree,pdegree,window=0,extra=0):
    n=degree+pdegree+2 if pdegree>=0 else degree+1
    end=len(records)-window
    count=n+extra
    if end<count:
        return None
    chosen=records[end-count:end]
    vals=values[end-count:end]
    rows=[]
    for r in chosen:
        G=complex_pair(r["kinematics"]["G"])
        rho=mp.mpf(r["kinematics"]["pT2"])/mp.mpf(r["kinematics"]["q2"])
        rows.append([G**j for j in range(degree+1)]+[rho*G**j for j in range(pdegree+1)])
    # Column scaling improves numerical conditioning without constraining coefficients.
    scale=[max(abs(row[j]) for row in rows) for j in range(n)]
    if any(x==0 for x in scale):
        raise ValueError("rank-deficient fit column")
    M=mp.matrix([[v/scale[j] for j,v in enumerate(row)] for row in rows])
    y=mp.matrix(vals)
    sol=mp.lu_solve(M,y) if count==n else mp.qr_solve(M,y)[0]
    coeff=[sol[j]/scale[j] for j in range(n)]
    residual=max(abs(sum(row[j]*coeff[j] for j in range(n))-value) for row,value in zip(rows,vals))
    residual/=max(mp.mpf(1),max(abs(v) for v in vals))
    return dict(coeff=coeff[:degree+1],nuisance=coeff[degree+1:],
                k=[r["k"] for r in chosen],points=count,parameters=n,
                scaled_residual=residual)


def analyse(records,helicity="++++"):
    if helicity not in LOG_HELICITIES:
        raise ValueError("This analysis targets the isotropic one-loop logarithm; select ++++, ++--, ---- or --++")
    if any(helicity not in r["M"] for r in records):
        raise ValueError("selected helicity not stored at every point")
    if len(records)<5:
        raise ValueError("need at least five distinct k values; nine are recommended")
    # Recompute G and rho from exact configurations; do not trust rounded fitted abscissae.
    checked=[]
    for r in records:
        kin=kinematics(r["configuration"],r["k"])
        tolerance=mp.power(10,-min(r["dps"]-15,mp.mp.dps-10))
        for name in ("t","u","pT2"):
            if abs(mp.mpf(r["kinematics"][name])-kin[name])>tolerance*max(1,abs(kin[name])):
                raise ValueError(f"inconsistent {name} at k={r['k']}")
        rr=dict(r,kinematics=dict(r["kinematics"],G=cpair(kin["G"]),q2=mp.nstr(kin["q2"]),pT2=mp.nstr(kin["pT2"])))
        # Use full analysis precision for dimensionful abscissae.
        rr["kinematics"]["q2"]=mp.nstr(kin["q2"],mp.mp.dps)
        rr["kinematics"]["pT2"]=mp.nstr(kin["pT2"],mp.mp.dps)
        checked.append(rr)
    records=checked
    values={p:[complex_pair(r["M"][helicity][p]) for r in records] for p in PARTS}
    one=fit(records,values["1l"],1,2)
    onew=fit(records,values["1l"],1,2,window=1)
    c1=one["coeff"][1]
    if abs(c1)<mp.power(10,-min(20,mp.mp.dps//2)):
        raise ValueError("one-loop logarithm too small for this normalization")
    pred=predictions(records[0]["configuration"])
    result=dict(configuration=records[0]["configuration"],helicity=helicity,distinct_points=len(records),
                c1=cpair(c1,35),c1_window_shift=mp.nstr(abs(c1-onew["coeff"][1]),15) if onew else None,
                assumption="r1=r2=0, tested rather than independently proved for unequal virtualities",
                fits=[])
    for P in (-1,0,1,2,3):
        for part in ("2l_CF","2l_CA"):
            f=fit(records,values[part],3,P)
            fw=fit(records,values[part],3,P,window=1)
            fl=fit(records,values[part],3,P,extra=1)
            if f is None:
                continue
            row=dict(part=part,power_correction_degree=P,k=f["k"],points=f["points"],parameters=f["parameters"],
                     scaled_fit_residual=mp.nstr(f["scaled_residual"],15),coefficients=[cpair(v,35) for v in f["coeff"]],comparison={})
            for j in (3,2):
                ratio=f["coeff"][j]/c1
                shift=abs(f["coeff"][j]-fw["coeff"][j])/abs(c1) if fw else None
                # Also vary the one-loop extraction independently.
                norm_shift=abs(f["coeff"][j]/onew["coeff"][1]-ratio) if onew else None
                row["comparison"][str(j)]=dict(extracted=cpair(ratio,35),prediction=mp.nstr(pred[part][j],35),
                    absolute_deviation=mp.nstr(abs(ratio-pred[part][j]),15),
                    window_shift=mp.nstr(shift,15) if shift is not None else None,
                    normalization_shift=mp.nstr(norm_shift,15) if norm_shift is not None else None,
                    least_squares_shift=mp.nstr(abs(fl["coeff"][j]/c1-ratio),15) if fl else None)
            result["fits"].append(row)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("points",type=Path)
    p.add_argument("--helicity",default="++++")
    p.add_argument("--output",type=Path)
    p.add_argument("--dps",type=int,default=110,help="fit arithmetic precision; does not increase input accuracy")
    args=p.parse_args()
    mp.mp.dps=args.dps
    records,controls=load_records(args.points)
    result=analyse(records,args.helicity)
    result["precision_controls"]=controls
    result["analysis_dps"]=args.dps
    result["input_sha256"]=provenance.file_hashes({str(f.name):f for f in sorted(args.points.glob("uv_*.json"))}) if args.points.is_dir() else provenance.file_hashes({args.points.name:args.points})
    result["analysis_source_sha256"]=provenance.file_hashes({n:HERE/n for n in ("analyse_scan.py","uv_core.py","stored_points.py")})
    if args.output:
        provenance.atomic_write(args.output,result)
    print(f"configuration={result['configuration']}, helicity={args.helicity}, c1={result['c1']}")
    for row in result["fits"]:
        print(f"P={row['power_correction_degree']:2} {row['part']} k={row['k']}")
        for j,v in row["comparison"].items():
            print(f"  G^{j}/c1={v['extracted']}, prediction={v['prediction']}, deviation={v['absolute_deviation']}, window={v['window_shift']}")
    if controls:
        print("Precision controls:",controls)


if __name__=="__main__":
    main()
