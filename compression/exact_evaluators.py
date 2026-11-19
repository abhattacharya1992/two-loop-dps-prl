#!/usr/bin/env python3
"""Time the ORIGINAL HH/ZZ full-amplitude evaluators at fixed scattering angle.

These are exact published-input backends, not new beta-series approximants.
External HH Gauge{1,2}Simp.in / ZZ translated VVamp inputs are mandatory.
The default is non-executing: pass --run-exact to invoke expensive evaluations.
Each process should run in its own OS process for peak RSS interpretation.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import importlib
import json
from pathlib import Path
import resource
import statistics
import sys
import time

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]


def fraction(text: str) -> mp.mpf:
    if "/" in str(text):
        a, b = str(text).split("/", 1)
        return mp.mpf(a) / mp.mpf(b)
    return mp.mpf(text)


def fixed_angle_kinematics(beta: str, costheta: str, *, digits: int) -> dict:
    """Equal external mass^2=1, s=4/(1-beta^2), t/u from fixed cos(theta)."""
    if not 20 <= digits <= 2000:
        raise ValueError("digits must be between 20 and 2000")
    with mp.workdps(digits):
        b, c = fraction(beta), fraction(costheta)
        if not 0 < b < 1 or not 0 <= c < 1:
            raise ValueError("require 0 < beta < 1 and 0 <= cos(theta) < 1")
        s = 4 / (1 - b*b)
        t = 1 - s/2 + s*b*c/2
        u = 1 - s/2 - s*b*c/2
        pt2 = (s/4 - 1)*(1 - c*c)
        x = (1-b)/(1+b)
        y = -t
        if abs(s+t+u-2) > mp.mpf(10)**(5-digits):
            raise AssertionError("Mandelstam identity failed")
        if pt2 <= 0:
            raise ValueError("Transverse momentum must be nonzero")
        return {k:mp.nstr(v,digits) for k,v in
                {"beta":b,"costheta":c,"s":s,"t":t,"u":u,
                 "pT2":pt2,"x":x,"y":y}.items()}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@contextmanager
def prepend_python(path: Path):
    p = str(path.resolve())
    sys.path.insert(0,p)
    try:
        yield
    finally:
        sys.path.remove(p)


def exact_hh(kin: dict, dps: int, repeats: int, directory: Path) -> dict:
    if not (directory/"anc/Gauge1Simp.in").is_file() or not (directory/"anc/Gauge2Simp.in").is_file():
        raise FileNotFoundError(
            "HH published inputs absent: obtain Gauge1Simp.in and Gauge2Simp.in "
            "through the documented paper2/hh fetch target")
    with prepend_python(directory):
        # Preserve the AST depth needed by the original large expressions.
        sys.setrecursionlimit(max(sys.getrecursionlimit(), 100000))
        ev=importlib.import_module("evaluate_hh")
        with mp.workdps(dps):
            started=time.perf_counter()
            codes={i:ev.compile_amp(i) for i in (1,2)}
            prep=time.perf_counter()-started
            x=mp.mpf(kin["x"]);y=mp.mpf(kin["y"])
            times=[];values={}
            for i in range(repeats):
                t0=time.perf_counter()
                values={j:mp.j*ev.eval_amp(codes[j],x,y,M2=1,
                                          muR2=mp.mpf("0.5"),N=3,Nf=5)/4
                        for j in (1,2)}
                times.append(time.perf_counter()-t0)
            return {
                "process":"HH","scheme":"published class-B, A_j=i Gauge_j/(4 M2), muR2/M2=1/2",
                "amplitudes":2,"input_sha256":{f"Gauge{i}Simp.in":sha256(directory/"anc"/f"Gauge{i}Simp.in") for i in (1,2)},
                "input_total_bytes":sum((directory/"anc"/f"Gauge{i}Simp.in").stat().st_size for i in (1,2)),
                "preparation_seconds":prep,"evaluation_seconds":times,
                "evaluation_median_seconds":statistics.median(times),
                "values":{f"A{j}":mp.nstr(values[j],min(dps,40)) for j in (1,2)},
            }


def exact_zz(kin: dict, dps: int, repeats: int, build: Path, scripts: Path) -> dict:
    missing=[p for p in (build/"p0002.py",build/"algebra.marshal",
                         build/"rationals.json",build/"coefficients.json") if not p.is_file()]
    if missing:
        raise FileNotFoundError("ZZ translated VVamp inputs absent: "
                                + ", ".join(str(p) for p in missing)
                                + "; follow paper2/zz/README.md")
    with prepend_python(scripts):
        zz=importlib.import_module("zz_points")
        gg=importlib.import_module("ggvv")
        with mp.workdps(dps):
            started=time.perf_counter()
            evaluator=gg.GGVV(build)
            prep=time.perf_counter()-started
            s=mp.mpf(kin["s"]);t=mp.mpf(kin["t"]);u=mp.mpf(kin["u"])
            kt={"s":s,"t":t,"u":u,"pT2":mp.mpf(kin["pT2"])}
            times=[];data={}
            for i in range(repeats):
                t0=time.perf_counter()
                fn=gg.Functions()
                At=evaluator.compute(s,t,1,1,fn)
                Au=evaluator.compute(s,u,1,1,fn)
                A={part:zz.all_coefficients({j:At[j][part] for j in zz.DIRECT},
                                            {j:Au[j][part] for j in zz.DIRECT})
                   for part in zz.PARTS}
                data=zz.helicity_amplitudes(A,kt)
                times.append(time.perf_counter()-t0)
            return {
                "process":"ZZ","scheme":"class-A massless-quark loops, qT scheme, mu^2=s",
                "helicity_count":len(data),"colour_parts":list(zz.PARTS),
                "input_sha256":{p.name:sha256(p) for p in missing if p.exists()} if missing else
                                {name:sha256(build/name) for name in
                                 ("p0002.py","algebra.marshal","rationals.json","coefficients.json")},
                "input_total_bytes":sum(p.stat().st_size for p in
                    (build/"p0002.py",build/"algebra.marshal",
                     build/"rationals.json",build/"coefficients.json")),
                "preparation_seconds":prep,"evaluation_seconds":times,
                "evaluation_median_seconds":statistics.median(times),
                "sample_++++":{part:mp.nstr(data["++++"][part],min(dps,40))
                               for part in zz.PARTS},
            }


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--process",required=True,choices=("HH","ZZ"))
    parser.add_argument("--run-exact",action="store_true",
                        help="required consent to expensive external-input evaluation")
    parser.add_argument("--beta",default="1/10")
    parser.add_argument("--costheta",default="0")
    parser.add_argument("--digits",type=int,default=50)
    parser.add_argument("--repeats",type=int,default=3)
    parser.add_argument("--hh-directory",type=Path,default=ROOT/"paper2/hh/hh_two_loop")
    parser.add_argument("--zz-build",type=Path,default=ROOT/"paper2/zz/build/ggv")
    parser.add_argument("--output",type=Path,default=ROOT/"compression-results/exact")
    args=parser.parse_args(argv)
    if not args.run_exact:
        parser.error("pass --run-exact to evaluate external HH/ZZ amplitudes")
    if not 1 <= args.repeats <= 10:
        parser.error("repeats must be 1 to 10")
    kin=fixed_angle_kinematics(args.beta,args.costheta,digits=args.digits)
    result=exact_hh(kin,args.digits,args.repeats,args.hh_directory) if args.process=="HH" else \
           exact_zz(kin,args.digits,args.repeats,args.zz_build,ROOT/"paper2/zz/scripts")
    report={"schema_version":1,"created_utc":datetime.now(timezone.utc).isoformat(),
            "scope":"original exact amplitude only, not threshold beta-series",
            "kinematics":kin,"digits":args.digits,"repeats":args.repeats,
            "peak_rss_kib_process":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            "result":result}
    args.output.mkdir(parents=True,exist_ok=True)
    out=args.output/(args.process.lower()+"-exact.json")
    out.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
