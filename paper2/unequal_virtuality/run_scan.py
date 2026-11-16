#!/usr/bin/env python3
"""Resumable, process-parallel unequal-virtuality VVamp forward scan."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
from uv_core import (HERE, GGVV, PARTS, ALL_HELICITIES, LOG_HELICITIES, SCHEMA,
                     configuration, kinematics, evaluate, cpair,
                     provenance, translated_inputs, sources)


def config_id(config):
    return hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()[:12]


def validate_record(r):
    if r.get("schema") != SCHEMA:
        raise ValueError("unsupported point schema")
    c = configuration(**r["configuration"])
    if c != r["configuration"] or r["dps"] < 30:
        raise ValueError("invalid configuration or precision")
    if r["diagnostics"]["li22_bad_order"] != 0:
        raise ValueError("unvalidated Li22 branch")
    if not r["M"] or not set(r["M"]) <= set(ALL_HELICITIES):
        raise ValueError("invalid/missing helicities")
    provenance.require_pair(r["kinematics"]["G"])
    for n in ("s","ma2","mb2","t","u","pT2","q2"):
        provenance.require_finite(r["kinematics"][n])
    if mp.mpf(r["kinematics"]["pT2"]) <= 0:
        raise ValueError("nonpositive pT2")
    for parts in r["M"].values():
        provenance.require_keys(parts,PARTS,exact=True)
        for value in parts.values():
            provenance.require_pair(value)
    for part in PARTS:
        provenance.require_keys(r["A"][part],[str(j) for j in range(1,21)],exact=True)
        for value in r["A"][part].values():
            provenance.require_pair(value)
    params=r["provenance"]["parameters"]
    if params["configuration"]!=c or params["k"]!=r["k"] or params["dps"]!=r["dps"]:
        raise ValueError("record/provenance parameter mismatch")
    if set(params["helicities"])!=set(r["M"]):
        raise ValueError("incomplete helicity set")


def run_point(g,config,k,dps,output,helicities=LOG_HELICITIES,recompute=False):
    if dps < max(30,5*k+10):
        raise ValueError("require dps >= max(30,5*k+10)")
    with mp.workdps(dps):
        c=configuration(**config)
        hs=tuple(sorted(set(helicities)))
        if not hs or not set(hs)<=set(ALL_HELICITIES):
            raise ValueError("unknown helicity requested")
        inputs=translated_inputs(g.build)
        if provenance.file_hashes(inputs)!=g.input_hashes:
            raise provenance.CacheMismatch("Translated inputs changed after evaluator was loaded")
        params=dict(configuration=c,k=k,dps=dps,helicities=hs,scheme="qT",mu2="s",loop_expansion="alpha_s/(2*pi)")
        meta=provenance.fingerprint("unequal-virtuality-forward",inputs=inputs,sources=sources(),parameters=params,libraries=("mpmath",))
        path=Path(output)/f"uv_{config_id(c)}_d{k}_p{dps}.json"
        old=provenance.load_result(path,meta,recompute=recompute,validator=validate_record)
        if old is not None:
            return str(path),"cached"
        start=time.monotonic()
        kin=kinematics(c,k)
        H,A,diag=evaluate(g,kin,hs)
        rec=dict(schema=SCHEMA,configuration=c,k=k,dps=dps,
                 kinematics={key:cpair(v) if key=="G" else mp.nstr(v,dps) for key,v in kin.items()},
                 M={h:{p:cpair(v) for p,v in vals.items()} for h,vals in H.items()},
                 A={p:{str(j):cpair(v) for j,v in vals.items()} for p,vals in A.items()},
                 diagnostics=diag,seconds=round(time.monotonic()-start,3),provenance=meta)
        if provenance.file_hashes(inputs)!=g.input_hashes:
            raise provenance.CacheMismatch("Translated inputs changed during evaluation")
        validate_record(rec)
        provenance.atomic_write(path,rec)
        return str(path),rec["seconds"]


_G=None

def init_worker(build):
    global _G
    _G=GGVV(build)


def worker(job):
    return run_point(_G,**job)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--s",default="15/2")
    p.add_argument("--ma2",default="1",help="first virtuality: mass SQUARED")
    p.add_argument("--mb2",default="2",help="second virtuality: mass SQUARED")
    p.add_argument("--k",type=int,nargs="+",default=list(range(4,13)))
    p.add_argument("--dps-base",type=int,default=30)
    p.add_argument("--dps-step",type=int,default=5)
    p.add_argument("--dps",type=int,help="override base+step*k, e.g. a precision control")
    p.add_argument("--workers",type=int,default=1,help="processes; start with 1 or 2")
    p.add_argument("--helicities",nargs="+",default=list(LOG_HELICITIES))
    p.add_argument("--all-helicities",action="store_true")
    p.add_argument("--build",type=Path,default=HERE/"build/ggv")
    p.add_argument("--output",type=Path,default=HERE/"runs/stage1")
    p.add_argument("--recompute",action="store_true")
    args=p.parse_args()
    if not 1<=args.workers<=32:
        p.error("workers must be between 1 and 32")
    c=configuration(args.s,args.ma2,args.mb2)
    hs=ALL_HELICITIES if args.all_helicities else args.helicities
    jobs=[dict(config=c,k=k,dps=args.dps if args.dps is not None else args.dps_base+args.dps_step*k,
               output=args.output,helicities=hs,recompute=args.recompute) for k in sorted(set(args.k))]
    for job in jobs:
        with mp.workdps(max(30,job["dps"])):
            kinematics(c,job["k"])
        if job["dps"]<max(30,5*job["k"]+10):
            p.error("require dps >= max(30,5*k+10)")
    if args.workers==1:
        init_worker(args.build)
        for job in jobs:
            print(worker(job),flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers,initializer=init_worker,initargs=(args.build,)) as pool:
            for f in as_completed([pool.submit(worker,j) for j in jobs]):
                print(f.result(),flush=True)


if __name__=="__main__":
    main()
