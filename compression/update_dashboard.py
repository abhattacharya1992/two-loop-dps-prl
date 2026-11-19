#!/usr/bin/env python3
"""Update tracked performance dashboard from completed, validated CI artifacts.

Run on GitHub Actions after successful compression or manual amplitude benchmark.
Does not perform the expensive amplitude calculations or access the network.
"""
import argparse
import copy
import json
import math
from pathlib import Path
import sys

DIR=Path(__file__).resolve().parent
SNAP=DIR/"benchmark_snapshot.json"
REPO="https://github.com/abhattacharya1992/two-loop-dps-prl"
LABELS=("HH","AA","ZZ")

def check_number(x, label, allow_zero=False):
    if not isinstance(x,(int,float)) or not math.isfinite(float(x)) or (x<0 if allow_zero else x<=0):
        raise ValueError("Invalid "+label)
    return float(x)

def find_reports(directory):
    found={}
    for file in directory.rglob("*.json"):
        if file.name not in ("report.json","stage2.json"):
            continue
        obj=json.loads(file.read_text())
        key=None
        if file.name=="report.json" and obj.get("scope")=="lossless stored amplitude-point JSON only":
            key="storage"
        elif file.name=="report.json" and obj.get("historical",{}).get("mode")=="published_historical_measurements":
            key="historical"
        elif file.name=="stage2.json" and obj.get("process") in ("hh","aa","zz"):
            key="exact"
        if key:
            if key in found:
                raise ValueError("Multiple "+key+" benchmark reports")
            found[key]=obj
    if not found:
        raise ValueError("No recognised verified benchmark reports")
    return found

def update(snapshot,reports,run_id):
    data=copy.deepcopy(snapshot)
    link=REPO+"/actions/runs/"+str(run_id)
    if "storage" in reports:
        rep=reports["storage"]
        if rep.get("compression_level")!=9:
            raise ValueError("Storage report must use zlib level 9")
        for name in LABELS:
            p=rep["processes"][name]
            if p.get("verified_byte_identical") is not True:
                raise ValueError(name+": SHA-256 lossless check missing")
            raw=check_number(p["raw_bytes"],"raw storage")
            packed=check_number(p["compressed_bytes"],"packed storage")
            if packed>=raw or abs(p["percent_smaller"]-100*(1-packed/raw))>.002:
                raise ValueError(name+": inconsistent compression totals")
            data["point_storage"]["items"][name]={
                "points":int(check_number(p["records"],"records")),
                "original_bytes":int(raw),"packed_bytes":int(packed),
                "raw_parse_ms":check_number(p["raw_json_parse_median_ms"],"parse time"),
                "verified_packed_parse_ms":check_number(p["compressed_verify_and_parse_median_ms"],"verified parse time")}
        data["point_storage"]["ci"]=link
    if "historical" in reports:
        rep=reports["historical"]
        if len(rep.get("upstream_git_revision",""))!=40:
            raise ValueError("AA upstream SHA is not pinned")
        hist=rep["historical"]
        row=next((x for x in hist["timing_rows"] if x["beta"]=="1/10" and x["cos_theta"]=="0" and x["working_digits"]==80),None)
        if row is None:
            raise ValueError("Expected AA historical timing point absent")
        exact=check_number(row["warm_median_seconds"]["exact"],"exact time")
        series=check_number(row["warm_median_seconds"]["beta6"],"beta6 time")
        if abs(exact/series-row["reported_warm_speedup"])>.0001:
            raise ValueError("AA reported warm speedup is inconsistent")
        aa=data["aa_analytic"]
        aa["original_bytes"]=int(check_number(hist["stored_expression_bytes"]["exact"],"AA original bytes"))
        aa["beta6_bytes"]=int(check_number(hist["stored_expression_bytes"]["beta6"],"AA beta6 bytes"))
        aa["upstream_source_commit"]=rep["upstream_git_revision"]
        aa["upstream_historical_ci"]=link
        aa["upstream_hashes"]=hist["input_sha256"]
        aa["historical_warm"].update({"exact_seconds":exact,"beta6_seconds":series,
            "beta":row["beta"],"cos_theta":row["cos_theta"],
            "digits":row["working_digits"],"repeats":row["warm_repeats"]})
    if "exact" in reports:
        obj=reports["exact"]
        kind=obj["process"].upper()
        if kind=="AA":
            gate=obj.get("accuracy_gate",{})
            if gate.get("passed") is not True:
                raise ValueError("AA numerical accuracy gate not passed")
            errs={}
            for name in ("M1","M2"):
                v=check_number(obj["relative_error_"+name],name+" relative error",True)
                if v>gate["max_relative_error"]:
                    raise ValueError(name+" accuracy failure")
                errs[name]=v
            ref=obj.get("independent_exact_reference",{})
            if ref.get("status")=="compared" and any(
                v>gate["max_relative_error"] for v in ref["exact_relative_errors"].values()):
                raise ValueError("AA archived exact reference failure")
            data["aa_analytic"]["fresh_accuracy"]={
                "beta":obj["configuration"]["beta"],
                "cos_theta":obj["configuration"]["cos_theta"],
                "digits":obj["configuration"]["digits"],
                "exact_seconds":check_number(obj["exact_median_seconds"],"AA exact time"),
                "beta6_seconds":check_number(obj["beta6_median_seconds"],"AA beta6 time"),
                "relative_errors":errs,
                "single_run_ci":link,
                "note":"Fresh one-shot Python run; not a warm repeated benchmark"}
        else:
            ref=obj.get("archived_reference",{})
            if (obj.get("precision_validation_status")!="passed" or
                    ref.get("pass") is not True or
                    obj.get("analytic_beta_series_status")!="not_derived"):
                raise ValueError(kind+" lacks certified numerical accuracy, or claims an unsupported expansion")
            err=check_number(obj["max_precision_convergence_error"],"convergence",True)
            if err>1e-10:
                raise ValueError(kind+" insufficient precision convergence")
            data["exact_precision"][kind].update({
                "digits":obj["precision_sweep"]["precision_levels"],
                "verified_min_digits":obj["accuracy_preserving_working_digits"],
                "low_precision_eval_seconds":check_number(obj["accuracy_preserving_evaluation_seconds"],"exact runtime"),
                "top_precision_convergence_error":err,
                "archived_reference_error":check_number(ref["max_relative_error"],"reference",True),
                "reference_is_independent_implementation":bool(ref["independent_algorithm"]),
                "ci":link})
    return data

def draw(data):
    aa=data["aa_analytic"]
    points=data["point_storage"]["items"]
    ratio=aa["original_bytes"]/aa["beta6_bytes"]
    timing=aa["historical_warm"]
    speed=timing["exact_seconds"]/timing["beta6_seconds"]
    from html import escape
    shapes=['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 830" width="1200" height="830">']
    def box(x,y,w,h,color,r=0):
        shapes.append('<rect x="%s" y="%s" width="%s" height="%s" rx="%s" fill="%s"/>'%(x,y,w,h,r,color))
    def text(x,y,value,size=16,color="#15314c",weight=400):
        shapes.append('<text x="%s" y="%s" fill="%s" font-size="%s" font-weight="%s" font-family="Arial,sans-serif">%s</text>'%(x,y,color,size,weight,escape(str(value))))
    navy,blue,green="#142943","#397bbc","#168b80"
    box(0,0,1200,830,"#f3f7fa");box(0,0,1200,100,navy)
    text(42,50,"COMPRESSION & PERFORMANCE",31,"white",700)
    text(42,80,"Validated AA analytical expansion and lossless HH / AA / ZZ storage",15,"#d9e5f0")
    box(38,123,546,413,"white",14);box(605,123,557,413,"white",14)
    text(58,165,"AA: analytic beta^6 threshold series",22,navy,700)
    text(58,194,"Matching Wolfram expression representations",14,"#65758a")
    text(58,239,"Expression size: %.2fx smaller"%ratio,21,green,700)
    text(58,266,"Original: %.2f MB"%(aa["original_bytes"]/1e6),14)
    box(58,276,493,16,blue,7)
    text(58,322,"Series: %.3f MB"%(aa["beta6_bytes"]/1e6),14)
    box(58,334,493,"16","#e5ecf2",7);box(58,334,max(6,493/ratio),16,green,7)
    text(58,388,"Repeated warm evaluation: %.2fx faster"%speed,20,green,700)
    text(58,417,"Original: %.2f ms"%(timing["exact_seconds"]*1000),14)
    box(58,430,493,16,blue,7)
    text(58,476,"Series: %.2f ms"%(timing["beta6_seconds"]*1000),14)
    box(58,487,493,16,"#e5ecf2",7);box(58,487,max(6,493/speed),16,green,7)
    text(627,165,"HH / AA / ZZ: stored-point compression",21,navy,700)
    text(627,194,"zlib level 9; SHA-256 byte-identical checks",14,"#65758a")
    for i,name in enumerate(LABELS):
        p=points[name]
        smaller=100*(1-p["packed_bytes"]/p["original_bytes"])
        y=246+i*92
        text(627,y,name+" (%s points)"%p["points"],19,navy,700)
        text(951,y,"%.1f%% smaller"%smaller,19,green,700)
        box(627,y+14,508,23,"#e5ecf2",8)
        rest=508*p["packed_bytes"]/p["original_bytes"]
        box(627,y+14,rest,23,blue,8);box(627+rest,y+14,508-rest,23,green)
        text(627,y+60,"%.1f KiB -> %.1f KiB"%(p["original_bytes"]/1024,p["packed_bytes"]/1024),13,"#65758a")
    box(38,558,1124,213,"white",14)
    text(59,602,"Numerical checks and scientific scope",22,navy,700)
    a=aa["fresh_accuracy"];h=data["exact_precision"]["HH"];z=data["exact_precision"]["ZZ"]
    text(61,646,"AA exact-to-beta^6 errors: M1 %.3g; M2 %.3g"%(a["relative_errors"]["M1"],a["relative_errors"]["M2"]),17)
    text(61,685,"HH: %d-digit exact evaluation passes stability (no analytic beta series)"%h["verified_min_digits"],16)
    text(61,724,"ZZ: %d-digit exact evaluation passes stability (no analytic beta series)"%z["verified_min_digits"],16)
    text(45,804,"Historical warm AA times are not equivalent to fresh one-shot timings. Compressed storage parsing is slower.",13,"#65758a")
    shapes.append("</svg>")
    return "\n".join(shapes)+"\n"

def table(data):
    aa=data["aa_analytic"];t=aa["historical_warm"];p=data["point_storage"]["items"]
    ratio=aa["original_bytes"]/aa["beta6_bytes"];speed=t["exact_seconds"]/t["beta6_seconds"]
    lines=["# Verified compression and performance results","",
           "![Performance overview](00-performance-overview.svg)","",
           "These benchmarks distinguish AA analytical beta-series reduction from lossless",
           "stored point-data compression; HH/ZZ analytical beta series are not derived.","",
           "## AA analytical beta^6 threshold series","",
           "| Metric | Original exact | Beta^6 series | Ratio |",
           "|---|---:|---:|---:|",
           "| Wolfram representation bytes | %s | %s | **%.2fx smaller** |"%(f"{aa['original_bytes']:,}",f"{aa['beta6_bytes']:,}",ratio),
           "| Historical warm median | %.3f ms | %.3f ms | **%.2fx faster** |"%(t["exact_seconds"]*1000,t["beta6_seconds"]*1000,speed),
           "",
           "Historical repeat timing from the pinned ggAA-threshold package (Wolfram 13.2,",
           "80 digits, beta=0.1, cos(theta)=0). Not a fresh one-shot Python speedup.",
           "[Upstream source](%s/tree/%s)"%(aa["upstream_source_repo"],aa["upstream_source_commit"]),
           "[Verification CI](%s)"%aa["upstream_historical_ci"],"",
           "## Lossless numerical-point storage (zlib + SHA-256)","",
           "| Process | Points | Raw bytes | Packed bytes | Smaller | Raw parse (ms) | Verified packed parse (ms) |",
           "|---|---:|---:|---:|---:|---:|---:|"]
    for name in LABELS:
        row=p[name];pct=100*(1-row["packed_bytes"]/row["original_bytes"])
        lines.append("| %s | %d | %d | %d | **%.2f%%** | %.4f | %.4f |"%(
            name,row["points"],row["original_bytes"],row["packed_bytes"],pct,
            row["raw_parse_ms"],row["verified_packed_parse_ms"]))
    lines+=["","Lossless byte-identical storage savings; decompression plus verified parsing",
            "is slower than raw JSON parsing. No HH/ZZ beta expansion is implied.",
            "[Storage CI](%s)"%data["point_storage"]["ci"],"",
            "## Fresh numerical checks",""]
    a=aa["fresh_accuracy"]
    lines+=["AA at beta=%s, cos(theta)=%s, %d digits: exact %.3f s; beta^6 %.4f s."%(
        a["beta"],a["cos_theta"],a["digits"],a["exact_seconds"],a["beta6_seconds"]),
        "Relative errors M1=%.3g, M2=%.3g. Fresh single-run timings are not historical warm speedups."%(
            a["relative_errors"]["M1"],a["relative_errors"]["M2"]),
        "[Accuracy-checked AA run](%s)"%a["single_run_ci"],""]
    for name in ("HH","ZZ"):
        row=data["exact_precision"][name]
        lines+=["**%s:** 40/60/80-digit original-amplitude sweep; lowest stable %d digits;"%(
            name,row["verified_min_digits"]),
            "accuracy-preserving runtime %.2f s; high-precision discrepancy %.3g;"%(
                row["low_precision_eval_seconds"],row["top_precision_convergence_error"]),
            "archived reference discrepancy %.3g. No analytic beta^6 series."%row["archived_reference_error"],
            "[Original-amplitude CI](%s)"%row["ci"],""]
    lines+=["HH archived reference uses the same analytical evaluator; ZZ FORM reference",
            "independently checks its 36 CF helicity components, not all colour parts.","",
            "## Automatic publication","",
            "Updated from successful benchmark GitHub Actions artifacts by",
            "the validation-aware update_dashboard.py script and a workflow_run publisher.",
            "Invalid AA accuracy, missing HH/ZZ stability checks, and unsupported",
            "analytical claims cause rejection. Original amplitude evaluations remain manual.",
            "The tracked benchmark_snapshot.json is the machine-readable source.",
            ""]
    return "\n".join(lines)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--artifact-dir",type=Path,required=True)
    parser.add_argument("--run-id",type=int,required=True)
    args=parser.parse_args()
    snapshot=json.loads(SNAP.read_text())
    revised=update(snapshot,find_reports(args.artifact_dir),args.run_id)
    SNAP.write_text(json.dumps(revised,indent=2,ensure_ascii=False)+"\n")
    (DIR/"00-performance-overview.svg").write_text(draw(revised))
    (DIR/"RESULTS.md").write_text(table(revised))
    print("Verified dashboard updated",args.run_id)

if __name__=="__main__":
    main()
