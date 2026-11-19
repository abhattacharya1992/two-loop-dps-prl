#!/usr/bin/env python3
"""Unified two-stage HH/AA/ZZ research pipeline.

Stage 1: evaluate the available AA analytic beta-series at fixed angle,
validate its provenance against independently archived results, and measure
lossless stored-point loading for HH/AA/ZZ. Optional original full-amplitude
evaluation requires authenticated third-party inputs.

Stage 2: consume exactly the Stage 1 record, verify its SHA-256, and report
runtime, storage, expression size, memory and numerical errors with explicit
scientific scope. No unsupported HH/ZZ beta coefficients are fabricated.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib
import json
from pathlib import Path
import resource
import statistics
import subprocess
import sys
import time

import mpmath as mp

import aa_threshold
import exact_evaluators
import point_benchmark

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "compression-results/unified"
DEFAULT_AA = ROOT / "ggAA-threshold"
SCHEMA = "dps-unified-v1"
PROCESSES = ("AA", "HH", "ZZ")


def canonical_bytes(value: dict) -> bytes:
    return json.dumps(value, sort_keys=True, ensure_ascii=True,
                      allow_nan=False, separators=(",", ":")).encode("ascii")


def seal(payload: dict) -> dict:
    return {"schema": SCHEMA, "sha256": hashlib.sha256(canonical_bytes(payload)).hexdigest(),
            "payload": payload}


def verify_record(record: dict) -> dict:
    if not isinstance(record, dict) or record.get("schema") != SCHEMA:
        raise ValueError("Wrong Stage 1 schema")
    payload = record.get("payload")
    if not isinstance(payload, dict) or set(payload.get("processes", {})) != set(PROCESSES):
        raise ValueError("Stage 1 missing a scientific process")
    expected = hashlib.sha256(canonical_bytes(payload)).hexdigest()
    if record.get("sha256") != expected:
        raise ValueError("Stage 1 contents have changed since evaluation")
    for process in PROCESSES:
        row = payload["processes"][process]
        if row.get("process") != process:
            raise ValueError(f"Mismatched process record: {process}")
        if not row.get("storage", {}).get("verified_byte_identical"):
            raise ValueError(f"Unverified numerical records: {process}")
    return payload


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)
                         + "\n", encoding="utf-8")
    temporary.replace(path)


def git_head(path: Path) -> str:
    result = subprocess.run(["git", "-C", str(path), "rev-parse", "--verify", "HEAD"],
                            capture_output=True, text=True, check=True)
    return result.stdout.strip()


def _aa_series(aa_repo: Path, *, beta: str, costheta: str,
               digits: int, repeats: int) -> dict:
    """Actual evaluation of M1/M2 at O(beta^4) and O(beta^6), same process."""
    from mpmath import mp as context
    source = (aa_repo / "src/python").resolve()
    files = sorted((aa_repo / "data/series").glob("M*_nf*.m"))
    if len(files) != 6 or not (source / "ggaa/series.py").is_file():
        raise FileNotFoundError("AA threshold source and six coefficient files are required")
    if not (0 < aa_threshold._number(beta) <= context.mpf("0.2")):
        raise ValueError("AA evaluated beta outside the release's validated beta <= 0.2")
    if not -1 < aa_threshold._number(costheta) < 1:
        raise ValueError("AA scattering angle must have -1 < cos(theta) < 1")
    aa_threshold._clear_foreign_modules()
    sys.path.insert(0, str(source))
    try:
        module = importlib.import_module("ggaa.series")
        with context.workdps(digits):
            t0 = time.perf_counter()
            series = module.ThresholdSeries()
            preparation_seconds = time.perf_counter() - t0
            times = {"beta4": [], "beta6": []}
            values = {}
            # A warm-up before repeated evaluations excludes initial import and
            # compilation; compilation is measured separately above.
            for power, label in ((4, "beta4"), (6, "beta6")):
                values[label] = {
                    f"M{i}": series(i, beta, costheta, order=power)
                    for i in (1, 2)
                }
            for iteration in range(repeats):
                order = ((4, "beta4"), (6, "beta6"))
                if iteration % 2:
                    order = tuple(reversed(order))
                for power, label in order:
                    started = time.perf_counter()
                    result = {f"M{i}": series(i, beta, costheta, order=power)
                              for i in (1, 2)}
                    times[label].append(time.perf_counter() - started)
                    for value in result.values():
                        if not context.isfinite(value):
                            raise ValueError("Non-finite AA series value")
                    values[label] = result
            encoded = {
                name: {key: context.nstr(value, digits) for key, value in row.items()}
                for name, row in values.items()
            }
        medians = {name: statistics.median(measurements)
                   for name, measurements in times.items()}
        return {
            "status": "evaluated_validated_analytic_series",
            "domain": "equal-mass AA threshold beta <= 0.2, fixed cos(theta)",
            "beta": beta, "cos_theta": costheta, "digits": digits,
            "repeats": repeats,
            "orders": [4, 6], "amplitude_labels": ["M1", "M2"],
            "outputs": encoded, "preparation_seconds": preparation_seconds,
            "evaluation_seconds": times,
            "median_seconds": medians,
            "series_input_bytes": sum(p.stat().st_size for p in files),
            "source_sha256": {p.name: aa_threshold._sha256(p) for p in files},
            "memory_note": "Only whole-process peak RSS recorded; not per evaluator",
        }
    finally:
        sys.path.remove(str(source))
        aa_threshold._clear_foreign_modules()


def _process_stored(process: str, repeats: int, compression_level: int) -> dict:
    paths = point_benchmark.files_for(process)
    measured = point_benchmark.measure(process, paths, repeats, compression_level)
    first = json.loads(paths[0].read_text(encoding="utf-8"))
    return {
        "process": process,
        "storage": measured,
        "stored_points": {
            "kind": "forward_limit at fixed beta and varying delta",
            "first_beta": first["beta"],
            "first_delta_exp": first.get("delta_exp"),
            "records_checked": len(paths),
        },
    }


def stage1(args) -> dict:
    if not 1 <= args.repeats <= 10:
        raise ValueError("repeats must be 1..10")
    if not 20 <= args.digits <= 500:
        raise ValueError("digits must be 20..500")
    aa_repo = args.aa_repo.resolve()
    if not aa_repo.is_dir():
        raise FileNotFoundError(
            f"AA analytical engine checkout missing: {aa_repo}. "
            "Clone abhattacharya1992/ggAA-threshold into the PRL repository.")
    rows = {p: _process_stored(p, args.repeats, args.level)
            for p in PROCESSES}
    old = aa_threshold.historical(aa_repo)
    rows["AA"]["analytic"] = _aa_series(
        aa_repo, beta=args.beta, costheta=args.costheta,
        digits=args.digits, repeats=args.repeats)
    rows["AA"]["historical_exact_comparison"] = old
    rows["HH"]["analytic"] = {
        "status": "not_derived",
        "reason": "No validated fixed-angle class-B HH beta expansion is supplied",
    }
    rows["ZZ"]["analytic"] = {
        "status": "not_derived",
        "reason": "No validated fixed-angle VVamp ZZ beta expansion is supplied",
    }
    for p in ("HH", "ZZ"):
        rows[p]["original_exact"] = {
            "status": "not_run", "reason": "External published inputs not bundled; "
                                          "use --run-exact after installation"}
    rows["AA"]["original_exact"] = {
        "status": "historical_comparison_only", "reason":
        "Fresh exact expression evaluation is opt-in with --run-exact"}
    if args.run_exact:
        # Intentional fail-closed behavior: absence of any one published
        # full-amplitude source aborts the stage, never uses stored forward
        # points as a replacement for the exact fixed-angle amplitude.
        kinematics = exact_evaluators.fixed_angle_kinematics(
            args.beta, args.costheta, digits=args.digits)
        fresh = aa_threshold.fresh(
            aa_repo, beta=args.beta, costheta=args.costheta,
            digits=args.digits, repeats=args.repeats, mA="125", muR="125",
            Nc=3, nf=5, order=6)
        rows["AA"]["original_exact"] = {"status": "fresh_evaluated",
                                        "comparison": fresh}
        hh = exact_evaluators.exact_hh(
            kinematics, args.digits, args.repeats,
            ROOT / "paper2/hh/hh_two_loop")
        zz = exact_evaluators.exact_zz(
            kinematics, args.digits, args.repeats,
            ROOT / "paper2/zz/build/ggv",
            ROOT / "paper2/zz/scripts")
        rows["HH"]["original_exact"] = {"status": "fresh_evaluated", "result": hh}
        rows["ZZ"]["original_exact"] = {"status": "fresh_evaluated", "result": zz}
    payload = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "repository_commit": git_head(ROOT),
        "aa_repository_commit": git_head(aa_repo),
        "configuration": {
            "beta": args.beta, "cos_theta": args.costheta,
            "digits": args.digits, "repeats": args.repeats,
            "compression_level": args.level, "run_exact": args.run_exact,
        },
        "processes": rows,
        "peak_rss_kib_process": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "memory_scope": "whole Stage 1 process high-water mark, not each evaluator",
        "data_integrity": "all stored point records individually decoded byte-identically",
    }
    path = args.output / "stage1.json"
    write_json(path, seal(payload))
    return payload


def summary(payload: dict) -> dict:
    rows = {}
    for process in PROCESSES:
        r = payload["processes"][process]
        storage = r["storage"]
        row = {
            "point_records": storage["records"],
            "point_bytes": storage["raw_bytes"],
            "lossless_packed_bytes": storage["compressed_bytes"],
            "storage_percent_smaller": storage["percent_smaller"],
            "raw_json_parse_median_ms": storage["raw_json_parse_median_ms"],
            "verified_packed_parse_median_ms":
                storage["compressed_verify_and_parse_median_ms"],
            "analytical_status": r["analytic"]["status"],
            "original_exact_status": r["original_exact"]["status"],
        }
        if process == "AA":
            ana = r["analytic"]
            old = r["historical_exact_comparison"]
            reported = old["timing_rows"]
            row["fresh_series_evaluation_median_seconds"] = ana["median_seconds"]
            row["fresh_series_preparation_seconds"] = ana["preparation_seconds"]
            row["aa_series_input_bytes"] = ana["series_input_bytes"]
            row["historical_exact_vs_beta6_size_ratio"] = (
                old["expression_size_ratio_exact_over_beta6"])
            row["historical_warm_speedup_range"] = [
                min(t["reported_warm_speedup"] for t in reported),
                max(t["reported_warm_speedup"] for t in reported),
            ]
            row["historical_beta6_component_max_relative_error"] = (
                old["validation"]["largest_relative_component_error_beta6"])
            if r["original_exact"]["status"] == "fresh_evaluated":
                fresh = r["original_exact"]["comparison"]
                row["fresh_exact_vs_beta6_speedup"] = fresh[
                    "median_speedup_exact_over_series"]
                row["fresh_exact_vs_beta6_relative_errors"] = (
                    fresh["relative_errors"])
        elif r["original_exact"]["status"] == "fresh_evaluated":
            row["fresh_full_amplitude_seconds"] = r[
                "original_exact"]["result"]["evaluation_median_seconds"]
        rows[process] = row
    return {
        "schema": SCHEMA,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_commit": payload["repository_commit"],
        "aa_source_commit": payload["aa_repository_commit"],
        "configuration": payload["configuration"],
        "processes": rows,
        "memory_peak_rss_kib": payload["peak_rss_kib_process"],
        "memory_scope": payload["memory_scope"],
        "science_limits": [
            "Only AA has a validated full analytic beta-series expansion",
            "HH and ZZ stored forward-limit values are not threshold expansions",
            "Stored JSON compression is lossless data compression, not analytic compression",
            "AA historical exact-versus-series timing is not a fresh same-run measurement",
            "Only an opt-in fresh exact run yields a same-machine exact-versus-series speedup",
        ],
    }


def markdown(report: dict) -> str:
    lines = [
        "# Unified HH/AA/ZZ compression and performance", "",
        "Stage 1 evaluates the available AA beta-series and authenticates HH/AA/ZZ",
        "numerical point records. Stage 2 analyses those observations without",
        "rerunning the amplitude evaluation.", "",
        f"Candidate: {report['source_commit']}",
        f"AA source: {report['aa_source_commit']}", "",
        "## Stored point compression and loading", "",
        "| Process | Points | Raw bytes | Packed bytes | Reduction | Raw parse (ms) | Checked packed parse (ms) |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for p in PROCESSES:
        row = report["processes"][p]
        lines.append(f"| {p} | {row['point_records']} | {row['point_bytes']:,} | "
                     f"{row['lossless_packed_bytes']:,} | "
                     f"{row['storage_percent_smaller']:.2f}% | "
                     f"{row['raw_json_parse_median_ms']:.4f} | "
                     f"{row['verified_packed_parse_median_ms']:.4f} |")
    a = report["processes"]["AA"]
    lines += [
        "", "## Analytic series: AA only", "",
        "The two complete AA remainders were numerically evaluated at orders",
        "beta^4 and beta^6 by the independent published series implementation.",
        f"Fresh beta4 evaluation median: "
        f"{a['fresh_series_evaluation_median_seconds']['beta4']:.5f} s",
        f"Fresh beta6 evaluation median: "
        f"{a['fresh_series_evaluation_median_seconds']['beta6']:.5f} s",
        f"AA upstream series coefficient source bytes: "
        f"{a['aa_series_input_bytes']:,}",
        f"Historical exact-vs-beta6 expression size ratio: "
        f"{a['historical_exact_vs_beta6_size_ratio']:.2f}x",
        f"Historical warm numerical speedup range: "
        f"{a['historical_warm_speedup_range'][0]:.2f}x to "
        f"{a['historical_warm_speedup_range'][1]:.2f}x",
        "The historical results were measured in an earlier independent run,",
        "not on this runner. Component-level historical maximum relative error",
        f"for beta<=0.2: {a['historical_beta6_component_max_relative_error']:.4g}.",
        "This component statistic is not the full-remainder accuracy bound.",
    ]
    if "fresh_exact_vs_beta6_speedup" in a:
        lines += ["",
                  f"Opt-in fresh exact-versus-beta6 speedup: "
                  f"{a['fresh_exact_vs_beta6_speedup']:.3f}x",
                  f"Fresh full-remainder relative errors: "
                  f"{a['fresh_exact_vs_beta6_relative_errors']}"]
    lines += ["", "## HH and ZZ status", "",
              "The original exact HH/ZZ evaluators can be invoked with external",
              "inputs. No complete analytic fixed-angle HH or ZZ threshold series",
              "is available or asserted. Their stored points approach the forward",
              "limit at fixed beta, which is a different asymptotic regime.",
              ""]
    for p in ("HH", "ZZ"):
        r = report["processes"][p]
        lines.append(f"- {p}: analytic = {r['analytical_status']}; "
                     f"exact = {r['original_exact_status']}")
    lines += ["", "## Memory and reproducibility", "",
              f"Stage 1 whole-process peak RSS: {report['memory_peak_rss_kib']:,} KiB.",
              "This cannot be attributed individually to HH/AA/ZZ or to beta4/beta6.",
              "Stage 2 requires a SHA-256-valid Stage 1 record. Neither stage",
              "overwrites published amplitude inputs, point records or manuscripts.", ""]
    return "\n".join(lines)


def stage2(input_file: Path, output: Path) -> dict:
    record = json.loads(input_file.read_text(encoding="utf-8"))
    payload = verify_record(record)
    report = summary(payload)
    write_json(output / "stage2.json", report)
    output.mkdir(parents=True, exist_ok=True)
    text = markdown(report)
    (output / "report.md").write_text(text, encoding="utf-8")
    print(text)
    return report


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=("stage1", "stage2", "run"))
    p.add_argument("--aa-repo", type=Path, default=DEFAULT_AA)
    p.add_argument("--beta", default="1/10")
    p.add_argument("--costheta", default="0")
    p.add_argument("--digits", type=int, default=50)
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--level", type=int, default=9)
    p.add_argument("--run-exact", action="store_true",
                   help="opt in to costly HH/AA/ZZ external exact amplitudes")
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--input", type=Path, default=None,
                   help="Stage 1 JSON for standalone Stage 2")
    return p


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    if args.stage in ("stage1", "run"):
        stage1(args)
        print(f"Stage 1 record: {args.output / 'stage1.json'}", flush=True)
    if args.stage in ("stage2", "run"):
        stage2(args.input or args.output / "stage1.json", args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
