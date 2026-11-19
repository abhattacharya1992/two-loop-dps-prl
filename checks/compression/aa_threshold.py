#!/usr/bin/env python3
"""Benchmark the independently published AA beta^6 series against exact data.

Requires an independent checkout of abhattacharya1992/ggAA-threshold. No
unpublished/full amplitude expressions are copied into this repository.
The default reports previously recorded timings and validation data (labelled
historical), NOT a new performance measurement. --fresh requires the exact
arXiv ancillary expressions and measures both evaluators on this machine.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib
import json
import os
import resource
from pathlib import Path
import subprocess
import sys
import time

from mpmath import mp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
AA_REPO_DEFAULT = ROOT.parent / "ggAA-threshold"


def _read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _number(number: str):
    if "/" in str(number):
        p, q = str(number).split("/", 1)
        return mp.mpf(p) / mp.mpf(q)
    return mp.mpf(number)


def complex_value(re: str, im: str):
    return mp.mpc(mp.mpf(re), mp.mpf(im))


def relative_error(exact, approximate) -> float:
    if exact == 0:
        return float(abs(approximate))
    return float(abs(exact - approximate) / abs(exact))


def historical(repo: Path) -> dict:
    timing_file = repo / "data/timing/timing_benchmark.csv"
    context_file = repo / "data/timing/timing_run.json"
    validation_file = repo / "data/validation/validation_1080.csv"
    for path in (timing_file, context_file, validation_file):
        if not path.is_file():
            raise FileNotFoundError(f"AA upstream checkout missing required file: {path}")
    context = _read_json(context_file)
    timing = []
    with timing_file.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            timing.append({
                "beta": row["Beta"], "cos_theta": row["CosTheta"],
                "working_digits": int(row["WorkingPrecision"]),
                "warm_repeats": int(row["WarmRepeats"]),
                "first_call_seconds": {
                    "exact": float(row["ExactFirstSeconds"]),
                    "beta6": float(row["SeriesFirstSeconds"]),
                },
                "warm_median_seconds": {
                    "exact": float(row["ExactWarmMedianSeconds"]),
                    "beta6": float(row["SeriesWarmMedianSeconds"]),
                },
                "reported_first_call_speedup": float(row["FirstCallSpeedup"]),
                "reported_warm_speedup": float(row["WarmSpeedup"]),
                "relative_error_M1": float(row["P1RelativeError"]),
                "relative_error_M2": float(row["P2RelativeError"]),
            })
    representations = context["representation_bytes"]
    exact_size = int(representations["unexpanded"])
    beta6_size = int(representations["beta6"])
    if exact_size <= 0 or beta6_size <= 0:
        raise ValueError("Malformed historical representation sizes")
    validation = []
    with mp.workdps(60), validation_file.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if not row["quantity"].startswith(("M1_", "M2_")):
                continue
            exact = complex_value(row["re_exact"], row["im_exact"])
            approx4 = complex_value(row["re_series4"], row["im_series4"])
            approx6 = complex_value(row["re_series6"], row["im_series6"])
            beta = _number(row["beta"])
            validation.append({
                "quantity": row["quantity"], "beta": float(beta),
                "cos_theta": row["costheta"],
                "error_beta4": relative_error(exact, approx4),
                "error_beta6": relative_error(exact, approx6),
            })
    if not timing or not validation:
        raise ValueError("No historical AA benchmark or validation records")
    checked = [v for v in validation if v["beta"] <= 0.2]
    return {
        "mode": "published_historical_measurements",
        "source_repository": "https://github.com/abhattacharya1992/ggAA-threshold",
        "input_sha256": {str(p.relative_to(repo)): _sha256(p) for p in
                         (timing_file, context_file, validation_file)},
        "environment": context.get("environment"),
        "stored_expression_bytes": {"exact": exact_size, "beta6": beta6_size},
        "expression_size_ratio_exact_over_beta6": round(exact_size / beta6_size, 4),
        "timing_rows": timing,
        "validation": {
            "total_component_rows": len(validation),
            "rows_beta_le_0_2": len(checked),
            "largest_relative_component_error_beta4": max(v["error_beta4"] for v in checked),
            "largest_relative_component_error_beta6": max(v["error_beta6"] for v in checked),
            "scope": "Stored individual flavour-component values, not all complete amplitudes",
        },
        "status": "published prior measurements, not freshly re-evaluated",
    }


def _clear_foreign_modules():
    for name in list(sys.modules):
        if name == "ggaa" or name.startswith("ggaa."):
            del sys.modules[name]


def fresh(repo: Path, *, beta: str, costheta: str, digits: int,
          repeats: int, mA: str, muR: str, Nc: int, nf: int, order: int) -> dict:
    if not (20 <= digits <= 2000 and 1 <= repeats <= 20 and order in (4, 6)):
        raise ValueError("digits, repeats or order outside supported range")
    aa_src = (repo / "src/python").resolve()
    if not (aa_src / "ggaa/twoloop.py").is_file():
        raise FileNotFoundError(f"Missing AA exact evaluator: {aa_src}")
    if not (repo / "data/series/M1_nf0.m").is_file():
        raise FileNotFoundError("Missing AA series coefficients")
    expected = [
        repo / "inputs/arXiv-1909.08993/anc/Finite_2Loop_P1",
        repo / "inputs/arXiv-1909.08993/anc/Finite_2Loop_P2",
    ]
    for path in expected:
        if not path.is_file():
            raise FileNotFoundError(
                f"Missing verified external amplitude input: {path}; "
                "run bash inputs/fetch_inputs.sh in the AA repository")
    _clear_foreign_modules()
    sys.path.insert(0, str(aa_src))
    try:
        series_mod = importlib.import_module("ggaa.series")
        exact_mod = importlib.import_module("ggaa.twoloop")
        if not (0 < _number(beta) <= mp.mpf("0.2")):
            raise ValueError("Fresh comparison limited to validated AA threshold beta <= 0.2")
        # Keep symbolic parsing/preparation distinct from repeated evaluations.
        start = time.perf_counter()
        compiled = series_mod.ThresholdSeries()
        prep_series = time.perf_counter() - start
        start = time.perf_counter()
        raw = {i: exact_mod.read_expression(expected[i - 1]) for i in (1, 2)}
        prep_exact = time.perf_counter() - start
        runs = {"exact": [], "beta_series": []}
        results = {}
        with mp.workdps(digits):
            symbols = exact_mod.symbols(_number(beta), _number(costheta), mA, muR, Nc, nf)
            # Alternate order to make the cached-special-function regime visible.
            for index in range(repeats):
                for mode in (("exact", "beta_series") if index % 2 == 0
                             else ("beta_series", "exact")):
                    start = time.perf_counter()
                    if mode == "exact":
                        evaluated = {i: exact_mod.Evaluator(symbols).evaluate(raw[i]) / 64
                                     for i in (1, 2)}
                    else:
                        evaluated = {i: compiled(i, beta, costheta, mA=mA, muR=muR,
                                                 Nc=Nc, nf=nf, order=order)
                                     for i in (1, 2)}
                    runs[mode].append(time.perf_counter() - start)
                    results[mode] = evaluated
            import statistics
            medians = {mode: statistics.median(times) for mode, times in runs.items()}
            errors = {
                f"M{i}": relative_error(results["exact"][i], results["beta_series"][i])
                for i in (1, 2)
            }
        return {
            "mode": "fresh_same_machine", "beta": beta, "cos_theta": costheta,
            "mA": mA, "muR": muR, "Nc": Nc, "nf": nf, "digits": digits,
            "order": order, "repeats": repeats,
            "prepare_seconds": {"exact_read": prep_exact, "series_parse": prep_series},
            "evaluation_seconds": runs, "evaluation_median_seconds": medians,
            "median_speedup_exact_over_series": medians["exact"] / medians["beta_series"],
            "relative_errors": errors,
            "exact_input_sha256": {p.name: _sha256(p) for p in expected},
            "peak_rss_kib_process": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            "rss_scope": "entire combined AA process; not isolated exact/series peak memory",
            "precision_note": ("Exact evaluator returns the published normalization "
                               "divided by 64; special functions use memoization "
                               "within each exact call."),
        }
    finally:
        sys.path.remove(str(aa_src))
        _clear_foreign_modules()


def markdown(report: dict) -> str:
    old = report["historical"]
    total = old["stored_expression_bytes"]
    lines = ["# AA: exact amplitude versus validated threshold beta series", "",
             "The upstream AA repository provides the only complete analytic",
             "beta series used here. HH and ZZ require independent derivations.",
             "",
             f"Historical exact expression bytes: {total['exact']:,}",
             f"Historical sixth-order expression bytes: {total['beta6']:,}",
             f"Historical expression size ratio: {old['expression_size_ratio_exact_over_beta6']:.2f}x",
             "",
             "## Published exact-versus-beta6 timings (historical; not run today)",
             "",
             "| Beta | cos theta | Exact warm median (s) | Beta6 warm median (s) | Historical speedup |",
             "|---|---:|---:|---:|---:|"]
    for t in old["timing_rows"]:
        lines.append(f"| {t['beta']} | {t['cos_theta']} | "
                     f"{t['warm_median_seconds']['exact']:.6f} | "
                     f"{t['warm_median_seconds']['beta6']:.6f} | "
                     f"{t['reported_warm_speedup']:.2f}x |")
    lines += ["", "## Validation", "",
              f"Checked {old['validation']['rows_beta_le_0_2']} stored flavour-component "
              "rows with beta <= 0.2. Largest relative deviations:",
              f"- beta4: {old['validation']['largest_relative_component_error_beta4']:.4g}",
              f"- beta6: {old['validation']['largest_relative_component_error_beta6']:.4g}",
              "",
              "Errors are for individual stored flavour components and must not",
              "be substituted for the published complete-remainder error bounds.",
              ""]
    if "fresh" in report:
        t = report["fresh"]
        lines += ["## Fresh exact-versus-series evaluation", "",
                  f"Beta: {t['beta']}; cos theta: {t['cos_theta']}; digits: {t['digits']}; "
                  f"series order: {t['order']}",
                  f"Median speedup: {t['median_speedup_exact_over_series']:.2f}x",
                  f"Relative errors: {t['relative_errors']}", ""]
    else:
        lines += ["No fresh exact-amplitude evaluation was requested or performed.", ""]
    lines += ["No HH/ZZ analytic beta-series acceleration is claimed.", ""]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aa-repo", type=Path, default=AA_REPO_DEFAULT)
    parser.add_argument("--fresh", action="store_true",
                        help="run expensive exact evaluator using upstream ancillary expressions")
    parser.add_argument("--beta", default="1/10")
    parser.add_argument("--costheta", default="0")
    parser.add_argument("--digits", type=int, default=80)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--order", type=int, choices=(4, 6), default=6)
    parser.add_argument("--mA", default="125")
    parser.add_argument("--muR", default="125")
    parser.add_argument("--Nc", type=int, default=3)
    parser.add_argument("--nf", type=int, default=5)
    parser.add_argument("--output", type=Path, default=ROOT / "compression-results/aa")
    args = parser.parse_args(argv)
    repo = args.aa_repo.resolve()
    report = {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "historical": historical(repo),
        "upstream_git_revision": subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True).stdout.strip(),
    }
    if args.fresh:
        report["fresh"] = fresh(
            repo, beta=args.beta, costheta=args.costheta, digits=args.digits,
            repeats=args.repeats, mA=args.mA, muR=args.muR, Nc=args.Nc,
            nf=args.nf, order=args.order,
        )
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    doc = markdown(report)
    (args.output / "report.md").write_text(doc)
    print(doc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
