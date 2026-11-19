#!/usr/bin/env python3
"""Benchmark retained scientific analysis targets across pinned Git revisions.

The two revisions are checked out in temporary detached worktrees. Nothing is
merged, pushed, or written into a caller's working tree by the measurements.
The tests use stored numerical records, not new external amplitude evaluations.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import platform
from pathlib import Path
import subprocess
import sys

from release_compare import benchmarks, git, worktrees

ORIGINAL_MAIN = "1259ee1fc9dedb42047cd3919f655713dd631f27"
FROZEN_V2 = "d362f551357f80125649e583b22b6dbf8aeccb4e"


def make_summary(results: dict, baseline: str) -> dict:
    """Make explicit which direction each timing ratio represents."""
    summary = {}
    for name, data in results.items():
        original = data["baseline_median_seconds"]
        current = data["candidate_median_seconds"]
        summary[name] = {
            "baseline_median_seconds": original,
            "candidate_median_seconds": current,
            "speedup_baseline_over_candidate": round(original / current, 4) if current else None,
            "candidate_percent_faster": round(100 * (original - current) / original, 3)
            if original else None,
            "peak_rss_kb": data["peak_rss_kb"],
            "failed": data["failed"],
            "possible_regression": data["possible_regression"],
        }
    return summary


def markdown(report: dict) -> str:
    lines = [
        "# DPS stored-data pipeline benchmark",
        "",
        f"Baseline type: {report['baseline_label']}",
        f"Baseline commit: {report['baseline']}",
        f"Candidate commit: {report['candidate']}",
        f"Python: {report['python']}",
        f"Environment: {report['platform']}",
        f"Rounds per workload: {report['runs']}",
        "",
        "Each paired measurement uses the same Make target, stored inputs and",
        "scientific dependency installation. Rounds alternate which revision",
        "runs first, to reduce systematic ordering effects.",
        "",
        "| Workload | Baseline median (s) | Candidate median (s) | Speedup | Candidate faster (%) | Baseline / candidate peak RSS (KiB) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for name, row in report.get("summary", {}).items():
        memory = row["peak_rss_kb"]
        ratio = row["speedup_baseline_over_candidate"]
        percent = row["candidate_percent_faster"]
        lines.append(
            f"| {name} | {row['baseline_median_seconds']:.3f} | "
            f"{row['candidate_median_seconds']:.3f} | "
            f"{ratio:.3f}x | {percent:+.3f}% | "
            f"{memory['baseline']} / {memory['candidate']} |"
        )
    lines += [
        "",
        "## Scope and caveats",
        "",
        "- The reported speedup equals baseline time divided by candidate time.",
        "- More tests may run in a newer revision, so a Make-target runtime is not",
        "  necessarily a matched microbenchmark of the same algorithm.",
        "- These tests analyse stored points, fit coefficients and validate",
        "  numerical records. They do NOT time complete HH, AA or ZZ amplitude",
        "  production from third-party formulas or beta-series evaluation.",
        "- Absolute timings and memory depend on hardware, Python, CPU load and",
        "  warm caches. Medians are indicative, not a universal performance claim.",
        "- Source and numerical validation are run separately by make reproduce;",
        "  this script never asserts mathematical equivalence from timing alone.",
        "",
    ]
    if report["errors"]:
        lines += ["## Failures", ""] + [f"- {e}" for e in report["errors"]] + [""]
    else:
        lines += ["All benchmark commands returned zero.", ""]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", default=ORIGINAL_MAIN)
    parser.add_argument("--candidate", default="HEAD")
    parser.add_argument("--baseline-label", choices=("original_main", "frozen_v2", "custom"),
                        default="original_main")
    parser.add_argument("--runs", type=int, default=5)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--output", default="benchmark-results/pipeline-original")
    args = parser.parse_args(argv)
    if not 1 <= args.runs <= 10:
        parser.error("--runs must be between 1 and 10")
    root = Path(git(Path.cwd(), "rev-parse", "--show-toplevel"))
    output = Path(args.output)
    if not output.is_absolute():
        output = root / output
    output.mkdir(parents=True, exist_ok=True)
    baseline = git(root, "rev-parse", "--verify", args.baseline + "^{commit}")
    candidate = git(root, "rev-parse", "--verify", args.candidate + "^{commit}")
    if baseline == candidate:
        parser.error("The baseline and candidate must be different commits")
    report = {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "baseline_label": args.baseline_label,
        "baseline": baseline, "candidate": candidate,
        "python": args.python,
        "platform": platform.platform(),
        "runs": args.runs,
        "raw": {}, "summary": {}, "errors": [],
    }
    try:
        with worktrees(root, baseline, candidate) as (reference, cleaned):
            report["raw"] = benchmarks(reference, cleaned, args.python, output, args.runs)
            report["summary"] = make_summary(report["raw"], args.baseline_label)
            for name, row in report["summary"].items():
                if row["failed"]:
                    report["errors"].append("Command failed in workload: " + name)
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        report["errors"].append("Benchmark interrupted: " + str(exc))
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    doc = markdown(report)
    (output / "report.md").write_text(doc, encoding="utf-8")
    print(doc)
    return 1 if report["errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
