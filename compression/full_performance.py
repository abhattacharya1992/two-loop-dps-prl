#!/usr/bin/env python3
"""Optional exact HH/AA/ZZ performance benchmark with two immutable stages.

Stage 1: run the published full evaluator(s) using SHA-256-pinned external
inputs. AA also runs its validated beta^6 series at the same phase-space point.
Stage 2: verify Stage 1's checksum and summarise CPU, wall time, peak RSS,
source bytes and (AA only) the numerical error and speedup.

Run each process in its OWN Python process. No fitting, invented HH/ZZ beta
series, or fallback to stored forward-point values is permitted.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import time

import aa_threshold
import exact_evaluators

ROOT = Path(__file__).resolve().parents[1]
MAGIC = "full-amplitude-benchmark-v1"
PROCESSES = ("aa", "hh", "zz")
ZZ_REQUIRED = (
    "ggvvamp-1.0.tar.gz",
    "gg_projectors-1.0.inc",
    "gg_Aj_A_qT-1.0.inc.gz",
    "vvpmasters_opt-1.0.m.gz",
)
AA_REQUIRED = ("Finite_2Loop_P1", "Finite_2Loop_P2")
HH_REQUIRED = ("Gauge1Simp.in", "Gauge2Simp.in")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _known_hashes(path: Path) -> dict[str, str]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing checksum manifest: {path}")
    result = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        pieces = line.split()
        if len(pieces) != 2 or len(pieces[0]) != 64:
            raise ValueError(f"Malformed checksum manifest line in {path}")
        name = pieces[1].removeprefix("*")
        if name in result:
            raise ValueError(f"Duplicate checksum entry: {name}")
        result[name] = pieces[0].lower()
    return result


def verified_sources(process: str, aa_repo: Path, root: Path = ROOT) -> dict:
    """Verify every required ORIGINAL input against an independent pinned hash."""
    if process == "aa":
        parent = aa_repo / "inputs/arXiv-1909.08993/anc"
        hashes = _known_hashes(aa_repo / "inputs/SHA256SUMS")
        names = AA_REQUIRED
    elif process == "hh":
        parent = root / "paper2/hh/hh_two_loop/anc"
        # This is the same fixed provenance as paper2/hh/fetch_inputs.py.
        hashes = {
            "Gauge1Simp.in": "3b2f6edede722228e8526a73b1aed1feaca6ec5bb834c5632c1995a1a444b7a3",
            "Gauge2Simp.in": "6fc90c7cc0bcb2fd8025aae8ca5279a4249519873a94168ff7373e0bc5ee7a44",
        }
        names = HH_REQUIRED
    elif process == "zz":
        parent = root / "paper2/zz/inputs"
        hashes = _known_hashes(parent / "SHA256SUMS")
        names = ZZ_REQUIRED
    else:
        raise ValueError(f"Unsupported process: {process}")
    verified = {}
    for name in names:
        if name not in hashes:
            raise ValueError(f"{process}: missing pinned SHA-256 for {name}")
        path = parent / name
        if not path.is_file():
            raise FileNotFoundError(f"{process}: missing published input: {path}")
        got = _sha256(path)
        if got != hashes[name]:
            raise ValueError(f"{process}: SHA-256 mismatch for {path}")
        verified[name] = {"sha256": got, "bytes": path.stat().st_size}
    return verified


def _git_head(path: Path) -> str:
    r = subprocess.run(
        ["git", "-C", str(path), "rev-parse", "--verify", "HEAD"],
        capture_output=True, text=True, check=True,
    )
    return r.stdout.strip()


def _canonical(value: dict) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      allow_nan=False, ensure_ascii=True).encode("ascii")


def seal(payload: dict) -> dict:
    return {"schema": MAGIC, "sha256": hashlib.sha256(_canonical(payload)).hexdigest(),
            "payload": payload}


def unseal(record: dict, expected_process: str) -> dict:
    if not isinstance(record, dict) or record.get("schema") != MAGIC:
        raise ValueError("Invalid benchmark schema")
    payload = record.get("payload")
    if not isinstance(payload, dict) or payload.get("process") != expected_process:
        raise ValueError("Wrong benchmark process")
    if record.get("sha256") != hashlib.sha256(_canonical(payload)).hexdigest():
        raise ValueError("Benchmark record has been modified")
    if payload.get("status") != "evaluated_published_full_amplitude":
        raise ValueError("Benchmark evaluation was not completed")
    if not payload.get("verified_external_inputs"):
        raise ValueError("Missing independently verified external inputs")
    return payload


def write(path: Path, value: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".temporary")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)
                    + "\n", encoding="utf-8")
    temp.replace(path)


def validated_args(process: str, beta: str, costheta: str,
                   digits: int, repeats: int) -> dict:
    if process not in PROCESSES:
        raise ValueError("process must be aa, hh, or zz")
    if not 20 <= digits <= 120:
        raise ValueError("precision must be 20..120 digits")
    if not 1 <= repeats <= 5:
        raise ValueError("repeat count must be 1..5")
    # AA's released numerical validation is limited to beta <= 0.2.
    b, c = aa_threshold._number(beta), aa_threshold._number(costheta)
    if not 0 < b <= aa_threshold._number("1/5"):
        raise ValueError("all-channel comparison requires 0 < beta <= 1/5")
    if not 0 <= c < 1:
        raise ValueError("require 0 <= cos(theta) < 1 for common HH/ZZ convention")
    return {"beta": beta, "cos_theta": costheta,
            "digits": digits, "repeats": repeats}


def stage1(process: str, args) -> dict:
    config = validated_args(process, args.beta, args.costheta,
                            args.digits, args.repeats)
    aa_repo = args.aa_repo.resolve()
    pinned = verified_sources(process, aa_repo, root=args.repo.resolve())
    kin = exact_evaluators.fixed_angle_kinematics(
        args.beta, args.costheta, digits=args.digits)
    pre_cpu = resource.getrusage(resource.RUSAGE_SELF)
    start = time.perf_counter()
    if process == "aa":
        # Both full and series coefficients are benchmarked in the same process.
        measurement = aa_threshold.fresh(
            aa_repo, beta=args.beta, costheta=args.costheta,
            digits=args.digits, repeats=args.repeats,
            mA="125", muR="125", Nc=3, nf=5, order=6,
        )
        series = sorted((aa_repo / "data/series").glob("M*_nf*.m"))
        if len(series) != 6:
            raise ValueError("AA six analytic coefficient files are missing")
        measurement["analytic_series_source_bytes"] = sum(p.stat().st_size for p in series)
        measurement["analytic_series_source_sha256"] = {p.name: _sha256(p) for p in series}
        measurement["analytic_series_source_kind"] = "published beta^6 coefficients"
    elif process == "hh":
        measurement = exact_evaluators.exact_hh(
            kin, args.digits, args.repeats,
            args.repo.resolve() / "paper2/hh/hh_two_loop",
        )
    else:
        measurement = exact_evaluators.exact_zz(
            kin, args.digits, args.repeats,
            args.repo.resolve() / "paper2/zz/build/ggv",
            args.repo.resolve() / "paper2/zz/scripts",
        )
    elapsed = time.perf_counter() - start
    post_cpu = resource.getrusage(resource.RUSAGE_SELF)
    payload = {
        "status": "evaluated_published_full_amplitude",
        "process": process, "configuration": config, "kinematics": kin,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "repo_commit": _git_head(args.repo.resolve()),
        "aa_repo_commit": _git_head(aa_repo) if process == "aa" else None,
        "verified_external_inputs": pinned,
        "measurement": measurement,
        "total_wall_seconds": elapsed,
        "cpu_user_seconds": post_cpu.ru_utime - pre_cpu.ru_utime,
        "cpu_system_seconds": post_cpu.ru_stime - pre_cpu.ru_stime,
        "whole_process_peak_rss_kib": post_cpu.ru_maxrss,
        "interpretation": ("Original complete AA Type-II two-loop coefficients with beta6 "
                           "series at a common point" if process == "aa" else
                           "Original published HH class-B or ZZ class-A subset only; "
                           "no HH/ZZ analytic beta expansion"),
    }
    write(args.output / "stage1.json", seal(payload))
    return payload


def stage2(path: Path, process: str, output: Path) -> dict:
    payload = unseal(json.loads(path.read_text(encoding="utf-8")), process)
    m = payload["measurement"]
    if process == "aa":
        median = m["evaluation_median_seconds"]
        exact_bytes = sum(v["bytes"] for v in payload["verified_external_inputs"].values())
        series_bytes = m["analytic_series_source_bytes"]
        result = {
            "exact_median_seconds": median["exact"],
            "beta6_median_seconds": median["beta_series"],
            "fresh_exact_over_beta6_speedup": m["median_speedup_exact_over_series"],
            "relative_error_M1": m["relative_errors"]["M1"],
            "relative_error_M2": m["relative_errors"]["M2"],
            "exact_published_input_bytes": exact_bytes,
            "saved_beta6_coefficient_file_bytes": series_bytes,
            "raw_input_over_coefficient_file_byte_ratio": exact_bytes / series_bytes,
            "representation_note": ("Raw third-party full files vs saved six beta-series "
                                    "coefficient files; the historic Mathematica "
                                    "expression-size comparison used different representations"),
            "preparation_seconds": m["prepare_seconds"],
            "timed_evaluation_seconds": m["evaluation_seconds"],
        }
    else:
        result = {
            "exact_median_seconds": m["evaluation_median_seconds"],
            "exact_preparation_seconds": m["preparation_seconds"],
            "input_bytes": m["input_total_bytes"],
            "timed_evaluation_seconds": m["evaluation_seconds"],
            "analytic_beta_series_status": "not_derived",
            "relative_error_status": "not_applicable_without_independent_beta_series",
            "speedup_status": "not_applicable_without_independent_beta_series",
            "comparison_note": ("Do not compare HH/ZZ exact evaluation times to AA: "
                                "different physical processes and conventions."),
        }
    result.update({
        "schema": MAGIC, "process": process,
        "source_commit": payload["repo_commit"],
        "aa_source_commit": payload["aa_repo_commit"],
        "configuration": payload["configuration"],
        "input_sha256": {name: value["sha256"] for name, value
                          in payload["verified_external_inputs"].items()},
        "whole_process_peak_rss_kib": payload["whole_process_peak_rss_kib"],
        "whole_process_total_wall_seconds": payload["total_wall_seconds"],
        "whole_process_user_cpu_seconds": payload["cpu_user_seconds"],
        "whole_process_system_cpu_seconds": payload["cpu_system_seconds"],
        "source_stage1_sha256": _sha256(path),
        "scientific_scope": payload["interpretation"],
    })
    write(output / "stage2.json", result)
    lines = [
        f"# {process.upper()} full-amplitude performance", "",
        f"Commit: {result['source_commit']}",
        f"Beta: {result['configuration']['beta']}; cos(theta): "
        f"{result['configuration']['cos_theta']}; "
        f"digits: {result['configuration']['digits']}; "
        f"repeats: {result['configuration']['repeats']}",
        f"Whole-process peak RSS: {result['whole_process_peak_rss_kib']:,} KiB",
        f"Total wall time: {result['whole_process_total_wall_seconds']:.3f} s",
        "", f"Exact evaluation median: {result['exact_median_seconds']:.4f} s",
    ]
    if process == "aa":
        lines += [
            f"Beta6 evaluation median: {result['beta6_median_seconds']:.4f} s",
            f"Fresh exact/beta6 speedup: {result['fresh_exact_over_beta6_speedup']:.3f}x",
            f"Relative error M1: {result['relative_error_M1']:.5g}",
            f"Relative error M2: {result['relative_error_M2']:.5g}",
            f"Published input bytes: {result['exact_published_input_bytes']:,}",
            f"Saved beta-series coefficient bytes: {result['saved_beta6_coefficient_file_bytes']:,}",
        ]
    else:
        lines += ["HH/ZZ analytic beta-series: NOT DERIVED.",
                  "HH/ZZ analytic speedup/relative errors: NOT AVAILABLE."]
    lines += ["", "The original published inputs were SHA-256 verified.",
              "Peak RSS is the entire process high-water mark (Linux),",
              "not an isolated per-evaluation measurement.", ""]
    output.mkdir(parents=True, exist_ok=True)
    (output / "report.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    return result


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=("run", "stage1", "stage2", "verify-inputs"))
    p.add_argument("--process", choices=PROCESSES, required=True)
    p.add_argument("--repo", type=Path, default=ROOT)
    p.add_argument("--aa-repo", type=Path, default=ROOT / "ggAA-threshold")
    p.add_argument("--beta", default="1/10")
    p.add_argument("--costheta", default="0")
    p.add_argument("--digits", type=int, default=40)
    p.add_argument("--repeats", type=int, default=1)
    p.add_argument("--output", type=Path, default=None)
    p.add_argument("--input", type=Path, default=None)
    return p


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    output = args.output or ROOT / "compression-results/full" / args.process
    if args.stage in ("run", "stage1", "verify-inputs"):
        validated_args(args.process, args.beta, args.costheta,
                       args.digits, args.repeats)
        if args.stage == "verify-inputs":
            print(json.dumps(verified_sources(args.process, args.aa_repo.resolve(),
                                             args.repo.resolve()), indent=2))
            return 0
        stage1(args.process, argparse.Namespace(
            **{**vars(args), "output": output}))
    if args.stage in ("run", "stage2"):
        stage2(args.input or output / "stage1.json", args.process, output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
