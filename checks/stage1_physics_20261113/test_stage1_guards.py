#!/usr/bin/env python3
"""Replay and failure-injection controls for the bounded Stage-1 checker.

Run with --output PATH. Only temporary copies of verify_stage1.py are mutated;
no source, production data, or retained result is overwritten.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import mpmath as mp

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "verify_stage1.py"
REFERENCE = HERE / "results.json"


def run_copy(script: Path, output: Path, optimized: bool = False) -> subprocess.CompletedProcess:
    command = [sys.executable] + (["-O"] if optimized else [])
    return subprocess.run(command + [str(script), "--output", str(output)],
                          text=True, capture_output=True, timeout=45)


def check_mutation(source: str, old: str, new: str, function: str,
                   expected_error: str, name: str) -> dict:
    if source.count(old) != 1:
        raise ValueError(f"{name}: mutation anchor must occur once")
    namespace = {"__name__": "stage1_failure_injection", "__file__": str(SOURCE)}
    exec(compile(source.replace(old, new), f"<mutation:{name}>", "exec"), namespace)
    mp.mp.dps = 50
    try:
        namespace[function]()
    except ArithmeticError as error:
        if expected_error not in str(error):
            raise ArithmeticError(f"{name}: failed for unexpected reason: {error}") from error
        return {"name": name, "function": function, "expected_error": expected_error,
                "rejected": True}
    raise ArithmeticError(f"{name}: deliberately incorrect calculation was accepted")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output must be new")
    source = SOURCE.read_text()
    source_sha = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    retained = json.loads(REFERENCE.read_text())
    if source_sha != retained["producer_sha256"]:
        raise ArithmeticError("checker source does not match the retained result")
    with tempfile.TemporaryDirectory(prefix="dps-stage1-guards-") as temporary:
        tmp = Path(temporary)
        normal_path, optimized_path = tmp / "normal.json", tmp / "optimized.json"
        for optimized, path in ((False, normal_path), (True, optimized_path)):
            result = run_copy(SOURCE, path, optimized)
            if result.returncode or not path.is_file():
                raise ArithmeticError(f"checker failed: {result.stdout}\n{result.stderr}")
        normal = json.loads(normal_path.read_text())
        optimized_result = json.loads(optimized_path.read_text())
        if normal != optimized_result:
            raise ArithmeticError("optimized and normal Python disagree")
        if normal["check_count"] != 66 or not normal["passed"]:
            raise ArithmeticError("unexpected baseline check count or pass flag")
        # Runtime strings are provenance, not physics. Compare the remainder
        # explicitly; also report whether the complete records match bytewise.
        no_runtime = lambda d: {k: v for k, v in d.items() if k != "runtime"}
        if no_runtime(normal) != no_runtime(retained):
            raise ArithmeticError("fresh physics result differs from retained result")
        prior = normal_path.read_bytes()
        refused = run_copy(SOURCE, normal_path)
        if refused.returncode == 0 or normal_path.read_bytes() != prior:
            raise ArithmeticError("existing-output protection failed")
        replay = {"normal_passed": True, "optimized_passed": True,
                  "normal_equals_optimized": True,
                  "physics_equals_reference": True,
                  "record_equals_reference_bytewise": prior == REFERENCE.read_bytes(),
                  "existing_output_rejected_and_preserved": True}
    mutations = [
        check_mutation(source, '"PP": -J01-J10+J11}', '"PP": sp.Integer(0)}',
                       "overlap_checks", "static overlap cancels at finite eta", "omit_parent_parent"),
        check_mutation(source, "matrix = -a.T if leg in (1, 2) else a",
                       "matrix = a.T if leg in (1, 2) else a", "colour_checks",
                       "explicit-generator colour check", "wrong_antiquark_sign"),
        check_mutation(source, "k1 = sp.expand(a1+b0*t+r0)",
                       "k1 = sp.expand(a1+r0)", "evolution_checks",
                       "two-loop nu0 G^2", "omit_running_coupling_of_leading_kernel"),
        check_mutation(source, "f2 = sp.series(exponential/(1-alpha), alpha, 0, 5)",
                       "f2 = sp.series(exponential/(1+alpha), alpha, 0, 5)",
                       "evolution_checks", "two-loop nu2 G^1", "wrong_spin_two_transform"),
        check_mutation(source, "boundary = sum(r[j]*X**j for j in range(5))",
                       "boundary = sum(r[j]*X**j for j in range(4))", "evolution_checks",
                       "quadratic boundary shifts nu2 third logarithm", "omit_quartic_boundary"),
    ]
    record = {"scope": "replay, output isolation, optimization, deliberate-error detection",
              "checker_sha256": source_sha,
              "guard_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "reference_sha256": hashlib.sha256(REFERENCE.read_bytes()).hexdigest(),
              "replay": replay, "failure_injections": mutations, "passed": True}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(record, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print("PASS: replay, Python -O, output protection, and five deliberate-error controls")


if __name__ == "__main__":
    main()
