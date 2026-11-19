#!/usr/bin/env python3
"""Accuracy-preserving performance of original HH/ZZ amplitudes.

Two *different* checks are performed:
1. At the requested fixed-angle point, recompute every available amplitude
   component at increasing decimal precision; test convergence against the
   highest-precision result (not an independent reference).
2. At a separately archived point, compare with pre-existing published-paper
   validation records. HH uses stored Bose-symmetry benchmark values from
   the original analytic gauge expressions (NOT an independent algorithm).
   ZZ uses separately obtained FORM/master-integral route-B helicity values,
   genuinely distinct from the translated route-A evaluator.

This does not construct HH/ZZ analytic beta series. The certified runtime is
the fastest working precision passing the stability gate, and is not an
analytic compression speedup.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib
import json
import math
from pathlib import Path
import statistics
import sys
import time

import mpmath as mp

import exact_evaluators as exact
import full_performance as full

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LEVELS = (40, 60, 80)
STABILITY_TOLERANCE = mp.mpf("1e-10")
REFERENCE_TOLERANCE = mp.mpf("1e-10")


def relative_difference(a, b):
    """Stable even if a component vanishes exactly by symmetry."""
    scale = max(abs(a), abs(b), mp.mpf("1e-20"))
    return abs(a-b)/scale


def parse_precisions(text: str) -> tuple[int, ...]:
    try:
        levels = tuple(int(p.strip()) for p in text.split(","))
    except (ValueError, TypeError) as error:
        raise ValueError("Precisions must be comma-separated integers") from error
    if (not 2 <= len(levels) <= 5
            or not all(20 <= n <= 120 for n in levels)
            or tuple(sorted(set(levels))) != levels):
        raise ValueError("Provide 2..5 increasing unique precisions in 20..120")
    return levels


def numeric_summary(samples: list[dict], timings: list[float],
                    precisions: tuple[int, ...], tolerance) -> dict:
    if len(samples) != len(timings) or len(samples) != len(precisions):
        raise ValueError("Precision sample and timing lengths differ")
    if any(set(row) != set(samples[-1]) for row in samples):
        raise ValueError("Precision sweeps compared different amplitude components")
    if not samples[-1]:
        raise ValueError("No numerical amplitude components evaluated")
    high = samples[-1]
    errors = []
    for row in samples:
        deviations = {k: float(relative_difference(row[k], high[k])) for k in high}
        if not all(math.isfinite(v) for v in deviations.values()):
            raise ValueError("Non-finite precision errors")
        errors.append(deviations)
    adjacent = [max(errors[i].values()) for i in range(len(errors)-1)]
    # Highest precision is a reference for convergence, not proof of truth.
    if not adjacent or adjacent[-1] > tolerance:
        raise ValueError(
            f"Highest two precisions do not converge to target {tolerance}: "
            f"observed {adjacent[-1] if adjacent else None}"
        )
    first_ok = next(
        (i for i in range(len(precisions)-1)
         if max(errors[i].values()) <= tolerance), len(precisions)-1
    )
    if not all(math.isfinite(t) and t > 0 for t in timings):
        raise ValueError("Invalid full-amplitude evaluation times")
    return {
        "precision_levels": list(precisions),
        "component_count": len(high),
        "accuracy_target": float(tolerance),
        "max_relative_difference_to_high_precision": [
            max(v.values()) for v in errors],
        "max_relative_difference_adjacent_top": adjacent[-1],
        "worst_components_at_lowest_precision": sorted(
            errors[0].items(), key=lambda v: v[1], reverse=True)[:8],
        "evaluation_seconds_by_precision": timings,
        "raw_baseline_digits": precisions[0],
        "raw_baseline_seconds": timings[0],
        "accuracy_preserving_digits": precisions[first_ok],
        "accuracy_preserving_seconds": timings[first_ok],
        "raw_over_accuracy_preserving_runtime_ratio": timings[0]/timings[first_ok],
        "precision_stability_passed": True,
    }


def _close_against_archive(calculated: dict, reference: dict, tolerance):
    if set(calculated) != set(reference):
        raise ValueError("Archived reference has a different component inventory")
    errors = {name: float(relative_difference(calculated[name], reference[name]))
              for name in reference}
    highest_error = max(errors.values())
    if not math.isfinite(highest_error) or highest_error > tolerance:
        raise ValueError(
            f"Archived reference check FAILED: {highest_error} > {tolerance}"
        )
    return {"pass": True, "component_count": len(errors),
            "max_relative_error": highest_error,
            "tolerance": float(tolerance),
            "worst_components": sorted(errors.items(),
                                       key=lambda v: v[1], reverse=True)[:8]}


def _prepare_hh(root: Path):
    scripts = root / "paper2/hh/hh_two_loop"
    with exact.prepend_python(scripts):
        ev = importlib.import_module("evaluate_hh")
        sys.setrecursionlimit(max(sys.getrecursionlimit(), 100000))
        started = time.perf_counter()
        codes = {i: ev.compile_amp(i) for i in (1, 2)}
        preparation = time.perf_counter() - started
    return ev, codes, preparation


def _hh_components(ev, codes, kin, digits: int) -> dict:
    with mp.workdps(digits):
        x = mp.mpf(kin["x"])
        y = mp.mpf(kin["y"])
        return {f"A{i}": mp.j * ev.eval_amp(codes[i], x, y, M2=1,
                                           muR2=mp.mpf("0.5"), N=3, Nf=5)/4
                for i in (1, 2)}


def measure_hh(root: Path, beta: str, costheta: str,
               precisions: tuple[int, ...], tolerance):
    ev, codes, preparation = _prepare_hh(root)
    timings, samples = [], []
    for n in precisions:
        kin = exact.fixed_angle_kinematics(beta, costheta, digits=n)
        started = time.perf_counter()
        samples.append(_hh_components(ev, codes, kin, n))
        timings.append(time.perf_counter()-started)
    convergence = numeric_summary(samples, timings, precisions, tolerance)
    # Archived Bose-symmetry input was produced at x=1/3, y=1/2.
    # It is an independent stored computation, not an independent method.
    path = root / "paper2/hh/hh_two_loop/bose_symmetry.json"
    with mp.workdps(precisions[-1]):
        raw = {f"Gauge{i}": ev.eval_amp(codes[i], mp.mpf(1)/3, mp.mpf(1)/2,
                    M2=1, muR2=mp.mpf("0.5"), N=3, Nf=5)
               for i in (1, 2)}
        record = json.loads(path.read_text(encoding="utf-8"))
        reference = {f"Gauge{i}": mp.mpmathify(record[str(i)]["at_y"])
                     for i in (1, 2)}
        archive = _close_against_archive(raw, reference, REFERENCE_TOLERANCE)
    archive.update({
        "provenance": "stored prior Bose-symmetry evaluation of same analytic expressions",
        "independent_algorithm": False,
        "reference_point": {"x": "1/3", "y": "1/2", "muR2": "1/2"},
        "reference_file": str(path.relative_to(root)),
        "reference_sha256": full._sha256(path),
        "note": "Archival consistency and Bose symmetry; not an independent amplitude implementation",
    })
    return convergence, archive, preparation


def _prepare_zz(root: Path):
    scripts = root / "paper2/zz/scripts"
    with exact.prepend_python(scripts):
        zz = importlib.import_module("zz_points")
        gg = importlib.import_module("ggvv")
        started = time.perf_counter()
        engine = gg.GGVV(root/"paper2/zz/build/ggv")
        preparation = time.perf_counter()-started
    return zz, gg, engine, preparation


def _zz_components(zz, gg, engine, kin: dict, digits: int) -> dict:
    with mp.workdps(digits):
        s, t, u = (mp.mpf(kin[k]) for k in ("s", "t", "u"))
        kin_num = {"s": s, "t": t, "u": u, "pT2": mp.mpf(kin["pT2"])}
        # The same function cache is used at t and u, as in zz_points.py.
        ff = gg.Functions()
        at = engine.compute(s, t, 1, 1, ff)
        au = engine.compute(s, u, 1, 1, ff)
        factors = {
            part: zz.all_coefficients({j: at[j][part] for j in zz.DIRECT},
                                      {j: au[j][part] for j in zz.DIRECT})
            for part in zz.PARTS
        }
        hel = zz.helicity_amplitudes(factors, kin_num)
        return {f"{label}:{part}": hel[label][part]
                for label in sorted(hel) for part in zz.PARTS}


def measure_zz(root: Path, beta: str, costheta: str,
               precisions: tuple[int, ...], tolerance):
    zz, gg, engine, preparation = _prepare_zz(root)
    timings, samples = [], []
    for n in precisions:
        kin = exact.fixed_angle_kinematics(beta, costheta, digits=n)
        started = time.perf_counter()
        samples.append(_zz_components(zz, gg, engine, kin, n))
        timings.append(time.perf_counter()-started)
    convergence = numeric_summary(samples, timings, precisions, tolerance)
    # The archived values are independently evaluated from the FORM amplitude
    # and Wolfram master integrals, as opposed to translated C++ route A.
    path = root/"paper2/zz/points_routeB/zzB_CF_b1o2_d4_p40.json"
    record = json.loads(path.read_text(encoding="utf-8"))
    with mp.workdps(precisions[-1]):
        k = zz.kinematics("1/2", 4)
        kin = {name: mp.nstr(k[name], precisions[-1])
               for name in ("s", "t", "u", "pT2")}
        values = _zz_components(zz, gg, engine, kin, precisions[-1])
        reference = {
            f"{helicity}:2l_CF": mp.mpc(parts["2l_CF"][0], parts["2l_CF"][1])
            for helicity, parts in record["M"].items()
        }
        selected = {label: values[label] for label in reference}
        archive = _close_against_archive(selected, reference, REFERENCE_TOLERANCE)
    archive.update({
        "provenance": "archived independent FORM coefficients and Wolfram master-integral route B",
        "independent_algorithm": True,
        "reference_point": {"beta": "1/2", "delta": "1e-4",
                            "scheme": "qT", "part": "2l_CF"},
        "reference_file": str(path.relative_to(root)),
        "reference_sha256": full._sha256(path),
        "note": "Checks CF helicities; other colour components tested by convergence only",
    })
    return convergence, archive, preparation


def verify_stage1(path: Path, process: str, config: dict) -> tuple[dict, str]:
    if not path.is_file():
        raise FileNotFoundError(f"Original exact benchmark is missing: {path}")
    original = full.unseal(json.loads(path.read_text(encoding="utf-8")), process)
    for key in ("beta", "cos_theta"):
        if full.aa_threshold._number(original["configuration"][key]) != full.aa_threshold._number(config[key]):
            raise ValueError(f"Precision sweep and Stage 1 disagree on {key}")
    if original["status"] != "evaluated_published_full_amplitude":
        raise ValueError("Original exact benchmark did not complete")
    return original, full._sha256(path)


def run(root: Path, process: str, beta: str, costheta: str,
        precisions: tuple[int, ...], stage1_file: Path, output: Path) -> dict:
    if process not in ("hh", "zz"):
        raise ValueError("Stability sweep applies only to HH and ZZ")
    config = full.validated_args(process, beta, costheta, precisions[0], 1)
    original, stage1_sha = verify_stage1(stage1_file, process, config)
    verified = full.verified_sources(process, root/"ggAA-threshold", root=root)
    if verified != original["verified_external_inputs"]:
        raise ValueError("Published input files changed since Stage 1")
    try:
        if process == "hh":
            stable, archive, setup = measure_hh(
                root, beta, costheta, precisions, STABILITY_TOLERANCE)
        else:
            stable, archive, setup = measure_zz(
                root, beta, costheta, precisions, STABILITY_TOLERANCE)
    except (ValueError, AssertionError) as error:
        full.write(output/"precision_failure.json", {
            "process": process, "precisions": list(precisions),
            "status": "failed", "reason": str(error),
        })
        raise
    payload = {
        "status": "precision_convergence_verified",
        "process": process,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "repo_commit": full._git_head(root),
        "original_stage1_sha256": stage1_sha,
        "configuration": config, "original_stage1_digits": original["configuration"]["digits"],
        "original_stage1_raw_seconds":
            original["measurement"]["evaluation_median_seconds"],
        "verified_external_inputs": verified,
        "setup_seconds": setup,
        "precision": stable, "archived_reference": archive,
        "interpretation": "self-convergence plus separately archived reference point; no beta expansion",
    }
    full.write(output/"precision.json", full.seal(payload))
    print(json.dumps({
        "process": process, "passed": True,
        "raw_seconds": stable["raw_baseline_seconds"],
        "accuracy_preserving_seconds": stable["accuracy_preserving_seconds"],
        "accuracy_preserving_digits": stable["accuracy_preserving_digits"],
        "archived_reference_max_error": archive["max_relative_error"],
    }, indent=2))
    return payload


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--process", choices=("hh","zz"), required=True)
    p.add_argument("--beta", default="1/10")
    p.add_argument("--costheta", default="1/3")
    p.add_argument("--precisions", default="40,60,80")
    p.add_argument("--repo", type=Path, default=ROOT)
    p.add_argument("--output", type=Path)
    args = p.parse_args(argv)
    levels = parse_precisions(args.precisions)
    out = args.output or args.repo/"compression-results/full"/args.process
    run(args.repo.resolve(), args.process, args.beta, args.costheta,
        levels, out/"stage1.json", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
