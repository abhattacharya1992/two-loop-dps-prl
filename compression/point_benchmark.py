#!/usr/bin/env python3
"""Lossless, authenticated compression and decode timing of HH/AA/ZZ point data.

This is STORAGE compression, not a threshold expansion or a replacement for
the external two-loop amplitude evaluators. Data are never modified in place.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import statistics
import struct
import time
import zlib

ROOT = Path(__file__).resolve().parents[1]
MAGIC = b"DPS-POINT-v1\0"
HEADER_LENGTH = len(MAGIC) + 8 + 32
MAX_POINT_BYTES = 16 * 1024 * 1024
POINT_FOLDERS = {
    "HH": ROOT / "paper2/hh/hh_two_loop/points",
    "AA": ROOT / "common/dps/points_aa",
    "ZZ": ROOT / "paper2/zz/points",
}
FILE_PREFIXES = {"HH": "g", "AA": "aa_g", "ZZ": "zz_b"}
AMPLITUDE_KEY = {"HH": "A", "AA": "A", "ZZ": "M"}


def compress_point(raw: bytes, level: int = 9) -> bytes:
    """Store exact bytes with a length and SHA-256, not rounded floats."""
    if not raw or len(raw) > MAX_POINT_BYTES:
        raise ValueError("Point size must be between 1 and 16777216 bytes")
    if level not in range(10):
        raise ValueError("Compression level must be 0 through 9")
    return (MAGIC + struct.pack(">Q", len(raw)) + hashlib.sha256(raw).digest()
            + zlib.compress(raw, level))


def decompress_point(data: bytes) -> bytes:
    """Reject corrupt, trailing, truncated or oversized data before JSON use."""
    if not data.startswith(MAGIC) or len(data) < HEADER_LENGTH:
        raise ValueError("Invalid point envelope")
    length = struct.unpack(">Q", data[len(MAGIC): len(MAGIC) + 8])[0]
    if not 0 < length <= MAX_POINT_BYTES:
        raise ValueError("Declared point size is invalid")
    checksum = data[len(MAGIC) + 8:HEADER_LENGTH]
    packed = data[HEADER_LENGTH:]
    try:
        obj = zlib.decompressobj()
        result = obj.decompress(packed, length + 1)
        result += obj.flush(length + 1 - len(result))
    except (zlib.error, ValueError) as error:
        raise ValueError("Invalid compressed data") from error
    if len(result) != length or not obj.eof or obj.unused_data or obj.unconsumed_tail:
        raise ValueError("Truncated, oversized or trailing compressed data")
    if hashlib.sha256(result).digest() != checksum:
        raise ValueError("Point checksum mismatch")
    return result


def files_for(process: str, folder: Path | None = None) -> list[Path]:
    path = folder or POINT_FOLDERS[process]
    if not path.is_dir():
        raise FileNotFoundError(f"Missing {process} stored-point directory: {path}")
    paths = [p for p in path.glob("*.json") if p.name.startswith(FILE_PREFIXES[process])]
    if not paths:
        raise FileNotFoundError(f"No {process} stored-point JSON inputs under {path}")
    return sorted(paths)


def read_verified_point(path: Path, process: str) -> bytes:
    raw = path.read_bytes()
    record = json.loads(raw)
    if (not isinstance(record, dict) or
            "beta" not in record or AMPLITUDE_KEY[process] not in record):
        raise ValueError(f"{path}: not a {process} amplitude-point record")
    if not isinstance(record[AMPLITUDE_KEY[process]], dict):
        raise ValueError(f"{path}: malformed amplitude mapping")
    return raw


def measure(process: str, paths: list[Path], repeats: int, level: int) -> dict:
    if not 1 <= repeats <= 100:
        raise ValueError("repeats must be 1 through 100")
    total_raw = total_compressed = 0
    cold_seconds = []
    original_parse = []
    compressed_parse = []
    historic_seconds = []
    for path in paths:
        raw = read_verified_point(path, process)
        packed = compress_point(raw, level=level)
        if decompress_point(packed) != raw:
            raise AssertionError(f"Lossless roundtrip failed for {path}")
        record = json.loads(raw)
        if isinstance(record.get("sec"), (int, float)) and record["sec"] >= 0:
            historic_seconds.append(float(record["sec"]))
        total_raw += len(raw)
        total_compressed += len(packed)
        start = time.perf_counter()
        zlib.decompress(packed[HEADER_LENGTH:])
        cold_seconds.append(time.perf_counter() - start)
        for i in range(repeats):
            # Alternate order to reduce systematic cache and order effects.
            if i % 2:
                start = time.perf_counter()
                json.loads(decompress_point(packed))
                compressed_parse.append(time.perf_counter() - start)
                start = time.perf_counter()
                json.loads(raw)
                original_parse.append(time.perf_counter() - start)
            else:
                start = time.perf_counter()
                json.loads(raw)
                original_parse.append(time.perf_counter() - start)
                start = time.perf_counter()
                json.loads(decompress_point(packed))
                compressed_parse.append(time.perf_counter() - start)
    return {
        "process": process, "records": len(paths), "raw_bytes": total_raw,
        "compressed_bytes": total_compressed,
        "bytes_saved": total_raw - total_compressed,
        "percent_smaller": round(100 * (total_raw - total_compressed) / total_raw, 3),
        "raw_json_parse_median_ms": round(statistics.median(original_parse) * 1000, 6),
        "compressed_verify_and_parse_median_ms": round(statistics.median(compressed_parse) * 1000, 6),
        "historical_full_evaluator_sec_median": (
            round(statistics.median(historic_seconds), 3) if historic_seconds else None),
        "historical_timing_count": len(historic_seconds),
        "verified_byte_identical": True,
    }


def markdown(report: dict) -> str:
    lines = [
        "# HH / AA / ZZ lossless point-data compression",
        "",
        "This reports *stored numerical-point compression*, not smaller analytic",
        "two-loop expressions or faster full-amplitude evaluation.",
        "",
        "| Process | Points | Original bytes | Compressed bytes | Smaller | JSON parse ms | Verified compressed parse ms |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for process, row in report["processes"].items():
        lines.append(f"| {process} | {row['records']} | {row['raw_bytes']:,} | "
                     f"{row['compressed_bytes']:,} | {row['percent_smaller']:.2f}% | "
                     f"{row['raw_json_parse_median_ms']:.4f} | "
                     f"{row['compressed_verify_and_parse_median_ms']:.4f} |")
    lines += [
        "", "Every compressed record was decoded and compared byte-for-byte with the",
        "original. The payloads are not committed or substituted for tracked inputs.",
        "Historical sec metadata was not measured on this machine and is not a",
        "fair cross-process timing comparison.", "",
        "HH/ZZ have no released, validated fixed-angle analytic beta series in",
        "this repository. For AA, use the separate exact-versus-beta6 benchmark",
        "with the published ggAA-threshold package and original inputs.", "",
    ]
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--process", choices=("all", "HH", "AA", "ZZ"), default="all")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--level", type=int, default=9)
    parser.add_argument("--output", type=Path, default=ROOT / "compression-results/points")
    args = parser.parse_args(argv)
    processes = ("HH", "AA", "ZZ") if args.process == "all" else (args.process,)
    results = {p: measure(p, files_for(p), args.repeats, args.level) for p in processes}
    report = {
        "schema_version": 1, "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "lossless stored amplitude-point JSON only", "algorithm": "zlib",
        "compression_level": args.level, "repeats": args.repeats, "processes": results,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    result = markdown(report)
    (args.output / "report.md").write_text(result, encoding="utf-8")
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
