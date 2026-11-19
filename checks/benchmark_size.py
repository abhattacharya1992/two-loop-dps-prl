#!/usr/bin/env python3
"""Measure tracked Git content, not the compressed .git object database.

Use immutable commit hashes so a prior main branch can be compared even after
the repository's default branch has changed.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

ORIGINAL_MAIN = "1259ee1fc9dedb42047cd3919f655713dd631f27"
FROZEN_V2 = "d362f551357f80125649e583b22b6dbf8aeccb4e"

CODE_SUFFIXES = (".py", ".c", ".cc", ".cpp", ".cxx", ".h", ".hpp",
                 ".f", ".f90", ".f95", ".sh", ".wl", ".wls", ".m")
DATA_SUFFIXES = (".json", ".jsonl", ".pkl", ".csv", ".dat", ".tsv", ".npy", ".npz")
PUBLICATION_SUFFIXES = (".tex", ".bib", ".bbl", ".pdf", ".eps", ".svg", ".png", ".jpg", ".tar.gz")


def git(root: Path, *args: str) -> bytes:
    return subprocess.run(["git", "-C", str(root), *args], capture_output=True,
                          check=True).stdout


def resolve_ref(root: Path, name: str) -> str:
    return git(root, "rev-parse", "--verify", name + "^{commit}").decode().strip()


def parse_tree(data: bytes) -> dict[str, int]:
    """Parse NUL-separated git ls-tree -r -l entries, including paths with spaces."""
    inventory = {}
    for row in data.split(b"\0"):
        if not row:
            continue
        meta, name = row.split(b"\t", 1)
        mode, kind, obj, size = meta.split()
        if kind != b"blob":
            continue
        if size == b"-":
            raise ValueError("Cannot measure a non-blob entry")
        inventory[name.decode("utf-8", errors="surrogateescape")] = int(size)
    return inventory


def classify(path: str) -> str:
    lower = path.lower()
    base = lower.rsplit("/", 1)[-1]
    if lower.startswith("paper2/r0/theory/secdec_cache_") and lower.endswith(".pkl"):
        return "optional_sector_caches"
    if base == "makefile" or base == "cmakelists.txt" or lower.endswith(CODE_SUFFIXES):
        return "source_code"
    if lower.endswith(".log"):
        return "provenance_logs"
    if lower.endswith(DATA_SUFFIXES):
        return "numerical_data"
    if lower.endswith(PUBLICATION_SUFFIXES):
        return "publication_material"
    return "other"


def describe(inventory: dict[str, int]) -> dict:
    categories = {kind: {"files": 0, "bytes": 0} for kind in
                  ("source_code", "numerical_data", "publication_material",
                   "provenance_logs", "optional_sector_caches", "other")}
    for path, size in inventory.items():
        category = categories[classify(path)]
        category["files"] += 1
        category["bytes"] += size
    return {
        "tracked_files": len(inventory),
        "tracked_bytes": sum(inventory.values()),
        "categories": categories,
    }


def reduction(reference: dict, candidate: dict) -> dict:
    initial, final = reference["tracked_bytes"], candidate["tracked_bytes"]
    removed = initial - final
    return {
        "bytes_removed_net": removed,
        "percent_smaller": round(100.0 * removed / initial, 3) if initial else None,
        "tracked_file_difference": reference["tracked_files"] - candidate["tracked_files"],
        "source_code_bytes_difference": (
            reference["categories"]["source_code"]["bytes"] -
            candidate["categories"]["source_code"]["bytes"]
        ),
    }


def markdown(report: dict) -> str:
    labels = ("original_main", "frozen_v2", "candidate")
    lines = ["# Tracked-content size comparison", "",
             "Sizes are raw Git-tree blob bytes, not packed history, clone size,",
             "filesystem allocation, or mathematical-expression length.", "",
             "| Revision | Files | Tracked bytes | MiB | Source bytes | Numerical-data bytes | Cache bytes |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for label in labels:
        row = report["revisions"][label]
        items = row["categories"]
        lines.append("| " + " | ".join([
            label, str(row["tracked_files"]), str(row["tracked_bytes"]),
            f"{row['tracked_bytes'] / 1048576:.3f}",
            str(items["source_code"]["bytes"]),
            str(items["numerical_data"]["bytes"]),
            str(items["optional_sector_caches"]["bytes"])
        ]) + " |")
    lines += ["", "## Net reduction from each baseline", ""]
    for label, values in report["comparisons"].items():
        lines.append(
            f"- {label}: {values['bytes_removed_net']:,} net bytes removed; "
            f"{values['percent_smaller']}% smaller; "
            f"source code bytes reduced by {values['source_code_bytes_difference']:,} "
            "(negative means code grew)."
        )
    lines += ["", "## Interpretation", "",
              "Historical blobs are intentionally retained. Deleting a performance cache",
              "from the current tree does not necessarily reduce full clone size.",
              "The optional sector caches are distinct from retained region integrands",
              "and scientific numerical inputs.",
              "These file sizes are not a measure of amplitude-expression complexity.", ""]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", default=ORIGINAL_MAIN)
    parser.add_argument("--v2", default=FROZEN_V2)
    parser.add_argument("--candidate", default="HEAD")
    parser.add_argument("--output", default="benchmark-results/size")
    args = parser.parse_args(argv)
    root = Path(git(Path.cwd(), "rev-parse", "--show-toplevel").decode().strip())
    revisions = {}
    for label, name in (("original_main", args.original),
                        ("frozen_v2", args.v2), ("candidate", args.candidate)):
        sha = resolve_ref(root, name)
        inventory = parse_tree(git(root, "ls-tree", "-r", "-l", "-z", sha))
        revisions[label] = {"commit": sha, **describe(inventory)}
    report = {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "metric": "uncompressed tracked Git blob bytes",
        "revisions": revisions,
        "comparisons": {
            label: reduction(revisions[label], revisions["candidate"])
            for label in ("original_main", "frozen_v2")
        },
    }
    output = Path(args.output)
    if not output.is_absolute():
        output = root / output
    output.mkdir(parents=True, exist_ok=True)
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    doc = markdown(report)
    (output / "report.md").write_text(doc, encoding="utf-8")
    print(doc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
