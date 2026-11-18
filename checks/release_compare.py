#!/usr/bin/env python3
"""Compare a candidate release with the frozen publication baseline.

All calculations run in detached temporary worktrees. This command never
updates a branch and does not modify the caller's working tree.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import sys
import tarfile
import tempfile
import time

BASELINE = "1259ee1fc9dedb42047cd3919f655713dd631f27"
DATA_SUFFIXES = (".json", ".jsonl", ".pkl", ".tex", ".bib", ".bbl", ".pdf")
HISTORICAL = "paper2/r0/theory/stored/history/2026-10-22-pr1-initial-fdh/"
PAPER = Path("paper2/manuscript_prl")
SCIENTIFIC_LINE = re.compile(
    r"^(?:PASS(?:\s|:)|passed$|All checks passed\.$|"
    r"gg -> (?:HH|AA|ZZ):|kinematics:|n4:|"
    r"Independent literature comparisons:|results_stage_b\.json reproduced|"
    r"r0_.*\.json reproduced|kappa1_n4\.json reproduced|"
    r"glauber\.json reproduced|three_loop\.json reproduced|"
    r"OK: |ok\s+test_tables_reproduced|"
    r"imaginary part of the kernel:|planar N = 4:|poles,|"
    r"\{'s':|== all scripts passed)"
)
TEMP_DIR = re.compile(r"/tmp/tmp\.[A-Za-z0-9]+")


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def run(args: list[str], *, cwd: Path | None = None, stdout=None,
        check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(args, cwd=cwd, text=True, stdout=stdout,
                          stderr=subprocess.STDOUT if stdout else subprocess.PIPE,
                          check=check)


def git(repo: Path, *args: str) -> str:
    result = subprocess.run(["git", "-C", str(repo), *args], capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f"git {' '.join(args)}: {result.stderr.strip()}")
    return result.stdout.strip()


def git_blobs(repo: Path, revision: str) -> dict[str, str]:
    record = subprocess.run(["git", "-C", str(repo), "ls-tree", "-rz", revision],
                            capture_output=True, check=True).stdout
    entries = {}
    for item in record.split(b"\0"):
        if not item:
            continue
        head, path = item.split(b"\t", 1)
        mode, kind, digest = head.decode("ascii").split()
        if kind == "blob":
            entries[path.decode("utf-8")] = digest
    return entries


def compare_inputs(reference: dict[str, str], candidate: dict[str, str]):
    checked, differences = [], []
    for path, digest in sorted(reference.items()):
        if path.startswith(HISTORICAL) or not path.endswith(DATA_SUFFIXES):
            continue
        checked.append(path)
        found = candidate.get(path)
        if found != digest:
            differences.append({"path": path, "expected": digest, "found": found})
    return checked, differences


def science_lines(content: str) -> list[str]:
    lines = []
    for line in content.splitlines():
        line = re.sub(r"\x1b\[[0-9;]*m", "", line).strip()
        line = TEMP_DIR.sub("/tmp/<temporary>", line)
        if SCIENTIFIC_LINE.search(line):
            lines.append(line)
    return lines


def compare_files(first: Path, second: Path, relpaths: list[Path]) -> list[str]:
    changed = []
    for path in relpaths:
        a, b = first / path, second / path
        if not a.is_file() or not b.is_file() or sha256(a.read_bytes()) != sha256(b.read_bytes()):
            changed.append(str(path))
    return changed


def publication_sources(root: Path) -> list[Path]:
    """Paths relative to the worktree root, including all generated PRL tables."""
    files = [PAPER / item for item in (
        "numbers.tex", "sm.bib", "sm_nocite.tex", "letter_labels.tex")]
    files += [PAPER / "tables" / f.name for f in sorted((root / PAPER / "tables").glob("*.tex"))]
    return files


def compare_publication_sources(reference: Path, candidate: Path) -> dict:
    """Compare the complete generated-file inventories and their exact contents."""
    main_files = set(publication_sources(reference))
    clean_files = set(publication_sources(candidate))
    paths = sorted(main_files | clean_files)
    differences = sorted({str(path) for path in main_files ^ clean_files}
                         | set(compare_files(reference, candidate, paths)))
    return {"counts": [len(main_files), len(clean_files)],
            "differences": differences, "match": not differences}


def archive_sources(path: Path) -> dict[str, str]:
    files = {}
    with tarfile.open(path, "r:gz") as bundle:
        for member in bundle.getmembers():
            if member.isfile() and not member.name.startswith("anc/"):
                files[member.name] = sha256(bundle.extractfile(member).read())
    return files


def pdf_pages(file: Path) -> int:
    info = subprocess.run(["pdfinfo", str(file)], capture_output=True, text=True, check=True).stdout
    found = re.search(r"^Pages:\s+(\d+)\s*$", info, re.MULTILINE)
    if not found:
        raise RuntimeError(f"pdfinfo did not return a page count for {file}")
    return int(found.group(1))


def pdf_text(file: Path) -> bytes:
    return subprocess.run(["pdftotext", "-layout", str(file), "-"],
                          capture_output=True, check=True).stdout


def rendered_page(file: Path, page: int, prefix: Path) -> bytes:
    subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page), "-singlefile",
                    "-r", "108", "-gray", str(file), str(prefix)],
                   check=True, capture_output=True)
    raster = prefix.with_suffix(".pgm")
    data = raster.read_bytes()
    raster.unlink()
    return data


def compare_pdfs(reference: Path, candidate: Path) -> dict:
    output = {}
    with tempfile.TemporaryDirectory(prefix="dps-pages-") as temp:
        tmp = Path(temp)
        for name in ("letter", "supplement"):
            a = reference / PAPER / f"{name}.pdf"
            b = candidate / PAPER / f"{name}.pdf"
            na, nb = pdf_pages(a), pdf_pages(b)
            mismatches = []
            if na == nb:
                for page in range(1, na + 1):
                    if rendered_page(a, page, tmp / "main") != rendered_page(b, page, tmp / "clean"):
                        mismatches.append(page)
            output[name] = {
                "baseline_pages": na, "candidate_pages": nb,
                "text_identical": pdf_text(a) == pdf_text(b),
                "different_rendered_pages": mismatches,
            }
    return output


@contextmanager
def worktrees(repo: Path, baseline: str, candidate: str):
    with tempfile.TemporaryDirectory(prefix="dps-release-check-") as dirname:
        directory = Path(dirname)
        a, b = directory / "main", directory / "clean"
        paths = []
        try:
            for path, revision in ((a, baseline), (b, candidate)):
                git(repo, "worktree", "add", "--detach", str(path), revision)
                paths.append(path)
            yield a, b
        finally:
            for path in reversed(paths):
                git(repo, "worktree", "remove", "--force", str(path))


def measured(command: list[str], cwd: Path, logfile: Path) -> dict:
    started = time.monotonic()
    metric = logfile.with_suffix(".memory")
    timed = Path("/usr/bin/time").is_file()
    invocation = (["/usr/bin/time", "-f", "%M", "-o", str(metric)] + command
                  if timed else command)
    with logfile.open("w") as output:
        result = run(invocation, cwd=cwd, stdout=output, check=False)
    rss = None
    if timed and metric.is_file():
        try:
            rss = int(metric.read_text().strip())
        except ValueError:
            pass
        metric.unlink()
    return {"command": command, "elapsed_seconds": round(time.monotonic() - started, 3),
            "peak_rss_kb": rss, "returncode": result.returncode, "log": logfile.name}


def reproduction(root: Path, label: str, python: str, output: Path) -> list[dict]:
    # The baseline predates the combined 'reproduce' Make target.
    if label == "baseline":
        commands = [
            ["make", "test", f"PY={python}"],
            ["make", "-C", "paper2/r0/derivation", "test", f"PYTHON={python}"],
        ]
        rows = []
        for i, command in enumerate(commands, 1):
            row = measured(command, root, output / f"{label}-{i}.log")
            rows.append(row)
            if row["returncode"]:
                return rows
        # These files are checked by the preceding tests, which can rewrite
        # their formatting. The original CI also restored them before archive.
        git(root, "restore", "--worktree", "--",
            "paper2/hh/rg_argument/results/verification.json",
            "paper2/hh/stage2a/results_stage2a.json")
        command = ["make", "archive", f"PY={python}"]
        rows.append(measured(command, root, output / f"{label}-3.log"))
        return rows
    return [measured(["make", "reproduce", f"PY={python}"],
                     root, output / f"{label}-1.log")]


def benchmarks(reference: Path, candidate: Path, python: str, output: Path, repeat: int) -> dict:
    commands = {
        "shared_dps": ["make", "-C", "common/dps", "test", f"PY={python}"],
        "vector_pairs": ["make", "-C", "paper2/zz", "test", f"PY={python}"],
        "matching": ["make", "-C", "paper2/r0", "test", f"PY={python}"],
    }
    result = {}
    for name, command in commands.items():
        measurements = {"baseline": [], "candidate": []}
        for i in range(repeat):
            # Reverse the order between rounds to reduce warm-cache bias.
            order = (("baseline", reference), ("candidate", candidate))
            if i % 2:
                order = tuple(reversed(order))
            for label, root in order:
                row = measured(command, root, output / f"benchmark-{name}-{label}-{i+1}.log")
                measurements[label].append(row)
        reference_times = [x["elapsed_seconds"] for x in measurements["baseline"]]
        candidate_times = [x["elapsed_seconds"] for x in measurements["candidate"]]
        base = statistics.median(reference_times)
        clean = statistics.median(candidate_times)
        result[name] = {
            "measurements": measurements, "baseline_median_seconds": base,
            "candidate_median_seconds": clean,
            "candidate_over_baseline": round(clean / base, 4) if base else None,
            "peak_rss_kb": {
                label: max((row["peak_rss_kb"] or 0) for row in records) or None
                for label, records in measurements.items()
            },
            "failed": any(row["returncode"] for values in measurements.values() for row in values),
            "possible_regression": bool(base and clean > base * 1.25 and clean - base > 1),
        }
    return result


def markdown(result: dict) -> str:
    failed = result.get("failures", [])
    status = ("FAILED" if failed else "SOURCE PRECHECK PASSED — calculations not run"
              if result["mode"] == "quick" else "PASSED for included scientific data")
    lines = ["# DPS release comparison", "", f"**{status}**", "",
             f"- Baseline: `{result['baseline']}`",
             f"- Candidate: `{result['candidate']}`",
             f"- Tracked scientific inputs checked: {result['input_files_checked']}",
             f"- Changed required inputs: {len(result['input_mismatches'])}",
             f"- Reproduction commands passed: {result.get('reproduction_passed', False)}",
             f"- Scientific output lines matched: {result.get('scientific_log_match', False)}",
             f"- Generated TeX files matched: {result.get('publication_sources_match', False)}",
             f"- Archive top-level sources matched: {result.get('archive_sources_match', False)}",
             ""]
    for doc, values in result.get("pdfs", {}).items():
        lines.append(f"- {doc.capitalize()}: {values['baseline_pages']} / {values['candidate_pages']} pages; "
                     f"text equal: {values['text_identical']}; rendered differences: "
                     f"{values['different_rendered_pages']}")
    lines += ["", "## Measured runtimes", ""]
    for kind, values in result.get("runs", {}).items():
        lines.append(f"- {kind}: " + ", ".join(f"{v['elapsed_seconds']:.2f}s ({v['returncode']})" for v in values))
    if result.get("benchmarks"):
        lines += ["", "## Repeated benchmarks", "",
                  "Runtimes are hardware- and cache-dependent; possible regressions require investigation.", ""]
        for name, b in result["benchmarks"].items():
            lines.append(f"- {name}: baseline {b['baseline_median_seconds']:.2f}s; "
                         f"candidate {b['candidate_median_seconds']:.2f}s; "
                         f"ratio {b['candidate_over_baseline']}; "
                         f"peak RSS (KB) main/clean: {b['peak_rss_kb']['baseline']}/"
                         f"{b['peak_rss_kb']['candidate']}; "
                         f"possible regression: {b['possible_regression']}")
    lines += ["", "## Remaining release condition", "",
              "The full original third-party amplitude evaluations and complete region-integral "
              "production were **not** rerun by this comparison. Their external inputs and HPC "
              "requirements are documented in `REPRODUCING.md`. A green result here does not "
              "certify full from-scratch production or automatically authorize a merge.", ""]
    if failed:
        lines += ["## Failures", ""] + [f"- {f}" for f in failed] + [""]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", default=BASELINE)
    parser.add_argument("--candidate", default="HEAD")
    parser.add_argument("--output", default="release-check")
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--quick", action="store_true", help="compare tracked inputs without running calculations")
    parser.add_argument("--benchmark-runs", type=int, default=0, metavar="N",
                        help="perform N alternating rounds of three representative scientific tests")
    args = parser.parse_args(argv)
    if not 0 <= args.benchmark_runs <= 10:
        parser.error("--benchmark-runs must be between 0 and 10")
    root = Path(git(Path.cwd(), "rev-parse", "--show-toplevel")).resolve()
    output = (Path(args.output) if Path(args.output).is_absolute()
              else root / args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    baseline = git(root, "rev-parse", f"{args.baseline}^{{commit}}")
    candidate = git(root, "rev-parse", f"{args.candidate}^{{commit}}")
    result = {
        "baseline": baseline, "candidate": candidate,
        "mode": "quick" if args.quick else "full", "failures": [],
    }
    try:
        first, second = git_blobs(root, baseline), git_blobs(root, candidate)
        checked, mismatches = compare_inputs(first, second)
        result["input_files_checked"] = len(checked)
        result["input_mismatches"] = mismatches
        if mismatches:
            result["failures"].append(f"{len(mismatches)} tracked scientific inputs differ")
        if not args.quick and not mismatches:
            with worktrees(root, baseline, candidate) as (reference, clean):
                result["runs"] = {
                    "baseline": reproduction(reference, "baseline", args.python, output),
                    "candidate": reproduction(clean, "candidate", args.python, output),
                }
                success = all(row["returncode"] == 0 for rows in result["runs"].values() for row in rows)
                result["reproduction_passed"] = success
                if not success:
                    result["failures"].append("Reproduction command failed; inspect the recorded logs")
                else:
                    main_log = "\n".join((output / row["log"]).read_text() for row in result["runs"]["baseline"])
                    clean_log = "\n".join((output / row["log"]).read_text() for row in result["runs"]["candidate"])
                    main_science, clean_science = science_lines(main_log), science_lines(clean_log)
                    result["scientific_log_lines"] = [len(main_science), len(clean_science)]
                    result["scientific_log_match"] = bool(main_science and main_science == clean_science)
                    if not result["scientific_log_match"]:
                        result["failures"].append("Scientific log lines differ")
                    publication = compare_publication_sources(reference, clean)
                    result["publication_source_counts"] = publication["counts"]
                    result["publication_source_differences"] = publication["differences"]
                    result["publication_sources_match"] = publication["match"]
                    if not result["publication_sources_match"]:
                        result["failures"].append("Generated manuscript source files differ")
                    archive = PAPER / "release" / "PRL_DPS_arxiv_sources.tar.gz"
                    main_archive = archive_sources(reference / archive)
                    clean_archive = archive_sources(clean / archive)
                    result["archive_top_level_file_counts"] = [len(main_archive), len(clean_archive)]
                    result["archive_source_differences"] = sorted(
                        name for name in set(main_archive) | set(clean_archive)
                        if main_archive.get(name) != clean_archive.get(name))
                    result["archive_sources_match"] = not result["archive_source_differences"]
                    if not result["archive_sources_match"]:
                        result["failures"].append("ArXiv top-level source files differ")
                    result["pdfs"] = compare_pdfs(reference, clean)
                    for doc, p in result["pdfs"].items():
                        if (p["baseline_pages"] != p["candidate_pages"] or
                                not p["text_identical"] or p["different_rendered_pages"]):
                            result["failures"].append(f"{doc} rendered or text content differs")
                    if args.benchmark_runs:
                        result["benchmarks"] = benchmarks(reference, clean, args.python, output, args.benchmark_runs)
                        for name, record in result["benchmarks"].items():
                            if record["failed"]:
                                result["failures"].append(f"Benchmark command failed: {name}")
                            elif record["possible_regression"]:
                                result["failures"].append(f"Possible performance regression in {name}; review measurements")
        result["complete"] = True
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        result["complete"] = False
        result["failures"].append(f"Comparison interrupted: {error}")
    finally:
        (output / "report.json").write_text(json.dumps(result, indent=2) + "\n")
        (output / "report.md").write_text(markdown(result))
        print((output / "report.md").read_text())
    return 1 if result["failures"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
