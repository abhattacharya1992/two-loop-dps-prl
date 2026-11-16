#!/usr/bin/env python3
"""Create and verify a self-contained arXiv source archive for the PRL Letter.

Run AFTER 'make' in paper2/manuscript_prl. Does not re-evaluate amplitudes.
Both letter.tex and supplement.tex are provided as independent source masters.

The Letter's data-availability sentence says that the scripts, the tables, the test suite and the
evaluated amplitude points are ancillary files of the submission. The archive therefore carries,
in anc/, the tracked public-facing scientific files from the current commit, excluding
files already packaged at the archive top level, the two generated PDFs, release/,
and internal development/CI files under revision/ and .github/.
The commit is recorded in anc/README_ANCILLARY.txt. A tracked
file outside those exclusions that differs from the commit stops the build, so that anc/ is exactly
that commit; --allow-dirty records the differing files instead. --no-anc builds the archive without
anc/.
"""
from __future__ import annotations

import argparse
import hashlib
import io
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
REQUIRED = (
    "letter.tex", "supplement.tex", "references.bib", "sm.bib",
    "letter.bbl", "supplement.bbl", "numbers.tex", "letter_labels.tex",
)
TEX_INPUT = re.compile(r"\\input\{([^}]+)\}")
GRAPHICS = re.compile(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}")
ANC = "anc"
ANC_README = "README_ANCILLARY.txt"
BUILD_OUTPUTS = frozenset({"letter.pdf", "supplement.pdf"})
RELEASE = "release"


def source_list() -> list[Path]:
    names = set(REQUIRED)
    for top in ("letter.tex", "supplement.tex"):
        src = (HERE / top).read_text()
        names.update(m.group(1) + ("" if m.group(1).endswith(".tex") else ".tex")
                     for m in TEX_INPUT.finditer(src))
        names.update(m.group(1) for m in GRAPHICS.finditer(src))
    # Includes in generated tables can themselves refer to other inputs.
    pending = list(names)
    while pending:
        name = pending.pop()
        if not name.endswith(".tex"):
            continue
        path = HERE / name
        if not path.is_file():
            raise FileNotFoundError(f"Required TeX source is missing: {name}")
        for match in TEX_INPUT.finditer(path.read_text()):
            child = match.group(1)
            child += "" if child.endswith(".tex") else ".tex"
            if child not in names:
                names.add(child)
                pending.append(child)
    files = []
    for name in sorted(names):
        path = (HERE / name).resolve()
        if not path.is_file() or not path.is_relative_to(HERE):
            raise FileNotFoundError(f"Missing or unsafe archive member: {name}")
        files.append(path)
    return files


def select_ancillary(tracked: list[str], changed: list[str], prefix: str,
                     top_level: set[str]) -> tuple[list[str], list[str]]:
    """Repository paths that go into anc/, and those of them that differ from the commit.

    prefix is the path of this folder in the repository; top_level holds the archive members
    relative to this folder. Pure function, so that the selection is testable without git.
    """
    def wanted(path: str) -> bool:
        # Exclude obsolete release helpers and historical superseded evaluator copies.
        if path == "compile_manuscipt.sh" or (path.startswith("paper2/kernel/")
            and path.endswith("_before_z6p_fix.py")):
            return False
        # Manuscript TeX files belong at archive top level, not in anc/.
        if path.endswith(".tex"):
            return False
        # Repository-only notes, internal reviews, and CI files are not arXiv ancillary data.
        if (path == 'revision' or path.startswith('revision/')
                or path.startswith('.github/') or path.startswith('derivations/')
                or path.startswith('paper2/diagrams/')):
            return False
        if not path.startswith(prefix + "/"):
            return True
        rel = path[len(prefix) + 1:]
        return not (rel in top_level or rel in BUILD_OUTPUTS or rel.startswith(RELEASE + "/"))

    paths = sorted(p for p in tracked if p and wanted(p))
    dirty = sorted(p for p in changed if p and wanted(p))
    return paths, dirty


def git(root: Path, *args: str) -> str:
    run = subprocess.run(["git", "-C", str(root), *args], capture_output=True, text=True)
    if run.returncode:
        raise RuntimeError(f"git {' '.join(args)} failed: {run.stderr.strip()}")
    return run.stdout


def ancillary(top_level: list[Path], allow_dirty: bool) -> tuple[Path, list[str], str]:
    """(repository root, paths for anc/, text of anc/README_ANCILLARY.txt)"""
    try:
        root = Path(git(HERE, "rev-parse", "--show-toplevel").strip()).resolve()
    except (OSError, RuntimeError) as error:
        raise RuntimeError("anc/ is taken from the git repository: run make_arxiv.py inside a "
                           "clone, or pass --no-anc") from error
    prefix = HERE.relative_to(root).as_posix()
    tracked = git(root, "ls-files", "-z").split("\0")
    changed = git(root, "diff", "--name-only", "-z", "HEAD").split("\0")
    members = {p.relative_to(HERE).as_posix() for p in top_level}
    paths, dirty = select_ancillary(tracked, changed, prefix, members)
    if dirty and not allow_dirty:
        listing = "\n  ".join(dirty[:20]) + ("\n  ..." if len(dirty) > 20 else "")
        raise RuntimeError("These tracked files differ from the commit; commit them first, so that "
                           f"anc/ equals the recorded commit (or pass --allow-dirty):\n  {listing}")
    commit = git(root, "rev-parse", "HEAD").strip()
    lines = [
        "Ancillary files of this arXiv submission (the Letter and its Supplemental Material).",
        "",
        f"They are the files of the repository two-loop-dps-prl at commit {commit}:",
        "the numerical inputs, scientific scripts, tables and reproducibility tests bundled",
        "with the publication; external published amplitude inputs and the original momentum",
        "region-selection scan are documented but not fully regenerated by these files.",
        "",
        "Not repeated here: the files at the top level of this arXiv source, which are the files of",
        f"{prefix}/ with the same relative paths (letter.tex, supplement.tex, the tables",
        "and figures they include, the .bib and .bbl files). To rebuild the manuscripts, copy them",
        f"into anc/{prefix}/. The PDFs are not included; arXiv builds them.",
        "",
        "Requirements: requirements-pinned.txt. Tests: run 'make test' in this directory (README.md).",
    ]
    if dirty:
        lines += ["", "Files that differed from that commit when this archive was made:"]
        lines += [f"  {p}" for p in dirty]
    return root, paths, "\n".join(lines) + "\n"


def add_member(bundle: tarfile.TarFile, name: str, data: bytes) -> None:
    entry = tarfile.TarInfo(name)
    entry.size = len(data)
    entry.mode = 0o644
    entry.mtime = 0
    bundle.addfile(entry, io.BytesIO(data))


def build_archive(dest: Path, with_anc: bool = True,
                  allow_dirty: bool = False) -> tuple[str, int, int]:
    files = source_list()
    anc = ancillary(files, allow_dirty) if with_anc else None
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(dest, "w:gz", compresslevel=9) as bundle:
        for path in files:
            add_member(bundle, path.relative_to(HERE).as_posix(), path.read_bytes())
        if anc:
            root, paths, readme = anc
            add_member(bundle, f"{ANC}/{ANC_README}", readme.encode())
            for rel in paths:
                add_member(bundle, f"{ANC}/{rel}", (root / rel).read_bytes())
    count = len(anc[1]) + 1 if anc else 0
    return hashlib.sha256(dest.read_bytes()).hexdigest(), len(files), count


def check_layout(dest: Path, with_anc: bool) -> None:
    """anc/ present (or absent) as requested, and no file both at the top level and in anc/"""
    with tarfile.open(dest, "r:gz") as bundle:
        names = bundle.getnames()
    top = {n for n in names if not n.startswith(ANC + "/")}
    anc = {n[len(ANC) + 1:] for n in names if n.startswith(ANC + "/")}
    if any(n.startswith('derivations/') for n in top | anc):
        raise RuntimeError('Repository-only derivations/ must not enter the arXiv source archive')
    if not with_anc:
        if anc:
            raise RuntimeError("anc/ found although --no-anc was given")
        return
    if ANC_README not in anc:
        raise RuntimeError(f"{ANC}/{ANC_README} is missing")
    tex_files = sorted(n for n in anc if n.endswith(".tex"))
    if tex_files:
        raise RuntimeError("TeX files are not allowed in anc/: " + ", ".join(tex_files))
    root = Path(git(HERE, "rev-parse", "--show-toplevel").strip()).resolve()
    prefix = HERE.relative_to(root).as_posix()
    repeated = sorted(t for t in top if f"{prefix}/{t}" in anc)
    if repeated:
        raise RuntimeError("Files both at the top level and in anc/: " + ", ".join(repeated))


def check_isolated(dest: Path, preview_dir: Path | None = None) -> None:
    with tempfile.TemporaryDirectory(prefix="dps-prl-arxiv-") as tmp:
        work = Path(tmp)
        with tarfile.open(dest, "r:gz") as bundle:
            if any(not (work / member.name).resolve().is_relative_to(work)
                   for member in bundle.getmembers()):
                raise RuntimeError("Unsafe archive member")
            bundle.extractall(work, filter="data")
        for master in ("letter", "supplement"):
            for _ in range(3):
                run = subprocess.run(
                    ["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                     "-file-line-error", master + ".tex"],
                    cwd=work, capture_output=True, text=True)
                if run.returncode:
                    raise RuntimeError(f"Clean {master} build failed:\n"
                                       + run.stdout[-3500:])
            log = (work / (master + ".log")).read_text(errors="replace")
            unresolved = re.findall(r"^.*(?:Citation .+ undefined|Reference .+ undefined|There were undefined references).*$", log, flags=re.M)
            if unresolved:
                raise RuntimeError(f"Unresolved citations/references in {master}:\\n" + "\\n".join(unresolved[:30]))
            if not (work / (master + ".pdf")).is_file():
                raise RuntimeError(f"No clean PDF for {master}")
            if preview_dir is not None:
                preview_dir.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(work / (master + ".pdf"), preview_dir / (master + ".pdf"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, default=HERE / RELEASE / "PRL_DPS_arxiv_sources.tar.gz")
    parser.add_argument("--check", action="store_true", help="Compile masters from clean extracted archive")
    parser.add_argument("--preview-dir", type=Path, default=None,
                        help="Save PDFs compiled from isolated archive here (requires --check)")
    parser.add_argument("--no-anc", action="store_true", help="Leave out the ancillary files (anc/)")
    parser.add_argument("--allow-dirty", action="store_true",
                        help="Build anc/ from a working tree that differs from the commit")
    args = parser.parse_args()
    if args.preview_dir is not None and not args.check:
        parser.error("--preview-dir requires --check")
    dest = args.output.resolve()
    digest, count, anc_count = build_archive(dest, not args.no_anc, args.allow_dirty)
    print(f"Archive: {dest} ({count} files at the top level, {anc_count} in {ANC}/; sha256 {digest})")
    check_layout(dest, not args.no_anc)
    if args.check:
        check_isolated(dest, args.preview_dir.resolve() if args.preview_dir else None)
        print("PASS: Letter and Supplement compile from isolated arXiv archive")


if __name__ == "__main__":
    main()
