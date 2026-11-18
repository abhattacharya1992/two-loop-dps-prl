# Comparing the cleaned release with the frozen baseline

The publication baseline is commit `1259ee1fc9dedb42047cd3919f655713dd631f27`.
The comparison is run from the cleaned branch. Neither branch is changed.

## On GitHub

1. Open **Actions** in `abhattacharya1992/two-loop-dps-prl`.
2. Select **Scientific and manuscript validation**.
3. Select **Run workflow** and choose branch `cleanup/reproducible-pipeline-20261118`.
4. Run it. The manually triggered workflow performs the frozen-baseline comparison.
5. Open the resulting run and download artifact `release-comparison`.
6. Read `report.md` first, then examine `report.json` and the individual logs for any failures.

The comparison checks every retained numerical dataset against its Git blob
identity, independently runs the published scientific tests in detached
worktrees, checks generated TeX and arXiv source files byte for byte, and
compares the Letter and Supplement as text and as rendered pages at 108 dpi.
Repeated timing benchmarks cover the shared DPS, ZZ and region-matching tests.
The benchmark measurements are not hardware-independent; investigate an
apparent regression rather than ignoring it or accepting a single noisy run.

## On a Linux workstation or HPC login node

From a full clone with Python 3.12, GNU Make, Git, the Python dependencies,
pdfLaTeX/BibTeX/REVTeX/TikZ-Feynman, and Poppler installed:

```bash
git clone https://github.com/abhattacharya1992/two-loop-dps-prl.git
cd two-loop-dps-prl
git fetch origin main cleanup/reproducible-pipeline-20261118
git switch cleanup/reproducible-pipeline-20261118
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-pinned.txt
make release-check PY=python BENCH_RUNS=3
```

This writes `release-check/report.md`, `release-check/report.json`, and the
complete command logs. For a quick source-data preflight without any numerical
runs or TeX, use `make release-check-quick`.

The comparison uses temporary Git worktrees and the same Python interpreter
and installed libraries for both versions. The baseline's older Makefile is
run using its original sequence of commands; the cleaned version uses
`make reproduce`. A working directory must contain the frozen baseline Git
commit: a shallow checkout may need `git fetch --unshallow` first.

## Remaining full-input stage

A passing `release-check` establishes equivalence for the included inputs,
analyses, reference derivations and generated publication. It does **not**
re-evaluate every original two-loop amplitude from external third-party files.
Those costly calculations must also be validated before claiming complete
from-scratch equivalence.

Use separate detached worktrees on the same HPC machine and preserve the
original inputs and code revisions:

```bash
git worktree add --detach ../dps-main 1259ee1fc9dedb42047cd3919f655713dd631f27
git worktree add --detach ../dps-clean 356e5285d6c6712ea88489d173856955082537e0
```

Install/verify the published input files **in each worktree** as documented in
`paper2/hh/README.md`, `paper2/zz/README.md` and `common/aa_source/README.md`.
For each worktree, run the same calculation commands, preserving the logs:

```bash
make -C paper2/hh fetch
make -C paper2/hh test-full
bash common/aa_source/inputs/fetch_inputs.sh
make -C common/dps aa-check
make -C paper2/zz inputs
make -C paper2/zz build
make -C paper2/zz checks
```

Regenerating entire scans and region-integral Laurent coefficients is a
separate, longer campaign. Follow the respective directory Makefiles for
`scan`, `rescan`, `aa-scan`, `quarks`, `gluons` and `n4`. Use isolated outputs and
record the third-party input SHA-256 sums, CPU configuration, precision, and
runtime in both worktrees. The produced JSON numerical records must be compared
by a precision-aware comparator and cross-checked against the paper's claims.

Several of these runs can take hours and require external amplitude sources
that are not distributed in the repository. Do not substitute existing JSON
records for a claimed fresh amplitude evaluation.

## Acceptance condition

A successful workflow does not merge code. The author decides whether the
archived historical helpers need to be present in the working branch and
whether the external amplitude-production comparisons are complete. Do not
replace `main` until these decisions and all numerical discrepancies are
resolved.
