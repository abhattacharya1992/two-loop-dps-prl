# Release comparison for the consolidated PRL v2 branch

The active working branch is `cleanup/reproducible-pipeline-20261118`.
Its scientific and publication baseline is the validated PRL v2 commit
`d362f551357f80125649e583b22b6dbf8aeccb4e`. The original `main` remains at
`1259ee1fc9dedb42047cd3919f655713dd631f27` and is never updated by these checks.

## GitHub Actions

Open the repository's **Actions → Scientific and manuscript validation**
workflow, click **Run workflow**, and select
`cleanup/reproducible-pipeline-20261118`. The comparison job runs
both revisions in temporary detached worktrees and uploads
`release-comparison`, containing `report.md`, `report.json`
and the command logs.

Normal pushes also run the scientific and manuscript CI checks.
The main/cleanup comparison from the original v1 is already recorded in the
GitHub Actions history; the present comparison concerns the updated PRL v2.

## Local release comparison

From a full clone, with Python 3.12, pinned scientific dependencies,
GNU Make, pdfLaTeX, BibTeX, REVTeX, TikZ-Feynman, Poppler and
`/usr/bin/time` available:

```bash
git clone https://github.com/abhattacharya1992/two-loop-dps-prl.git
cd two-loop-dps-prl
git switch cleanup/reproducible-pipeline-20261118
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-pinned.txt
make release-check PY=python BENCH_RUNS=3
```

The result is written to `release-check/report.md` and
`release-check/report.json`. For the source and cache-provenance
preflight only, use `make release-check-quick`.

The comparison rejects changes to required scientific inputs. Exactly six
source-fingerprinted sector-decomposition **performance** caches are excluded
from that byte-for-byte requirement. Their original Git blobs and SHA-256
checksums are recorded in `checks/region_cache_manifest.json` and verified
before running the comparison. They can be restored for expensive region
integrations using:

```bash
python checks/region_cache_archive.py --restore
```

The full release comparison regenerates numerical analyses, published
tables, Letter and Supplement sources, arXiv top-level sources and PDFs.
It compares the PDFs as extracted text and rendered pages, and measures
the runtime and peak resident memory of representative scientific tests.
An unexplained difference or performance regression fails the comparison.

## What these checks do not establish

The entire original two-loop amplitude-production chain and expensive
region-integral generation from external third-party input expressions
are outside the standard workflow. The commands for those calculations
are documented in `REPRODUCING.md` and the relevant process-specific
Makefiles. A successful comparison demonstrates equivalence for the
included published calculations, not a universal guarantee that every
external calculation has been rerun.

Nothing in this workflow merges any branch. Replacing `main` requires
the author's explicit approval.
