# Release comparison for consolidated PRL v2 on main

The permanent working and release branch is `main`. Its scientific and
publication comparison baseline is the validated PRL v2 commit
`d362f551357f80125649e583b22b6dbf8aeccb4e`. The original, pre-v2
`main` commit, `1259ee1fc9dedb42047cd3919f655713dd631f27`, remains in
Git history. The comparison does not modify either revision or update `main`.

## GitHub Actions

Open the repository's **Actions → Scientific and manuscript validation**
workflow, click **Run workflow**, and select **main**. This manual dispatch
runs the `release_compare` job alongside normal scientific validation.
The comparison checks both revisions in temporary detached worktrees and uploads
`release-comparison`, containing `report.md`, `report.json`
and the command logs.

Normal pushes to `main` run scientific and manuscript validation; the full
frozen-v2 release comparison runs only when manually dispatched. Historical
comparisons remain available in the GitHub Actions history.

## Local release comparison

From a full clone, with Python 3.12, pinned scientific dependencies,
GNU Make, pdfLaTeX, BibTeX, REVTeX, TikZ-Feynman, Poppler and
`/usr/bin/time` available:

```bash
git clone https://github.com/abhattacharya1992/two-loop-dps-prl.git
cd two-loop-dps-prl
# The clone already checks out main.
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
before running the comparison. Their historical bytes can be recovered for
inspection or provenance using:

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
region-integral generation are outside the standard workflow. Full-amplitude
runs require published third-party inputs; the independent region generators
use the included source code. The commands for both are documented in
`REPRODUCING.md` and the relevant process-specific Makefiles. A successful
comparison demonstrates equivalence for the included published calculations,
not a guarantee that every original calculation has been rerun.

The release comparison runs in temporary detached worktrees. It does not
merge branches, modify the checked-out scientific sources, or push to `main`.
