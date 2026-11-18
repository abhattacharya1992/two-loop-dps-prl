# Reproducing the DPS calculations and manuscript

The canonical development branch is `cleanup/reproducible-pipeline-20261118`.
It contains the PRL v2 manuscripts and both stages of tested cleanup.
Use `make release-check` to compare this branch with frozen PRL v2
(`d362f551357f80125649e583b22b6dbf8aeccb4e`). The original `main` commit
(`1259ee1fc9dedb42047cd3919f655713dd631f27`) is unchanged.


The default reproduction uses the evaluated amplitude points and region Laurent
coefficients included in this repository. It reruns the scientific comparisons,
checks the independent analytic derivation, rebuilds the two manuscripts, and
compiles the arXiv source package in isolation.

## Clean clone

Use Python 3.12 on Linux. The GitHub Actions workflow runs these steps on
Ubuntu with Python 3.12 and the dependencies in `requirements-pinned.txt`.

```bash
git clone https://github.com/abhattacharya1992/two-loop-dps-prl.git
cd two-loop-dps-prl
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-pinned.txt
```

For a full manuscript build, install GNU Make, Bash, Git, pdfLaTeX, BibTeX,
REVTeX 4.2, TikZ-Feynman, and the TeX fonts. On Ubuntu, the publication CI
uses the following packages:

```bash
sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  texlive-latex-base texlive-latex-recommended texlive-latex-extra \
  texlive-pictures texlive-publishers texlive-science \
  texlive-plain-generic texlive-fonts-recommended lmodern \
  poppler-utils libginac-dev
```

After installation, the complete included reproduction is:

```bash
make reproduce PY=python
```

This executes the scientific tests, the diagram and analytic matching tests,
generates the manuscript tables and figures, compiles `letter.pdf` and
`supplement.pdf`, builds the arXiv source archive, and checks that the archive
compiles outside the repository. The two PDFs are in
`paper2/manuscript_prl/`; the archive is in
`paper2/manuscript_prl/release/`.

The arXiv archive records the exact Git commit used for ancillary files. Run
this command on a committed checkout: changing tracked scientific inputs
without committing them causes the archive check to stop.

## Running parts independently

```bash
make test-science PY=python
make -C paper2/r0/derivation test PYTHON=python
make test PY=python
make build PY=python
make archive PY=python
```

`make test-science` does not need LaTeX. `make test` adds figure compilation.
`make build` rebuilds the final paper from the stored numerical data.
`make archive` builds the manuscripts and then tests the isolated source
package.

The stored-data tests and the analytic derivation are distinct checks.
The tests compare coefficients and generated numbers with recorded results;
they do not silently fetch or recompute all published two-loop amplitudes.
The HH bookkeeping and RG tests restore their original tracked output files
after comparing freshly generated values.

## Regenerating the underlying inputs

The original numerical amplitude evaluation is a separate, substantially
more expensive task and needs published input files that cannot be bundled
in this repository. The relevant entry points are:

- Scalar pairs: `make -C paper2/hh fetch`, then
  `make -C paper2/hh test-full`. Full scans are `scan`/`rescan` in
  `paper2/hh/`; see its README and Makefile.
- Vector pairs: follow `paper2/zz/README.md` to obtain the VVamp inputs.
  The directory Makefile provides `inputs`, `build`, `checks`, `scan`
  and `analyse` targets.
- Pseudoscalar pairs: `common/dps/README.md` documents the published
  amplitude inputs; `make -C common/dps aa-scan` re-evaluates the scan
  after those inputs have been installed.
- One-loop finite matching from regions:
  `make -C paper2/r0/theory quarks` and
  `make -C paper2/r0/theory gluons` recompute the stored Laurent
  coefficients. These take hours rather than seconds. Other detailed
  checks are documented in that directory's Makefile.

These routes are intentionally not part of `make reproduce`: they require
external published expressions, substantial computing time or both.
The default reproduction is complete for the publicly included numerical
records and their interpretation, not a claim to regenerate external
third-party amplitudes from first principles.

Before committing manually regenerated data, inspect `git status` and the
numerical differences. The cache checks reject partial results and mismatched
source fingerprints; do not bypass those safeguards.

## Optional region-sector performance caches

The six large `paper2/r0/theory/secdec_cache_*.pkl` files were removed from
this distribution because they contain cached sector decompositions, not
amplitude samples, Laurent coefficients, or matching results. The loader
rebuilds their entries on cache misses or source changes. Fresh integrand and
region production still requires the external inputs, precision controls and
computing time described above. The default `make reproduce` uses neither
these caches nor their archived versions.

The exact six original files remain in Git history at `d362f551357f80125649e583b22b6dbf8aeccb4e`.
To verify the Git blobs and their SHA-256 checksums:

```bash
make cache-archive-check PY=python
```

To restore optional acceleration data for the full original region runs:

```bash
python checks/region_cache_archive.py --restore
```

Restoration requires a clone containing the frozen v2 commit. Each recovered
file is checked against its pinned Git blob and SHA-256 hash. Avoid executing
untrusted pickle files from other sources.
