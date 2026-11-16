# Tests with off-shell vector bosons of different invariant masses

A separate extension of the massless-quark `gg -> V1* V2*` DPS comparison,
built on `main` **3d59dd953735bfbe2e3f8ca48dfd12b9e7b86195**.
Branch: **`unequal-virtuality`**. The original amplitude points and equal-mass scientific datasets are preserved.

Two independent tests with **different invariant masses of the off-shell vector bosons**
have been completed: `(ma2,mb2,s)=(1,2,15/2)` and `(2,3,17)`, in
consistent mass-squared units. Each includes nine forward distances and
an 80/95-digit control at k=10. The `++++` and `++--` transverse
helicity projections test the logarithms; their parity partners
are consistency checks, not separate processes.

**Both benchmarks agree with the conditional predictions** for the
leading two logarithms, separately for C_F and C_A; see
[RESULTS.md](RESULTS.md). The original first-run CSV and provenance
remain in `results/`, and full independently recomputed desktop
amplitudes for both configurations are archived under
[local_reproductions/](results/local_reproductions/).

Replay the result without downloading or reevaluating VVamp:

```bash
python paper2/unequal_virtuality/verify_stage1.py --retained-records
python paper2/unequal_virtuality/analyse_scan.py \
  paper2/unequal_virtuality/results/stage1_amplitudes.csv.gz
```

Install the two dependencies below first when working in a fresh environment.
The compressed CSV retains every unrounded helicity amplitude used in the fits;
these analysis inputs are not resumable point caches. The historical source-run full tensor archive is unavailable, so its
original recontraction check is not reproduced by the CSV alone.
However, the independent desktop reproductions of **both** benchmarks
now include the complete raw tensor and helicity JSON records.
Fresh `run_scan.py` evaluations also save all twenty coefficients.
The `--raw-archive` mode of `verify_stage1.py` checks only the original
historical archive; for desktop records, unpack the archives and rerun
`analyse_scan.py` as shown below.

Read [DERIVATION.md](DERIVATION.md) for the prediction, its boundary assumption,
normalizations, crossing map and scope. The prediction was fixed in
[predictions_stage1.json](predictions_stage1.json) before numerical extraction.
The massless-quark class-A qT finite remainder at `mu2=s` is used. Electroweak
couplings, top-mass effects and rates are not included. The finite matching
comparison is tested for the two stated invariant-mass configurations.
The `--retained-records` option verifies the historical point fingerprints,
precision controls, predictions and fits against their stored records; it
preserves their original source hashes rather than equating them to the
current evaluator. The original evaluator source would be required for a
strict source-byte replay.

## Reproduce both desktop validations

From the repository root, verify the recorded checksums:

```bash
sha256sum -c paper2/unequal_virtuality/results/local_reproductions/lhcpheno48_20261115.sha256
sha256sum -c paper2/unequal_virtuality/results/local_reproductions/lhcpheno48_benchmark2_20261115.sha256
```

To reanalyse the second configuration without expensive amplitude evaluation:

```bash
mkdir -p paper2/unequal_virtuality/runs/replay_benchmark2
tar -xzf paper2/unequal_virtuality/results/local_reproductions/lhcpheno48_benchmark2_20261115.tar.gz \
  -C paper2/unequal_virtuality/runs/replay_benchmark2
python paper2/unequal_virtuality/analyse_scan.py \
  paper2/unequal_virtuality/runs/replay_benchmark2/a2_b3_s17 \
  --output paper2/unequal_virtuality/runs/replay_benchmark2/fits_pppp.json
python paper2/unequal_virtuality/analyse_scan.py \
  paper2/unequal_virtuality/runs/replay_benchmark2/a2_b3_s17 \
  --helicity '++--' \
  --output paper2/unequal_virtuality/runs/replay_benchmark2/fits_ppmm.json
```

For the first benchmark, extract `lhcpheno48_20261115.tar.gz`
and analyse its `stage1/` directory similarly. Each archive contains
ten points (nine distances plus the high-precision repeat) and two
fit-result JSON files. These are **research validation data**, separate
from the PRL arXiv ancillary package. GitHub CI checks the archive
checksums and contents but does not rerun the full VVamp calculation.

## Install and validate

From the repository root, using Python 3.10 or later:

```bash
python3 -m venv .venv-uv
source .venv-uv/bin/activate
python -m pip install -r paper2/unequal_virtuality/requirements.txt
make -C paper2/unequal_virtuality test PY=python
python paper2/unequal_virtuality/prepare.py
```

The preparation command downloads the published VVamp archive from CERN's
source mirror, verifies its SHA-256, unpacks it safely and runs the existing
translator. Generated input and build files stay in this directory and are
ignored by Git. A previous verified archive can be used offline:

```bash
python paper2/unequal_virtuality/prepare.py --archive /path/to/ggvvamp-1.0.tar.gz
```

Only `mpmath` is needed for amplitude evaluation and fitting. `sympy` is used
by the independent algebra/unit tests. The Python-version-specific
`algebra.marshal` must be regenerated after changing Python versions; do not
copy it from a different interpreter. The code uses the existing evaluator
unchanged and makes four crossed evaluations when the virtualities differ.

## A single point, the benchmark scan, and precision control

```bash
# A single initial point
python paper2/unequal_virtuality/run_scan.py --k 4

# Nine distinct distances; processes, not threads. Start with one or two.
python paper2/unequal_virtuality/run_scan.py --workers 2

# Repeat k=10 at 95 rather than 80 digits in the same output directory.
python paper2/unequal_virtuality/run_scan.py --k 10 --dps 95

# Fit the primary helicity and write coefficients, residuals and window shifts.
python paper2/unequal_virtuality/analyse_scan.py \
  paper2/unequal_virtuality/runs/stage1 \
  --output paper2/unequal_virtuality/runs/stage1_analysis.json

# The second helicity is already stored; no new amplitude evaluation is needed.
python paper2/unequal_virtuality/analyse_scan.py \
  paper2/unequal_virtuality/runs/stage1 --helicity '++--' \
  --output paper2/unequal_virtuality/runs/stage1_analysis_ppmm.json
```

Point `k` means `t=t_plus-sqrt(ma2*mb2)*10^-k`, not `pT2=10^-k`.
Default working precision is `30+5*k` decimal digits. `--dps` overrides it;
`--dps-base` and `--dps-step` configure a schedule. Fitting uses higher arithmetic
precision, which does not increase the accuracy of stored amplitude inputs.

Each point is saved atomically as soon as it finishes. Resume only reuses a
file when its inputs, code, parameters and runtime agree exactly. A partial,
legacy or incompatible file is rejected rather than silently relabelled.
Preserve such records and use a fresh output directory, or explicitly request
`--recompute`. Failed recomputations do not overwrite an existing valid file.
No command modifies the original equal-mass points or the manuscripts.

## General configurations

The input flags are **squared** virtualities, not masses. For example:

```bash
python paper2/unequal_virtuality/run_scan.py \
  --ma2 2 --mb2 3 --s 17 --k 4 5 6 7 8 9 10 11 12 \
  --workers 2 --output paper2/unequal_virtuality/runs/a2_b3_s17
python paper2/unequal_virtuality/analyse_scan.py \
  paper2/unequal_virtuality/runs/a2_b3_s17 \
  --output paper2/unequal_virtuality/runs/a2_b3_s17_analysis.json
```

This is the **second validated configuration**, independently rerun on
`lhcpheno48`; its fitted results and precision control are in
[RESULTS.md](RESULTS.md). Supply rational or decimal strings. All three invariants
must use the same units. Require positive timelike virtualities and energy above
threshold; the DPS expansion is intended away from threshold and fraction
endpoints. Very extreme scale hierarchies need their own numerical study.

`--all-helicities` stores all 36 projections, while `--helicities` selects them.
This does not save much evaluation time: the expensive tensor coefficients
are shared. `analyse_scan.py` deliberately only tests the isotropic logarithmic
channels `++++`, `++--`, `----`, `--++`. A longitudinal or angular-constant
channel requires its appropriate theoretical prediction, not division by a
vanishing c1. A run directory must contain one physical configuration only.

## Interpretation

All coefficients are fitted freely. The predictions are used only afterward.
The output compares no power correction and correction polynomials P=0,1,2,3,
central and shifted windows, an extra-point least-squares fit, and independent
one-loop normalization shifts. Nine distinct points are recommended: P=3 needs
eight parameters and a ninth point for the shifted window.

A small interpolation residual alone is not evidence of the theory: square
interpolation has no remaining degrees of freedom. Inspect coefficient stability
and the overdetermined comparison. Window shifts are not statistical or rigorous
error bounds. A second-precision point controls roundoff at that point, not every
point or every possible mass ratio. Two successful configurations with different virtualities
do not prove general factorisation.

The complete input identity is pinned to
`1c69a6beb3be85b07225af867586f242e999337a1ccf1fb940016a5341dfc091`
for `ggvvamp-1.0.tar.gz`. The primary reference is arXiv:1503.08835; see
[DERIVATION.md](DERIVATION.md) and `../zz/README.md` for attribution and known
shared ingredients of the amplitude checks.
