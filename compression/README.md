# Compression and Performance — research branch

## Distinct scientific and storage questions

This is a self-contained HH/AA/ZZ research module inside the PRL repository. Its tests and outputs remain separate from the canonical Letter and Supplemental Material.

1. HH, AA, ZZ: point_benchmark.py implements byte-identical zlib storage
   and checks SHA-256, size and payload integrity on every decoded JSON point.
   It compares load-and-parse timings, NOT exact-amplitude evaluation.
2. AA: aa_threshold.py reads the independently published and validated
   ggAA-threshold beta^6 results, size and timing measurements. The historical
   timings are labelled as historical. An optional --fresh mode compares
   the exact and series evaluations with independently downloaded arXiv inputs.
3. HH/ZZ: no validated complete fixed-angle beta-series amplitude expressions
   are currently present. The included fixed-beta forward-limit fits and point
   records are not analytic beta series. Do not claim a numerical speedup
   from a missing beta expansion.

## Usage

Install requirements-pinned.txt, then from the PRL repository root:

    make compression-tests PY=python
    make compression-points PY=python

For the AA source data, checkout the independent project first:

    git clone https://github.com/abhattacharya1992/ggAA-threshold.git ggAA-threshold
    make compression-aa-recorded PY=python AA_REPO=ggAA-threshold

This reads the published timing_run.json, timing_benchmark.csv and
validation_1080.csv records, checking input SHA-256 and upstream Git SHA.
No additional amplitude files are tracked or copied into this repository.

For a *fresh* AA exact-vs-beta6 benchmark (requires full external inputs and
considerably longer evaluation time):

    (cd ggAA-threshold && bash inputs/fetch_inputs.sh)
    python compression/aa_threshold.py --aa-repo ggAA-threshold --fresh --beta 1/10 --costheta 0 --order 6 --digits 80 --repeats 3

The two AA amplitudes are compared in the same convention (exact divided by
64). The published historical timings also preserve the original CPU,
software and representations; do not substitute them for new measurements.

## HH and ZZ: still-open physics work

HH: derive a fixed-angle expansion of the relevant class-B full two-loop
remainders, keeping helicity and colour decomposition, before numerical
checks against the external Gauge1/2 expressions.

ZZ: derive and validate the corresponding threshold expansion of the full
massless-quark two-loop amplitudes, including all helicity and colour parts;
compare with the VVamp evaluator, not the fixed-beta small-pT DPS fits.

Both processes require independent validation over a bounded beta domain,
including withheld points, precision sensitivity, relative errors,
expression bytes, operation counts, preparation/evaluation time and memory.
This branch does not fabricate those analytic coefficients.

The optional CI only runs lossless HH/AA/ZZ and published AA record checks.
It cannot claim that original third-party full amplitudes were recomputed.
All benchmark reports are generated outputs and ignored by Git. No existing
scientific numerical records or manuscript files are deleted or modified.

## Optional original exact-amplitude evaluations at fixed angle

The original HH and ZZ evaluator backends can be timed, but only after the
external published expressions are installed:

    make -C paper2/hh fetch
    make compression-exact-hh PY=python

For ZZ, first follow paper2/zz/README.md to install VVamp expressions and
translate them into paper2/zz/build/ggv, then:

    make compression-exact-zz PY=python

These commands use s=4/(1-beta^2), equal external mass squared 1, and
t/u=1-s/2 +/- s beta cos(theta)/2. Positive beta and 0 <= cos(theta) < 1
are required by the current ZZ helicity reconstruction convention.
They report preparation time, repeated evaluation time, expression/input
bytes, normalized amplitudes, source SHA-256 and process high-water RSS.

Run HH and ZZ in separate OS processes. RSS is the whole process high-water
mark, not the incremental cost of individual function calls. Results are
stored under compression-results/exact and ignored by Git. These exact
backend calls are opt-in and excluded from CI because the external expressions
are not distributed in this repository.

Do not compare HH and ZZ exact backend timings against the AA beta6 timing:
they are different physical processes, conventions and numerical workloads.
Direct HH/ZZ exact-vs-beta benchmarking requires their own validated series.


## Unified two-stage implementation

The single entry point combines the established HH/AA/ZZ adapters into one
experiment, rather than maintaining three unrelated reporting paths.

    make compression-unified PY=python AA_REPO=ggAA-threshold

Alternatively, run each stage independently:

    make compression-stage1 PY=python AA_REPO=ggAA-threshold
    make compression-stage2 PY=python

Stage 1 writes compression-results/unified/stage1.json. It really evaluates
the independently validated AA beta^4 and beta^6 series at fixed angle. It
also authenticates and measures every stored HH/AA/ZZ forward-point record.
For HH and ZZ the default Stage 1 status is not_derived for analytic beta
series, and not_run for original full evaluators. This is intentional and
must not be interpreted as an analytical result.

Use the explicit --run-exact flag on compression/unified_pipeline.py
for the optional full-amplitude backend, only after ALL required external
HH, AA and ZZ inputs are installed. A missing input is an error, never a
fallback to stored points.

Stage 2 checks the SHA-256 sealed Stage 1 JSON file and creates
compression-results/unified/{stage2.json,report.md}. It reports
point-data compression, loading runtime, actual AA series timing,
historical AA exact/series speedups, size measurements, and the entire
Stage 1 process peak RSS. It never reruns an evaluator. Exact HH/ZZ
versus their missing analytic beta expansions cannot yet be reported.

The engine is a unified driver and interface to existing evaluators, not
a new, validated multi-process symbolic beta-series construction engine.
A process-independent rational DAG/MPL series constructor and new HH/ZZ
analytic coefficients remain separate scientific development.


## Optional full HH/AA/ZZ amplitude performance workflow

The independent .github/workflows/full-performance.yml runs by
workflow_dispatch only: it will NEVER launch expensive external evaluations
on routine commits. When the workflow file is available on the default
GitHub branch, select Actions > Full HH AA ZZ amplitude performance (manual)
> Run workflow; choose aa, hh, zz or all. Each process runs separately,
with a two-hour safety timeout and separately uploaded JSON/Markdown reports.

AA: pin ggAA-threshold to the previously verified commit; fetch its published
two-loop exact expressions and evaluate the complete M1/M2 Type-II finite
remainders together with the validated beta^6 series. Report exact/series
runtime, input sizes and relative numerical errors on the same runner.

HH: fetch both published class-B Gauge{1,2} expressions, check their SHA-256,
compile and evaluate them at fixed angle. Report setup, evaluation time, input
bytes and whole-process memory. No full HH beta expansion currently exists.

ZZ: fetch the original VVamp archive and three auxiliary files, check each
against the tracked SHA-256 manifest, translate the original expression
program and evaluate the complete massless-quark-loop class-A helicities.
Report time, input bytes and whole-process memory. No full ZZ beta expansion
currently exists.

All processes use fixed-angle equal-mass kinematics at beta <= 1/5,
0 <= cos(theta) < 1, with 20..120 working digits and 1..5 repeats.
The Stage 1 JSON record is SHA-256 sealed; Stage 2 independently verifies
the seal before reporting results. No stored forward-point data are accepted
as a substitute for complete original amplitudes.

Use locally, after acquiring all required external inputs:

    python compression/full_performance.py run --process aa --beta 1/10 --costheta 0 --digits 40 --repeats 1
    python compression/full_performance.py run --process hh --beta 1/10 --costheta 0 --digits 40 --repeats 1
    python compression/full_performance.py run --process zz --beta 1/10 --costheta 0 --digits 40 --repeats 1

Reports are generated under compression-results/full/PROCESS/ and ignored by Git.
The manual GitHub Actions workflow will not appear in the Actions picker
until its workflow file also exists on the repository's default branch.
This research workflow is distinct from the normal release-check.
