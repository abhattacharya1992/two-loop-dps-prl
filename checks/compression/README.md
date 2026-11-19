# Compression and Performance — research branch

## Distinct scientific and storage questions

This is an independent HH/AA/ZZ experiment, not the PRL release branch.

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
    python checks/compression/aa_threshold.py --aa-repo ggAA-threshold --fresh --beta 1/10 --costheta 0 --order 6 --digits 80 --repeats 3

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
