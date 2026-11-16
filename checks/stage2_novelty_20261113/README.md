# Stage 2 novelty and finite-boundary checks

These standalone checks support Stages 2.1–2.3 at baseline
`da14291dc48904a77ffb93b98cb625d704d2c607`. Requirements: Python >=3.10,
SymPy and mpmath. No production DPS module is imported.

```sh
python verify_stage2.py --output /tmp/dps-stage2-fresh.json
python test_stage2_guards.py --output /tmp/dps-stage2-guards-fresh.json
```

Output paths must be new. The first script checks the endpoint helicity relation,
known Sudakov constants and amplitude coefficients, exact rational-helicity
projections for arbitrary x, and the soft/self-energy algebra. It identifies
remaining common collinear integration targets without declaring them integrated.
The supplied numerical quark comparison uses retained data, not a fresh integral.

`inputs.json` records the source blob identifiers and labels its transcriptions.
`results.json` contains 127 named bounded checks, source/input hashes and runtime
versions. The guard requires fresh normal and optimized-Python replays to agree
with the retained record apart from runtime labels. It also rejects overwriting
an output and five deliberately incorrect physical/algebraic variants.

This does not prove complete DPS factorisation, exhaustive originality or a new
all-order NMHV law. It does not evaluate the full common collinear functions at
symbolic x, regenerate raw amplitudes or rerun the primary NMHV ancillary data.
See `../../paper2/derivations/stage2_novelty_r0_20261113.md` and
`../../docs/STAGE2_NOVELTY_STATUS_20261113.md` for the detailed conclusions.
