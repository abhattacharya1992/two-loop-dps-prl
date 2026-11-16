# Numerical tests with off-shell vector bosons of different invariant masses

**Outcome:** two independently evaluated configurations with different invariant masses of the off-shell vector bosons support the conditional predictions for the two highest
DPS logarithms. The tests do not establish a general factorisation theorem or
determine the finite matching for arbitrary virtualities. The extension also tests the one-loop normalisation and finite matching
without changing existing equal-mass amplitude records.

## Benchmark 1: (m1^2,m2^2,s) = (1,2,15/2)

The configuration is `ma2=1`, `mb2=2`, `s=15/2` (all mass-squared units).
At the forward endpoint, `x1=2/3`, `x2=1/5`, `t_plus=-1/2`, `t_minus=-4`.
Nine distances `d=sqrt(2)*10^(-k)`, k=4,...,12, were freshly evaluated with
50,...,90 decimal working digits. The k=10 point was repeated at 95 digits.
Four transverse helicities and all twenty tensor coefficients were evaluated
for each point using all four unequal-mass crossing configurations.

Predictions and the fit plan were committed **before amplitude extraction**:
`cb37bb3cb7b3a632818f9afacba348cffa1f4f20`. The evaluator implementation is
`8067769dbc689f15c00febfdedbfcaebffff669a`. The numerical input is the
published VVamp amplitude, not amplitudes synthesized from the DPS kernel.
The translation and special-function implementation are shared with the existing
ZZ implementation, so this is not an independent new implementation of VVamp.

Use `G=0.5*log(s*pT2/(ma2*mb2))+i*pi/2` and
`L=log(sqrt(ma2*mb2)/s)=-1.6683294302622921...`. Predictions in VVamp's
alpha_s/(2*pi) loop convention are

    CF: f3/c1=-4/3, f2/c1=-2*L-3,
    CA: f3/c1= 2/3, f2/c1= 2*L.

Here c1 is extracted independently from the one-loop amplitude of the same
helicity. All two-loop coefficients are fitted freely. No predicted coefficient
is imposed. The same r1=r2=0 boundary assumption as the Letter is being tested;
unequal-virtuality finite matching r0 is not assumed or derived.

## First benchmark: primary helicity ++++

The primary fit has four leading polynomial coefficients and four power-correction
coefficients (P=3). It uses k=5,...,12, compared with k=4,...,11 and an
all-nine-point overdetermined fit. Deviations include the small imaginary part
of each extracted ratio, not only differences of real parts.

| Ratio | Prediction | Extracted real part | Absolute complex deviation | Shifted-window change |
|---|---:|---:|---:|---:|
| CF G^3/c1 | -1.3333333333333333 | -1.3333333333335162 | 2.017e-13 | 1.282e-11 |
| CF G^2/c1 | 0.3366588605245842 | 0.3366588605179940 | 7.766e-12 | 4.490e-10 |
| CA G^3/c1 | 0.6666666666666667 | 0.6666666666668394 | 1.744e-13 | 1.284e-11 |
| CA G^2/c1 | -3.3366588605245842 | -3.3366588605180945 | 6.715e-12 | 4.497e-10 |

The second helicity ++-- gives maximum deviations 1.355e-13 (cubic) and
5.216e-12 (quadratic). Its parity partners are consistency checks, not extra
independent processes. Complete complex ratios, one-loop normalization shifts
and all model variations are in [stage1_summary.json](results/stage1_summary.json).

The all-nine-point least-squares result changes the cubic coefficients by
0.500 and the quadratic coefficients by 0.527 times the interpolation window
shift. Lower-order power corrections give progressively larger deviations:
for ++++, the largest quadratic deviation is 1.673e-6 without power corrections,
7.456e-7 at P=0, 2.373e-7 at P=1, 4.144e-8 at P=2, and 7.766e-12 at P=3.
This model dependence is retained rather than hidden by selecting only the best
agreement. A tiny square-interpolation residual is not itself a physics test.
Window shifts are stability diagnostics, not statistical confidence intervals
or rigorous bounds on an error shared by neighboring windows.

## Benchmark 2: (m1^2,m2^2,s) = (2,3,17)

A second independent `lhcpheno48` run used nine forward distances
`k=4,...,12` with 50,...,90 decimal digits and a further `k=10`
evaluation at 95 digits. Both `++++` and `++--` were fitted freely,
with the C_F and C_A terms separated. All momenta are expressed in
consistent mass-squared units. The same four unequal-mass crossings
and published VVamp amplitude are used; this is a distinct kinematic
test, not an independent amplitude program.

The theoretical predictions are unchanged:
`G=log(pT/sqrt(ma2*mb2/s))+i*pi/2`,
`L=log(sqrt(6)/17)`, and, with c1 the fitted one-loop G coefficient,

    CF: f3/c1=-4/3, f2/c1=-2*L-3
    CA: f3/c1= 2/3, f2/c1= 2*L

For the `++++` component, the central P=3 fits (k=5,...,12) give:

| Ratio | Prediction | Fitted real part | Absolute complex deviation | Shifted-window change |
|---|---:|---:|---:|---:|
| CF G^3/c1 | -1.3333333333333333 | -1.3333333333335049 | 1.949e-13 | 1.194e-11 |
| CF G^2/c1 | 0.8746672188843772 | 0.8746672188784078 | 7.332e-12 | 4.077e-10 |
| CA G^3/c1 | 0.6666666666666667 | 0.6666666666668518 | 1.851e-13 | 1.358e-11 |
| CA G^2/c1 | -3.8746672188843772 | -3.8746672188774572 | 6.966e-12 | 4.638e-10 |

For `++--` the maximum deviations are 1.329e-13 (cubic)
and 4.998e-12 (quadratic). Across both helicities, the
largest P=2 quadratic deviation is 2.241e-8, versus
7.332e-12 at P=3. The observed improvement with P,
and the shifted-window comparison, test stability rather
than furnishing statistically calibrated errors.

The k=10 80/95-digit control gives maximum scaled amplitude
differences of 5.943e-59 at one loop, 5.592e-56
for the C_F part and 2.609e-55 for C_A, with scaling
`max(1,abs(amplitude))`. This controls roundoff at that
point, not a systematic extrapolation uncertainty.

Both full desktop reproduction archives, with ten raw
tensor/helicity point files and two helicity-fit outputs each,
are retained in [local_reproductions/](results/local_reproductions/):
`lhcpheno48_20261115.tar.gz` (benchmark 1) and
`lhcpheno48_benchmark2_20261115.tar.gz` (benchmark 2).
Both have SHA-256 manifests. They allow the separate desktop
fits to be rerun from full amplitude data.

## Validation completed

- Thirteen offline tests passed: invariant identities, timelike polarizations,
  mass/scale symmetries, four-crossing bookkeeping, invalid-input rejection,
  cache integrity and unconstrained extraction of deliberately wrong synthetic
  coefficients. The equal-mass limit reproduces all 36 stored helicity amplitudes
  from their tensor coefficients; it does not freshly evaluate those old points.
- The fresh unequal-virtuality scan uses four distinct crossings at every point,
  with zero Li22 convergence-order failures.
- The 80/95-digit k=10 repeat has maximum scaled differences 2.542e-59 at one
  loop, 5.992e-57 for CF and 6.940e-57 for CA, across the four stored helicities.
  The scale in this comparison is max(1,|amplitude|). It controls roundoff at
  that point, not every possible kinematic configuration.
- Parity partners agree within a maximum scaled difference 5.795e-51.
- Recontracting all twenty saved tensor coefficients in all ten raw records at
  110 digits reproduced the saved helicity amplitudes within 2.021e-49 locally.
  This was a tensor-to-helicity consistency check, not a new amplitude calculation.
  The original tensor archive is not distributed; see data availability below.
- Exporting the unrounded helicities to CSV reproduces all fitted coefficients,
  window shifts and precision comparisons identically to the raw JSON analysis.

The [provenance manifest](results/stage1_provenance.json) records the original
input and source hashes, runtime, individual raw-record hashes and the exact CSV
checksum. [tensor_recontraction.json](results/tensor_recontraction.json) records
the historical local full-record check. Branch CI reruns the thirteen tests, replays the retained first-benchmark
analysis, and validates the SHA-256 checksums and member layout of both
desktop archives. CI does not recompute the expensive VVamp amplitudes.

## Reproduce and extend

From the repository root, after installing `requirements.txt` in this directory:

```bash
python paper2/unequal_virtuality/verify_stage1.py
python paper2/unequal_virtuality/analyse_scan.py \
  paper2/unequal_virtuality/results/stage1_amplitudes.csv.gz --helicity '++--'
```

To reevaluate all amplitudes or change virtualities/energy, follow
[README.md](README.md). Fresh scans save full tensor and helicity records.
The committed CSV contains every unrounded helicity amplitude needed for the
logarithmic fits and the precision comparison. The historical first-run tensor archive
`unequal_virtuality_stage1_raw_records.json.gz` was not distributed. Its
hash and historical check report remain in the manifest, but hashes alone
cannot reproduce that historical calculation. The two subsequent *desktop*
reproductions do include full tensor JSON data in the archives above. Reevaluate
with `run_scan.py` to generate full tensor records. The optional `--raw-archive`
mode verifies the original archive only, when that hash-matching file is available.

Neither third-party amplitude sources nor Python-specific compiled caches are
committed. `prepare.py` downloads and verifies the original source reproducibly.
The general code is ready for desktop scans. The two specified kinematic configurations have been tested. Other mass
ratios and energies are supported inputs, not automatically validated.
No publication claim or manuscript text changes automatically from these tests.
