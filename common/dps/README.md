# One code for gg → HH, gg → AA and gg → ZZ

The DPS analysis of Paper 2 is the same chain for every process: evaluate the published amplitudes
on the forward line at fixed velocity, fit the powers of G = ½ ln[(m⁴ − tu)/m⁴] with power
corrections in p_T², and compare the coefficients with the factorised formula. This package does
the chain once, with one backend per process:

| module | what it does |
|---|---|
| `dps/kinematics.py` | the forward line in units of m: v = (1−β)/(1+β), y = v + δ, z, p_T², G, cos θ; the momentum fraction x = (1+β)/2 |
| `dps/fit.py` | interpolation fits in G with p_T² corrections; the reported window variation is an empirical stability diagnostic, not a statistical uncertainty |
| `dps/theory.py` | F_A, F_S, W; Δ_C, B; the two-loop predictions (eq. 5.5); the known parts of eq. (7.1) with the hard functions and scheme conversions of the three processes; the closed forms (7.2), (7.3) of r₀ |
| `dps/points.py` | the point records of the three processes in one form (loaders for `hh/hh_two_loop/points`, `zz/points` and `points_aa`) |
| `dps/colour.py` | Laurent polynomials in N and n_f: the colour decomposition of the gg → AA files in one pass |
| `dps/backends/hh.py`, `zz.py` | gg → HH and gg → ZZ: one point through the scan drivers of `hh/hh_two_loop/scan_colour_scale.py` and `zz/scripts/zz_points.py` (they need the third-party inputs of those folders) |
| `dps/backends/aa.py` | gg → AA: the one- and two-loop Type II remainders of arXiv:1909.08993, through the evaluator of the P1 package (`ggaa.twoloop`), colour term by colour term |
| `dps/analysis.py` | the leading logarithms against eq. (5.5) and r₀ against eqs. (7.2), (7.3), for gluon and quark splittings |
| `dps/report.py` | `results_hh.json`, `results_zz.json`, `results_aa.json` and a summary |
| `tests/test_dps.py` | the shared code against the stored results of `hh/`, `zz/` and `r0/` |
| `dps/provenance.py`, `tests/test_provenance.py` | actual input/source hashes, run parameters, runtime records, atomic writes and regressions for safe cache reuse |

This package is `common/dps/` of the master folder, next to `common/aa_source/` (the published AA evaluator, whose
evaluator it uses for gg → AA) and `paper2/` (whose `hh/` and `zz/` points it reads); it also runs
inside `paper2/` or inside the task folder `L_unified_code_ggAA/` (`dps/layout.py` finds the files in
each layout).

The gg → HH and gg → ZZ points are those of `hh/` and `zz/`, evaluated by the scan drivers there
(`hh/hh_two_loop/scan_colour_scale.py`, `zz/scripts/zz_points.py`); the
analysis scripts of those folders remain the record of how the amplitude tables were made, and
`make test` here shows that this package reproduces their fitted coefficients to the last digit.
The gg → AA points are new (Sec. 7.1 of the paper, Table 9) and are produced by this package.

The supplied point records remain analysis inputs. Resuming their numerical
evaluation now requires matching provenance and a complete finite payload;
untagged legacy files are not silently reused. Use a fresh output directory to
preserve those files, or explicitly request `--recompute`. See the
[reproduction guide](../../REPRODUCING.md) for commands and the distinction
between reanalysis and fresh amplitude evaluation.

## The gg → AA result (scan of 18–19 October 2026 on lhcpheno48, 80 points at β = 1/10, 1/4, 1/2, 2/3, 4/5)

Both components give the constants and logarithms of eq. (7.2) with the pseudoscalar hard function
h = C_A(−L² + 4 + ζ₂), and the helicity term follows the parity of the vertex in each component in its
own way: for ν = 0 the overlap changes, F_S → F_A (predicted before the evaluation, held to 2 × 10⁻¹⁰);
for ν = 2 the helicity-flip term changes sign, r₀^{AA} = r₀^{HH} − (4/3)(C_A − 2T_Fn_f)/W, at all five
velocities with the 1/W dependence (read off the amplitudes: the flip puts the two daughters of one
parent into the two vertices with the helicity opposite to that of the tree configuration, so its sign
relative to the tree term is M₋₋/M₊₊ of the vertex, −1 for the pseudoscalar). `dps/theory.py`
carries this sign as `Process.rho2`; `make_tables.py` asserts both differences from gg → HH at build time.

## Running

    make test          # seconds: HH/ZZ stored results; fresh AA fits and complete point inventory
    make analyse       # seconds: results_*.json for every process with points, and a summary

The gg → AA scan needs the published evaluator inputs. From `common/dps/`, prepare
the bundled evaluator and verify the downloaded expressions with:

    (cd ../aa_source && bash inputs/fetch_inputs.sh)

The evaluator is in `common/aa_source/`; an existing compatible `paper1/` tree
or a directory specified with `P1_ROOT=/path/to/AA_source` can also be used.

    make aa-check                 # about 10 min: one point at two precisions, agreeing digits
    make aa-scan WORKERS=4        # the 48 points of the scan and 2 control points, resumable
    make aa-blind WORKERS=4       # the velocities beta = 1/10, 4/5 not used to find r0

Each point is one JSON file in `points_aa/`; a run that is stopped continues where it left off.
The scan costs about 7 CPU hours (4 to 14 minutes per point, growing with the precision), the
blind velocities about 5.

## The gg → AA evaluation

The published files are polynomials in N, 1/N and n_f with kinematic coefficients. Replacing N and
n_f by the elementary Laurent polynomials of `dps/colour.py` inside the evaluator of the P1 package
gives every colour monomial in one pass, at the cost of one numerical evaluation (the arithmetic on
the polynomials is a small part of the work, which is the token stream of the 71 MB and 133 MB
files). The two-loop monomials are N², N n_f, n_f² and n_f/N, as for gg → HH.

Units m_A = μ_R = 1, for which the published one-loop files are the remainders (P1, eq. 2.18) and
the two-loop files are 64 times the remainders of P1's normalisation; the records hold M⁽²⁾ = file/64.
The files cancel strongly near the forward limit: 57–98 digits at 1 − cos θ = 10⁻¹⁰ in P1's stored
points, more at small β. The scan uses 50 + 9k digits at δ = 10⁻ᵏ, with 15 additional digits
when β < 1/2 and another 15 when β < 1/4. It repeats the point δ = 10⁻¹⁰, β = 1/2 with 15 more
digits; the regular scan itself extends to δ = 10⁻¹¹.

`make test` requires all 80 regular AA records at the five velocities above and both precision
controls. It recomputes the one-loop coefficients, two-loop logarithms, absent colour terms and
r₀ from those records. Missing records fail the test, and cached `results_aa.json` metrics are not
used to certify the two-loop data.
