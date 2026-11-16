# The DPS singularity of two-loop gluon-fusion amplitudes: code and data

This package reproduces the numbers, tables and figures of

> A. Bhattacharya, *The double parton scattering singularity of two-loop gluon-fusion
> amplitudes*, arXiv:XXXX.XXXXX.

It complements the package of the gg → AA threshold paper,
https://github.com/abhattacharya1992/ggAA-threshold.

## Quick start

    pip install -r requirements.txt
    make test          # about a minute, no downloads
    make fetch         # the two-loop gg -> HH expressions of arXiv:1809.05388, checked by SHA-256
    make test-full     # also evaluates those expressions directly (2-3 minutes more)

`make all` reruns every analysis and check, rebuilds the tables and the figures, and tests.

## Contents

| directory | paper | what it does |
|---|---|---|
| `one_loop/` | Sec. 2 | the one-loop DPS term from the Feynman rules (`verify_dps_one_loop.py`), the opposite-helicity numerator for CP-even vertices (`hh_spin2_numerator.py`), the Landau determinants of the one-loop topologies (`landau_boxes.py`), and the exact one-loop gg → AA amplitudes of P1 from 90° to the forward limit (`aa_full_range.py`, Fig. 2; `make aa`, which needs the P1 package with its inputs) |
| `hh_two_loop/` | Secs. 2.5, 3 | evaluation of the two-loop gg → HH remainders, the colour and scale scan, its 50 evaluated points, and the analyses: colour structure (`analyse_colour.py`), the one-loop amplitude from the μ_R dependence (`analyse_one_loop.py`), scale dependence (`analyse_scale.py`), the forward limit of the opposite-helicity amplitude (`forward_approach.py`, Sec. 2.5), the dependence of the fits on the window and the model (`analyse_windows.py`, Fig. 5), the choice of the model by the window shift, which does not use the predictions, and least-squares fits through more points (`check_fit_models.py`, App. B), and Li₂,₂ against an independent representation for every argument of the scan and of the blind test (`check_li22_values.py`, App. A; `make li22` after `make fetch`) |
| `stage2a/` | Sec. 4 | what UV renormalisation and the infrared subtraction contribute to the leading logarithms |
| `rg_argument/` | Sec. 5 | the renormalisation-group structure and its boundary term |
| `orientation/` | Sec. 2.4 | the two quark-loop orientations of the massless crossed box |
| `../manuscript/` | the paper: Tables 1–8, Figs. 1–5, every number quoted in the text | `make_tables.py` writes `numbers.tex` and the tables from the result files here, in `../zz` and in `../n4` (and from `p1_numbers.json` for the few numbers quoted from P1); `make_figure.py` draws Figs. 2–5 from the same files; `make figures` here runs both (Fig. 1 is a TikZ drawing, compiled by `make` in `../manuscript`) |
| `tests/` | | `test_release.py` |

The check scripts in `one_loop/`, `stage2a/`, `rg_argument/` and `orientation/` stop with
an error if a check fails. The analyses in `hh_two_loop/` print their comparisons, and
`tests/test_release.py` turns those into pass or fail. The logs of the runs used for the
paper are in the `logs/` directories.

## Third-party input

The two-loop amplitudes are those of P. Banerjee, S. Borowka, P. K. Dhani, T. Gehrmann and
V. Ravindran, JHEP 11 (2018) 130 [arXiv:1809.05388], ancillary files `Gauge1Simp.in` and
`Gauge2Simp.in`. The active input directory is `hh_two_loop/anc/`; `fetch_inputs.py` downloads
them there and checks

    Gauge1Simp.in  3b2f6edede722228e8526a73b1aed1feaca6ec5bb834c5632c1995a1a444b7a3
    Gauge2Simp.in  6fc90c7cc0bcb2fd8025aae8ca5279a4249519873a94168ff7373e0bc5ee7a44

The files are not redistributed here. From this directory:

    make fetch
    make test-full

The full test checks the input hashes before evaluating them. Its separate AA curve check needs
the P1 one-loop inputs and is reported as skipped when those inputs are absent.

They are read without Mathematica. The only function beyond classical polylogarithms,
Li₂,₂(a, b), is evaluated from its integral representation, −a∫₀¹dt ln t Li₂(abt)/(1 − at), with
mpmath. Towards the forward limit the branch point t = 1/(ab) of Li₂ approaches the end of the path
from beyond: arguments with ab = 1 − ε, b > 1 and ε of the order of δ occur. mpmath's tanh-sinh
quadrature then reports convergence before it resolves the branch point. The integral is therefore
split at points 1 − d, d = 1/2, 1/8, 1/32, …, that accumulate at t = 1 until d is no larger than the
distance of the nearest singular point, and evaluated with 15 guard digits (`evaluate_hh.li22`).
An earlier version, split at t = 1/2 only, lost up to 21 digits for such arguments (a relative error
of 10⁻²⁹ at 51 digits). All points, those of the scan and those of the blind test in
`../r0/data/points_hh`, were recomputed with the present one on 3 November 2026 (`make rescan`,
`run_hh_new_velocities.py --recompute`). The colour terms changed by up to 3 × 10⁻¹⁶ (n_f/N, gauge 2,
β = 1/2, δ = 10⁻⁷), the N² terms by up to 2 × 10⁻¹⁷ (gauge 2) and 6 × 10⁻²⁷ (gauge 1); n_f² did not
change. Of the numbers quoted in the papers only the precision control moved: δ = 10⁻¹⁰ with 15 more
digits now changes the physical colour sum by less than 8 × 10⁻³⁵ instead of 2 × 10⁻²⁶, since the earlier
routine lost different digits at the two precisions; every fit, deviation, window shift, r₀ value and
blind-test result agrees with the earlier one at the precision quoted. The direct evaluation of the
complete expressions in `hh_two_loop/stored_forward_results.jsonl` was made with the earlier routine;
its agreement with the sum of the colour terms, 4 × 10⁻¹⁸, is set by that routine (gauge 2, β = 1/4,
δ = 10⁻⁸).

`make li22` (after `make fetch`) compares the evaluation, for every argument with which the
expressions call it on the lines of the scan, at its precision control and on the two lines of the
blind test of r₀ (`../r0/data/points_hh`), with the representation obtained by summing over the
smaller index first, ∫₀¹du (−ln u)[b Li₂(a) − Li₂(abu)/u]/(1 − bu), split in the same way and
evaluated with twice the digits (logs/check_li22_values.log). Run over all 41 points on 3 November 2026
(10 086 argument evaluations, the branch point as close as 1.1 × 10⁻¹³ to the end of the path; about 11
CPU hours), it finds a largest relative difference of 8.6 × 10⁻⁴⁴, and at every point at most 0.13 of
the working precision; App. A of Paper II quotes it. Where that comparison shows a
difference far above the working precision, `check_li22_impact.py BETA:K ...` names the arguments, says
which of the two evaluations loses the digits (against `evaluate_hh.li22` at twice the digits), checks that it
reproduces the stored point, and gives the largest relative change of each colour term when the scan's
Li₂,₂ values are replaced by those at twice the digits: the effect on the stored amplitude values.

The gg → AA amplitudes of Fig. 2 are the one-loop finite remainders of arXiv:1909.08993
(ancillary files `Finite_1Loop_P1` and `Finite_1Loop_P2`), evaluated with the P1 package,
in which this one lives and which downloads and checks them (`make inputs` there). Their
values along the curves of the figure are in `one_loop/results_aa_curve.json`.

The evaluated points in `hh_two_loop/points/` are included, so the analyses, tables, figures
and fast tests need no download. Recomputing them needs the ancillary files: `make rescan
WORKERS=N` replaces every point (about 3 CPU hours), `make scan` computes only missing points and
keeps those computed with the present code (a point computed with other code, inputs or
parameters stops the scan until it is recomputed).

## Requirements

Python ≥ 3.10 with mpmath, sympy, numpy, scipy and matplotlib (`requirements.txt`). No
Mathematica, no compiled code.

## Licence and citation

Licence: to be added (`LICENSE`). If you use the package, please cite the paper
(`CITATION.cff`).
