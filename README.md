# Double Parton Scattering Singularities in Two-Loop Gluon-Fusion Amplitudes

This repository contains the source code, numerical results and validation tools accompanying the manuscript:

**Double parton scattering singularities in two-loop gluon-fusion amplitudes**

The study investigates the leading double parton scattering (DPS) singularity in two-loop gluon-fusion amplitudes with massive colour-singlet final states. It analyses the logarithmic structure, finite matching coefficients and their dependence on colour, helicity and parity.

The calculations cover three QCD processes:

- $gg \to HH$: scalar pair production
- $gg \to AA$: pseudoscalar pair production
- $gg \to ZZ$: vector boson pair production through massless quark loops

The analysis also includes comparisons with planar $\mathcal N=4$ super Yang–Mills amplitudes and conditional predictions for the leading three-loop logarithms.

## Scientific results

The main results are:

1. The two-loop logarithmic structure of the leading DPS contribution, obtained using renormalisation group evolution and an analysis of the relevant momentum regions.
2. Finite matching coefficients $r_0(x)$ for gluon and quark splittings, tested against two-loop amplitudes.
3. Analytic helicity and parity relations for the rational contributions to the gluon-splitting coefficient.
4. Numerical tests of the predicted logarithms and finite terms in $gg\to HH$, $gg\to AA$ and $gg\to ZZ$.
5. Comparisons with known planar $\mathcal N=4$ amplitude results through five loops.
6. Conditional constraints on the leading logarithmic coefficients of three-loop QCD amplitudes.

The complete classification of pinched surfaces and matching onto renormalised operators remain open. The three-loop results depend on the additional boundary conditions specified in the manuscript.

## Repository structure

| Directory | Contents |
|---|---|
| `paper2/manuscript_prl/` | Letter, Supplemental Material, figures, tables, bibliography and arXiv source package |
| `paper2/_build/` | Scripts required to generate manuscript tables and figures |
| `common/dps/` | Shared DPS calculations, numerical data and fitting procedures |
| `paper2/hh/` | Scalar pair amplitude analysis and validation |
| `paper2/zz/` | Vector boson pair amplitude analysis and validation |
| `paper2/r0/` | Finite matching coefficients and region calculations |
| `paper2/kernel/` | One-loop DPS kernel and boundary analysis |
| `paper2/n4/` | Planar $\mathcal N=4$ amplitude comparisons |
| `checks/` | Independent analytical and numerical validation |
| `checks/threshold_ggAA/` | Comparison with the published pseudoscalar threshold expansion |
| `checks/pinch_regions/` | Exact scalar-pinch identities and bounded checks on supplied momentum-region data |
| `checks/dps_links/` | Exact checks of the coupling structure of the massless quark loop and of the relation of the one-loop DPS term to the lowest-order splitting kernels of double parton distributions |
| `common/aa_source/` | Published pseudoscalar amplitude evaluator and input checksums, retained as a shared dependency |

The repository contains one manuscript and its Supplemental Material. Earlier manuscript versions and unrelated development files are not included.

## Requirements

The Python dependencies are provided in `requirements-pinned.txt`.

The package uses Python 3.12, NumPy, SciPy, SymPy, mpmath and Matplotlib.

Manuscript compilation additionally requires:

- GNU Make and Bash
- pdfLaTeX and BibTeX
- REVTeX 4.2
- TikZ-Feynman

### Installation

```bash
git clone https://github.com/abhattacharya1992/two-loop-dps-prl.git

cd two-loop-dps-prl

python3 -m venv .venv
source .venv/bin/activate

python -m pip install -r requirements-pinned.txt
```

## Reproducing the results

A fresh clone can run the included scientific tests, regenerate the Letter and
Supplement, and check the arXiv archive with one command:

```bash
make reproduce
```

This requires the pinned Python dependencies and the TeX packages listed above.
For installation commands, the individual reproduction routes, and a precise
account of inputs not redistributed here, see [REPRODUCING.md](REPRODUCING.md).
Use `make test-science` to run the scientific checks without TeX.

### 1. Scientific validation

Run the complete included test suite:

```bash
make test
```

This runs the shared QCD analyses, independent Stage 1 and Stage 2 checks, threshold comparisons, planar $\mathcal N=4$ checks, and the kernel and finite matching tests.

The default tests use the numerical data included in this repository. They verify the published analysis, the algebraic relations, the fitting procedures and the consistency of the stored results.

Some tests requiring externally supplied full-amplitude inputs are excluded or skipped. A successful default test does not constitute a new evaluation of every original two-loop amplitude.

### 2. Independent threshold validation

The DPS predictions are checked against the earlier pseudoscalar threshold calculation:

A. Bhattacharya, *Two-loop $gg\to AA$ amplitudes near threshold and the double parton scattering singularity*, [arXiv:2610.19024](https://arxiv.org/abs/2610.19024).

Run:

```bash
make threshold
```

The checker verifies seven coefficients symbolically at $x=1/2$, including the one-loop DPS terms, the two-loop logarithmic coefficients and the single logarithm after the infrared subtraction scheme conversion.

It also verifies the scalar-to-pseudoscalar helicity-weight ratio $9/7$.

This provides a consistency check between the threshold and forward-limit analyses of the same underlying amplitude. It does not reproduce the complete threshold expansion through $\beta_A^6$, which belongs to the earlier paper.

### 3. Manuscript reproduction

Regenerate the tables, figures and PDFs:

```bash
make build
```

The resulting documents are:

- `paper2/manuscript_prl/letter.pdf`
- `paper2/manuscript_prl/supplement.pdf`

To generate and independently compile the arXiv source package:

```bash
make archive
```

The archive is written to:

`paper2/manuscript_prl/release/PRL_DPS_arxiv_sources.tar.gz`

The **arXiv archive includes the Letter and Supplemental Material source files at its top level**, together with an `anc/` directory containing the tracked scientific code, stored amplitude points, numerical results, tests and provenance documentation. The originating Git commit is recorded in `anc/README_ANCILLARY.txt`. Published third-party full-amplitude expressions are not redistributed here.

Ancillary packaging requires tracked inputs to match the current Git commit. During editing, use `make build`, or check compilation without ancillary data using `python3 paper2/manuscript_prl/make_arxiv.py --check --no-anc`. Commit the final changes before running `make archive`. The standard test suite preserves the tracked HH bookkeeping and RG reference outputs.
If running optional numerical regeneration commands, inspect `git status` and
commit or restore only independently checked by-products before packaging.

The manuscript's canonical LaTeX sources are `letter.tex` and `supplement.tex` in `paper2/manuscript_prl/`.

## Historical records

The initial FDH comparison and superseded kernel evaluators remain available
in Git history at commit `1259ee1fc9dedb42047cd3919f655713dd631f27`. They are not used
by the current publication or validation workflows. The original FDH checksums
and the recovery command are documented in
[the historical FDH record](paper2/r0/theory/stored/history/2026-10-22-pr1-initial-fdh/README.md).
The current numerical checkpoints and independent reference implementations
remain in the active repository.

## PRL manuscript v2

This branch starts from the validated cleanup commit `6d54247406102e201bb09e0d365ab1c1d3a0d940`.
It contains revised manuscripts and figures, independent coupling and DPD
splitting checks, and an additional near-endpoint validation at `beta=19/20`.
The original `main` and frozen cleanup branch are unchanged.
The endpoint summary in `paper2/r0/data/endpoint_test.json` is regenerated
by the manuscript table target and checked by the r0 regression tests.

The frozen-baseline `make release-check` tests identical publication sources;
this new manuscript intentionally differs. Validate v2 with `make reproduce`
and independently compare the unchanged scientific inputs and calculations.

## Frozen-baseline comparison

The scientific equivalence check against the original `main` commit runs in
separate temporary worktrees. It compares retained numerical inputs, re-runs
the included analyses, checks the generated Letter and Supplement, and measures
representative runtimes. The command does not update either branch:

```bash
make release-check
```

For complete instructions and the distinct external-input stage, see
[release comparison instructions](checks/RELEASE_CHECK.md).

## Blind-point provenance

The HH and ZZ independent-velocity point data are copied without modification from their source repository. Their exact Git blob identities can be checked using:

    python3 checks/blind_point_provenance.py

This verifies the source data rather than recalculating the external two-loop amplitudes. See [blind-point provenance](checks/BLIND_POINTS.md).

## Further validation

The repository retains the relevant amplitude evaluation and region calculation tools. More extensive tests require the original external inputs.

**Pseudoscalar amplitudes:** The published two-loop inputs are available through [arXiv:1909.08993](https://arxiv.org/abs/1909.08993). Their input manifests and relevant evaluation code are included under `common/aa_source/`.

**Scalar amplitudes:** The two-loop $gg\to HH$ expressions are associated with [arXiv:1809.05388](https://arxiv.org/abs/1809.05388). To obtain and test them:

```bash
make -C paper2/hh fetch
make -C paper2/hh test-full
```

**Vector boson amplitudes:** Fresh $gg\to ZZ$ evaluations require the original VVamp inputs. Instructions for obtaining, checking and compiling them are provided in `paper2/zz/README.md`.

**Finite matching and regions:** The stored Laurent coefficients and related checks are included. Recomputing the underlying region integrals is a separate calculation documented in `paper2/r0/theory/Makefile`.

**Planar $\mathcal N=4$:** The repository includes comparisons against published amplitude coefficients. Instructions for checks using the primary multi-Regge input data are provided in `paper2/n4/README.md`.

## Reproducibility and scope

The included tests distinguish analytical identities, numerical comparisons with retained data, and calculations requiring external amplitude inputs.

Agreement with the retained numerical results supports the stated relations within their tested kinematic ranges. It does not independently establish the completeness of the region analysis or an all-order factorisation theorem.

The three-loop logarithmic coefficients are predictions within the boundary assumptions described in the Letter and Supplemental Material.

Some validation scripts regenerate local numerical records. Inspect `git status` after running tests before committing generated changes.

## License

The repository is distributed under the terms specified in `LICENSE` and `common/aa_source/LICENSE`. Third-party inputs remain subject to their respective licenses.

## Acknowledgments

See the manuscript.
