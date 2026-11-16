# Scalar DPS pinch cross-checks

This publication component constructs scalar Symanzik polynomials and checks the local DPS pinch. It does **not** recompute QCD amplitudes or enumerate all momentum-space regions.

```bash
make -C checks/pinch_regions check
```

The calculation needs Python and SymPy from the repository requirements. It has no network access or private repository dependency. The XZ-compressed JSON is data, not executable or pickle content.

## What is computed

`exact_checks.py` independently constructs the tree/forest and loop-quadratic-form polynomials for **54 canonical graphs**, with `x`, `y` and `pT^2/s` symbolic. It checks the one-loop box identity and proves positive-interior stationary-point exclusion for **16 interjet exchanges** at equal virtualities, at arbitrary interior fractions. These certificates do not exclude boundary pinches or regulator-defined Glauber terms.

`param_space.py` reconstructs leading polynomials at **2823 supplied rows / 2363 distinct graph-parameter vectors**, in the 70 parametrisations retained by the original survey. The supplied labels distinguish DPS-type and SPS-type reduced diagrams; pointwise powers $P$ include diagram numerators whereas $P_{den}$ counts only measure and denominators. These powers are inputs of the momentum-space survey, not outputs of the scalar program. All 1709 hidden rows / 1305 distinct hidden vectors have two mixed-sign factors, one per incoming family. The same algebra occurs in 678 SPS-labelled rows; momentum labels are explicitly not independent output. All 142 selected hard/collinear/ultrasoft scalar-power comparisons agree with the supplied measure-and-denominator powers; 127 of these have $P_{den}\le 0$. The extra 15 are ultrasoft and power suppressed ($P_{den}=2$). None is a test of the QCD numerator powers. The code retains the 422/85 direction of the other local-power discrepancies and the 140 cases with a logarithmic normal integral (`b=1`). These diagnostics are not a global region count or integrated contributions.

Sign witnesses are constructive and rational. Local saddle powers use `(A D-B C)/A`, not a numerical search for a zero of the leading factors. The result is deterministic, including under a different Python hash seed and `python -O`.

## Provenance

The input enumeration originates from an exploratory momentum-space survey separate from this release. Its full scan programs are not distributed here, so the released scalar checks are independently reproducible *conditional on the supplied virtuality vectors*, but the original momentum-space selection is not independently regenerated. The graph routing is adapted from BSD-3-Clause source code.

`vectors.json.xz` is a pinned compact transcription of the supplied virtuality records. Its source replay file SHA256 is `41c4747fba61fbb274ab89c7cfce00f973eb2483cf3512f290ce26d05780d207`. The JSON records its schema, graph/mode dictionaries, source commit and label warning; the compressed-byte SHA256 is enforced by the checker. The historical surface labels are *input*, not a region-by-region classification independently derived by this package. Alternate loop routings are not regenerated from the original momentum-space survey.

## Controls and scope

Tests cover exact sign witnesses, inappropriate use on nonmultiaffine polynomials, shifted saddles, false jet labels, an interior stationary positive control, a deliberately wrong polynomial, corrupt input bytes, wrong scalar powers and an always-true jet classifier. The positive stationary control is a bilinear pinch, not a spectator-Glauber graph; that separate control remains outside this package.

No full dissection, numerator scan, axial-pole prescription, regulated overlap integration or operator matching is supplied here. A nonisolated homogeneous scaling is not automatically a vanishing regulated contribution. The 223 source rows with tangent sign changes require their own analysis; the local chart test is not that analysis. No source research files are silently repaired or overwritten.
