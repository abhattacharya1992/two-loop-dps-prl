# Benchmarked results (frozen snapshot: 19 November 2026)

![Results overview](00-performance-overview.svg)

**Scope:** This page separates AA **analytic threshold-series compression** from HH/AA/ZZ **lossless storage compression** and HH/ZZ **original-amplitude precision convergence**. None of the HH/ZZ results is an analytic beta-series speedup. The plot is a dated snapshot, not a live dashboard. All values below are based on completed GitHub Actions results and the pinned upstream AA data.

## 1. AA analytical series: expression size and repeated evaluation

| Metric | Original exact expression | β⁶ threshold series | Improvement |
|---|---:|---:|---:|
| Wolfram expression representation | 44,844,413 bytes | 637,268 bytes | **70.37× smaller** |
| Repeated evaluation, median | 342.209 ms | 17.842 ms | **19.18× faster** |

The warm timing was measured previously in **Wolfram 13.2** on AMD Ryzen 9 7950X, at β=1/10, cosθ=0 and **80-digit precision** (seven repeats). Expression byte counts compare matching Wolfram representations, *not* the original published file byte counts against the compact six-coefficient files. The warm timings exclude initial preparation, and are not the new one-shot Python benchmark.

Source: [pinned AA upstream commit](https://github.com/abhattacharya1992/ggAA-threshold/tree/76491fd105e6f498be51affa213fe6839f879ce6), using `data/timing/timing_run.json` and `data/timing/timing_benchmark.csv`. Historical records reverified by [compression CI #44588650024](https://github.com/abhattacharya1992/two-loop-dps-prl/actions/runs/44588650024).

### Separate fresh AA accuracy-checked comparison

At β=1/10, cosθ=1/3 and 120-digit Python precision, the exact calculation took **237.11 s** versus **0.0274 s** for β⁶, with relative errors **M1 = 6.38 × 10⁻⁹** and **M2 = 8.34 × 10⁻⁸**. This is a *single-run* comparison affected by the cost of evaluating a large original expression; **do not interpret its ~8,643× ratio as a repeatable warm-evaluation speedup**. The original 40-digit Python exact run failed accuracy, motivating the automatic accuracy gate. [Validated AA run](https://github.com/abhattacharya1992/two-loop-dps-prl/actions/runs/44556434508) · [Invalid 40-digit run](https://github.com/abhattacharya1992/two-loop-dps-prl/actions/runs/44564737847).

## 2. HH / AA / ZZ stored numerical-point compression

| Process | Verified points | Raw JSON | zlib + checksum | Saved | Raw parse | Verified packed parse |
|---|---:|---:|---:|---:|---:|---:|
| HH | 50 | 157,276 B | 75,139 B | **52.22%** | 0.0150 ms | 0.0313 ms |
| AA | 82 | 236,949 B | 121,491 B | **48.73%** | 0.0108 ms | 0.0269 ms |
| ZZ | 28 | 804,266 B | 225,887 B | **71.91%** | 0.0864 ms | 0.1806 ms |

Every packed record was decompressed and byte-for-byte verified using SHA-256. **This is lossless stored-data compression, not a beta expansion.** Parsing a verified compressed record is slower than raw JSON parsing; the gain here is disk/storage size, not loading speed. [Source: compression CI #44588650024](https://github.com/abhattacharya1992/two-loop-dps-prl/actions/runs/44588650024), artifact `compression-performance-reports`, `points/report.json`.

## 3. HH and ZZ original amplitudes: accuracy-preserving performance

| Property | HH | ZZ |
|---|---:|---:|
| Checked original-amplitude components | 2 | 108 (36 helicities × 3 colour parts) |
| Precision levels | 40, 60, 80 digits | 40, 60, 80 digits |
| Lowest passing precision at β=0.1, cosθ=1/3 | **40 digits** | **40 digits** |
| Accuracy-preserving evaluation time | **66.85 s** | **27.63 s** |
| Maximum discrepancy, 60 vs 80 digits | 1.78 × 10⁻⁵¹ | 4.34 × 10⁻⁵⁰ |
| Archived reference discrepancy | 2.81 × 10⁻²⁷ | 2.14 × 10⁻²¹ |
| Full analytic β⁶ threshold series | **Not derived** | **Not derived** |

The HH archived Bose-symmetry benchmark was generated using the same analytic evaluator (archival consistency, **not** an independent implementation). The ZZ FORM/master-integral reference independently verifies the 36 two-loop `C_F` helicity components. Precision convergence of the other ZZ components does not establish independent correctness. These are different physical processes and evaluator implementations; their absolute runtimes are not a direct complexity ranking. [HH precision run](https://github.com/abhattacharya1992/two-loop-dps-prl/actions/runs/44573748819) · [ZZ precision run](https://github.com/abhattacharya1992/two-loop-dps-prl/actions/runs/44578459993).

## Data and updating this snapshot

The adjacent [`benchmark_snapshot.json`](benchmark_snapshot.json) contains the underlying numbers, source-run links, input provenance and scientific scope. The plot [`00-performance-overview.svg`](00-performance-overview.svg) is a compact, static presentation of those values. When repeating benchmarks, update the snapshot, review scientific comparability, and regenerate the plot. Generated detailed reports are also available as CI artifacts, but these may expire; this tracked snapshot remains in Git history.
