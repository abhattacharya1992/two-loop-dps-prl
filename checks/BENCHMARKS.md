# Historical size and performance benchmarks

## Two independent questions

**Repository compactness**: raw tracked Git blob sizes, counted for a specified
revision. This includes scientific inputs, publication material and source code.
It is NOT a measure of packed .git history, fresh clone size, or the number of
algebraic terms in an amplitude expression.

**Analysis pipeline speed**: running the same three Make targets in detached
worktrees at pinned Git revisions, on the same runner and dependencies. The
reference commands exercise stored point analysis, not fresh expensive
third-party two-loop amplitude evaluations or a beta expansion.

The original reference is pinned at
1259ee1fc9dedb42047cd3919f655713dd631f27, and frozen PRL v2 at
d362f551357f80125649e583b22b6dbf8aeccb4e.
These remain meaningful after the repository main branch was merged or renamed.

## Commands

Run from the repository root after installing requirements-pinned.txt:

~~~~bash
make benchmark-size PY=python
make benchmark-pipeline PY=python PERF_RUNS=5
make benchmark-v2 PY=python PERF_V2_RUNS=3
~~~~

Size outputs:
- benchmark-results/size/report.json
- benchmark-results/size/report.md

Pipeline outputs:
- benchmark-results/pipeline-original/report.json and report.md
- benchmark-results/pipeline-v2/report.json and report.md
- Individual command logs for each test and repetition

The report filenames above are local generated outputs and should not be
committed as scientific source inputs. GitHub Actions uploads the benchmark
results as downloadable artifacts with the exact baseline and candidate SHA.

The scientific targets are:
- make -C common/dps test
- make -C paper2/zz test
- make -C paper2/r0 test

Three types of result are reported: median wall-clock time and the per-run
distributions, maximum peak RSS in KiB if /usr/bin/time is available, and any
nonzero command status. The timing ratio is baseline seconds divided by
candidate seconds; greater than one means the candidate was faster.

The implementation alternates baseline and candidate first within successive
rounds. Both worktrees are removed afterward; no branch is updated. The
benchmark always records the original and exact candidate commits.

## Historical tracked-content snapshot before T07

Previously measured with full recursive Git tree listings:

| Revision | Tracked files | Tracked bytes | MiB |
|---|---:|---:|---:|
| Original main (1259ee1) | 772 | 12,233,790 | 11.67 |
| Frozen PRL v2 (d362f55) | 818 | 11,288,288 | 10.77 |
| Cleanup before Phase 1 (d9052d6) | 812 | 7,213,931 | 6.88 |

The latest cleanup includes the subsequent documentation and benchmark code,
so its live size is recalculated instead of copied from these historical values.

Do not call this a 41% reduction in executable code or a 41% smaller Git clone.
The main reduction is optional performance cache data, retained in Git history.

## Scientific equivalence and limitations

Normal scientific reproducibility is independently checked by make reproduce.
Frozen PRL v2 equivalence is independently checked by make release-check.
The benchmark targets contain different amounts of unit-test code across
revisions; therefore an end-to-end speed comparison of a Make target is not
a matched microbenchmark of identical executable instructions.

A future full-amplitude versus beta-series comparison requires the actual
validated expressions in each process, a shared test domain, matched arithmetic
precision, warm-up control, separate load/preparation and evaluation timings,
operation counts, expression serialization sizes, and numerical-error bounds.
That larger mathematical benchmark is NOT implemented or implied by these
stored-data pipeline results.


## Consolidated main (19 November 2026)
The original main tree (1259ee1) contained 772 files / 12,233,790 raw
tracked bytes (11.667 MiB). The first consolidated main tree (42d808a)
contained 829 files / 7,332,078 bytes (6.992 MiB): 40.07% fewer
tracked bytes than the original. Later code/CI additions change the
current count; regenerate it with make benchmark-size.
These numbers are NOT compressed Git history or clone size.
