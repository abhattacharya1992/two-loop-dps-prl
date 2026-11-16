# Published gg -> AA amplitude evaluator

This directory retains the amplitude evaluator and input checksums from the earlier pseudoscalar calculation (arXiv:1909.08993). The current two-loop DPS study uses this code for optional fresh AA amplitude comparisons and one-loop checks. It is **not** a second manuscript.

The Python module `ggaa` is in `src/python/ggaa/`. Obtain the published input expressions with `bash inputs/fetch_inputs.sh`; `inputs/SHA256SUMS` verifies the downloaded files. Those external expressions are not redistributed.

The code retains its BSD 3-Clause license in `LICENSE`. The DPS analysis finds this directory automatically. For an external compatible source tree, set `P1_ROOT` to its root directory.
