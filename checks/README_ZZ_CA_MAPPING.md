# Regenerating the ZZ crossed-master mapping check

`check_zz_ca_mapping.py` identifies the two crossed nonplanar master families
from six weight-four coefficients of the published bare and finite FORM exports.
It does not use a finite C++ amplitude value to choose the mapping.

## Retained evidence and fresh evaluation

The compact result is
`paper2/zz/scripts/fixtures/nonplanar_weight4_residues.json`.
The default test suite compares its six retained coefficients, rejects the
interchanged mapping, and rejects incomplete or nonfinite records:

```bash
python3 paper2/zz/scripts/test_master_families.py
```

That test rechecks the retained result. To regenerate the coefficients from the
published expressions, obtain these two checksum-pinned files in
`paper2/zz/inputs/`:

- `gg_Aj_A-1.0.inc.gz`
- `gg_Aj_A_qT-1.0.inc.gz`

The checker validates both files against `paper2/zz/inputs/SHA256SUMS` before
evaluating them. The files are external inputs; they are not redistributed here.
See the input instructions in `paper2/zz/README.md`.

The restoration was also validated on 15 November 2026 by downloading those
published expressions, checking both hashes, and evaluating all six comparisons
with the restored reader. All comparisons passed and the interchanged mapping
was rejected; the regenerated coefficients agree with the retained fixture at
every serialized digit. The fresh result and compiler/library provenance are in
[`results/zz_ca_mapping_20261115.json`](results/zz_ca_mapping_20261115.json).
Its local directory prefixes are normalized, as documented in that record.

The native calculation additionally requires a C++17 compiler, GiNaC and CLN
development files, and Linux `ldd`. On Debian or Ubuntu, `libginac-dev` supplies
the GiNaC headers and pulls in the CLN development dependency. Python requirements
are listed in the repository's `requirements-pinned.txt`.

From the repository root, use a fresh working directory:

```bash
validation_dir=$(mktemp -d)
python3 checks/check_zz_ca_mapping.py \
  --finite paper2/zz/inputs/gg_Aj_A_qT-1.0.inc.gz \
  --bare paper2/zz/inputs/gg_Aj_A-1.0.inc.gz \
  --work "$validation_dir/run" \
  --output "$validation_dir/zz_ca_mapping.json"
```

The default installation prefix is `/usr`, with libraries under
`lib/x86_64-linux-gnu`; `--prefix` selects an equivalent installation in another
prefix. The output records source, executable and shared-library hashes,
compiler details, precision settings, all six comparisons, and rejection of the
interchanged mapping. The `--work` path must not already exist.

## Native parser integrity checks

The parser checks can run without the large published amplitude expressions:

```bash
make -C paper2/zz native-parser-test
```

This target compiles the reader in a temporary directory and removes the build
after the checks. `CXX`, `CPPFLAGS`, `CXXFLAGS` and `LDFLAGS` select a compiler or
library installation when needed.

They test operator precedence, exact division at the requested precision,
complex values, inverse colour powers, function arities, malformed input,
unsupported functions, and singular denominators.

## Restored source provenance

The native parser and its integrity tests are restored byte-for-byte from
`abhattacharya1992/DPS_papers` at commit
`65efa9a70c73d56aa3f55cc612e3d412603d63f5`:

| File | Git blob SHA-1 |
| --- | --- |
| `checks/native_form.cpp` | `b16c9138cc64034b954e2ea591a64e050ce2f126` |
| `checks/test_native_form.py` | `93ba6e0c3b8b11239ea334ff319a7bfb6a23052a` |

The small `libraries()` function from that commit's
`checks/regenerate_native_pairs.py` is retained directly in
`check_zz_ca_mapping.py`. The full HH/AA regeneration driver is not a dependency
of this ZZ check. This restoration adds the missing execution dependencies;
it does not change the stored residue values or the master-family assignment.
