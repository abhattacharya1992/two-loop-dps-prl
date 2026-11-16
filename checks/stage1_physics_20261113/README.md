# Bounded Stage 1 physics checks

Start with [the status summary](../../docs/STAGE1_PHYSICS_STATUS_20261113.md)
and [the derivation](../../paper2/derivations/stage1_physics_20261113.md).

The checker independently integrates the RG equation and constructs the
Fourier–Bessel transforms, evaluates a static rapidity-overlap integral by a
separate quadrature route, verifies colour matrices against explicit generators,
and tests the stated Landau and rapidity-separation identities. It does not
import production DPS modules or claim a complete factorisation theorem.

Run from the repository root with Python, SymPy, mpmath and NumPy:

```sh
python checks/stage1_physics_20261113/verify_stage1.py --output /tmp/stage1-fresh.json
python checks/stage1_physics_20261113/test_stage1_guards.py --output /tmp/stage1-guards.json
```

Use new output paths. `results.json` retains the 66 named checks, source hash,
runtime versions, results and exclusions. `guard_results.json` records replay,
optimized-Python, output-preservation and deliberate-error controls. Neither
program modifies the production data or the retained records. The scripts are
standalone: the root `make test` does not invoke them automatically.

The prescribed sign and homogeneous-integral convention are documented in the
note. A passing result verifies these mathematical statements and inputs; it
does not determine the unknown coefficients of possible three-loop boundary
terms or prove that every physically allowed region has been included.
