# Initial FDH comparison snapshot

The seven records listed in SHA256SUMS document the initial FDH check from
22 October 2026. Their contents and checksums remain available in Git history
at commit `1259ee1fc9dedb42047cd3919f655713dd631f27`. They are not active inputs for the
Letter, Supplement, or current numerical tests.

To retrieve any record from a clone containing that commit, use:

```bash
git show 1259ee1fc9dedb42047cd3919f655713dd631f27:paper2/r0/theory/stored/history/2026-10-22-pr1-initial-fdh/FILE > FILE
```

Replace `FILE` with the filename in `SHA256SUMS`. The current FDH results
and checkpoints remain in `paper2/r0/theory/`. The archival checksum manifest
is retained here so the earlier comparison can still be identified and recovered.

The initial check reported a residual of `1.2779728801307832e-10`;
another older stored record reported `2.1717718860454835e-13`.
Both passed the existing tolerance, but the residuals do not establish
the original producer provenance.

The active FDH files were regenerated from fresh symbolic integrands and fresh
sector decompositions using the reviewed scripts. Their metadata records actual
source/input hashes, parameters and runtime. The CDR and scalar Laurent inputs
to the scheme comparison remain explicitly identified stored inputs; this was
not a new full quark/gluon/scalar production campaign. For scope and evidence,
see the [follow-up audit at the preserved revision](https://github.com/abhattacharya1992/DPS_papers/blob/550b060f7aee329bdab183de3aa7919165ea2a2a/docs/AUDIT_FOLLOWUP_2026-10-22.md)
The original `checks/reproduce_fdh.py` is not included in this
release; the retained records and active `paper2/r0/theory/check_fdh.py`
provide the available verification evidence.
