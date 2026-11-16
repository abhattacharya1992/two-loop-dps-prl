# FDH records before the follow-up regeneration

These seven files were preserved byte for byte from the initial draft PR head,
`2b5cbeb6b6dd7979b1a39de925a2dc43791b9f35`, before the 22 October follow-up.
`SHA256SUMS` records their exact contents. The old check quotes a residual of
`1.2779728801307832e-10`; the separate older `stored/` record quotes
`2.1717718860454835e-13`. Both passed the existing numerical tolerance, but their
original producer provenance is not established by those numbers.

The active FDH files were regenerated from fresh symbolic integrands and fresh
sector decompositions using the reviewed scripts. Their metadata records actual
source/input hashes, parameters and runtime. The CDR and scalar Laurent inputs
to the scheme comparison remain explicitly identified stored inputs; this was
not a new full quark/gluon/scalar production campaign. For scope and evidence,
see the [follow-up audit at the preserved revision](https://github.com/abhattacharya1992/DPS_papers/blob/550b060f7aee329bdab183de3aa7919165ea2a2a/docs/AUDIT_FOLLOWUP_2026-10-22.md)
The original `checks/reproduce_fdh.py` is not included in this
release; the retained records and active `paper2/r0/theory/check_fdh.py`
provide the available verification evidence.
