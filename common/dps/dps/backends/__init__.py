"""The amplitude backends.  Each one evaluates the published amplitudes of its process at a point
of the forward line and writes a point record (points.py); the analysis never sees the differences.

    hh   the two-loop class-B gg -> HH remainders of Banerjee et al. (hh/hh_two_loop/evaluate_hh.py)
    zz   the one- and two-loop gg -> ZZ helicity amplitudes of VVamp (zz/scripts/zz_points.py)
    aa   the one- and two-loop Type II gg -> AA remainders of arXiv:1909.08993, through the P1
         package (ggaa.twoloop), with the colour decomposition done in one pass (colour.py)

The hh and zz backends call the scan drivers of those folders unchanged; this package adds the
aa backend and the common interface `evaluate(component, beta, k, dps, out_dir)`.
"""
