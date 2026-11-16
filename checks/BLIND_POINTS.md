# Blind-velocity data provenance

The two-loop DPS finite matching coefficient was identified using the main scan velocities. Its independent blind tests use beta = 1/10 and 4/5 (gg -> HH), and beta = 4/5 (gg -> ZZ).

## Source

The 32 HH input files in paper2/r0/data/points_hh and the 9 ZZ input files in paper2/r0/data/points_zz are byte-identical to the files in abhattacharya1992/DPS_papers, branch review/prl-physics-revision, commit b56dcd8f092f6b721bced03d3bcfbb617c14dab1. The two saved results and three extraction/source scripts are also identical. The Git object SHA-1 digests for all 46 files are in blind_points_manifest.json.

From the project root, check the recorded sources and rerun the fits:

    python3 checks/blind_point_provenance.py
    make -C paper2/r0 test

These operations verify byte-level origin and re-extract r0 from stored points. They do not independently regenerate the original loop amplitudes.

## Fresh input generation

HH points can be regenerated with paper2/r0/data/run_hh_new_velocities.py, which uses paper2/hh/hh_two_loop/scan_colour_scale.py, after obtaining/checking published input expressions from arXiv:1809.05388 (see paper2/hh/README.md). Work in a fresh copy to preserve committed points; the driver accepts --recompute to bypass its cache.

ZZ points require the external VVamp inputs and their translation/compilation (see paper2/zz/README.md). The generator paper2/zz/scripts/zz_points.py was used at beta = 4/5, delta = 10^(-k) with k = 4,...,12 and 30+5k digits. Extract using paper2/r0/data/blind_test_zz.py. HH extracts with paper2/r0/data/blind_test_hh.py.

No expensive fresh original-amplitude evaluation was performed solely to create this provenance record.
