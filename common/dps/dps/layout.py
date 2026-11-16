"""Where the stored points and results of the other processes are, in the three layouts this package
lives in: the master folder (common/dps next to paper1/ and paper2/), the Paper 2 folder (paper2/dps
next to paper2/hh, paper2/zz, paper2/r0) and the task folders (L_unified_code_ggAA/dps next to
B_colour_test_HH, J_ggZZ_quark_splitting, K_kernel_boundary_term); and where the P1 package is:
`common/aa_source/` of the master folder, a legacy `paper1/` folder, or the
directory named by the environment variable P1_ROOT."""
import os
import pathlib

PKG = pathlib.Path(__file__).resolve().parents[1]          # the dps/ directory
P2 = PKG.parents[1] / 'paper2'                             # the master layout: <master>/common/dps
MASTER = {
    'hh_points': P2 / 'hh' / 'hh_two_loop' / 'points',
    'hh_colour': P2 / 'hh' / 'hh_two_loop' / 'results_colour.json',
    'zz_points': P2 / 'zz' / 'points',
    'zz_results': P2 / 'zz' / 'results_zz.json',
    'r0_hh': P2 / 'r0' / 'data' / 'r0_from_data.json',
    'r0_zz': P2 / 'r0' / 'data' / 'r0_from_zz.json',
}
PACKAGE = {
    'hh_points': PKG.parent / 'hh' / 'hh_two_loop' / 'points',
    'hh_colour': PKG.parent / 'hh' / 'hh_two_loop' / 'results_colour.json',
    'zz_points': PKG.parent / 'zz' / 'points',
    'zz_results': PKG.parent / 'zz' / 'results_zz.json',
    'r0_hh': PKG.parent / 'r0' / 'data' / 'r0_from_data.json',
    'r0_zz': PKG.parent / 'r0' / 'data' / 'r0_from_zz.json',
}
TASKS = {
    'hh_points': PKG.parents[1] / 'B_colour_test_HH' / 'points',
    'hh_colour': PKG.parents[1] / 'B_colour_test_HH' / 'results_colour.json',
    'zz_points': PKG.parents[1] / 'J_ggZZ_quark_splitting' / 'points',
    'zz_results': PKG.parents[1] / 'J_ggZZ_quark_splitting' / 'results_zz.json',
    'r0_hh': PKG.parents[1] / 'K_kernel_boundary_term' / 'stageC' / 'r0_from_data.json',
    'r0_zz': PKG.parents[1] / 'K_kernel_boundary_term' / 'stageC' / 'r0_from_zz.json',
}
AA_POINTS = PKG / 'points_aa'


def path(key):
    """The path for `key` in whichever layout has it (the master layout first)."""
    for layout in (MASTER, PACKAGE, TASKS):
        if layout[key].exists():
            return layout[key]
    return MASTER[key]


def p1_root():
    """AA amplitude sources: prefer P1_ROOT, then common/aa_source, then a legacy source tree."""
    env = os.environ.get('P1_ROOT')
    if env:
        return pathlib.Path(env).resolve()
    master = PKG.parents[1]                                  # <master> for common/dps, the P1 folder for paper2/dps
    for cand in (master / 'common' / 'aa_source', master / 'paper1', master):
        if (cand / 'src' / 'python' / 'ggaa').is_dir():
            return cand.resolve()
    return (master / 'common' / 'aa_source').resolve()
