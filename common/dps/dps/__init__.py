"""One code for the DPS analysis of gg -> HH, gg -> AA and gg -> ZZ (Paper 2).

The chain is the same for every process: evaluate the published amplitudes on the forward line
at fixed velocity (backends), fit the powers of G = (1/2) ln[(m^4 - t u)/m^4] with power
corrections in p_T^2 (fit), and compare the coefficients with the factorised DPS formula
(theory).  The modules:

    kinematics   the forward line in units of m: x, y, z, p_T^2, G, cos(theta)
    fit          interpolation fits in G with p_T^2 corrections, uncertainty from the window shift
    theory       F_A, F_S, W, the anomalous dimensions, the two-loop predictions, r_0
    points       point records: one JSON file per point, the same fields for every process
    backends     hh (Banerjee et al.), zz (VVamp) and aa (arXiv:1909.08993 through the P1 package)
    analysis     the fits of the leading logarithms and of r_0 for any process
"""
