#!/usr/bin/env python3
"""Independent, bounded Stage-1 checks. No production DPS module is imported.

This verifies specified mathematical identities, not completeness of a region list
or an all-order factorisation theorem. Requirements: sympy, mpmath, numpy.
Run from any directory with --output PATH; existing output is not overwritten.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

BASELINE = "98cd34fc7e5c707fabdebeac7401398daea4a532"
CHECKS: list[str] = []


def zero(name: str, expression) -> None:
    entries = list(expression) if isinstance(expression, sp.MatrixBase) else [expression]
    if any(sp.simplify(e) != 0 for e in entries):
        raise ArithmeticError(f"{name}: {expression}")
    CHECKS.append(name)


def overlap_checks() -> dict:
    eta = sp.Symbol("eta", real=True)
    u, v = sp.symbols("u v", positive=True)
    F = -2 * sp.pi**2 * (1 - sp.I * sp.tan(sp.pi * eta / 2))
    J10, J01, J11 = F*u**(-eta), F*v**(-eta), F*(u+v)**(-eta)
    parts = {"DD": J11, "PD": J01-J11, "DP": J10-J11,
             "PP": -J01-J10+J11}
    zero("static overlap cancels at finite eta", sum(parts.values()))
    limits = {k: sp.limit(w, eta, 0) for k, w in parts.items()}
    expected = {"DD": -2*sp.pi**2, "PD": 0, "DP": 0, "PP": 2*sp.pi**2}
    for k in parts:
        zero(f"static {k} eta limit", limits[k]-expected[k])
    if sp.simplify(sum(limits.values())-limits["PP"]) == 0:
        raise ArithmeticError("omitting parent-parent incorrectly passed")
    CHECKS.append("omitted-parent negative control is nonzero")
    zero("O(eta) PD term retained", sp.limit(parts["PD"]/eta, eta, 0)
         + 2*sp.pi**2*sp.log((u+v)/v))

    # Independent route: d=c-a makes the regulator independent of c.
    # Contour integration in c leaves a PV integral in d. Pair its positive
    # and negative halves, invert the tail, and obtain the finite interval below.
    # This is different from the baseline's nested Hilbert-transform quadrature.
    samples = []
    for es in ("-0.6", "-0.25", "0.2", "0.55"):
        e = mp.mpf(es)
        def integrand(t):
            if t == 1:
                return -e
            return (t**e-t**(-e))/(1-t*t)
        # t=z^4 resolves the integrable endpoint at t=0.
        def mapped(z):
            if z == 0:
                return mp.mpf(0)
            return 4*z**3*integrand(z**4)
        val = mp.quad(mapped, [0, mp.mpf("0.5"), 1])
        target = -mp.pi/2*mp.tan(mp.pi*e/2)
        err = abs(val-target)
        if err > mp.mpf("1e-35"):
            raise ArithmeticError(f"independent overlap quadrature: {err}")
        samples.append({"eta": es, "absolute_error": mp.nstr(err, 8)})
    CHECKS.append("four independent finite-eta PV quadratures")

    # Keep an arbitrary common transverse factor and d-dimensional Born tensor.
    e = sp.Symbol("eps")
    b0, b1, b2, z2, z1, z0 = sp.symbols("b0 b1 b2 z2 z1 z0")
    born = b0 + e*b1 + e**2*b2
    pole = z2/e**2+z1/e+z0
    zero("common factor cancels before eps expansion", sp.expand(born*pole-born*pole))
    finite_mismatch = sp.expand((born-b0)*pole).coeff(e, 0)
    zero("four-dimensional-only subtraction finite defect", finite_mismatch-b1*z1-b2*z2)
    return {"J00_prescription": "homogeneous longitudinal integral set to zero",
            "closed_parts": {k: str(vv) for k, vv in parts.items()},
            "eta_zero_limits": {k: str(vv) for k, vv in limits.items()},
            "quadrature": samples,
            "finite_defect_if_only_counterterm_uses_4D_Born": str(finite_mismatch),
            "scope": "complete static eikonal attachment sum; not unexpanded Glauber graphs"}


def transverse_checks() -> dict:
    z, delta, c1, c2, ca = sp.symbols("z Delta C1 C2 Ca", real=True)
    full = -c1-c2+delta*sp.cos(z)
    zero("coincident four-line factor is parent factor",
         (full+ca).subs(z, 0).subs(delta, c1+c2-ca))
    zero("subtracted soft factor", full-full.subs(z, 0)-delta*(sp.cos(z)-1))
    zero("ultrasoft separation suppression", sp.limit((sp.cos(z)-1)/z**2, z, 0)+sp.Rational(1, 2))
    records = []
    # Auxiliary mass only in the common scalar transverse integral, not a new
    # massive-gauge-theory prescription for the production calculation.
    for es in ("0.2", "0.35"):
        e, b, m = mp.mpf(es), mp.mpf("1.7"), mp.mpf("1e-15")
        ib = (b/m)**e*mp.besselk(e, m*b)/(2*mp.pi)**(1-e)
        i0 = m**(-2*e)*mp.gamma(e)/(4*mp.pi)**(1-e)
        limit = mp.gamma(-e)*(b*b/4)**e/(4*mp.pi)**(1-e)
        rel = abs((ib-i0)/limit-1)
        if rel > mp.mpf("1e-18"):
            raise ArithmeticError(f"mass-regulated transverse subtraction: {rel}")
        records.append({"epsilon": es, "mass": "1e-15", "relative_error": mp.nstr(rel, 8)})
    CHECKS.append("two mass-regulated transverse limits")
    e = sp.Symbol("eps")
    bracket = sp.series(sp.exp(-sp.EulerGamma*e)*sp.gamma(-e), e, 0, 3).removeO().expand()
    zero("rapidity coefficient has no finite constant at mu_b", bracket.coeff(e, 0))
    zero("rapidity coefficient residue", bracket.coeff(e, -1)+1)
    return {"subtracted_colour_factor": "Delta*(cos(ell_perp.b)-1)",
            "mass_regulated_checks": records,
            "rapidity_coefficient_series": str(bracket),
            "scope": "one-exchange soft/coincident overlap, after colour projection"}


def generators(n: int) -> list[np.ndarray]:
    out = []
    for i in range(n):
        for j in range(i+1, n):
            a = np.zeros((n, n), complex); a[i, j] = a[j, i] = 0.5; out.append(a)
            a = np.zeros((n, n), complex); a[i, j] = -0.5j; a[j, i] = 0.5j; out.append(a)
    for k in range(1, n):
        d = np.zeros(n); d[:k] = 1; d[k] = -k
        out.append(np.diag(d/np.sqrt(2*k*(k+1))).astype(complex))
    return out


def colour_checks() -> dict:
    n = sp.Symbol("N", positive=True)
    cf = (n*n-1)/(2*n)
    metric = sp.Matrix([[n*n, n], [n, n*n]])
    singlet, hard = sp.Matrix([1, 0]), sp.Matrix([0, 1])
    psi = (hard-singlet/n)/2
    bra = hard.T*metric
    norm = (bra*psi)[0]
    mean = lambda op: sp.factor((bra*op*psi)[0]/norm)
    jet = sp.Matrix([[-cf, -sp.Rational(1, 2)], [0, 1/(2*n)]])
    same = sp.Matrix([[1/(2*n), 0], [-sp.Rational(1, 2), -cf]])
    cross = sp.Matrix([[-1/(2*n), sp.Rational(1, 2)], [sp.Rational(1, 2), -1/(2*n)]])
    for label, op in (("jet", jet), ("hard", same), ("cross", cross)):
        zero(label+" is Hermitian in colour metric", op.T*metric-metric*op)
    zero("four-leg colour conservation", jet+same+cross+cf*sp.eye(2))
    zero("parent-adjoint splitting tensor eigenvalue", jet*psi-psi/(2*n))
    d = 2*cross
    zero("one-exchange colour projection", mean(d)+1/n)
    zero("two-exchange colour projection", mean(d*d)-1-1/n**2)
    zero("nonzero projected commutator", mean(same*cross-cross*same)-sp.Rational(1, 4))
    f = d+sp.eye(2)/n
    zero("zero projected one-loop operator", mean(f))
    zero("its square does not project to zero", mean(f*f)-1)
    zero("shifted dipole squares to identity", f*f-sp.eye(2))
    t = sp.Symbol("t", real=True)
    exp_d = sp.exp(-t/n)*(sp.cosh(t)*sp.eye(2)+sp.sinh(t)*f)
    zero("exact projected matrix exponential", mean(exp_d)-sp.exp(-t/n)*sp.cosh(t))
    zero("matrix exponential initial value", exp_d.subs(t, 0)-sp.eye(2))
    zero("matrix exponential evolution", sp.diff(exp_d, t)-d*exp_d)

    # Explicit generators independently verify the Fierz-derived matrices.
    max_error = 0.0
    pairs = {(0, 1): jet, (2, 3): jet, (0, 2): same, (1, 3): same,
             (0, 3): cross, (1, 2): cross}
    for dim in (2, 3, 4):
        gs = generators(dim)
        eye = np.eye(dim)
        s = np.einsum("ij,lk->ijkl", eye, eye)
        h = np.einsum("ik,lj->ijkl", eye, eye)
        basis = [s, h]
        def act(a, leg, tensor):
            matrix = -a.T if leg in (1, 2) else a
            return np.moveaxis(np.tensordot(matrix, tensor, axes=(1, leg)), 0, leg)
        for (i, j), op in pairs.items():
            m = np.array(op.subs(n, dim)).astype(complex)
            for col, tensor in enumerate(basis):
                actual = sum((act(a, i, act(a, j, tensor)) for a in gs), np.zeros_like(s, complex))
                expected = m[0, col]*s+m[1, col]*h
                max_error = max(max_error, float(np.max(abs(actual-expected))))
    if max_error > 1e-12:
        raise ArithmeticError(f"explicit-generator colour check: {max_error}")
    CHECKS.append("36 explicit-generator dipole actions in SU(2), SU(3), SU(4)")
    return {"basis": "S=delta_ij delta_lk, H=delta_ik delta_lj; legs q,qbar,qbar,q",
            "metric": str(metric), "jet_dipole": str(jet), "hard_dipole": str(same),
            "cross_dipole": str(cross), "D_mean": str(mean(d)), "D_squared_mean": str(mean(d*d)),
            "projected_hard_cross_commutator": str(mean(same*cross-cross*same)),
            "projected_exponential": "exp(-t/N)*cosh(t)",
            "explicit_generator_max_abs_error": max_error,
            "scope": "colour obstruction to projection-before-evolution; not a calculated DPS anomalous dimension"}


def evolution_checks() -> dict:
    X, G, L, t, y, alpha = sp.symbols("X G L t y alpha")
    D, B, b0, b1, gam1, dg1, hs, r0 = sp.symbols("Delta B beta0 beta1 Gamma1 Dgamma1 Hsum r0")
    dg0 = B-b0
    q0 = -2*D*(X-y)+dg0
    a1 = sp.integrate(q0, (y, 0, t))
    a2 = sp.integrate(b0*(t-y)*q0-D*gam1*(X-y)/2+dg1, (y, 0, t))
    k1 = sp.expand(a1+b0*t+r0)
    k2 = sp.expand(a2+a1*a1/2+(b0*t+r0)*a1+b0*b0*t*t+b1*t+2*b0*t*r0+hs*k1)
    k1 = sp.expand(k1.subs(t, X+L))
    k2 = sp.expand(k2.subs(t, X+L))
    z3, z5 = sp.symbols("zeta3 zeta5")
    exponential = sp.series(sp.exp(-2*G*alpha-2*z3*alpha**3/3-2*z5*alpha**5/5), alpha, 0, 6).removeO()
    f0 = sp.expand(exponential/(2*alpha))
    f2 = sp.series(exponential/(1-alpha), alpha, 0, 5).removeO().expand()
    transforms0 = [-sp.factorial(j)*f0.coeff(alpha, j) for j in range(5)]
    transforms2 = [sp.factorial(j)*f2.coeff(alpha, j) for j in range(5)]
    def transform(poly, ts):
        return sp.expand(sum(poly.coeff(X, j)*ts[j] for j in range(5)))
    one0, one2 = transform(k1, transforms0), transform(k1, transforms2)
    for power, expected in ((3, -4*D/3), (2, -B)):
        zero(f"two-loop nu0 G^{power}", one0.coeff(G, power)-expected)
    for power, expected in ((2, -4*D), (1, 4*D-2*B)):
        zero(f"two-loop nu2 G^{power}", one2.coeff(G, power)-expected)
    two0, two2 = transform(k2, transforms0), transform(k2, transforms2)
    c4 = 2*D*(B+2*b0/3)
    c3 = 2*B*(B+b0)/3-D*gam1/3-4*D*hs/3-4*D*r0/3-4*D*L*(B+b0+D*L)/3
    expected0 = {5: sp.Rational(8, 5)*D**2, 4: c4, 3: c3}
    expected2 = {4: 8*D**2, 3: 4*c4-16*D**2,
                 2: 3*c3-4*D*(3*B+2*b0-6*D)}
    for harmonic, actual, expected in ((0, two0, expected0), (2, two2, expected2)):
        for power, value in expected.items():
            zero(f"three-loop nu{harmonic} G^{power}", actual.coeff(G, power)-value)
    # A homogeneous boundary term at relative order a_s^2 is invisible to
    # the mu RGE through total order a_s^3, but its Fourier transform is not zero.
    r = sp.symbols("r0_3 r1_3 r2_3 r3_3 r4_3")
    boundary = sum(r[j]*X**j for j in range(5))
    shifts = {"nu0": transform(boundary, transforms0), "nu2": transform(boundary, transforms2)}
    zero("quadratic boundary shifts nu0 third logarithm", shifts["nu0"].coeff(G, 3)-4*r[2]/3)
    zero("quadratic boundary shifts nu2 third logarithm", shifts["nu2"].coeff(G, 2)-4*r[2]-12*r[3]-48*r[4])
    zero("quartic boundary shifts highest nu0 logarithm", shifts["nu0"].coeff(G, 5)-sp.Rational(16, 5)*r[4])
    zero("cubic boundary shifts next nu0 logarithm", shifts["nu0"].coeff(G, 4)+2*r[3])
    a = sp.Symbol("a")
    beta = -2*b0*a*a-2*b1*a**3
    residual = sp.expand(beta*sp.diff(a**3*boundary, a)-a*(-4*D*(X-t)+2*dg0)*a**3*boundary)
    for order in range(4):
        zero(f"homogeneous boundary RGE at order a^{order}", residual.coeff(a, order))
    return {"integrated_relative_kernel_order1": str(k1),
            "integrated_relative_kernel_order2_without_new_boundary": str(k2),
            "transforms_nu0": [str(sp.expand(e)) for e in transforms0],
            "transforms_nu2": [str(sp.expand(e)) for e in transforms2],
            "general_three_loop_boundary_shifts": {k: str(sp.collect(v, G)) for k, v in shifts.items()},
            "conclusion": "all six quoted coefficients agree within the stated boundary ansatz; RGE alone does not impose that ansatz"}


def separation_checks() -> dict:
    u, v, w, X, a = sp.symbols("u v w X a")
    c = sp.symbols("c0:9")
    phi1 = sum(c[j]*u**j for j in range(3))
    phi2 = sum(c[3+j]*v**j for j in range(3))
    sigma = sum(c[6+j]*w**j for j in range(3))
    equation = sp.Poly(-sp.diff(phi1, u)-sp.diff(phi2, v)+sp.diff(sigma, w), u, v, w)
    solution = sp.solve(equation.coeffs()+[c[1]-c[4]], c, dict=True)[0]
    zero("separation forbids u squared", solution[c[2]])
    zero("separation forbids v squared", solution[c[5]])
    zero("separation forbids w squared", solution[c[8]])
    gamma1, gamma2, r0 = sp.symbols("gamma1 gamma2 r0")
    scalar = sp.series((1+a*r0)*sp.exp(X*(a*gamma1+a*a*gamma2)), a, 0, 3).removeO().expand()
    zero("scalar rapidity boundary quadratic coefficient", scalar.coeff(a, 2).coeff(X, 2)-gamma1**2/2)
    zero("scalar zero gamma1 implies zero quadratic boundary", scalar.subs(gamma1, 0).coeff(a, 2).coeff(X, 2))
    # Broad hard-product factorisation permits a nonseparable extra factor.
    lam = sp.Symbol("lambda")
    broad_log = phi1+phi2+sigma+a*a*lam*(u+v+2*w)**2
    zero("joint extra boundary term is regulator independent",
         (-sp.diff(broad_log, u)-sp.diff(broad_log, v)+sp.diff(broad_log, w)
          ).subs(solution))
    zero("joint extra term is not separately factorised",
         sp.diff(broad_log, u, v)-2*a*a*lam)
    return {"separation_solution": {str(k): str(vv) for k, vv in solution.items()},
            "counterexample": "exp[a_s^2*lambda*(u+v+2w)^2] is boost and rapidity-regulator independent but not J1(u)J2(v)S(w)",
            "conclusion": "the separation proof is valid under J1 J2 S closure and argument locality, not from A=H1 H2 K alone"}


def pinch_checks() -> dict:
    x, y, A, B, q = sp.symbols("x y A B q", real=True)
    p1, p2 = sp.Matrix([1, 0, 0, 1])/2, sp.Matrix([1, 0, 0, -1])/2
    metric = sp.diag(1, -1, -1, -1)
    dot = lambda u, v: (u.T*metric*v)[0]
    k, Q1 = x*p1, x*p1+y*p2
    lines = [k, k-p1, k-Q1, k+p2-Q1]
    weights = [A*(1-x), A*x, B*(1-y), B*y]
    for i, line in enumerate(lines):
        zero(f"crossed-box active line {i} is on shell", dot(line, line))
    zero("crossed-box stationary Landau equation",
         sum((w*line for w, line in zip(weights, lines)), sp.zeros(4, 1)))
    Q1 += sp.Matrix([0, q, 0, 0])
    Q2 = p1+p2-Q1
    m1, m2 = dot(Q1, Q1), dot(Q2, Q2)
    t, u = dot(p1-Q1, p1-Q1), dot(p2-Q1, p2-Q1)
    zero("DPS determinant invariant", t*u-m1*m2-q*q)
    gram = sp.Matrix([[dot(a, b) for b in [p1, p2, Q1]] for a in [p1, p2, Q1]])
    zero("external Gram determinant", gram.det()-q*q/4)
    cs = sp.symbols("c0:8")
    ws = sp.symbols("alpha0:4", nonnegative=True)
    rows = [sp.Matrix([cs[2*i], cs[2*i+1]]) for i in range(4)]
    mat = sum((w*row*row.T for w, row in zip(ws, rows)), sp.zeros(2))
    squares = sum(ws[i]*ws[j]*(rows[i][0]*rows[j][1]-rows[j][0]*rows[i][1])**2
                  for i in range(4) for j in range(i+1, 4))
    zero("two-loop positive-parameter Cauchy-Binet identity", mat.det()-squares)
    # Relative order a^2 in a complete local MATRIX rapidity evolution:
    # U=I+a X Gamma1+a^2 (X Gamma2+X^2 Gamma1^2/2).
    # Gamma1=0 as an operator, not only after projection, is sufficient.
    X, a = sp.symbols("X a")
    entries = sp.symbols("g0:8")
    g1, g2 = sp.Matrix(2, 2, entries[:4]), sp.Matrix(2, 2, entries[4:])
    evolution = sp.eye(2)+a*X*g1+a*a*(X*g2+X*X*g1*g1/2)
    quadratic = evolution.applyfunc(lambda z: sp.expand(z).coeff(a, 2).coeff(X, 2))
    zero("matrix rapidity evolution quadratic boundary", quadratic-g1*g1/2)
    zero("operator-zero one-loop generator suffices", quadratic.subs(dict.fromkeys(entries[:4], 0)))
    return {"positive_stationary_parameters": [str(z) for z in weights],
            "domain": "0<x,y<1 and A,B>0", "cycle_identity": "det(C^T diag(alpha) C)=sum alpha_e alpha_f minor(C_e,C_f)^2",
            "matrix_sufficient_condition": "local rapidity evolution in a closed operator basis; Gamma1(mu_b)=0 as an operator; X-independent finite matching",
            "not_proved": "classification of contracted Landau graphs, all field assignments, or local rapidity closure"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output must be new; preserve the previous validation record")
    mp.mp.dps = 50
    result = {"baseline_commit": BASELINE,
              "scope": "static overlap; common transverse subtraction; colour algebra; conditional RG/boundary analysis",
              "not_established": ["complete pinched-surface classification", "complete regulated DPS operator matching", "three-loop factorisation or rapidity closure"],
              "overlap": overlap_checks(), "transverse": transverse_checks(),
              "colour": colour_checks(), "evolution": evolution_checks(),
              "separation": separation_checks(), "pinch_and_matrix": pinch_checks()}
    result["checks"] = CHECKS
    result["check_count"] = len(CHECKS)
    result["runtime"] = {"python": platform.python_version(), "sympy": sp.__version__,
                         "mpmath": mp.__version__, "numpy": np.__version__, "mp_dps": mp.mp.dps}
    result["producer_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result["passed"] = True
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(f"PASS: {len(CHECKS)} bounded identity/numerical checks. See scope and exclusions in {args.output}.")


if __name__ == "__main__":
    main()
