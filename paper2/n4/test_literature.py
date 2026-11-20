#!/usr/bin/env python3
"""Compare seven stored N=4 constants with independent published expressions.

Source: Caron-Huot et al., arXiv:1903.10890v2, equations (2.20), (5.4),
(5.10), and (5.11). The published input is transcribed below, without using
the repository's MHV coefficient-fitting code. In the repository's 2016
BDS-like convention, a = 2 g^2, gamma_K = 2 Gamma_cusp,
B = -Gamma_virt, and 1 + sum(kappa_n a^n) = exp(zeta_2 gamma_K/4 + 2 Gamma_3).
The additional cosmic-normalization factor of the 2019 paper is removed.

This script only reads repository data. Optional --output writes an audit
JSON result to the caller's chosen location.
"""

import argparse
import json
from pathlib import Path

import sympy as sp


def compare(repo):
    a = sp.Symbol("a")
    z = {k: sp.Symbol(f"zeta_{k}") for k in range(2, 9)}
    R = sp.Rational
    # Published cusp, converted to gamma_K = 2 Gamma_cusp and a = 2 g^2.
    gamma_k = (
        4*a - 4*z[2]*a**2 + 22*z[4]*a**3
        + (-R(219, 2)*z[6] - 4*z[3]**2)*a**4
    )
    # Gamma_virt: Eq. (5.10); Gamma_3: Eq. (5.11), with g^(2n)=(a/2)^n.
    gamma_virt = (
        -12*z[3]*(a/2)**2
        + (80*z[5] + 16*z[2]*z[3])*(a/2)**3
        - (700*z[7] + 80*z[2]*z[5] + 168*z[3]*z[4])*(a/2)**4
    )
    gamma_3 = (
        R(7, 2)*z[4]*(a/2)**2
        - (R(115, 3)*z[6] + R(56, 3)*z[3]**2)*(a/2)**3
        + (R(1701, 4)*z[8] + 260*z[3]*z[5]
           + 28*z[2]*z[3]**2)*(a/2)**4
    )
    prefactor = sp.series(
        sp.exp(z[2]*gamma_k/4 + 2*gamma_3), a, 0, 5
    ).removeO().expand()
    data_path = Path(repo)/"paper2/n4/results_stage_d.json"
    stored = json.loads(data_path.read_text())["solved_constants"]
    symbols = {str(v): v for v in z.values()}
    even_zeta = {z[k]: sp.zeta(k) for k in (2, 4, 6, 8)}
    rows = []
    for n, key in ((1, "r0"), (2, "r0b"), (3, "r0c"), (4, "r0d")):
        target = sp.sympify(stored[key], locals=symbols)
        difference = sp.simplify(
            (prefactor.coeff(a, n) - target).subs(even_zeta)
        )
        rows.append({
            "constant": f"kappa_{n}", "repository_key": key,
            "repository_expression": stored[key],
            "symbolic_difference": str(difference), "passed": difference == 0,
        })
    for n in (2, 3, 4):
        key = f"B{n}"
        target = sp.sympify(stored[key], locals=symbols)
        difference = sp.simplify(
            (-gamma_virt.coeff(a, n) - target).subs(even_zeta)
        )
        rows.append({
            "constant": f"B_{n}", "repository_key": key,
            "repository_expression": stored[key],
            "symbolic_difference": str(difference), "passed": difference == 0,
        })
    return {
        "source": "https://arxiv.org/pdf/1903.10890v2",
        "source_equations": ["2.20", "5.4", "5.10", "5.11"],
        "independence_note": (
            "The B=-Gamma_virt relation is already a published alternative "
            "representation of the same self-crossing singular terms, tested "
            "against Dixon-Esterlis through five loops in that publication. "
            "Matching B2/B3 is a literature consistency check, not two "
            "new independent self-crossing predictions."
        ),
        "doi": "10.1007/JHEP08(2019)016",
        "normalization": {
            "coupling": "a = 2 g^2",
            "cusp": "gamma_K = 2 Gamma_cusp",
            "linear_exponent": "B = -Gamma_virt",
            "boundary_prefactor": "exp(zeta_2 gamma_K/4 + 2 Gamma_3)",
            "amplitude": "2016 BDS-like normalization; 2019 cosmic factor removed",
        },
        "repository_data": "paper2/n4/results_stage_d.json",
        "sympy_version": sp.__version__,
        "comparisons": rows,
        "passed": all(row["passed"] for row in rows),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo", type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = compare(args.repo)
    rendered = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered)
    for row in result["comparisons"]:
        status = "PASS" if row["passed"] else "FAIL"
        print(f'{status} {row["constant"]}: symbolic difference = {row["symbolic_difference"]}')
    print(f'Independent literature comparisons: {sum(r["passed"] for r in result["comparisons"])}/7 passed')
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
