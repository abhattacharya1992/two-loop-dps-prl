"""Laurent expansion of the scalar-loop graphs of J1 for g -> gg at a fixed x (n4_diagrams.py), graph by graph, b space
at L_b = 0, from the integrands of run_n4.py; as assemble_gg.py (assemble_lib.py).  Colour weight {'ns': .} per complex
adjoint scalar, in units of C_A.

    python3 assemble_n4.py 3/10 [names]        (reads n4_integrands_x3o10.pkl, writes n4_laurent_x3o10.json)
    python3 assemble_n4.py 3/10 --fdh          (the gluon-loop graphs with g^mu_mu = 4: fdh_laurent_x3o10.json)
    NQ=96,72,40 / NEPS=16 as for assemble_gg.py
"""
import sys, pickle, os
import sympy as sp
import gg_diagrams
import n4_diagrams
from assemble_lib import assemble
from cache_utils import arguments, load_integrands, quadrature

args = arguments(__doc__, fdh=True, assembly=True)
fdh = args.fdh
xval = sp.Rational(args.x)
tag = str(xval).replace('/', 'o')
only = args.only
nq, neps, suffix = quadrature()
kind = 'fdh' if fdh else 'n4'
input_path = f'{kind}_integrands_x{tag}.pkl'
data = load_integrands(input_path, kind, xval, allow_legacy=args.allow_legacy_input)
if fdh:
    assemble(data, gg_diagrams.DIAGRAMS, gg_diagrams.colour, f'fdh_laurent_x{tag}{suffix}.json',
             f'secdec_cache_fdh_x{tag}.pkl', only=only, nq=nq, neps=neps,
             input_path=input_path, recompute=args.recompute)
else:
    assemble(data, n4_diagrams.DIAGRAMS, n4_diagrams.colour, f'n4_laurent_x{tag}{suffix}.json',
             f'secdec_cache_n4_x{tag}.pkl', only=only, nq=nq, neps=neps,
             input_path=input_path, recompute=args.recompute)
