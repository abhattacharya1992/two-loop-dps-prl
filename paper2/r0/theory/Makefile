# r0 of the one-loop DPS kernel from the regions (Paper 2, Sec. 7.2 and appendix F).  Python >= 3.10 with numpy,
# scipy, mpmath and sympy; no input files.
#
#   make test      recombine the stored Laurent coefficients (gqq_laurent_x*.json, gg_laurent_x*.json) into r0 and
#                  check that r0_theory_quarks.json and r0_theory_gluons.json are reproduced and that the statements
#                  of the paper hold (a few seconds)
#   make results   rewrite r0_theory_quarks.json and r0_theory_gluons.json from the stored Laurent coefficients
#   make quarks    recompute the collinear functions for g -> q qbar from scratch at x = 1/10, 3/10, 1/2, 7/10, 9/10,
#                  with the precision runs at x = 1/2 (about half an hour on one core); the stored Laurent coefficients,
#                  results and sector caches are moved to stored/ first and compared with the new ones at the end
#   make gluons    the same for g -> gg, with the same precision runs (about two hours on one core)
#   make n4        planar N = 4 (Sec. 7.3): the scalar-loop graphs of n4_diagrams.py at the five x (run_n4.py,
#                  assemble_n4.py, about ten minutes on one core), the scheme check at x = 1/2 (the gluon-loop graphs
#                  with g^mu_mu = 4, check_fdh.py), the precision runs, and kappa_1 from the regions (combine_n4.py,
#                  n4_hexagon.py -> kappa1_n4.json); the stored files are moved to stored/ first
#   make threeloop the three-loop DPS logarithms from the renormalisation group (Paper II, Sec. 7.2): three_loop.py -> three_loop.json
#                  (seconds; sympy); make test reruns it and compares with the stored file
#   make glauber   the imaginary part of the kernel (task 5, workbook VII.15): the Coulomb phase of the soft function and
#                  the Glauber region, glauber.py -> glauber.json (about ten seconds; no amplitude evaluations; reads the
#                  stored fits of r0/data and common/dps, the model graph of ../../kernel and n4/de_data.json); part 2b
#                  runs glauber_regions.py, the four eta-regulated regions of the off-shell exchange in closed form
#   make checks    the pipeline checks: sector decomposition against known integrals, the normal form against the
#                  integrands, charge conjugation, sectors against direct quadrature, the origin of the
#                  helicity-dependent term of r0 for g -> gg (helicity_terms.json), the analytic steps of
#                  Sec. 7 and appendix F as the workbook (Part VII) states them, and the one-loop amplitudes of
#                  planar N = 4 from the box functions (check_bdslike.py) (about three minutes)
#
#   make test also reruns glauber.py and compares with the stored glauber.json.

SHELL := /bin/bash
.SHELLFLAGS := -e -o pipefail -c
PY ?= python3
XS = 3/10 7/10 1/10 9/10 1/2

.PHONY: test results quarks gluons n4 glauber threeloop checks

test:
	$(PY) test_theory.py
	$(PY) test_cache.py
	$(PY) test_laurent_inputs.py
	$(PY) test_numerical_checks.py
	$(PY) test_colour.py

glauber:
	$(PY) glauber.py | tee glauber.log

threeloop:
	$(PY) three_loop.py | tee three_loop.log

results:
	$(PY) r0_theory_quarks.py
	$(PY) r0_theory_gluons.py

quarks:
	$(PY) archive_run.py gqq
	for x in $(XS); do t=$$(echo $$x | tr '/' 'o'); \
	  $(PY) run_gqq.py $$x > run_gqq_x$$t.log 2>&1 && $(PY) assemble_gqq.py $$x > assemble_x$$t.log 2>&1 || exit 1; done
	NQ=96,72,40 $(PY) assemble_gqq.py 1/2 > assemble_x1o2_nq96.log 2>&1
	NEPS=16 $(PY) assemble_gqq.py 1/2 > assemble_x1o2_neps16.log 2>&1
	$(PY) r0_theory_quarks.py
	$(PY) compare_stored.py gqq

gluons:
	$(PY) archive_run.py gg
	for x in $(XS); do t=$$(echo $$x | tr '/' 'o'); \
	  $(PY) run_gg.py $$x > run_gg_x$$t.log 2>&1 && $(PY) assemble_gg.py $$x > assemble_gg_x$$t.log 2>&1 || exit 1; done
	NQ=96,72,40 $(PY) assemble_gg.py 1/2 > assemble_gg_x1o2_nq96.log 2>&1
	NEPS=16 $(PY) assemble_gg.py 1/2 > assemble_gg_x1o2_neps16.log 2>&1
	$(PY) r0_theory_gluons.py
	$(PY) compare_stored.py gg

n4:
	$(PY) archive_run.py n4
	for x in $(XS); do t=$$(echo $$x | tr '/' 'o'); \
	  $(PY) run_n4.py $$x > run_n4_x$$t.log 2>&1 && $(PY) assemble_n4.py $$x > assemble_n4_x$$t.log 2>&1 || exit 1; done
	NQ=96,72,40 $(PY) assemble_n4.py 1/2 > assemble_n4_x1o2_nq96.log 2>&1
	NEPS=16 $(PY) assemble_n4.py 1/2 > assemble_n4_x1o2_neps16.log 2>&1
	$(PY) run_n4.py 1/2 --fdh > run_fdh_x1o2.log 2>&1 && $(PY) assemble_n4.py 1/2 --fdh > assemble_fdh_x1o2.log 2>&1
	$(PY) check_fdh.py 1/2 > check_fdh.log 2>&1
	@tail -n 2 check_fdh.log
	$(PY) n4_hexagon.py
	$(PY) compare_stored.py n4

checks:
	$(PY) test_secdec.py > test_secdec.log 2>&1
	$(PY) test_gg_tree.py > test_gg_tree.log 2>&1
	$(PY) check_steps.py 3/10 Vqq ek -3/10 --allow-legacy-input > check_steps.log 2>&1
	$(PY) check_csym.py Vqq 1/10 -3/10 > check_csym.log 2>&1
	$(PY) check_steps_gg.py 3/10 Vgh1,Vbuba1,Vbubb2,W1g,W1kb,W1k,W2g,W2k,W2kb > check_steps_gg.log 2>&1
	@tail -n 2 check_steps_gg.log
	$(PY) check_helicity_terms.py > check_helicity_terms.log 2>&1
	@tail -n 2 check_helicity_terms.log
	$(PY) check_formulas.py > check_formulas.log 2>&1
	@tail -n 1 check_formulas.log
	$(PY) check_bdslike.py > check_bdslike.log 2>&1
	@tail -n 1 check_bdslike.log
