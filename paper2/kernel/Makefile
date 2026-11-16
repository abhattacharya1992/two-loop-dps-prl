# The boundary term of the one-loop DPS kernel (Paper 2, appendix E): r1 = r2 = 0.  Every target runs
# from this directory; Python >= 3.10 with numpy, mpmath and sympy, no input files.
#
#   make test      recompute results_stage_b.json (the numbers of appendix E) and check that the stored
#                  file is reproduced and that the statements of the appendix hold (about 10 s)
#   make results   rewrite results_stage_b.json
#   make checks    rerun every check script into its .log (about 5 minutes): the one-loop soft
#                  function, the off-shell triangle, the residues of the DPS loop, the ultrasoft
#                  cancellation, the bookkeeping of the region scales, the explicit two-loop graph and
#                  the two evaluators of it against each other

PY ?= python3

.PHONY: test results checks

test:
	$(PY) test_kernel.py
	$(PY) test_model_integral.py
	$(PY) model_fast.py

results:
	$(PY) stage_b_results.py

checks:
	$(PY) soft_function.py > soft_function.log
	$(PY) lab_triangle.py > lab_triangle.log
	$(PY) model_kres.py > model_kres.log
	$(PY) us_telescoping.py > us_telescoping.log
	$(PY) kernel_logs.py > kernel_logs.log
	$(PY) model_numerics.py > model_numerics.log
	$(PY) model_ref_check.py > model_ref_check.log
	@grep -l "passed" *.log
