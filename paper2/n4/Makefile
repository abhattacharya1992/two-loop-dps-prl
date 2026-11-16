# The b-space structure of Paper 2, Sec. 5.1, against planar N = 4 super Yang-Mills (the
# self-crossing expansions of Dixon and Esterlis). Every target runs from this directory; only
# sympy is needed, no input files.
#
#   make test      rerun the comparison and check that it reproduces results_stage_d.json, the
#                  file the manuscript reads for Tables 5 and 6, then compare seven constants
#                  with the independent analytic expressions of arXiv:1903.10890 (about 10 s)
#   make analyse   the comparison itself: results_stage_d.json and stage_d.log
#   make blind     the two independent reimplementations, MHV and NMHV (about 6 minutes, most
#                  of it in the sympy series of the Gamma functions); they print one PASS/FAIL
#                  line per coefficient into stage_d_blind.log and nmhv_blind.log

PY ?= python3

.PHONY: test analyse blind primary-inputs primary-check

test:
	$(PY) test_n4.py
	$(PY) test_literature.py
	$(PY) test_failures.py

analyse:
	$(PY) stage_d.py de_data.json results_stage_d.json > stage_d.log
	@grep -c PASS stage_d.log | sed 's/^/coefficients that agree: /'

blind:
	$(PY) stage_d_blind.py de_data.json > stage_d_blind.log
	$(PY) nmhv_blind.py de_data.json > nmhv_blind.log
	@grep -c "PASS" stage_d_blind.log nmhv_blind.log

primary-inputs:
	$(PY) fetch_primary.py

primary-check:
	$(PY) check_primary_nmhv.py
