# One code for the DPS analysis of gg -> HH, gg -> AA and gg -> ZZ (Paper 2).
#
#   make test        the shared code against the stored results of hh/, zz/ and r0/ (seconds)
#   make aa-check    digit check of the gg -> AA evaluation at one point (about 10 minutes)
#   make aa-scan     the gg -> AA scan: 48 points + 2 control points, resumable, WORKERS processes
#                    (about 6 CPU hours; needs the P1 package with its inputs, see README.md)
#   make aa-blind    the two velocities not used to find the closed forms of r_0 (about 4 CPU hours)
#   make analyse     fits and comparisons for every process with points: results_*.json
#
# The published AA evaluator is under <master>/common/aa_source (legacy: paper1/) or the one in
# which paper2/ sits; another one: make P1_ROOT=/path/to/P1 ...
PY ?= python3
WORKERS ?= 2
export P1_ROOT

.PHONY: test aa-check aa-scan aa-blind analyse

test:
	$(PY) tests/test_dps.py
	$(PY) tests/test_provenance.py
	$(PY) tests/test_numerical_guards.py

aa-check:
	$(PY) -m dps.backends.aa --check --gauge 1 --beta 1/2 --k 4

aa-scan:
	$(PY) -m dps.backends.aa --jobs --workers $(WORKERS)

aa-blind:
	$(PY) -m dps.backends.aa --blind --workers $(WORKERS)

analyse:
	$(PY) -m dps.report
