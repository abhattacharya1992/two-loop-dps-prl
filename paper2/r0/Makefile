# The constant r0 of the one-loop DPS kernel (Paper 2, Sec. 7, appendix F, Table 9).  Every target runs from
# this directory; Python >= 3.10 with numpy, scipy, mpmath and sympy, no input files.
#
#   make test      data/: rerun the extraction of r0 from the stored gg -> HH points (../hh) and gg -> ZZ fits
#                  (../zz) and the blind tests, and check the stored results (a few seconds);
#                  theory/: recombine the stored Laurent coefficients of the collinear functions into r0 and check
#                  the stored results and the statements of the paper (a few seconds)
#   make theory    recompute the collinear functions for g -> q qbar and g -> gg from scratch (about 2.5 hours on
#                  one core; see theory/Makefile)

PY ?= python3

.PHONY: test theory

test:
	cd data && $(PY) test_data.py
	$(MAKE) -C theory test

theory:
	cd theory && $(MAKE) quarks gluons
