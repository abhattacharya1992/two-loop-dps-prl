# Paper 2 package. Every target runs from this directory.
#
#   make test        fast tests, no third-party input needed (about a minute)
#   make fetch       download the two ancillary files of arXiv:1809.05388 and check them
#   make test-full   the tests plus the published expressions themselves (after make fetch;
#                    2-3 minutes more)
#   make analyse     rerun the gg -> HH analyses on the stored points (seconds), including the choice
#                    of the fit model by the window shift and the least-squares fits of App. B
#   make checks      rerun the one-loop, subtraction, RG and orientation checks (about a minute)
#   make figures     tables and Figs. 2-5 of the paper from the analysis results (written
#                    into ../manuscript_prl, which also uses the gg -> ZZ results in ../zz)
#   make li22        Li_{2,2} of evaluate_hh.py against an independent representation, for every
#                    argument of the scan, its precision control and the blind test of ../r0
#                    (after make fetch; logs/check_li22_values.log, about 11 CPU hours, run on
#                    WORKERS cores)
#   make scan        compute missing scan points (after make fetch; existing points are kept if
#                    they were computed with the present code; WORKERS cores)
#   make rescan      recompute every scan point (after a change of the evaluation code; about 3 CPU
#                    hours, run on WORKERS cores); the blind-test points of ../r0/data are
#                    recomputed with run_hh_new_velocities.py --recompute there
#   make aa          the exact one-loop gg -> AA amplitudes of P1 from 90 degrees to the forward
#                    limit (Fig. 2; about a minute). Needs the AA evaluator in common/aa_source of the master
#                    folder or the one in which paper2/ lives, with its inputs downloaded (make inputs
#                    there); for a P1 package elsewhere, make aa P1=<its top directory>

# A failed scan must propagate through tee to make.
SHELL := /bin/bash
.SHELLFLAGS := -o pipefail -c
PY ?= python3
WORKERS ?= 2
P1 ?=

.PHONY: all test test-full fetch analyse checks figures scan rescan aa li22

all: analyse checks figures test

test:
	$(PY) tests/test_release.py

test-full:
	$(PY) tests/test_release.py --full

fetch:
	$(PY) fetch_inputs.py

analyse:
	cd hh_two_loop && $(PY) analyse_colour.py > logs/analyse_colour.log
	cd hh_two_loop && $(PY) analyse_one_loop.py > logs/analyse_one_loop.log
	cd hh_two_loop && $(PY) analyse_scale.py > logs/analyse_scale.log
	cd hh_two_loop && $(PY) forward_approach.py > logs/forward_approach.log
	cd hh_two_loop && $(PY) analyse_windows.py > logs/analyse_windows.log
	cd hh_two_loop && $(PY) check_fit_models.py > logs/check_fit_models.log

checks:
	mkdir -p one_loop/logs stage2a/logs rg_argument/logs orientation/logs
	cd one_loop && $(PY) verify_dps_one_loop.py > logs/verify_dps_one_loop.log
	cd one_loop && $(PY) hh_spin2_numerator.py > logs/hh_spin2_numerator.log
	cd one_loop && $(PY) landau_boxes.py > logs/landau_boxes.log
	cd stage2a && $(PY) stage2a_bookkeeping.py > logs/stage2a_bookkeeping.log
	cd rg_argument && $(PY) src/verify_integrated_setup.py > logs/verify_integrated_setup.log
	cd orientation && $(PY) orientation_check.py > logs/orientation_check.log

figures:
	$(MAKE) -C ../manuscript_prl tables figure PY="$(PY)"

scan:
	cd hh_two_loop && $(PY) scan_colour_scale.py --workers $(WORKERS) | tee -a logs/scan_rerun.log

rescan:
	cd hh_two_loop && $(PY) scan_colour_scale.py --recompute --workers $(WORKERS) | tee -a logs/scan_rerun.log

li22:
	cd hh_two_loop && $(PY) check_li22_values.py --workers $(WORKERS) | tee logs/check_li22_values.log

aa:
	mkdir -p one_loop/logs
	cd one_loop && $(PY) aa_full_range.py $(abspath $(P1)) > logs/aa_full_range.log
