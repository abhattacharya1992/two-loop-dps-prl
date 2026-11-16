# gg -> ZZ test of Paper 2, Sec. 6. Every target runs from this directory.
#
#   make test      re-analyse the stored points and logs and compare with results_zz.json,
#                  results_windows.json, results_forward.json, results_wide.json, summary.json
#                  and summary.md
#   make analyse   fits and tables from the stored points and logs: results_zz.json,
#                  results_windows.json (fit stability), results_forward.json (approach of the
#                  one-loop amplitude to its forward limit), results_wide.json (the one-loop
#                  amplitudes from 90 degrees to the forward limit), summary.md, summary.json,
#                  results_fit_models.json (the fit model chosen by the window shift, without the
#                  predictions, and least-squares fits)
#   make inputs    check the VVamp files placed in inputs/ and unpack them into build/
#   make build     translate the ggvvamp code into Python (build/ggv, about 30 s)
#   make checks    tensor basis, Li22 argument order, p0002 statement by statement, and the
#                  comparison with the FORM file at the benchmark point of ggvvamp (one loop, and
#                  the two-loop C_F and C_A parts of A1, A2, A3, A10, A17); writes logs/ (about
#                  2 minutes)
#   make routeb    the FORM file on the forward line at beta = 1/2: the C_F part at delta =
#                  10^-4 and the C_A part at 10^-3, 10^-4, 10^-5 (points_routeB/, about 10 min)
#   make cxx       compile the generated C++ files of ggvvamp (algebraic part, quad precision,
#                  needs g++ and libquadmath) and compare them with the Python translation at the
#                  benchmark point and at beta = 1/2, delta = 10^-4 (logs/check_cxx.log, about
#                  15 min, most of it compiling; the objects stay in build/cxx)
#   make native-parser-test  compile and test the independent FORM reader; needs g++,
#                  GiNaC and CLN development files, but no published amplitude inputs
#   make li22      Li22 of the translation against an independent integral representation, for
#                  all arguments at the benchmark point and at delta = 10^-12 for the three
#                  velocities (logs/check_li22_values.log, about 45 min)
#   make li22-scan the integral branch of that Li22 at the working precision of every scan point,
#                  its control and the blind test of ../r0, against the integral split at points
#                  accumulating at t = 1 (logs/check_li22_scan.log, a few CPU hours on WORKERS cores)
#   make point BETA=1/2 K=6    compute one forward point (2-10 minutes), 30 + 5K digits
#   make scan      compute all missing points (about 2.5 CPU hours; existing points are kept)
#   make wide      compute the points away from the forward limit, from 90 degrees down to
#                  sin^2(theta) = 10^-3.5, next to the first scan point (points_wide/, about a CPU
#                  hour; existing points are kept)
#
# inputs, build, checks, routeb, cxx, li22, li22-scan, point, scan and wide need the VVamp files (see
# README.md, Inputs). After checks, routeb, cxx or li22, run make analyse to update summary.md and
# summary.json.

# Stop a multi-point recipe at the first failed evaluation, including inside loops.
SHELL := /bin/bash
.SHELLFLAGS := -e -o pipefail -c
PY ?= python3
CXX ?= g++
WORKERS ?= 2
BETA ?= 1/2
K ?= 6
DPS = $(shell echo $$((30 + 5 * $(K))))

.PHONY: test analyse inputs build checks routeb cxx native-parser-test li22 li22-scan point scan wide

test:
	$(PY) scripts/test_cache.py
	$(PY) scripts/test_check_cxx.py
	$(PY) scripts/test_make.py
	$(PY) scripts/test_numerical_checks.py
	$(PY) scripts/test_p0002.py
	$(PY) scripts/test_master_families.py
	$(PY) scripts/test_zz.py points points_routeB logs results_zz.json summary.md

analyse:
	$(PY) scripts/analyse_zz.py points results_zz.json > logs/analyse_zz.log
	$(PY) scripts/make_summary.py results_zz.json points points_routeB logs summary.md > /dev/null
	$(PY) scripts/analyse_windows.py points results_windows.json > logs/analyse_windows.log
	$(PY) scripts/check_fit_models.py points results_fit_models.json > logs/check_fit_models.log
	$(PY) scripts/forward_curve.py points results_zz.json results_forward.json > logs/forward_curve.log
	$(PY) scripts/wide_curve.py points points_wide results_wide.json > logs/wide_curve.log
	@echo "results_zz.json, results_windows.json, results_fit_models.json, results_forward.json, results_wide.json, summary.md and summary.json written"

inputs:
	cd inputs && sha256sum --check --ignore-missing SHA256SUMS
	mkdir -p build/inputs
	tar -xzf inputs/ggvvamp-1.0.tar.gz -C build
	cp inputs/gg_projectors-1.0.inc build/inputs/
	gunzip -c inputs/gg_Aj_A_qT-1.0.inc.gz > build/inputs/gg_Aj_A_qT-1.0.inc
	gunzip -c inputs/vvpmasters_opt-1.0.m.gz > build/inputs/vvpmasters_opt-1.0.m

build:
	$(PY) scripts/translate_ggvvamp.py build/ggvvamp-1.0 build/ggv

checks:
	$(PY) scripts/check_tensors.py build/inputs > logs/check_tensors.log
	$(PY) scripts/check_li22_order.py build/ggv > logs/check_li22_order.log
	$(PY) scripts/check_p0002.py build/ggvvamp-1.0 build/ggv > logs/check_p0002.log
	$(PY) scripts/crosscheck_form.py build/ggv build/inputs 15625 -2384.89727896833 6467.215561 625 30 1l \
	  > logs/crosscheck_1l_benchmark.log
	$(PY) scripts/crosscheck_form.py build/ggv build/inputs 15625 -2384.89727896833 6467.215561 625 30 1 3 10 17 \
	  > logs/crosscheck_2l_benchmark.log
	$(PY) scripts/crosscheck_form.py build/ggv build/inputs 15625 -6147.88716003167 6467.215561 625 30 2 \
	  > logs/crosscheck_2l_benchmark_x12.log
	@tail -qn1 logs/check_tensors.log logs/check_p0002.log logs/crosscheck_1l_benchmark.log
	@cat logs/check_li22_order.log
	@grep -h "2l_C" logs/crosscheck_2l_benchmark*.log

cxx:
	@cxx_log=$$(mktemp logs/check_cxx.XXXXXX); trap 'rm -f "$$cxx_log"' EXIT HUP INT TERM; \
	if $(PY) scripts/check_cxx.py build/ggv build/ggvvamp-1.0 build/cxx bench 1/2:4 > "$$cxx_log"; then \
	  mv "$$cxx_log" logs/check_cxx.log; cat logs/check_cxx.log; \
	else status=$$?; cat "$$cxx_log" >&2; exit $$status; fi

native-parser-test:
	@validation_dir=$$(mktemp -d); trap 'rm -rf "$$validation_dir"' EXIT; \
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -O2 -std=c++17 ../../checks/native_form.cpp \
	  $(LDFLAGS) -lginac -lcln -o "$$validation_dir/native_form"; \
	$(PY) ../../checks/test_native_form.py "$$validation_dir/native_form"

li22:
	$(PY) scripts/check_li22_values.py build/ggv bench 1/4:12 1/2:12 2/3:12 > logs/check_li22_values.log
	@cat logs/check_li22_values.log

li22-scan:
	$(PY) scripts/check_li22_scan.py build/ggv --workers $(WORKERS) | tee logs/check_li22_scan.log

routeb:
	for k in 3 4 5; do $(PY) scripts/zz_routeB.py build/inputs build/routeB 1/2 $$k $$((20 + 5 * k)) 2l_CA; done
	$(PY) scripts/zz_routeB.py build/inputs build/routeB_CF 1/2 4 40 2l_CF
	cp build/routeB/zzB_b1o2_d*.json points_routeB/
	cp build/routeB_CF/zzB_b1o2_d4_p40.json points_routeB/zzB_CF_b1o2_d4_p40.json

point:
	$(PY) scripts/zz_points.py build/ggv points $(BETA) $(K) $(DPS)

scan:
	for b in 1/2 1/4 2/3; do for k in 4 5 6 7 8 9 10 11 12; do \
	  $(PY) scripts/zz_points.py build/ggv points $$b $$k $$((30 + 5 * k)); done; done
	$(PY) scripts/zz_points.py build/ggv points 1/2 10 95

wide:
	for b in 1/4 1/2 2/3; do $(PY) scripts/zz_wide.py build/ggv points_wide $$b 0 1 2 3 4 5 6 7; done
