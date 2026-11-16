# One PRL Letter and one Supplemental Material; tests run sequentially.
SHELL := /bin/bash
.SHELLFLAGS := -e -o pipefail -c
PY ?= python3
.NOTPARALLEL:
.PHONY: help threshold independent-checks test build archive

help:
	@printf '%s\n' 'make threshold: earlier AA threshold coefficients vs PRL at x=1/2' \
	 'make test: all included stored-data tests and bounded Stage 1/2 checks' \
	 'make build: regenerate tables, figures, Letter and Supplement' \
	 'make archive: build and validate isolated arXiv source package'

threshold:
	$(PY) checks/threshold_ggAA/check.py

independent-checks:
	@set -eu; d=$$(mktemp -d); trap 'rm -rf "$$d"' EXIT; \
	$(PY) checks/stage1_physics_20261113/verify_stage1.py --output "$$d/stage1.json"; \
	$(PY) checks/stage1_physics_20261113/test_stage1_guards.py --output "$$d/stage1-guards.json"; \
	$(PY) checks/stage2_novelty_20261113/verify_stage2.py --output "$$d/stage2.json"; \
	$(PY) checks/stage2_novelty_20261113/test_stage2_guards.py --output "$$d/stage2-guards.json"

test: threshold independent-checks pinch-checks
	$(PY) checks/blind_point_provenance.py
	$(PY) checks/test_boundary_coefficients.py
	$(MAKE) -C paper2/manuscript_prl test PY="$(PY)"
	$(MAKE) -C common/dps test PY="$(PY)"
	$(MAKE) -C paper2/hh test PY="$(PY)"
	$(MAKE) -C paper2/zz test PY="$(PY)"
	$(MAKE) -C paper2/n4 test PY="$(PY)"
	$(MAKE) -C paper2/kernel test PY="$(PY)"
	$(MAKE) -C paper2/r0 test PY="$(PY)"
	$(MAKE) -C paper2/unequal_virtuality test PY="$(PY)"
	$(MAKE) -C paper2/unequal_virtuality verify PY="$(PY)"

build:
	$(MAKE) -C paper2/manuscript_prl PY="$(PY)"

archive: build
	$(PY) paper2/manuscript_prl/make_arxiv.py --check

.PHONY: pinch-checks
pinch-checks:
	$(MAKE) -C paper2/manuscript_prl figure PY="$(PY)"
	$(PY) -m unittest discover -s paper2/manuscript_prl -p 'test_make_tables.py' -v
	$(MAKE) -C checks/pinch_regions check PY="$(PY)"
