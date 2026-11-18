# One PRL Letter and one Supplemental Material; tests run sequentially.
SHELL := /bin/bash
.SHELLFLAGS := -e -o pipefail -c
PY ?= python3
BENCH_RUNS ?= 3
.NOTPARALLEL:
.PHONY: cache-archive-check phase2-check help threshold independent-checks test-science test pinch-checks pinch-figures manuscript-tests derivation-tests build archive reproduce release-check release-check-quick release-check-tests

help:
	@printf '%s\n' 'make threshold: earlier AA threshold coefficients vs PRL at x=1/2' \
	 'make test-science: scientific tests without TeX figure compilation' \
	 'make test: scientific tests plus figure generation checks' \
	 'make build: regenerate tables, figures, Letter and Supplement' \
	 'make archive: build and validate isolated arXiv source package' \
	 'make reproduce: test, analytic derivations, PDFs and source archive' \
	 'make release-check: compare this branch against frozen main and benchmark'

threshold:
	$(PY) checks/threshold_ggAA/check.py

independent-checks:
	@set -eu; d=$$(mktemp -d); trap 'rm -rf "$$d"' EXIT; \
	$(PY) checks/stage1_physics_20261113/verify_stage1.py --output "$$d/stage1.json"; \
	$(PY) checks/stage1_physics_20261113/test_stage1_guards.py --output "$$d/stage1-guards.json"; \
	$(PY) checks/stage2_novelty_20261113/verify_stage2.py --output "$$d/stage2.json"; \
	$(PY) checks/stage2_novelty_20261113/test_stage2_guards.py --output "$$d/stage2-guards.json"

test: pinch-figures manuscript-tests test-science

test-science: threshold independent-checks pinch-checks dps-links
	$(PY) checks/blind_point_provenance.py
	$(PY) checks/test_boundary_coefficients.py
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

derivation-tests:
	$(MAKE) -C paper2/r0/derivation test PYTHON="$(PY)"

reproduce: test derivation-tests archive

pinch-checks:
	$(MAKE) -C checks/pinch_regions check PY="$(PY)"

manuscript-tests:
	$(MAKE) -C paper2/manuscript_prl tables PY="$(PY)"
	$(MAKE) -C paper2/manuscript_prl pdf PY="$(PY)"
	$(PY) -m unittest discover -s paper2/manuscript_prl -p 'test_make_tables.py' -v
	$(MAKE) -C paper2/manuscript_prl test PY="$(PY)"

pinch-figures:
	$(MAKE) -C paper2/manuscript_prl figure PY="$(PY)"

# Run both revisions in detached worktrees. Does not update a branch.
release-check-tests:
	$(PY) -m unittest discover -s checks -p 'test_release_compare.py' -v

release-check-quick: release-check-tests
	$(PY) checks/release_compare.py --quick --python "$(PY)" --output release-check

release-check: release-check-tests
	$(PY) checks/release_compare.py --python "$(PY)" --output release-check --benchmark-runs "$(BENCH_RUNS)"

.PHONY: dps-links
dps-links:
	$(PY) checks/dps_links/check_couplings.py
	$(PY) checks/dps_links/check_dpd_kernels.py

# Cache snapshots are optional; this check only reads the frozen v2 Git history.
cache-archive-check:
	$(PY) -m unittest discover -s checks -p 'test_region_cache_archive.py' -v
	$(PY) checks/region_cache_archive.py --verify
