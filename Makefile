# One PRL Letter and one Supplemental Material; tests run sequentially.
SHELL := /bin/bash
.SHELLFLAGS := -e -o pipefail -c
PY ?= python3
BENCH_RUNS ?= 3
PERF_RUNS ?= 5
PERF_V2_RUNS ?= 3
.NOTPARALLEL:
.PHONY: benchmark-size benchmark-pipeline benchmark-v2 benchmark-tests cache-archive-check phase2-check help threshold independent-checks test-science test pinch-checks pinch-figures manuscript-tests derivation-tests build archive archive-check reproduce release-check release-check-quick release-check-tests

help:
	@printf '%s\n' 'make threshold: earlier AA threshold coefficients vs PRL at x=1/2' \
	 'make test-science: scientific tests without TeX figure compilation' \
	 'make test: scientific tests plus figure generation checks' \
	 'make build: regenerate tables, figures, Letter and Supplement' \
	 'make archive: build and validate isolated arXiv source package' \
	 'make reproduce: test, analytic derivations, PDFs and source archive' \
	 'make release-check: compare publication and calculations against frozen PRL v2' \
	 'make benchmark-size: tracked Git tree size for original, v2, and candidate' \
	 'make benchmark-pipeline: original main versus candidate stored-data runtime' \
	 'make benchmark-v2: frozen v2 versus candidate stored-data runtime'

threshold:
	$(PY) checks/threshold_ggAA/check.py

independent-checks:
	@set -eu; d=$$(mktemp -d); trap 'rm -rf "$$d"' EXIT; \
	$(PY) checks/stage1_physics_20261113/verify_stage1.py --output "$$d/stage1.json"; \
	$(PY) checks/stage1_physics_20261113/test_stage1_guards.py --output "$$d/stage1-guards.json"; \
	$(PY) checks/stage2_novelty_20261113/verify_stage2.py --output "$$d/stage2.json"; \
	$(PY) checks/stage2_novelty_20261113/test_stage2_guards.py --output "$$d/stage2-guards.json"

test: manuscript-tests test-science

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

archive: build archive-check

archive-check:
	$(PY) paper2/manuscript_prl/make_arxiv.py --check

derivation-tests:
	$(MAKE) -C paper2/r0/derivation test PYTHON="$(PY)"

reproduce: test derivation-tests archive-check

pinch-checks:
	$(MAKE) -C checks/pinch_regions check PY="$(PY)"

manuscript-tests:
	$(MAKE) -C paper2/manuscript_prl all PY="$(PY)"
	$(PY) -m unittest discover -s paper2/manuscript_prl -p 'test_make_tables.py' -v
	$(MAKE) -C paper2/manuscript_prl test PY="$(PY)"

pinch-figures:
	$(MAKE) -C paper2/manuscript_prl figure PY="$(PY)"

# Run both revisions in detached worktrees. Does not update a branch.
release-check-tests:
	$(PY) -m unittest discover -s checks -p 'test_release_compare.py' -v

release-check-quick: cache-archive-check release-check-tests
	$(PY) checks/release_compare.py --quick --archived-cache-manifest checks/region_cache_manifest.json --python "$(PY)" --output release-check

release-check: cache-archive-check release-check-tests
	$(PY) checks/release_compare.py --archived-cache-manifest checks/region_cache_manifest.json --python "$(PY)" --output release-check --benchmark-runs "$(BENCH_RUNS)"

.PHONY: dps-links
dps-links:
	$(PY) checks/dps_links/check_couplings.py
	$(PY) checks/dps_links/check_dpd_kernels.py

# Reproduce v2 and the candidate in detached worktrees, then compare their outputs.
phase2-check: release-check

# Cache snapshots are optional; this check only reads the frozen v2 Git history.
cache-archive-check:
	$(PY) -m unittest discover -s checks -p 'test_region_cache_archive.py' -v
	$(PY) checks/region_cache_archive.py --verify

# Optional performance reporting: no scientific input is changed and no ref is moved.
benchmark-tests:
	$(PY) -m unittest discover -s checks -p 'test_benchmark_tools.py' -v

benchmark-size: benchmark-tests
	$(PY) checks/benchmark_size.py --output benchmark-results/size

benchmark-pipeline: benchmark-tests
	$(PY) checks/benchmark_pipeline.py --baseline-label original_main --runs $(PERF_RUNS) --python "$(PY)" --output benchmark-results/pipeline-original

benchmark-v2: benchmark-tests
	$(PY) checks/benchmark_pipeline.py --baseline d362f551357f80125649e583b22b6dbf8aeccb4e --baseline-label frozen_v2 --runs $(PERF_V2_RUNS) --python "$(PY)" --output benchmark-results/pipeline-v2

# Dedicated compression/ module; never part of the default PRL reproduction.
AA_REPO ?= ggAA-threshold
.PHONY: compression-tests compression-points compression-aa-recorded
compression-tests:
	$(PY) -m unittest discover -s compression -p 'test_*.py' -v

compression-points: compression-tests
	$(PY) compression/point_benchmark.py --process all --repeats 5

compression-aa-recorded: compression-tests
	$(PY) compression/aa_threshold.py --aa-repo "$(AA_REPO)"

# These two targets require the published large third-party evaluator inputs.
.PHONY: compression-exact-hh compression-exact-zz
compression-exact-hh:
	$(PY) compression/exact_evaluators.py --process HH --run-exact --beta 1/10 --costheta 0 --digits 50 --repeats 3

compression-exact-zz:
	$(PY) compression/exact_evaluators.py --process ZZ --run-exact --beta 1/10 --costheta 0 --digits 50 --repeats 3

# Two-stage unified HH/AA/ZZ runner: shared evaluation stage and separate report stage.
.PHONY: compression-unified compression-stage1 compression-stage2
COMPRESSION_BETA ?= 1/10
COMPRESSION_COSTHETA ?= 0
COMPRESSION_DIGITS ?= 50
COMPRESSION_REPEATS ?= 3

compression-unified: compression-tests
	$(PY) compression/unified_pipeline.py run --aa-repo "$(AA_REPO)" --beta $(COMPRESSION_BETA) --costheta $(COMPRESSION_COSTHETA) --digits $(COMPRESSION_DIGITS) --repeats $(COMPRESSION_REPEATS)

compression-stage1: compression-tests
	$(PY) compression/unified_pipeline.py stage1 --aa-repo "$(AA_REPO)" --beta $(COMPRESSION_BETA) --costheta $(COMPRESSION_COSTHETA) --digits $(COMPRESSION_DIGITS) --repeats $(COMPRESSION_REPEATS)

compression-stage2: compression-tests
	$(PY) compression/unified_pipeline.py stage2

# Check the opt-in original-amplitude benchmark guards; no expensive input needed.
.PHONY: compression-full-tests
compression-full-tests: compression-tests
	$(PY) -m unittest discover -s compression -p 'test_full_performance.py' -v
