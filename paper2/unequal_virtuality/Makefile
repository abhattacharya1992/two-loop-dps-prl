PY ?= python3
WORKERS ?= 1
.PHONY: prepare test verify scan analyse control
prepare:
	$(PY) prepare.py
test:
	$(PY) test_unequal_virtuality.py
verify:
	$(PY) verify_stage1.py --retained-records
	$(PY) verify_matching.py --check
scan:
	$(PY) run_scan.py --workers $(WORKERS)
analyse:
	$(PY) analyse_scan.py runs/stage1 --output runs/stage1_analysis.json
control:
	$(PY) run_scan.py --k 10 --dps 95 --output runs/stage1
