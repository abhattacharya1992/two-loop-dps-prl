PY ?= python3
.PHONY: test check
# A fresh polynomial check, not a full amplitude or momentum-region rescan.
test:
	$(PY) test_checks.py
check: test
	$(PY) param_space.py
	$(PY) exact_checks.py
