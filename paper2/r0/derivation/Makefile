# Analytic derivation of r0: run all checks, or regenerate the files in outputs/.
#   make -C paper2/r0/derivation test
#   make -C paper2/r0/derivation outputs
PYTHON ?= python3
ROOT   := ../../..
CHECKS := check_qq_W1q_W2qb check_qq_W1qb_W2q check_qq_Vqq check_qq_Vggq_W1g_W2g check_qq_sum check_gg

.PHONY: test outputs

test:
	@for c in $(CHECKS); do echo "== $$c.py"; $(PYTHON) -B $$c.py $(ROOT) || exit 1; done
	@echo "== all scripts passed"

outputs:
	@mkdir -p outputs
	@for c in $(CHECKS); do echo "== $$c.py"; $(PYTHON) -B $$c.py $(ROOT) > outputs/$$c.txt || exit 1; done
