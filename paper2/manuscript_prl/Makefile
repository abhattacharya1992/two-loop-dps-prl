# Builds the single-paper layout: the Letter (letter.tex, REVTeX 4.2, PRL, with End Matter) and its
# Supplemental Material (supplement.tex), which holds the derivations and the QCD amplitudes. The
# numbers, tables and figures come from the results in ../hh, ../zz, ../n4, ../kernel, ../r0 and
# ../../common/dps through the generators of ../_build. Every target runs from this
# directory. REVTeX 4.2 must be installed (kpsewhich revtex4-2.cls).
#
#   make          write numbers.tex, the tables and sm.bib, draw the figures, build the two PDFs and
#                 check the length of the Letter against the PRL limits (count_words.py)
#   make tables   only numbers.tex, tables/ and sm.bib (make_tables.py)
#   make figure   only the figures (make_figure.py; Fig. 1 of the Letter needs TikZ-Feynman)
#   make pdf      only the LaTeX runs: letter.pdf, then letter_labels.tex (the numbers of the Letter
#                 that the Supplement cites, letter_labels.py), then supplement.pdf
#   make words    only the APS length estimate of the Letter
#   make test     word-counter and ancillary archive regression tests
#   make clean    remove the LaTeX by-products (the .bbl files are kept, since arXiv needs them)

PY ?= python3

.PHONY: all tables figure pdf words test clean

all: tables figure pdf words

tables:
	$(PY) ../r0/data/endpoint_test.py > endpoint_test.log
	$(PY) make_tables.py > tables.log
	@tail -1 tables.log

figure:
	$(PY) make_figure.py

define latex_runs
	@pdflatex -interaction=nonstopmode -halt-on-error $(1).tex > $(1)_latex.log || { tail -30 $(1)_latex.log; exit 1; }
	@bibtex $(1) > $(1)_bibtex.log || { cat $(1)_bibtex.log; exit 1; }
	@for i in 1 2 3 4; do \
	  pdflatex -interaction=nonstopmode -halt-on-error $(1).tex > $(1)_latex.log || { tail -30 $(1)_latex.log; exit 1; }; \
	  if [ $$i -ge 2 ] && ! grep -q "Rerun to get" $(1)_latex.log; then break; fi; \
	done
	@grep -E "LaTeX Warning|Overfull|Underfull" $(1)_latex.log || echo "$(1).pdf built without LaTeX warnings"
	@grep -iE "warning|error" $(1)_bibtex.log | grep -vE "missing number|empty year in SM|There w(as|ere)" || echo "bibtex ($(1)): no warnings"   # apsrev4-2 warns for entries without an issue number, and for the Supplemental Material entry, which has no year
endef

pdf:
	$(call latex_runs,letter)
	$(PY) letter_labels.py
	$(call latex_runs,supplement)

words:
	$(PY) count_words.py

test:
	$(PY) -m unittest discover -s . -p 'test_count_words.py' -v
	$(PY) -m unittest discover -s . -p 'test_make_arxiv.py' -v

clean:
	rm -f *.aux *.log *.out *.blg *Notes.bib
