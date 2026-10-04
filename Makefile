LATEXMK ?= latexmk
PYTHON ?= python3

.PHONY: all details check
.NOTPARALLEL:
all: details
details:
	mkdir -p build
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error -cd -outdir=../build tex/compton_calculation_details.tex
check:
	$(PYTHON) scripts/run_compton_checks.py

.PHONY: check-list
check-list:
	$(PYTHON) scripts/run_compton_checks.py --list
