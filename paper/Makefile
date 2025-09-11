SHELL := /bin/bash

.PHONY: pdf fast clean distclean open

MAIN := main.tex

# Full build with bibtex, multiple passes, using latexmk
pdf:
	latexmk -pdf -interaction=nonstopmode -halt-on-error $(MAIN)
	@echo "✅ PDF ready: build/main.pdf"

# Quick pass without bibtex (useful while drafting)
fast:
	pdflatex -interaction=nonstopmode -output-directory=build $(MAIN) || true

open:
	@if command -v open >/dev/null 2>&1; then open build/main.pdf; \
	elif command -v xdg-open >/dev/null 2>&1; then xdg-open build/main.pdf; \
	else echo "👉 Open build/main.pdf manually"; fi

clean:
	latexmk -C
	rm -rf build

distclean: clean
	rm -f *.synctex.gz
