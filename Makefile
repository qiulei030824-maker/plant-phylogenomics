# Plant Phylogenomics Makefile

.PHONY: help install test clean lint check

help:
	@echo "Targets:"
	@echo "  install   Install Python dependencies"
	@echo "  test      Run tests"
	@echo "  clean     Remove __pycache__ and .pyc files"
	@echo "  lint      Run basic syntax check"
	@echo "  check     Check external tool dependencies"

install:
	pip install -r requirements.txt
	pip install -e .

test:
	python -m pytest tests/ -v

clean:
	find . -type d -name '__pycache__' -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name '*.pyc' -delete

lint:
	python -m py_compile scripts/01_download_pfam_hmm.py
	python -m py_compile config/species_config.py

check:
	@echo "=== Tool availability ==="
	@which hmmsearch 2>/dev/null && echo "hmmsearch: OK" || echo "hmmsearch: NOT FOUND"
	@which mafft 2>/dev/null && echo "mafft: OK" || echo "mafft: NOT FOUND"
	@which raxml-ng 2>/dev/null && echo "raxml-ng: OK" || echo "raxml-ng: NOT FOUND"
	@which FastTree 2>/dev/null && echo "FastTree: OK" || echo "FastTree: NOT FOUND"
	@which trimal 2>/dev/null && echo "trimal: OK" || echo "trimal: NOT FOUND"
