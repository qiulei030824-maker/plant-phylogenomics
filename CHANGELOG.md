# Changelog

All notable changes to this project will be documented in this file.

## [0.4.0] - 2026-05-06

### Added
- `scripts/04_build_domain_tree.py`: Build phylogenetic tree of domain-containing proteins
  - Strategy: "all" (default, every sequence) or "longest" (one per species)
  - MAFFT alignment with full parameter control (--mafft-custom, --mafft-args, --mafft-threads)
  - IQ-TREE2 tree inference with full parameter control (--iqtree-custom, --iqtree-args)
  - {input} and {prefix} placeholders for custom command substitution
  - --skip-align / --skip-tree for step-by-step execution
  - --force to re-run existing outputs
  - Full argparse CLI (-i/--input-dir, -o/--output-dir, --align-dir, --tree-dir)
  - Semi-automated: auto-resolves paths from species_config if no explicit paths provided

## [0.3.0] - 2026-05-06

### Added
- `scripts/03_extract_domain_seqs.py`: Extract domain-containing PEP/CDS sequences per species
  - Task 1: Auto-scan species directories for PEP and CDS files
  - Task 2: Parse HMMER domtblout targets, extract matching sequences via seqkit grep
  - Task 3: CDS protein-to-transcript ID mapping (handles _P→_t, -P suffix removal, etc.)
  - Full argparse CLI (positional pfam_id, -i/--input, -o/--output-dir, --skip-cds, --force, --pep-suffix/--cds-suffix)
  - Semi-automated: auto-resolves paths from species_config if no explicit paths provided

## [0.2.0] - 2026-05-06

### Added
- `scripts/02_analyze_domtblout.py`: Analyze HMMER domtblout output
  - Task 1: Count domains per gene -> gene_domain_count.tsv
  - Task 2: Count unique genes per species -> species_gene_count.tsv
  - Auto-generates PF{ID}_summary.txt report
  - Full argparse CLI (positional pfam_id, -i/--input, -o/--output-dir, --delimiter)
  - Semi-automated: auto-resolves paths from species_config if no explicit paths provided

### Fixed
- Column index mapping in parse_domtblout (domtblout v3 22-column format, 0-based)

## [0.1.0] - 2026-05-06

### Added
- Initial repository structure (scripts/, config/, docs/, tests/)
- `scripts/01_download_pfam_hmm.py`: Download Pfam HMM + hmmsearch with full CLI
- `config/species_config.py`: Shared species configuration with env var overrides
- `docs/usage.md`: Detailed usage instructions
- `docs/installation.md`: Installation and dependency guide
- `Makefile`: Build automation targets
- `requirements.txt` / `setup.py`: Python packaging
- `CHANGELOG.md`: Version tracking
