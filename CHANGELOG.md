# Changelog

All notable changes to this project will be documented in this file.

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
