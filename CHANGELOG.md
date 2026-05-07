# Changelog

All notable changes to this project will be documented in this file.

## [1.0.0] - 2026-05-07

### Added
- `scripts/04_build_domain_tree.py`: **Industrial-grade rewrite** (~1229 lines, up from ~400)
  - **Profile system**: 4 presets (`fast`/`standard`/`accurate`/`ultra`) + `auto` (AutoStrategy)
    - fast: fftns2 + trimal-gappyout + IQ-TREE LG + 100 BS
    - standard: localpair + MFP + 1000 UFBoot (best for most use cases)
    - accurate: linsi + MFP+MERGE + 2000 UFBoot
    - ultra: einsi + MFP+MERGE + 5000 UFBoot
    - AutoStrategy: auto-decides based on n_seqs (<50→linsi, 50-200→einsi, 200-1000→fftns2, >1000→fftns2+LG)
  - **6 sequence selection strategies**: `all`, `longest`, `canonical`, `domain_best` (HMMER domtblout parsing), `longest_isoform`, `representative` (cd-hit/mmseqs clustering)
  - **MAFFT alignment**: auto/linsi/einsi/fftns2 modes with auto fallback logic
  - **Trimming engines**: clipkit (kpi-smart-gap/smart-gap/strict/adaptive) + trimal (gappyout/automated1/strictplus)
  - **IQ-TREE2**: model selection MFP/MFP+MERGE/LG, conditional -b (standard) vs -B -bnni (UFBoot2) based on sequence count
  - **Checkpoint/resume system**: 5 checkpoints (01_input→02_filter→03_align→04_trim→05_tree), auto-skip completed steps
  - **QC report**: alignment statistics (sequences, sites, gaps, identity, occupancy)
  - **Metadata**: `run_metadata.json` with full provenance tracking
  - **Auto-generated methods paragraph**: ready for publication methods section
  - **Tree visualization**: ete3 → toytree → matplotlib fallback, PDF/SVG/PNG output
  - **Taxonomy filtering**: --clade / --family / --genus / --species (48 species, 6 clades, 8 groups)
  - **Large family optimization**: auto cd-hit when >2000 sequences
  - **YAML config support**: `--config config.yaml` for reproducible runs
  - **Full logging**: simultaneous file + console output
  - **20+ CLI parameters** with argparse, sensible defaults
- `config/species_config.py`: Shared configuration with TAXONOMY_HIERARCHY (6 clades: Rosids, Asterids, Monocots, Basal_Angiosperms, Bryophytes, Lycophytes), GROUP_COLORS, 48 species in 8 groups, new helper functions (get_genus, get_family, get_group_color, parse_newick_tips, resolve_pfam_paths)

### Changed
- `scripts/04_build_domain_tree.py`: Complete rewrite from basic script to industrial-grade engine (v0.4 → v1.0.0)
- `config/species_config.py`: Expanded from basic path config to full taxonomy/group management

### Fixed
- argparse: Strategy default value not propagating correctly
- `__main__` guard: Missing `if __name__ == "__main__": main()` protection
- clipkit CLI: Incompatible command-line parameter format
- IQ-TREE2: Bootstrap parameter logic for large sequence sets

### Tested
- C2 domain (all species, ~1200 sequences): profile=accurate, strategy=domain_best — **PASS**
- PF00168 (Rosids clade, 19 species): profile=standard, strategy=longest — **PASS**
- Resume/checkpoint mechanism — **PASS**

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
