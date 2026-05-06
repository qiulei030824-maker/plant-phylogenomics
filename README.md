# Plant Phylogenomics

A modular, production-grade pipeline for plant phylogenomics analysis — from Pfam domain identification to phylogenetic tree construction and visualization.

## Repository Structure

```
plant-phylogenomics/
├── scripts/          # Core analysis scripts (step-numbered)
│   ├── 01_download_pfam_hmm.py
│   ├── 02_extract_domain_seqs.py
│   ├── 03_build_domain_tree.py
│   ├── 04_dnds_analysis.py
│   └── 04_run_last_alignment.py
├── config/           # Shared configuration files
│   ├── species_config.py       # Species groups, colors, path resolution
│   └── default.yaml            # YAML-based configuration template
├── docs/             # Documentation
│   ├── usage.md                # Detailed usage instructions
│   └── installation.md         # Installation & dependency guide
├── tests/            # Unit tests
├── workflows/        # Pipeline workflow definitions
│   └── steps/                  # Individual workflow steps
├── data/             # Data directory (gitignored)
├── Makefile          # Build automation
├── requirements.txt  # Python dependencies
├── setup.py          # Python package setup
└── CHANGELOG.md      # Version history
```

## Quick Start

```bash
# Clone the repository
git clone https://github.com/qiulei030824-maker/plant-phylogenomics.git
cd plant-phylogenomics

# Install dependencies
pip install -r requirements.txt

# Run a Pfam domain analysis
python scripts/01_download_pfam_hmm.py PF00168
```

## Key Features

- **Modular step-by-step design**: Each step is an independent script with a well-defined interface
- **CLI-first**: All scripts accept command-line arguments for flexible usage
- **Configurable**: Species groups, colors, and paths are centralized in `config/species_config.py`
- **Reproducible**: Complete provenance from HMM download to phylogenetic tree
- **Extensible**: Add new analysis scripts following the existing conventions

## Pipeline Overview

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_download_pfam_hmm.py` | Download Pfam HMM → run hmmsearch on proteomes |
| 2 | `02_extract_domain_seqs.py` | Extract domain sequences from hmmsearch results |
| 3 | `03_build_domain_tree.py` | Multiple sequence alignment → phylogenetic tree |
| 4 | `04_dnds_analysis.py` | dN/dS ratio analysis |
| 4 | `04_run_last_alignment.py` | LAST alignment for synteny analysis |
| R1 | `R/01_classic_tree.R` | Tree visualization with ggtree |

## License

MIT
