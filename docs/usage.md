# Usage Guide

## Installation

See [installation.md](installation.md) for full dependency setup.

## Quick Start

### 1. Download Pfam HMM and run hmmsearch

```bash
# Basic usage (uses default paths from config)
python scripts/01_download_pfam_hmm.py PF00168

# Custom input/output directories
python scripts/01_download_pfam_hmm.py PF00168 \
    -i /path/to/genome/pep/files \
    -o /path/to/output/root

# Specify hmmsearch binary path
python scripts/01_download_pfam_hmm.py PF00168 \
    --hmmsearch /usr/local/bin/hmmsearch

# Download HMM only (skip hmmsearch)
python scripts/01_download_pfam_hmm.py PF00168 --skip-hmmsearch

# Force re-download HMM
python scripts/01_download_pfam_hmm.py PF00168 --force-redownload

# Custom PEP file suffix
python scripts/01_download_pfam_hmm.py PF00168 --pep-suffix ".proteome.fa"

# Skip additional subdirectories
python scripts/01_download_pfam_hmm.py PF00168 --skip-dirs "extra_data" "backup"
```

### 2. Extract domain sequences

```bash
# Coming soon
# python scripts/02_extract_domain_seqs.py PF00168
```

### 3. Build phylogenetic tree

```bash
# Coming soon
# python scripts/03_build_domain_tree.py PF00168
```

## Configuration

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `PFAM_DATA_DIR` | `/data5/qiulei/PPI/data` | Directory with species PEP subdirectories |
| `PFAM_BASE_ROOT` | `/data5/qiulei/PPI/03.OrthoFinder_MCScanX/0430` | Output root directory |
| `PFAM_SKIP_DIRS` | `logs,pfamdb,__pycache__,...` | Comma-separated dirs to skip |

Example:
```bash
export PFAM_DATA_DIR=/custom/genome/path
python scripts/01_download_pfam_hmm.py PF00168
```

### Species Configuration

Edit `config/species_config.py` to:
- Add/remove species in species groups
- Change group colors for tree visualization
- Modify species lists for different plant clades

## Pipeline Steps

```
Step 1: Pfam HMM Download + hmmsearch
    ↓
Step 2: Domain Sequence Extraction
    ↓
Step 3: MSA + Phylogenetic Tree
    ↓
Step 4: Downstream Analysis (dN/dS, synteny, etc.)
    ↓
Step R1: Tree Visualization (R/ggtree)
```

## Output Structure

```
<output_root>/
└── <PFAM_ID>/
    └── hmmer/
        ├── <PFAM_ID>.hmm              # Pfam HMM file
        ├── all_<PFAM_ID>.domtblout     # Combined hmmsearch results
        ├── all_proteomes_<PFAM_ID>.fa  # Concatenated proteomes
        └── hmmsearch.log               # hmmsearch log
    ├── pep/                            # Extracted domain sequences
    ├── algin/                          # Multiple sequence alignments
    └── tree/                           # Phylogenetic trees
```

## Adding New Scripts

1. Place the script in `scripts/` with a step-numbered prefix
2. If it needs shared config, import from `config.species_config`
3. Use argparse for CLI parameters (see `01_download_pfam_hmm.py` as template)
4. Update the pipeline table in `README.md`
5. Add entry to `CHANGELOG.md`
