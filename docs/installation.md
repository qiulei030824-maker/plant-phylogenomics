# Installation Guide

## Prerequisites

### System Requirements
- **OS**: Linux / macOS (Windows via WSL2)
- **Python**: 3.8+
- **R**: 4.0+ (for visualization scripts)
- **Disk**: 10+ GB for genome/proteome data

### External Tools

| Tool | Version tested | Purpose |
|------|---------------|---------|
| HMMER (hmmsearch) | ≥3.3 | Profile HMM search |
| MAFFT | ≥7.4 | Multiple sequence alignment |
| RAxML-NG / FastTree | latest | Maximum likelihood phylogenetic trees |
| trimAl | ≥1.4 | Alignment trimming |

## Installation Steps

### 1. Clone the repository

```bash
git clone https://github.com/qiulei030824-maker/plant-phylogenomics.git
cd plant-phylogenomics
```

### 2. Install Python dependencies

```bash
pip install -r requirements.txt
```

Or using setup.py:

```bash
pip install -e .
```

### 3. Install R packages (for visualization)

```r
install.packages(c("ggplot2", "treeio", "ggtree", "ape"))
# ggtree is from Bioconductor:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree")
```

### 4. Install external bioinformatics tools

#### HMMER
```bash
# conda
conda install -c bioconda hmmer

# apt (Ubuntu/Debian)
sudo apt install hmmer

# source
wget http://eddylab.org/software/hmmer/hmmer.tar.gz
tar -xzf hmmer.tar.gz && cd hmmer-*
./configure && make && sudo make install
```

#### MAFFT
```bash
conda install -c bioconda mafft
# or
sudo apt install mafft
```

#### RAxML-NG
```bash
conda install -c bioconda raxml-ng
```

#### FastTree
```bash
conda install -c bioconda fasttree
```

## Verification

Verify installation:

```bash
# Check Python imports
python -c "from scripts._version import __version__; print('OK')"

# Run help
python scripts/01_download_pfam_hmm.py --help

# Check external tools
which hmmsearch && hmmsearch -h | head -3
```

## Troubleshooting

### "hmmsearch not found"
Set the path explicitly:
```bash
python scripts/01_download_pfam_hmm.py PF00168 --hmmsearch /custom/path/hmmsearch
```

### "No module named 'config'"
Run from the repository root:
```bash
cd /path/to/plant-phylogenomics
python -m scripts.01_download_pfam_hmm PF00168
```

### Permission denied for output directory
Use a custom output dir:
```bash
python scripts/01_download_pfam_hmm.py PF00168 -o /tmp/my_analysis
```
