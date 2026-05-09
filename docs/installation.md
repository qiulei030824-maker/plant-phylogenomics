# Installation Guide

## Prerequisites

### System Requirements
- **OS**: Linux / macOS (Windows via WSL2)
- **Python**: 3.8+
- **R**: 4.0+ (for R visualization scripts)
- **Disk**: 10+ GB for genome/proteome data

### External Tools

| Tool | Version tested | Used by |
|------|---------------|---------|
| HMMER (`hmmsearch`) | ≥3.3 | Core step 01 |
| MAFFT | ≥7.4 | Core step 04, N2 |
| IQ-TREE2 | ≥2.2 | Core step 04 |
| ClipKIT | ≥1.3 | Core step 04 |
| trimAl | ≥1.4 | Core step 04, N2 |
| RAxML-NG | ≥1.1 | N2 |
| FastTree | ≥2.1 | N2 |
| OrthoFinder | ≥2.5 | N1 |
| LAST (`lastdb`, `lastal`) | ≥1027 | N3 |
| seqkit | ≥2.0 | Core step 03 |
| jcvi (Python package) | ≥1.3 | N4–N7 |

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

### 3. Install R packages (for R1, R2 visualization)

```r
install.packages(c("ggplot2", "treeio", "ggtree", "ape", "optparse",
                   "jsonlite", "scales", "RColorBrewer", "patchwork", "reshape2"))
# ggtree is from Bioconductor:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree")
```

### 4. Install external bioinformatics tools

All tools can be installed via conda:

```bash
# Core pipeline tools
conda install -c bioconda hmmer mafft iqtree clipkit trimal seqkit

# OrthoFinder pipeline tools
conda install -c bioconda orthofinder raxml-ng fasttree last
```

Or install individually:

```bash
# HMMER (Ubuntu/Debian)
sudo apt install hmmer

# MAFFT
sudo apt install mafft

# Others via conda
conda install -c bioconda iqtree raxml-ng fasttree last
```

## Verification

```bash
# Run help (core pipeline scripts)
python /data5/qiulei/script/0506/01_download_hmm_search_domain.py --help

# Check external tools
which hmmsearch && hmmsearch -h | head -3
which mafft && mafft --help 2>&1 | head -3
which iqtree2 && iqtree2 --help 2>&1 | head -3
```

## Troubleshooting

### "hmmsearch not found"
Set the path explicitly:
```bash
python 01_download_hmm_search_domain.py PF00168 --hmmsearch /custom/path/hmmsearch
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
python 01_download_hmm_search_domain.py PF00168 -o /tmp/my_analysis
```
