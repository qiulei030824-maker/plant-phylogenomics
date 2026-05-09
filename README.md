# Plant Phylogenomics

A modular pipeline for plant phylogenomics analysis — from Pfam domain identification to phylogenetic tree construction, synteny visualization, and publication-ready figures.

> **Note on repository scope**: This repository currently hosts the shared configuration layer (`config/`), documentation, and packaging infrastructure. The production pipeline scripts are maintained at a separate path (see [Execution Guide](docs/usage.md) for details). A future release will integrate all scripts directly into this repository.

---

## Pipeline Architecture

The full analysis system consists of **13 scripts** organized into two independent pipelines:

### Core Domain Pipeline (4 scripts)

```
01_download_hmm_search_domain.py    Download Pfam HMM → hmmsearch on proteomes
02_parse_domain_table.py             Parse domtblout → gene/species counts
03_extract_domain_sequences.py      Extract domain PEP/CDS per species
04_build_phylogenetic_tree.py       MSA (MAFFT) → trimming → IQ-TREE2 → QC
```

### OrthoFinder Synteny Pipeline (7 Python + 2 R scripts)

```
N1  Species tree + orthogroups (OrthoFinder)
N2  Ortholog tree (MAFFT → trimAl → FastTree → RAxML-NG)
N3  LAST alignment between species
N4  Anchor conversion (→ jcvi SimpleFile)
N5  BED generation (from GFF3)
N6  Macrosynteny chain plot (jcvi karyotype)
N7  Microsynteny plot with gene structure (jcvi synteny)
R1  Classic phylogenetic tree (ggtree)
R2  Tree + gene count heatmap (patchwork)
```

### Dependencies

```
01_download_hmm_search_domain.py  ─┐
02_parse_domain_table.py           ├─ Core pipeline (independent)
03_extract_domain_sequences.py    ─┘
04_build_phylogenetic_tree.py     ─── MSA → tree (final step)

N1 → N2 ──┬── R1, R2 (visualization)
           ├── N3 → N4 ──┬── N6 (macrosynteny)
           │              └── N7 (microsynteny)
           └── N5 ───────→ N6, N7 (BED files)
```

---

## Quick Start

```bash
# Clone the repository
git clone https://github.com/qiulei030824-maker/plant-phylogenomics.git
cd plant-phylogenomics

# Install Python dependencies
pip install -r requirements.txt

# Setup species groups and paths (edit config/species_config.py)
# Then run the core pipeline:

python 01_download_hmm_search_domain.py PF00168 \
    -i /path/to/proteomes \
    -o /path/to/output

python 02_parse_domain_table.py PF00168 \
    -i /path/to/output

python 03_extract_domain_sequences.py PF00168 \
    -i /path/to/output

python 04_build_phylogenetic_tree.py PF00168 \
    -i /path/to/output \
    --profile standard --bootstrap 1000
```

See [docs/usage.md](docs/usage.md) for complete parameter tables and examples.

---

## Repository Structure

```
plant-phylogenomics/
├── config/
│   └── species_config.py      # Species groups, colors, path resolution
├── docs/
│   ├── installation.md        # Dependency installation guide
│   └── usage.md               # Full pipeline usage documentation
├── .gitignore
├── CHANGELOG.md
├── Makefile
├── README.md
├── requirements.txt
└── setup.py
```

---

## Key Features

| Feature | Description |
|---------|-------------|
| **Modular** | Each step is an independent script with well-defined inputs/outputs |
| **CLI-first** | All scripts accept command-line arguments; YAML config also supported |
| **Checkpoint/resume** | Persistent QC reports and run metadata enable safe restarts |
| **Auto-strategy** | Sequence thresholds (<50, 50–200, 200–1000, >1000) auto-select MSA/tree parameters |
| **Profiles** | `fast` / `standard` / `accurate` / `ultra` presets for different scale analyses |
| **Species-aware** | Centralized species groups and colors for consistent visualization |

---

## External Tools Required

| Tool | Version | Used by |
|------|---------|---------|
| HMMER (`hmmsearch`) | ≥3.3 | Core step 1 |
| MAFFT | ≥7.4 | Core step 4, N2 |
| IQ-TREE2 | ≥2.2 | Core step 4 |
| ClipKIT / trimAl | latest | Core step 4, N2 |
| RAxML-NG | ≥1.1 | N2 |
| FastTree | ≥2.1 | N2 |
| OrthoFinder | ≥2.5 | N1 |
| LAST | ≥1027 | N3 |
| jcvi (Python) | ≥1.3 | N4–N7 |

See [docs/installation.md](docs/installation.md) for installation instructions.

---

## License

MIT
