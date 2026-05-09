# Usage Guide

## Pipeline Locations

Production scripts are maintained at two locations:

| Pipeline | Path | Scripts |
|----------|------|---------|
| Core domain | `/data5/qiulei/script/0506/` | 01–04 |
| OrthoFinder synteny | `/data5/qiulei/script/0506/OrthoFinder`$'\n'`/` | N1–N7, R1–R2 |

> The OrthoFinder directory name contains a trailing newline character. Use a glob to access: `"OrthoFinder"*`

---

## Core Domain Pipeline (4 scripts)

### 01 — Pfam HMM Download & hmmsearch

**File:** `01_download_hmm_search_domain.py`

Downloads a Pfam HMM profile, builds a combined proteome database, and runs `hmmsearch`.

**Inputs:** Species proteome PEP files in a directory (one `.pep.fa` per species, or custom suffix via `--pep-suffix`)

**Outputs:**
```
<output_root>/hmmer/
├── <PFAM_ID>.hmm                    # Pfam HMM
├── <PFAM_ID>.hmm.h3[fimop]          # HMMER binary indices
├── all_<PFAM_ID>.domtblout          # Concatenated hmmsearch results (domtblout)
├── all_proteomes_<PFAM_ID>.fa       # Concatenated proteome sequences
└── hmmsearch.log
```

**Parameters**

| Arg | Required | Default | Description |
|-----|----------|---------|-------------|
| `pfam_id` | ✓ | — | Pfam accession (e.g., PF00168) |
| `-i` | ✗ | `config.species_config.DATA_DIR` | Input proteomes directory |
| `-o` | ✗ | `config.species_config.BASE_ROOT` | Output root directory |
| `--hmmsearch` | ✗ | `hmmsearch` | hmmsearch binary path |
| `--skip-dirs` | ✗ | from config | Additional directories to skip |
| `--skip-hmmsearch` | ✗ | `False` | Only download HMM, skip search |
| `--force-redownload` | ✗ | `False` | Re-download HMM even if exists |
| `--pep-suffix` | ✗ | `.pep.fa` | Proteome file suffix |

**Examples**

```bash
# Minimal
python 01_download_hmm_search_domain.py PF00168

# Custom paths
python 01_download_hmm_search_domain.py PF00168 \
    -i /path/to/proteomes \
    -o /path/to/output

# Specify hmmsearch binary
python 01_download_hmm_search_domain.py PF00168 \
    --hmmsearch /usr/local/bin/hmmsearch

# Download only (skip hmmsearch)
python 01_download_hmm_search_domain.py PF00168 --skip-hmmsearch
```

---

### 02 — Parse Domain Table

**File:** `02_parse_domain_table.py`

Parses the concatenated hmmsearch domtblout into structured TSV reports.

**Inputs:** `all_<PFAM_ID>.domtblout` (from step 01)

**Outputs:**
```
<output_root>/
├── gene_domain_count.tsv         # Per-gene domain count
├── species_gene_count.tsv        # Per-species gene/repeat summary
└── summary.txt                   # Text summary report
```

**Parameters**

| Arg | Required | Default | Description |
|-----|----------|---------|-------------|
| `pfam_id` | ✓ | — | Pfam accession |
| `-i` | ✗ | `config.species_config.DATA_DIR` | Input proteomes directory |
| `-o` | ✗ | `config.species_config.BASE_ROOT` | Output root directory |
| `--delimiter` | ✗ | `|` | Protein ID delimiter |

**Examples**

```bash
python 02_parse_domain_table.py PF00168 \
    -i /path/to/output \
    -o /path/to/output
```

**Common pitfalls**
- The output root should match step 01's `-o` path so the script can find `hmmer/all_<PFAM_ID>.domtblout`

---

### 03 — Extract Domain Sequences

**File:** `03_extract_domain_sequences.py`

Extracts PEP and CDS sequences for each species using `seqkit grep`, based on gene IDs identified by step 02. Each species gets its own output file.

**Inputs:**
- `gene_domain_count.tsv` (from step 02)
- Species PEP files in `<data_dir>/<species>/`
- Species CDS files in `<data_dir>/<species>/`

**Outputs:**
```
<output_root>/pep/
├── <species_1>_<PFAM_ID>.pep.fa
├── <species_2>_<PFAM_ID>.pep.fa
└── ...

<output_root>/cds/
├── <species_1>_<PFAM_ID>.cds.fa
├── <species_2>_<PFAM_ID>.cds.fa
└── ...
```

**Parameters**

| Arg | Required | Default | Description |
|-----|----------|---------|-------------|
| `pfam_id` | ✓ | — | Pfam accession |
| `-i` | ✓ | — | Output root from step 01/02 |
| `-o` | ✓ | — | Output directory |
| `--data-dir` | ✗ | `config.species_config.DATA_DIR` | Species genome data directory |
| `--skip-dirs` | ✗ | from config | Directories to skip |
| `--delimiter` | ✗ | `|` | Protein ID delimiter |
| `--skip-cds` | ✗ | `False` | Skip CDS extraction |
| `--force` | ✗ | `False` | Overwrite existing output |
| `--pep-suffix` | ✗ | `.pep.fa` | PEP file suffix |
| `--cds-suffix` | ✗ | `.cds.fa` | CDS file suffix |

**Examples**

```bash
python 03_extract_domain_sequences.py PF00168 \
    -i /path/to/output \
    -o /path/to/output \
    --data-dir /path/to/genome_data

# Skip CDS extraction
python 03_extract_domain_sequences.py PF00168 \
    -i /path/to/output \
    -o /path/to/output --skip-cds
```

---

### 04 — Build Phylogenetic Tree

**File:** `04_build_phylogenetic_tree.py`

The main phylogenetic pipeline: reads domain PEP sequences → MAFFT alignment → ClipKIT/trimAl trimming → IQ-TREE2 tree inference → QC report + methods text generation.

**Auto-strategy engine** selects parameters based on sequence count:

| Sequences | MSA | Trimming | Tree |
|-----------|-----|----------|------|
| <50 | `--localpair --maxiterate 1000` | ClipKIT `-m gappy` | IQ-TREE2 `-m MFP -B 1000` |
| 50–200 | `--genafpair --maxiterate 1000` | ClipKIT `-m gappy` | IQ-TREE2 `-m MFP -B 1000` |
| 200–1000 | `--auto` | trimAl `-gt 0.8` | IQ-TREE2 `-m MFP -B 1000` |
| >1000 | `--auto` | trimAl `-gt 0.8` | IQ-TREE2 `-m MFP` (no bootstrap) |

**Preset profiles** override auto-strategy:

| Profile | MSA | Trimming | Bootstrap | Notes |
|---------|-----|----------|-----------|-------|
| `fast` | `--auto` | trimAl `-gt 0.8` | 100 | Quick exploration |
| `standard` | `--localpair --maxiterate 1000` | ClipKIT `-m gappy` | 1000 | Default |
| `accurate` | `--localpair --maxiterate 1000` | ClipKIT `-m gappy-kpi` | 1000 | Higher precision |
| `ultra` | `--localpair --maxiterate 1000` | ClipKIT `-m kpi` | 5000 | Publication-grade |

**Inputs:**
- `<output_root>/pep/<species>_<PFAM_ID>.pep.fa` (from step 03)
- Optionally: `<output_root>/cds/` for codon alignment

**Outputs:**
```
<output_root>/algin/
├── combined_unaligned.fa          # Concatenated unaligned sequences
├── combined.aln.fa                # MAFFT alignment
├── combined.trim.fa               # Trimmed alignment
├── combined.aln.html              # Alignment visualization (PDF/SVG/PNG)

<output_root>/tree/
├── <PFAM_ID>_tree.contree         # IQ-TREE2 consensus tree
├── <PFAM_ID>_tree.treefile        # IQ-TREE2 best tree
├── <PFAM_ID>_tree.iqtree          # Full IQ-TREE2 report
├── <PFAM_ID>_tree.log             # IQ-TREE2 log
├── <PFAM_ID>_tree.pdf             # Tree visualization (PDF)
├── <PFAM_ID>_tree.svg             # (SVG)
├── <PFAM_ID>_tree.png             # (PNG)
├── qc_report.txt                  # Quality control metrics
├── run_metadata.json              # Pipeline provenance
└── methods.txt                    # Auto-generated methods text
```

**Parameters**

| Arg | Required | Default | Description |
|-----|----------|---------|-------------|
| `pfam_id` | ✓ | — | Pfam accession |
| `-i` | ✗ | `config.species_config.BASE_ROOT` | Output root |
| `-o` | ✗ | `config.species_config.BASE_ROOT` | Output root |
| `--profile` | ✗ | `standard` | Preset: fast/standard/accurate/ultra |
| `--strategy` | ✗ | `auto` | Strategy: auto/fast/standard/accurate/ultra |
| `--clade` / `--family` / `--genus` / `--species` | ✗ | — | Taxonomic filtering |
| `--align-mode` | ✗ | `pep` | pep/cds/codon |
| `--trim` | ✗ | `auto` | Trimming tool: auto/clipkit/trimal/none |
| `--bootstrap` | ✗ | `auto` | Bootstrap replicates (auto/0–10000) |
| `--force` | ✗ | `False` | Overwrite existing |
| `--resume` | ✗ | `False` | Resume from checkpoint |
| `--config` | ✗ | None | YAML config file |

**Examples**

```bash
# Default (standard profile)
python 04_build_phylogenetic_tree.py PF00168

# Fast exploration
python 04_build_phylogenetic_tree.py PF00168 --profile fast

# Publication-grade
python 04_build_phylogenetic_tree.py PF00168 \
    --profile accurate --bootstrap 5000

# Filter to specific clade
python 04_build_phylogenetic_tree.py PF00168 \
    --clade Cucurbits Brassicaceae

# Resume interrupted run
python 04_build_phylogenetic_tree.py PF00168 --resume
```

**Common pitfalls**
- Large alignments (>1000 sequences) skip bootstrap by default for performance
- The `--resume` flag checks for existing checkpoint files; re-runs only incomplete steps

---

## OrthoFinder Synteny Pipeline (N1–N7 + R1–R2)

### N1 — Species Tree & Orthogroups

**File:** `stepN1_build_species_tree.py`

Runs OrthoFinder to build a species tree and identify orthogroups from proteome FASTA files.

**Inputs:** Species proteome FASTA files in a directory (`.fa`, `.faa`, or `.fasta`)

**Outputs:**
```
<output_dir>/proteomes_<output_name>/         # Symlinked inputs
<output_dir>/Results_<output_name>/
├── Species_Tree/SpeciesTree_rooted.txt
├── Single_Copy_Orthologue_Sequences/
├── Orthogroups/Orthogroups.tsv
└── Gene_Trees/
```

**Parameters**

| Arg | Required | Default | Description |
|-----|----------|---------|-------------|
| `--proteome-dir` | ✓ | — | Proteome FASTA directory |
| `--output-name` | ✗ | `my_species` | Output name prefix |
| `--threads` | ✗ | 32 | OrthoFinder search threads |
| `--ortho-threads` | ✗ | 8 | OrthoFinder analysis threads |
| `--orthofinder` | ✗ | `orthofinder` | OrthoFinder binary |
| `--output-dir` | ✗ | `./data/orthofinder` | Output root |
| `--config` | ✗ | None | YAML config |

**Idempotency:** Skips if `SpeciesTree_rooted.txt` exists.

**Example**

```bash
python stepN1_build_species_tree.py \
    --proteome-dir /path/to/proteomes \
    --output-name my_species \
    --output-dir /path/to/orthofinder \
    --threads 16
```

---

### N2 — Ortholog Tree

**File:** `stepN2_build_ortholog_tree.py`

Finds the orthogroup containing a target gene → MAFFT alignment → trimAl → FastTree → RAxML-NG.

**Inputs:**
- Orthogroups.tsv (from N1)
- Single copy orthologue sequences (from N1)

**Outputs:**
```
<output_root>/<pfam_id>/tree/ortholog_tree_<target_gene>/
├── <target>.aln.fa
├── <target>.trim.fa
├── <target>.fasttree.nwk
├── <target>.raxml.bestTree
└── <target>.raxml.support
```

**Parameters**

| Arg | Required | Default | Description |
|-----|----------|---------|-------------|
| `--target-gene` | ✓ | — | Target gene ID (e.g., AT1G09070) |
| `--orthofinder-dir` | ✗ | (hardcoded) | OrthoFinder results directory |
| `--pfam-id` | ✗ | PF00168 | Pfam ID |
| `--pipeline-root` | ✗ | (hardcoded) | Pipeline data root |
| `--threads` | ✗ | 8 | CPU threads |
| `--bootstrap` | ✗ | 100 | RAxML-NG bootstrap replicates |

**Example**

```bash
python stepN2_build_ortholog_tree.py \
    --target-gene AT1G09070 \
    --orthofinder-dir /path/to/orthofinder/Results_my_species \
    --pfam-id PF00168 \
    --pipeline-root /path/to/data
```

---

### N3 — LAST Alignment

**File:** `stepN3_run_last_alignment.py`

Builds LAST databases and runs pairwise alignment between adjacent species in a chain.

**Inputs:** `<pipeline_root>/<pfam_id>/cds/<species>.cds.fa`

**Outputs:**
```
<pipeline_root>/<pfam_id>/jcvi_chain/
├── <sp_a>.<sp_b>.last
└── <sp_a>.<sp_b>.last.filtered
```

**Parameters**

| Arg | Required | Default | Description |
|-----|----------|---------|-------------|
| `--pfam-id` | ✓ | — | Pfam ID |
| `--species` | ✓ | — | Species list (ordered, 2+) |
| `--pipeline-root` | ✗ | (hardcoded) | Pipeline data root |
| `--threads` | ✗ | 4 | LAST alignment threads |

**Example**

```bash
python stepN3_run_last_alignment.py \
    --pfam-id PF00168 \
    --species ChineseLong arabidopsis_thaliana solanum_lycopersicum \
    --pipeline-root /path/to/data
```

---

### N4 — Anchor Conversion

**File:** `stepN4_convert_anchors.py`

Converts LAST alignment data to jcvi SimpleFile (6-column tab-separated) format.

**Inputs:** `.last.filtered` or `.last` (from N3)

**Outputs:** `<pipeline_root>/<pfam_id>/jcvi_chain/<pair>.simple`

**Parameters**

| Arg | Required | Default | Description |
|-----|----------|---------|-------------|
| `--pfam-id` | ✓ | — | Pfam ID |
| `--species` | ✓ | — | Species list (matches N3 order) |
| `--pipeline-root` | ✗ | (hardcoded) | Pipeline data root |

**Example**

```bash
python stepN4_convert_anchors.py \
    --pfam-id PF00168 \
    --species ChineseLong arabidopsis_thaliana solanum_lycopersicum
```

---

### N5 — BED Generation

**File:** `stepN5_generate_bed.py`

Generates BED files for each species (from existing `.gene.bed` or by converting GFF3 via jcvi).

**Inputs:** GFF3 or existing `.gene.bed` in `<data_dir>/<species>/`

**Outputs:**
```
<pipeline_root>/<pfam_id>/jcvi_chain/
├── <species>.bed
└── seqids.txt
```

**Parameters**

| Arg | Required | Default | Description |
|-----|----------|---------|-------------|
| `--pfam-id` | ✓ | — | Pfam ID |
| `--species` | ✓ | — | Species list |
| `--pipeline-root` | ✗ | (hardcoded) | Pipeline data root |
| `--data-dir` | ✗ | None | Species data directory |

**Example**

```bash
python stepN5_generate_bed.py \
    --pfam-id PF00168 \
    --species ChineseLong arabidopsis_thaliana solanum_lycopersicum \
    --data-dir /path/to/genome_data
```

---

### N6 — Macrosynteny Chain Plot

**File:** `stepN6_highlight_chain.py`

Generates a chromosome-level synteny plot with highlighted gene pairs using jcvi karyotype.

**Inputs:** `.simple` (N4), `.bed` (N5), `seqids.txt` (N5)

**Outputs:** `<pipeline_root>/<pfam_id>/jcvi_chain/chain_synteny_fancy.pdf`

**Parameters**

| Arg | Required | Default | Description |
|-----|----------|---------|-------------|
| `--pfam-id` | ✓ | — | Pfam ID |
| `--species` | ✓ | — | Species list (chain order) |
| `--pipeline-root` | ✗ | (hardcoded) | Pipeline data root |
| `--highlight-pairs-json` | ✗ | None | JSON with gene pairs to highlight |

**Highlight JSON format**
```json
{
    "ChineseLong.arabidopsis_thaliana": [
        ["CsaV3_2G006560", "AT1G09070"]
    ]
}
```

**Example**

```bash
python stepN6_highlight_chain.py \
    --pfam-id PF00168 \
    --species ChineseLong arabidopsis_thaliana solanum_lycopersicum \
    --highlight-pairs-json highlight.json
```

---

### N7 — Microsynteny Plot

**File:** `stepN7_microsynteny.py`

Creates a gene-level microsynteny plot showing gene neighborhood around target genes, with exon structure.

**Inputs:** `.stripped.simple` (N6), `.bed` (N5), GFF3, config JSON

**Outputs:** `<pipeline_root>/<pfam_id>/jcvi_chain/microsynteny.pdf`

**Parameters**

| Arg | Required | Default | Description |
|-----|----------|---------|-------------|
| `--pfam-id` | ✓ | — | Pfam ID |
| `--config-json` | ✓ | — | JSON config with targets/tracks/edges |
| `--pipeline-root` | ✗ | (hardcoded) | Pipeline data root |
| `--neighborhood` | ✗ | 10 | Genes flanking each target |
| `--target-labels` | ✗ | "" | Comma-separated genes to label |

**Config JSON structure**
```json
{
    "targets": {
        "ChineseLong": ["CsaV3_1G038880"],
        "arabidopsis_thaliana": ["AT1G09070"]
    },
    "tracks_def": [
        ["ChineseLong", "CsaV3_1G038880"],
        ["arabidopsis_thaliana", "AT1G09070"]
    ],
    "track_edges": [[0, 1]]
}
```

**Example**

```bash
python stepN7_microsynteny.py \
    --pfam-id PF00168 \
    --config-json microsynteny_config.json \
    --pipeline-root /path/to/data \
    --neighborhood 10
```

---

### R1 — Classic Tree Visualization

**File:** `stepR1_orthofinder_tree.R`

Reads a RAxML-NG support tree and produces a publication-grade ggtree figure with group colors, bootstrap labels, and species labels.

**Inputs:** `.raxml.support` (from N2)

**Outputs:** `<output_dir>/<og_name>_tree_classic.pdf`

**Parameters**

| Arg | Required | Default | Description |
|-----|----------|---------|-------------|
| `--og` | ✓ | — | Orthogroup name (e.g., AT1G09070_OG0000001) |
| `--tree` | ✓ | — | Newick tree file path |
| `--output-dir` | ✗ | (hardcoded) | Output directory |
| `--width` | ✗ | 16 | PDF width (inches) |
| `--height` | ✗ | 12 | PDF height |
| `--tip-size` | ✗ | 2.5 | Tip label font size |

**Example**

```bash
Rscript stepR1_orthofinder_tree.R \
    --og AT1G09070_OG0000001 \
    --tree /path/to/tree.raxml.support \
    --output-dir /path/to/output \
    --width 16 --height 12
```

---

### R2 — Tree + Heatmap

**File:** `stepR2_orthofinder_heatmap.R`

Combines a ggtree phylogeny with a gene-count heatmap using patchwork (left panel: tree, right panel: YlOrRd heatmap).

**Outputs:** `<output_dir>/<og_name>_tree_heatmap.pdf`

**Parameters:** Same as R1.

**Example**

```bash
Rscript stepR2_orthofinder_heatmap.R \
    --og AT1G09070_OG0000001 \
    --tree /path/to/tree.raxml.support
```

---

## Output Structure Summary

```
<project_root>/
├── hmmer/                           # Core step 01: HMM + hmmsearch
│   ├── <PFAM_ID>.hmm
│   └── all_<PFAM_ID>.domtblout
├── pep/                             # Core step 03: domain PEP sequences
│   └── <species>_<PFAM_ID>.pep.fa
├── cds/                             # Core step 03: domain CDS sequences
│   └── <species>_<PFAM_ID>.cds.fa
├── algin/                           # Core step 04: alignments
│   ├── combined.aln.fa
│   └── combined.trim.fa
├── tree/                            # Core step 04 + N2: phylogenetic trees
│   ├── <PFAM_ID>_tree.treefile
│   ├── <PFAM_ID>_tree.pdf
│   └── ortholog_tree_<gene>/
├── jcvi_chain/                      # N3–N7: synteny analysis
│   ├── *.last / *.last.filtered
│   ├── *.simple
│   ├── *.bed
│   ├── chain_synteny_fancy.pdf
│   └── microsynteny.pdf
├── orthofinder/                     # N1: OrthoFinder results
├── gene_domain_count.tsv            # Core step 02
└── species_gene_count.tsv           # Core step 02
```

## Configuration

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `PFAM_DATA_DIR` | `/data5/qiulei/PPI/data` | Species PEP subdirectories |
| `PFAM_BASE_ROOT` | `/data5/qiulei/PPI/03.OrthoFinder_MCScanX/0430` | Output root |
| `PFAM_SKIP_DIRS` | `logs,pfamdb,...` | Comma-separated dirs to skip |

### Species Config

Edit `config/species_config.py` to customize species groups, colors, or path defaults.

---

## Reproducibility

- Step 04 generates `run_metadata.json` with all parameters and versions
- Step 04 generates `methods.txt` for direct use in publications
- YAML config files can capture full pipeline runs (see `configs/` in the OrthoFinder pipeline directory)

## Naming Conventions

- Script prefix: `N` for core pipeline, `R` for R visualization
- Step numbers: 01–04 for core pipeline; N1–N7 for synteny pipeline
- Output tree files: `<pfam_id>_tree.treefile` (core) or `<target_gene>_raxml.support` (OrthoFinder)
