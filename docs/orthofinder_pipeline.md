# OrthoFinder Pipeline Documentation

## Overview

This pipeline integrates OrthoFinder-based orthology inference, LAST alignment for macrosynteny analysis, and gene structure visualization. It consists of three sub-pipelines:

1. **N-Series (OrthoFinder Analysis)**: Species tree, ortholog trees, synteny analysis
2. **R-Series (Visualization)**: Tree visualization and heatmap generation
3. **Domain Pipeline**: HMM-based domain analysis (see scripts/domain_pipeline/)

## Pipeline Steps

### N1: Build Species Tree (OrthoFinder)
- Symlinks proteome files and runs OrthoFinder
- Uses MAFFT + FastTree for MSA mode species tree
- Output: Rooted species tree + Orthogroups.tsv

### N2: Build Ortholog Tree
- Finds orthogroup containing target gene from Orthogroups.tsv
- Extracts sequences, MAFFT alignment, trimAl trimming
- FastTree preview + RAxML-NG final tree with bootstrapping

### N3: LAST Alignment
- Builds LAST database for each species' CDS
- Runs pairwise lastal between consecutive species
- Filters best hit per query

### N4: Convert Anchors
- Reads jcvi anchors + LAST alignment scores
- Merges into jcvi SimpleFile (.simple) format

### N5: Generate BED files
- Converts GFF3 -> BED via jcvi, or finds existing BEDs
- Creates seqids.txt for chromosome ordering

### N6: Macro-synteny Highlight Chain
- Runs jcvi.graphics.karyotype for chain layout
- Highlights specific gene pairs

### N7: Microsynteny Plot
- Extracts neighborhoods around target genes
- Builds blocks, merged BED, extra BED (exons)
- Creates layout file + runs jcvi.graphics.synteny

### R1: OrthoFinder Tree Visualization
- R script using ggtree for species tree visualization
- Group coloring with taxonomic groups

### R2: Tree + Gene Count Heatmap
- R script combining ggtree tree + heatmap of gene counts
- Uses patchwork for combined layout

## Quick Start

```bash
# 1. Build species tree
python scripts/orthofinder_pipeline/stepN1_build_species_tree.py \
    --config configs/orthofinder_pipeline/stepN1_species_tree.yaml

# 2. Build ortholog tree for a target gene
python scripts/orthofinder_pipeline/stepN2_build_ortholog_tree.py \
    --target-gene AT1G09070 --threads 8

# 3-7. Synteny analysis
python scripts/orthofinder_pipeline/stepN3_run_last_alignment.py --config-json configs/orthofinder_pipeline/stepN3_last_alignment.yaml
python scripts/orthofinder_pipeline/stepN4_convert_anchors.py --config-json configs/orthofinder_pipeline/stepN4_convert_anchors.yaml
python scripts/orthofinder_pipeline/stepN5_generate_bed.py --config-json configs/orthofinder_pipeline/stepN5_generate_bed.yaml
python scripts/orthofinder_pipeline/stepN6_highlight_chain.py --pfam-id PF00168 --config-json configs/orthofinder_pipeline/stepN6_highlight_chain.yaml
python scripts/orthofinder_pipeline/stepN7_microsynteny.py --pfam-id PF00168 --config-json configs/orthofinder_pipeline/stepN7_microsynteny.yaml

# R1-R2: Visualize results
Rscript scripts/visualization/stepR1_orthofinder_tree.R --tree results/species_tree.nwk --output-dir results/
Rscript scripts/visualization/stepR2_orthofinder_heatmap.R --tree results/species_tree.nwk --output-dir results/
```

## Output Structure

```
<pfam_id>/jcvi_chain/
  ├── *.bed                      # Species BED files
  ├── *.simple                   # Pairwise .simple files
  ├── seqids.txt                 # Chromosome ordering
  ├── layout.csv                 # Karyotype layout
  ├── microsynteny.blocks        # Synteny blocks
  ├── microsynteny.bed           # Merged BED
  ├── microsynteny.extra.bed     # Exon features
  ├── microsynteny.layout        # Synteny layout
  └── microsynteny.pdf           # Final microsynteny plot

<pfam_id>/ortholog_tree/
  ├── OG*_raw.fa                 # Extracted sequences
  ├── OG*_aligned.fa             # MAFFT alignment
  ├── OG*_trimmed.fa             # trimAl-trimmed
  ├── OG*_quick.tre              # FastTree preview
  └── raxml/OG*.*               # RAxML-NG results
```

## Dependencies

- Python 3.8+: BioPython, PyYAML, jcvi
- External tools: OrthoFinder, MAFFT, trimAl, FastTree, RAxML-NG, LAST
- R 4.0+: ape, ggtree, ggplot2, RColorBrewer, patchwork, reshape2

## Configuration

Each step has a corresponding YAML config in `configs/orthofinder_pipeline/`. Edit paths in these files to match your data locations before running.

## Notes

- The N1 pipeline (OrthoFinder) requires significant memory for large proteome sets
- N3-N7 work with a specific Pfam ID and require prior domain analysis
- R scripts expect a species prefix naming convention for tip labels
- All paths default to `/data5/qiulei/pfam_pipeline/data/` — update per your environment
