#!/usr/bin/env Rscript
# [R2] Phylogenetic tree + orthogroup presence/absence heatmap
#
# Generates a combined plot: species tree + heatmap of OG presence/absence.
#
# Usage:
#   Rscript stepR2_orthofinder_heatmap.R [--config config.yaml]
#
# Dependencies:
#   - ape, ggtree, ggplot2, reshape2, optparse, yaml, RColorBrewer

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(reshape2)
  library(optparse)
  library(yaml)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
})

option_list <- list(
  make_option(c('--config', '-c'), help='YAML config file'),
  make_option(c('--tree'), help='Newick species tree file'),
  make_option(c('--orthogroups'), help='Orthogroups.tsv/csv file'),
  make_option(c('--outdir', '-o'), default='output/heatmap',
              help='Output directory'),
  make_option(c('--format'), default='pdf',
              help='Output format (pdf, svg, png)'),
  make_option(c('--width'), type='integer', default=18, help='Figure width'),
  make_option(c('--height'), type='integer', default=12, help='Figure height'),
  make_option(c('--top-ogs'), type='integer', default=50,
              help='Number of top OGs to display'),
  make_option(c('--min-species'), type='integer', default=2,
              help='Minimum species an OG must appear in'),
  make_option(c('--colors'), default='YlOrRd',
              help='RColorBrewer palette name')
)

opt <- parse_args(OptionParser(option_list = option_list))

# Load config
cfg <- list()
if (!is.null(opt$config) && file.exists(opt$config)) {
  cfg <- yaml.load_file(opt$config)
}

# Resolve parameters
tree_file <- opt$tree %||% cfg$tree_file %||% stop('Missing --tree')
og_file <- opt$orthogroups %||% cfg$orthogroups %||% stop('Missing --orthogroups')
outdir <- opt$outdir %||% cfg$outdir %||% 'output/heatmap'
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
top_ogs <- opt$top_ogs %||% cfg$top_ogs %||% 50
min_species <- opt$min_species %||% cfg$min_species %||% 2

# 1. Read tree
tr <- read.tree(tree_file)
if (!is.rooted(tr)) {
  tr <- midpoint.root(tr)
}
cat(sprintf('[R2] Tree: %d tips\n', Ntip(tr)))

# 2. Read orthogroups (presence/absence matrix)
raw <- read.delim(og_file, header = TRUE, check.names = FALSE,
                  stringsAsFactors = FALSE, sep = '\t')
cat(sprintf('[R2] Orthogroups: %d rows x %d cols\n', nrow(raw), ncol(raw)))

# Create binary presence/absence matrix
og_names <- raw[, 1]
mat <- raw[, -1, drop = FALSE]
mat[] <- lapply(mat, function(x) as.integer(!is.na(x) & nchar(trimws(x)) > 0))
rownames(mat) <- og_names

# Filter: keep OGs present in at least min_species species
keep <- rowSums(mat) >= min_species
mat <- mat[keep, , drop = FALSE]
cat(sprintf('[R2] After filtering: %d OGs\n', nrow(mat)))

# Select top OGs by total presence
if (nrow(mat) > top_ogs) {
  presence_order <- order(rowSums(mat), decreasing = TRUE)
  mat <- mat[presence_order[1:top_ogs], , drop = FALSE]
}

# Ensure columns match tree tip labels
common_species <- intersect(colnames(mat), tr$tip.label)
if (length(common_species) == 0) {
  stop('No matching species between OG matrix and tree tip labels')
}
mat <- mat[, common_species, drop = FALSE]

# Reorder matrix by tree tip order
tip_order <- rev(tr$tip.label[tr$tip.label %in% common_species])
mat <- mat[, tip_order, drop = FALSE]

# 3. Create ggtree plot with heatmap
p <- ggtree(tr, layout = 'rectangular', size = 0.5) +
  geom_tiplab(size = 3, offset = 0.02) +
  xlim(NA, max(tr$edge.length) * 1.3) +
  theme(plot.margin = margin(5, 5, 5, 5))

# Melt matrix for heatmap
mat_long <- melt(as.matrix(mat), varnames = c('OG', 'Species'),
                 value.name = 'Presence')

# Add heatmap
p2 <- gheatmap(p, mat_long, offset = max(tr$edge.length) * 0.1,
               width = 0.8, colnames = FALSE,
               color = 'gray90') +
  scale_fill_gradientn(
    colors = brewer.pal(9, opt$colors),
    name = 'Present',
    breaks = c(0, 1),
    labels = c('Absent', 'Present')
  ) +
  theme(legend.position = 'right',
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        plot.margin = margin(10, 10, 10, 10))

# Save
ext <- opt$format
out_path <- file.path(outdir, paste0('orthogroup_heatmap.', ext))
ggsave(out_path, p2, width = opt$width, height = opt$height,
       dpi = 300, limitsize = FALSE)
cat(sprintf('[R2] Heatmap saved: %s\n', out_path))

# Write filtered matrix
write.csv(mat, file.path(outdir, 'orthogroup_presence_matrix.csv'))
cat('[R2] Done\n')
