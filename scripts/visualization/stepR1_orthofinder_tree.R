#!/usr/bin/env Rscript
# [R1] Classic phylogenetic tree visualization from OrthoFinder results
#
# Usage:
#   Rscript stepR1_orthofinder_tree.R [--config config.yaml]
#
# Dependencies:
#   - ape, phytools, ggtree, ggplot2, optparse, yaml

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(ggtree)
  library(ggplot2)
  library(optparse)
  library(yaml)
  library(dplyr)
})

option_list <- list(
  make_option(c('--config', '-c'), help='YAML config file'),
  make_option(c('--tree'), help='Newick tree file'),
  make_option(c('--outdir', '-o'), default='output/tree_viz',
              help='Output directory'),
  make_option(c('--format'), default='pdf',
              help='Output format (pdf, svg, png)'),
  make_option(c('--width'), type='integer', default=12, help='Figure width'),
  make_option(c('--height'), type='integer', default=10, help='Figure height'),
  make_option(c('--layout'), default='rectangular',
              help='Tree layout (rectangular, circular, fan, slanted)'),
  make_option(c('--tip-labels'), type='logical', default=TRUE,
              help='Show tip labels'),
  make_option(c('--bootstrap-threshold'), type='double', default=70,
              help='Minimum bootstrap value to display')
)

opt <- parse_args(OptionParser(option_list = option_list))

# Load config
cfg <- list()
if (!is.null(opt$config) && file.exists(opt$config)) {
  cfg <- yaml.load_file(opt$config)
}

# Resolve parameters
tree_file <- opt$tree %||% cfg$tree_file %||% stop('Missing --tree')
outdir <- opt$outdir %||% cfg$outdir %||% 'output/tree_viz'
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Read tree
tr <- read.tree(tree_file)
cat(sprintf('[R1] Tree loaded: %d tips, %d nodes\n',
            Ntip(tr), Nnode(tr)))

# Root at midpoint if unrooted
if (!is.rooted(tr)) {
  tr <- midpoint.root(tr)
}

# Create base plot
p <- ggtree(tr, layout = opt$layout, size = 0.5)

# Add bootstrap support (node labels)
if (!is.null(tr$node.label)) {
  bs <- as.numeric(tr$node.label)
  bs[is.na(bs)] <- 0
  # Only show values above threshold
  bs_label <- ifelse(bs >= opt$bootstrap_threshold,
                     round(bs), '')
  # Label internal nodes
  internal_nodes <- (Ntip(tr) + 1):(Ntip(tr) + Nnode(tr))
  bs_df <- data.frame(node = internal_nodes, label = bs_label)
  bs_df <- bs_df[bs_df$label != '', ]
  if (nrow(bs_df) > 0) {
    p <- p %<+% bs_df +
      geom_label(aes(label = label), fill = 'white', size = 2.5,
                 label.size = 0.2, na.rm = TRUE)
  }
}

# Add tip labels
if (opt$tip_labels) {
  p <- p + geom_tiplab(size = 3, align = TRUE, linesize = 0.3)
}

# Theme
p <- p + theme_tree2() +
  theme(plot.margin = margin(10, 10, 10, 10))

# Save
ext <- opt$format
out_path <- file.path(outdir, paste0('species_tree.', ext))
ggsave(out_path, p, width = opt$width, height = opt$height,
       dpi = 300, limitsize = FALSE)
cat(sprintf('[R1] Tree saved: %s\n', out_path))

# Also write a basic summary
tree_stats <- data.frame(
  Metric = c('N_tips', 'N_nodes', 'Tree_length', 'Is_rooted',
             'Tip_labels'),
  Value = c(Ntip(tr), Nnode(tr),
            round(sum(tr$edge.length), 4),
            is.rooted(tr),
            paste(tr$tip.label[1:min(5, Ntip(tr))], collapse = ', '))
)
write.csv(tree_stats, file.path(outdir, 'tree_stats.csv'),
          row.names = FALSE)
cat('[R1] Done\n')
