#!/usr/bin/env Rscript
"""
[R2] OrthoFinder Tree + Gene Count Heatmap

Usage:
    Rscript stepR2_orthofinder_heatmap.R --tree <newick_file> --output-dir <dir> [options]

Requires:
    R packages: ape, ggtree, ggplot2, RColorBrewer, patchwork, reshape2
"""

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-t", "--tree"), type="character", help="Input newick tree file"),
  make_option(c("-o", "--output-dir"), type="character", default=".", help="Output directory"),
  make_option(c("--og-name"), type="character", default="orthogroup", help="Orthogroup name"),
  make_option(c("--width"), type="numeric", default=14, help="Plot width in inches"),
  make_option(c("--height"), type="numeric", default=10, help="Plot height in inches"),
  make_option(c("--tip-size"), type="numeric", default=2.5, help="Tip label font size")
)

parser <- OptionParser(option_list=option_list, description="R2: Tree + Heatmap")
args <- parse_args(parser)

if (is.null(args$tree)) {
  print_help(parser)
  quit(status=1)
}

tree_path <- args$tree
output_dir <- args$output_dir
og_name <- args$`og-name`
width <- args$width
height <- args$height
tip_size <- args$`tip-size`

# ── Species prefix mapping ────────────────────────────────────────────────
default_prefix_map <- list(
  "A_fistulosum" = "Fistulosum",
  "A_trichopoda" = "Amborella",
  "A_thaliana" = "Arabidopsis",
  "B_oleracea" = "Cabbage",
  "B_rapa" = "ChineseCabbage",
  "B_vulgaris" = "Beet",
  "C_arietinum" = "Chickpea",
  "C_clementina" = "Clementine",
  "C_lanatus" = "Watermelon",
  "C_melo" = "Melon",
  "C_papaya" = "Papaya",
  "C_pepo" = "Pumpkin",
  "C_rubella" = "ShepherdsPurse",
  "C_sativus" = "ChineseLong",
  "C_sinensis" = "Orange",
  "C_tinctorius" = "Safflower",
  "D_carinota" = "Carrot",
  "D_glomerata" = "Orchardgrass",
  "E_curvula" = "WeepingLovegrass",
  "E_grandis" = "Eucalyptus",
  "E_salsugineum" = "SaltCress",
  "G_arboreum" = "TreeCotton",
  "G_barbadense" = "PimaCotton",
  "G_hirsutum" = "UplandCotton",
  "G_max" = "Soybean",
  "G_raimondii" = "Raimondii",
  "H_annus" = "Sunflower",
  "H_brachypoda" = "AquaticXY",
  "H_sapiens" = "Human",
  "J_curcas" = "Jatropha",
  "L_angustifolius" = "Lupin",
  "L_japonicus" = "Lotus",
  "L_sativa" = "Lettuce",
  "M_acuminata" = "Banana",
  "M_domestica" = "Apple",
  "M_esculenta" = "Cassava",
  "M_guttatus" = "Monkeyflower",
  "M_polymorpha" = "Liverwort",
  "M_truncatula" = "Medicago",
  "N_attenuata" = "Tobacco",
  "N_nucifera" = "SacredLotus",
  "O_bartlettii" = "BartlettOak",
  "O_europaea" = "Olive",
  "O_sativa" = "Rice",
  "P_abies" = "Spruce",
  "P_americana" = "Avocado",
  "P_dactylifera" = "DatePalm",
  "P_deltoides" = "Cottonwood",
  "P_granatum" = "Pomegranate",
  "P_hallii" = "HallPanicgrass",
  "P_patens" = "Moss",
  "P_persica" = "Peach",
  "P_trichocarpa" = "Poplar",
  "P_virgatum" = "Switchgrass",
  "Q_lobata" = "ValleyOak",
  "Q_rubra" = "RedOak",
  "R_communis" = "CastorBean",
  "S_bicolor" = "Sorghum",
  "S_fallax" = "Timothy",
  "S_italica" = "FoxtailMillet",
  "S_lycopersicum" = "Tomato",
  "S_moellendorffii" = "Spikemoss",
  "S_parvula" = "Parvula",
  "S_pinnata" = "Pinnata",
  "S_polyrhiza" = "Duckweed",
  "S_tuberosum" = "Potato",
  "S_viridis" = "GreenMillet",
  "T_cacao" = "Cocoa",
  "T_hassleriana" = "SpiderFlower",
  "T_palmeri" = "PalmerAmaranth",
  "T_pratense" = "RedClover",
  "T_salsuginea" = "SaltCress",
  "U_gibba" = "Bladderwort",
  "V_angularis" = "AdzukiBean",
  "V_radiata" = "MungBean",
  "V_unguiculata" = "Cowpea",
  "V_vinifera" = "Grape",
  "X_irriguum" = "Iris",
  "Z_jujuba" = "Jujube",
  "Z_mays" = "Maize",
  "Z_mays_AGPv4" = "Maize_AGPv4",
  "e_grandis" = "Eucalyptus"
)

build_species_map <- function(tip_labels, prefix_map) {
  sapply(tip_labels, function(label) {
    for (prefix in names(prefix_map)) {
      if (grepl(paste0("^", prefix, "[|_]"), label) || label == prefix) {
        return(prefix_map[[prefix]])
      }
    }
    parts <- strsplit(label, "[|_]")[[1]]
    return(parts[1])
  })
}

melon_cuke_group <- c("ChineseLong", "Watermelon", "Melon", "Pumpkin", "AquaticXY", "Bladderwort")
brassica_group <- c("Arabidopsis", "Cabbage", "ChineseCabbage", "ShepherdsPurse", "SaltCress", "Parvula", "SpiderFlower")
legume_group <- c("Soybean", "Lotus", "Medicago", "Chickpea", "Lupin", "Cowpea", "AdzukiBean", "MungBean", "RedClover")
solanaceae_group <- c("Tobacco", "Tomato", "Potato")
grass_group <- c("Rice", "Maize", "Sorghum", "FoxtailMillet", "GreenMillet", "Switchgrass", "HallPanicgrass", "Timothy", "Orchardgrass", "WeepingLovegrass")
rosid_group <- c("Peach", "Apple", "Jujube", "Poplar", "Cottonwood", "Cassava", "CastorBean", "Jatropha", "Cocoa", "Eucalyptus", "Grape", "Orange", "Clementine", "Papaya", "Avocado", "Olive", "Pomegranate", "ValleyOak", "RedOak", "BartlettOak", "Pinnata")

group_colors <- c(
  "Cucurbitaceae" = "#E41A1C",
  "Brassicaceae" = "#377EB8",
  "Fabaceae" = "#4DAF4A",
  "Solanaceae" = "#984EA3",
  "Poaceae" = "#FF7F00",
  "Rosid" = "#FFFF33",
  "Outgroup" = "#A65628",
  "Other" = "#999999"
)

get_group <- function(sp) {
  for (g in names(group_defs)) {
    if (sp %in% group_defs[[g]]) return(g)
  }
  return("Other")
}

# ── Load libraries ─────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)
})

# ── Read tree ─────────────────────────────────────────────────────────────
tree <- read.tree(tree_path)
cat(sprintf("  Tree loaded: %d tips\n", length(tree$tip.label)))

# ── Map tips ───────────────────────────────────────────────────────────────
tip_species <- build_species_map(tree$tip.label, default_prefix_map)
tip_groups <- sapply(tip_species, get_group)

# ── Count genes per species ────────────────────────────────────────────────
species_counts <- table(tip_species)
species_stats <- data.frame(
  species = names(species_counts),
  gene_count = as.integer(species_counts),
  stringsAsFactors = FALSE
)
cat(sprintf("  Species: %d\n", nrow(species_stats)))

# ── Build ordered data ─────────────────────────────────────────────────────
tree_plot <- ggtree(tree, layout = "rectangular", branch.length = "branch.length")
tree_data <- tree_plot$data
tip_order <- rev(tree_data$label[tree_data$isTip])

heatmap_df <- data.frame(
  label = tip_order,
  species = tip_species[tip_order],
  stringsAsFactors = FALSE
)
heatmap_df <- merge(heatmap_df, species_stats, by = "species", all.x = TRUE, sort = FALSE)
heatmap_df$gene_count[is.na(heatmap_df$gene_count)] <- 0

heatmap_long <- reshape2::melt(heatmap_df,
                                id.vars = c("label", "species"),
                                measure.vars = c("gene_count"),
                                variable.name = "metric", value.name = "value")
heatmap_long$label <- factor(heatmap_long$label, levels = rev(tip_order))

group_df <- data.frame(
  label = tree$tip.label,
  Group = factor(tip_groups, levels = names(group_colors)),
  stringsAsFactors = FALSE
)

# ── Tree plot ──────────────────────────────────────────────────────────────
tip_labels_formatted <- paste0(tip_species, " | ", tree$tip.label)
names(tip_labels_formatted) <- tree$tip.label

p_tree <- ggtree(tree, layout = "rectangular", branch.length = "branch.length") %<+% group_df +
  geom_tiplab(aes(color = Group, label = tip_labels_formatted[label]),
              size = tip_size, fontface = "italic") +
  scale_color_manual(values = group_colors, name = "Group") +
  geom_nodepoint(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) >= 50),
                 size = 1.5, shape = 21, fill = "white", color = "grey30") +
  labs(title = paste(og_name, "OrthoFinder Tree with Gene Counts")) +
  theme_tree(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  geom_treescale(x = 0, y = 0, width = 0.1, offset = 2)

# ── Heatmap plot ───────────────────────────────────────────────────────────
heatmap_colors <- colorRampPalette(brewer.pal(9, "YlOrRd"))(100)

p_heat <- ggplot(heatmap_long, aes(x = metric, y = label, fill = value)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradientn(colors = heatmap_colors, name = "Genes") +
  scale_x_discrete(position = "top", labels = c("Gene Count")) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 8),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

# ── Save ───────────────────────────────────────────────────────────────────
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_pdf <- file.path(output_dir, paste0(og_name, "_tree_heatmap.pdf"))

combined <- p_tree + p_heat + plot_layout(widths = c(3, 1))
ggsave(output_pdf, plot = combined, width = width, height = height, dpi = 300)
cat(sprintf("  Output: %s\n", output_pdf))
cat("R2 complete!\n")