#!/usr/bin/env Rscript
"""
[R1] OrthoFinder Tree Visualization with Group Coloring

Usage:
    Rscript stepR1_orthofinder_tree.R --tree <newick_file> --output-dir <dir> [options]

Requires:
    R packages: ape, ggtree, ggplot2, RColorBrewer
"""

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-t", "--tree"), type="character", help="Input newick tree file"),
  make_option(c("-o", "--output-dir"), type="character", default=".", help="Output directory"),
  make_option(c("--og-name"), type="character", default="orthogroup", help="Orthogroup name for plot title"),
  make_option(c("--width"), type="numeric", default=12, help="Plot width in inches"),
  make_option(c("--height"), type="numeric", default=8, help="Plot height in inches"),
  make_option(c("--tip-size"), type="numeric", default=3, help="Tip label font size")
)

parser <- OptionParser(option_list=option_list, description="R1: OrthoFinder Tree Visualization")
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
  "N_nucifera" = "LotusSacred",
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
    matched <- FALSE
    for (prefix in names(prefix_map)) {
      if (startsWith(label, prefix) || grepl(paste0("^", prefix, "|"), label)) {
        return(prefix_map[[prefix]])
      }
    }
    # Try to extract the first part before | or _
    parts <- strsplit(label, "[|_]")
    if (length(parts[[1]]) > 0) {
      return(parts[[1]][1])
    }
    return("Unknown")
  })
}

# Group definitions
melon_cuke_group <- c("Cucumber", "Watermelon", "Melon", "Pumpkin", "AquaticXY", "Bladderwort")
brassica_group <- c("Arabidopsis", "Cabbage", "ChineseCabbage", "ShepherdsPurse", "SaltCress", "Parvula", "SpiderFlower")
legume_group <- c("Soybean", "Lotus", "Medicago", "Chickpea", "Lupin", "Cowpea", "AdzukiBean", "MungBean", "RedClover")
solanaceae_group <- c("Tobacco", "Tomato", "Potato", "Pepper")
grass_group <- c("Rice", "Maize", "Sorghum", "FoxtailMillet", "GreenMillet", "Switchgrass", "HallPanicgrass", "Timothy", "Orchardgrass", "WeepingLovegrass")
rosid_group <- c("Peach", "Apple", "Strawberry", "Jujube", "Poplar", "Cottonwood", "Cassava", "CastorBean", "Jatropha", "Cocoa", "Cotton", "Eucalyptus", "Grape", "Orange", "Clementine", "Papaya", "Avocado", "Olive", "Pomegranate", "ValleyOak", "RedOak", "BartlettOak", "SpiderFlower", "Pinnata")

# Color palette
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
  if (sp %in% melon_cuke_group) return("Cucurbitaceae")
  if (sp %in% brassica_group) return("Brassicaceae")
  if (sp %in% legume_group) return("Fabaceae")
  if (sp %in% solanaceae_group) return("Solanaceae")
  if (sp %in% grass_group) return("Poaceae")
  if (sp %in% rosid_group) return("Rosid")
  return("Other")
}

# ── Load libraries ─────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(RColorBrewer)
})

# ── Read tree ─────────────────────────────────────────────────────────────
tree <- read.tree(tree_path)
cat(sprintf("  Tree loaded: %d tips\n", length(tree$tip.label)))

# ── Map tips ───────────────────────────────────────────────────────────────
tip_species <- build_species_map(tree$tip.label, default_prefix_map)
tip_groups <- sapply(tip_species, get_group)

# ── Build annotation data ─────────────────────────────────────────────────
group_df <- data.frame(
  label = tree$tip.label,
  Species = tip_species,
  Group = factor(tip_groups, levels = names(group_colors)),
  stringsAsFactors = FALSE
)

# ── Plot ───────────────────────────────────────────────────────────────────
tip_labels_formatted <- paste0(tip_species, " | ", tree$tip.label)
names(tip_labels_formatted) <- tree$tip.label

p <- ggtree(tree, layout = "rectangular", branch.length = "branch.length") %<+% group_df +
  geom_tiplab(aes(color = Group, label = tip_labels_formatted[label]),
              size = tip_size, fontface = "italic") +
  scale_color_manual(values = group_colors, name = "Group") +
  geom_treescale(x = 0, y = 0, width = 0.1, offset = 2) +
  labs(title = paste(og_name, "- OrthoFinder Tree")) +
  theme_tree(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# ── Save ───────────────────────────────────────────────────────────────────
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
output_pdf <- file.path(output_dir, paste0(og_name, "_orthofinder_tree.pdf"))
ggsave(output_pdf, plot = p, width = width, height = height, dpi = 300)
cat(sprintf("  Output: %s\n", output_pdf))
cat("R1 complete!\n")