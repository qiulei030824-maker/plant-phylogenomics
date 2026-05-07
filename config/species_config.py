#!/usr/bin/env python3
"""
Shared configuration for C2 domain / Pfam domain analysis pipeline.

This file centralizes:
  - Species groups and colors
  - Base paths
  - Utility functions (get_species, get_protein_id, get_group, parse_newick_tips)

All scripts should import from here instead of duplicating these definitions.
"""

import os
import re
from pathlib import Path

# ── Base paths ─────────────────────────────────────────────────────────────
DATA_DIR = "/data5/qiulei/PPI/data"
BASE_ROOT = Path("/data5/qiulei/PPI/03.OrthoFinder_MCScanX/0430")
SKIP_DIRS = {"logs", "pfamdb", "__pycache__", "taxonomy_classification"}

# ── Species groups (major plant clades) ────────────────────────────────────
SPECIES_GROUP = {
    "Cucurbits": [
        "97103", "Cargyrosperma", "Cmaxima_Rimu", "Cmoschata_Rifu",
        "Cpepo", "ChineseLong", "DHL92", "USVL1VR-Ls",
    ],
    "Brassicaceae": [
        "arabidopsis_halleri", "arabidopsis_lyrata", "arabidopsis_thaliana",
        "brassica_napus", "brassica_oleracea", "brassica_rapa",
        "eutrema_salsugineum", "arabis_alpina",
    ],
    "Solanaceae": [
        "capsicum_annuum", "nicotiana_attenuata",
        "solanum_lycopersicum", "solanum_tuberosum",
    ],
    "Fabids": [
        "corylus_avellana", "quercus_lobata", "eucalyptus_grandis",
        "pistacia_vera", "gossypium_raimondii", "theobroma_cacao",
        "ficus_carica", "malus_domestica_golden", "rosa_chinensis",
        "glycine_max", "medicago_truncatula", "phaseolus_vulgaris",
        "citrus_clementina", "manihot_esculenta", "populus_trichocarpa",
        "prunus_persica",
    ],
    "Vitales": ["vitis_vinifera"],
    "Monocots": [
        "asparagus_officinalis", "ananas_comosus",
        "brachypodium_distachyon", "sorghum_bicolor", "zea_mays",
        "oryza_sativa", "musa_acuminata",
    ],
    "Basal_Angiosperms": ["amborella_trichopoda", "nymphaea_colorata"],
    "Non_Vascular": [
        "marchantia_polymorpha", "physcomitrium_patens",
        "selaginella_moellendorffii",
    ],
}

# ── Group colors (for plotting) ────────────────────────────────────────────
GROUP_COLORS = {
    "Cucurbits": "#E41A1C",
    "Brassicaceae": "#377EB8",
    "Solanaceae": "#4DAF4A",
    "Fabids": "#984EA3",
    "Vitales": "#FF7F00",
    "Monocots": "#FFFF33",
    "Basal_Angiosperms": "#A65628",
    "Non_Vascular": "#F781BF",
    "Other": "#999999",
}

# ── C2-specific subdirectory names (used when PFAM_ID == "C2") ────────────
C2_PEP_DIRNAME = "C2domain_pep"
C2_CDS_DIRNAME = "C2domain_cds"
C2_HMMER_DIRNAME = "hmmer"
C2_TREE_DIRNAME = "tree"
C2_ALIGN_DIRNAME = "algin"
C2_HMMER_FILENAME = "all_C2.domtblout"
C2_TREE_FILENAME = "C2_tree.treefile"


# ── Utility functions ──────────────────────────────────────────────────────

# ── Taxonomy hierarchy (built from taxonomy_classification dir) ────────────
TAXONOMY_HIERARCHY = {
    "Rosids": {
        "Cucurbitaceae": {
            "Cucumis": ["ChineseLong", "DHL92"],
            "Lagenaria": ["USVL1VR-Ls"],
            "Citrullus": ["97103"],
            "Cucurbita": ["Cargyrosperma", "Cmaxima_Rimu", "Cmoschata_Rifu", "Cpepo"],
        },
        "Brassicaceae": {
            "Arabidopsis": ["arabidopsis_halleri", "arabidopsis_lyrata", "arabidopsis_thaliana"],
            "Brassica": ["brassica_napus", "brassica_oleracea", "brassica_rapa"],
            "Eutrema": ["eutrema_salsugineum"],
            "Arabis": ["arabis_alpina"],
        },
        "Myrtaceae": {"Eucalyptus": ["eucalyptus_grandis"]},
        "Euphorbiaceae": {"Manihot": ["manihot_esculenta"]},
        "Salicaceae": {"Populus": ["populus_trichocarpa"]},
        "Fagaceae": {"Quercus": ["quercus_lobata"]},
        "Betulaceae": {"Corylus": ["corylus_avellana"]},
        "Malvaceae": {"Theobroma": ["theobroma_cacao"], "Gossypium": ["gossypium_raimondii"]},
        "Rutaceae": {"Citrus": ["citrus_clementina"]},
        "Anacardiaceae": {"Pistacia": ["pistacia_vera"]},
        "Vitaceae": {"Vitis": ["vitis_vinifera"]},
        "Rosaceae": {"Malus": ["malus_domestica_golden"], "Rosa": ["rosa_chinensis"], "Prunus": ["prunus_persica"]},
        "Moraceae": {"Ficus": ["ficus_carica"]},
        "Fabaceae": {"Glycine": ["glycine_max"], "Phaseolus": ["phaseolus_vulgaris"], "Medicago": ["medicago_truncatula"]},
    },
    "Asterids": {
        "Solanaceae": {
            "Solanum": ["solanum_lycopersicum", "solanum_tuberosum"],
            "Nicotiana": ["nicotiana_attenuata"],
            "Capsicum": ["capsicum_annuum"],
        },
        "Asteraceae": {},
    },
    "Monocots": {
        "Asparagaceae": {"Asparagus": ["asparagus_officinalis"]},
        "Musaceae": {"Musa": ["musa_acuminata"]},
        "Poaceae": {
            "Oryza": ["oryza_sativa"],
            "Zea": ["zea_mays"],
            "Brachypodium": ["brachypodium_distachyon"],
            "Sorghum": ["sorghum_bicolor"],
        },
        "Bromeliaceae": {"Ananas": ["ananas_comosus"]},
    },
    "Basal_Angiosperms": {
        "Nymphaeaceae": {"Nymphaea": ["nymphaea_colorata"]},
        "Amborellaceae": {"Amborella": ["amborella_trichopoda"]},
    },
    "Bryophytes": {
        "Funariaceae": {"Physcomitrium": ["physcomitrium_patens"]},
        "Marchantiaceae": {"Marchantia": ["marchantia_polymorpha"]},
    },
    "Lycophytes": {
        "Selaginellaceae": {"Selaginella": ["selaginella_moellendorffii"]},
    },
}


def get_all_species():
    """Return a flat list of all species names."""
    result = []
    for clade in TAXONOMY_HIERARCHY.values():
        for family in clade.values():
            for genus in family.values():
                result.extend(genus)
    return result


def get_species(tip_label):
    """Extract species name from a tree tip label (format: 'species|protein_id')."""
    return tip_label.split("|")[0]


def get_protein_id(tip_label):
    """Extract protein ID from a tree tip label (format: 'species|protein_id')."""
    return tip_label.split("|")[1] if "|" in tip_label else tip_label


def get_group(species_name):
    """Return the group/clade name for a given species name."""
    for clade_name, clade_data in TAXONOMY_HIERARCHY.items():
        for family_data in clade_data.values():
            for genus_species in family_data.values():
                if species_name in genus_species:
                    return clade_name
    return "Other"


def get_genus(species_name):
    """Return the genus name for a given species name."""
    for clade_data in TAXONOMY_HIERARCHY.values():
        for family_data in clade_data.values():
            for genus, species_list in family_data.items():
                if species_name in species_list:
                    return genus
    return None


def get_family(species_name):
    """Return the family name for a given species name."""
    for clade_data in TAXONOMY_HIERARCHY.values():
        for family, family_data in clade_data.items():
            for genus_species in family_data.values():
                if species_name in genus_species:
                    return family
    return None


def get_group_color(species_name):
    """Return the color for a given species name."""
    return GROUP_COLORS.get(get_group(species_name), GROUP_COLORS["Other"])


def parse_newick_tips(filepath):
    """Extract tip labels from a Newick tree file, preserving order."""
    with open(filepath) as f:
        newick = f.read().strip()
    cleaned = re.sub(r'\)\d+:', '):', newick)
    cleaned = re.sub(r':[\d.eE+-]+', '', cleaned)
    matches = re.findall(r'[,(](\S+?)[,)]', cleaned)
    seen = set()
    tips = []
    for m in matches:
        m = m.strip()
        if m and not m.startswith('(') and m not in seen:
            seen.add(m)
            tips.append(m)
    return tips


def resolve_pfam_paths(pfam_id):
    """
    Resolve input/output paths for a given PFAM_ID.

    For C2 domain, data is stored directly under BASE_ROOT.
    For other PFAM IDs, data is stored under BASE_ROOT/{pfam_id}/.

    Returns a dict with keys:
        base_dir, pep_dir, cds_dir, hmmer_dir, tree_dir, align_dir,
        hmmer_file, tree_file
    """
    is_c2 = pfam_id.upper() == "C2"

    if is_c2:
        base_dir = BASE_ROOT
        pep_dir = BASE_ROOT / C2_PEP_DIRNAME
        cds_dir = BASE_ROOT / C2_CDS_DIRNAME
        hmmer_dir = BASE_ROOT / C2_HMMER_DIRNAME
        tree_dir = BASE_ROOT / C2_TREE_DIRNAME
        align_dir = BASE_ROOT / C2_ALIGN_DIRNAME
        hmmer_file = hmmer_dir / C2_HMMER_FILENAME
        tree_file = tree_dir / C2_TREE_FILENAME
    else:
        base_dir = BASE_ROOT / pfam_id
        pep_dir = base_dir / "pep"
        cds_dir = base_dir / "cds"
        hmmer_dir = base_dir / "hmmer"
        tree_dir = base_dir / "tree"
        align_dir = base_dir / "algin"
        hmmer_file = hmmer_dir / f"all_{pfam_id}.domtblout"
        tree_file = tree_dir / f"{pfam_id}_tree.treefile"

    return {
        "base_dir": base_dir,
        "pep_dir": pep_dir,
        "cds_dir": cds_dir,
        "hmmer_dir": hmmer_dir,
        "tree_dir": tree_dir,
        "align_dir": align_dir,
        "hmmer_file": hmmer_file,
        "tree_file": tree_file,
        "is_c2": is_c2,
    }
