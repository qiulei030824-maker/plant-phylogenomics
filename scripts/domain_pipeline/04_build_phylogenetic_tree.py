#!/usr/bin/env python3
"""
[Step 4] Build phylogenetic tree of domain-containing proteins — Industrial-grade.

Features:
  - Profile system (fast|standard|accurate|ultra)
  - Sequence selection strategies (all|longest|canonical|domain_best|longest_isoform|representative)
  - Auto strategy engine (auto-detects <50 / 50-200 / 200-1000 / >1000 sequences)
  - MAFFT alignment (auto|linsi|einsi|fftns2) + trimming (clipkit|trimal|none)
  - IQ-TREE2 with MFP, UFBoot2, or fast LG
  - Taxonomy filtering via species_config (--clade / --family / --genus / --species)
  - Resume/checkpoint system
  - Metadata + QC report + auto-generated methods paragraph
  - Tree visualization (PDF/SVG/PNG)
  - Config file support (--config config.yaml)
  - Full logging

Usage:
    python 04_build_phylogenetic_tree.py <PFAM_ID>
    python 04_build_phylogenetic_tree.py <PFAM_ID> --profile accurate
    python 04_build_phylogenetic_tree.py <PFAM_ID> --strategy longest --species ChineseLong,DHL92
    python 04_build_phylogenetic_tree.py <PFAM_ID> --config config.yaml
"""

import os, sys, re, json, time, shutil, textwrap, subprocess, argparse
from pathlib import Path
from collections import defaultdict
from datetime import datetime

# ── Path setup ──────────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)

try:
    from config.species_config import resolve_pfam_paths, TAXONOMY_HIERARCHY, get_all_species, get_group
except ModuleNotFoundError:
    from species_config import resolve_pfam_paths, TAXONOMY_HIERARCHY, get_all_species, get_group

# ── Logging ─────────────────────────────────────────────────────────────────
import logging

def setup_logger(log_dir, name="04_build_domain_tree", level=logging.INFO):
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.handlers.clear()

    fh = logging.FileHandler(log_dir / "run.log", mode="a")
    fh.setLevel(level)
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logger.addHandler(fh)

    ch = logging.StreamHandler(sys.stdout)