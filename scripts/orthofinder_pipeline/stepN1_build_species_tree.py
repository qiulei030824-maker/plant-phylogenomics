#!/usr/bin/env python3
"""[N1] Build Species Tree — Industrial-grade.

Steps:
  1. Gathers single-copy orthogroups from OrthoFinder output
  2. Concatenates protein alignments
  3. Builds species tree with IQ-TREE2

Usage:
    python stepN1_build_species_tree.py [--config config.yaml]
"""

import os
import sys
import json
import yaml
import shutil
import argparse
import logging
import subprocess
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)-7s | %(message)s')
logger = logging.getLogger(__name__)


def load_config(path: str) -> dict:
    with open(path) as f:
        return yaml.safe_load(f) or {}


def build_parser():
    p = argparse.ArgumentParser(description='[N1] Build Species Tree from OrthoFinder results')
    p.add_argument('--config', '-c', help='YAML config')
    p.add_argument('--orthofinder-dir', help='OrthoFinder results directory (OrthoFinder/Results_XXX)')
    p.add_argument('--outdir', '-o', default='output/species_tree')
    p.add_argument('--threads', '-t', type=int, default=0)
    return p


def find_single_copy_ogs(og_dir: Path) -> list:
    """Identify single-copy orthogroups from Orthogroups.csv or Orthogroups.tsv"""
    # Try Orthogroups.csv first (newer OrthoFinder), else .tsv
    csv = og_dir / 'Orthogroups.csv'
    tsv = og_dir / 'Orthogroups.tsv'
    path = csv if csv.exists() else tsv
    if not path.exists():
        logger.error(f'No orthogroups file found in {og_dir}')
        sys.exit(1)

    single_copies = []
    with open(path) as f:
        header = f.readline().strip().split('\t')
        n_species = len(header) - 1
        for line in f:
            parts = line.strip().split('\t')
            genes_per_sp = [len(p.split()) if p else 0 for p in parts[1:]]
            if all(c == 1 for c in genes_per_sp):
                single_copies.append(parts[0])
    logger.info(f'Found {len(single_copies)} single-copy orthogroups')
    return single_copies


def main():
    parser = build_parser()
    args = parser.parse_args()

    cfg = load_config(args.config) if args.config else {}
    orthofinder_dir = Path(args.orthofinder_dir or cfg.get('orthofinder_dir', 'OrthoFinder'))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    threads = args.threads or os.cpu_count() or 4

    # 1. Find single-copy OGs
    og_dir = orthofinder_dir / 'Orthogroups'
    if not og_dir.exists():
        # Try Results_<date>
        results = sorted(orthofinder_dir.glob('Results_*'))
        if results:
            og_dir = results[-1] / 'Orthogroups'
    single_copies = find_single_copy_ogs(og_dir)

    # 2. Concatenate alignments
    aln_dir = orthofinder_dir / 'MultipleSequenceAlignments'
    concat_path = outdir / 'concatenated_aln.fasta'
    with open(concat_path, 'w') as out:
        for og in single_copies:
            fa = aln_dir / f'{og}.fa'
            if fa.exists():
                with open(fa) as f:
                    out.write(f.read())
    logger.info(f'Concatenated alignment: {concat_path}')

    # 3. Build tree with IQ-TREE2
    tree_path = outdir / 'species_tree.treefile'
    cmd = ['iqtree2', '-s', str(concat_path), '--prefix', str(outdir / 'species_tree'),
           '-T', str(threads), '-m', 'MFP', '-B', '1000']
    logger.info(f'Running IQ-TREE2: {" ".join(cmd)}')
    subprocess.run(cmd, check=True)
    logger.info(f'Species tree: {tree_path}')


if __name__ == '__main__':
    main()
