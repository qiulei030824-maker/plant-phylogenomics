#!/usr/bin/env python3
"""[N2] Build Ortholog Tree — Industrial-grade.

Builds a gene tree for a specified orthogroup from OrthoFinder results,
using IQ-TREE2 with appropriate model selection.

Usage:
    python stepN2_build_ortholog_tree.py --og OG0000001
"""

import os
import sys
import json
import yaml
import argparse
import logging
import subprocess
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)-7s | %(message)s')
logger = logging.getLogger(__name__)


def build_parser():
    p = argparse.ArgumentParser(description='[N2] Build Ortholog Tree')
    p.add_argument('--config', '-c', help='YAML config')
    p.add_argument('--og', required=True, help='Orthogroup ID (e.g., OG0000001)')
    p.add_argument('--orthofinder-dir', help='OrthoFinder results directory')
    p.add_argument('--outdir', '-o', default='output/ortholog_tree')
    p.add_argument('--threads', '-t', type=int, default=0)
    p.add_argument('--model', default='MFP', help='IQ-TREE2 model selection (MFP, LG, WAG)')
    p.add_argument('--bootstrap', type=int, default=1000, help='Number of bootstrap replicates')
    return p


def main():
    parser = build_parser()
    args = parser.parse_args()

    cfg = {}
    if args.config:
        with open(args.config) as f:
            cfg = yaml.safe_load(f) or {}

    orthofinder_dir = Path(args.orthofinder_dir or cfg.get('orthofinder_dir', 'OrthoFinder'))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    threads = args.threads or os.cpu_count() or 4

    # Locate the alignment file for the orthogroup
    aln_dir = orthofinder_dir / 'MultipleSequenceAlignments'
    if not aln_dir.exists():
        results = sorted(orthofinder_dir.glob('Results_*'))
        if results:
            aln_dir = results[-1] / 'MultipleSequenceAlignments'

    fa_path = aln_dir / f'{args.og}.fa'
    if not fa_path.exists():
        logger.error(f'Alignment not found: {fa_path}')
        sys.exit(1)

    # Build tree
    prefix = outdir / f'{args.og}_tree'
    cmd = ['iqtree2', '-s', str(fa_path), '--prefix', str(prefix),
           '-T', str(threads), '-m', args.model, '-B', str(args.bootstrap)]
    logger.info(f'Running IQ-TREE2: {" ".join(cmd)}')
    subprocess.run(cmd, check=True)
    logger.info(f'Tree saved: {prefix}.treefile')


if __name__ == '__main__':
    main()
