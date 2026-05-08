#!/usr/bin/env python3
"""[N3] Run LAST alignment — Industrial-grade.

Runs LASTZ or LAST alignments for synteny detection using
paired genomes derived from OrthoFinder results.

Usage:
    python stepN3_run_last_alignment.py [--config config.yaml]
"""

import os
import sys
import json
import yaml
import glob
import argparse
import logging
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)-7s | %(message)s')
logger = logging.getLogger(__name__)


def build_parser():
    p = argparse.ArgumentParser(description='[N3] Run LAST alignment')
    p.add_argument('--config', '-c', help='YAML config')
    p.add_argument('--genome-dir', help='Directory with genome FASTA files')
    p.add_argument('--outdir', '-o', default='output/last_alignment')
    p.add_argument('--threads', '-t', type=int, default=4)
    p.add_argument('--pairs', nargs='+', help='Genome pairs (ref query ref query ...)')
    return p


def run_last(ref: Path, query: Path, outdir: Path) -> bool:
    """Run LAST alignment for a pair of genomes."""
    ref_name = ref.stem
    query_name = query.stem
    pair_dir = outdir / f'{ref_name}_vs_{query_name}'
    pair_dir.mkdir(parents=True, exist_ok=True)

    # 1. Index reference
    if not (pair_dir / f'{ref_name}.prj').exists():
        logger.info(f'Indexing reference: {ref_name}')
        subprocess.run(['lastdb', '-P1', str(pair_dir / ref_name), str(ref)], check=True)

    # 2. Run lastal
    maf_path = pair_dir / f'{ref_name}_vs_{query_name}.maf'
    if not maf_path.exists():
        logger.info(f'Aligning {query_name} -> {ref_name}')
        with open(maf_path, 'w') as f:
            subprocess.run(['lastal', '-P1', str(pair_dir / ref_name), str(query),
                           '-f', 'MAF'], stdout=f, check=True)

    # 3. Convert to anchors
    anchors_path = pair_dir / f'{ref_name}_vs_{query_name}.anchors'
    if not anchors_path.exists():
        logger.info(f'Converting to anchors: {query_name}')
        subprocess.run(['maf-convert', 'tab', str(maf_path)],
                       stdout=open(anchors_path, 'w'), check=True)

    return True


def main():
    parser = build_parser()
    args = parser.parse_args()

    cfg = {}
    if args.config:
        with open(args.config) as f:
            cfg = yaml.safe_load(f) or {}

    genome_dir = Path(args.genome_dir or cfg.get('genome_dir', 'genomes'))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    threads = args.threads

    # Collect genome files
    genomes = sorted(genome_dir.glob('*.fa')) + sorted(genome_dir.glob('*.fasta'))
    if not genomes:
        logger.error(f'No genome FASTA files found in {genome_dir}')
        sys.exit(1)

    # Generate pairs
    if args.pairs:
        pairs = list(zip(args.pairs[::2], args.pairs[1::2]))
    else:
        # All-vs-all: first genome as reference, rest as queries
        ref = genomes[0]
        pairs = [(ref, q) for q in genomes[1:]]

    logger.info(f'Running {len(pairs)} genome pairs with {threads} threads')

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(run_last, ref, query, outdir): (ref.name, query.name)
                   for ref, query in pairs}
        for fut in as_completed(futures):
            ref_n, q_n = futures[fut]
            try:
                fut.result()
                logger.info(f'Done: {ref_n} vs {q_n}')
            except Exception as e:
                logger.error(f'Failed: {ref_n} vs {q_n}: {e}')

    logger.info('All LAST alignments complete')


if __name__ == '__main__':
    main()
