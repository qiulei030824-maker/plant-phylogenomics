#!/usr/bin/env python3
"""[N4] Convert anchors + filter — Industrial-grade.

Converts LAST/LASTZ alignment anchors into filtered
synteny anchor files for MCscanX or similar tools.

Usage:
    python stepN4_convert_anchors.py [--config config.yaml]
"""

import os
import sys
import json
import yaml
import argparse
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)-7s | %(message)s')
logger = logging.getLogger(__name__)


def build_parser():
    p = argparse.ArgumentParser(description='[N4] Convert anchors + filter')
    p.add_argument('--config', '-c', help='YAML config')
    p.add_argument('--input-dir', help='Directory with anchor files from step N3')
    p.add_argument('--outdir', '-o', default='output/anchors')
    p.add_argument('--min-length', type=int, default=50, help='Minimum alignment length')
    p.add_argument('--min-identity', type=float, default=0.3, help='Minimum identity fraction')
    return p


def convert_anchors(input_dir: Path, outdir: Path, min_len: int = 50, min_idy: float = 0.3):
    """Convert tab-delimited anchors to filtered format."""
    outdir.mkdir(parents=True, exist_ok=True)

    anchor_files = list(input_dir.glob('*.anchors'))
    if not anchor_files:
        logger.error(f'No .anchors files found in {input_dir}')
        sys.exit(1)

    for af in anchor_files:
        pair_name = af.stem
        out_path = outdir / f'{pair_name}.filtered.anchors'
        total = 0
        kept = 0

        with open(af) as fin, open(out_path, 'w') as fout:
            for line in fin:
                parts = line.strip().split('\t')
                if len(parts) < 12:
                    continue
                total += 1
                # Parse MAF tab format
                qry_id, qry_len, qry_start, qry_end, qry_strand = parts[0], parts[1], parts[2], parts[3], parts[4]
                ref_id, ref_len, ref_start, ref_end, ref_strand = parts[5], parts[6], parts[7], parts[8], parts[9]
                score = float(parts[10])
                identity = float(parts[11]) if len(parts) > 11 else 0

                aln_len = min(int(qry_end) - int(qry_start), int(ref_end) - int(ref_start))
                if aln_len >= min_len and identity >= min_idy:
                    fout.write(f'{ref_id}\t{ref_start}\t{ref_end}\t{qry_id}\t{qry_start}\t{qry_end}\t{identity:.4f}\n')
                    kept += 1

        logger.info(f'{pair_name}: {kept}/{total} anchors kept')


def main():
    parser = build_parser()
    args = parser.parse_args()

    cfg = {}
    if args.config:
        with open(args.config) as f:
            cfg = yaml.safe_load(f) or {}

    input_dir = Path(args.input_dir or cfg.get('input_dir', 'output/last_alignment'))
    outdir = Path(args.outdir)

    convert_anchors(input_dir, outdir, args.min_length, args.min_identity)
    logger.info(f'Filtered anchors written to {outdir}')


if __name__ == '__main__':
    main()
