#!/usr/bin/env python3
"""[N5] Generate BED files from anchor-filtered synteny — Industrial-grade.

Converts filtered synteny anchors to BED format for
visualization with pyCirgo or similar tools.

Usage:
    python stepN5_generate_bed.py [--config config.yaml]
"""

import os
import sys
import json
import yaml
import argparse
import logging
from pathlib import Path
from collections import defaultdict

logging.basicConfig(level=logging.INFO, format='%(asctime)s | %(levelname)-7s | %(message)s')
logger = logging.getLogger(__name__)


def build_parser():
    p = argparse.ArgumentParser(description='[N5] Generate BED files')
    p.add_argument('--config', '-c', help='YAML config')
    p.add_argument('--input-dir', help='Filtered anchors directory')
    p.add_argument('--outdir', '-o', default='output/bed')
    p.add_argument('--flank', type=int, default=1000, help='Flank size around anchors')
    return p


def anchors_to_bed(input_dir: Path, outdir: Path, flank: int = 1000):
    """Convert filtered anchors to BED format."""
    outdir.mkdir(parents=True, exist_ok=True)

    anchor_files = list(input_dir.glob('*.filtered.anchors'))
    if not anchor_files:
        logger.error(f'No .filtered.anchors files in {input_dir}')
        sys.exit(1)

    for af in anchor_files:
        pair_name = af.stem.replace('.filtered', '')
        ref_bed = outdir / f'{pair_name}_ref.bed'
        qry_bed = outdir / f'{pair_name}_qry.bed'

        ref_records = defaultdict(list)
        qry_records = defaultdict(list)

        with open(af) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 7:
                    continue
                ref_id, ref_start, ref_end = parts[0], int(parts[1]), int(parts[2])
                qry_id, qry_start, qry_end = parts[3], int(parts[4]), int(parts[5])
                identity = parts[6]

                ref_records[ref_id].append((ref_start, ref_end, qry_id, qry_start, qry_end, identity))
                qry_records[qry_id].append((qry_start, qry_end, ref_id, ref_start, ref_end, identity))

        def write_bed(path, records, is_ref=True):
            with open(path, 'w') as out:
                for sid in sorted(records.keys()):
                    blocks = sorted(records[sid])
                    for b in blocks:
                        start = max(0, b[0] - flank)
                        end = b[1] + flank
                        name = f'{sid}|{b[2]}|{b[5]}'
                        out.write(f'{sid}\t{start}\t{end}\t{name}\t0\t+\n')

        write_bed(ref_bed, ref_records)
        write_bed(qry_bed, qry_records)
        logger.info(f'Wrote {pair_name} BED files')


def main():
    parser = build_parser()
    args = parser.parse_args()

    cfg = {}
    if args.config:
        with open(args.config) as f:
            cfg = yaml.safe_load(f) or {}

    input_dir = Path(args.input_dir or cfg.get('input_dir', 'output/anchors'))
    outdir = Path(args.outdir)

    anchors_to_bed(input_dir, outdir, args.flank)
    logger.info(f'BED files written to {outdir}')


if __name__ == '__main__':
    main()
