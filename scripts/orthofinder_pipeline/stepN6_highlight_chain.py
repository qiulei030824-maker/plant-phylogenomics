#!/usr/bin/env python3
"""[N6] Generate highlight chain plot data — Industrial-grade.

Prepares chain-level synteny visualization data from
anchor-filtered LAST alignments.

Usage:
    python stepN6_highlight_chain.py [--config config.yaml]
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
    p = argparse.ArgumentParser(description='[N6] Generate highlight chain data')
    p.add_argument('--config', '-c', help='YAML config')
    p.add_argument('--input-dir', help='Filtered anchors directory')
    p.add_argument('--outdir', '-o', default='output/chains')
    p.add_argument('--min-anchors', type=int, default=5, help='Minimum anchors per chain')
    p.add_argument('--max-gap', type=int, default=50000, help='Max gap between anchors to connect chain')
    return p


def build_chains(input_dir: Path, outdir: Path, min_anchors: int = 5, max_gap: int = 50000):
    """Build chain structures from anchors."""
    outdir.mkdir(parents=True, exist_ok=True)

    anchor_files = list(input_dir.glob('*.filtered.anchors'))
    if not anchor_files:
        logger.error(f'No .filtered.anchors files in {input_dir}')
        sys.exit(1)

    for af in anchor_files:
        pair_name = af.stem.replace('.filtered', '')
        # Group anchors by ref_qry pair
        pairs = defaultdict(list)
        with open(af) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 7:
                    continue
                ref_id, ref_s, ref_e = parts[0], int(parts[1]), int(parts[2])
                qry_id, qry_s, qry_e = parts[3], int(parts[4]), int(parts[5])
                identity = float(parts[6])
                key = (ref_id, qry_id)
                pairs[key].append({
                    'ref_start': ref_s, 'ref_end': ref_e,
                    'qry_start': qry_s, 'qry_end': qry_e,
                    'identity': identity
                })

        # Build chains
        chains = []
        for (ref_id, qry_id), anchors in pairs.items():
            anchors.sort(key=lambda x: (x['ref_start'], x['qry_start']))
            current_chain = [anchors[0]]
            for a in anchors[1:]:
                last = current_chain[-1]
                gap_ref = a['ref_start'] - last['ref_end']
                gap_qry = a['qry_start'] - last['qry_end']
                if gap_ref <= max_gap and gap_qry <= max_gap:
                    current_chain.append(a)
                else:
                    if len(current_chain) >= min_anchors:
                        chains.append((ref_id, qry_id, current_chain[:]))
                    current_chain = [a]
            if len(current_chain) >= min_anchors:
                chains.append((ref_id, qry_id, current_chain[:]))

        # Write output
        chain_path = outdir / f'{pair_name}_chains.json'
        chain_data = []
        for ref_id, qry_id, anchors in chains:
            chain_data.append({
                'ref_id': ref_id,
                'qry_id': qry_id,
                'n_anchors': len(anchors),
                'anchors': anchors
            })

        with open(chain_path, 'w') as f:
            json.dump(chain_data, f, indent=2)
        logger.info(f'{pair_name}: {len(chain_data)} chains written')


def main():
    parser = build_parser()
    args = parser.parse_args()

    cfg = {}
    if args.config:
        with open(args.config) as f:
            cfg = yaml.safe_load(f) or {}

    input_dir = Path(args.input_dir or cfg.get('input_dir', 'output/anchors'))
    outdir = Path(args.outdir)

    build_chains(input_dir, outdir, args.min_anchors, args.max_gap)
    logger.info(f'Chain data written to {outdir}')


if __name__ == '__main__':
    main()
