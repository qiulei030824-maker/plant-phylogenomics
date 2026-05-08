#!/usr/bin/env python3
"""[N7] Create microsynteny visualization — Industrial-grade.

Generates publication-ready microsynteny plots showing
gene order conservation across species using anchor/chain data.

Usage:
    python stepN7_microsynteny.py [--config config.yaml]
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

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.collections import LineCollection
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    logger.warning('matplotlib not available; will output JSON only')


def build_parser():
    p = argparse.ArgumentParser(description='[N7] Create microsynteny visualization')
    p.add_argument('--config', '-c', help='YAML config')
    p.add_argument('--chain-dir', help='Directory with chain JSON files')
    p.add_argument('--outdir', '-o', default='output/microsynteny')
    p.add_argument('--species', nargs='+', help='Species order (ref qry ...)')
    p.add_argument('--region', help='Region to visualize (e.g., Chr1:1000000-2000000)')
    p.add_argument('--format', choices=['pdf', 'svg', 'png'], default='pdf')
    return p


def plot_microsynteny(chain_data: dict, species_order: list, out_path: Path, region: str = None):
    """Generate microsynteny plot."""
    if not HAS_MPL:
        # Save JSON only
        with open(out_path.with_suffix('.json'), 'w') as f:
            json.dump(chain_data, f, indent=2)
        return

    fig, axes = plt.subplots(len(species_order) - 1, 1,
                              figsize=(12, 2 * (len(species_order) - 1)),
                              sharex=True)
    if len(species_order) == 2:
        axes = [axes]

    colors = plt.cm.tab20.colors

    for idx, (ref_name) in enumerate(species_order[:-1]):
        ax = axes[idx]
        qry_name = species_order[idx + 1]
        pair_key = f'{ref_name}_vs_{qry_name}'
        pair_chains = chain_data.get(pair_key, [])

        # Draw chromosomes
        ax.axhline(0, color='gray', linewidth=2, alpha=0.5)
        ax.axhline(-1, color='gray', linewidth=2, alpha=0.5)

        # Draw chains
        for ci, chain in enumerate(pair_chains[:50]):  # Limit to 50 chains
            anchors = chain.get('anchors', [])
            color = colors[ci % len(colors)]
            for a in anchors:
                # Ref segment
                ax.plot([a['ref_start'], a['ref_end']], [0, 0],
                        color=color, linewidth=2, alpha=0.7)
                # Query segment (offset at -1)
                ax.plot([a['qry_start'], a['qry_end']], [-1, -1],
                        color=color, linewidth=2, alpha=0.7)
                # Connection line
                ax.plot([a['ref_start'] + (a['ref_end'] - a['ref_start']) / 2,
                         a['qry_start'] + (a['qry_end'] - a['qry_start']) / 2],
                        [0, -1], color=color, linewidth=0.5, alpha=0.3)

        # Labels
        ax.set_ylabel(f'{ref_name}\nvs\n{qry_name}', fontsize=8)
        ax.set_yticks([0, -1])
        ax.set_yticklabels([ref_name, qry_name], fontsize=7)
        ax.set_xlim(auto=True)

    ax.set_xlabel('Genomic Position (bp)')
    plt.tight_layout()
    plt.savefig(str(out_path), dpi=300, bbox_inches='tight')
    plt.close()
    logger.info(f'Microsynteny plot: {out_path}')


def main():
    parser = build_parser()
    args = parser.parse_args()

    cfg = {}
    if args.config:
        with open(args.config) as f:
            cfg = yaml.safe_load(f) or {}

    chain_dir = Path(args.chain_dir or cfg.get('chain_dir', 'output/chains'))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load chain data
    chain_files = sorted(chain_dir.glob('*_chains.json'))
    if not chain_files:
        logger.error(f'No chain files found in {chain_dir}')
        sys.exit(1)

    all_chain_data = {}
    for cf in chain_files:
        pair_name = cf.stem.replace('_chains', '')
        with open(cf) as f:
            all_chain_data[pair_name] = json.load(f)

    # Determine species order
    if args.species:
        species_order = args.species
    else:
        # Extract from chain file names
        species_order = list(set(
            pair.split('_vs_')[0] for pair in all_chain_data
        ) | set(
            pair.split('_vs_')[1] for pair in all_chain_data
        ))
        species_order.sort()

    # Plot
    out_path = outdir / f'microsynteny.{args.format}'
    plot_microsynteny(all_chain_data, species_order, out_path, args.region)
    logger.info('Microsynteny visualization complete')


if __name__ == '__main__':
    main()
