#!/usr/bin/env python3
"""
[N3] Run LAST alignment between consecutive species pairs.

Pipeline:
  1. Build LAST DB for each species
  2. Run lastal for each consecutive pair
  3. Filter to keep best hit per query

Usage:
    python stepN3_run_last_alignment.py [--config config.yaml]
"""
import os, sys, json, argparse, subprocess, logging
from pathlib import Path
from datetime import datetime


def setup_logging(log_dir=None, name="stepN3"):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    if logger.handlers:
        return logger
    fmt = logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s")
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(ch)
    if log_dir:
        log_dir = Path(log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_dir / f"stepN3_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger


def run(species_list, threads, cds_dir, output_dir, logger=None):
    log = logger or setup_logging()
    log.info(f"{'='*60}")
    log.info(f"N3: LAST alignment for {len(species_list)} species")
    log.info(f"{'='*60}")
    
    output_dir = Path(output_dir) / "last"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Build LAST databases
    log.info("\n[Step 1] Building LAST databases...")
    for sp in species_list:
        cds_fa = Path(cds_dir) / f"{sp}.cds.fa"
        if not cds_fa.exists():
            log.warning(f"  CDS not found: {cds_fa}, skipping {sp}")
            continue
        db_prefix = output_dir / f"{sp}"
        log.info(f"  Building DB for {sp}...")
        subprocess.run(f"lastdb -P{threads} {db_prefix} {cds_fa}", shell=True, check=True,
                      stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    # Run pairwise alignments
    log.info("\n[Step 2] Running pairwise alignments...")
    for i in range(len(species_list) - 1):
        q = species_list[i]
        t = species_list[i + 1]
        db_prefix = output_dir / t
        cds_fa = Path(cds_dir) / f"{q}.cds.fa"
        if not db_prefix.with_suffix(".prj").exists() or not cds_fa.exists():
            log.warning(f"  Skipping {q} -> {t}")
            continue
        raw_maf = output_dir / f"{q}_vs_{t}.maf"
        best_maf = output_dir / f"{q}_vs_{t}.best.maf"
        log.info(f"  Aligning {q} -> {t}...")
        cmd = f"lastal -P{threads} -m50 -E0.05 {db_prefix} {cds_fa} > {raw_maf} 2>/dev/null"
        subprocess.run(cmd, shell=True, check=True)
        # Filter best hit
        cmd = f"last-split -m1 {raw_maf} | last-maplet | maf-convert html > {best_maf} 2>/dev/null"
        subprocess.run(cmd, shell=True, check=True)
        log.info(f"    -> {best_maf}")
    
    log.info(f"\n{'='*60}")
    log.info("N3 complete!")
    log.info(f"{'='*60}")
    return True


def main():
    parser = argparse.ArgumentParser(description="N3: LAST alignment")
    parser.add_argument("--config-json", help="JSON config")
    parser.add_argument("--species", help="Comma-separated species list")
    parser.add_argument("--cds-dir", default="/data5/qiulei/pfam_pipeline/data/PF00168/cds")
    parser.add_argument("--output-dir", default="/data5/qiulei/pfam_pipeline/data/PF00168")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--log-dir")
    args = parser.parse_args()
    
    logger = setup_logging(args.log_dir)
    
    if args.config_json:
        with open(args.config_json) as f:
            config = json.load(f)
        species = config.get("stepN3", {}).get("species_list", args.species.split(",") if args.species else [])
        cds_dir = config.get("stepN3", {}).get("cds_dir", args.cds_dir)
        threads = config.get("stepN3", {}).get("threads", args.threads)
        output_dir = config.get("stepN3", {}).get("output_dir", args.output_dir)
    elif args.species:
        species = args.species.split(",")
        cds_dir = args.cds_dir
        threads = args.threads
        output_dir = args.output_dir
    else:
        logger.error("Provide --config-json or --species")
        sys.exit(1)
    
    success = run(species, threads, cds_dir, output_dir, logger)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()