#!/usr/bin/env python3
"""
[N5] Generate BED files for jcvi synteny visualization.

Pipeline:
  1. Find existing BED files or convert GFF3 -> BED via jcvi
  2. Create seqids.txt for chromosome ordering

Usage:
    python stepN5_generate_bed.py [--config config.yaml]
"""
import os, sys, json, argparse, subprocess, logging
from pathlib import Path
from datetime import datetime


def setup_logging(log_dir=None, name="stepN5"):
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
        fh = logging.FileHandler(log_dir / f"stepN5_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger


def find_gff3(data_dir, species):
    sp_dir = Path(data_dir) / species
    if not sp_dir.is_dir():
        return None
    for f in sp_dir.iterdir():
        if f.name.endswith(".gff3") or f.name.endswith(".gff3.gz"):
            return f
    return None


def find_bed(data_dir, species):
    sp_dir = Path(data_dir) / species
    if not sp_dir.is_dir():
        return None
    for f in sp_dir.iterdir():
        if f.name.endswith(".bed"):
            return f
    return None


def run(species_list, data_dir, output_dir, logger=None):
    log = logger or setup_logging()
    log.info(f"{'='*60}")
    log.info(f"N5: Generate BED files for {len(species_list)} species")
    log.info(f"{'='*60}")
    
    output_dir = Path(output_dir) / "jcvi_chain"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    for sp in species_list:
        log.info(f"\n  Processing {sp}")
        
        # Check for existing BED
        existing_bed = find_bed(data_dir, sp)
        if existing_bed:
            log.info(f"    Found existing BED: {existing_bed.name}")
            os.system(f"cp {existing_bed} {output_dir / f'{sp}.bed'}")
            continue
        
        # Convert GFF3 -> BED via jcvi
        gff3_file = find_gff3(data_dir, sp)
        if gff3_file is None:
            log.warning(f"    No GFF3 found for {sp}")
            continue
        
        log.info(f"    Converting {gff3_file.name} -> BED...")
        os.chdir(str(output_dir))
        cmd = f"python -m jcvi.formats.gff bed --type=mRNA --key=ID {gff3_file} -o {sp}.bed"
        subprocess.run(cmd, shell=True, check=True, capture_output=True)
        log.info(f"    -> {output_dir / f'{sp}.bed'}")
    
    # Generate seqids
    log.info(f"\n  Generating seqids.txt...")
    seqids_file = output_dir / "seqids.txt"
    with open(seqids_file, "w") as out:
        for sp in species_list:
            bed_file = output_dir / f"{sp}.bed"
            if bed_file.exists():
                with open(bed_file) as f:
                    chrs = sorted(set(line.split()[0] for line in f if line.strip()))
                for chr_name in chrs:
                    out.write(f"{sp}:{chr_name}\n")
    log.info(f"    -> {seqids_file}")
    
    log.info(f"\n{'='*60}")
    log.info("N5 complete!")
    log.info(f"{'='*60}")
    return True


def main():
    parser = argparse.ArgumentParser(description="N5: Generate BED files")
    parser.add_argument("--config-json")
    parser.add_argument("--species", help="Comma-separated species")
    parser.add_argument("--data-dir", default="/data5/qiulei/pfam_pipeline/data")
    parser.add_argument("--output-dir")
    parser.add_argument("--log-dir")
    args = parser.parse_args()
    
    logger = setup_logging(args.log_dir)
    species = []
    if args.config_json:
        with open(args.config_json) as f:
            config = json.load(f)
        species = config.get("stepN5", {}).get("species_list", [])
        args.data_dir = config.get("stepN5", {}).get("data_dir", args.data_dir)
    if args.species:
        species = args.species.split(",")
    if not species:
        logger.error("No species provided")
        sys.exit(1)
    
    output_dir = args.output_dir or args.data_dir
    success = run(species, args.data_dir, output_dir, logger)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()