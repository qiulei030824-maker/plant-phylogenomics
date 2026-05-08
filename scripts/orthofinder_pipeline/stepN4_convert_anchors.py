#!/usr/bin/env python3
"""
[N4] Convert anchors + LAST alignment to jcvi SimpleFile format (.simple).

Pipeline:
  1. Read anchor file (jcvi pairwise alignment anchors)
  2. Read LAST alignment
  3. Merge into .simple format

Usage:
    python stepN4_convert_anchors.py [--config config.yaml]
"""
import os, sys, json, argparse, logging
from pathlib import Path
from datetime import datetime


def setup_logging(log_dir=None, name="stepN4"):
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
        fh = logging.FileHandler(log_dir / f"stepN4_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger


def load_last_maf(maf_file):
    """Load LAST MAF and return dict of alignments."""
    alignments = {}
    with open(maf_file) as f:
        current = {}
        for line in f:
            if line.startswith("a"):
                if current:
                    pass
                current = {"score": 0, "pairs": []}
                for part in line.strip().split():
                    if part.startswith("score="):
                        current["score"] = float(part.split("=")[1])
            elif line.startswith("s"):
                parts = line.strip().split()
                if len(parts) >= 7:
                    seqid = parts[1]
                    start = int(parts[2])
                    alen = int(parts[3])
                    strand = parts[4]
                    src_size = int(parts[5])
                    current.setdefault("pairs", []).append({"seqid": seqid, "start": start, "alen": alen, "strand": strand, "src_size": src_size})
            elif line.strip() == "":
                if current and len(current.get("pairs", [])) == 2:
                    q = current["pairs"][0]
                    s = current["pairs"][1]
                    alignments[(q["seqid"], s["seqid"])] = {"score": current["score"], "strand": s["strand"]}
                current = {}
    return alignments


def run(species_list, anchor_dir, last_dir, output_dir, logger=None):
    log = logger or setup_logging()
    log.info(f"{'='*60}")
    log.info(f"N4: Convert anchors to SimpleFile format")
    log.info(f"{'='*60}")
    
    output_dir = Path(output_dir) / "jcvi_chain"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    for i in range(len(species_list) - 1):
        q = species_list[i]
        t = species_list[i + 1]
        log.info(f"\n  Processing {q} -> {t}")
        
        # Read anchors
        anchor_file = Path(anchor_dir) / f"{q}_{t}.anchors"
        if not anchor_file.exists():
            log.warning(f"    Anchors not found: {anchor_file}")
            continue
        anchors = []
        with open(anchor_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    anchors.append((parts[0], parts[1]))
        log.info(f"    {len(anchors)} anchors loaded")
        
        # Load LAST alignment
        maf_best = Path(last_dir) / f"{q}_vs_{t}.best.maf"
        last_align = {}
        if maf_best.exists():
            last_align = load_last_maf(maf_best)
            log.info(f"    {len(last_align)} LAST alignments loaded")
        
        # Write .simple output
        simple_file = output_dir / f"{q}_{t}.simple"
        written = 0
        with open(simple_file, "w") as f:
            for ga, gb in anchors:
                key = (ga, gb)
                aln = last_align.get(key, {})
                score = aln.get("score", 0)
                strand = aln.get("strand", "+")
                f.write(f"{ga}\t{ga}\t0\t0\t{score}\t{strand}\t{gb}\t{gb}\t0\t0\n")
                written += 1
                # Also write reversed
                f.write(f"{gb}\t{gb}\t0\t0\t{score}\t-\t{ga}\t{ga}\t0\t0\n")
                written += 1
        log.info(f"    {written} lines -> {simple_file}")
    
    log.info(f"\n{'='*60}")
    log.info("N4 complete!")
    log.info(f"{'='*60}")
    return True


def main():
    parser = argparse.ArgumentParser(description="N4: Convert anchors")
    parser.add_argument("--config-json")
    parser.add_argument("--species", help="Comma-separated species")
    parser.add_argument("--anchor-dir")
    parser.add_argument("--last-dir")
    parser.add_argument("--output-dir")
    parser.add_argument("--log-dir")
    args = parser.parse_args()
    
    logger = setup_logging(args.log_dir)
    species = args.species.split(",") if args.species else []
    
    if args.config_json:
        with open(args.config_json) as f:
            config = json.load(f)
        species = config.get("stepN4", {}).get("species_list", species)
        args.anchor_dir = config.get("stepN4", {}).get("anchor_dir", args.anchor_dir)
    
    last_dir = args.last_dir or (Path(args.anchor_dir).parent / "last" if args.anchor_dir else None)
    output_dir = args.output_dir or (Path(last_dir).parent if last_dir else Path.cwd())
    
    if not species:
        logger.error("No species provided")
        sys.exit(1)
    
    success = run(species, args.anchor_dir, str(last_dir), str(output_dir), logger)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()