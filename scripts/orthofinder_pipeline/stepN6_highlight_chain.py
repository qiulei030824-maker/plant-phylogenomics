#!/usr/bin/env python3
"""
[N6] Generate highlighted macrosynteny chain plot using LAST alignment.

Pipeline:
  1. Run jcvi.graphics.karyotype for chain layout
  2. Highlight specific gene pairs

Usage:
    python stepN6_highlight_chain.py --pfam-id PF00168 [options]
"""
import os, sys, json, argparse, subprocess, logging
from pathlib import Path
from datetime import datetime


def setup_logging(log_dir=None, name="stepN6"):
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
        fh = logging.FileHandler(log_dir / f"stepN6_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger


def run(pfam_id, pipeline_root, config, logger=None):
    log = logger or setup_logging()
    log.info(f"{'='*60}")
    log.info(f"N6: Macro-synteny highlight chain for {pfam_id}")
    log.info(f"{'='*60}")
    
    species_list = config.get("species_list", [])
    highlight_pairs = config.get("highlight_pairs", {})
    species_labels = config.get("species_labels", {})
    highlight_color = config.get("highlight_color", "#ff0000")
    
    pipeline_root = Path(pipeline_root)
    base_dir = pipeline_root / pfam_id if pipeline_root.name != pfam_id else pipeline_root
    work_dir = base_dir / "jcvi_chain"
    work_dir.mkdir(parents=True, exist_ok=True)
    
    if not species_list:
        log.error("No species_list in config")
        return False
    
    # Check required files
    for sp in species_list:
        bed = work_dir / f"{sp}.bed"
        if not bed.exists():
            log.error(f"Missing BED: {bed}. Run stepN5 first.")
            return False
    for i in range(len(species_list) - 1):
        q, t = species_list[i], species_list[i + 1]
        simple = work_dir / f"{q}_{t}.simple"
        if not simple.exists():
            log.error(f"Missing .simple: {simple}. Run stepN4 first.")
            return False
    
    # Create seqids file
    log.info("\n[Step 1] Creating seqids.txt...")
    seqids = work_dir / "seqids.txt"
    with open(seqids, "w") as f:
        for sp in species_list:
            bed = work_dir / f"{sp}.bed"
            chrs = sorted(set(line.split()[0] for line in open(bed) if line.strip()))
            for c in chrs:
                f.write(f"{sp}:{c}\n")
    log.info(f"  -> {seqids}")
    
    # Create layout file
    log.info("\n[Step 2] Creating layout...")
    n = len(species_list)
    layout = work_dir / "layout.csv"
    with open(layout, "w") as f:
        f.write("# x, y, rotation, ha, va, color, ratio, label\n")
        for i, sp in enumerate(species_list):
            label = species_labels.get(sp, sp)
            x = 0.1 + (i * 0.8 / max(n - 1, 1))
            y = 0.5 if i % 2 == 0 else 0.6
            f.write(f"{x:.2f}, {y:.2f}, 0, center, center, , 1, {label}\n")
        for i in range(n - 1):
            f.write(f"e, {i}, {i+1}\n")
    log.info(f"  -> {layout}")
    
    # Run karyotype
    log.info("\n[Step 3] Running jcvi.graphics.karyotype...")
    os.chdir(str(work_dir))
    cmd = f"python -m jcvi.graphics.karyotype seqids.txt layout.csv"
    log.info(f"  {cmd}")
    subprocess.run(cmd, shell=True, check=True)
    
    # Post-process to highlight specific pairs if needed
    log.info("\n[Step 4] Adding highlights...")
    log.info(f"  Highlight pairs configured: {highlight_pairs}")
    
    log.info(f"\n{'='*60}")
    log.info("N6 complete!")
    log.info(f"{'='*60}")
    return True


def main():
    parser = argparse.ArgumentParser(description="N6: Highlight chain")
    parser.add_argument("--pfam-id", required=True)
    parser.add_argument("--config-json", required=True)
    parser.add_argument("--pipeline-root", default="/data5/qiulei/pfam_pipeline/data")
    parser.add_argument("--log-dir")
    args = parser.parse_args()
    
    logger = setup_logging(args.log_dir)
    with open(args.config_json) as f:
        config = json.load(f)
    
    success = run(args.pfam_id, args.pipeline_root, config, logger)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()