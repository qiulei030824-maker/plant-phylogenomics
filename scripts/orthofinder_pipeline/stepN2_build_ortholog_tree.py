#!/usr/bin/env python3
"""
[N2] Build ortholog tree for a target gene using OrthoFinder results.

Pipeline:
  1. Find orthogroup for target gene from Orthogroups.tsv
  2. Extract sequences from species proteomes
  3. MAFFT alignment
  4. trimAl trimming
  5. FastTree quick preview
  6. RAxML-NG final tree with bootstrapping

Usage:
    python stepN2_build_ortholog_tree.py --target-gene AT1G09070 [options]
"""
import os, sys, re, json, shutil, argparse, subprocess, logging, textwrap
from pathlib import Path
from datetime import datetime
from collections import defaultdict


def setup_logging(log_dir=None, name="stepN2", level=logging.INFO):
    logger = logging.getLogger(name)
    logger.setLevel(level)
    if logger.handlers:
        return logger
    fmt = logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s")
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(fmt)
    logger.addHandler(ch)
    if log_dir:
        log_dir = Path(log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(log_dir / f"{name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger


def find_orthogroup(target_gene, orthogroups_tsv):
    """Find which orthogroup contains the target gene."""
    with open(orthogroups_tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            og_name = parts[0]
            for cell in parts[1:]:
                genes = [g.strip() for g in cell.split(", ") if g.strip()]
                if target_gene in genes:
                    return og_name, dict(zip(orthogroups_tsv.parent.parent.joinpath("proteomes").iterdir() if False else [], []))
    return None, None


def run(target_gene, orthofinder_dir, pfam_id, output_dir, threads, bootstrap, logger=None):
    log = logger or setup_logging()
    log.info(f"{'='*60}")
    log.info(f"N2: Building ortholog tree for {target_gene}")
    log.info(f"{'='*60}")
    
    orthofinder_dir = Path(orthofinder_dir)
    results_dir = orthofinder_dir / "Results_" + max(
        (d for d in orthofinder_dir.iterdir() if d.name.startswith("Results_")),
        key=lambda d: d.stat().st_mtime, default=None
    ).name if any(d.name.startswith("Results_") for d in orthofinder_dir.iterdir()) else None
    
    if results_dir is None or not results_dir.exists():
        log.error(f"No OrthoFinder Results_* directory found in {orthofinder_dir}")
        return False
    
    orthogroups_tsv = results_dir / "Orthogroups.tsv"
    if not orthogroups_tsv.exists():
        log.error(f"Orthogroups.tsv not found: {orthogroups_tsv}")
        return False
    
    # Find orthogroup
    og_name = None
    species_genes = {}
    with open(orthogroups_tsv) as f:
        header = f.readline().strip().split("\t")
        species_list = header[1:]
        for line in f:
            parts = line.strip().split("\t")
            name = parts[0]
            for i, cell in enumerate(parts[1:]):
                if i < len(species_list):
                    genes = [g.strip() for g in cell.split(", ") if g.strip()]
                    if target_gene in genes:
                        og_name = name
                    if og_name == name and genes:
                        species_genes[species_list[i]] = genes
    
    if og_name is None:
        log.error(f"Target gene {target_gene} not found in any orthogroup")
        return False
    log.info(f"  Orthogroup: {og_name}")
    
    pfam_dir = Path(pfam_id) if not Path(pfam_id).exists() else Path(pfam_id)
    output_dir = Path(output_dir) if output_dir else pfam_dir / "ortholog_tree"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract sequences
    proteome_dir = orthofinder_dir / "proteomes"
    seqs_fa = output_dir / f"{og_name}_raw.fa"
    
    log.info(f"\n[Step 1] Extracting sequences...")
    from Bio import SeqIO
    count = 0
    with open(seqs_fa, "w") as out:
        for sp, genes in species_genes.items():
            proteome_file = proteome_dir / f"{sp}.fa"
            if not proteome_file.exists():
                proteome_file = proteome_dir / f"{sp}.fasta"
            if not proteome_file.exists():
                continue
            for record in SeqIO.parse(proteome_file, "fasta"):
                rec_id = record.id.split()[0]
                if rec_id in genes or record.id in genes:
                    record.id = f"{sp}|{rec_id}"
                    record.description = ""
                    SeqIO.write(record, out, "fasta")
                    count += 1
    log.info(f"  {count} sequences -> {seqs_fa}")
    
    # MAFFT alignment
    aligned_fa = output_dir / f"{og_name}_aligned.fa"
    log.info(f"\n[Step 2] MAFFT alignment...")
    cmd = f"mafft --auto --thread {threads} {seqs_fa} > {aligned_fa} 2>/dev/null"
    subprocess.run(cmd, shell=True, check=True)
    log.info(f"  -> {aligned_fa}")
    
    # trimAl
    trimmed_fa = output_dir / f"{og_name}_trimmed.fa"
    log.info(f"\n[Step 3] trimAl trimming...")
    cmd = f"trimal -in {aligned_fa} -out {trimmed_fa} -automated1"
    subprocess.run(cmd, shell=True, check=True)
    log.info(f"  -> {trimmed_fa}")
    
    # FastTree preview
    quick_tree = output_dir / f"{og_name}_quick.tre"
    log.info(f"\n[Step 4] FastTree quick tree...")
    cmd = f"FastTree -lg {trimmed_fa} > {quick_tree} 2>/dev/null"
    subprocess.run(cmd, shell=True, check=True)
    
    # RAxML-NG
    log.info(f"\n[Step 5] RAxML-NG with {bootstrap} bootstraps...")
    raxml_dir = output_dir / "raxml"
    raxml_dir.mkdir(exist_ok=True)
    cmd = (f"raxml-ng --msa {trimmed_fa} --model LG+G4 --prefix {raxml_dir / og_name} "
           f"--threads {threads} --seed 42 --bsconverge {bootstrap} --bs-tree {bootstrap}")
    subprocess.run(cmd, shell=True, check=True)
    
    log.info(f"\n{'='*60}")
    log.info(f"N2 complete! Results in {output_dir}")
    log.info(f"{'='*60}")
    return True


def main():
    parser = argparse.ArgumentParser(description="N2: Build ortholog tree")
    parser.add_argument("--target-gene", required=True)
    parser.add_argument("--orthofinder-dir", default="/data5/qiulei/pfam_pipeline/data/orthofinder")
    parser.add_argument("--pfam-id", default="PF00168")
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--bootstrap", type=int, default=100)
    args = parser.parse_args()
    
    success = run(args.target_gene, args.orthofinder_dir, args.pfam_id,
                  args.output_dir, args.threads, args.bootstrap)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()