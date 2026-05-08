#!/usr/bin/env python3
"""
[Step 2] Analyze HMMER domtblout output from Step 1.

Usage:
    python scripts/02_parse_domain_table.py <PFAM_ID>
    python scripts/02_parse_domain_table.py <PFAM_ID> -i /path/to/domtblout -o /path/to/output
    python scripts/02_parse_domain_table.py PF00168
    python scripts/02_parse_domain_table.py PF00168 -i /data/xxx/hmmer/all_PF00168.domtblout

This script:
  1. Counts how many domains each gene (target name) has
  2. Counts how many unique genes each species has
  3. Outputs TSV files and a summary report
"""

import os
import sys
import argparse
import csv
from collections import Counter, defaultdict

# ── Path setup: allow import from config/ or local species_config ──────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)

try:
    from config.species_config import BASE_ROOT, resolve_pfam_paths
except ModuleNotFoundError:
    from species_config import BASE_ROOT, resolve_pfam_paths


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Analyze HMMER domtblout: gene domain counts + species gene counts",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s PF00168
  %(prog)s PF00168 -i /path/to/all_PF00168.domtblout -o /path/to/output
  %(prog)s PF00168 --delimiter "|"
        """,
    )
    parser.add_argument("pfam_id", type=str, help="Pfam accession (e.g. PF00168)")
    parser.add_argument("-i", "--input", type=str, default=None,
                        help="Path to domtblout file (default: auto-resolve from PFAM_ID)")
    parser.add_argument("-o", "--output-dir", type=str, default=None,
                        help="Output directory for stats (default: <domtblout_dir>/stats)")
    parser.add_argument("--delimiter", type=str, default="|",
                        help="Delimiter between species and gene ID (default: '|')")
    return parser.parse_args()


def parse_domtblout(filepath):
    """
    Parse HMMER domtblout file.
    
    Returns list of dicts with keys: gene, species, domain_idx, domain_count, 
    full_seq_evalue, dom_ievalue, ali_from, ali_to, env_from, env_to.
    """
    records = []
    skipped_comment = 0
    skipped_short = 0
    skipped_gene = 0

    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                skipped_comment += 1
                continue
            cols = line.strip().split()
            if len(cols) < 19:
                skipped_short += 1
                continue

            gene = cols[0]
            if not gene or gene == '-':
                skipped_gene += 1
                continue

            record = {
                "gene": gene,
                "seq_evalue": float(cols[6]) if cols[6] not in ('-', '') else None,
                "score": float(cols[7]) if cols[7] not in ('-', '') else None,
                "bias": float(cols[8]) if cols[8] not in ('-', '') else None,
                "domain_idx": int(cols[9]),
                "domain_count": int(cols[10]),
                "dom_c_evalue": float(cols[11]) if cols[11] not in ('-', '') else None,
                "dom_i_evalue": float(cols[12]) if cols[12] not in ('-', '') else None,
                "dom_score": float(cols[13]) if cols[13] not in ('-', '') else None,
                "ali_from": int(cols[17]),
                "ali_to": int(cols[18]),
            }
            records.append(record)

    info = {
        "total_lines": len(records),
        "skipped_comment": skipped_comment,
        "skipped_short": skipped_short,
        "skipped_gene": skipped_gene,
    }
    return records, info


def resolve_species(gene, delimiter):
    """Extract species name from gene ID using the delimiter."""
    if delimiter and delimiter in gene:
        return gene.split(delimiter, 1)[0]
    return "UNKNOWN"


def task1_gene_domain_count(records, output_path, delimiter=None):
    """
    Task 1: Count how many domains each gene (target name) has.
    Output TSV: Gene_ID | Domain_Count | Domain_Indices | Species
    """
    gene_domains = defaultdict(list)
    for r in records:
        gene_domains[r["gene"]].append(r["domain_idx"])

    species = {}
    for gene in gene_domains:
        species[gene] = resolve_species(gene, delimiter)

    rows = []
    for gene, domains in gene_domains.items():
        rows.append({
            "Gene_ID": gene,
            "Domain_Count": len(domains),
            "Domain_Indices": ",".join(str(d) for d in sorted(domains)),
            "Species": species[gene],
        })
    rows.sort(key=lambda x: -x["Domain_Count"])

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["Gene_ID", "Domain_Count", "Domain_Indices", "Species"], delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    return rows


def task2_species_gene_count(records, output_path, delimiter=None):
    """
    Task 2: Count how many unique genes each species has.
    Output TSV: Species | Gene_Count
    """
    species_genes = defaultdict(set)
    for r in records:
        sp = resolve_species(r["gene"], delimiter)
        species_genes[sp].add(r["gene"])

    rows = [{"Species": sp, "Gene_Count": len(genes)}
            for sp, genes in species_genes.items()]
    rows.sort(key=lambda x: -x["Gene_Count"])

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["Species", "Gene_Count"], delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main():
    args = parse_args()
    pfam_id = args.pfam_id.upper()

    # ── Resolve input domtblout path ────────────────────────────────────────
    if args.input is not None:
        domtblout = os.path.abspath(args.input)
    else:
        paths = resolve_pfam_paths(pfam_id)
        domtblout = str(paths["hmmer_file"])

    if not os.path.exists(domtblout):
        print(f"ERROR: domtblout not found: {domtblout}")
        sys.exit(1)

    # ── Resolve output dir ──────────────────────────────────────────────────
    if args.output_dir is not None:
        outdir = os.path.abspath(args.output_dir)
    else:
        outdir = os.path.join(os.path.dirname(domtblout), "stats")
    os.makedirs(outdir, exist_ok=True)

    print("=" * 60)
    print(f"Step 2: Analyze HMMER domtblout for {pfam_id}")
    print(f"  Input:   {domtblout}")
    print(f"  Output:  {outdir}/")
    print(f"  Delimiter: '{args.delimiter}'")
    print("=" * 60)

    # ── Parse domtblout ─────────────────────────────────────────────────────
    print(f"\n[Step 1] Parsing domtblout...")
    records, info = parse_domtblout(domtblout)
    print(f"  Total domain hits:     {info['total_lines']}")
    print(f"  Comments skipped:      {info['skipped_comment']}")
    print(f"  Short lines skipped:   {info['skipped_short']}")
    print(f"  Empty genes skipped:   {info['skipped_gene']}")

    # ── Task 1: Gene domain count ───────────────────────────────────────────
    print(f"\n[Task 1] Counting domains per gene...")
    task1_out = os.path.join(outdir, f"{pfam_id}_gene_domain_count.tsv")
    gene_rows = task1_gene_domain_count(records, task1_out, args.delimiter)
    total_genes = len(gene_rows)
    print(f"  Total genes with domain hits: {total_genes}")
    print(f"  Output: {task1_out}")
    print(f"  Top 5 multi-domain genes:")
    multi = [r for r in gene_rows if r["Domain_Count"] > 1]
    for r in gene_rows[:5]:
        print(f"    {r['Gene_ID']}: {r['Domain_Count']} domains")

    # ── Task 2: Species gene count ──────────────────────────────────────────
    print(f"\n[Task 2] Counting unique genes per species...")
    task2_out = os.path.join(outdir, f"{pfam_id}_species_gene_count.tsv")
    species_rows = task2_species_gene_count(records, species_rows_path := task2_out, args.delimiter)
    print(f"  Total species: {len(species_rows)}")
    print(f"  Output: {task2_out}")
    print(f"  Top 10 species by gene count:")
    for r in species_rows[:10]:
        print(f"    {r['Species']}: {r['Gene_Count']} genes")

    # ── Summary ─────────────────────────────────────────────────────────────
    summary_path = os.path.join(outdir, f"{pfam_id}_summary.txt")
    with open(summary_path, 'w') as f:
        f.write(f"# Summary: {pfam_id} domtblout analysis\n")
        f.write(f"# Input: {domtblout}\n")
        f.write(f"# Date: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"Total domain hits:     {info['total_lines']}\n")
        f.write(f"Total genes:           {total_genes}\n")
        f.write(f"Multi-domain genes:    {len(multi)}\n")
        f.write(f"Single-domain genes:   {total_genes - len(multi)}\n")
        f.write(f"Total species:         {len(species_rows)}\n\n")
        f.write("Files:\n")
        f.write(f"  gene_domain_count:    {task1_out}\n")
        f.write(f"  species_gene_count:   {task2_out}\n")
        f.write(f"  this summary:         {summary_path}\n")

    print(f"\n  Summary: {summary_path}")

    # ── Print report ───────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("REPORT")
    print("=" * 60)
    print(f"  Input:            {domtblout}")
    print(f"  Total hits:       {info['total_lines']}")
    print(f"  Total genes:      {total_genes}")
    print(f"  Multi-domain:     {len(multi)} ({100*len(multi)//max(total_genes,1)}%)")
    print(f"  Single-domain:    {total_genes - len(multi)} ({100*(total_genes-len(multi))//max(total_genes,1)}%)")
    print(f"  Total species:    {len(species_rows)}")
    print(f"\n  Gene domain stats:    {task1_out}")
    print(f"  Species gene stats:  {task2_out}")
    print(f"  Summary:             {summary_path}")
    print("=" * 60)


if __name__ == "__main__":
    main()
