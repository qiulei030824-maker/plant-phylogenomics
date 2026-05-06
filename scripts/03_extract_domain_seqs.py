#!/usr/bin/env python3
"""
[Step 3] Extract domain-containing CDS and PEP sequences per species
from the HMMER domtblout result.

Usage:
    python 03_extract_domain_seqs.py <PFAM_ID>
    python 03_extract_domain_seqs.py <PFAM_ID> -i /path/to/domtblout -o /path/to/output_dir
    python 03_extract_domain_seqs.py PF00168
    python 03_extract_domain_seqs.py PF00168 -i /data/xxx/hmmer/all_PF00168.domtblout --skip-cds
    python 03_extract_domain_seqs.py PF00168 --data-dir /custom/data/dir

Output:
    {output_dir}/pep/{species}.pep.fa
    {output_dir}/cds/{species}.cds.fa
"""

import os
import re
import sys
import argparse
import subprocess
from collections import defaultdict

# ── Path setup: allow import from config/ or local species_config ──────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)

try:
    from config.species_config import DATA_DIR as _DEF_DATA_DIR, SKIP_DIRS as _DEF_SKIP_DIRS, resolve_pfam_paths
except ModuleNotFoundError:
    from species_config import DATA_DIR as _DEF_DATA_DIR, SKIP_DIRS as _DEF_SKIP_DIRS, resolve_pfam_paths


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract domain-containing PEP/CDS sequences per species from HMMER domtblout",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s PF00168
  %(prog)s PF00168 -i /path/to/all_PF00168.domtblout -o /path/to/output_dir
  %(prog)s PF00168 --data-dir /custom/data --skip-cds
  %(prog)s PF00168 --delimiter "@" --force
        """,
    )
    parser.add_argument("pfam_id", type=str, help="Pfam accession (e.g. PF00168)")
    parser.add_argument("-i", "--input", type=str, default=None,
                        help="Path to domtblout file (default: auto-resolve from PFAM_ID)")
    parser.add_argument("-o", "--output-dir", type=str, default=None,
                        help="Output base directory for pep/ and cds/ subdirs (default: auto-resolve)")
    parser.add_argument("--data-dir", type=str, default=None,
                        help=f"Species data directory (default: {_DEF_DATA_DIR})")
    parser.add_argument("--skip-dirs", type=str, nargs="*", default=None,
                        help="Directories to skip when scanning species")
    parser.add_argument("--delimiter", type=str, default="|",
                        help="Delimiter between species and gene ID in domtblout (default: '|')")
    parser.add_argument("--skip-cds", action="store_true",
                        help="Skip CDS extraction (PEP only)")
    parser.add_argument("--force", action="store_true",
                        help="Re-extract even if output files already exist")
    parser.add_argument("--pep-suffix", type=str, default=None,
                        help="PEP file suffix pattern (auto-detect if not set)")
    parser.add_argument("--cds-suffix", type=str, default=None,
                        help="CDS file suffix pattern (auto-detect if not set)")
    return parser.parse_args()


def find_pep_cds_files(species_dir, pep_suffix=None, cds_suffix=None):
    """Find pep and cds files in a species directory."""
    pep_file = None
    cds_file = None
    for f in os.listdir(species_dir):
        fpath = os.path.join(species_dir, f)
        if os.path.isdir(fpath):
            continue
        fl = f.lower()

        if pep_suffix:
            if pep_suffix in f:
                pep_file = fpath
        elif 'pep' in fl and (f.endswith('.fa') or f.endswith('.faa')):
            if pep_file is None or '_pep_' in fl or '.pep.' in fl:
                pep_file = fpath

        if cds_suffix:
            if cds_suffix in f:
                cds_file = fpath
        elif ('cds' in fl or '_CDS_' in f) and (f.endswith('.fa') or f.endswith('.faa')):
            if cds_file is None or '_CDS_' in f or '.cds.' in fl:
                cds_file = fpath
    return pep_file, cds_file


def build_species_file_map(data_dir, skip_dirs, pep_suffix=None, cds_suffix=None):
    """Build mapping of species name -> (pep_path, cds_path)."""
    species_map = {}
    for entry in sorted(os.listdir(data_dir)):
        if entry in skip_dirs:
            continue
        sp_dir = os.path.join(data_dir, entry)
        if not os.path.isdir(sp_dir):
            continue
        pep_file, cds_file = find_pep_cds_files(sp_dir, pep_suffix, cds_suffix)
        if pep_file is None:
            print(f"WARNING: No pep file for {entry}")
            continue
        species_map[entry] = (pep_file, cds_file)
    return species_map


def parse_hmmer_targets(hmmer_file):
    """Parse HMMER domtblout and return sorted target IDs."""
    targets = set()
    with open(hmmer_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split()
            if cols:
                full_id = cols[0]
                if '|' in full_id:
                    _, _, target_id = full_id.partition('|')
                    targets.add(target_id)
                else:
                    targets.add(full_id)
    return sorted(targets)


def protein_to_transcript_candidates(protein_id):
    """Convert a protein ID to possible transcript IDs using common patterns."""
    candidates = []

    m = re.match(r'^(.+)_P(\d+)$', protein_id)
    if m:
        candidates.append(f"{m.group(1)}_t{m.group(2)}")
        candidates.append(f"{m.group(1)}_T{m.group(2)}")

    m = re.match(r'^(.+)_P(\d+)$', protein_id)
    if m:
        candidates.append(f"{m.group(1)}.t{m.group(2)}")
        candidates.append(f"{m.group(1)}.T{m.group(2)}")

    m = re.match(r'^(.+)\.\d+-P$', protein_id)
    if m:
        candidates.append(m.group(1))

    if protein_id.endswith(':cds'):
        candidates.append(protein_id[:-4])

    if protein_id.startswith('cds-'):
        candidates.append(protein_id[4:])

    m = re.match(r'^(.+)_p(\d+\.\d+)$', protein_id)
    if m:
        candidates.append(f"{m.group(1)}_t{m.group(2)}")
        candidates.append(f"{m.group(1)}_T{m.group(2)}")

    if protein_id.endswith('-P'):
        candidates.append(protein_id[:-2])

    m = re.match(r'^(.+)_P\d+$', protein_id)
    if m:
        candidates.append(m.group(1))

    return candidates


def build_cds_target_map(cds_file, hmmer_targets):
    """Build a mapping from HMMER protein target IDs to transcript IDs in CDS file."""
    transcript_ids = set()
    with open(cds_file) as f:
        for line in f:
            if line.startswith('>'):
                tid = line[1:].strip().split()[0]
                transcript_ids.add(tid)

    mapping = {}
    unmatched = []

    for pid in hmmer_targets:
        if pid in transcript_ids:
            mapping[pid] = pid
            continue

        found = False
        for candidate in protein_to_transcript_candidates(pid):
            if candidate in transcript_ids:
                mapping[pid] = candidate
                found = True
                break

        if not found:
            unmatched.append(pid)

    return mapping, unmatched


def main():
    args = parse_args()
    pfam_id = args.pfam_id.upper()

    # ── Resolve species data dir and skip dirs ──────────────────────────
    data_dir = args.data_dir or _DEF_DATA_DIR
    skip_dirs = set(args.skip_dirs) if args.skip_dirs else _DEF_SKIP_DIRS

    # ── Resolve input domtblout path ───────────────────────────────────
    if args.input is not None:
        hmmer_file = os.path.abspath(args.input)
    else:
        paths = resolve_pfam_paths(pfam_id)
        hmmer_file = str(paths["hmmer_file"])

    if not os.path.exists(hmmer_file):
        print(f"ERROR: domtblout not found: {hmmer_file}")
        print("  Use -i to specify a custom path, or run Step 1 first.")
        sys.exit(1)

    # ── Resolve output dirs ────────────────────────────────────────────
    if args.output_dir is not None:
        out_base = os.path.abspath(args.output_dir)
        out_pep_dir = os.path.join(out_base, "pep")
        out_cds_dir = os.path.join(out_base, "cds")
    else:
        paths = resolve_pfam_paths(pfam_id)
        out_pep_dir = str(paths["pep_dir"])
        out_cds_dir = str(paths["cds_dir"])
        out_base = str(paths["base_dir"])

    print("=" * 60)
    print(f"Extracting {pfam_id} domain sequences per species")
    print(f"  Input:   {hmmer_file}")
    print(f"  Output:  {out_base}/")
    print(f"  Data:    {data_dir}")
    print(f"  CDS:     {'enabled' if not args.skip_cds else 'skipped'}")
    print("=" * 60)

    # ── Step 1: Build species file map ─────────────────────────────────
    print("\n[Step 1] Building species file map...")
    species_map = build_species_file_map(
        data_dir, skip_dirs,
        pep_suffix=args.pep_suffix, cds_suffix=args.cds_suffix
    )
    print(f"  Found {len(species_map)} species with PEP files")

    # ── Step 2: Parse HMMER targets ────────────────────────────────────
    print("\n[Step 2] Parsing HMMER targets...")
    targets = parse_hmmer_targets(hmmer_file)
    print(f"  Found {len(targets)} unique target IDs")

    target_file = f"/tmp/all_{pfam_id}_target_ids.txt"
    with open(target_file, 'w') as f:
        for t in targets:
            f.write(t + '\n')

    # ── Step 3: Extract sequences per species ─────────────────────────
    print("\n[Step 3] Extracting sequences per species...")
    os.makedirs(out_cds_dir, exist_ok=True)
    os.makedirs(out_pep_dir, exist_ok=True)

    total_extracted = 0
    species_counts = {}

    for sp_name, (pep_file, cds_file) in sorted(species_map.items()):
        out_pep = os.path.join(out_pep_dir, f"{sp_name}.pep.fa")

        # Skip if output already exists and not --force
        if os.path.exists(out_pep) and os.path.getsize(out_pep) > 0 and not args.force:
            r = subprocess.run(["grep", "-c", "^>", out_pep], capture_output=True, text=True)
            pep_count = int(r.stdout.strip() or 0)
            print(f"\n  [{sp_name}]")
            print(f"    PEP: {out_pep} (exists, {pep_count} seqs, skip)")
        else:
            print(f"\n  [{sp_name}]")
            print(f"    PEP: {pep_file}")

            result = subprocess.run(
                ["seqkit", "grep", "-f", target_file, "-o", out_pep, pep_file],
                capture_output=True, text=True,
            )
            if result.returncode != 0:
                print(f"    ERROR (PEP): {result.stderr.strip()}")
                continue

            pep_count = 0
            if os.path.exists(out_pep) and os.path.getsize(out_pep) > 0:
                r = subprocess.run(["grep", "-c", "^>", out_pep], capture_output=True, text=True)
                pep_count = int(r.stdout.strip() or 0)
            print(f"    PEP extracted: {pep_count}")

        # ── CDS extraction (optional) ─────────────────────────────────
        cds_count = 0
        if not args.skip_cds and cds_file and os.path.exists(cds_file):
            out_cds = os.path.join(out_cds_dir, f"{sp_name}.cds.fa")

            if os.path.exists(out_cds) and os.path.getsize(out_cds) > 0 and not args.force:
                r = subprocess.run(["grep", "-c", "^>", out_cds], capture_output=True, text=True)
                cds_count = int(r.stdout.strip() or 0)
                print(f"    CDS: {out_cds} (exists, {cds_count} seqs, skip)")
            else:
                print(f"    CDS: {cds_file}")

                mapping, unmatched = build_cds_target_map(cds_file, targets)

                if mapping:
                    cds_target_file = f"/tmp/{pfam_id}_cds_targets_{sp_name}.txt"
                    with open(cds_target_file, 'w') as f:
                        for tid in mapping.values():
                            f.write(tid + '\n')

                    result = subprocess.run(
                        ["seqkit", "grep", "-f", cds_target_file, "-o", out_cds, cds_file],
                        capture_output=True, text=True,
                    )
                    if result.returncode == 0 and os.path.exists(out_cds) and os.path.getsize(out_cds) > 0:
                        r = subprocess.run(["grep", "-c", "^>", out_cds], capture_output=True, text=True)
                        cds_count = int(r.stdout.strip() or 0)

                    os.remove(cds_target_file)

                if unmatched:
                    print(f"    CDS unmatched IDs: {len(unmatched)} (e.g., {unmatched[:3]})")
                print(f"    CDS extracted: {cds_count} (mapped {len(mapping)} protein IDs)")
        elif not args.skip_cds:
            print(f"    CDS: not available")

        species_counts[sp_name] = (pep_count, cds_count)
        total_extracted += pep_count

    # ── Summary ────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Summary:")
    for sp_name, (pc, cc) in sorted(species_counts.items()):
        if pc > 0:
            print(f"  {sp_name}: PEP={pc}, CDS={cc}")
    print(f"\nTotal sequences extracted: {total_extracted}")
    print(f"Output PEP: {out_pep_dir}/")
    print(f"Output CDS: {out_cds_dir}/")
    print(f"\nNext step: python 03_build_domain_tree.py {pfam_id}")
    print("=" * 60)


if __name__ == "__main__":
    main()
