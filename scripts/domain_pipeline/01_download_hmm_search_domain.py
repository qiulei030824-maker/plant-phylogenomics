#!/usr/bin/env python3
"""
[Step 1] Download Pfam HMM for a given PFAM ID and run hmmsearch on all species.

Usage:
    python 01_download_hmm_search_domain.py <PFAM_ID>
    python 01_download_hmm_search_domain.py <PFAM_ID> -i /path/to/genomes -o /path/to/output
    python 01_download_hmm_search_domain.py <PFAM_ID> --hmmsearch /usr/bin/hmmsearch --force-redownload

Examples:
    python 01_download_hmm_search_domain.py PF00168                          # default paths
    python 01_download_hmm_search_domain.py PF00168 -i /data/genomes -o /data/results
    python 01_download_hmm_search_domain.py PF00168 --skip-hmmsearch         # HMM only
    python 01_download_hmm_search_domain.py C2

This script:
  1. Downloads the Pfam HMM from the Pfam website
  2. Runs hmmsearch on all species' proteome (PEP) files
  3. Outputs a combined domtblout file
"""

import gzip
import os
import sys
import argparse
import urllib.request
import shutil

# Import shared config (used as defaults)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from species_config import DATA_DIR, BASE_ROOT, SKIP_DIRS, resolve_pfam_paths


# ── Constants (mirrored from species_config for custom output paths) ──────────
C2_PEP_DIRNAME = "C2domain_pep"
C2_CDS_DIRNAME = "C2domain_cds"
C2_HMMER_DIRNAME = "hmmer"
C2_TREE_DIRNAME = "tree"
C2_ALIGN_DIRNAME = "algin"
C2_HMMER_FILENAME = "all_C2.domtblout"
C2_TREE_FILENAME = "C2_tree.treefile"


# ── Helper functions ──────────────────────────────────────────────────────────

def find_pep_file(species_dir):
    """Find the PEP (proteome) file in a species directory."""
    for f in os.listdir(species_dir):
        fpath = os.path.join(species_dir, f)
        if os.path.isdir(fpath):
            continue
        fl = f.lower()
        if 'pep' in fl and (f.endswith('.fa') or f.endswith('.faa')):
            return fpath
    return None


def resolve_output_paths(pfam_id, output_root):
    """
    Resolve output paths for a given PFAM_ID under a custom output_root.
    Mirrors the logic from species_config.resolve_pfam_paths.
    """
    is_c2 = pfam_id.upper() == "C2"
    if is_c2:
        base_dir = output_root
        hmmer_dir = os.path.join(output_root, C2_HMMER_DIRNAME)
        pep_dir = os.path.join(output_root, C2_PEP_DIRNAME)
        cds_dir = os.path.join(output_root, C2_CDS_DIRNAME)
        tree_dir = os.path.join(output_root, C2_TREE_DIRNAME)
        align_dir = os.path.join(output_root, C2_ALIGN_DIRNAME)
        hmmer_file = os.path.join(hmmer_dir, C2_HMMER_FILENAME)
        tree_file = os.path.join(tree_dir, C2_TREE_FILENAME)
    else:
        base_dir = os.path.join(output_root, pfam_id)
        hmmer_dir = os.path.join(base_dir, "hmmer")
        pep_dir = os.path.join(base_dir, "pep")
        cds_dir = os.path.join(base_dir, "cds")
        tree_dir = os.path.join(base_dir, "tree")
        align_dir = os.path.join(base_dir, "algin")
        hmmer_file = os.path.join(hmmer_dir, f"all_{pfam_id}.domtblout")
        tree_file = os.path.join(tree_dir, f"{pfam_id}_tree.treefile")
    return {
        "base_dir": base_dir,
        "hmmer_dir": hmmer_dir,
        "pep_dir": pep_dir,
        "cds_dir": cds_dir,
        "tree_dir": tree_dir,
        "align_dir": align_dir,
        "hmmer_file": hmmer_file,
        "tree_file": tree_file,
        "is_c2": is_c2,
    }


def download_pfam_hmm(pfam_id, output_path):
    """Download Pfam HMM from the Pfam website."""
    url = f"https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{pfam_id}?annotation=hmm"
    print(f"  Downloading HMM from: {url}")
    try:
        tmp_path = output_path + ".tmp"
        urllib.request.urlretrieve(url, tmp_path)
        if not os.path.exists(tmp_path) or os.path.getsize(tmp_path) == 0:
            print(f"  ERROR: Downloaded file is empty")
            return False

        with open(tmp_path, 'rb') as f:
            magic = f.read(2)

        if magic == b'\x1f\x8b':
            print(f"  Detected gzip compression, decompressing...")
            with gzip.open(tmp_path, 'rt', encoding='utf-8', errors='replace') as gz_in:
                content = gz_in.read()
            with open(output_path, 'w') as f_out:
                f_out.write(content)
            os.remove(tmp_path)
            print(f"  Decompressed and saved to: {output_path}")
        else:
            os.rename(tmp_path, output_path)
            print(f"  Saved to: {output_path}")

        return True
    except Exception as e:
        print(f"  ERROR downloading HMM: {e}")
        return False


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Download Pfam HMM and run hmmsearch on species proteomes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s PF00168
  %(prog)s PF00168 -i /custom/genome/dir -o /custom/output/root
  %(prog)s PF00168 -i /custom/genome/dir -o /custom/output/root --hmmsearch /usr/local/bin/hmmsearch
  %(prog)s PF00168 --skip-hmmsearch                # download HMM only
        """,
    )
    parser.add_argument("pfam_id", type=str, help="Pfam accession (e.g. PF00168, C2)")
    parser.add_argument("-i", "--input-dir", type=str, default=None,
                        help="Directory containing species PEP subdirectories")
    parser.add_argument("-o", "--output-dir", type=str, default=None,
                        help="Output root directory")
    parser.add_argument("--hmmsearch", type=str, default=None,
                        help="Path to hmmsearch binary (default: auto-detect from PATH)")
    parser.add_argument("--skip-dirs", type=str, nargs="*", default=None,
                        help="Additional directory names to skip (space separated)")
    parser.add_argument("--skip-hmmsearch", action="store_true",
                        help="Only download HMM, skip hmmsearch step")
    parser.add_argument("--force-redownload", action="store_true",
                        help="Re-download HMM even if it already exists")
    parser.add_argument("--pep-suffix", type=str, default=None,
                        help="Custom substring to identify PEP files (default: 'pep')")
    return parser.parse_args()


def find_species_pep_files(input_dir, skip_dirs, pep_suffix=None):
    """
    Scan input_dir for species subdirectories containing PEP files.
    Returns dict of {species_name: pep_file_path}.
    """
    species_pep = {}
    for entry in sorted(os.listdir(input_dir)):
        if entry in skip_dirs:
            continue
        sp_dir = os.path.join(input_dir, entry)
        if not os.path.isdir(sp_dir):
            continue

        if pep_suffix:
            for f in os.listdir(sp_dir):
                fpath = os.path.join(sp_dir, f)
                if os.path.isfile(fpath) and pep_suffix in f:
                    species_pep[entry] = fpath
                    break
        else:
            pep_file = find_pep_file(sp_dir)
            if pep_file is not None:
                species_pep[entry] = pep_file
    return species_pep


def main():
    args = parse_args()

    pfam_id = args.pfam_id.upper()

    # ── Resolve input directory ─────────────────────────────────────────────
    input_dir = os.path.abspath(args.input_dir) if args.input_dir else DATA_DIR

    # ── Resolve output paths ────────────────────────────────────────────────
    if args.output_dir is not None:
        output_root = os.path.abspath(args.output_dir)
        paths = resolve_output_paths(pfam_id, output_root)
    else:
        output_root = str(BASE_ROOT)
        paths = resolve_pfam_paths(pfam_id)

    HMMER_DIR = paths["hmmer_dir"]

    # ── Resolve hmmsearch binary ────────────────────────────────────────────
    if args.hmmsearch is not None:
        HMMSEARCH = args.hmmsearch
    else:
        HMMSEARCH = shutil.which("hmmsearch")
        if HMMSEARCH is None:
            HMMSEARCH = "/usr/bin/hmmsearch"

    # ── Resolve skip dirs ────────────────────────────────────────────────────
    skip_dirs = set(SKIP_DIRS)
    if args.skip_dirs:
        skip_dirs.update(args.skip_dirs)

    print("=" * 60)
    print(f"Step 1: Download Pfam HMM and run hmmsearch for {pfam_id}")
    print(f"  Input dir:   {input_dir}")
    print(f"  Output dir:  {paths['base_dir']}/")
    print(f"  hmmsearch:   {HMMSEARCH}")
    print("=" * 60)

    # ── Step 1: Create output directories ──────────────────────────────────
    os.makedirs(HMMER_DIR, exist_ok=True)

    # ── Step 2: Download Pfam HMM ──────────────────────────────────────────
    print(f"\n[Step 2] Downloading Pfam HMM for {pfam_id}...")
    hmm_file = os.path.join(HMMER_DIR, f"{pfam_id}.hmm")
    if os.path.exists(hmm_file) and os.path.getsize(hmm_file) > 0 and not args.force_redownload:
        print(f"  HMM file already exists: {hmm_file}")
    else:
        if os.path.exists(hmm_file) and args.force_redownload:
            print(f"  Force re-download: removing existing {hmm_file}")
            os.remove(hmm_file)
        success = download_pfam_hmm(pfam_id, hmm_file)
        if not success:
            url2 = f"http://pfam.xfam.org/family/{pfam_id}/hmm"
            print(f"  Trying alternative URL: {url2}")
            try:
                urllib.request.urlretrieve(url2, hmm_file)
                if os.path.getsize(hmm_file) > 0:
                    print(f"  Saved to: {hmm_file}")
                else:
                    print(f"  ERROR: Downloaded file is empty")
                    sys.exit(1)
            except Exception as e:
                print(f"  ERROR downloading HMM: {e}")
                sys.exit(1)

    # ── Skip hmmsearch if requested ────────────────────────────────────────
    if args.skip_hmmsearch:
        print("\n  --skip-hmmsearch set, skipping hmmsearch.")
        print(f"\n  Pfam HMM: {hmm_file}")
        print("Done!")
        return

    # ── Step 3: Find all species with PEP files ────────────────────────────
    print(f"\n[Step 3] Finding species with PEP files in {input_dir}...")
    if args.pep_suffix:
        species_pep = find_species_pep_files(input_dir, skip_dirs, pep_suffix=args.pep_suffix)
    else:
        species_pep = find_species_pep_files(input_dir, skip_dirs)
    print(f"  Found {len(species_pep)} species with PEP files")

    # ── Step 4: Run hmmsearch on concatenated PEPs ──────────────────────────
    print(f"\n[Step 4] Running hmmsearch on all species...")

    combined_domtblout = os.path.join(HMMER_DIR, f"all_{pfam_id}.domtblout")
    log_file = os.path.join(HMMER_DIR, "hmmsearch.log")

    combined_pep = os.path.join(HMMER_DIR, f"all_proteomes_{pfam_id}.fa")
    if not os.path.exists(combined_pep):
        print(f"  Concatenating all PEP files into {combined_pep}...")
        with open(combined_pep, 'w') as out:
            for sp_name, pep_file in sorted(species_pep.items()):
                with open(pep_file) as f:
                    for line in f:
                        if line.startswith('>'):
                            out.write(f">{sp_name}|{line[1:]}")
                        else:
                            out.write(line)
        print(f"  Combined PEP file created: {combined_pep}")

    if not os.path.exists(combined_domtblout) or os.path.getsize(combined_domtblout) == 0:
        cmd = [
            HMMSEARCH,
            "--domtblout", combined_domtblout,
            "-o", "/dev/null",
            hmm_file,
            combined_pep
        ]
        print(f"  Command: {' '.join(cmd)}")
        print("  Running hmmsearch...")
        cmd_str = " ".join(cmd) + f" > {log_file} 2>&1"
        retcode = os.system(cmd_str)
        if retcode != 0:
            print(f"  ERROR: hmmsearch failed with return code {retcode}")
            sys.exit(1)
    else:
        print(f"  hmmsearch results already exist: {combined_domtblout}")

    # ── Step 5: Count results ──────────────────────────────────────────────
    print(f"\n[Step 5] Counting results...")
    target_count = 0
    species_hits = {}
    with open(combined_domtblout) as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split()
            if len(cols) < 5:
                continue
            target_id = cols[0]
            target_count += 1
            if '|' in target_id:
                sp = target_id.split('|')[0]
                species_hits[sp] = species_hits.get(sp, 0) + 1

    print(f"  Total domain hits: {target_count}")
    print(f"  Species with hits: {len(species_hits)}")
    for sp, count in sorted(species_hits.items(), key=lambda x: -x[1]):
        print(f"    {sp}: {count}")

    # ── Summary ────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Pfam HMM: {hmm_file}")
    print(f"  Combined domtblout: {combined_domtblout}")
    print(f"  Total hits: {target_count}")
    print(f"\nNext step: python 02_extract_domain_seqs.py {pfam_id}")
    print("=" * 60)


if __name__ == "__main__":
    main()
