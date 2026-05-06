#!/usr/bin/env python3
"""
[Step 4] Build a phylogenetic tree of domain-containing proteins across species.

By default, uses ALL domain-containing protein sequences from all species.
Optionally, use --strategy longest to pick one representative per species.

Usage:
    python 04_build_domain_tree.py <PFAM_ID>
    python 04_build_domain_tree.py <PFAM_ID> -i /path/to/pep_dir -o /path/to/output
    python 04_build_domain_tree.py PF00168 --strategy longest
    python 04_build_domain_tree.py PF00168 --mafft-custom "mafft --localpair --maxiterate 1000 --thread 32"
    python 04_build_domain_tree.py PF00168 --iqtree-custom "iqtree2 -s {input} --prefix {prefix} -m MFP -B 1000 -T AUTO"

Input:
    {pep_dir}/*.pep.fa           (from 03_extract_domain_seqs.py)
Output:
    {align_dir}/{pfam_id}_aligned.fa
    {tree_dir}/{pfam_id}_tree.treefile
    {tree_dir}/{pfam_id}_tree.contree
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path

# ── Path setup ──────────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)

try:
    from config.species_config import resolve_pfam_paths
except ModuleNotFoundError:
    from species_config import resolve_pfam_paths


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build phylogenetic tree of domain-containing proteins across species",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  %(prog)s PF00168
  %(prog)s PF00168 -i /data/my_project/pep -o /data/my_project/tree
  %(prog)s PF00168 --strategy longest
  %(prog)s PF00168 --skip-align --skip-tree   (just collect sequences)
  %(prog)s PF00168 --mafft-custom "mafft --localpair --maxiterate 1000 --thread 32"
  %(prog)s PF00168 --iqtree-custom "iqtree2 -s {input} --prefix {prefix} -m MFP -B 1000 -T AUTO"
  %(prog)s PF00168 --iqtree-args "-m MFP -B 1000 -T AUTO --no-terrace"
        """,
    )
    parser.add_argument("pfam_id", type=str, help="Pfam accession (e.g. PF00168)")
    parser.add_argument("-i", "--input-dir", type=str, default=None,
                        help="PEP directory containing *.pep.fa (default: auto-resolve)")
    parser.add_argument("-o", "--output-dir", type=str, default=None,
                        help="Output base directory (default: auto-resolve)")
    parser.add_argument("--strategy", type=str, default="all", choices=["all", "longest"],
                        help="'all' = use every sequence; 'longest' = one per species (default: all)")
    parser.add_argument("--align-dir", type=str, default=None,
                        help="Custom alignment output directory (default: {output_dir}/alignment)")
    parser.add_argument("--tree-dir", type=str, default=None,
                        help="Custom tree output directory (default: {output_dir}/tree)")
    parser.add_argument("--skip-align", action="store_true",
                        help="Skip MAFFT alignment (requires existing aligned FASTA)")
    parser.add_argument("--skip-tree", action="store_true",
                        help="Skip IQ-TREE2 (alignment only)")

    # ── MAFFT control ─────────────────────────────────────────────────────
    parser.add_argument("--mafft-custom", type=str, default=None,
                        help="Full MAFFT command. Use {input} for rep FASTA path. "
                             "E.g.: 'mafft --localpair --maxiterate 1000 --thread 32 {input}'")
    parser.add_argument("--mafft-args", type=str, default="--auto",
                        help="MAFFT arguments (ignored if --mafft-custom set). Default: '--auto'")
    parser.add_argument("--mafft-threads", type=str, default="-1",
                        help="MAFFT threads (ignored if --mafft-custom set). Default: -1 = all")

    # ── IQ-TREE2 control ──────────────────────────────────────────────────
    parser.add_argument("--iqtree-custom", type=str, default=None,
                        help="Full IQ-TREE2 command. Use {input} for aligned FASTA, {prefix} for tree prefix. "
                             "E.g.: 'iqtree2 -s {input} --prefix {prefix} -m MFP -B 1000 -T AUTO'")
    parser.add_argument("--iqtree-args", type=str, default="",
                        help="Extra IQ-TREE2 arguments (appended to default cmd, ignored if --iqtree-custom set)")

    parser.add_argument("--force", action="store_true",
                        help="Re-run even if output files exist")
    parser.add_argument("--delimiter", type=str, default="|",
                        help="Delimiter in FASTA header between species and gene (default: '|')")
    return parser.parse_args()


def collect_sequences(pep_dir, strategy):
    """Collect domain sequences from pep_dir. Returns list of (species, gene_id, seq, length)."""
    pep_files = sorted(Path(pep_dir).glob("*.pep.fa"))
    if not pep_files:
        print(f"ERROR: No *.pep.fa files found in {pep_dir}")
        sys.exit(1)

    print(f"  Found {len(pep_files)} PEP files")

    records = []

    if strategy == "longest":
        for fpath in pep_files:
            species = fpath.stem.replace(".pep", "")
            longest_id, longest_seq, longest_len = None, None, -1
            cur_id, cur_lines = None, []

            with open(fpath) as f:
                for line in f:
                    line = line.rstrip("\n")
                    if line.startswith(">"):
                        if cur_id is not None and cur_lines:
                            seq = "".join(cur_lines)
                            if len(seq) > longest_len:
                                longest_len = len(seq)
                                longest_id = cur_id
                                longest_seq = seq
                        cur_id = line[1:].split()[0]
                        cur_lines = []
                    else:
                        cur_lines.append(line)
                if cur_id is not None and cur_lines:
                    seq = "".join(cur_lines)
                    if len(seq) > longest_len:
                        longest_len = len(seq)
                        longest_id = cur_id
                        longest_seq = seq

            if longest_id is not None:
                records.append((species, longest_id, longest_seq, longest_len))
                print(f"  {species:30s} \u2192 {longest_id:40s}  ({longest_len} aa) [LONGEST]")
            else:
                print(f"  {species:30s} \u2192 NO SEQUENCE", file=sys.stderr)
    else:
        for fpath in pep_files:
            species = fpath.stem.replace(".pep", "")
            count = 0
            with open(fpath) as f:
                for line in f:
                    if line.startswith(">"):
                        header = line[1:].strip().split()[0]
                        records.append((species, header, None, None))
                        count += 1
            print(f"  {species:30s} \u2192 {count} sequences")

    return records


def write_representative_fasta(records, out_path, delimiter="|"):
    """Write multi-FASTA for 'longest' strategy (seqs already loaded)."""
    count = 0
    with open(out_path, "w") as out:
        for species, gid, seq, _slen in records:
            out.write(f">{species}{delimiter}{gid}\n{seq}\n")
            count += 1
    return count


def load_all_sequences(pep_dir, records, delimiter="|"):
    """Load actual sequences for all records (used in 'all' strategy)."""
    from collections import defaultdict
    species_genes = defaultdict(set)
    for sp, gid, _, _ in records:
        species_genes[sp].add(gid)

    pep_files = sorted(Path(pep_dir).glob("*.pep.fa"))
    loaded = []

    for fpath in pep_files:
        species = fpath.stem.replace(".pep", "")
        wanted = species_genes.get(species, set())
        if not wanted:
            continue
        cur_id, cur_lines = None, []
        with open(fpath) as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if cur_id is not None and cur_lines and cur_id in wanted:
                        loaded.append((species, cur_id, "".join(cur_lines)))
                    cur_id = line[1:].split()[0]
                    cur_lines = []
                else:
                    cur_lines.append(line)
            if cur_id is not None and cur_lines and cur_id in wanted:
                loaded.append((species, cur_id, "".join(cur_lines)))

    return loaded


def shell_split(s):
    """Split a command string into argv list, respecting quotes."""
    import shlex
    return shlex.split(s)


def main():
    args = parse_args()
    pfam_id = args.pfam_id.upper()

    # ── Resolve paths ─────────────────────────────────────────────────────
    if args.input_dir is not None:
        pep_dir = Path(os.path.abspath(args.input_dir))
    else:
        paths = resolve_pfam_paths(pfam_id)
        pep_dir = paths["pep_dir"]

    if args.output_dir is not None:
        out_base = Path(os.path.abspath(args.output_dir))
        align_dir = Path(os.path.abspath(args.align_dir)) if args.align_dir else (out_base / "alignment")
        tree_dir = Path(os.path.abspath(args.tree_dir)) if args.tree_dir else (out_base / "tree")
    else:
        paths = resolve_pfam_paths(pfam_id)
        align_dir = paths["align_dir"]
        tree_dir = paths["tree_dir"]
        out_base = paths["base_dir"]

    rep_fasta = tree_dir / f"{pfam_id}_representatives.fa"
    aligned_fasta = align_dir / f"{pfam_id}_aligned.fa"
    tree_prefix = tree_dir / f"{pfam_id}_tree"

    align_dir.mkdir(parents=True, exist_ok=True)
    tree_dir.mkdir(parents=True, exist_ok=True)

    if not pep_dir.exists():
        print(f"ERROR: PEP directory not found: {pep_dir}")
        print(f"  Run: python 03_extract_domain_seqs.py {pfam_id}")
        sys.exit(1)

    # ── Step 1: Collect sequences ─────────────────────────────────────────
    print("=" * 60)
    print(f"[Step 4] Building tree for {pfam_id} (strategy: {args.strategy})")
    print("=" * 60)

    print(f"\n[1/4] Collecting sequences from {pep_dir} ...")
    records = collect_sequences(pep_dir, args.strategy)
    print(f"  Total records: {len(records)}")

    # ── Step 2: Write representative FASTA ────────────────────────────────
    print(f"\n[2/4] Writing representative FASTA to {rep_fasta} ...")

    if args.strategy == "longest":
        n = write_representative_fasta(records, rep_fasta, args.delimiter)
        print(f"  Written: {n} sequences")
    else:
        loaded = load_all_sequences(pep_dir, records, args.delimiter)
        with open(rep_fasta, "w") as out:
            for species, gid, seq in loaded:
                out.write(f">{species}{args.delimiter}{gid}\n{seq}\n")
        print(f"  Written: {len(loaded)} sequences")

    if args.skip_align and args.skip_tree:
        print("\nDone! (--skip-align and --skip-tree, representative FASTA only)")
        return

    # ── Step 3: MAFFT alignment ──────────────────────────────────────────
    if not args.skip_align:
        print(f"\n[3/4] Aligning with MAFFT ...")

        if aligned_fasta.exists() and aligned_fasta.stat().st_size > 0 and not args.force:
            print(f"  Alignment already exists: {aligned_fasta} (use --force to redo)")
        else:
            if args.mafft_custom:
                cmd_str = args.mafft_custom.replace("{input}", str(rep_fasta))
                mafft_cmd = shell_split(cmd_str)
                if "{input}" not in args.mafft_custom and str(rep_fasta) not in mafft_cmd:
                    mafft_cmd.append(str(rep_fasta))
                print(f"  Custom command: {' '.join(mafft_cmd)}")
            else:
                mafft_cmd = ["mafft"] + shell_split(args.mafft_args) + ["--thread", args.mafft_threads, str(rep_fasta)]
                print(f"  Command: {' '.join(mafft_cmd)}")

            with open(aligned_fasta, "w") as out:
                result = subprocess.run(mafft_cmd, stdout=out, stderr=subprocess.PIPE, text=True)

            if result.returncode != 0:
                print(f"  MAFFT failed (rc={result.returncode})", file=sys.stderr)
                print(f"  stderr: {result.stderr[:500]}", file=sys.stderr)
                sys.exit(1)

            n_seqs = 0
            with open(aligned_fasta) as f:
                for line in f:
                    if line.startswith(">"):
                        n_seqs += 1
            print(f"  Alignment: {aligned_fasta} ({n_seqs} sequences)")
    else:
        print(f"\n[3/4] Skipping MAFFT (--skip-align)")
        if not aligned_fasta.exists():
            print(f"  WARNING: aligned FASTA not found: {aligned_fasta}")

    # ── Step 4: IQ-TREE2 ────────────────────────────────────────────────
    if not args.skip_tree:
        print(f"\n[4/4] Running IQ-TREE2 ...")

        if not aligned_fasta.exists():
            print(f"  ERROR: aligned FASTA not found: {aligned_fasta}")
            print(f"  Remove --skip-align or run without it.")
            sys.exit(1)

        treefile = tree_prefix.with_suffix(".treefile")
        if treefile.exists() and treefile.stat().st_size > 0 and not args.force:
            print(f"  Tree already exists: {treefile} (use --force to redo)")
        else:
            if args.iqtree_custom:
                cmd_str = args.iqtree_custom
                cmd_str = cmd_str.replace("{input}", str(aligned_fasta))
                cmd_str = cmd_str.replace("{prefix}", str(tree_prefix))
                iqtree_cmd = shell_split(cmd_str)
                if "-s" not in iqtree_cmd:
                    iqtree_cmd += ["-s", str(aligned_fasta)]
                print(f"  Custom command: {' '.join(iqtree_cmd)}")
            else:
                iqtree_cmd = [
                    "/usr/bin/iqtree2",
                    "-s", str(aligned_fasta),
                    "--prefix", str(tree_prefix),
                    "-m", "MFP",
                    "-B", "1000",
                    "-T", "AUTO",
                    "--redo",
                ]
                if args.iqtree_args:
                    iqtree_cmd += shell_split(args.iqtree_args)
                print(f"  Command: {' '.join(iqtree_cmd)}")

            result = subprocess.run(iqtree_cmd, capture_output=True, text=True)

            if result.returncode != 0:
                print(f"  IQ-TREE2 failed (rc={result.returncode})", file=sys.stderr)
                print(f"  stderr: {result.stderr[:500]}", file=sys.stderr)
                sys.exit(1)

            lines = result.stdout.strip().split("\n")
            print("  IQ-TREE2 output (last 15 lines):")
            for line in lines[-15:]:
                print(f"    {line}")
    else:
        print(f"\n[4/4] Skipping IQ-TREE2 (--skip-tree)")

    # ── Summary ──────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print(f"SUMMARY \u2014 {pfam_id} ({args.strategy})")
    print("=" * 60)
    print(f"  Representative FASTA: {rep_fasta}")
    print(f"  Alignment:           {aligned_fasta}")
    print(f"  Tree file:           {tree_prefix}.treefile")
    print(f"  Consensus tree:      {tree_prefix}.contree")
    print(f"  IQ-TREE2 report:     {tree_prefix}.iqtree")
    print("=" * 60)
    print("Done!")


if __name__ == "__main__":
    main()
