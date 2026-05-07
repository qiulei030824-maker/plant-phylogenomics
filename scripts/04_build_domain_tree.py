#!/usr/bin/env python3
"""
[Step 4] Build phylogenetic tree of domain-containing proteins — Industrial-grade.

Features:
  - Profile system (fast|standard|accurate|ultra)
  - Sequence selection strategies (all|longest|canonical|domain_best|longest_isoform|representative)
  - Auto strategy engine (auto-detects <50 / 50-200 / 200-1000 / >1000 sequences)
  - MAFFT alignment (auto|linsi|einsi|fftns2) + trimming (clipkit|trimal|none)
  - IQ-TREE2 with MFP, UFBoot2, or fast LG
  - Taxonomy filtering via species_config (--clade / --family / --genus / --species)
  - Resume/checkpoint system
  - Metadata + QC report + auto-generated methods paragraph
  - Tree visualization (PDF/SVG/PNG)
  - Config file support (--config config.yaml)
  - Full logging

Usage:
    python 04_build_domain_tree.py <PFAM_ID>
    python 04_build_domain_tree.py <PFAM_ID> --profile accurate
    python 04_build_domain_tree.py <PFAM_ID> --strategy longest --species ChineseLong,DHL92
    python 04_build_domain_tree.py <PFAM_ID> --config config.yaml
"""

import os, sys, re, json, time, shutil, textwrap, subprocess, argparse
from pathlib import Path
from collections import defaultdict
from datetime import datetime

# ── Path setup ──────────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _SCRIPT_DIR)

try:
    from config.species_config import resolve_pfam_paths, TAXONOMY_HIERARCHY, get_all_species, get_group
except ModuleNotFoundError:
    from species_config import resolve_pfam_paths, TAXONOMY_HIERARCHY, get_all_species, get_group

# ── Logging ─────────────────────────────────────────────────────────────────
import logging

def setup_logger(log_dir, name="04_build_domain_tree", level=logging.INFO):
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.handlers.clear()

    fh = logging.FileHandler(log_dir / "run.log", mode="a")
    fh.setLevel(level)
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logger.addHandler(fh)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    logger.addHandler(ch)
    return logger

# ── Profiles ────────────────────────────────────────────────────────────────
PROFILES = {
    "fast":     {"mafft": "auto",        "trim": "relaxed",  "bootstrap": 100,  "model": "LG",     "threads": "auto"},
    "standard": {"mafft": "localpair",   "trim": "moderate", "bootstrap": 1000, "model": "MFP",    "threads": "auto"},
    "accurate": {"mafft": "linsi",       "trim": "strict",   "bootstrap": 2000, "model": "MFP+MERGE", "threads": "auto"},
    "ultra":    {"mafft": "einsi",       "trim": "adaptive", "bootstrap": 5000, "model": "MFP+MERGE", "threads": "auto"},
}

TRIM_PROGRAMS = {
    "clipkit": {"relaxed":  "-m kpi-smart-gap",  "moderate": "-m smart-gap",   "strict": "-m strict",   "adaptive": "-m adaptive"},
    "trimal":  {"relaxed":  "-gappyout",          "moderate": "-automated1",   "strict": "-strictplus", "adaptive": "-automated1"},
}

# ── Species resolution (shared with 03) ─────────────────────────────────────
def list_species():
    print("=" * 70)
    print("AVAILABLE TAXONOMY HIERARCHY FOR SPECIES SELECTION")
    print("=" * 70)
    for clade, families in TAXONOMY_HIERARCHY.items():
        print(f"\n  [{clade}]")
        for family, genera in families.items():
            if not genera:
                continue
            print(f"    Family: {family}")
            for genus, species_list in genera.items():
                for sp in species_list:
                    print(f"      {genus:15s} → {sp}")
    print("\n" + "=" * 70)
    print(f"Total species: {len(get_all_species())}")

def resolve_species(args):
    all_set = set(get_all_species())
    if args.species:
        selected = set(s.strip() for s in args.species.split(","))
        unknown = selected - all_set
        if unknown:
            print(f"WARNING: Unknown species: {', '.join(sorted(unknown))}")
        return selected & all_set

    selected = set()
    if args.clade:
        if args.clade not in TAXONOMY_HIERARCHY:
            print(f"ERROR: Unknown clade '{args.clade}'. Use --list to see options.")
            sys.exit(1)
        for families in TAXONOMY_HIERARCHY[args.clade].values():
            for genus_species in families.values():
                selected.update(genus_species)
        print(f"  Selected clade '{args.clade}': {len(selected)} species")
    if args.family:
        found = False
        for clade_data in TAXONOMY_HIERARCHY.values():
            if args.family in clade_data:
                for genus_species in clade_data[args.family].values():
                    selected.update(genus_species)
                found = True
                print(f"  Selected family '{args.family}': {len(selected)} species")
                break
        if not found:
            print(f"ERROR: Unknown family '{args.family}'. Use --list.")
            sys.exit(1)
    if args.genus:
        found = False
        for clade_data in TAXONOMY_HIERARCHY.values():
            for family_data in clade_data.values():
                if args.genus in family_data:
                    selected.update(family_data[args.genus])
                    found = True
        if found:
            print(f"  Selected genus '{args.genus}': {len(selected)} species")
        else:
            print(f"ERROR: Unknown genus '{args.genus}'. Use --list.")
            sys.exit(1)
    return selected

# ── Shell helper ────────────────────────────────────────────────────────────
def shell_split(s):
    import shlex; return shlex.split(s)

def run_cmd(cmd, logger=None, desc=None, cwd=None, check=True, timeout=None):
    s = " ".join(str(c) for c in cmd)
    if logger:
        logger.info(f"Running: {s}")
    if desc:
        print(f"  {desc} ...")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd, timeout=timeout)
    except subprocess.TimeoutExpired:
        msg = f"Command timed out ({timeout}s): {s}"
        if logger: logger.error(msg)
        print(f"  ERROR: {msg}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        msg = f"Command not found: {cmd[0]}. Is it installed?"
        if logger: logger.error(msg)
        print(f"  ERROR: {msg}", file=sys.stderr)
        sys.exit(1)
    if result.returncode != 0:
        msg = f"Command failed (rc={result.returncode}): {s}"
        if logger:
            logger.error(msg)
            logger.error(f"stderr: {result.stderr[:1000]}")
        print(f"  ERROR: {msg}", file=sys.stderr)
        if result.stderr:
            print(f"  stderr: {result.stderr[:500]}", file=sys.stderr)
        print(f"  Suggestion: Check inputs and tool installation.", file=sys.stderr)
        if check:
            sys.exit(1)
    return result

# ── Auto strategy engine ────────────────────────────────────────────────────
class AutoStrategy:
    def __init__(self, n_seqs, avg_len, domain_len=None, gap_ratio=None, redundancy=None, n_species=None):
        self.n_seqs = n_seqs
        self.avg_len = avg_len
        self.domain_len = domain_len
        self.gap_ratio = gap_ratio
        self.redundancy = redundancy
        self.n_species = n_species

    def decide(self, profile="standard"):
        n = self.n_seqs
        if n < 50:
            mafft_mode = "linsi"
            trim_level = "moderate"
            model = "MFP"
            bootstrap = 1000
        elif n < 200:
            mafft_mode = "einsi"
            trim_level = "moderate"
            model = "MFP"
            bootstrap = 1000
        elif n < 1000:
            mafft_mode = "fftns2"
            trim_level = "relaxed"
            model = "MFP+MERGE"
            bootstrap = 100
        else:
            mafft_mode = "fftns2"
            trim_level = "relaxed"
            model = "LG"
            bootstrap = 100

        if profile and profile in PROFILES:
            p = PROFILES[profile]
            bootstrap = p["bootstrap"]
            model = p["model"]

        return {
            "mafft_mode": mafft_mode,
            "trim_level": trim_level,
            "model": model,
            "bootstrap": bootstrap,
            "recommend_representative": n > 2000,
            "recommend_fasttree": n > 5000,
        }

# ── Sequence collection ─────────────────────────────────────────────────────
def collect_sequences(pep_dir, species_filter=None):
    pep_files = sorted(Path(pep_dir).glob("*.pep.fa"))
    if not pep_files:
        print(f"ERROR: No *.pep.fa files found in {pep_dir}")
        sys.exit(1)

    if species_filter:
        filtered = []
        for fpath in pep_files:
            species = fpath.stem.replace(".pep", "")
            if species in species_filter:
                filtered.append(fpath)
        pep_files = filtered
        print(f"  Species filter active: {len(pep_files)} species selected")
    else:
        print(f"  Found {len(pep_files)} PEP files (all species)")

    records = []
    for fpath in pep_files:
        species = fpath.stem.replace(".pep", "")
        with open(fpath) as f:
            cur_id, cur_lines = None, []
            for line in f:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if cur_id is not None and cur_lines:
                        records.append((species, cur_id, "".join(cur_lines)))
                    cur_id = line[1:].split()[0]
                    cur_lines = []
                else:
                    cur_lines.append(line)
            if cur_id is not None and cur_lines:
                records.append((species, cur_id, "".join(cur_lines)))
    return records

def select_longest_per_species(records):
    best = {}
    for species, gid, seq in records:
        slen = len(seq)
        if species not in best or slen > best[species][2]:
            best[species] = (gid, seq, slen)
    result = [(sp, gid, seq) for sp, (gid, seq, _) in best.items()]
    return result

def select_domain_best(records, domtblout_path, logger=None):
    best_per_target = {}
    try:
        with open(domtblout_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split()
                if len(parts) < 23:
                    continue
                target = parts[0]
                evalue = float(parts[12])
                ali_from = int(parts[17])
                ali_to = int(parts[18])
                env_from = int(parts[19])
                env_to = int(parts[20])
                seq_len = int(parts[2])
                dom_bitscore = float(parts[13])
                dom_coverage = (ali_to - ali_from + 1) / max(seq_len, 1)
                is_partial = dom_coverage < 0.5

                key = target
                if key not in best_per_target or evalue < best_per_target[key][0]:
                    best_per_target[key] = (evalue, dom_coverage, dom_bitscore, f">{target}", is_partial)
    except FileNotFoundError:
        if logger:
            logger.warning(f"Domain table not found: {domtblout_path}, falling back to longest")
        return select_longest_per_species(records)

    seq_map = {}
    for species, gid, seq in records:
        seq_map[gid] = (species, seq)

    species_best = {}
    for target, (evalue, coverage, bitscore, header, is_partial) in best_per_target.items():
        if target not in seq_map:
            continue
        species, seq = seq_map[target]
        if species not in species_best:
            species_best[species] = (target, seq, evalue, coverage, is_partial)
        else:
            _, _, cur_e, cur_cov, cur_partial = species_best[species]
            if not is_partial and cur_partial:
                species_best[species] = (target, seq, evalue, coverage, is_partial)
            elif is_partial == cur_partial and evalue < cur_e:
                species_best[species] = (target, seq, evalue, coverage, is_partial)

    result = []
    for species, (gid, seq, evalue, coverage, partial) in species_best.items():
        if partial:
            pass
        result.append((species, gid, seq))

    species_with_best = {s for s, _, _ in result}
    fallback = select_longest_per_species([r for r in records if r[0] not in species_with_best])
    result.extend(fallback)
    return result

def select_representative(records, seq_path, logger=None):
    import tempfile
    if len(records) < 2:
        return records
    try:
        tmp_fa = Path(tempfile.mkstemp(suffix=".fa")[1])
        with open(tmp_fa, "w") as f:
            for species, gid, seq in records:
                f.write(f">{species}|{gid}\n{seq}\n")

        out_cluster = tmp_fa.with_suffix(".cdhit")
        cmd = ["cd-hit", "-i", str(tmp_fa), "-o", str(out_cluster), "-c", "0.9", "-M", "0", "-T", "0"]
        result = run_cmd(cmd, logger=logger, desc="cd-hit clustering", check=False)
        if result.returncode != 0:
            out_cluster2 = tmp_fa.with_suffix(".mmseqs")
            cmd2 = ["mmseqs", "easy-cluster", str(tmp_fa), str(out_cluster2), str(tmp_fa.parent / "tmp")]
            result2 = run_cmd(cmd2, logger=logger, desc="mmseqs clustering", check=False)
            if result2.returncode != 0:
                if logger: logger.warning("Clustering failed, using all sequences")
                return records
            out_file = out_cluster2.with_suffix("_rep_seq.fasta")
        else:
            out_file = out_cluster

        selected = []
        seen_species = set()
        if out_file.exists():
            with open(out_file) as f:
                cur_header = None
                for line in f:
                    if line.startswith(">"):
                        cur_header = line[1:].strip().split()[0]
                    elif cur_header:
                        species = cur_header.split("|")[0] if "|" in cur_header else cur_header
                        if species not in seen_species:
                            selected.append((species, cur_header, line.strip()))
                            seen_species.add(species)

        if not selected:
            if logger: logger.warning("No representatives found, using all sequences")
            return records
        return selected
    except Exception as e:
        if logger: logger.warning(f"Clustering failed: {e}, using all sequences")
        return records
    finally:
        try:
            for f in Path(tmp_fa.parent).glob("tmp*"):
                if f.is_file(): f.unlink()
        except: pass

# ── Filtering ───────────────────────────────────────────────────────────────
def filter_sequences(records, min_len=20, remove_duplicates=True):
    stats = {"total": len(records), "removed_short": 0, "removed_duplicates": 0, "kept": 0}
    filtered = []
    for species, gid, seq in records:
        if len(seq) >= min_len:
            filtered.append((species, gid, seq))
        else:
            stats["removed_short"] += 1

    if remove_duplicates:
        seen_seqs = set()
        deduped = []
        for species, gid, seq in filtered:
            fp = seq[:50]
            if fp not in seen_seqs:
                seen_seqs.add(fp)
                deduped.append((species, gid, seq))
            else:
                stats["removed_duplicates"] += 1
        filtered = deduped

    stats["kept"] = len(filtered)
    return filtered, stats

# ── MAFFT ───────────────────────────────────────────────────────────────────
def build_mafft_cmd(mafft_mode, input_fasta, threads="-1"):
    mode_flags = {
        "auto":     ["--auto", "--maxiterate", "1000"],
        "linsi":    ["--localpair", "--maxiterate", "1000"],
        "einsi":    ["--genafpair", "--maxiterate", "1000"],
        "fftns2":   ["--retree", "2", "--maxiterate", "2"],
        "localpair":["--localpair", "--maxiterate", "1000"],
    }
    flags = mode_flags.get(mafft_mode, mode_flags["auto"])
    cmd = ["mafft"] + flags + ["--thread", str(threads), str(input_fasta)]
    return cmd

def run_mafft(mafft_mode, input_fasta, output_fasta, threads="-1", logger=None, force=False):
    if output_fasta.exists() and output_fasta.stat().st_size > 0 and not force:
        if logger: logger.info(f"Alignment exists, skipping: {output_fasta}")
        n = sum(1 for l in open(output_fasta) if l.startswith(">"))
        return n

    cmd = build_mafft_cmd(mafft_mode, input_fasta, threads)
    with open(output_fasta, "w") as out:
        result = run_cmd(cmd, logger=logger, desc=f"MAFFT ({mafft_mode})", check=False)
        if result.returncode == 0:
            out.write(result.stdout)
        else:
            if logger: logger.warning("MAFFT failed, retrying with --auto")
            cmd2 = build_mafft_cmd("auto", input_fasta, threads)
            with open(output_fasta, "w") as out2:
                result2 = run_cmd(cmd2, logger=logger, desc="MAFFT (auto fallback)")
                out2.write(result2.stdout)

    n = sum(1 for l in open(output_fasta) if l.startswith(">"))
    return n

# ── Trimming ────────────────────────────────────────────────────────────────
def run_trim(trim_program, trim_level, aligned_fasta, output_fasta, logger=None, force=False):
    if output_fasta.exists() and output_fasta.stat().st_size > 0 and not force:
        if logger: logger.info(f"Trimmed alignment exists, skipping: {output_fasta}")
        return

    if trim_program == "clipkit":
        flag = TRIM_PROGRAMS["clipkit"].get(trim_level, "-m smart-gap")
        cmd = ["clipkit", str(aligned_fasta)] + shell_split(flag) + ["-o", str(output_fasta)]
        run_cmd(cmd, logger=logger, desc=f"ClipKIT ({trim_level})")

    elif trim_program == "trimal":
        flag = TRIM_PROGRAMS["trimal"].get(trim_level, "-automated1")
        cmd = ["trimal", "-in", str(aligned_fasta), "-out", str(output_fasta)] + shell_split(flag)
        run_cmd(cmd, logger=logger, desc=f"trimAl ({trim_level})")

    else:
        shutil.copy2(aligned_fasta, output_fasta)
        if logger: logger.info("No trimming (copying alignment as-is)")

# ── IQ-TREE2 ────────────────────────────────────────────────────────────────
def build_iqtree_cmd(seq_fasta, prefix, model="MFP", bootstrap=1000, threads="AUTO"):
    cmd = ["iqtree2", "-s", str(seq_fasta), "--prefix", str(prefix), "-T", str(threads), "--redo"]
    if model == "MFP":
        cmd += ["-m", "MFP"]
    elif model == "MFP+MERGE":
        cmd += ["-m", "MFP+MERGE"]
    elif model == "LG":
        cmd += ["-m", "LG+G4"]
    else:
        cmd += ["-m", model]

    if bootstrap >= 1000:
        cmd += ["-B", str(bootstrap), "-bnni"]
    else:
        cmd += ["-b", str(bootstrap)]
    return cmd

def run_iqtree(seq_fasta, prefix, model="MFP", bootstrap=1000, threads="AUTO", logger=None, force=False):
    treefile = prefix.with_suffix(".treefile")
    if treefile.exists() and treefile.stat().st_size > 0 and not force:
        if logger: logger.info(f"Tree exists, skipping: {treefile}")
        return treefile

    cmd = build_iqtree_cmd(seq_fasta, prefix, model, bootstrap, threads)
    run_cmd(cmd, logger=logger, desc=f"IQ-TREE2 (model={model}, B={bootstrap})")

    if not treefile.exists():
        print(f"  ERROR: IQ-TREE2 did not produce {treefile}")
        sys.exit(1)
    return treefile

# ── QC ──────────────────────────────────────────────────────────────────────
def compute_qc(records, aligned_fasta, trimmed_fasta):
    stats = {}
    stats["total_species"] = len(set(r[0] for r in records))
    stats["total_sequences"] = len(records)
    avg_len = sum(len(r[2]) for r in records) / max(len(records), 1)
    stats["avg_seq_length"] = round(avg_len, 1)
    stats["min_seq_length"] = min(len(r[2]) for r in records)
    stats["max_seq_length"] = max(len(r[2]) for r in records)

    if aligned_fasta and aligned_fasta.exists():
        with open(aligned_fasta) as f:
            lines = f.readlines()
        aln_seqs = []
        cur = ""
        for line in lines:
            if line.startswith(">"):
                if cur: aln_seqs.append(cur)
                cur = ""
            else:
                cur += line.strip()
        if cur: aln_seqs.append(cur)
        if aln_seqs:
            aln_len = len(aln_seqs[0])
            stats["alignment_length"] = aln_len
            gaps = sum(s.count("-") + s.count(".") for s in aln_seqs)
            total = aln_len * len(aln_seqs)
            stats["gap_ratio"] = round(gaps / max(total, 1), 3)

            if aln_len > 0:
                pi_count = 0
                for pos in range(aln_len):
                    chars = set()
                    for s in aln_seqs:
                        c = s[pos].upper()
                        if c not in ("-", ".", "X", "?"):
                            chars.add(c)
                            if len(chars) >= 2:
                                pi_count += 1
                                break
                stats["parsimony_informative_sites"] = pi_count
                stats["pi_ratio"] = round(pi_count / aln_len, 3)

    if trimmed_fasta and trimmed_fasta.exists():
        with open(trimmed_fasta) as f:
            lines = f.readlines()
        trim_seqs = []
        cur = ""
        for line in lines:
            if line.startswith(">"):
                if cur: trim_seqs.append(cur)
                cur = ""
            else:
                cur += line.strip()
        if cur: trim_seqs.append(cur)
        if trim_seqs:
            stats["trimmed_alignment_length"] = len(trim_seqs[0])
            trim_len = len(trim_seqs[0])
            if trim_len > 0:
                pi_count = 0
                for pos in range(trim_len):
                    chars = set()
                    for s in trim_seqs:
                        c = s[pos].upper()
                        if c not in ("-", ".", "X", "?"):
                            chars.add(c)
                            if len(chars) >= 2:
                                pi_count += 1
                                break
                stats["trimmed_pi_sites"] = pi_count
                stats["trimmed_pi_ratio"] = round(pi_count / trim_len, 3)

    return stats

def write_qc_report(stats, out_path):
    lines = [
        "=" * 55,
        "QC REPORT — Domain Tree Pipeline",
        "=" * 55,
        f"Total species:           {stats.get('total_species', '?')}",
        f"Total sequences:         {stats.get('total_sequences', '?')}",
        f"Avg sequence length:     {stats.get('avg_seq_length', '?')} aa",
        f"Min length:              {stats.get('min_seq_length', '?')} aa",
        f"Max length:              {stats.get('max_seq_length', '?')} aa",
        "",
        f"Alignment length:        {stats.get('alignment_length', '?')}",
        f"Gap ratio:               {stats.get('gap_ratio', '?')}",
        f"Parsimony informative:   {stats.get('parsimony_informative_sites', '?')}",
        f"PI ratio:                {stats.get('pi_ratio', '?')}",
        "",
        f"Trimmed alignment len:   {stats.get('trimmed_alignment_length', '?')}",
        f"Trimmed PI sites:        {stats.get('trimmed_pi_sites', '?')}",
        f"Trimmed PI ratio:        {stats.get('trimmed_pi_ratio', '?')}",
        "=" * 55,
    ]
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return lines

# ── Metadata ────────────────────────────────────────────────────────────────
def generate_metadata(pfam_id, args, strategy, stats, runtime, out_path):
    meta = {
        "pfam_id": pfam_id,
        "date": datetime.now().isoformat(),
        "profile": args.profile,
        "strategy": args.strategy,
        "align_mode": strategy.get("mafft_mode", args.align_mode),
        "trim_program": args.trim_program,
        "trim_level": strategy.get("trim_level", ""),
        "iqtree_model": strategy.get("model", "MFP"),
        "bootstrap": strategy.get("bootstrap", 1000),
        "species_count": stats.get("total_species", 0),
        "sequence_count": stats.get("total_sequences", 0),
        "alignment_length": stats.get("alignment_length", 0),
        "trimmed_alignment_length": stats.get("trimmed_alignment_length", 0),
        "parsimony_informative_sites": stats.get("parsimony_informative_sites", 0),
        "runtime_seconds": round(runtime, 1),
        "command": " ".join(sys.argv),
        "cwd": os.getcwd(),
        "python_version": sys.version,
        "clade": args.clade,
        "family": args.family,
        "genus": args.genus,
        "species": args.species,
    }
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(meta, f, indent=2)
    return meta

# ── Methods auto-generation ─────────────────────────────────────────────────
def generate_methods_text(pfam_id, stats, strategy, args, model_found=None):
    profile = args.profile
    mafft_names = {"linsi": "L-INS-i", "einsi": "E-INS-i", "fftns2": "FFT-NS-2", "auto": "auto"}
    mafft_name = mafft_names.get(strategy.get("mafft_mode", args.align_mode), strategy.get("mafft_mode", args.align_mode))

    n_seq = stats.get("total_sequences", "?")
    n_sp = stats.get("total_species", "?")
    aln_len = stats.get("alignment_length", "?")
    trim_len = stats.get("trimmed_alignment_length", "?")
    pi = stats.get("trimmed_pi_sites", stats.get("parsimony_informative_sites", "?"))
    boot = strategy.get("bootstrap", 1000)
    model = model_found or strategy.get("model", "MFP")

    trim_name = args.trim_program
    if trim_name == "clipkit":
        trim_name = "ClipKIT"
    elif trim_name == "trimal":
        trim_name = "trimAl"

    seq_strategy = args.strategy
    if seq_strategy == "longest":
        seq_desc = "the longest protein per species"
    elif seq_strategy == "domain_best":
        seq_desc = "the best-scoring domain-containing protein per species"
    elif seq_strategy == "all":
        seq_desc = "all domain-containing proteins"
    else:
        seq_desc = f"{seq_strategy} sequences"

    methods = textwrap.dedent(f"""\
    ## Methods: Phylogenetic Analysis of {pfam_id}

    Protein sequences containing {pfam_id} domains were retrieved from {n_sp} species
    using HMMER3. For phylogenetic reconstruction, {n_seq} {seq_desc} were selected.
    Multiple sequence alignment was performed using MAFFT ({mafft_name}).
    The alignment spanned {aln_len} positions (trimmed to {trim_len} positions using
    {trim_name}). A total of {pi} parsimony-informative sites were retained.
    Maximum likelihood phylogeny was inferred using IQ-TREE2 under the {model} model
    with {boot} ultrafast bootstrap replicates.
    """).strip()
    return methods

# ── Visualization ───────────────────────────────────────────────────────────
def visualize_tree(treefile, output_stem, logger=None):
    try:
        from ete3 import Tree, TreeStyle, NodeStyle
        t = Tree(str(treefile))
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.branch_vertical_margin = 2

        for leaf in t.iter_leaves():
            sp = leaf.name.split("|")[0] if "|" in leaf.name else leaf.name
            grp = get_group(sp)
            nstyle = NodeStyle()
            leaf.set_style(nstyle)

        for ext in ["pdf", "svg", "png"]:
            out = f"{output_stem}.{ext}"
            try:
                t.render(out, tree_style=ts)
                if logger: logger.info(f"Tree visualization saved: {out}")
                print(f"  Tree visualization: {out}")
            except Exception as e:
                if logger: logger.warning(f"Failed to render {ext}: {e}")

    except ImportError:
        try:
            import matplotlib.pyplot as plt
            try:
                from toytree import tree as toytree_mod
                toyt = toytree_mod(str(treefile), tree_format="newick")
                for ext in ["pdf", "svg", "png"]:
                    canvas = toyt.draw(width=800, height=600)
                    canvas.figure.savefig(f"{output_stem}.{ext}", dpi=150)
                    plt.close()
                    if logger: logger.info(f"Tree visualization saved: {output_stem}.{ext}")
            except ImportError:
                with open(treefile) as f:
                    newick = f.read().strip()
                n_tips = newick.count(",") + 1
                print(f"  Tree has {n_tips} tips. Install ete3 or toytree for visualization.")
                if logger: logger.info("ete3/toytree not installed, skipping visualization")
        except Exception as e:
            if logger: logger.warning(f"Visualization failed: {e}")
            print(f"  WARNING: Tree visualization failed: {e}")

# ── Config file ─────────────────────────────────────────────────────────────
def load_config(config_path):
    try:
        import yaml
        with open(config_path) as f:
            cfg = yaml.safe_load(f)
        if not isinstance(cfg, dict):
            print(f"  WARNING: Config file is empty or invalid")
            return {}
        return cfg
    except ImportError:
        print(f"  WARNING: PyYAML not installed, ignoring --config")
        return {}
    except FileNotFoundError:
        print(f"  ERROR: Config file not found: {config_path}")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"  WARNING: Invalid YAML in config: {e}")
        return {}

# ── Main ────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Industrial-grade domain tree builder",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  %(prog)s PF00168                              # standard profile, auto strategy
  %(prog)s PF00168 --profile accurate            # L-INS-i + strict trim + 2000 BS
  %(prog)s PF00168 --strategy longest --clade Rosids  # one per species, Rosids only
  %(prog)s PF00168 --species ChineseLong,DHL92   # custom species
  %(prog)s PF00168 --align-mode linsi --trim clipkit
  %(prog)s PF00168 --tree-mode ultrafast          # UFBoot2
  %(prog)s PF00168 --config config.yaml           # YAML config
  %(prog)s PF00168 --resume                       # resume from last checkpoint
  %(prog)s PF00168 --list                         # show available species
        """,
    )

    parser.add_argument("pfam_id", type=str, help="Pfam accession (e.g. PF00168)")
    parser.add_argument("--clade", type=str, help="Clade filter (Rosids, Asterids, Monocots...)")
    parser.add_argument("--family", type=str, help="Family filter (Cucurbitaceae, Brassicaceae...)")
    parser.add_argument("--genus", type=str, help="Genus filter (Cucumis, Arabidopsis...)")
    parser.add_argument("--species", type=str, help="Comma-separated species names")
    parser.add_argument("--list", action="store_true", help="List available taxonomic groups")
    parser.add_argument("--profile", type=str, default="standard", choices=list(PROFILES.keys()),
                        help="Execution profile (default: standard)")
    parser.add_argument("--strategy", type=str, default="all",
                        choices=["all", "longest", "canonical", "domain_best", "longest_isoform", "representative"],
                        help="Sequence selection strategy (default: all)")
    parser.add_argument("--align-mode", type=str, default="auto",
                        choices=["auto", "linsi", "einsi", "fftns2", "localpair"],
                        help="MAFFT alignment mode (default: auto)")
    parser.add_argument("--mafft-threads", type=str, default="-1",
                        help="MAFFT threads (default: -1 = all)")
    parser.add_argument("--trim", dest="trim_program", type=str, default="clipkit",
                        choices=["none", "clipkit", "trimal"],
                        help="Alignment trimming program (default: clipkit)")
    parser.add_argument("--trim-level", type=str, default=None,
                        choices=["relaxed", "moderate", "strict", "adaptive"],
                        help="Trim stringency (overrides profile default)")
    parser.add_argument("--tree-mode", type=str, default="standard",
                        choices=["fast", "standard", "ultrafast", "bayesian-ready"],
                        help="Tree building mode (default: standard)")
    parser.add_argument("--iqtree-model", type=str, default=None,
                        help="Override IQ-TREE2 model (e.g. LG+G4, WAG+I+G4)")
    parser.add_argument("--bootstrap", type=int, default=None,
                        help="Override bootstrap replicates")
    parser.add_argument("--iqtree-args", type=str, default="",
                        help="Extra IQ-TREE2 arguments")

    parser.add_argument("-i", "--input-dir", type=str, default=None,
                        help="PEP directory (default: auto-resolve)")
    parser.add_argument("-o", "--output-dir", type=str, default=None,
                        help="Output base directory (default: auto-resolve)")
    parser.add_argument("--delimiter", type=str, default="|",
                        help="FASTA header delimiter (default: '|')")

    parser.add_argument("--force", action="store_true", help="Force redo all steps")
    parser.add_argument("--resume", action="store_true", help="Resume from last checkpoint")
    parser.add_argument("--skip-align", action="store_true", help="Skip MAFFT")
    parser.add_argument("--skip-trim", action="store_true", help="Skip trimming")
    parser.add_argument("--skip-tree", action="store_true", help="Skip IQ-TREE2")
    parser.add_argument("--skip-viz", action="store_true", help="Skip visualization")
    parser.add_argument("--config", type=str, default=None, help="YAML config file")
    parser.add_argument("--verbose", action="store_true", help="Verbose logging")
    parser.add_argument("--debug", action="store_true", help="Debug logging")
    parser.add_argument("--threads", type=str, default="auto",
                        help="Computation threads (default: auto = all CPUs)")
    parser.add_argument("--min-seq-len", type=int, default=20, help="Minimum sequence length (default: 20)")

    args = parser.parse_args()

    if args.list:
        list_species()
        return

    cfg = {}
    if args.config:
        cfg = load_config(args.config)
        for key in ("profile", "strategy", "align_mode", "trim_program", "trim_level",
                    "bootstrap", "iqtree_model", "threads"):
            if key in cfg and getattr(args, key) == parser.get_default(key):
                setattr(args, key, cfg[key])

    pfam_id = args.pfam_id.upper()

    if args.input_dir is not None:
        pep_dir = Path(os.path.abspath(args.input_dir))
    else:
        paths = resolve_pfam_paths(pfam_id)
        pep_dir = paths["pep_dir"]

    if args.output_dir is not None:
        tree_base = Path(os.path.abspath(args.output_dir))
    else:
        if args.input_dir is not None:
            tree_base = Path(os.path.abspath(args.input_dir)).parent / "tree"
        else:
            tree_base = paths["tree_dir"]

    dirs = {
        "input":    tree_base / "input",
        "filtered": tree_base / "filtered",
        "align":    tree_base / "align",
        "trim":     tree_base / "trim",
        "iqtree":   tree_base / "iqtree",
        "logs":     tree_base / "logs",
        "metadata": tree_base / "metadata",
        "figures":  tree_base / "figures",
        "checkpoints": tree_base / "checkpoints",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)

    log_level = logging.DEBUG if args.debug else logging.INFO
    logger = setup_logger(dirs["logs"], level=log_level)
    logger.info("=" * 60)
    logger.info(f"Starting 04_build_domain_tree.py for {pfam_id}")
    logger.info(f"Args: {' '.join(sys.argv)}")
    t_start = time.time()

    if not pep_dir.exists():
        logger.error(f"PEP directory not found: {pep_dir}")
        print(f"  ERROR: PEP directory not found: {pep_dir}")
        print(f"  Run 03 first to generate domain sequences.")
        sys.exit(1)

    selected_species = resolve_species(args)
    species_tag = ""
    if selected_species:
        slist = sorted(selected_species)
        species_tag = "_".join(slist[:3])
        if len(slist) > 3:
            species_tag += f"_etc{len(slist)}"

    ckpt_input = dirs["checkpoints"] / "01_input.done"
    if ckpt_input.exists() and args.resume and not args.force:
        logger.info("Checkpoint 01_input found, resuming...")
        records = []
        input_fasta = dirs["input"] / f"{pfam_id}_{species_tag}_all.fa" if species_tag else dirs["input"] / f"{pfam_id}_all.fa"
        if input_fasta.exists():
            with open(input_fasta) as f:
                cur_species, cur_id, cur_seq = None, None, []
                for line in f:
                    if line.startswith(">"):
                        if cur_id and cur_seq:
                            records.append((cur_species, cur_id, "".join(cur_seq)))
                        h = line[1:].strip().split()[0]
                        parts = h.split("|", 1)
                        cur_species = parts[0]
                        cur_id = parts[1] if len(parts) > 1 else h
                        cur_seq = []
                    else:
                        cur_seq.append(line.strip())
                if cur_id and cur_seq:
                    records.append((cur_species, cur_id, "".join(cur_seq)))
        else:
            logger.warning("Checkpoint exists but input FASTA missing, re-collecting")
            records = collect_sequences(pep_dir, selected_species)
    else:
        records = collect_sequences(pep_dir, selected_species)

    logger.info(f"Collected {len(records)} sequences total")
    print(f"\n  Total sequences: {len(records)}")

    if len(records) == 0:
        logger.error("No sequences found. Exiting.")
        print("  ERROR: No sequences found.")
        sys.exit(1)

    input_fasta = dirs["input"] / f"{pfam_id}_{species_tag}_all.fa" if species_tag else dirs["input"] / f"{pfam_id}_all.fa"
    with open(input_fasta, "w") as f:
        for species, gid, seq in records:
            f.write(f">{species}{args.delimiter}{gid}\n{seq}\n")
    ckpt_input.touch()

    avg_len = sum(len(r[2]) for r in records) / max(len(records), 1)
    n_species = len(set(r[0] for r in records))
    auto = AutoStrategy(n_seqs=len(records), avg_len=avg_len, n_species=n_species)
    auto_decisions = auto.decide(args.profile)
    if args.profile and args.profile in PROFILES:
        p = PROFILES[args.profile]
        if not args.align_mode or args.align_mode == "auto":
            args.align_mode = p["mafft"]
        trim_level = args.trim_level or p["trim"]
    else:
        trim_level = args.trim_level or auto_decisions["trim_level"]
        if args.align_mode == "auto":
            args.align_mode = auto_decisions["mafft_mode"]

    model = args.iqtree_model or auto_decisions["model"]
    bootstrap = args.bootstrap or auto_decisions["bootstrap"]

    logger.info(f"Strategy: align={args.align_mode}, trim={trim_level}, model={model}, B={bootstrap}")

    print(f"\n  Strategy:")
    print(f"    Profile:         {args.profile}")
    print(f"    Strategy:        {args.strategy}")
    print(f"    MAFFT mode:      {args.align_mode}")
    print(f"    Trim:            {args.trim_program} ({trim_level})")
    print(f"    Model:           {model}")
    print(f"    Bootstrap:       {bootstrap}")
    print(f"    Species count:   {n_species}")
    print(f"    Sequence count:  {len(records)}")

    n_orig = len(records)
    if len(records) > 2000 and args.strategy == "all":
        logger.info(f"Large family detected ({len(records)} seqs), auto-enabling representative selection")
        print(f"\n  Large family: {len(records)} sequences > 2000")
        print(f"  Auto-enabling representative-based selection")
        args.strategy = "representative"

    strategy_method = "all"
    ckpt_filter = dirs["checkpoints"] / "02_filter.done"
    if ckpt_filter.exists() and args.resume and not args.force:
        logger.info("Checkpoint 02_filter found, loading filtered sequences")
        filtered_fasta = dirs["filtered"] / f"{pfam_id}_{species_tag}_filtered.fa"
        if filtered_fasta.exists():
            records_filtered = []
            with open(filtered_fasta) as f:
                cur_species, cur_id, cur_seq = None, None, []
                for line in f:
                    if line.startswith(">"):
                        if cur_id and cur_seq:
                            records_filtered.append((cur_species, cur_id, "".join(cur_seq)))
                        h = line[1:].strip().split()[0]
                        parts = h.split("|", 1)
                        cur_species = parts[0]
                        cur_id = parts[1] if len(parts) > 1 else h
                        cur_seq = []
                    else:
                        cur_seq.append(line.strip())
                if cur_id and cur_seq:
                    records_filtered.append((cur_species, cur_id, "".join(cur_seq)))
            records = records_filtered
        else:
            logger.warning("Checkpoint exists but filtered FASTA missing, re-filtering")
    else:
        print(f"\n[Filtering] Strategy: {args.strategy}")
        if args.strategy == "longest":
            records = select_longest_per_species(records)
            strategy_method = "longest per species"
            logger.info(f"After longest selection: {len(records)} sequences")
        elif args.strategy == "domain_best":
            domtblout = paths["hmmer_file"] if not args.input_dir else None
            if domtblout and domtblout.exists():
                records = select_domain_best(records, domtblout, logger)
                strategy_method = "domain best per species"
                logger.info(f"After domain_best selection: {len(records)} sequences")
            else:
                logger.warning(f"Domain table not found, falling back to longest")
                records = select_longest_per_species(records)
                strategy_method = "longest per species (domain_best fallback)"
        elif args.strategy == "representative":
            tmp_fa = dirs["input"] / f"{pfam_id}_for_cluster.fa"
            with open(tmp_fa, "w") as f:
                for species, gid, seq in records:
                    f.write(f">{species}|{gid}\n{seq}\n")
            records = select_representative(records, tmp_fa, logger)
            strategy_method = "representative (clustered)"
            logger.info(f"After representative selection: {len(records)} sequences")
        elif args.strategy == "canonical":
            seen = set()
            canonical = []
            for species, gid, seq in records:
                if species not in seen:
                    canonical.append((species, gid, seq))
                    seen.add(species)
            records = canonical
            strategy_method = "canonical (1st per species)"
            logger.info(f"After canonical selection: {len(records)} sequences")
        elif args.strategy == "longest_isoform":
            gene_best = {}
            for species, gid, seq in records:
                gene_key = f"{species}|{gid.split('.')[0] if '.' in gid else gid}"
                if gene_key not in gene_best or len(seq) > len(gene_best[gene_key][1]):
                    gene_best[gene_key] = (species, gid, seq)
            records = list(gene_best.values())
            strategy_method = "longest isoform"
            logger.info(f"After longest_isoform selection: {len(records)} sequences")
        else:
            strategy_method = "all sequences"

        records, filter_stats = filter_sequences(records, min_len=args.min_seq_len)
        logger.info(f"After filter: {filter_stats['kept']} sequences "
                    f"(removed {filter_stats['removed_short']} short, {filter_stats['removed_duplicates']} duplicates)")

        filtered_fasta = dirs["filtered"] / f"{pfam_id}_{species_tag}_filtered.fa"
        with open(filtered_fasta, "w") as f:
            for species, gid, seq in records:
                f.write(f">{species}{args.delimiter}{gid}\n{seq}\n")
        logger.info(f"Filtered FASTA: {filtered_fasta} ({len(records)} seqs)")
        ckpt_filter.touch()

    if len(records) < 4:
        logger.error(f"Only {len(records)} sequences after filtering. Need >=4 for tree building.")
        print(f"  ERROR: Only {len(records)} sequences. Need >=4 for tree building.")
        sys.exit(1)

    aligned_fasta = dirs["align"] / f"{pfam_id}_{species_tag}_aligned.fa"
    ckpt_align = dirs["checkpoints"] / "03_align.done"

    if ckpt_align.exists() and args.resume and not args.force:
        logger.info("Checkpoint 03_align found, skipping MAFFT")
        print(f"  [Resume] Using existing alignment: {aligned_fasta}")
    elif args.skip_align:
        if aligned_fasta.exists() and aligned_fasta.stat().st_size > 0:
            logger.info("--skip-align: using existing alignment")
        else:
            logger.warning("--skip-align but alignment not found!")
            print(f"  WARNING: --skip-align but {aligned_fasta} not found")
    else:
        print(f"\n[MAFFT] Aligning {len(records)} sequences ({args.align_mode}) ...")
        n_aln = run_mafft(args.align_mode, filtered_fasta, aligned_fasta,
                          threads=args.mafft_threads, logger=logger, force=args.force)
        print(f"  Alignment: {aligned_fasta} ({n_aln} sequences)")
        ckpt_align.touch()

    if not aligned_fasta.exists():
        logger.error("Alignment not available")
        print(f"  ERROR: Alignment not available")
        sys.exit(1)

    trimmed_fasta = dirs["trim"] / f"{pfam_id}_{species_tag}_trimmed.fa"
    ckpt_trim = dirs["checkpoints"] / "04_trim.done"

    if ckpt_trim.exists() and args.resume and not args.force:
        logger.info("Checkpoint 04_trim found, skipping trim")
        print(f"  [Resume] Using existing trimmed alignment: {trimmed_fasta}")
    elif args.skip_trim:
        shutil.copy2(aligned_fasta, trimmed_fasta)
        logger.info("--skip-trim: copying alignment as-is")
    else:
        print(f"\n[Trim] Running {args.trim_program} ({trim_level}) ...")
        run_trim(args.trim_program, trim_level, aligned_fasta, trimmed_fasta, logger=logger, force=args.force)
        ckpt_trim.touch()

    print(f"\n[QC] Computing alignment statistics ...")
    qc_stats = compute_qc(records, aligned_fasta, trimmed_fasta)
    qc_path = dirs["metadata"] / "qc_report.txt"
    qc_lines = write_qc_report(qc_stats, qc_path)
    for line in qc_lines:
        print(f"  {line}")
    logger.info(f"QC report saved to {qc_path}")

    iqtree_prefix = dirs["iqtree"] / f"{pfam_id}_{species_tag}_tree"
    treefile = None
    ckpt_tree = dirs["checkpoints"] / "05_tree.done"

    if ckpt_tree.exists() and args.resume and not args.force:
        logger.info("Checkpoint 05_tree found, IQ-TREE2 already completed")
        treefile = iqtree_prefix.with_suffix(".treefile")
        if not treefile.exists():
            logger.warning("Checkpoint exists but treefile missing")
    elif args.skip_tree:
        logger.info("--skip-tree: IQ-TREE2 skipped")
    else:
        print(f"\n[IQ-TREE2] Running (model={model}, B={bootstrap}) ...")
        treefile = run_iqtree(trimmed_fasta, iqtree_prefix,
                             model=model, bootstrap=bootstrap,
                             threads=args.threads if args.threads != "auto" else "AUTO",
                             logger=logger, force=args.force)
        ckpt_tree.touch()

    model_found = None
    iqtree_report = iqtree_prefix.with_suffix(".iqtree")
    if iqtree_report.exists():
        with open(iqtree_report) as f:
            for line in f:
                m = re.search(r'Model of substitution:\s+(\S+)', line)
                if m:
                    model_found = m.group(1)
                    break

    if treefile and treefile.exists() and not args.skip_viz:
        print(f"\n[Visualization] Generating tree figures ...")
        viz_stem = dirs["figures"] / f"{pfam_id}_{species_tag}_tree"
        visualize_tree(treefile, viz_stem, logger)

    runtime = time.time() - t_start
    qc_stats.update({"total_species": n_species, "total_sequences": len(records)})
    strategy_dict = {
        "mafft_mode": args.align_mode,
        "trim_level": trim_level,
        "model": model_found or model,
        "bootstrap": bootstrap,
    }
    meta = generate_metadata(pfam_id, args, strategy_dict, qc_stats, runtime,
                            dirs["metadata"] / "run_metadata.json")
    logger.info(f"Metadata saved")

    methods_text = generate_methods_text(pfam_id, qc_stats, strategy_dict, args, model_found)
    methods_path = dirs["metadata"] / "methods.txt"
    with open(methods_path, "w") as f:
        f.write(methods_text + "\n")
    logger.info(f"Methods saved to {methods_path}")

    print("\n" + "=" * 60)
    print(f"SUMMARY — {pfam_id} (profile={args.profile})")
    print("=" * 60)
    print(f"  Input sequences:      {qc_stats.get('total_sequences', '?')}")
    print(f"  Alignment:            {aligned_fasta}")
    print(f"  Trimmed:              {trimmed_fasta}")
    if treefile:
        print(f"  Tree file:            {treefile}")
        print(f"  Consensus:            {iqtree_prefix}.contree")
        print(f"  IQ-TREE2 report:      {iqtree_report}")
    print(f"  Model:                {model_found or model}")
    print(f"  Bootstraps:           {bootstrap}")
    print(f"  Runtime:              {runtime:.1f}s")
    print(f"  Log file:             {dirs['logs'] / 'run.log'}")
    print(f"  Metadata:             {dirs['metadata'] / 'run_metadata.json'}")
    print(f"  Methods:              {methods_path}")
    print("=" * 60)
    print(f"Done! (profile={args.profile}, strategy={args.strategy}, align={args.align_mode}, trim={args.trim_program})")
    logger.info(f"Completed in {runtime:.1f}s")

if __name__ == "__main__":
    main()
