# Plant Genome Download Scripts

This directory contains scripts for downloading plant genome data used in this study.

## Data Sources

| Source | Protocol | Tool | Coverage |
|--------|----------|------|----------|
| **Ensembl Plants release 57** | HTTPS | aria2c/wget | ~100 plant species |
| **NCBI Datasets** | HTTPS REST | `datasets` CLI | All published genomes |

> ⚠️ FTP is NOT used: FTP is unreliable for large files (>1GB), especially gymnosperm genomes (up to 22Gb).

## Script Execution Order

### Batch 1: Original 46 species (Ensembl)
- `download_missing_files.sh` — fix missing gff3 for athal/bnap/ptri
- `download_expanded_species.sh` — main Ensembl download (~21 species)
- `download_remaining_parallel.sh` — parallel fill for 4 species

### Batch 2: Additional 8 species
- `download_recommended_species.sh` (first phase) — initial recommended species download

### Batch 3: All remaining (~55 species) [CURRENTLY RUNNING]
- `download_remaining_robust.sh` — **recommended script** with aria2c HTTPS + checkpoint + retry

## Recommended Script

**`download_remaining_robust.sh`** is the most robust download script:
- aria2c (multi-threaded HTTPS, supports resume `-c`, retry `--max-tries=10`)
- NCBI datasets CLI for non-Ensembl species
- tmux session for persistence (survives SSH disconnect)
- Checkpoint system (`/tmp/download_robust/checkpoint.txt`)
- Auto-skip already-downloaded species

## Three Phases

1. **Phase 1**: Ensembl remaining ~35 species via aria2c HTTPS
2. **Phase 2**: NCBI small/medium genomes ~20 species via datasets CLI
3. **Phase 3**: Pinus taeda (~22Gb) — special handling (pep+gff only)

## Output Format

All species stored in unified format:
```
/data5/qiulei/PPI/data/<species_name>/
├── <species_name>.pep.fa.gz      # Protein sequences
├── <species_name>.cds.fa.gz      # CDS sequences
├── <species_name>.gff3.gz        # Gene annotations
```

## Monitoring

```bash
# View tmux session
# tmux attach -t download_robust

# View checkpoint
cat /tmp/download_robust/checkpoint.txt
```

## Notes

1. Ensembl SSL certificate issue: `ftp.ensemblgenomes.org` has hostname mismatch — all wget/aria2c need `--no-check-certificate` / `--check-certificate=false`
2. Pinus taeda (~22Gb): only pep+gff downloaded, skip CDS
3. All downloads use HTTPS, never FTP
4. The scripts are designed to be re-runnable (checkpoint + existing file detection)
