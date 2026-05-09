# Quality Assurance Report

## Summary

| Check | Status | Notes |
|-------|--------|-------|
| **Scripts documented** | 13/13 (100%) | 4 core (01–04) + 7 N-step + 2 R scripts |
| **Scripts with parameter tables** | 13/13 (100%) | Each script has complete parameter table |
| **Scripts with minimal example** | 13/13 (100%) | Each script has at least one runnable example |
| **Scripts with full example** | 13/13 (100%) | Each script has real-parameter example |
| **Scripts with common pitfalls** | 8/13 (62%) | Steps 01, 03, 04, N1–N7 covered |
| **Scripts with output structure** | 13/13 (100%) | Inputs/outputs specified for all |
| **README reflects actual repo** | ✅ | No fictional content; repository scope noted |
| **Hardcoded paths flagged** | ✅ | docs/usage.md notes defaults are hardcoded |
| **Unresolved TODO/Coming soon** | 0 | Removed all placeholders |
| **Broken commands in docs** | ✅ | All commands verified against actual script signatures |
| **Directory naming issues noted** | ✅ | OrthoFinder trailing newline documented |

## Detailed Findings

### A. Repository Structure Issues

| Issue | Severity | Recommendation |
|-------|----------|----------------|
| Repo has scripts/ empty (no actual scripts) | **HIGH** | Only `config/species_config.py` exists; production scripts are at `/data5/qiulei/script/0506/`. README now transparently notes this. |
| `scripts/_version.py` referenced in old docs doesn't exist | **LOW** | Removed from installation.md verification section |
| `config/default.yaml` listed in old README doesn't exist | **LOW** | Removed from listing |
| `tests/`, `workflows/`, `data/` directories don't exist | **LOW** | Removed from repo structure; these were fictional |

### B. Script Naming & Config Inconsistencies

| Issue | Severity | Recommendation |
|-------|----------|----------------|
| `config/stepN1_species_tree.yaml` vs `stepN1_build_species_tree.yaml` — duplicate | **MEDIUM** | Remove one; `build_species_tree` is the canonical script name |
| `stepN3_last_alignment.yaml` vs script `stepN3_run_last_alignment.py` — missing `run_` | **LOW** | Rename YAML to `stepN3_run_last_alignment.yaml` |
| N1 YAML has orphan keys `msa_program` and `tree_method` never read by script | **LOW** | Either remove or implement in script |
| CLI param `--tree` vs YAML key `tree_file` (R1, R2) | **LOW** | Inconsistent but not broken since YAML mapped manually |

### C. Path Hardcoding Issues

Hardcoded `/data5/qiulei/pfam_pipeline/` paths found in:

| File | Count | Example |
|------|-------|---------|
| `stepN1_species_tree.yaml` | 2 | `proteome_dir`, `output_dir` |
| `stepN1_build_species_tree.yaml` | 2 | Same paths |
| `stepN2_ortholog_tree.yaml` | 1 | `orthofinder_dir` |
| `stepN3_last_alignment.yaml` | 2 | `cds_dir` + implicit pipeline-root |
| `stepN4_convert_anchors.yaml` | 0 | (species only) |
| `stepN5_generate_bed.yaml` | 1 | `data_dir` |
| `stepN6_highlight_chain.yaml` | 0 | (species only) |
| `stepN7_microsynteny.yaml` | 1 | `data_dir` |
| `stepR1_tree.yaml` | 2 | `tree_file`, `output_dir` |
| `stepR2_heatmap.yaml` | 2 | `tree_file`, `output_dir` |
| N2–N7 script defaults | 6 | `--pipeline-root`, `--orthofinder-dir`, etc. |
| R1–R2 script defaults | 2 | `--output-dir` fallback |

**Total: ~21 hardcoded paths** — Migration required for any new deployment.

### D. Documentation Completeness

| Document | Word Count | Scripts Covered | Status |
|----------|-----------|-----------------|--------|
| `README.md` | ~600 | Pipeline overview | ✅ Rewritten |
| `docs/README.md` | ~100 | Index + diagram | ✅ New |
| `docs/usage.md` | ~3200 | 13 scripts | ✅ Rewritten |
| `docs/installation.md` | ~500 | N/A | ✅ Updated |

### E. Reproducibility

| Feature | Status |
|---------|--------|
| Auto-generated methods.txt for publications | ✅ (step 04 only) |
| Run metadata JSON | ✅ (step 04 only) |
| YAML config support | ✅ (N1–N7, R1–R2) |
| Checkpoint/resume | ✅ (step 04 only) |
| Seed specification for stochastic steps | ❌ (not implemented) |
| Container/Docker support | ❌ (not implemented) |

## Recommendations

1. **Integrate production scripts** into this repository (`scripts/` directory) so the repo is self-contained
2. **Remove duplicate YAML** `stepN1_species_tree.yaml` — keep only `stepN1_build_species_tree.yaml`
3. **Resolve orphan keys** in N1 YAML (`msa_program`, `tree_method`)
4. **Rename** `stepN3_last_alignment.yaml` → `stepN3_run_last_alignment.yaml`
5. **Parameterize all hardcoded paths** — either via env vars or a centralized config
6. **Add `_version.py`** or a version module if packaging is intended
7. **Remove fictional directories** (`tests/`, `workflows/`, `data/`, `config/default.yaml`) from any future repo descriptions

---

*Generated: 2026-05-09*
