# Documentation Index

| Document | Description |
|----------|-------------|
| [installation.md](installation.md) | System requirements, external tool installation, Python/R dependencies |
| [usage.md](usage.md) | Full pipeline usage: parameter tables, examples, output structure for all 13 scripts |

## Quick Links

- **Core domain pipeline** (4 scripts): Pfam HMM → hmmsearch → domain parsing → sequence extraction → phylogenetic tree
- **OrthoFinder synteny pipeline** (7 Python + 2 R scripts): Species tree → ortholog tree → synteny analysis → visualization

## Pipeline Overview

```
                        ┌──────────────────────────────────────┐
                        │         Core Domain Pipeline         │
                        │                                      │
                        │  01_download_hmm_search_domain.py    │
                        │         ↓                            │
                        │  02_parse_domain_table.py            │
                        │         ↓                            │
                        │  03_extract_domain_sequences.py      │
                        │         ↓                            │
                        │  04_build_phylogenetic_tree.py       │
                        └──────────────────────────────────────┘

                        ┌──────────────────────────────────────┐
                        │     OrthoFinder Synteny Pipeline     │
                        │                                      │
                        │  N1 → N2 ──┬── R1 (classic tree)     │
                        │            ├── R2 (tree + heatmap)   │
                        │            │                          │
                        │            ├── N3 → N4 ──┬── N6      │
                        │            │              └── N7      │
                        │            └── N5 ───────→ N6, N7    │
                        └──────────────────────────────────────┘
```

## Related Documents

- `config/species_config.py` — Central species groups, colors, and path resolution
- `CHANGELOG.md` — Version history and updates
