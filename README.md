# MEM-finding

Code for papers on finding maximal exact matches (MEMs) with indexes based on the run-length compressed Burrows–Wheeler transform (RLBWT).

| Paper | Code | Reproduction |
|---|---|---|
| *Reportable Shortcuts with RLBWTs*, arXiv:2609.06596 (v3) | `chr19/shortcuts/` (including its chr19 MEM experiment) | `chr19/shortcuts/README.md`, "Reproducing the experiments" |
| *More RLBWT Tricks*, arXiv:2609.40299 (v2): two-level indexing at closed syncmers with blocks of continuations, the STPD comparison, two-level locating | `tool/` (branch `tool`); `chr19/twolevel/` is the earlier implementation, on forward strands only | `tool/CHR19.md` |
| *KATKA* (preliminary taxonomic classification), arXiv:2610.10500 | `katka/` (the original code) and `katka-readable/`; superseded by `tool/`'s `zymo` and `silva` presets | `katka/README.md`, `katka-readable/README.md` |
| In progress: the KeBaB revision experiments | `tool/results/kebab_revision/`, `tool/results/zymo_kebab/` | their `README.md` files |
| In progress: classifying SILVA 16S reads, against Cliffy | `tool/` (`silva` preset) | `tool/results/silva/README.md`; `tool/NOTES.md`, "SILVA queue" |

The codebases have separate build systems. `tool/` builds with `make -C tool` (see `tool/README.md`).

## License

GPL-3.0 (see `LICENSE` and `katka/LICENSE`).
