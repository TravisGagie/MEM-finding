# Block length as genomes are added (2026-11-05, cloud container)

`experiments/block_scaling.sh` on a synthetic pangenome from `scripts/pangenome_sim.py --haplotypes 1024 --queries 20
--length 300000` (seed 1; parameters and population statistics in `pangenome_stats.txt`). Collections are the reference
plus the first N haplotypes, N = 16 … 1024 (n = N + 1 sequences). The patterns are the 20 held-out haplotypes (6.0 million
characters), MEMs of length at least 1000, k = 51, s = 12.

- **Correctness:** at every n, the MEMs and intervals were byte-identical under both rules, with and without block
  steps, and identical to plain search (`CHECK_PLAIN=1`).
- **Files:** `table.txt` (the script's table), `blocks.tsv` (the same, machine-readable), `logs/` (index and search logs).
