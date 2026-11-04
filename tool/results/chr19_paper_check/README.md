# CHR19.md's pipeline on real data, at small scale (cloud container, 2026-11-04)

A check that `experiments/chr19_paper.sh all` runs end to end on real HPRC data, not a reproduction of the
paper's numbers.
- **The collection:** CHM13 and both haplotypes of HG00097 (3 sequences, `manifest.tsv`). The queries are the
  paper's two haplotypes. Everything was downloaded with `experiments/hprc_chr19.py --manifest manifest.tsv`
  (23 s), with every length and N count checked.
- **The text:** S is 367.8 million characters (both strands).
- **The BWT:** built with `katka bwt`.
- **The machine:** 4 cores, 15 GB, no cgroup.

Every configuration gave the same MEMs as plain search: 33,452 at L = 1000 and 97,959 at L = 100.
- **Locating:** the occurrence counts equal ep − sp (48,942 and 465,458 occurrences).
- **Positions:** two-level locating and STPD over the parse wrote byte-identical position lists for every MEM
  (`locating.md`).

The times are this container's, with only one individual. That collection is far less repetitive than the
paper's: 6.7 million distinct phrases out of 19.9 million, against 8.5 of 220.5 with 16 individuals. Its MEMs
are shorter and more numerous. So the times are not comparable to the paper's tables; they are here only as
a record. The whole run took 24 minutes.

Files: `table1.md`, `table2.md`, `locating.md` (the script's outputs); `logs/` (the build steps).
