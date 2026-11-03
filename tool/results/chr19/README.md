# Chromosome 19 (Section 5), new tool, desktop

Desktop (i7-8700, 14 GiB), 2026-11-02/03. The tool was at commits 5a3d161 (blocks runs) and c09732c (STPD).
Data: `~/shortcuts/chr19_scaling/hprc_chr19` (CHM13 plus the HPRC haplotypes; the first 33 or 65 lines are 16 or
32 individuals) and its `queries.txt`: 4 patterns, 122,434,624 characters once split at runs of N.
Each step ran on its own, in a cgroup capped at 11 GB with no swap. Queries ran single-threaded, pinned to CPU 9
(as the old `chr19/twolevel/run_final.sh` runs were). Times are one run each; the old runs used the median of 5.

## Two-level search with blocks (`run_scripts/run_n.sh`, as in `tool/NOTES.md`)

| | 16 individuals (33 seq.) | 32 individuals (65 seq.) |
|---|---|---|
| text S (both strands) | 4,065,485,530 | 8,004,217,316 |
| `prep` | 18 s, 0.13 GB | 38 s, 0.13 GB |
| Big-BWT | 345 s, 2.6 GB | 599 s, 3.5 GB |
| `index --preset chr19` | 214 s, 4.6 GB | 368 s, 8.8 GB (9.04 GB as the tool reports it) |
| RLBWT / RLCSA / two-level (`.lvl`) | 191 / 579 / 209 MB | 218 / 1,193 / 250 MB |
| parse / blocks in `.lvl` | 166.4 / 42.2 MB | 194.3 / 55.6 MB |
| L = 1000: MEMs | 22,711 | 18,479 |
| two-level with blocks, no jumps, plain (RLCSA) | 5.66 s, 8.46 s, 40.4 s | 5.54 s, 8.95 s, 50.7 s |
| forward-backward (two-level with blocks) | 5.65 s | 5.36 s |
| L = 100: MEMs | 38,984 | 30,809 |
| two-level with blocks, no jumps, plain | 5.84 s, 8.57 s, 39.2 s | 5.65 s, 9.03 s, 49.8 s |
| query peak (two-level / plain) | 1.36 / 1.03 GB | 2.02 / 1.65 GB |

Every two-level and forward-backward output is byte-identical to plain search.

The MEMs are on both strands (one text), so the counts differ from the old two-index code, which had 38,890 at
L = 100 for 16 individuals. For comparison, the old code (`~/twolevel_runs/twolevel/final/summary.txt`,
16 individuals, L = 1000) took 290.5 s for plain search (RLBWT), 12.3 s for two-level and 8.05 s for two-level
with blocks.

64 individuals were not run: the 32-individual index already peaks at 9.0 GB, and the 64-individual one needs
about 40 GB of disk.

## STPD: katka against the published STPD-index (16 individuals)

Their code: github.com/regindex/STPD-index at 8a432a4. It was built with CMake 3.31.8 and two configure-time
workarounds for this machine; none of their files was changed:
- `-DCMAKE_PROJECT_INCLUDE=run_scripts/linker_lang_shim.cmake`, which gives their header-only OBJECT libraries
  LINKER_LANGUAGE CXX, without which CMake refuses them;
- `-Wno-template-body`, because GCC 15 rejects an uninstantiated template in sdsl's `louds_tree.hpp`;
- plus the divsufsort include directory, as `NOTES.md` says.

The comparison is `tool/experiments/stpd_index_chr19.sh`.
- **The script:** it needs an absolute WORK directory. With a relative one, as in the `NOTES.md` example, it
  fails at once: it `cd`s into WORK, then writes `WORK/summary.txt`.
- **Their forward-backward search (`stpd-index-mems --fb`):** it was stopped after 54 minutes at 100% CPU with
  no output. Their BML takes 45 s on the same index.
- **The remaining runs:** `run_scripts/stpd33_rest.sh` runs them with the script's commands, plus its comparison
  code verbatim (`stpd33_compare.py`), without their forward-backward search.

`katka index --stpd c33` took 400 s with a peak of 8.9 GB. It added 44,361,753 sampled positions (0.648 r):
177 MB on top of the text, plus 350 MB for φ⁻¹.

Their index was built from katka's files (`stpd33/build_st.log`):

| step | time | peak |
|---|---|---|
| `stpd-prep` | 31 s | 0.14 GB |
| Big-BWT | 359 s | 2.8 GB |
| `index --stpd` | 431 s | 8.9 GB |
| `stpd-export` | 2,590 s | 6.4 GB |
| their `stpd-build` | 507 s | 4.9 GB |

Their index is 584.3 MiB (612.7 MB) and includes their RLZ text oracle.

### Queries (L = 1000 / L = 100), single-threaded

| configuration | search time (s) | wall incl. loading (s) | query peak | files needed (MB) |
|---|---|---|---|---|
| katka two-level with blocks (`--csa --twolevel`) | 5.96 / 6.13 | 6.66 / 6.70 | 1.36 GB | 979: RLBWT 191 + RLCSA 579 + `.lvl` 209 |
| katka STPD, MEMs only, no intervals (`--stpd --no-intervals`) | 3.35 / 3.64 | 9.48 / 9.42 | 4.95 GB | 4,784: text 4,065 + `.stpd` 528 + RLBWT 191 |
| katka STPD with intervals by plain search (`--stpd`) | 132.4 / 136.1 | 138.2 / 142.0 | 4.95 GB | 4,784 |
| katka STPD, intervals by two-level search (`--csa --stpd --twolevel`) | 6.47 / 7.26 | 13.6 / 14.3 | 5.84 GB | 5,571: the above + RLCSA + `.lvl` |
| katka STPD with locating (`--stpd --locate`) | 132.0 / 135.8 | 138.0 / 141.8 | 4.95 GB | 4,784 |
| theirs, BML (MEMs only) | 44.6 / 54.8 | — / 72.4 | 0.96 GB (L = 100) | 613: `st.stpdi`, RLZ oracle included |
| theirs, BML with locating | 42.4 + 11.4 / 54.8 + 15.4 | — / 86.1 | 0.96 GB (L = 100) | 613 |

- **Locating:** 72,876 occurrences at L = 1000 and 384,341 at L = 100.
- **Peaks:** the script didn't record peak memory for their L = 1000 runs. The L = 100 runs use the same index.
- **Their driver's check rejected** 70,241 candidates at L = 1000 (of 1,285,721 jumps) and 101,235 at L = 100
  (of 1,834,289).
- **Text:** katka's STPD keeps the whole text in memory (4.07 GB), which dominates its space. Their RLZ oracle
  is about a seventh of that.
- **Locating time:** katka's `--stpd --locate` time is mostly the plain-search intervals (438 M backward steps;
  `--stpd` alone takes as long). `--stpd --locate --twolevel` was not in the script.

### Comparison (reference: katka's plain BML)

At L = 1000 and at L = 100, every configuration, katka's and theirs, gives exactly the same MEMs: 22,711 and
38,984. Theirs flags 4 and 7 MEMs, respectively, whose first occurrence crosses a junction of their separator-free
text. These MEMs have the same coordinates as katka's. Leaving them out, as asked, the comparison is 22,707 and
38,977 MEMs, all identical. Both locating runs give occurrence counts equal to ep − sp for every MEM.

32 individuals were not tried for STPD. At 16 individuals, `index --stpd` already peaks at 8.9 GB and
`stpd-export` at 6.4 GB, both with a 4 GB text; with the 8 GB text they would exceed the 11 GB cap.
