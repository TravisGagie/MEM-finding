# Zymo results (Section 3), new tool, desktop

Machine: Travis's desktop (Intel i7-8700, 6 cores / 12 threads, 14 GiB RAM, NVMe disk, Linux). Run 2026-11-02.
The old numbers in the paper were measured on the laptop; KATKA was rerun here for the comparison.
Each step ran on its own, with nothing else running, in a cgroup capped at 11 GB with no swap.
Timing: `/usr/bin/time` (wall, user CPU, peak RSS). Searches use 12 threads.
The scripts that produced every number are in `run_scripts/`, and their logs are in `logs/`.

Data: the 7-species reference of SPUMONI 2 (14,280 genomes) and its 581,802 ONT reads (`all.fq`, `all.truth`),
plus the 720-base prefix test (581,802 bacterial and 10,245 yeast prefixes). The index is homopolymer-compressed.

## Index build

The BWT is KATKA's laptop `zy.bwt` (one string per species, Big-BWT and pfp-merge). It was checked against a fresh
`prep`: the same length, 7 datasets, the same table and the same character counts. The tag runs are KATKA's laptop
`zy.sp.tagruns`, converted by `scripts/convert_katka_tagruns.py` and checked against the BWT by
`tests/check_tagruns.cpp`. Their construction took KATKA 7.2 h on the laptop; that is KATKA's figure, not
measured here.

| step | wall | peak RSS | output | size on disk |
|---|---|---|---|---|
| `prep --preset zymo` (`-d -s -U`) | 311 s | 23.6 MB | `zy.S` (39,277,433,146 chars; deleted after the build), `zy.tbl` | 39.28 GB; 1.4 MB |
| `index --csa` | 263 s | 6.59 GB | `zy.rlbwt` (n = 39,277,433,153, r = 338,097,187), `zy.csa` | 1.024 GB; 5.563 GB |
| tag runs: conversion | 17 s | 4.70 GB | `zy.sp.tagruns` (448,723,822 runs) | 5.385 GB (build cache only) |
| tag runs: check against the BWT | 119 s | 7.01 GB | 35,000,000 rows from the 7 terminators, in 31.5 M runs: 0 mismatches | — |
| `index --tags --by-tag --no-rmq` | 152 s | 6.97 GB | `zy.tag` (from the cached runs: 11.7 s; the rest is the RLBWT, 140 s) | 0.687 GB |
| `index --kebab 20` | 2,569 s | 6.06 GB | `zy.k20.kbf` (filter alone: 2,428 s), byte-identical to KATKA's | 0.135 GB |

The index used by the searches is `zy.rlbwt`, `zy.csa`, `zy.tag` and `zy.k20.kbf`: 7.41 GB in total, and 1.85 GB
without the RLCSA. The BWT file (`zy.bwt`, 39.28 GB) is an input only.

The steps were run separately because `index --preset zymo` in one process keeps the RLCSA in memory during the
tag step (see `NOTES.md`).

## Classification, all 581,802 reads

Accuracy uses KATKA's scripts and the laptop's protocol: `katka-readable/scripts/credit_schemes.py` and
`run_scripts/score_spu2.py` for L = 30. For the longest MEMs it uses `katka-readable/scripts/maxmem_roc.py`,
which keeps only each read's longest MEMs, by homopolymer-compressed length. Here it runs as
`scripts/maxmem_roc_stream.py`, identical output with the reads streamed.

| run | options | wall | user CPU | peak RSS | MEMs | accuracy (prop / even) | paper (laptop) |
|---|---|---|---|---|---|---|---|
| threshold L = 30 | `-U -L 30 -D 14 -c --csa` | 566 s | 6,372 s | 7.21 GB | 89,640,946; 2,212 reads without | 99.433% / 99.538% | 99.43%, ~780 s |
| longest MEMs, no KeBaB | `-U -X -L 20 -D 14 -c --csa` | 126 s | 1,344 s | 7.20 GB | 3,993,026 | 99.61% / 99.27% | 211 s |
| longest MEMs, KeBaB k = 20 | `--preset zymo` (adds `-K 20`) | 90 s | 919 s | 7.34 GB | 814,201 | 99.61% / 99.27% | 99.61%, 118 s |

- With and without KeBaB, every read has the same longest MEMs: the same spans and `{tag:count}`
  (`scripts/same_longest_mems.py`: 0 of 581,802 reads differ). The no-KeBaB run lists more MEMs because adaptive
  BML reports every MEM it finds as L rises, while KeBaB's longest-first search stops sooner.
- With L = 30 the other credit schemes give: one 99.438%, cliffs 99.576%, range 99.538% (as in the paper).
- On the laptop, KATKA's no-KeBaB run used `-L 1` instead of `-L 20`.

## Prefix test (SPUMONI 2's Table 2), 720 bases

| run | options | wall | user CPU | peak RSS | sensitivity | specificity | paper |
|---|---|---|---|---|---|---|---|
| L = 30, θ = 14, RLBWT and tags (no RLCSA) | `-U -L 30 -D 14 -c` | 90 s | 1,032 s | 1.68 GB | 96.88% | 97.90% | 96.88% / 97.90%, 109 s |

A prefix is accepted if it has at least one MEM of at least L compressed characters after the −D 14 filter.

## Regression against KATKA, same machine

KATKA is katka-readable's `katka-classify` on KATKA's own laptop index files (`zy.rz`, `zy.csa`, `zy.s1.tag`, `zy.k20.kbf`).

| run | KATKA wall / user / peak | new wall / user / peak | output |
|---|---|---|---|
| `-U -L 30 -D 14 -c`, RLCSA | 571 s / 6,384 s / 7.95 GB | 566 s / 6,372 s / 7.21 GB | byte-identical (1,837,518,643 bytes) |
| longest MEMs, KeBaB, `-D 14`, RLCSA | 89 s / 893 s / 8.04 GB | 90 s / 919 s / 7.34 GB | byte-identical (28,090,295 bytes) |

Wall times include loading the index (about 9 s for the new tool). The read file (23 GB) does not fit in the page cache,
so the order of the runs was alternated.

## Grammar-compressed tags (2026-11-05)

`grammar_run.sh`: the per-species tag runs B_λ verified (`zy.sp.tagruns`, 448,723,822 runs), `bigrepair -i -r -m 0`, then `index --tags --by-tag --no-rmq --bigrepair`. The built-in RePair wasn't run: it needs about 25 bytes per run, about 11.4 GB here, over the 11 GB cap.

| | run-length tags | grammar tags |
|---|---:|---:|
| tag array | 687,420,273 bytes | 750,537,217 bytes |
| tag sequence | 168.3 MB (3 bits per run) | 231.4 MB (46,052,325 rules, top level 3,539,044 symbols) |
| run starts, unchanged | about 519 MB | about 519 MB |
| time per read (one run, 100,000 reads, 12 threads) | 200.8 µs | 202.9 µs |

- **Ratio 1.09: the grammar makes the tag array 9% larger.** The prediction of little gain holds, and goes further: there is a loss. On SILVA bigrepair gives 0.61 (836.4 to 509.5 MB) and the built-in RePair 0.54 (449.7 MB).
- **Why:** two limits apply. Tag runs per BWT run are 1.33, against SILVA's 5.3. And only the tag sequence can shrink, which is 24% of Zymo's array, because 7 tags need 3 bits; on SILVA, with 9,118 tags at 14 bits, it is 75%. bigrepair did compress the sequence: 1,795 MB of 32-bit integers became 383 MB, or 173 MB as it reports its own compressed size. But 3 bits per run is already below the 4.1 bits per run the grammar costs here.
- **The built-in RePair** gave a 12% smaller array than bigrepair on SILVA. To gain here its grammar would have to fall below 168 MB from 231 MB, a 27% cut. That is not measured.
- **Build:** `tagruns_to_int32.py` 10 s, 0.11 GB; bigrepair 197 s, peak 5.56 GB; `index` 2,057 s, peak 8.18 GB. Of the `index` time, 2,013 s is rebuilding the RLBWT from `zy.bwt` on the external drive, and the tag step takes 44 s. All within the cap.
- **Check:** on the first 100,000 reads of `all.fq`, `find -U -X -L 20 -D 14 -c --csa` gives byte-identical outputs, 14,099,474 bytes. `zy.rlbwt`, rewritten by `index`, kept its SHA-256.
