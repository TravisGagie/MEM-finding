# STPD over the parse on chr19, 16 individuals (c33)

Desktop (i7-8700, 14 GiB), 2026-11-03. The tool was at 4700b46: `dc63572`'s packed dictionary and `aba178a`'s
RLZAP oracle, one binary for every run.
- **Runs:** `experiments/stpd_parse_chr19.sh`, wrapped in `run33.sh`. The wrapper adds `index --parse-oracle`
  first, for the dictionary's size, and katka's STPD over the text.
- **Conditions:** queries single-threaded, pinned to CPU 9; each step in an 11 GB cgroup with no swap.
- **Logs:** `summary33.txt`, `run33.log`, `logs/`.
- **The STPD sample:** rebuilt through the parse oracle, without the text. It is byte-identical to the one built
  from the text the night before.

## Building the components

| step | wall | peak | output |
|---|---|---|---|
| `index --parse-oracle` | 104.7 s | 1.52 GB | `c33.ptxt` 1,036.4 MB, see below |
| `index --stpd --parse-oracle` (character-level sample, text read through the oracle) | 664.2 s | 5.97 GB | `c33.stpd` 527.6 MB (sample 177.4 + φ⁻¹ 350.2) |
| `index --stpd-parse` | 57.2 s | 3.39 GB | `c33.pstpd` 116.9 MB |
| `index --rlzap` (J, D, M = 4, 4, 2) | 39.9 s | 2.78 GB | `c33.rptxt` 316.6 MB |

**Parse oracle with the 2-bit dictionary:**
- **Dictionary:** **151.2 MB**, down from 630.5: 598,602,502 bases at 2 bits. Only 89 phrases contain other
  characters, with 89 runs of X. The dictionary's offsets take 31.9 MB.
- **Rest of the oracle:** parse 661.4 MB, phrase starts 191.8 MB.
- **Total:** **1,036.4 MB**, against 1,483.7 MB before: 2.04 bits per character.

**The RLZAP parse:** 12,352,342 factors (17.31 phrases each), 18,644,761 literal phrases, 1,837,448 pointer changes
and 1,611,960 explicit pointers. These are the counts of `results/chr19_rlz/`. Its parse part takes 133.4 MB,
against 853.2 MB for the plain parse and its starts:
- reference ids 20.1 and starts 10.0;
- factor starts 10.7 (parse) and 17.2 (text);
- literal counts 5.3;
- literals 55.9;
- pointers 14.2.

## Queries (search time, wall including loading, peak; L = 1000 / L = 100)

| configuration | search (s) | wall (s) | peak | files read (MB) |
|---|---|---|---|---|
| plain search (RLCSA) | 40.66 / 39.18 | 41.7 / 39.6 | 1.05 GB | 770: RLBWT 191 + RLCSA 579 |
| **two-level with blocks** | **5.62 / 5.72** | 6.2 / 6.3 | 1.39 GB | 979: + `.lvl` 209 |
| STPD over the parse, plain parse | 5.69 / 6.40 | 7.7 / 7.9 | 2.98 GB | 1,982: see below |
| STPD over the parse, plain parse, `--locate` | 5.73 / 9.26 | 7.2 / 10.7 | 2.98 GB | 2,378: + parse φ⁻¹ 45.5 + character φ⁻¹ 350.2 |
| STPD over the parse, RLZAP parse | 8.65 / 10.17 | 9.7 / 11.2 | 2.26 GB | 1,263: the plain parse and starts (853.2) replaced by 133.4 |
| STPD over the parse, RLZAP parse, `--locate` | 8.98 / 17.96 | 10.0 / 19.0 | 2.26 GB | 1,659 |
| character-level STPD through the parse oracle (`--stpd --parse-oracle`) | 16.19 / 19.25 | 17.2 / 20.3 | 2.04 GB | 1,755: `.ptxt` 1,036 + `.stpd` 528 + RLBWT 191 |
| katka's STPD over the text (`--stpd --no-intervals`) | 3.25 / 3.57 | 8.5 / 6.3 | 4.95 GB | 4,784: text 4,065 + `.stpd` 528 + RLBWT 191 |
| katka's STPD over the text, `--locate` | 3.23 / 3.76 | 6.1 / 6.5 | 4.95 GB | 4,784 |

How the table was read:
- **"Search"** is the tool's own timing, which excludes loading the index.
- **Correctness:** every configuration gives the same MEMs as plain search, 22,711 at L = 1000 and 38,984 at
  L = 100. Every `--locate` run gives occurrence counts equal to ep − sp. The RLZAP runs give byte-identical
  output to the plain parse's.

STPD over the parse (plain parse) reads these components, in MB:
- **the parse:** plain parse 661.4, starts 191.8, dictionary 151.2 plus its offsets 31.9;
- **the samples:** parse sample 33.5 (9,571,667 of r_Z = 10,332,018), entry table 29.8, C 8.1;
- **the character-level sample:** 177.4;
- **the right end:** RLCSA 579.3, B 12.8, dictionary hash 105.2;
- **for `--locate` only:** parse φ⁻¹ 45.5 and character φ⁻¹ 350.2.

Its step counts at L = 1000:
- right-end character steps: 4.24 M;
- middle phrase comparisons: 13.9 M;
- parse-sample binary searches: 157,005, with 3.81 M comparisons;
- left-end character steps: 1.11 M;
- character-level STPD jumps: 65,987.

Plain search takes 259.8 M character steps.

## Reading

- **STPD over the parse is not faster than two-level search with blocks.** At L = 1000 they tie (5.69 against
  5.62 s); at L = 100 blocks are 12% faster (5.72 against 6.40 s). STPD over the parse also uses more than twice
  the memory (2.98 against 1.39 GB) and reads twice the files (1,982 against 979 MB). With locating it stays the
  same speed at L = 1000 (5.73 s), but slows to 9.26 s at L = 100.
- **RLZAP saves space but costs time.** It cuts STPD over the parse's files by 720 MB, to 1,263 MB, and its peak to
  2.26 GB. It slows the search by 52% at L = 1000 and 59% at L = 100, and locating at L = 100 by 94%.
- **Character-level STPD through the parse oracle** saves 2.9 GB against the text, but is 5 times slower.
- **katka's STPD over the text is still the fastest search** (3.25 s; the same with `--locate`, since `5f1c76c` no
  longer computes intervals for it). It has to hold the 4 GB text, though.

# 32 individuals (c65)

The same run (`run65.sh`, `summary65.txt`, `run65.log`, `logs65/`) on c65 (n = 8.0 × 10⁹, r = 72.6 M).

## Building the components

| step | wall | peak | output |
|---|---|---|---|
| `index --parse-oracle` | 199.0 s | 2.16 GB | `c65.ptxt` 1,886.2 MB, see below |
| `index --stpd --parse-oracle` | 1,230.4 s | 7.10 GB | `c65.stpd` 581.3 MB (sample 193.4 + φ⁻¹ 387.8) |
| `index --stpd-parse` | 119.1 s | 5.98 GB | `c65.pstpd` 145.8 MB (11,521,370 sampled of r_Z = 12,954,298) |
| `index --rlzap` | 78.7 s | 5.21 GB | `c65.rptxt` 443.8 MB |

- **Parse oracle:** 1,886.2 MB in all, 1.89 bits per character. Its parts:
  - dictionary 170.3 MB: 673,944,330 bases at 2 bits; 177 phrases with other characters;
  - dictionary offsets 35.7 MB;
  - parse 1,302.5 MB (434,172,121 phrases, 9,522,052 distinct);
  - phrase starts 377.7 MB.
- **RLZAP parse part:** 237.9 MB, against 1,680.2 MB for the plain parse and its starts.

## Queries (search time, wall including loading, peak; L = 1000 / L = 100)

| configuration | search (s) | wall (s) | peak | files read (MB) |
|---|---|---|---|---|
| plain search (RLCSA) | 50.30 / 49.29 | 52.1 / 50.1 | 1.69 GB | 1,410 |
| **two-level with blocks** | **5.39 / 5.57** | 6.3 / 6.5 | 2.07 GB | 1,660 |
| STPD over the parse, plain parse | 5.85 / 6.50 | 8.5 / 8.8 | 4.56 GB | 3,490 |
| STPD over the parse, plain parse, `--locate` | 5.97 / 11.38 | 8.9 / 14.2 | 4.57 GB | 3,937 |
| STPD over the parse, RLZAP parse | 8.96 / 10.40 | 10.9 / 12.1 | 3.12 GB | 2,048 |
| STPD over the parse, RLZAP parse, `--locate` | 9.31 / 22.85 | 10.9 / 24.4 | 3.13 GB | 2,495 |
| character-level STPD through the parse oracle | 15.93 / 18.19 | 17.7 / 19.8 | 2.97 GB | 2,685 |
| katka's STPD over the text | 3.18 / 3.43 | 13.7 / 14.5 | 8.87 GB | 8,803 |
| katka's STPD over the text, `--locate` | 3.24 / 3.65 | 14.2 / 15.0 | 8.87 GB | 8,803 |

STPD over the parse reads these components, in MB:
- **the parse:** plain parse 1,302.5, starts 377.7, dictionary 170.3 plus offsets 35.7;
- **the samples:** parse sample 41.8, entry 34.5, C 9.9;
- **the character-level sample:** 193.4;
- **the right end:** RLCSA 1,192.7, B 14.7, dictionary hash 116.7;
- **for `--locate` only:** parse φ⁻¹ 59.6 and character φ⁻¹ 387.8.

**Correctness:** every configuration gives the same MEMs as plain search, 18,479 and 30,809. Every `--locate` run
gives occurrence counts equal to ep − sp. The RLZAP runs give byte-identical output to the plain parse's.

**From 16 to 32 individuals:**
- **Search times:** two-level with blocks (5.62 → 5.39 s at L = 1000), STPD over the parse (5.69 → 5.85 s) and
  STPD over the text (3.25 → 3.18 s) all stay flat. Plain search grows (40.7 → 50.3 s).
- **Peak memory:** STPD over the parse grows from 2.98 to 4.56 GB, two-level with blocks from 1.39 to 2.07 GB,
  and STPD over the text from 4.95 to 8.87 GB.

So the c33 conclusion holds: STPD over the parse is not faster than blocks, and needs more than twice the memory.

# 16 individuals with a sparser parse: closed syncmers k = 200, s = 20 (2026-11-04, a22b892)

`run33_k200s20.sh` rebuilt c33's two-level index with `--syncmer 200,20`, in a separate directory with links to
c33's text, BWT and character-level STPD sample. Then it ran `stpd_parse_chr19.sh` as above
(`summary33_k200s20.txt`, `run33_k200s20.log`, `logs_sparse/`).

## The two parses

| | k, s = 51, 12 | k, s = 200, 20 |
|---|---|---|
| phrases m / distinct σ | 220.5 M / 8.51 M | **48.4 M** / 3.80 M |
| `index --preset chr19` | 214 s, 4.6 GB | 172 s, 3.28 GB |
| `.lvl`: parse + blocks | 166.4 + 42.2 MB | **69.1 + 19.4 MB** (88.5 MB in all) |
| parse oracle: parse + starts + dictionary (+ offsets) | 661.4 + 191.8 + 151.2 (+ 31.9) MB | 133.1 + 55.5 + 278.3 (+ 14.7) MB |
| RLZAP parse part (instead of parse + starts) | 133.4 MB, 17.3 phrases per factor | 68.9 MB, 5.13 phrases per factor |
| parse STPD (sample + entry + C; + φ⁻¹) | 71.4 (+ 45.5) MB | 29.3 (+ 18.6) MB |

## Queries (search time, single-threaded, CPU 9; L = 1000 / L = 100)

| configuration | k, s = 51, 12 | k, s = 200, 20 | peak, 200,20 |
|---|---|---|---|
| plain search (RLCSA) | 40.66 / 39.18 s | 47.40 / 42.24 s (the same search; the difference is run-to-run noise) | 1.05 GB |
| **two-level with blocks** | **5.62 / 5.72 s** | 6.52 / 5.95 s | 1.14 GB |
| STPD over the parse, plain parse | 5.69 / 6.40 s | 7.00 / 12.77 s | 2.17 GB |
| … with `--locate` | 5.73 / 9.26 s | 7.07 / 15.14 s | 2.17 GB |
| STPD over the parse, RLZAP parse | 8.65 / 10.17 s | 8.89 / 26.64 s | 2.05 GB |
| … with `--locate` | 8.98 / 17.96 s | 9.05 / 31.62 s | 2.05 GB |
| character-level STPD through the parse oracle | 16.19 / 19.25 s | 13.48 / 17.13 s | 1.48 GB |

**Correctness:** every configuration gives the same MEMs as plain search (22,711 and 38,984). Every `--locate` run
gives counts equal to ep − sp. The RLZAP runs give byte-identical output to the plain parse's.

## The parse-level part of the search (from `--stats`; there are no dedicated counters)

For STPD over the parse:
- **"Right-end" and "left-end" character steps** are the characters crossed at character level, after the pattern's
  last syncmer and before its first. Per search they are an upper bound on the length added at each end to turn
  the parse-level match into an exact MEM, including one failing step.
- **"Parse-sample binary searches"** count the entries into the parse level.

| per search (STPD over the parse) | 51,12, L = 1000 | 200,20, L = 1000 | 51,12, L = 100 | 200,20, L = 100 |
|---|---|---|---|---|
| searches | 70,246 | 70,246 | 101,241 | 101,241 |
| parse-level entries (parse-sample binary searches) | 157,005 | 90,276 | 181,242 | 86,456 |
| right-end characters | 60.3 | **238.2** | 59.7 | **96.2** |
| left-end characters | 15.9 | **75.3** | 18.8 | **119.4** |
| middle phrase comparisons | 198.4 | 41.6 | 132.0 | 27.8 |

For two-level search with blocks at L = 1000, the sparse parse raises character steps from 6.80 M to 28.23 M, and
lowers parse steps from 2.55 M to 0.85 M and block steps from 1.81 M to 0.42 M. At L = 100 the character steps go
from 10.26 M to 25.73 M.

Without timers, the time of the end steps can't be separated from the parse-level time. An upper estimate: a
character step costs about plain search's time per step, 40.7 s / 259.8 M = 0.16 µs. At L = 1000 the sparse parse's
22.0 M end steps of STPD over the parse would then take about 3.4 s of its 7.0 s, against 0.8 s of 5.7 s at 51,12.

## Reading

**Space:** the sparser parse shrinks the two-level structures to 42% (209 → 88.5 MB) and the parse to 4.6 times fewer
phrases. But the dictionary of longer phrases grows from 151 to 278 MB, so the parse oracle only drops from 1,036
to 482 MB, and RLZAP's factors get shorter (5.1 against 17.3 phrases).

**Time:**
- **Two-level with blocks:** a little slower (6.52 against 5.62 s at L = 1000; 5.95 against 5.72 s at L = 100).
  Longer phrases mean longer character-level ends.
- **STPD over the parse:** slower too, sharply at L = 100 (12.8 against 6.4 s). Its ends grow from about 76 to 314
  characters per search at L = 1000.
- **Character-level STPD through the parse oracle:** the only configuration that gets faster, since its seeks find
  longer phrases (13.5 against 16.2 s).

# Locating with two-level search and blocks (bb42c5a, 65e46c4; 2026-11-04)

`index --twolevel-locate PREFIX` builds `PREFIX.tloc`, then `find --format mems --csa --twolevel --locate` runs at
L = 1000 and 100. Both run in the 11 GB cgroup, and searches are pinned to CPU 9 (`tloc.sh`, `locpos.sh`,
`logs_tloc/`).

## Build

| | c33 (16 individuals) | c65 (32 individuals) |
|---|---|---|
| `index --twolevel-locate` | 432 s, peak 4.21 GB | 828 s, peak 7.20 GB |
| `.tloc` in all | **888.7 MB** (the NOTES estimate: about 880) | 1,180.8 MB |
| characters: run samples + φ⁻¹ | 273.8 + 350.2 = 624.0 MB | 299.4 + 387.8 = 687.3 MB |
| parse: run samples + phrase → text | 45.5 + 191.8 = 237.3 MB | 59.6 + 377.7 = 437.3 MB |
| blocks: run samples + block → phrase | 10.1 + 17.4 = 27.5 MB | 18.4 + 37.9 = 56.3 MB |

## Search (single-threaded; L = 1000 / L = 100)

| | c33 | c65 |
|---|---|---|
| two-level with blocks, no locating | 5.62 / 5.72 s (1.39 GB) | 5.39 / 5.57 s (2.07 GB) |
| **two-level with blocks, `--locate`** | **7.50 / 8.10 s** (2.22 GB) | **7.13 / 7.72 s** (3.17 GB) |
| STPD over the parse, `--locate` | 5.73 / 9.26 s (2.98 GB) | 5.97 / 11.38 s (4.57 GB) |
| occurrences located | 72,876 / 384,341 | 64,945 / 362,636 |
| MEMs in more than one sequence | 9,629 / 19,653 | 7,147 / 14,082 |

The two-level times are from the first runs. In the `--locate-out` runs below, two-level locating took the same
time: 7.50 / 8.25 s on c33 and 7.05 / 7.83 s on c65.

Files read:
- **Two-level locating:** about 1,868 MB for c33 (RLBWT 191 + RLCSA 579 + `.lvl` 209 + `.tloc` 889) and about
  2,841 MB for c65.
- **STPD over the parse with `--locate`:** 2,378 MB and 3,937 MB.

## Checks

- **Against plain search:** columns 1–5 equal plain search's.
- **Counts:** column 6, the occurrences, equals ep − sp for every MEM.
- **Against STPD over the parse with `--locate`:** the occurrence and distinct-sequence counts are equal, MEM by MEM.
- **The occurrence totals** at L = 1000 and 100 on c33 are also those of the published STPD-index on the same
  queries (above).
- **Positions:** `--locate-out` writes each MEM's sorted positions, those of its reverse complement in S, for BML.
  The two-level locating files and STPD over the parse's files are **byte-identical** at both L on both collections:

| | MEMs | file |
|---|---|---|
| c33, L = 1000 | 22,711 | 1,317,329 bytes |
| c33, L = 100 | 38,984 | 5,036,107 bytes |
| c65, L = 1000 | 18,479 | 1,140,044 bytes |
| c65, L = 100 | 30,809 | 4,660,604 bytes |

## Reading

- **Locating with two-level search:** it adds about 2 s to the 5.6–5.7 s search (34–42% on c33, 32–39% on c65),
  and 0.9–1.2 GB of files.
- **Against STPD over the parse with `--locate`:** at L = 1000 it is slower (7.5 against 5.7 s on c33). At L = 100
  it is faster (8.1 against 9.3 s on c33, 7.7 against 11.4 s on c65). It needs 25–30% less memory, and reads about
  21% (c33) to 28% (c65) fewer bytes.
