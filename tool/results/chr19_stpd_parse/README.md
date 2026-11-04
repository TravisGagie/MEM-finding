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
