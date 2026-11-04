# Reproducing the chromosome 19 experiments of *More RLBWT Tricks* (v2)

This document is for reproducing Section 5 of *More RLBWT Tricks* (Travis Gagie, arXiv:2609.40299v2) without asking us: Table 1, Table 2 and the paragraph on locating. Everything runs with one tool, `katka`, in this directory. One script, `experiments/chr19_paper.sh`, builds the indexes, runs every configuration, checks every output against plain backward search and writes the tables in the paper's layout. A synthetic stand-in for the data checks the whole pipeline in under a minute (see the last section).

## The code

- **Repository:** https://github.com/TravisGagie/MEM-finding, branch `tool`, directory `tool/`. Use the commit that last changed this file (`git log -1 -- tool/CHR19.md`).
- **The paper's own runs** used earlier commits of the same branch: `5a3d161` for Table 1 (16 and 32 individuals) and `9550a7a` for 64 individuals; `4700b46` for Table 2; `bb42c5a` and `65e46c4` for locating. Their logs are in `results/chr19/`, `results/chr19_stpd_parse/` and `results/chr19_rlz/`. Since then nothing has changed in what the searches compute; the self-test checks that every configuration still gives the same MEMs, intervals and occurrences as plain search.
- **Requirements:**
  - a C++17 compiler (GCC 11 or later; the paper used GCC 15), make, zlib and Python 3;
  - [Big-BWT](https://github.com/alshai/Big-BWT) for the BWTs of the real collections (`make` in its directory; tested at commit `5b00650`). Without it, `katka bwt` builds the BWT with libsais, but needs about 6 bytes per character: 24 GB at 16 individuals, 96 GB at 64;
  - for the last two rows of Table 2 only, the published STPD-index (below).
- **Hardware:** the paper used an Intel Core i7-8700 with 14 GiB of RAM. Every step ran alone under an 11 GB memory cap with no swap (`systemd-run --user --scope -p MemoryMax=11G -p MemorySwapMax=0`; the script uses it if it is available), and searches ran single-threaded, pinned to one core. The 64-individual index build peaks at 8.3 GB.
- **Disk:** about 12 GB per 16 individuals for the text, its BWT and the indexes (the 64-individual text and BWT are 16 GB each). Keep at least 60 GB free for all three sizes.

```
git clone -b tool https://github.com/TravisGagie/MEM-finding.git && cd MEM-finding/tool
make                       # bin/katka, bin/katka-test, bin/katka-blocks; third_party/libsais is included
./selftest.sh              # about 10 minutes; ends with "self-test passed"
git clone https://github.com/alshai/Big-BWT.git ~/Big-BWT && make -C ~/Big-BWT
```

## Data

The collections are chromosome 19 of CHM13 and of both haplotypes of the first 16, 32 or 64 individuals of HPRC release 2 whose haplotypes each have it as one sequence of at least 90% of CHM13's length (33, 65 and 129 sequences), all in CHM13's orientation. The queries are haplotype 1 of HG00423 and HG00438, which are left out of the collections. The tool indexes each sequence and its reverse complement.

`experiments/hprc_chr19.py` downloads them and writes the two files the script reads: `collection.txt` (one sequence per line, CHM13 first, then both haplotypes of each individual) and `queries.txt` (one per line). Every character other than A, C, G and T becomes N, and each sequence is checked against its length and number of Ns. It fetches only the bases of chromosome 19, block by block from the bgzipped assemblies on S3. With 16 individuals that's 35 sequences, about 2.1 GB of bases, a few minutes:

```
python3 experiments/hprc_chr19.py ~/chr19data                  # 16 individuals + the 2 queries
python3 experiments/hprc_chr19.py ~/chr19data --plan           # only list the sequences
```

**32 and 64 individuals.** The script's built-in list (from `chr19/twolevel/download_hprc_chr19.py`) has the first 16 individuals. `experiments/hprc_chr19_64.tsv` lists the 129 sequences of the paper's 64-individual collection and the two queries: role, sample, haplotype, sequence name, length and number of N, then the source assembly, with the individuals in order. Its first 16 individuals are the built-in list. For all three sizes, download the 64 individuals (129 sequences, about 7.9 GB of bases):

```
python3 experiments/hprc_chr19.py ~/chr19data --manifest experiments/hprc_chr19_64.tsv
```

`--individuals N` takes only the first N individuals of a manifest.

Use the 64-individual `collection.txt` for every size: the script takes its first 33, 65 or 129 lines.

## Running

```
export BIGBWT=~/Big-BWT/bigbwt PIN=9               # PIN: a core to pin the searches to
T=~/MEM-finding/tool/experiments/chr19_paper.sh
$T build    ~/chr19data ~/chr19work               # indexes for 16, 32 and 64 individuals (33 65 129 sequences)
$T table1   ~/chr19data ~/chr19work               # Table 1        -> ~/chr19work/table1.md
$T table2   ~/chr19data ~/chr19work               # Table 2        -> ~/chr19work/table2.md
$T locating ~/chr19data ~/chr19work               # the locating paragraph -> ~/chr19work/locating.md
```

`$T all` does the four in order. A last argument restricts the sizes, for example `"33"` for 16 individuals only. Table 2 and the locating paragraph use the smallest size. Each command is restartable: the build skips indexes that exist. Every search's MEMs are compared with plain backward search's, and a difference stops the script.

### The build

For each size n, in `~/chr19work`:

| step | command | output | the paper's time and peak, 16 / 32 / 64 individuals |
|---|---|---|---|
| prep | `katka prep listN.tsv cN` | `cN.S` (both strands, 4.07 / 8.00 / 15.9 × 10⁹ characters), `cN.tbl` | 18 / 38 / — s, 0.13 GB |
| BWT | `bigbwt cN.S` | `cN.bwt` | 345 / 599 / — s, 2.6 / 3.5 / — GB |
| index | `katka index --preset chr19 cN` | `cN.rlbwt`, `cN.csa` (RLCSA), `cN.lvl` (two-level index with blocks) | 4.6 / 5.2 / 8.3 GB |

For the smallest size, the build also runs:
- `katka index --stpd --parse-oracle`: `cN.ptxt` (the parse as a text oracle, with the dictionary at 2 bits per base) and `cN.stpd` (the character-level STPD sample, built reading the text through the oracle);
- `katka index --stpd-parse`: `cN.pstpd`;
- `katka index --rlzap`: `cN.rptxt`, the RLZ parse oracle;
- `katka index --twolevel-locate`: `cN.tloc`.

The log of each step ends with `wall W s, peak P GB, exit 0`. With 16 individuals the index log reports 220.5 million phrases, 8.51 million distinct, and 23.2 million blocks of 9.52 phrases on average; `c33.lvl` is 209 MB.

### Table 1

For each size, `table1` runs `katka find --format mems --stats -L L c$n queries.fa out` with:
- `--csa` (plain backward search on the RLCSA);
- `--csa --twolevel --no-jumps` (two-level, no block steps);
- `--csa --twolevel` (two-level with blocks);

each for L = 1000 and 100, and `--csa --twolevel --fb` (forward–backward with blocks) for L = 1000. `queries.fa` is the two queries split at runs of N: 4 patterns of 122.4 million characters in all. The time reported is the tool's own search time (`wall … s` on its `reads=` line), which excludes loading the index. The paper's Table 1:

| | 16 individuals | 32 individuals | 64 individuals |
|---|---|---|---|
| plain backward search (RLCSA) | 40.4 / 39.2 | 50.7 / 49.8 | 59.6 / 56.9 |
| two-level, no block steps | 8.46 / 8.57 | 8.95 / 9.03 | 9.45 / 9.74 |
| two-level with blocks | 5.66 / 5.84 | 5.54 / 5.65 | 5.65 / 5.80 |
| forward–backward with blocks | 5.65 | 5.36 | 5.60 |
| two-level index (MB) | 209 | 250 | 318 |
| peak memory of the index build (GB) | 4.6 | 5.2 | 8.3 |

`table1.md` has the same rows and columns. With 16 individuals there are 22,711 MEMs at L = 1000 and 38,984 at L = 100. Times on other machines will differ; the MEMs must not. With `--stats` the logs give the step counts the text quotes: at L = 1000 with 16 individuals, plain search takes 260 million character steps, and two-level search with blocks 6.8 million character steps, 2.5 million parse steps and 1.8 million block steps.

### Table 2

With 16 individuals, `table2` runs the same searches for these configurations (all `katka find --format mems -L L`):

| configuration | options | the paper's time (s), peak (GB), files read (MB) |
|---|---|---|
| plain backward search (RLCSA) | `--csa` | 40.7 / 39.2, 1.05, 770 |
| two-level with blocks | `--csa --twolevel` | 5.62 / 5.72, 1.39, 979 |
| with locating | `--csa --twolevel --locate` | 7.50 / 8.10, 2.22, 1,868 |
| STPD over the text | `--stpd --no-intervals` | 3.25 / 3.57, 4.95, 4,784 |
| with locating | `--stpd --locate` | 3.23 / 3.76, 4.95, 4,784 |
| STPD over the text, read through the parse | `--stpd --parse-oracle --no-intervals` | 16.2 / 19.3, 2.04, 1,755 |
| STPD over the parse, plain parse | `--csa --stpd-parse --no-intervals` | 5.69 / 6.40, 2.98, 1,982 |
| with locating | `--csa --stpd-parse --locate` | 5.73 / 9.26, 2.98, 2,378 |
| STPD over the parse, RLZ parse | `--csa --stpd-parse --rlzap-oracle --no-intervals` | 8.65 / 10.2, 2.26, 1,263 |
| published STPD-index (RLZ text) | `experiments/stpd_index_chr19.sh` | 44.6 / 54.8, 0.96, 613 |
| with locating | (the same, `--locate`) | 53.8 / 70.2, 0.96, 613 |

- **The columns:** "files" is the size of what the search reads. It counts whole files (RLBWT, RLCSA, `.lvl`, `.tloc`, the text `.S`, `.stpd`, `.ptxt`), and, for STPD over the parse, the components the tool lists on its `STPD over the parse, components (MB)` line. The peak is that of the search at L = 1000.
- **Other numbers in the text:** with 32 individuals, STPD over the text takes 3.18 s and 8.87 GB, and STPD over the parse 5.85 / 6.50 s and 4.56 GB. For those, run `table2` with the size `"65"`. The sparser parse (k = 200, s = 20) is `katka index --preset chr19 --syncmer 200,20` on a copy of `c33`; the desktop's run is `results/chr19_stpd_parse/run33_k200s20.sh`.
- **The RLZ factorization:** `katka index --rlzap` prints its factors. With 16 individuals: 12.35 million factors of 17.3 phrases, and the parse part shrinks from 853 MB (ids and phrase starts) to 133 MB. Greedy RLZ, 9.4 phrases per factor, is in `results/chr19_rlz/` (`rlz_parse.cpp`, an analysis program).

**The published STPD-index** (github.com/regindex/STPD-index, commit `8a432a4`; its rows need `STPD=path/to/the/built/clone`):
- **Building it:** its CMake fetches dependencies from `chromium.googlesource.com`, and its first configure needs the generated libdivsufsort include directory added:
  ```
  git clone https://github.com/regindex/STPD-index && cd STPD-index && git checkout 8a432a4 && mkdir build && cd build
  cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-I$PWD/_deps/sdsl-build/external/libdivsufsort/include -Wno-template-body" \
           -DCMAKE_PROJECT_INCLUDE=$HOME/MEM-finding/tool/results/chr19/run_scripts/linker_lang_shim.cmake
  make -j4
  ```
  The shim and `-Wno-template-body` were needed with CMake 3.31 and GCC 15 (`results/chr19/README.md`); older versions may not need them.
- **Its index:** its own construction builds a suffix tree, which does not fit in memory at this scale. So `stpd_index_chr19.sh` builds its input files with `katka stpd-prep` and `katka stpd-export`, which are identical to those its construction writes on smaller texts. Then their `stpd-build` builds the index. With 16 individuals that took about 65 minutes in all, most of it in `stpd-export`.
- **Searching:** `experiments/stpd_index_mems.cpp` is a BML driver around its longest-prefix search, built automatically. `table2` runs only their BML, with and without locating: their forward–backward search did not finish in 54 minutes.

### The locating paragraph

`locating` reports:
- the size of `c33.tloc`: the samples at the run starts of the three BWTs, φ⁻¹, the block → phrase and phrase → text lists. The paper reports 889 MB with 16 individuals and 1,181 MB with 32;
- the times with and without locating (from Table 2's runs);
- the number of occurrences of all MEMs: 72,876 at L = 1000 and 384,341 at L = 100 with 16 individuals. It checks these equal the sums of ep − sp of plain search;
- whether the positions of every MEM are byte-identical with two-level locating and with STPD over the parse with locating. It writes both with `--locate-out`: per MEM, name, start, end and the sorted positions of its reverse complement in S. They must be.

Expected `locating.md`, with 16 individuals:

```
- the samples, phi^-1 and the two lists (PREFIX.tloc): 889 MB
- two-level with blocks: 5.62 / 5.72 s without locating, 7.50 / 8.10 s with (L = 1000 / 100; Table 2)
- STPD over the parse with locating: 5.73 / 9.26 s
- occurrences of all MEMs: 72876 (L = 1000), 384341 (L = 100); equal to the sums of ep - sp of plain search: True
- positions (--locate-out) byte-identical to STPD over the parse, MEM by MEM: L = 1000 True, L = 100 True
```

With 32 individuals the paper reports 7.13 / 7.72 s with locating instead of 5.39 / 5.57 s, against 7.72 / 11.4 s for STPD over the parse: run `locating` with the size `"65"`.

## Using two-level search with blocks from code

`experiments/twolevel_example.cpp` is a complete program, about 90 lines. It loads an index built by `katka index --preset chr19 PREFIX` (and optionally `--twolevel-locate`) and, for each pattern of a FASTA file, does two things:
- runs one longest-suffix query;
- lists the MEMs of length at least L, by BML, with the BWT interval of each MEM's reverse complement and, optionally, its occurrences.

```
g++ -std=c++17 -O3 -march=native -Isrc -Ithird_party/libsais/include experiments/twolevel_example.cpp \
    build/libsais.o build/libsais64.o -o bin/twolevel-example -lz -pthread
bin/twolevel-example ~/chr19work/c33 ~/chr19work/queries.fa 1000 --locate
```

Its MEMs and intervals are byte-identical to `katka find`'s. Everything is header-only, in `src/`. The pieces:
- **The index:** `rlbwt` (`.rlbwt`), `rlcsa` (`.csa`, for character steps), `two_level` (`.lvl`: B, the parse BWT, the dictionary, the blocks), and optionally `twolevel_locate` (`.tloc`). A `stepper` bundles the RLBWT and RLCSA.
- **`segment_search`** holds a pattern segment, its parse into phrases (`prepare`), and scratch space: one per thread. Set `Q.TL` to search over the two-level index, and `Q.jumps = true` for block steps.
- **The primitive:** `Q.back(dir, lo, end, I)` returns the start j of the longest suffix of X[lo..end) that occurs in the text, X being the segment (dir 0) or its reverse complement (dir 1), and its BWT interval I. That's the same j and I as plain backward search (`search2` in `src/twolevel.hpp`). With a `twolevel_locate`, `search2` also returns SA at the interval's first row (a toehold).
- **MEMs:** `find_mems(Q, read, options, report)` runs BML (`options.L`) or forward–backward search (`options.fb`) over the A/C/G/T runs of a read. It calls `report(mem)` with [a, e) and the interval of the MEM's reverse complement, which occurs in the text since the text holds both strands. With `Q.locate` and `Q.TLC` set, `Q.occbuf` holds the occurrences when `report` is called.

## Checking it end to end on synthetic data

`experiments/chr19_synthetic.py` writes a small stand-in for the data: a random reference and haplotypes derived from four founders, with SNPs, short indels and runs of N. The whole pipeline then runs in about 10 seconds:

```
python3 experiments/chr19_synthetic.py /tmp/c19/data --length 300000 --individuals 8
experiments/chr19_paper.sh all /tmp/c19/data /tmp/c19/work "5 9 17"       # 2, 4 and 8 "individuals"
```

It ends with `done: all`, after printing the three results. Every search's MEMs equal plain search's, the occurrence counts equal ep − sp, and the positions of two-level locating and STPD over the parse are byte-identical (`True` on the last line). The times are meaningless at this size. We also ran it on real data at small scale: CHM13 and both haplotypes of HG00097 (3 sequences) with the two real query haplotypes, downloaded with `hprc_chr19.py --manifest` (see `results/chr19_paper_check/`).
