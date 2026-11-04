# KATKA, the tool

One index and one search engine for the three experiments of *Three Examples of MEM-Finding*: long noisy reads with few taxa (Zymo, Section 3), short reads with many taxa (SILVA, Section 4) and long accurate matches (chromosome 19, Section 5). It merges the code of KATKA (`../katka-readable/`) and of *More RLBWT Tricks* (`../chr19/twolevel/`), on one text holding every sequence and its reverse complement. The old directories are unchanged, and they still reproduce the numbers in the current paper.

Every technique is an option, and a preset chooses the options for each use case. Every search method returns the same longest suffix and BWT interval for every query, so any combination of them gives the same MEMs. `selftest.sh` checks this against brute force.

## Building

You need a C++17 compiler, make and zlib. libsais is in `third_party/` (Apache 2.0). sdsl-lite is not needed.

```
make            # bin/katka and bin/katka-test
./selftest.sh   # builds, then checks everything on small synthetic data (about 2 minutes)
```

The build uses `-march=native`. Run `make CXXFLAGS="-std=c++17 -O3 -DNDEBUG"` for a portable binary.

## Quick start

```
katka prep list.tsv ref            # ref.S (both strands), ref.tbl, ref.info
katka bwt ref                      # ref.bwt = BWT(ref.S $), with libsais (or bring one from Big-BWT / pfp-merge)
katka index --preset silva ref     # ref.rlbwt, ref.csa, ref.tag, ref.tag.rmq
katka find --preset silva ref reads.fq out.listings
```

`list.tsv` has one line per FASTA file (plain or gzipped), in the order the genomes should have in the text: `<tag>\t<path>`. For SILVA, that is tree order.

## Commands

### `katka prep [--preset NAME] [-d] [-s] [-U] [--digest k,w] list.tsv PREFIX`

Builds the text `S = G_1 X rc(G_1) X G_2 X rc(G_2) X ...`, where every character other than A, C, G or T becomes `X`. It writes `PREFIX.S`, `PREFIX.tbl` (index, tag, name, length, start in S) and `PREFIX.info`. The text and table are KATKA's `katka-prep` output, byte for byte.

| option | meaning |
|---|---|
| `-s` | a FASTA file may hold several genomes, a genome being a run of consecutive records with the same sample ID (the header up to the first `.` or blank) |
| `-d` | also write one text per dataset (consecutive genomes with the same tag), `PREFIX.ds<j>.S`, listed in `PREFIX.datasets`, for Big-BWT and pfp-merge |
| `-U` | homopolymer compression of every genome; `katka find` then compresses the reads too (recorded in `PREFIX.info`) |
| `--digest k,w` | index symmetric minimizer digests instead of the bases (App. B; the paper uses 4,11): a canonical hash, all tied minimizers in each window, equal neighbours collapsed, the rarest codes sharing one byte if more than 252 occur. Writes the byte map `PREFIX.map`; `katka find` then digests the reads too. The text, table and map are the same as KATKA's `rz-prep -B` |

### `katka bwt [-d] PREFIX`

Builds `PREFIX.bwt` with libsais: one byte per row, with the terminator as byte 0. With `-d`, the BWT has one string per dataset, as pfp-merge builds it, for at most 60 datasets and under 2^31 characters. For large texts, use Big-BWT (`bigbwt PREFIX.S && mv PREFIX.S.bwt PREFIX.bwt`) or Big-BWT plus pfp-merge, as KATKA did. Terminators may be bytes 0, 1 or 2.

### `katka index [--preset NAME] [options] PREFIX`

Reads `PREFIX.S`, `PREFIX.tbl` and the BWT, and always writes the RLBWT `PREFIX.rlbwt`.

| option | writes | meaning |
|---|---|---|
| `--bwt FILE` | | the BWT, if not `PREFIX.bwt` |
| `--csa` | `PREFIX.csa` | explicit RLCSA: per-character run lists with jump tables |
| `--efcsa` | `PREFIX.efcsa` | Elias–Fano RLCSA: run starts shifted by c·n, so one sequence serves every character |
| `--tags` | `PREFIX.tag`, `.tag.rmq` | the complete, run-length compressed tag array, built without a suffix array (LF pass, then φ) |
| `--by-tag` | | tags are the tag column of the table, numbered by first appearance (species for Zymo); by default each genome is its own tag (one genus per file for SILVA) |
| `--sample s` | | store the tag of about one run in s, recover the others by LF (one sample per cycle) |
| `--no-rmq` | | no RMQ: counting only, with no listing |
| `--grammar` | | grammar-compress the runs' tags with the built-in approximate RePair |
| `--bigrepair P` | | grammar-compress them with `bigrepair -i` output `P.R`, `P.C` (as KATKA did) |
| `--tag-out FILE` | | write the tag array somewhere other than `PREFIX.tag` |
| `--tag-lcp l1,l2,…` | | also report B_λ for each λ: the tag boundaries whose two suffixes share more than λ characters. It compares the text at each boundary while φ streams the suffix array, and prints the bound (r + σ + 2d)(λ + 1) + B_λ + 1 of Lemma 11 against the number of tag runs |
| `--kebab k` | `PREFIX.k<k>.kbf` | KeBaB: a blocked Bloom filter of the canonical k-mers of S |
| `--kebab-eps e` | | its false-positive rate (default 0.1) |
| `--twolevel` | `PREFIX.lvl` | two-level index with blocks; needs a BWT with one terminator |
| `--stpd` | `PREFIX.stpd` | suffix-tree path decomposition sample for `find --stpd`: the needed BWT run starts, from the LF pass and text comparisons; needs one terminator |
| `--stpd-all` | | keep every BWT run start (no pruning) |
| `--stpd-parse` | `PREFIX.pstpd` | STPD over the parse: the parse sample, entry table, C and parse φ⁻¹, from `PREFIX.ptxt` (built first if missing) and `PREFIX.lvl`'s k, s; alone, it reads no BWT |
| `--twolevel-locate` | `PREFIX.tloc` | samples at the run starts of the BWT, the parse BWT and the block BWT, block → phrase and phrase → text maps, and φ⁻¹, for `find --twolevel --locate`; needs `PREFIX.lvl` with blocks and the BWT |
| `--rlzap[=J,D,M]` | `PREFIX.rptxt` | the parse oracle with the parse as an RLZAP-style factorization against the first sequence's parse, both strands (default 4,4,2); from `PREFIX.ptxt`; alone, it reads no BWT |
| `--parse-oracle` | `PREFIX.ptxt` | the text read through the two-level parse (phrase starts, parse, dictionary), for `find --stpd --parse-oracle`; uses `--syncmer`; run alone, it reads only `PREFIX.S` |
| `--syncmer k,s` | | closed syncmers (default 51,12), canonical under reverse complement |
| `--no-blocks` | | the parse only |
| `--continuation smallest\|frequent` | | which successor continues a block (default smallest; frequent stores the order bit ω) |

Tag runs are cached in `PREFIX.tagruns` (or `PREFIX.sp.tagruns` with `--by-tag`), so other kinds of tag array reuse them.

### `katka find [--preset NAME] [options] PREFIX reads out`

Finds the MEMs of each read and answers each one from the tags. Reads may be FASTQ or FASTA, plain or gzipped. Write `-` as `out` for standard output.

| option | meaning |
|---|---|
| `-L n` | MEMs of at least n bases; with `-X` or `--top`, this is the starting value (default 1) |
| `-X`, `--adaptive` | BML with adaptive L: keep L equal to the length of the longest MEM found so far. All the longest MEMs, ties included, are reported, plus the shorter ones found before L grew (as KATKA does) |
| `--top t` | report only the t longest MEMs: keep the t longest found so far and set L to one more than the shortest of them |
| `--longest t` | the t longest MEMs with ties: every MEM at least as long as the t-th longest (L is raised to the t-th longest found so far). With `-K`, the streaming rule: pseudo-MEMs longest first, until the next is shorter than L (exact if L ≥ k). Overrides `-X` |
| `--kebab-top t` | with `-K`: `--longest t` over the t longest pseudo-MEMs only (KeBaB's heuristic; no guarantee) |
| `--kebab-batch t` | with `-K`: `--longest t` over batches of pseudo-MEMs (the 2t longest, then 4t, 8t, ...), stopping once t MEMs passing `-D` are longer than the next pseudo-MEM; if they run out with fewer than t MEMs of length ≥ k and L < k, BML on the whole read (exact) |
| `--kebab-check T t` | with `-K`: `--longest t` over the T longest pseudo-MEMs (T ≥ t); accepted if there are at least T + 1 pseudo-MEMs and t MEMs at least as long as the (T+1)st (`--strict`: longer), else BML on the whole read (exact) |
| `--check-seed` | with `--kebab-check`: the fallback starts with L = the t-th longest MEM found in the scanned pseudo-MEMs, if at least t were found (same output) |
| `--fixed-L` | with `--kebab-top t`: every MEM of length ≥ max(L, k) in the t longest pseudo-MEMs (KeBaB's pipeline) |
| `--read-stats FILE` | with the three above: per read, pseudo-MEMs, those searched, batches, fallback, backward steps, MEMs (TSV) |
| `--fb` | forward–backward search (all MEMs of at least L) instead of BML |
| `-D d` | low-complexity filter: ignore MEMs with fewer than d distinct 3-mers |
| `-U`, `--no-hpc` | homopolymer-compress the reads, or don't (default: as the index was built). Spans are reported in the original read |
| `-K k` | KeBaB pseudo-MEMs with `PREFIX.k<k>.kbf`. With `-X` or `--top`, searched longest first, stopping at the first one shorter than L |
| `--inline-kebab k` | KeBaB window check inside BML: test the window's k-mers lazily from its right end leftward (the rolling codes restarted at the window's last k-mer) and, at the first that fails, jump past it, before any backward step |
| `--kbf FILE` | the filter file, if not `PREFIX.k<k>.kbf` |
| `--csa`, `--efcsa` | backward steps on the explicit or Elias–Fano RLCSA (default: the RLBWT) |
| `-F t` | lookup table of the intervals of all t-mers, t ≤ 13 (plain search only; over digests, of all t-symbol strings, e.g. `-F 2`) |
| `--twolevel` | two-level search with block jumps (`PREFIX.lvl`); character steps use the RLCSA if `--csa` is given |
| `--no-jumps` | two-level search without block steps |
| `--stpd` | STPD search (`PREFIX.stpd`, with the text `PREFIX.S` as the oracle): binary searches over a sample of at most r text positions instead of backward steps. The intervals of the reported MEMs are then computed by backward search, so the output is the same as plain search |
| `--no-intervals` | with `--stpd` and no tags: skip the intervals (the mems format prints `-` for them) |
| `--stpd-parse` | STPD over the parse (`PREFIX.lvl`'s B and dictionary, `PREFIX.pstpd`, `PREFIX.ptxt`, `PREFIX.stpd`; the RLCSA with `--csa`, else the RLBWT, at the right end). Same MEMs as plain search. With `--locate`, φ⁻¹ over the parse; the mems format then prints the occurrences and the distinct sequences among them |
| `--rlzap-oracle` | with `--stpd-parse` or `--stpd --parse-oracle`: read the text through `PREFIX.rptxt` |
| `--parse-oracle` | with `--stpd`: read the text through `PREFIX.ptxt` instead of holding `PREFIX.S` (same output; see NOTES.md) |
| `--twolevel --locate` | two-level search with blocks keeps a toehold (SA at the interval's first row) through character, parse and block steps (`PREFIX.tloc`), and lists every occurrence by φ⁻¹; the mems format prints the occurrences and the distinct sequences |
| `--locate` | with `--stpd`: locate every occurrence of each MEM (of its reverse complement, for BML and FB) by φ⁻¹ steps from the first one, checked against the text at the end of each batch of doubling size, then a binary search in the batch that overshot, as Becker et al. locate; the mems format gets a sixth column, the number of occurrences |
| `--stpd --twolevel` | STPD finds the MEMs, and two-level search (with blocks unless `--no-jumps`) computes their intervals |
| `-c`, `--count` | each tag with its number of occurrences (proportional credit) |
| `-l`, `--list` | the distinct tags (even credit), by Muthukrishnan–Sadakane listing |
| `--range` | the smallest and largest tags (range scheme, LCA). This is the default when there is a tag array |
| `--no-tags` | no answers |
| `-T FILE` | the tag array, if not `PREFIX.tag` |
| `--format listing\|mems` | Cliffy's listing format extended with counts (default), or one line per MEM, `name start end sp ep` (0-based, half-open), with no tags |
| `-j p` | p threads, which take chunks of 256 reads; the output is the same as with one thread |
| `--stats` | numbers of searches and of character, parse and block steps |

Over a digested index (`prep --digest`), `find` runs BML on the digests with L in bases (L ≥ w): a failed window jumps to the first window that selects only later symbols, and each digested MEM is reported with its span in bases, keeping only the spans not inside another's (App. A, "BML with digests"). Answers, the RLCSAs, the lookup table and threads work as usual; `-X`, `--top`, `--fb`, KeBaB, two-level search and `-U` do not apply.

In the listing format, each read gets `>name`, then `[a,b] {answer}` for each MEM, where `[a,b]` is its span in the read (0-based, inclusive). This is the format of `katka-classify`, so the scorers in `scripts/` (copied from `katka-readable/scripts/`) work unchanged:

```
python3 scripts/credit_schemes.py out.listings truth.txt PREFIX.tbl      # Zymo: proportional, even, range, ...
python3 scripts/norm_credit.py ...                                       # normalised credit
python3 scripts/maxmem_roc.py out.listings truth.txt PREFIX.tbl reads    # -X runs
python3 scripts/score_silva.py --mode freq ...                           # SILVA (freq = proportional)
```

The interval reported for each BML MEM is that of its reverse complement, as in KATKA: the forward search ends there. Its tags and counts are the same as the MEM's. Forward–backward search reports the same interval.

### `katka stpd-prep [-s] list.tsv OUT`, `katka stpd-export OUT`

These write input files for the published STPD-index (github.com/regindex/STPD-index), so that its `stpd-build` can index whole chromosome collections. Its own sampling step builds a suffix tree, at about 68 bytes per character.
- `stpd-prep` writes `OUT.T`: their text, both strands, split at non-A/C/G/T characters and concatenated without separators. It also writes `OUT.S` = reverse(`OUT.T`), the piece starts and `OUT.reflen`, the first sequence's length on both strands, for their `-l`.
- After the BWT of `OUT.S` and `katka index --stpd OUT`, `stpd-export` writes `OUT.T.{colex_m,rbwt,pa,lcs}`, byte-identical to their `stpd-sampling -c -P` on the texts tested. It needs O(r) memory plus the text: φ streams the prefix array, and the LCS array comes from PLCP. `stpd-build -i OUT.T` then finds these files and skips the suffix tree. `experiments/stpd_index_chr19.sh` runs the whole comparison.

### `katka info PREFIX`, `katka presets`

`info` shows the sizes of the index files; `presets` lists each preset's options.

## Presets

Each preset is one flag (`--preset NAME`) for `prep`, `index` and `find`. Options given after it override its settings.

| preset | use case | prep | index | find |
|---|---|---|---|---|
| `zymo` | noisy long reads, few taxa | `-d -s -U` | `--csa --tags --by-tag --no-rmq --kebab 20` | `-U -X -L 20 -K 20 -D 14 --count --csa`: homopolymer compression, θ = 14, adaptive L, KeBaB k = 20 longest first, proportional credit |
| `silva` | short reads, many taxa | | `--csa --tags` | `-L 30 --count --csa -F 10`: explicit RLCSA, 10-mer table, run-length tags, proportional credit |
| `chr19` | long accurate matches | | `--csa --twolevel` | `-L 1000 --twolevel --format mems --csa`: BML, two-level with blocks, no tags (add `-L 100` for L = 100) |
| `novel` | a few long MEMs in mostly novel sequence | | `--csa --kebab 20` | `-L 1000 --inline-kebab 20 --format mems --csa`: fixed large L, inline KeBaB check |

For `zymo`, build the BWT with one string per dataset: run Big-BWT on each `PREFIX.ds<j>.S`, then pfp-merge them, or use `katka bwt -d` for small data.

## Files

| file | contents |
|---|---|
| `PREFIX.S`, `.tbl`, `.info` | the text, the genome table, whether homopolymer-compressed or digested |
| `PREFIX.map` | the digest's byte map (KATKA's format, `RZDGMAP1`) |
| `PREFIX.bwt` | the BWT (input) |
| `PREFIX.rlbwt` | RLBWT: Elias–Fano run starts, a wavelet matrix over the run heads, per-character Elias–Fano cumulative run lengths |
| `PREFIX.csa`, `.efcsa` | RLCSAs |
| `PREFIX.tag`, `.tag.rmq` | tag array (run starts; tags, sampled tags or grammar); RMQ for listing |
| `PREFIX.tagruns`, `.sp.tagruns` | cache of the tag runs |
| `PREFIX.k<k>.kbf` | KeBaB filter; the same format as KATKA's `.kbf` (`RZKBLM01`) |
| `PREFIX.stpd` | STPD sample (sorted text positions) |
| `PREFIX.lvl` | two-level index: B, the parse BWT, the phrase dictionary, and the blocks (succ, BS, block BWT, symbols) |

Apart from the `.kbf`, the formats are new. The tool doesn't read KATKA's `.rz`, `.tag` or `.csa` files, because it doesn't use sdsl.

## Map from the paper to the code

Section numbers are those of *Three Examples of MEM-Finding* (`../paper/`).

| Paper | File | What is there |
|---|---|---|
| §2 "RLBWTs and RLCSAs" | `src/bits.hpp` | Elias–Fano, packed integers, rank/select bitvectors, run-length bitvectors, `EFRank` (the Elias–Fano RLCSA, also used for the BWTs of the parse and of the blocks), wavelet matrix |
| §2 "RLBWTs and RLCSAs" | `src/rlbwt.hpp`, `src/rlcsa.hpp` | the RLBWT (rank, LF, backward step); explicit and Elias–Fano RLCSAs |
| §2 "BWTs"; App. A "Lookup tables" | `src/search.hpp` | intervals, backward steps on either structure, the t-mer lookup table |
| §2 "Both strands"; App. A "Reference text" | `src/prep.hpp` | the text of both strands; datasets; homopolymer compression |
| §2 "MEMs", "Boyer–Moore–Li", Fig. 1, "Longest MEMs", "KeBaB"; App. A "BML without digests" | `src/mems.hpp` | `segment_search` (longest-suffix searches, plain or two-level); BML with fixed or adaptive L or the t longest; forward–backward search; KeBaB pseudo-MEMs and the inline window check |
| §2 "KeBaB", §3 | `src/kmer_filter.hpp` | the Bloom filter (from KATKA) |
| §3 "Long, noisy reads" | `src/mems.hpp` (`distinct_3mers`), `src/main.cpp` (`finder::run`) | the low-complexity filter; homopolymer compression of reads, with spans mapped back |
| §4 "Counting with complete tag arrays" | `src/tag_array.hpp` | counting (proportional credit); listing (even credit, range/LCA); sampled tags via LF |
| §4 "Counting" (listing without the previous-run array) | `src/rmq.hpp` | a range-minimum structure that does not store the array (balanced parentheses of a min-tree, about 6 bits per run) |
| §4 "Construction" | `src/tag_build.hpp` | LF pass for the r-index samples; φ streaming of the suffix array into tag runs; sampling |
| §4 "Grammar-compressing the runs" | `src/grammar.hpp` | approximate RePair, bigrepair input, the encoding of Gagie et al. |
| §5 "Two-level indexing with a syncmer parse", "Blocks" (Lemmas 6–14) | `src/twolevel.hpp` | canonical closed syncmers, the phrase dictionary of both strands, the pattern parse, `search2` (block steps and the back-up rule) |
| §5 "Implementation and setup" | `src/twolevel_build.hpp` | streaming parse of S, B, the parse BWT (libsais), successors, blocks, symbols, block BWT |
| App. B "Minimizer digests"; App. A "Digests", "BML with digests" | `src/digest.hpp`, `src/mems.hpp` (`bml_digest_run`) | symmetric digests, the byte map, the windows' digest ranges and spans; BML on digests |
| Becker et al., "Compressing Suffix Trees by Path Decompositions" (ICALP 2026) | `src/stpd.hpp` | an STPD search backend: the mirror image of their Algorithm 1 (right to left, π = ISA, i.e. their st-colex⁻ of the reversed text); the sample built from the BWT runs |
| App. A "Parallelism"; output | `src/main.cpp` | the command line, presets, threads writing in input order |
| §3, §4, Table 1 | `scripts/*.py` | the credit schemes and scorers (from `katka-readable/scripts/`) |

Reading order: `bits.hpp` → `rlbwt.hpp` → `search.hpp` → `mems.hpp` → `tag_array.hpp` → `twolevel.hpp` → `main.cpp`.

## Experiments

- `bin/katka-blocks PREFIX reads listing [-k 8,12,16] [-L 30] [-n 5000] [--check]` (`experiments/blocks_feasibility.cpp`): a feasibility study of order-k character-level blocks of continuations (a character continues a block if it is the most frequent one after its preceding k-mer). For each k it reports:
  - the number of distinct k-mers;
  - the distribution of block lengths in S;
  - for the listing's MEMs (the reads give their bases), the steps to cross each MEM under the blocks model, against its length;
  - for BML with threshold L run on the same reads, the fraction of its backward steps that fall in failed windows, at the ends of searches, or inside long crossings that block steps would replace.

  The model is at the top of the source file.
- `experiments/stpd_index_chr19.sh STPDIR list.tsv patterns.fa WORK [L]`: the published STPD-index against katka, on the same patterns, with L = 1000 by default.
  - It builds their index through `stpd-prep`/`stpd-export`, with their RLZ text oracle and r-index φ.
  - It runs BML and forward–backward search on their longest-prefix matching (`experiments/stpd_index_mems.cpp`, compiled against their code by `experiments/build_stpd_index_mems.sh`), and their `locate_all_occurrences` for every MEM.
  - It runs katka's plain, two-level and STPD searches (with and without intervals, with two-level intervals, with `--locate`).
  - It compares the MEM sets and the occurrence counts, and writes times, sizes and peak memory to `WORK/summary.txt`.
- `tests/exp_kebab.sh PREFIX queries.fa [L] [k] [R]`: KeBaB then BML, plain BML, and KeBaB inside BML, on queries mostly absent from the text (Section 3.2). It reports the median time of R runs, k-mer tests, window restarts, backward steps and MEMs, and checks that all three give the same MEMs. `tests/gen_novel.py` makes synthetic queries of that kind.

## Testing

`selftest.sh` runs everything below on data from `tests/gen_data.py` (4 species of 5 genomes each, 300 reads with errors and reverse complements, and a chromosome-19-like collection):

- **Unit tests** (`bin/katka-test unit`): the bitvectors, Elias–Fano, the wavelet matrix, `EFRank`, the RMQ and the grammar, on random inputs.
- **Brute force** (`bin/katka-test mems`), using a suffix array of `S$`:
  - BWT characters and ranks;
  - all MEMs of length ≥ L (L = 1, 12, 25, 40) by BML and by forward–backward search, with their intervals;
  - adaptive L (every reported MEM is a MEM, and all the longest ones are reported);
  - the t longest (t = 1, 3, 10; their lengths are the top t);
  - KeBaB pseudo-MEMs and the inline check (exactly the MEMs of length ≥ L, for L ≥ k);
  - counting and listing for every MEM, with full, sampled (s = 4), grammar-compressed and per-species tag arrays;
  - reads plus random patterns from both strands, with and without homopolymer compression;
  - B_λ of the tag runs (`--tag-lcp`) against a suffix array and LCP array (`katka-test blam`);
  - STPD search: the same MEMs and intervals as plain search with the pruned and the full samples, every longest-suffix search ending at the lexicographically smallest occurrence (t = SA[sp]), and the sample's size against the paper's definition of st-colex⁻ (`katka-test stpd-size`);
  - BML over digests: every digested MEM by brute force, with its span, for L = 11, 15, 20, 30 (as KATKA's `RZ_BML_CHECK`), with and without a lookup table.
- **All search configurations** give identical output: RLBWT, explicit and Elias–Fano RLCSA, lookup table, two-level without and with block jumps, with either continuation rule. Two-level search is also checked query by query against plain backward search, on both strands.
- **Command line**: threads; two-level search with tags; gzipped and FASTA reads; a BWT with one string per dataset; each preset.
- **Regressions**, if `KATKA_OLD` (katka-readable's `bin/`), `RZ_OLD` (a build of `katka/rz-index`) or `CHR19_OLD` (`chr19/twolevel/build/`) is set:
  - `tests/regress_katka.sh`: the same text, table, KeBaB filter and listings, byte for byte, as `katka-classify`. It covers counting, listing and range answers, sampled, per-species and grammar tags, `-X`, `-K`, `-D`, `-U`, a per-dataset BWT and 4 threads;
  - `tests/regress_digest.sh`: the same digested text, table and byte map as `rz-prep -B`, and the same listings, byte for byte, as `rz-classify -B` with tags, for L = 11 to 30, with and without a lookup table;
  - `tests/regress_stpd.sh` (with `STPD_OLD`, a built clone of https://github.com/regindex/STPD-index): our pruned STPD sample of rev(T) is, position for position, the paper's st-colex⁻ sampling of T (`stpd-sampling -c`), except for their position 0. `stpd-export`'s four files are byte-identical to `stpd-sampling -c -P`'s, on their yeast text and on a chromosome-like text. Then `experiments/stpd_index_chr19.sh` on the chromosome-like collection must give the same MEMs, by every method of both tools, and the same occurrence counts;
  - `tests/regress_chr19.sh`: the same MEMs, with the same numbers of occurrences, as `query_twolevel`, for BML and FB with L = 100 and 1000.

```
make -C ../katka-readable SDSL=/path/to/sdsl
(cd ../chr19/twolevel && cmake -S . -B build && cmake --build build --target build_twolevel query_twolevel build_index)
make -C ../katka/rz-index SDSL=/path/to/sdsl rz-prep rz-build rz-classify rz-csabuild rz-tagbuild
git clone https://github.com/regindex/STPD-index ../../STPD-index    # builds with CMake; see NOTES.md
KATKA_OLD=../katka-readable/bin RZ_OLD=../katka/rz-index CHR19_OLD=../chr19/twolevel/build STPD_OLD=../../STPD-index ./selftest.sh
```
