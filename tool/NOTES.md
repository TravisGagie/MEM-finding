# Notes on the tool

Written 2026-11-01 for Travis, by the session that built `tool/` from `BRIEF_tool.md`.

## Status

All tested with `selftest.sh` (brute force, plus both regressions):

- **Index and search core (priority 1): done.**
  - RLBWT (wavelet matrix plus Elias–Fano); explicit and Elias–Fano RLCSAs; the t-mer lookup table.
  - Homopolymer compression, with spans mapped back to the read; the low-complexity filter θ.
  - BML with fixed L, adaptive L, or the t longest MEMs; forward–backward search.
  - KeBaB: the Bloom filter, pseudo-MEMs searched longest first, and the new inline window check.
- **Answering each MEM (priority 2): done.**
  - The complete run-length tag array, built by LF and φ without a suffix array.
  - Counting; Muthukrishnan–Sadakane listing, with an RMQ that doesn't store the previous-run array; range/LCA answers.
  - Sampled tags via LF; grammar-compressed tag runs (built-in approximate RePair, or `bigrepair -i` output).
  - Threads, writing in input order; Cliffy's listing format with counts.
- **Two-level indexing with blocks (priority 3): done**, on the text of both strands:
  - closed syncmers canonical under reverse complement, and one dictionary for both strands;
  - B, the parse BWT, blocks with the smallest or the most frequent successor (with ω);
  - it plugs into BML and forward–backward search, and it works with tags. Two-level search with tags gives byte-identical listings to plain search.
- **Minimizer digests (priority 4): done** (added 2026-11-01, after the first push).
  - `prep --digest 4,11` builds symmetric digests: a canonical hash, all tied minimizers per window, equal neighbours collapsed, and the rarest codes (in reverse-complement pairs) sharing one byte. It checks that every genome's digest is strand-symmetric.
  - `find` runs BML on digests with L in bases, the failure jump and spans as in App. A, and keeps only the span-maximal MEMs. The lookup table works over digest symbols (`-F 2`).
  - It is checked against brute force, and against KATKA's original `rz-prep -B` and `rz-classify -B` (identical text, map and listings).
  - Over digests, only BML with a fixed L is supported: no adaptive L, top t, FB, KeBaB, two-level search or homopolymer compression, as in KATKA.
- **STPD search (added 2026-11-02, at Travis's request):** a backend from *Compressing Suffix Trees by Path Decompositions* (Becker et al., ICALP 2026); see the next section.
- **Fixed 2026-11-03: `index` kept the RLCSA in memory through the tag step.** Now each step frees what it built after saving it. Only the RLBWT stays throughout, and the two-level step reloads `PREFIX.csa`, or builds it if it is missing or belongs to another BWT. Each step prints the peak memory so far. On the 64-Mb synthetic collection, `index --csa --tags --by-tag --no-rmq --kebab 20` peaked at 0.11 GB instead of 0.14 GB (the RLCSA's 32 MB), with a byte-identical tag array. On Zymo it should peak near the tag step alone (about 7 GB), not 7 GB plus 5.6 GB. The original report:
  - **(was) Known issue (2026-11-02, desktop): `index` keeps the RLCSA in memory through the tag step.** In `cmd_index` (`src/main.cpp`), `rlcsa X` lives until the function returns, though only the two-level step uses it after `X.save`. So `index --preset zymo` needs the RLCSA (5.6 GB on Zymo) and the tag-run construction (move structure, φ samples, run vectors: several GB more) at once, which exceeds 14 GiB. Fix: free `X` after saving it when `--twolevel` is off (or build the tags first). Until then, build the steps separately: `index --csa`, then `index --tags --by-tag --no-rmq`, then `index --kebab 20` (peaks 6.6, 7.0 and 6.1 GB on Zymo; see `results/zymo/`).
- **The inline KeBaB check is now lazy (2026-11-02, at the paper task's request):** each window's k-mers are tested from its right end leftward, with the rolling codes restarted at the window's last k-mer, stopping at the first failure. Nothing is precomputed over the pattern. `--stats` counts k-mer tests and restarts. The output is unchanged (it jumps to the same place: the first failure from the right is the rightmost one).
- **Not in the tool, as the brief says:** the rz-index, shortcuts, the frequency-rank encoding of the parse, RLZ verification, whole-phrase matching, and the combined LCA-and-listing rule.

Also new: `katka bwt` (libsais) and `katka bwt -d` (one string per dataset, like pfp-merge, for small data). The test data can therefore be built entirely here.

### Added 2026-11-02 for the paper's TODOs

- **STPD with two-level intervals** (`find --stpd --twolevel`): STPD finds the MEMs, and two-level search (with blocks) computes the interval of each reported MEM instead of plain backward search. The output is identical (brute force and command-line tests). On the 2 Mbp × 16 synthetic collection with L = 1000:

  | configuration | time |
  |---|---|
  | STPD with plain intervals | 0.98 s |
  | STPD with two-level intervals | 0.40 s |
  | two-level search alone | 0.34 s |

  Here, recovering each MEM's interval costs about as much as finding it.
- **B_λ** (`index --tags --tag-lcp 0,10,30,100,1000`): the tag boundaries whose suffixes share more than λ characters, with Lemma 11's bound. It is exact against brute force (`katka-test blam`, in the self-test). For the paper's TODO, run it on the Zymo and SILVA references on the desktop, for example:

  ```
  $T index --tags --by-tag --no-rmq --tag-out /dev/null --tag-lcp 0,10,30,100,1000 zy
  ```

  It loads the text into memory for the comparisons, and stops each comparison after max λ + 1 characters. On the synthetic data, per-species tags have B_30 = 29,501 of 267,599 boundaries, against 706,851 of 1,004,033 boundaries with B_10 for per-genome tags.
- **The KeBaB experiment** (`tests/exp_kebab.sh`, Section 3.2's TODO). On 20 synthetic 100-kb queries, each with three 3-kb segments copied from the 64-Mb synthetic collection, with L = 1000 and k = 20 (median of 5 runs, same 71 MEMs):

  | method | time | k-mer tests | backward steps |
  |---|---|---|---|
  | KeBaB then BML | 0.140 s | 1,999,620 | 212,678 |
  | plain BML | 0.100 s | 0 | 277,758 |
  | KeBaB inside BML | 0.110 s | 134,640 (2,093 window restarts) | 212,803 |

  Inside BML it tests 15 times fewer k-mers than KeBaB's separate pass. With an index this small, though, a failed window costs plain BML only about 13 cheap steps, and most of the time goes into the shared segments, so it is not faster. The case the paper describes needs a large index (cache misses per step), i.e. the desktop: for example the chromosome 19 index with queries from `tests/gen_novel.py` built on other sequence.

### The published STPD-index on chromosome 19 (2026-11-03)

`experiments/stpd_index_chr19.sh` compares the published code (github.com/regindex/STPD-index, Becker et al.) with katka on the chr19 question: both strands, L = 1000 (or 100), the same patterns.

**Their index.** Their `stpd-build` first runs `stpd-sampling`, which builds the suffix array, the LCP array and an explicit suffix tree: about 68 bytes per character on their yeast text (0.83 GB for 12.2 million characters). Their README says it has been tested only on files of a few gigabytes. Chromosome 19 with 16 individuals on both strands is 4.07 × 10⁹ characters, so it can't run on the desktop. But `stpd-build` skips that step when `T.colex_m` exists, and builds its index by reading `T.rbwt`, `T.pa` and `T.lcs` sequentially. katka writes all four in O(r) memory plus the text:
- `katka stpd-prep list.tsv st`: their text `st.T` (both strands, split at runs of N into A/C/G/T pieces, concatenated: their oracle accepts only A/C/G/T), `st.S` = reverse(`st.T`), the piece starts, and `st.reflen`;
- Big-BWT on `st.S`, then `katka index --stpd st`: our sample of S is their st-colex⁻ sampling of T;
- `katka stpd-export st`: `st.T.colex_m`, `.rbwt` (our BWT file of S$), `.pa` (SA of S$ streamed by φ) and `.lcs` (PLCP by φ in text order, kept as a 2N-bit vector).

On their yeast text and on a chromosome-like text, these files are **byte-identical** to `stpd-sampling -c -P`'s (`tests/regress_stpd.sh`). Given them, `stpd-build` produces a byte-identical index. The script then runs `stpd-build -i st.T -o st.stpdi -l $(cat st.reflen)`, with their RLZ oracle (the first sequence, both strands, as the reference) and their r-index φ.

**Their MEMs.** Their tools locate whole patterns; MEM finding needs only the longest prefix of a string that occurs, which their Algorithm 1 computes on the way. `experiments/stpd_index_mems.cpp` runs BML and forward–backward search on that primitive over their index of both strands, like katka's, and locates every MEM with their `locate_all_occurrences` (φ⁻¹, exponential search, LCS checks). It is compiled against their code with their flags by `experiments/build_stpd_index_mems.sh`.

One addition to their algorithm was necessary: their binary search returns a candidate sample without checking it, and they check only once, at the end of a whole pattern. For MEMs, the step that fails at every MEM's end would otherwise continue from a wrong sample. So the driver checks each candidate with their `O.LCS`; when the check fails, the extended prefix doesn't occur. The number of rejected candidates is reported (2,265 of 21,398 jumps on the test data, L = 100).

Their text has no separators, so a match can run across the junction between two pieces; the script flags MEMs whose first occurrence crosses one.

**katka's locating.** `find --stpd --locate` lists every occurrence of each MEM by φ⁻¹ from the first, checking the text at the end of each doubling batch, then binary-searching the batch that overshot, as in their paper. φ⁻¹ adds 2 numbers per BWT run to `PREFIX.stpd`. Tested against the suffix array: 347 million occurrences in the self-test, all identical. Counts equal ep − sp for every MEM.

**Tested here** (the self-test, with `STPD_OLD` set): on the chromosome-like collection with L = 100, their BML, their forward–backward search, and katka's plain, two-level, STPD, STPD with two-level intervals and STPD with locate searches all give the same 556 MEMs. Their occurrence counts and katka's equal ep − sp.

**For the desktop:**

1. Build their code (pinned; its CMake fetches sdsl, malloc_count and gsacak, and needs `chromium.googlesource.com`). The first `make` fails with "divsufsort.h: No such file or directory"; reconfiguring with the generated include directory fixes it:

   ```
   git clone https://github.com/regindex/STPD-index ~/rz/STPD-index && cd ~/rz/STPD-index && git checkout 8a432a4
   mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && (make -j4 || true)
   cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-I$PWD/_deps/sdsl-build/external/libdivsufsort/include" && make -j4
   ```

2. Add the STPD sample (with φ⁻¹) to katka's chr19 index of the same collection, e.g. 16 individuals (`c33` from "Chromosome 19" above). It rebuilds the RLBWT too, which is quick:

   ```
   cd ~/tool/chr19 && $T index --stpd c33
   ```

3. Run the comparison (single-threaded, pinned; L = 1000, then 100):

   ```
   cd ~/tool/chr19
   BIGBWT=~/rz/MEM-finding/katka/deps/Big-BWT/bigbwt OURS=$PWD/c33 PIN=0 \
     ~/rz/MEM-finding/tool/experiments/stpd_index_chr19.sh ~/rz/STPD-index list33.tsv queries.fa stpd33 1000
   OURS=$PWD/c33 PIN=0 ~/rz/MEM-finding/tool/experiments/stpd_index_chr19.sh ~/rz/STPD-index list33.tsv queries.fa stpd33 100
   ```

   The second run reuses their index in `stpd33`. Results are in `stpd33/summary.txt`, which overwrites itself, so copy it between runs.

**Resources.**
- **16 individuals (N ≈ 4.07 × 10⁹):**
  - Disk: about 60 GB in `stpd33` while building (`st.T`, `st.S` and the BWT, 4 GB each; `st.T.pa` and `st.T.lcs`, 20 GB each, which `stpd-build` deletes afterwards).
  - Memory: `katka index --stpd` and `stpd-export` hold the text (4 GB) plus O(r) structures (a few GB); their `stpd-build` holds the text for RLZ (4 GB), a suffix array of the reference (about 1 GB) and O(r) for φ. Every step's peak is recorded in the summary.
- **32 individuals:** twice that: the 8 GB text plus the rest may not fit in 16 GB.
- **64 individuals:** not possible on the desktop: both katka's STPD steps and their RLZ build hold the 16 GB text.
- **Queries:** katka `--stpd` holds the text (4 GB for 16 individuals); theirs holds their RLZ oracle.

### Order-k character-level blocks: feasibility (`katka-blocks`, 2026-11-03)

`bin/katka-blocks` (`experiments/blocks_feasibility.cpp`) measures whether blocks of continuations at the character level would pay on SILVA, before anyone builds them. The model (also at the top of the source):

- **Continuations:** a character continues a block if its preceding k-mer is over A/C/G/T and it is that k-mer's most frequent next character in S.
- **Crossing a match:** with block starts s0 < s1 < … < sz whose contexts lie inside the match, crossing X[a..e) costs z block steps and (e − a) − (sz − s1) character steps on T's BWT. The character steps cover the right end, the back-up (the leftmost complete block redone, as in the back-up lemma) and the left end. With fewer than two complete blocks it costs e − a character steps, as now.

It reports, for each k:
- (a) distinct k-mers;
- (b) block lengths in S;
- (c) the steps to cross the listing's MEMs, against their lengths;
- (d) BML's backward steps (L = 30, run by the tool itself on the same reads, over PREFIX.rlbwt or PREFIX.csa), split into failed windows, ends of searches, and inside long crossings. The split accounts for every step BML takes.

**Tested here** (in `selftest.sh`): the hash table against direct counts, and the block count and distinct k-mers against a Python brute force. The program checks that BML's MEMs equal the listing's.

On the synthetic reads (40-bp MEMs, k = 8), crossings would cost 0.68 steps per character and BML 1.14 times fewer steps, with 54% of its steps in failed windows. With k = 12 and 16, these short MEMs contain too few block starts to save anything, although blocks in the text are long. On the chromosome-like queries (k = 8), BML would take 2.0 times fewer steps.

**For the desktop (SILVA):** it needs the tool's SILVA index (`bac.S`, `bac.rlbwt`, and `bac.csa` for speed), a few thousand reads, and the listing of `katka find -L 30` on exactly those reads:

```
cd ~/tool/silva                         # where katka index --preset silva built bac.*
head -n 20000 READS.fq > blocks_reads.fq      # 5,000 reads (FASTQ: 4 lines each), e.g. one region's mate-1 reads
$T find -L 30 -c --csa -F 10 bac blocks_reads.fq blocks_L30.listings
~/rz/MEM-finding/tool/bin/katka-blocks bac blocks_reads.fq blocks_L30.listings -k 8,12,16 -L 30 > blocks_silva.txt 2>&1
```

For several regions, run it once per region, each with its own reads and listing. Memory:
- the text (n bytes, about 1.2 GB);
- the index, during the BML pass only;
- for k = 16, a hash table of about 12 bytes × (4/3) × the number of distinct 16-mers, which briefly needs 1.5 times that while it doubles. For a few hundred million distinct 16-mers, that is up to about 10 GB.

If k = 16 runs out of memory, run `-k 8,12` and report that. k = 8 and 12 use arrays of 1 MB and 268 MB. It checks that BML's MEMs equal the listing's, and says so if the listing was made with other options.

### STPD search

`katka index --stpd` and `katka find --stpd` (`src/stpd.hpp`). BML and forward–backward search only need longest-suffix searches, so this is the mirror image of the paper's Algorithm 1.

- **The search:** it matches right to left with π = ISA, keeping the occurrence of the matched suffix whose suffix of S is smallest (t = SA[sp]). If S[t−1] is the next pattern character, t decreases by one. Otherwise a binary search over the sample, comparing the text directly, finds the smallest suffix prefixed by the extended string.
- **The sample:** it consists of BWT run starts x (as text positions SA[x] − 1), sorted by LF(x). A run start is kept only if lcp(SA[y], SA[x]) < lcp(SA[x−1], SA[x]), where y is the previous row with the same character. The sample is built from one LF pass, as for the tag array, plus text comparisons: no suffix array or suffix tree, unlike the paper's builder (which computes the explicit suffix tree).
- **Size:** on every test collection, the pruned sample is exactly 2 smaller than |st-colex⁻| of the reversed text computed by the paper's definition (`katka-test stpd-size`). Those 2 are boundary entries, for row 0 and for text position 0, which a search never jumps to. It holds about 0.66 r positions on the test collections.
- **Intervals:** the search gives positions, not BWT intervals. The interval of each reported MEM is computed afterwards by backward search over the MEM, so tags and the mems format are unchanged; `--no-intervals` skips this.
- **Memory:** the text S$ is held in memory as the oracle (1 byte per character), plus the sample. A compressed oracle (RLZ or 2-bit packing) is not implemented.

On the 2 Mbp × 16 chromosome-like collection, BML with L = 1000 (same MEMs in every case):

| configuration | time |
|---|---|
| plain search | 2.1 s |
| STPD, with intervals | 1.2 s |
| STPD, `--no-intervals` | 0.20 s |
| two-level with blocks | 0.40 s |

The sample takes 8.2 MB (2.5 million positions, against r = 3.8 million) on top of the 64 MB text.

**Against the paper's code** (https://github.com/regindex/STPD-index, `stpd-sampling -c`, `tests/regress_stpd.sh`): our pruned sample of rev(T) is exactly their st-colex⁻ sampling of T, position for position. The only difference is their position 0 (the first prefix), which a search never needs. This holds on their `data/yeast.txt` (12.2 million characters: 5,479,329 positions, 0.654 r) and on a chromosome-like text. Building the sample took 4.7 s here (BWT 0.7 s, then the LF pass and comparisons), against 28 s for their `stpd-sampling`, which builds the suffix tree.

To build their code, allow `chromium.googlesource.com` (its CMake fetches dependencies there), and add the generated libdivsufsort include directory to the compiler flags:

```
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-I$PWD/_deps/sdsl-build/external/libdivsufsort/include"
```

The first `cmake` run creates that directory. The plain `cmake ..` failed with "divsufsort.h: No such file or directory".

### sdsl-lite

The tool has no sdsl dependency, only zlib and the bundled libsais (the auto-mode permission check refused to build sdsl at first). Its own versions replace sdsl's `sd_vector` (Elias–Fano), `wt_huff` (a wavelet matrix), `rank_support_v5` and `rmq_succinct_sct` (`src/rmq.hpp`, a min-tree in balanced parentheses). They are checked against brute force in `katka-test unit`.

sdsl was then built here, only for the regressions. Its consequences:

- The tool's file formats are new (only the KeBaB `.kbf` format is KATKA's), so the old `.rz`, `.tag` and `.csa` files can't be reused. Rebuilding is fast (see below).
- On 12,000 short reads here, against `katka-classify` (sdsl):
  - the RLBWT backend is about 10% faster (3.8 s against 4.3 s);
  - the explicit RLCSA, which is the same code, takes about the same time or up to 15% longer (1.2–1.5 s against 1.2–1.3 s; noisy).

  Run-to-run noise here is about 10%. Please time it on the desktop.

## Where both strands (and other changes) make the tool differ from the old code

1. **The chr19 text.** The old code indexes T, the forward strand, and the reversed collection T′. The tool indexes S = every sequence and its reverse complement, and does all searches on it. So:
   - **MEMs are with respect to both strands.** Matching the reverse strand can make a MEM longer or add MEMs, as for ropebwt3 and Movi 2. The paper's comparison section notes these differences.
   - **A MEM's interval** is in the BWT of S. For BML and FB it is the interval of the MEM's reverse complement (the forward search ends there), as in KATKA. Its size counts occurrences on both strands.
   - **The index covers twice as many characters.** n = 4.07 × 10⁹ for 16 individuals, just under 2³²; with 32 and 64 individuals the explicit RLCSA switches to 64-bit lists.
   - **Measured sizes.** On the 2 Mbp × 16 synthetic collection:
     - the RLBWT of S takes 8.7 MB, against 16.0 MB for the old RLBWTs of T and T′. The encodings differ, so this isn't a like-for-like comparison;
     - the parse index with blocks takes 9.6 MB, against 6.9 MB for the old two levels. The parse of both strands has 3.2 million phrases and 456,000 distinct phrases.

     The real sizes need measuring on the desktop.

   On the synthetic chromosome, `tests/regress_chr19.sh` gets the same MEMs and counts as `query_twolevel` for BML and FB with L = 100 and 1000, because random sequence has no long reverse-strand matches. On real chromosome 19, expect the same small differences the paper reports for ropebwt3 (e.g., 234 reverse-strand-only MEMs with L = 100).
2. **The syncmer order g (decision 2).** g(γ) is now a hash (splitmix64) of min(code(γ), code(rc(γ))), not a 32-bit hash of min(code(γ), code(γ^rev)). I checked the two claims in the brief:
   - **Prefix-freeness and the order lemma (Lemma 6) still hold.** The proof uses only that whether a k-mer is a closed syncmer doesn't depend on its context, which holds for any order g on s-mers. Nothing in it uses symmetry.
   - **One dictionary serves both strands.** The s-mers of rc(κ) are the reverse complements of κ's, in reverse order, with the same values of g, so κ is a closed syncmer if and only if rc(κ) is, with first and last swapped. So the syncmer starts of rc(P) are |P| − k − q for the starts q of P, and its phrases are the reverse complements of P's, in reverse order. The dictionary stores, for each phrase of S over A/C/G/T, its id and the id of its reverse complement, so one parse of a pattern serves backward search on P and forward search, which is backward search on rc(P). The build checks that every such phrase has its reverse complement in S. The only exception is the first phrase of S, which starts at position 0 rather than at a syncmer, so no pattern phrase can equal it. The builder dies if more than two are missing. The tests compare two-level against plain search query by query on P and on rc(P).

   Phrases containing `X` or the terminator are left out of the dictionary, since no pattern phrase can contain them. The old dictionary held all phrases.
3. **BML variant.** The tool uses KATKA's BML everywhere: after a window matches, it extends left to a, then forward to e, and continues at x = e − L + 1. The chr19 code extends forward from x, then backward from e, and continues at the start of the longest suffix of P[0..e+1). Both find exactly the MEMs of length at least L, but the step counts differ. On a 2 Mbp × 16 synthetic collection, plain BML with L = 1000 took 7.3 million character steps in the tool against 11.5 million in the old code. With two-level search, extending left after the window is done by searching again from the window's end (`segment_search::extend_left`), which costs extra parse steps on the window.
4. **Two-level indexing needs one terminator**, as before. Zymo's per-dataset BWT can't be used with `--twolevel`; `katka index` says so.
5. **N.** The old chr19 code indexed N as a character; S has `X` wherever the input has anything other than A, C, G or T. Patterns are split at non-ACGT characters in both, so MEMs over A/C/G/T are unaffected.
6. **The t longest MEMs.** When several MEMs tie for t-th place, the tool keeps the one found first.
7. **`-K` with `-X`.** KATKA required the filter's k to be at most L. The tool requires it only with a fixed L (pseudo-MEMs), since adaptive L starts low and the longest-first order handles it, as KATKA's `-X -L 20 -K` runs did.
8. **Tags** may be any 32-bit integers (KATKA stored 14 bits).
9. **Grammar.** The built-in RePair is approximate: each round replaces every pair whose frequency is at least half the largest. The grammar differs from bigrepair's, so the size reported in the paper (0.23 GB for SILVA) should be remeasured, or the tool given bigrepair's output with `--bigrepair`. Answers are identical either way (tested). The built-in one needs about 16 bytes per run, roughly 6 GB for SILVA's 359 million runs, so on the desktop bigrepair may be the safer choice.
10. **Digests of reads collapse equal neighbouring codes, as the reference's do.** KATKA collapsed equal neighbouring bytes in reads but codes in the reference. These differ only when two different codes that share the merged byte are adjacent. In that case KATKA's read digest could miss a match, which the tool's cannot. The test data includes merged codes, and the listings were identical.
11. **Output order with `--top`**: the t longest are printed at the end of the read, sorted by position.

## What was tested here (small data only)

`./selftest.sh` takes about 2.5 minutes here. It runs:

- the unit tests;
- the brute-force checks, on 340 patterns (300 reads plus 40 random patterns from both strands). These cover about 260,000 MEMs from plain BML and FB, every search configuration, KeBaB in both modes, adaptive L, the t longest, and counting and listing with four kinds of tag array. They also run on the homopolymer-compressed text, and on the chromosome-like collection with both continuation rules;
- command-line checks of each preset and option combination;
- the three regressions, with `KATKA_OLD`, `RZ_OLD` and `CHR19_OLD` set. All pass:
  - **KATKA:** byte-identical listings for 12 option combinations on each of three texts (plain, homopolymer-compressed, and one BWT string per dataset);
  - **digests:** the same digested text, byte map and listings as KATKA's original `rz-classify -B`, for L = 11 to 30 (`tests/regress_digest.sh`);
  - **chr19:** identical MEMs and counts.

KATKA's 40-genome `test-small` isn't in the repository, so the KATKA regression ran on synthetic data from `tests/gen_data.py`. To run it on `test-small`, put `list.tsv` and `reads.fq` in a directory DIR and run `tests/regress_katka.sh OLD_BIN DIR WORKDIR`.

Scale check here: on a chromosome-like collection of 16 mutated copies of 2 Mbp (64 million characters with both strands), building took 6 s for the BWT and 9 s for the rest. BML with L = 1000 on two 2-Mbp queries took 1.78 s with plain search, 0.39 s two-level and 0.38 s with blocks, with identical output. The old code took 6.87 s, 0.77 s and 0.72 s on the same data, with T and T′.

## Commands for the desktop

Build (no sdsl needed) and test:

```
cd MEM-finding/tool && make && ./selftest.sh
T=$PWD/bin/katka
```

Optionally, also run the regressions (needs sdsl-lite, as for the old code):

```
make -C ../katka-readable SDSL=../katka/deps/sdsl
(cd ../chr19/twolevel && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j4 --target build_twolevel query_twolevel build_index)
KATKA_OLD=$PWD/../katka-readable/bin CHR19_OLD=$PWD/../chr19/twolevel/build ./selftest.sh
```

Big-BWT and pfp-merge are where `katka/setup_deps.sh` put them: `katka/deps/Big-BWT/bigbwt` and `katka/pfp-merge/`.

### SILVA

The reference is the one built by the first step of `katka/silva/build_katka.sh`: `katka/silva/work/ref/dna/doc_<i>_seq.fa`, one file per genus, in tree order.

```
S=$HOME/MEM-finding/katka/silva; mkdir -p ~/tool/silva && cd ~/tool/silva
n=$(cut -d' ' -f1 $S/work/ref/done)
for i in $(seq 1 $n); do printf 'silva\t%s\n' "$S/work/ref/dna/doc_${i}_seq.fa"; done > list.tsv
$T prep list.tsv bac
~/MEM-finding/katka/deps/Big-BWT/bigbwt bac.S && mv bac.S.bwt bac.bwt    # or: $T bwt bac (libsais, about 6 bytes per character)
$T index --preset silva bac                  # explicit RLCSA, tags (genus = genome index) with the RMQ
$T index --tags --sample 4 --tag-out bac.s4.tag bac
$T index --tags --grammar --tag-out bac.g.tag bac   # or --bigrepair PREFIX, for bigrepair's grammar
```

Classification uses the same command lines as before, with `katka find` in place of `katka-classify`. For example, `katka-classify -j 8 -L 30 -c -C bac.csa -T bac.s1.tag -F 10 bac.rlbwt R out` becomes:

```
$T find --preset silva -j 8 bac R out        # = -L 30 -c --csa -F 10
$T find -j 8 -L 30 -l --csa -F 10 bac R out  # even credit (listing)
```

The digested index of App. B (k = 4, w = 11), with an explicit RLCSA and tags:

```
$T prep --digest 4,11 list.tsv dg
~/MEM-finding/katka/deps/Big-BWT/bigbwt dg.S && mv dg.S.bwt dg.bwt
$T index --csa --efcsa --tags dg
$T find -j 8 -L 30 -c --csa -F 2 dg R out      # L = 15 ... 50 for the table in App. B; --efcsa for its Elias-Fano row
```

Then score with `scripts/score_silva.py` (or the `compare_lib.sh` pipeline, giving it the listings) as before.

### Zymo

Use the list and data `katka-readable/README.md` describes: one line per genome FASTA, tagged by species.

```
mkdir -p ~/tool/zymo && cd ~/tool/zymo
$T prep --preset zymo list.tsv zy            # -d -s -U
for f in $(cat zy.datasets); do ~/MEM-finding/katka/deps/Big-BWT/bigbwt $f; done
(cd ~/MEM-finding/katka/pfp-merge && ./PFPmerge.py -p 4 -o ~/tool/zymo/zy.bwt $(cat ~/tool/zymo/zy.datasets))
$T index --preset zymo zy                    # explicit RLCSA, species tags without the RMQ, KeBaB k = 20
$T find --preset zymo -j 6 zy reads.fa longest.listings         # -U -X -L 20 -K 20 -D 14 -c --csa
$T find -U -L 30 -D 14 -c --csa -j 6 zy reads.fa L30.listings     # threshold L = 30
python3 scripts/credit_schemes.py L30.listings truth.txt zy.tbl
python3 scripts/maxmem_roc.py longest.listings truth.txt zy.tbl reads.fa
```

The tag build streams the suffix array once by φ. KATKA's construction took 7.2 hours on the laptop for this index. This one keeps the runs in memory: about 12 bytes per run, roughly 5.4 GB for 449 million runs, plus about 10 bytes per BWT run for the samples.

### Chromosome 19

`$H` is the directory `chr19/twolevel/run_final.sh` uses (`~/shortcuts/chr19_scaling/hprc_chr19`), with `collection.txt` (one sequence per line, CHM13 first, 129 lines for 64 individuals) and `queries.txt`. Take the first 33, 65 or 129 lines for 16, 32 or 64 individuals:

```
H=$HOME/shortcuts/chr19_scaling/hprc_chr19; mkdir -p ~/tool/chr19 && cd ~/tool/chr19
for n in 33 65 129; do
  mkdir -p fa$n; head -n $n $H/collection.txt | awk -v d=fa$n '{f = sprintf("%s/s%03d.fa", d, NR); print ">s" NR > f; print > f; close(f)}'
  ls fa$n/*.fa | sed 's/^/-\t/' > list$n.tsv
  $T prep list$n.tsv c$n
  ~/MEM-finding/katka/deps/Big-BWT/bigbwt c$n.S && mv c$n.S.bwt c$n.bwt
  $T index --preset chr19 c$n                # RLBWT, explicit RLCSA, two-level with blocks
done
tr -c 'ACGT\n' '\n' < $H/queries.txt | awk 'length($0) > 0 { print ">p" NR - 1; print }' > queries.fa
for n in 33 65 129; do
  for L in 1000 100; do
    taskset -c 0 $T find --preset chr19 -L $L --stats c$n queries.fa bml$L.$n.tsv             # two-level, blocks
    taskset -c 0 $T find --preset chr19 -L $L --no-jumps --stats c$n queries.fa nj$L.$n.tsv  # two-level, no blocks
    taskset -c 0 $T find -L $L --format mems --csa --stats c$n queries.fa plain$L.$n.tsv      # plain search
    cmp bml$L.$n.tsv plain$L.$n.tsv && cmp nj$L.$n.tsv plain$L.$n.tsv
  done
  taskset -c 0 $T find --preset chr19 --fb --stats c$n queries.fa fb1000.$n.tsv
done
```

`--continuation frequent` at index time gives the most-frequent-successor blocks. For STPD search, add `--stpd` at index time and use `find -L 1000 --format mems --stpd [--no-intervals]`. It needs the text in memory: 4 GB for 16 individuals (both strands).

Memory: the parse of both strands has about twice as many phrases as T's. With 64 individuals, if that is about 800 million, building blocks needs roughly 13 bytes per phrase, about 10 GB. That should fit in 16 GB, but tightly. If it doesn't, `--no-blocks` needs about 12 bytes per phrase.

## Open questions for Travis

- **test-small:** the 40-genome collection, to repeat the KATKA regression on real genomes. Ideally added to the repository, or attached, as the brief suggests.
- **Preset `novel`:** I chose L = 1000 and k = 20; the brief only says "fixed large L". Change it in `PRESETS` in `src/main.cpp` if you want something else.
- **Presets `zymo` and `chr19`:** the `zymo` preset includes `-D 14` for the longest MEMs, as the brief's table says; the paper's KATKA `-X` command line in `katka-readable/README.md` doesn't. The `chr19` preset has L = 1000; give `-L 100` for the other column.
- **Commits:** the commits on branch `tool` are authored as "Claude"; change the author if you prefer.
