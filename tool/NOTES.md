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

  On the synthetic data, per-species tags have B_30 = 29,501 of 267,599 boundaries, against 706,851 of 1,004,033 boundaries with B_10 for per-genome tags.

  **Memory, rewritten 2026-11-04 so that Zymo fits in 11 GB.** The old code read the whole text into memory (39 GB on Zymo), and the tag-run vectors grew by doubling during the φ pass. That pass alone held about 11 GB on Zymo: the RLBWT, the φ samples, and the vectors. Now (`tag_runs::compute`, `src/tag_build.hpp`):
  - **The tag runs** are streamed to a temporary file during the φ pass. Once the φ samples are freed, they are read back into vectors of exactly their size: 12 bytes per run, 5.4 GB on Zymo.
  - **A text within `--blam-mem`** (default 6 GB; SILVA) is mapped with `mmap`, and the two suffixes at each boundary are compared directly, as before.
  - **A larger text** (Zymo) is never read at random:
    - During the pass, the two text positions of every tag boundary (and how many characters each has before its dataset's end) are written to 256 files by position.
    - Then, a group of λ values at a time, the files are read in position order with the text mapped and read sequentially. Each suffix's first λ + 1 characters get a Karp–Rabin fingerprint (mod 2⁶¹ − 1).
    - A boundary counts if both suffixes have λ + 1 characters before their dataset's end and their fingerprints are equal. A count can only be too high, by a collision, with probability at most (λ + 1)/2⁶¹ per boundary.
    - The fingerprint arrays take 8 bytes per boundary per λ (3.6 GB on Zymo), and the λ values are grouped to fit `--blam-mem`. On Zymo with the default, that means one λ per pass. Each pass reads the 39 GB text and the 10.8 GB of records once, sequentially.
  - **`--tmp DIR`** puts the large temporary files there: the φ samples' spill (about 8 GB on Zymo), the tag runs (5.4 GB) and the boundary records (10.8 GB). The desktop's internal disk has about 24 GB free, so use the external drive.
  - **`--blam-external`** forces the records and fingerprints, for testing.

  **Checks:**
  - Against brute force: the default and the external path, each with one λ per pass and with `--tmp`, on the self-test's reads collection.
  - On a collection with one terminator per dataset (`z`, as Zymo is built), against a new brute force over T₀ $₀ T₁ $₁ ⋯, the string `katka bwt -d` sorts. The old brute force used S$ with one terminator, whose LCPs run across dataset ends, so it disagreed slightly with both the old and the new implementation there (217,323 against 217,319 at λ = 0). The old and new implementations agree with each other and with the new brute force.
  - The tag arrays are byte-identical to the old code's.
  - On the 64 Mbp synthetic collection with per-genome tags (63.7 million boundaries), both paths give the same B_λ and tag array as the old code: 0.89 GB and 62 s direct (old: 0.92 GB, 66 s), 2.19 GB and 183 s external. There, boundaries outnumber the text's bytes.

  **Expected on Zymo:**
  - the φ pass: the RLBWT (1.0 GB) and the φ samples (about 3 GB), after `sample_rindex`'s own peak (its LF structure, about 5–6 GB as before);
  - B_λ: the RLBWT and 3.6 GB per λ;
  - then the tag runs (5.4 GB) and the tag array build, as before (about 7 GB).

  **For the desktop.** Zymo, with the temporary files on the external drive (`/mnt/ext/tmp` stands for a directory there). The φ pass rereads the BWT file (39 GB), and each λ pass the text (39 GB) and the records. Run it as a `systemd-run --user --unit=...` service:

  ```
  mkdir -p /mnt/ext/tmp
  systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "blam: %e s, peak %M KB" \
    ~/rz/MEM-finding/tool/bin/katka index --tags --by-tag --no-rmq --tag-out /dev/null --tag-lcp 0,10,30,100,1000 --tmp /mnt/ext/tmp ~/zymo/hpc/zy > blam.zymo.err 2>&1
  ```

  For SILVA, the same command on its index without `--by-tag` if its tags are per genome, as in `results/silva_profile/`. Its text fits the budget, so it is compared directly. Report the `tag runs rho` and `lambda = ...` lines (B_λ and the bound), the time and peak, and the run count against the cached `.tagruns` (it must be equal). Put the results in `tool/results/blam/`.
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

**Their forward–backward search** (2026-11-03): the desktop stopped it after 54 minutes on chromosome 19, while their BML took 45 s. It runs correctly; the slowness is in their code, not in the driver.
- **Where the time goes:** on the 2 Mbp × 16 synthetic collection (L = 1000; the same 1,996 MEMs as BML and as katka), their BML takes 11.6 s and their FB 74.6 s. Of that, 99.6% is in their `binary_search_with_tabluation`: the driver's check takes 0.1 s, their LCP extensions 0.09 s, and the string copies are negligible.
- **Why per jump:** each probe of their binary search compares the pattern with a sample backwards from the end, over the whole matched prefix (`LCS_char`). So a jump costs O(matched length × log samples), about 110 µs per jump here, for BML and FB alike.
- **Why FB is slower:** FB makes 7.8 times as many longest-prefix searches as BML (77,720 against 9,949) and 7 times as many jumps (709,296 against 101,433). Each search starts again from the left and re-crosses the overlapping parts of consecutive MEMs.
- **Why chromosome 19 is worse:** MEMs are much longer there, so each jump costs more as well.
- **katka:** its STPD search keeps the common prefix shared with each bound of the binary search (the lcp-skip), so a probe doesn't recompare what is known to match. Its FB takes 0.61 s on the same data (729,026 jumps, about 0.8 µs each), and its BML 0.21 s.

**katka's locating.** `find --stpd --locate` (since 2026-11-03 it computes no intervals unless tags are asked for: the desktop found 132 of its 135 s going to plain-search intervals) lists every occurrence of each MEM by φ⁻¹ from the first, checking the text at the end of each doubling batch, then binary-searching the batch that overshot, as in their paper. φ⁻¹ adds 2 numbers per BWT run to `PREFIX.stpd`. Tested against the suffix array: 347 million occurrences in the self-test, all identical. Counts equal ep − sp for every MEM.

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

**Fixed 2026-11-03:** listing records are now paired with reads by position (with a check that the names agree), not by name. On SILVA, names repeat across regions and mates, and the desktop had to make them unique (`results/silva_blocks/README.md`).

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

### STPD with the text read through the parse (`--parse-oracle`, 2026-11-03)

`katka index --parse-oracle [--syncmer k,s] PREFIX` writes `PREFIX.ptxt` (`src/parse_text.hpp`). It contains the two-level index's parse of S$ (the same closed syncmers, the same lexicographic phrase ids) stored as a text oracle:
- the phrase starts: the boundary list, in Elias–Fano;
- the parse: one packed id per phrase;
- the phrase dictionary: the distinct phrases concatenated in id order, with packed offsets.

Run alone, it reads only `PREFIX.S`, not the BWT. `katka find --stpd --parse-oracle` loads it instead of `PREFIX.S`. Every text read of STPD search (the comparisons of the binary search, the S[t−1] test of the LF-like step, and the checks of φ⁻¹ locating) goes through it:
- **Position to phrase:** a predecessor query on the starts, then the phrase's characters from the dictionary. Consecutive phrases overlap by k characters, so either phrase can serve a position in the overlap.
- **The current phrase:** each thread keeps a cursor on its current phrase. A read inside it costs no query. A read just past either end moves to the next or previous phrase using only the phrase lengths (the next phrase starts k characters before the current one ends). Only a read anywhere else costs a predecessor query.
- `--stats` prints how the seeks were served (same phrase, neighbour, query).

**Tests:** `katka-test mems --stpd` builds an oracle with k,s = 12,3 (short phrases, so many boundary crossings). It checks every position read forward and backward, plus 20,000 random reads and LCEs, against the text. It also runs the STPD MEM checks and the φ⁻¹ locating checks (against the suffix array) through the oracle. `selftest.sh` compares `find --stpd` and `find --stpd --locate` outputs with and without the oracle, byte for byte.

**On the synthetic collection** (16 × 2 Mbp, both strands, n = 64.0 × 10⁶; 2 patterns of 4.0 Mbp; search time, the minimum of 3 runs, since this container's timings vary by up to 30%). The outputs were byte-identical to the plain text in every case, with the same numbers of jumps and comparisons:

| | plain text | oracle 51,12 | oracle 31,8 |
|---|---|---|---|
| size | 64.0 MB | 45.0 MB (5.62 bits/char) | 43.9 MB (5.71 bits/char) |
| (dictionary / parse / starts) | | 33.1 / 7.6 / 2.8 MB | 26.2 / 13.4 / 4.3 MB |
| L = 1000, `--no-intervals` | 0.20 s | 0.84 s | 1.05 s |
| L = 1000, `--locate` | 0.20 s | 0.89 s | 1.10 s |
| L = 1000, FB, `--no-intervals` | 0.74 s | 4.64 s | 4.80 s |
| L = 100, `--no-intervals` | 0.29 s | 1.53 s | 1.86 s |
| L = 100, `--locate` | 0.33 s | 1.64 s | 1.83 s |

At L = 1000 (k,s = 51,12) the seeks split into 8.36 million served by the current phrase, 1.01 million by a neighbour and 1.88 million predecessor queries. The queries are almost exactly the binary search's comparisons (1.81 million): each comparison starts at a sampled position anywhere in the text. The LF-like steps and the comparisons' rightward reads almost never need a query. So the slowdown (4 to 6.5 times) is the cost of a random seek: about ten dependent memory accesses (two Elias–Fano selects, the packed id, two dictionary offsets, the bytes) against one for the plain text. On this collection the dictionary is half the text, because the synthetic genomes differ more than human haplotypes. On chr19 the dictionary should be a small fraction, and the parse and starts should take about (log σ + 2 + log(n/m)) bits per phrase, roughly 30 bits per 18 characters. That is 1.7 bits per character against 8.

A cheap fix for the random seeks, not implemented: store each sampled position as (phrase index, offset) instead of a text position. That costs about the same bits, and a comparison would then start without a query. Only φ⁻¹ positions would still need one. That is also the first step of STPD over the parse.

**For the desktop:** run `experiments/parse_oracle_chr19.sh` on the 16-individual index, which has `c33.S` and `c33.stpd`, with the chromosome-19 preset's syncmers and one smaller setting:

```
cd ~/rz/MEM-finding && git pull origin tool && make -C tool
PIN=9 ~/rz/MEM-finding/tool/experiments/parse_oracle_chr19.sh ~/tool/chr19/c33 ~/tool/chr19/queries.fa ~/tool/chr19/poracle33 51,12 31,8
```

It needs `c33.rlbwt` too (for `find`), but not `c33.bwt`. Each oracle is built in `poracle33/k51s12/` and `poracle33/k31s8/` next to links to c33's files, so c33 is not changed. Single-threaded, it runs `--stpd --no-intervals` and `--stpd --locate` at L = 1000 and 100, and FB at L = 1000, each with the text and with each oracle. It cmps every oracle output against the text's and writes `poracle33/summary.txt` (sizes, times, peaks, seek counts). The plain-text runs hold the 4.07 GB text. The oracle runs should take about 1 GB. Copy `summary.txt` and the `*.log` files to `tool/results/chr19/poracle33/` and push.

### STPD over the parse (sketch; implemented below)

The idea is to run STPD search over the parse Z (phrase ids, the oracle's `ids`) for the whole phrases of the pattern, and over characters only at its two ends. The setting, as in two-level search: a pattern piece X = X_l · P_1 ⋯ P_h · X_r, where:
- X_r starts at X's last closed syncmer, so it is a prefix of a phrase;
- P_1 ⋯ P_h are the whole phrases between X's first and last syncmers (consecutive ones overlapping by k);
- X_l is the part before X's first syncmer.

Three facts make this work:
1. A closed syncmer depends only on its k-mer. So every occurrence of X[f..] (f = X's first syncmer) in S starts at a phrase start and has the same phrase boundaries as X. (A syncmer at position 0 is the start of phrase 0 anyway.)
2. The dictionary is prefix-free (the build checks this), and the ids are lexicographic. So for suffixes of S that start at phrase starts, the character order equals the order of their parse suffixes compared as id strings: equal ids mean equal phrases, and the first differing id decides.
3. Hence the smallest suffix of S prefixed by X[f..] starts at start(u), where Z[u..] is the smallest parse suffix that matches P_1 ⋯ P_h as ids and then X_r as a character prefix of phrase Z[u + h].

**The per-phrase table.** F[a], for each id a, is the smallest parse suffix starting with a, i.e. SA_Z at the first row of a's block. It takes σ log m bits. Optionally, G[a] is the range of the parse sample's entries whose suffix starts with a, which takes σ log r_Z bits. F[a] is always a parse-sample entry, because a's first occurrence in the parse BWT is a run start. So F only saves a binary search, and G shortens the others to the entries starting with a.

**Entering from the character-level interval.** The right end X_r is matched by characters: character-level STPD (or backward search) gives the smallest occurrence t, or the interval [sp, ep). Then:
- From t: t is a phrase start by fact 1, and u = pred(t) on the boundary list (check start(u) = t).
- From [sp, ep), without a toehold: ρ = rank_B(sp) is the first parse row prefixed by X_r. The ids whose phrases have prefix X_r are contiguous (lexicographic ids), and ρ is the first row of the smallest of them, a₁ (ρ = C_Z[a₁]). So a₁ comes from C_Z, and u = F[a₁].

Without either: a₁ is the smallest id whose phrase has prefix X_r, found by a binary search over the dictionary, which is sorted.

**Parse-level steps.** These go right to left over P_h, …, P_1, with state u = the smallest parse suffix prefixed by the pattern's ids matched so far plus X_r. To prepend id a:
- If Z[u − 1] = a, then u − 1. This is one read of `ids`, next to the cursor.
- Otherwise, binary search over the parse sample for the smallest parse suffix prefixed by a·(ids)·X_r. The parse sample is the STPD sample of Z: run starts of the parse BWT, pruned with the LCP of Z in phrases, built by the same LF pass over the parse BWT, which the `.lvl` already has. With G it searches only among the entries starting with a. Each comparison reads ids of Z sequentially from a phrase index (no predecessor query), using the lcp-skip in ids. Only the last phrase's comparison reads characters from the dictionary, for X_r.
- If no parse suffix is prefixed by it, the match ends inside phrase a. Convert u to t = start(u) (one select) and continue by characters into a's characters, right to left from its last uncovered one, with character-level STPD through the same oracle. That gives the exact match start for BML, and the next search restarts there.

A parse-level jump replaces up to |P_i| − k character-level jumps, each of which would cost log r comparisons with random seeks. The parse sample has at most r_Z entries (the parse BWT runs: 514 thousand for 3.2 million phrases on the synthetic data, against r = 3.8 million).

**Finishing at the left end.** After P_1, t = start(u) (one Elias–Fano select) is the smallest suffix of S prefixed by X[f..] (fact 3). Character-level STPD continues over X_l, testing S[t − 1] through the oracle: the cursor is already on phrase u − 1, whose last characters precede start(u), so this costs no query. When it jumps, it uses the character-level sample. The final t is the smallest occurrence of the whole match, so locating (φ⁻¹ with oracle checks) and the BML/FB drivers are unchanged.

**What it needs on top of the oracle:**
- F (and optionally G);
- the parse sample, as phrase indices;
- C_Z (the parse BWT already has its counts);
- pattern parsing (`parse_pattern`, as in two-level search), plus a map from each pattern phrase to its id. That is the dictionary's hash table, `two_level::dict`; a phrase that is not in the dictionary ends the parse-level match.

Open points: an id-level pruning rule equivalent to the character-level one (the same proof should go through with Z's LCP in phrases); and how often Z[u − 1] ≠ a on human haplotypes, i.e. whether parse-level jumps are rare enough. `--stats` on the two-level index already counts parse steps, which bounds this.

### KeBaB: `--kebab-check T t` and `--kebab-top t --fixed-L` (2026-11-03)

- **`--kebab-check T t`** (T ≥ t ≥ 1, with `-K k`):
  - Search the T longest pseudo-MEMs, longest first, with the t-longest selection with ties (`--longest t`). L is carried across them: it rises to the t-th longest MEM so far, and inside pseudo-MEMs it is at least k.
  - Accept if there are at least T + 1 pseudo-MEMs and at least t MEMs are at least as long as the (T+1)st pseudo-MEM. With `--strict`, they must be longer.
  - Otherwise, run BML on the whole read with the same selection from L. That gives the t longest MEMs, or all of them if there are fewer than t.
  - *Why it is exact:* every MEM not yet found lies in a pseudo-MEM beyond the T-th, so it is at most λ = |p_{T+1}| long.
    - Without `--strict`, accepting therefore gives the right lengths for the t longest. A MEM of length exactly λ that ties the t-th place may be left out.
    - With `--strict`, the output is exactly every MEM at least as long as the t-th longest, ties included.
  - `--read-stats` has a new last column: the fallback's reason, `few_pseudo_mems` (at most T pseudo-MEMs) or `few_long_mems` (fewer than t MEMs reach the bound). `--stats` prints the counts of each.
- **`--kebab-top t --fixed-L`** (KeBaB's own pipeline): report every MEM of length at least max(L, k) inside the t longest pseudo-MEMs, and nothing else. There is no longest-t selection, and `-X` is ignored. Inside a pseudo-MEM a match shorter than k need not be a MEM, hence max(L, k).

- **`--check-seed`** (with `--kebab-check`): if the check fails but at least t MEMs were found in the scanned pseudo-MEMs, the whole-read BML starts with L equal to the t-th longest of them, instead of from L.
  - *Why it is exact:* every MEM found there is a MEM of the read of length at least k (inside pseudo-MEMs only such matches are reported), so the t-th longest found is at most the true t-th longest, and no MEM of the output is shorter.
  - The fallback re-finds the MEMs it had, so it starts with an empty selection and only L carried over.
  - Brute force: the output is identical to the unseeded fallback in every run. Seeded fallbacks: 934 of 6,426 and 922 of 6,458. The backward steps drop by 3.0% and 4.4% in all, but some runs take more steps (377 and 208, against 543 and 671 that take fewer): a larger L means longer windows, and a longer window can cost more steps before it fails.
  - `--stats` counts the seeded fallbacks.

**Tests** (`katka-test mems`, brute force; t = 1, 3, 5; T = t, 2t, t + 2; `-D` 0 and 14; L = k and L = 1; plain and homopolymer-compressed):
- `--strict` gives exactly the t longest with ties.
- Without it, the output contains only MEMs from that set and has the same t longest lengths.
- `--fixed-L` gives exactly the MEMs of length ≥ max(L, k) contained in the t longest pseudo-MEMs.
- Every configuration (RLCSA, RLBWT, two-level, STPD, STPD over the parse) gives the same output.

On the self-test reads: 24,480 runs (6,426 fallbacks) and 2,040 runs; on the homopolymer-compressed copy, 21,600 runs (6,458 fallbacks) and 1,800 runs. All pass.

**For the desktop (Zymo, as in `results/zymo_kebab/`):** for t in 1, 3 and 5, and T in t and 2t:

```
K=~/rz/MEM-finding/tool/bin/katka; H=~/zymo/hpc/zy; R=~/zymo/all.fq; Z="-U -D 14 -c --csa -j 12 --stats"
$K find $Z -L 20 -K 20 --kebab-check $T $t          --read-stats e.$T.$t.tsv  $H $R e.$T.$t     # exact check rule
$K find $Z -L 20 -K 20 --kebab-check $T $t --strict --read-stats es.$T.$t.tsv $H $R es.$T.$t    # strict (ties)
```

Report for each run:
- time and peak memory;
- the `kebab check` line (fallbacks and their reasons);
- `scripts/compare_longest_t.py` against `c.$t` (the streaming rule, exact at L = 20) and against `d.$t`.

The lengths should never differ from `c.$t`'s, and with `--strict` the spans shouldn't either.

**After the current queue:** rerun the t = 3 and t = 5 check runs, with and without `--strict`, adding `--check-seed`, to files `ec.$T.$t` and `ecs.$T.$t`. `cmp` each with its unseeded run: they must be byte-identical. Report the times and the `kebab check` line, whose "seeded" count is the number of fallbacks that started from the seed.

For the paper's pipeline, run `$K find $Z -L 20 -K 20 --kebab-top $t --fixed-L $H $R f.$t` for t = 1, 3, 5. Score it with every MEM voting, as the L = 30 threshold runs in `results/zymo/` were (`credit_schemes.py`, `score_spu2.py`), against the same scoring of `-L 20 -K 20` without `--kebab-top`. The difference is the cost of searching only the t longest pseudo-MEMs.

### KeBaB stopping rules for the journal revision (2026-11-03)

Three rules, all with the t longest MEMs *with ties* (`--longest t`: every MEM at least as long as the t-th longest; L is raised to the t-th longest found so far, not one more, so ties are found). Pseudo-MEMs are sorted longest first. Inside them only matches of length ≥ k are reported: a shorter match that ends at a pseudo-MEM's boundary need not be a MEM. MEMs dropped by `-D` are not kept and do not count anywhere.
- **Streaming** (`-K k --longest t`): search until the next pseudo-MEM is shorter than L. Exact when L ≥ k. `-K --top t` also streams, but `--top` keeps exactly t (the first found among ties), so it is not tie-complete.
- **Heuristic** (`-K k --kebab-top t`): only the t longest pseudo-MEMs, then stop. Pseudo-MEMs are at least k long by construction, so a read whose longest MEM is shorter than k gets nothing.
- **Batch** (`-K k --kebab-batch t`): search the 2t longest pseudo-MEMs, then batches of 4t, 8t, .... After each batch, λ is the length of the next pseudo-MEM (k − 1 if none is left). Stop once at least t MEMs are strictly longer than λ. If the pseudo-MEMs run out first (fewer than t MEMs of length ≥ k), run BML on the whole read, but only if L < k: with L ≥ k, every MEM the fallback could find lies in a pseudo-MEM already searched.

A pseudo-MEM shorter than the current L is skipped (no backward steps) and not counted as searched. So with the zymo preset the batch rule does essentially the streaming rule's work and differs only in when it checks.

**Tests** (`katka-test mems`, brute force on the self-test data; t = 1, 3, 5; `-D` 0 and 14; plain and homopolymer-compressed):
- `--longest` without a filter at L = 1 is exact.
- KeBaB streaming at L = k is exact.
- KeBaB batch at L = k and at L = 1 is exact. The L = 1 case covers 447 whole-read fallbacks and 326 second batches.
- The heuristic reports only true MEMs, but differs from the t longest in 613 of 4,080 runs (synthetic reads, many shorter than k).
- On the self-test reads, `--kebab-batch 3 -L 1` gives byte-identical output to the exact reference `--longest 3 -L 1`.

`scripts/compare_longest_t.py run.txt ref.txt reads.fq t [k]` compares two listings read by read, using the t longest with ties at compressed lengths. It reports the reads whose sorted lengths differ (a longest MEM missed), the reads whose spans differ, and the reads whose reference longest MEM is shorter than k.

**For the desktop (Zymo, `~/zymo/hpc/zy`, `~/zymo/all.fq`; one run at a time, 12 threads, in the 11 GB cgroup as in `results/zymo/run_scripts/runs2.sh`).** After `git pull origin tool && make -C tool` and `tool/selftest.sh`, for t in 1, 3 and 5:

```
K=~/rz/MEM-finding/tool/bin/katka; H=~/zymo/hpc/zy; R=~/zymo/all.fq; Z="-U -D 14 -c --csa -j 12 --stats"
$K find $Z -L 20 -K 20 --kebab-top $t   --read-stats a.$t.tsv $H $R a.$t      # (a) heuristic
$K find $Z -L 20 -K 20 --kebab-batch $t --read-stats b.$t.tsv $H $R b.$t      # (b) batch rule
$K find $Z -L 1  -K 20 --kebab-batch $t --read-stats b1.$t.tsv $H $R b1.$t    # (b') batch rule with L = 1 (the fallback can run)
$K find $Z -L 20 -K 20 --longest $t     --read-stats c.$t.tsv $H $R c.$t      # (c) streaming rule
$K find $Z -L 1        --longest $t                         $H $R d.$t      # (d) exact reference, no filter
```

These answer the desktop's points:
- (d) uses `--longest t` for every t. For t = 1 its longest MEMs are those of `-X -L 1`, and for t > 1 it is tie-complete, unlike `--top`.
- (c) runs for every t.
- `-D 14` is part of every rule's guarantee.
- The pseudo-MEMs are all at least k long.

Report for each run:
- time and peak memory;
- accuracy with `maxmem_roc.py`, which votes with the longest MEMs only, for every t (a t-longest scorer would be an analysis script under `results/`);
- the `kebab (...)` line of `--stats`: pseudo-MEMs and searched per read, batches, reads with a second batch, fallbacks, backward steps per read.

For item 3, run:
- `python3 tool/scripts/compare_longest_t.py a.$t c.$t $R $t`: the heuristic against the exact rule at the same L = 20. This is the number that says how often the heuristic misses.
- The same comparison for `b.$t` against `c.$t`, which should give 0.
- `a.$t`, `b1.$t` and `c.$t` against `d.$t`. Against (d), runs at L = 20 also differ on reads whose t longest include MEMs shorter than 20; the script counts those reads separately.

Put the logs, the `--stats` lines and the script outputs in `tool/results/zymo_kebab/`, then commit and push. Compress the listings with zstd and don't commit them.

### The parse oracle's dictionary at 2 bits per base (2026-11-03)

`PREFIX.ptxt` (format `KTPTXT02`) stores the phrase dictionary packed (`packed_dict`, `src/parse_text.hpp`):
- **the bases:** the phrases concatenated in id order at 2 bits per base, with any character other than A/C/G/T stored as A;
- **the exceptions:** a bit vector over the phrases marking those with other characters. For each marked phrase, its runs of X as (offset, length) and the offset of $. The text has no other characters: `katka prep` writes X for every non-A/C/G/T character, and $ is the terminator.

Decoding is exact (X and $ come back as themselves), so every comparison gives the same order as with the bytes. That matters: the binary searches decide less/greater at X and $ too, not only equal/unequal. Each cursor decodes a window of its current phrase into its own buffer: 64 characters to the left of the position read and 192 to its right, 4 bases per byte by table lookup, with the exception runs found by binary search. Decoding whole phrases was far too slow on phrases spanning N-separated scaffolds: one phrase of the self-test's reads collection holds thousands of X runs. The parse sample's comparisons decode only the prefix they compare. Old `.ptxt` files (`KTPTXT01`) are rebuilt by `index --stpd-parse`; `index --parse-oracle` overwrites them.

**On the synthetic collection:**
- **Size:** the dictionary dropped from 33.1 MB to 8.4 MB (33.07 million bases; 7 phrases with exceptions, 7 X runs), and the oracle from 45.0 to 20.3 MB (2.53 bits per character).
- **Outputs:** the parse sample built from it is byte-identical. Every output of STPD over the parse (with and without `--locate`) and of `--stpd --parse-oracle`, at L = 1000 and 100, is byte-identical to the previous run and has the same MEMs as plain search. `katka-test` reads every position forward and backward, plus random reads and LCEs, against the text.
- **Time** (minimum of 3 runs; was → now): STPD over the parse 0.31 → 0.34 s (L = 1000) and 0.53 → 0.50 s (L = 100); character-level STPD through the oracle 0.72 → 0.94 s and 1.28 → 1.39 s. The character-level search seeks at random once per comparison, and each seek decodes a window.

**Expected on c33:**
- **Dictionary:** the desktop measured 630.5 MB of phrase bytes, so about 158 MB packed, plus about 1–2 MB for the bit vector over 8.5 million phrases and the few exceptions.
- **Offsets:** about 32 MB, unchanged.
- **The whole oracle:** about 1.04 GB instead of 1.48 GB, with ids 661 MB and starts 192 MB. With RLZ for the parse part (the desktop's estimate), about 360 MB.

**For the desktop (c33):**

```
cd ~/rz/MEM-finding && git pull origin tool && make -C tool
cd ~/tool/chr19 && ~/rz/MEM-finding/tool/bin/katka index --parse-oracle c33 2>&1 | grep "parse oracle"
```

That rebuilds `c33.ptxt` (reading only `c33.S`). Its line gives the new dictionary size, the number of phrases with exceptions and the total. (`stpd_parse_chr19.sh` also rebuilds a `.ptxt` of the old format by itself. The `.pstpd` does not depend on the dictionary's format.)

### Locating with two-level search and blocks (`index --twolevel-locate`, `find --twolevel --locate`, 2026-11-04)

Two-level search with blocks now keeps a toehold, SA at the first row of the current interval, through every kind of step (`twolevel_locate` and `search2` in `src/twolevel.hpp`, the builder in `src/tloc_build.hpp`). Every occurrence then comes from φ⁻¹, and each is mapped to its sequence with the sequence starts.

**`PREFIX.tloc`:**

| part | what |
|---|---|
| character run samples | SA[x] − 1 at the start x of every BWT run, in each character's run order (`csamp`, r entries) |
| φ⁻¹ | SA at the end of every run → SA at the start of the next (as `index --stpd`'s) |
| parse run samples | SA_Z at the start of every run of the parse BWT, in row order (`prow`, `psamp`) |
| block run samples | the same for the block BWT, with block indices (`brow`, `bsamp`) |
| block → phrase | Elias–Fano over the first phrase of each block (`bstarts`) |
| phrase → text | Elias–Fano over the start of each phrase (`pstarts`) |

The sequence starts come from `PREFIX.tbl` at query time. `index --twolevel-locate` builds `PREFIX.tloc` with:
- one LF pass over the BWT (characters and φ⁻¹, as `index --stpd`);
- a suffix array of the parse (from `PREFIX.ptxt`, built if missing);
- the run lists of the parse BWT and block BWT, read from the two-level index, with each run's symbol checked against the parse.

**The toehold**, t = SA[sp] (as the r-index keeps one at the last row, mirrored):
- **Character step** c on [sp, ep): if BWT[sp] = c, the new toehold is t − 1. Otherwise the first c in the range is at a row x that starts a BWT run (BWT[x − 1] ≠ c, and x − 1 is in the range). LF maps x to the new sp, so the new toehold is SA[x] − 1, a sample. Its index is c's run number, found from the RLBWT's count of c's runs before the new sp.
- **Parse step** a on parse rows [plo, phi] with toehold u (a phrase): the same over BWT_Z, with the parse samples. Whether BWT_Z[plo] = a, and the first a-run starting after plo, both come from the parse BWT's runs.
- **Block step:** the same over the block BWT, with block indices.
- **Into the blocks.** Every row of [plo, phi] is a suffix starting with the pattern phrase p = pid[i]. A row is a block row exactly when the phrase before it, BWT_Z, is not a phrase z with succ(z) = p. So the first block row x ≥ plo (`BS`) is plo itself, or a row whose BWT symbol is not in {z : succ(z) = p} while the row before it is. That makes x the start of a run, and SA_Z[x] is u or a parse sample. Its block is found with block → phrase.
- **Out of the blocks:** both exits of the block stage (exact, and the backup to the interval before the last block) set plo to a block row. The toehold becomes the first phrase of the block held there (or before the last block step, for the backup).
- **Out of the parse:** the toehold becomes the text start of the phrase, via phrase → text.
- **Locating:** from SA[sp], apply φ⁻¹ ep − sp − 1 times. Locating runs one more two-level search, with the toehold, per reported MEM; the window searches of BML don't pay for it.

**Tests** (`katka-test mems --twolevel`):
- After random two-level searches on both strands, with and without blocks, the toehold equals SA[sp]: 3,120 checks on the chromosome-like collection.
- Every occurrence of every MEM, for BML and FB at L = 1 and 25, equals the suffix array's (110,456 occurrences).
- `selftest.sh` checks `find --twolevel --locate` against plain search: the same MEMs and intervals, and occurrence counts equal to ep − sp.

**On the synthetic collection** (search time, minimum of 3 runs):
- **Time:** two-level with blocks 0.32 → 0.49 s with `--locate` at L = 1000 (1,996 MEMs, 4,410 occurrences), and 0.44 → 0.77 s at L = 100 (4,663 MEMs, 29,460 occurrences). The difference includes the second search per MEM.
- **Agreement:** the MEMs and intervals equal plain search's. Every MEM's number of occurrences and of distinct sequences equals that of character-level STPD with `--locate`, which is checked separately against the suffix array. 3,487 of those MEMs occur in more than one sequence.

**Extra space** on the synthetic collection (n = 64.0 × 10⁶, r = 3.82 million, r_Z = 514 thousand, r_B = 154 thousand), 33.5 MB in all:
- characters 28.1 MB: run samples 12.4, φ⁻¹ 15.7;
- parse 4.6 MB: run samples 1.8, phrase → text 2.8;
- blocks 0.8 MB: run samples 0.45, block → phrase 0.32.

On c33 (r = 68.4 million, r_Z = 10.3 million, m = 220 million), I estimate about 880 MB:
- characters about 625 MB: run samples about 274 MB, φ⁻¹ 350 MB (as in `c33.stpd`);
- parse about 240 MB: samples about 45 MB, phrase → text about 192 MB;
- blocks about 20 MB.

For comparison: two-level search with blocks reads 979 MB (RLBWT, RLCSA, `.lvl`), and STPD over the parse with `--locate` reads 2,378 MB.

**For the desktop (c33, then c65 if it fits):**

```
cd ~/rz/MEM-finding && git pull origin tool && make -C tool && tool/selftest.sh
cd ~/tool/chr19
systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "tloc: %e s, peak %M KB" ~/rz/MEM-finding/tool/bin/katka index --twolevel-locate c33 > tloc33.err 2>&1
for L in 1000 100; do
  taskset -c 9 ~/rz/MEM-finding/tool/bin/katka find -L $L --format mems --csa --twolevel --locate --stats c33 queries.fa tlloc$L.33.tsv 2> tlloc$L.33.err
done
```

`index --twolevel-locate` needs `c33.bwt`, `c33.lvl` and `c33.ptxt`; the last is rebuilt if it is missing or of the old format. Check each output against plain search (`plain$L.33.tsv` from `run_n.sh`): columns 1–5 must be equal, and column 6 (occurrences) must equal ep − sp. Report:
- the build's time and peak, and the `two-level locate` lines (component sizes);
- the search times next to two-level with blocks without `--locate` and STPD over the parse with `--locate`, from `results/chr19_stpd_parse/`.

Put the results in `tool/results/chr19_stpd_parse/`.

### The RLZAP parse oracle (`index --rlzap`, `find --rlzap-oracle`, 2026-11-03)

`PREFIX.rptxt` is the parse oracle with the parse (ids and phrase starts) replaced by an RLZAP-style factorization (`rlz_parse` and `parse_text::make_rlz` in `src/parse_text.hpp`). The dictionary (2 bits per base) and its offsets are unchanged.
- **The reference:** the phrases starting in S[0 .. 2(|G₁| + 1)): the first sequence of `PREFIX.tbl` and its reverse complement (CHM13 on chr19), as packed ids with their starts in Elias–Fano.
- **The factorization:** exactly that of `results/chr19_rlz/rlz_parse.cpp`, J, D, M = 4, 4, 2 by default (`--rlzap=J,D,M`). On the self-test's chromosome-like collection it gives the same counts as that program: 5,505 factors, 11,291 literals (2,994 factors only a literal), 1,162 pointer changes, 61 explicit pointers.
- **Factor f:** n_f literal ids, then a copy of reference phrases at a relative pointer (source = target − rel_f).
  - Stored: the factors' first phrases (Elias–Fano over parse positions) and their text starts (Elias–Fano over S); the literal counts (Elias–Fano over their running sum); the literal ids, packed; the pointers, run-length coded (Elias–Fano over the factors where the pointer changes, and the values there).
  - A copied phrase's start follows from the reference's starts: the phrase spacing depends only on the ids.
- **Access:**
  - Phrase id and text start of phrase i: the factor holding i (one predecessor query, cached per thread; consecutive phrases usually share a factor), then a literal or the reference.
  - The phrase at text position x: the factor by its text start, then its few literals by their lengths, then a predecessor query on the reference's starts, shifted by the pointer.
  - All of STPD over the parse (the middle stage's Z[u − 1], conversions between parse and text positions, the parse sample's comparisons, φ⁻¹ locating) and the character-level oracle go through these three accessors.
- **Check at build time:** every phrase's id and start, and many position lookups, against the plain parse.
- **Tests:**
  - `katka-test` builds an RLZAP oracle and reads every position forward and backward, plus random reads and both kinds of LCE, against the text. It runs STPD over the parse through it in every brute-force MEM check.
  - `selftest.sh` compares `find --stpd-parse --rlzap-oracle` with plain search (L = 1, 100, 1000), and its `--locate` output with the plain parse's.

**On the synthetic collection** (reference: the first of 16 genomes, both strands, 200,801 phrases):
- **The factorization:** the target's 3,012,064 phrases form 198,682 factors (15.2 phrases each), with 464,340 literals.
- **The parse part:** 2.59 MB instead of 10.48 MB (ids 7.63 + starts 2.85). The whole oracle is 12.4 MB instead of 20.3 MB.

Search time (minimum of 3 runs; every output had the same MEMs as plain search, and `--locate` the same output as with the plain parse):

| | L = 1000 | L = 100 |
|---|---|---|
| STPD over the parse, plain parse | 0.32 s | 0.47 s |
| STPD over the parse, RLZAP parse | 0.44 s | 0.62 s |
| STPD over the parse, RLZAP parse, `--locate` | 0.47 s | 0.66 s |
| character-level STPD through the RLZAP oracle | 1.46 s | 2.85 s |

**For the desktop:** `experiments/stpd_parse_chr19.sh` now also builds `PREFIX.rptxt` (`index --rlzap`) and runs STPD over the parse through it, with and without `--locate`, at each L. It checks the MEMs against plain search and the `--locate` output against the plain parse's. Its summary gives the RLZAP oracle's components.

### STPD over the parse (implemented 2026-11-03)

`katka index --stpd-parse PREFIX` and `katka find --stpd-parse` (`src/pstpd.hpp`, and `ps_back`, `ps_extend`, `ps_locate` in `src/mems.hpp`). It is built on an existing two-level index; only its B and its dictionary's hash table are loaded.

**Components:**

| component | file | what it is |
|---|---|---|
| plain parse | `.ptxt` | Z[0..m): phrase ids packed at ⌈log₂(σ+1)⌉ bits |
| boundaries | `.ptxt` | Elias–Fano over the phrase starts b_0 < ⋯ < b_{m−1} (parse position → text position, and back by predecessor) |
| dictionary strings | `.ptxt` | the phrases concatenated in id order, with packed offsets |
| parse sample | `.pstpd` | the right-to-left STPD sample of Z$ (holding lexicographically smallest occurrences): parse positions SA_Z[x] − 1 for the needed run starts x of BWT_Z, sorted by parse suffix |
| entry table | `.pstpd` | entry[y] = SA_Z[C[y]], for the σ ids, and C (the first parse row of each id) |
| parse φ⁻¹ | `.pstpd` | for `--locate`: SA_Z at the ends of BWT_Z's runs → SA_Z at the next run's start |
| character-level sample | `.stpd` | `index --stpd`; for 32 individuals or more, `index --stpd --parse-oracle` builds it reading the text through `.ptxt` (byte-identical, about 3.5 times slower, without holding the text) |
| right end | `.csa` or `.rlbwt`, `.lvl` | backward steps (RLCSA with `--csa`), then B (rank) |
| pattern phrases | `.lvl` | the dictionary's hash table (phrase → id, and its reverse complement's id) |

`index --stpd-parse` builds the parse sample with a suffix array of the parse (libsais_int, 8 bytes per phrase) and one pass over its rows. `find --stpd-parse` prints every component's size.

**Search** (the longest suffix of X[lo..e) that occurs, X = P or rc(P)). The state is a match X[j..e) and its smallest occurrence, held either as a text position t (characters) or as a parse position u (the parse):
1. **Right end.** Let q = S[r] be the last syncmer start of X with q + k ≤ e. If q ≥ lo, take backward steps from e down to q. Then ρ = rank_B(sp) + 1 is the first parse row; y = the id with C[y] ≤ ρ < C[y+1], which in fact satisfies C[y] = ρ (checked at run time); and u = entry[y]. If there is no such q, start with characters from t = SA[0].
2. **Middle.** While the pattern phrase to the left, a = pid[i−1] = id(X[S[i−1]..S[i]+k)), exists, is in the dictionary, and S[i−1] ≥ lo:
   - if Z[u−1] = a, then u − 1 (one read of the packed parse);
   - otherwise binary-search the parse sample for the smallest sampled suffix in the range of (a, pid[i], …, pid[r−1]; X[q..e)). Comparisons read Z by ids (with the lcp skip), then the next phrase's string against X[q..e).
3. **Left end.** Set t = b_u (one select) and run character-level STPD with the character-level sample. Its text reads go through `.ptxt`.
4. **Back into the parse.** In the character stage, whenever j reaches a syncmer start S[i′] whose left phrase can be matched, set u = pred(t). So BML's extension to the left (after a window search) runs over phrases again.

`--locate`: start φ⁻¹ over the parse from the smallest occurrence of the MEM's core X[S[i₀]..e), where i₀ is its first syncmer. Use doubling batches with checks on phrase ids and the right end's characters, then a binary search in the last batch. Then check the left end's characters X[j..S[i₀]) through the oracle, and map each occurrence to its sequence with the sequence starts from `PREFIX.tbl`. A MEM with no core falls back to character-level locating.

Counters (`--stats`): right-end character steps, middle phrase comparisons, parse-sample binary searches and their comparisons, left-end character steps (and character-level STPD jumps).

**Correctness.** Let b_0 = 0 < b_1 < ⋯ < b_{m−1} be the boundaries of S$: 0 and the starts of closed syncmers. Phrase u is S[b_u .. b_{u+1} + k), and the last phrase runs to the end. Three facts are used, the first two already in the two-level argument:
- (F1) Whether S[p..p+k) is a closed syncmer depends only on those k characters (canonical under reverse complement, too).
- (F2) The dictionary is prefix-free.
- (F3) Lemma order: for boundaries b_u, b_v, S[b_u..] < S[b_v..] if and only if Z[u..] < Z[v..], ids being lexicographic.

*Lemma A (alignment).* Let q be a syncmer start of X with q + k ≤ e, and Y = X[q..e). In every occurrence S[p..p+|Y|) = Y:
- p is a boundary;
- for each syncmer start q′ of X with q ≤ q′ and q′ + k ≤ e, p + q′ − q is a boundary;
- no other position in [p, p + |Y| − k] is a boundary.

Proof: by F1, a k-mer inside the occurrence is a closed syncmer in S if and only if it is one in X.

*Entry.* Let q = S[r] be the last syncmer start of X with q + k ≤ e, and Y = X[q..e).
- (i) For every boundary b_u with S[b_u..] prefixed by Y, Y is a prefix of phrase Z[u]. By Lemma A, the next boundary after b_u is beyond b_u + |Y| − k, since otherwise X would have a syncmer start in (q, e − k]. So phrase u, which ends k characters after the next boundary, is longer than |Y|.
- (ii) Conversely, if phrase Z[u] has prefix Y, so does S[b_u..]. So the parse rows whose suffixes are prefixed by Y (as characters) are those whose first id's phrase has prefix Y.
- (iii) Those ids form a contiguous range [y₁, y₂], because the ids are lexicographic. The parse rows are sorted by first id, so the rows form [C[y₁], C[y₂ + 1]).
- (iv) B maps the character interval of Y to these parse rows, because every occurrence of Y is at a boundary (Lemma A) and the order agrees (F3). So rank_B(sp) + 1 = C[y₁], and the parse suffix at that row, entry[y₁], is the smallest suffix of S$ prefixed by Y, at text position b_{entry[y₁]}.

*Middle step.* Let R_i be the parse rows whose suffixes (as characters) are prefixed by X_i = X[S[i]..e): a contiguous range, by F3. Let u = SA_Z[min R_i], and a = pid[i−1].
- By Lemma A, every occurrence of X_{i−1} = X[S[i−1]..e) is a boundary v whose phrase is exactly X[S[i−1]..S[i]+k), so Z[v] = a, and the next boundary is v + 1, with Z[v+1..] in R_i. Conversely, any such v is an occurrence. So R_{i−1} = LF_a(R_i), an ordinary backward step over BWT_Z on a contiguous range. If a is not in the dictionary, X_{i−1} has no occurrence.
- The character-level STPD lemma (`stpd.hpp`) then applies to BWT_Z word for word: if BWT_Z[min R_i] = Z[u−1] = a, then u − 1 holds the smallest occurrence. Otherwise the first row x of R_i with BWT_Z[x] = a starts a run, and SA_Z[x] − 1 holds it.
- *Pruning.* R_i is not a range of an id string. It is the set of rows whose expansions share a prefix of |X_i| characters, so it is an interval of the character LCP between boundary suffixes. Hence if lcp(SA_Z[y], SA_Z[x]) ≥ lcp(SA_Z[x−1], SA_Z[x]) in characters, then any R_i containing x − 1 and x also contains y, and x is not the first a in R_i. So the pruning rule is applied with character LCPs: the shared phrases contribute len − k each, then the LCP of the first two different phrases, which ends inside both by F2. (Id LCPs would be wrong here: two rows can share a ids and differ in the next id while both lie in R_i.)
- *The binary search's comparisons.* Comparing Z[v..] with (a, pid[i..r−1]) id by id decides the order at the first differing id, by F2 and F3. Equality on all of them, followed by Y being a prefix of phrase Z[v + r − i + 1], is membership in R_{i−1}, by Entry (i)–(ii).

*Exit.* When the middle stage stops at i with state u, every occurrence of X_i is a boundary (Lemma A). The smallest in parse order is the smallest in character order (F3). So t = b_u is the smallest suffix of S$ prefixed by X_i. That is exactly the invariant of character-level STPD for the match X[S[i]..e), so it continues correctly to the left. If the stop is a failure, X_{i−1} does not occur, so the longest suffix starts in (S[i−1], S[i]], inside the phrase, where character-level STPD finds it. Re-entry (step 4) is the same argument read backwards: at a syncmer start S[i′], the smallest occurrence t of X_{i′} is a boundary b_u (Lemma A), and u = pred(t) is the smallest in parse order (F3).

*Locating.* Every occurrence of the MEM X[j..e) contains an occurrence of its core X_{i₀} at offset S[i₀] − j (Lemma A). The core's occurrences are the contiguous parse rows from its smallest one, so φ⁻¹ over the parse lists them in row order; membership is the comparison above. The left-end check keeps exactly those preceded by X[j..S[i₀]).

**Tests:**
- `katka-test mems --twolevel --stpd` adds two configurations, with the RLCSA and with the RLBWT at the right end, to every brute-force MEM check. These cover BML and FB on both strands, at L = 1, 12, 25, 40, 100 and 1000, with intervals, and the KeBaB, top-t and longest-t modes.
- It also checks:
  - 1,040 (chromosome-like) and 20,920 (self-test reads) longest-suffix searches against plain search, including BML's extend-left resumption;
  - every occurrence from φ⁻¹ over the parse against the suffix array: 1.6 million and 259 million occurrences.
- `selftest.sh` compares `find --stpd-parse` with plain search byte for byte (BML and FB, L = 1, 25, 100, 1000), and the `--locate` counts with ep − sp.

**On the synthetic collection** (16 × 2 Mbp, both strands, n = 64.0 × 10⁶; m = 3.21 million phrases, σ = 456 thousand, r_Z = 514 thousand against r = 3.82 million; search time, the minimum of 3 runs; every output had the same MEMs as plain search):

| | L = 1000 | L = 100 |
|---|---|---|
| plain search (RLCSA) | 1.38 s | 1.61 s |
| two-level with blocks | 0.33 s | 0.45 s |
| STPD over the parse (RLCSA at the right end) | 0.31 s | 0.53 s |
| STPD over the parse, `--locate` | 0.42 s | 0.65 s |
| character-level STPD, parse oracle | 0.72 s | 1.28 s |
| character-level STPD, plain text (holds the text) | 0.18 s | 0.28 s |

At L = 1000, STPD over the parse made:
- 435 thousand right-end character steps;
- 343 thousand middle phrase comparisons;
- 8,293 parse-sample binary searches (165 thousand comparisons);
- 113 thousand left-end character steps, with 6,362 character-level jumps.

The character-level STPD through the oracle made 7.28 million character steps and 84,619 jumps. The right end dominates: about 55 backward steps per search, since each search starts with up to a phrase of characters.

Components (MB): plain parse 7.6; boundaries 2.8; dictionary strings 33.1 (+ offsets 1.4); parse sample 1.3 (485,007 entries, 0.94 r_Z); entry table 1.3 (+ C 0.3); parse φ⁻¹ 1.8; character-level sample 8.2 (+ φ⁻¹ 15.7); RLCSA 31.6 (or RLBWT 8.7); B 0.6; dictionary hash 5.2. On this collection the dictionary dominates, because the synthetic genomes differ more than human haplotypes do. The pruning keeps 0.94 of r_Z, against 0.66 of r at the character level.

**For the desktop.** Run on the 16-individual index:

```
cd ~/rz/MEM-finding && git pull origin tool && make -C tool && tool/selftest.sh
PIN=9 ~/rz/MEM-finding/tool/experiments/stpd_parse_chr19.sh ~/tool/chr19/c33 ~/tool/chr19/queries.fa ~/tool/chr19/stpd_parse33
```

The script adds `c33.ptxt` and `c33.pstpd` next to c33 (`index --stpd-parse`, which reads only `c33.S` and `c33.lvl`), and `c33.stpd` if it is missing. It runs each step in an 11 GB cgroup. Then, single-threaded and pinned, for L = 1000 and 100, it runs:
- plain search (the reference);
- STPD over the parse, with `--no-intervals` and with `--locate`;
- character-level `--stpd` with the parse oracle;
- two-level with blocks.

Each output is compared with plain search (MEMs, and the `--locate` counts against ep − sp). Times, peaks, counters and component sizes go to `stpd_parse33/summary.txt`.

Then try 32 individuals: run the same command on `~/tool/chr19/c65`, with the work directory `stpd_parse65`. It needs `c65.bwt` once, to build `c65.stpd` through the parse oracle. That build holds O(r) and the oracle, not the 8 GB text; if it exceeds 11 GB, the cgroup stops it, and the summary says so.

Copy the summaries and logs (not the TSVs) to `tool/results/chr19_stpd_parse/`, then commit and push.

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

Memory (rewritten 2026-11-03, after the desktop's 32-individual build peaked at 9.04 GB with m = 434 million phrases, about 21 bytes per phrase): the parse vector grew by doubling (2.1 GB of capacity), libsais_int got a copy of it, the suffix array and the parse-BWT vector were held together with it, and the RLCSA (1.2 GB) and RLBWT stayed loaded throughout. `build_twolevel` (`src/twolevel_build.hpp`) now:
- writes the parse to `PREFIX.lvl.tmp` while scanning S, and reads it back at its exact size;
- calls back to `index` once B is built, which frees the RLCSA and RLBWT;
- derives successors, block starts (one bit per phrase) and block symbols (4 bytes per block; distinct keys through a hash map) from the parse before the suffix array;
- runs libsais_int with the parse itself as input (it is destroyed), so the parse and its suffix array take 8 bytes per phrase in all;
- takes the block BWT from the suffix array;
- reads the parse back and overwrites the suffix array with the parse BWT.

The peak is about 8.6 bytes per phrase plus the dictionary and the RLCSA/RLBWT (until B). Output is byte-identical (smallest and frequent continuation, `--no-blocks`, selftest). On the synthetic collection with `--syncmer 8,2` (25.4 million phrases, 2.9 phrases per block) the step peaks at 0.47 GB instead of 0.85 GB. The step now prints its own peak (`peak memory of this step`, VmHWM reset through `/proc/self/clear_refs`). Expected on chr19: about 4 GB at 32 individuals, and about 8 GB at 64 individuals (m ≈ 870 million), which should fit in 11 GB. `--continuation frequent` counts distinct successor pairs in a hash map (about 40 bytes per distinct pair) instead of sorting all m pairs (8 bytes per phrase).

Not changed: `index --stpd` holds the text (1 byte per character) plus O(r). At 16 individuals it peaked at 8.9 GB (n = 4.07 × 10⁹, r = 68 million), so about 4.8 GB is something else. I haven't measured what. That is the next thing to look at if STPD at 32 individuals is wanted.

**For the desktop (2026-11-03): the new two-level build, at 32 and then 64 individuals.** After `git pull origin tool && make -C tool`, rebuild the 32-individual index in the same 11 GB cgroup as `run_n.sh`, keeping the old `.lvl` to compare:

```
cd ~/tool/chr19 && mv c65.lvl c65.lvl.old
systemd-run --user --scope -q -p MemoryMax=11G -p MemorySwapMax=0 /usr/bin/time -f "index: %e s wall, peak %M KB" ~/rz/MEM-finding/tool/bin/katka index --preset chr19 c65 > index.65.new.err 2>&1
tail -5 index.65.new.err; cmp c65.lvl c65.lvl.old && echo "same .lvl"
```

The log line `two-level: parse ... (peak memory of this step X GB)` gives the two-level step alone; the expected value is about 4 GB. If the `.lvl` files are the same and the peak is well under 11 GB, run `run_scripts/run_n.sh 129` (64 individuals; about 40 GB of disk, so delete `c65.bwt` first if space is short). Put the logs in `tool/results/chr19/` and push.

## Open questions for Travis

- **test-small:** the 40-genome collection, to repeat the KATKA regression on real genomes. Ideally added to the repository, or attached, as the brief suggests.
- **Preset `novel`:** I chose L = 1000 and k = 20; the brief only says "fixed large L". Change it in `PRESETS` in `src/main.cpp` if you want something else.
- **Presets `zymo` and `chr19`:** the `zymo` preset includes `-D 14` for the longest MEMs, as the brief's table says; the paper's KATKA `-X` command line in `katka-readable/README.md` doesn't. The `chr19` preset has L = 1000; give `-L 100` for the other column.
- **Commits:** the commits on branch `tool` are authored as "Claude"; change the author if you prefer.
