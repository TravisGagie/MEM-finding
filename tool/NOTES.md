# Notes on the tool

**Overnight batch (2026-11-07), for the morning.** Each item is in `results/silva/README.md`; this is the short version.
1. **Cliffy's tie rule.** Done. Cliffy counts a read correct if the true genus is among those tied at the top; `cliffy_listing_score.py` reports that as `theirs` and also `strict`. No Cliffy listings are here (they are Gonzalo's), so they couldn't be rescored: he needs to run `score_cliffy_listing.sh`. **Found on the way:** Cliffy's `classify_method2` adds the whole match length to every listed genus, not |M|/|D′| as `silva.tex` lines 38 and 42 say. With that credit on exact lists the tie rule adds 3.47 points (92.50 strict, 95.97 under its rule, 3.6% of reads tied). Naive Bayes, EM and the EM prior gain nothing from it.

**For the chat session (2026-11-06).**
- **New results:** `results/silva/README.md`, last sections, with `narrow_tables.txt`. Lookup table at every step: neither (a) nor (b) speeds up search. RLCSA with sampled tags and no RLBWT: full 1,401.6 MB and 45.69 µs; s = 4 995.5 MB and 59.33 µs; built-in grammar 1,015.0 MB and 53.03 µs. The default's footprint is 1.40 GB of files, peak 1.43 GB.
- **Baseline:** 45.69 µs, on the fixed build. A refactor had slowed search by 3.5% (47.13 µs); it is undone, and those two sections say which numbers came before the fix. The lookup-table table is the one measured before it.
- **The index-size note:** the paper's 1.44 GB and the new 1.40 GB count nearly the same things: RLCSA and tags, with the table adding 18 MB to the paper's. So the "1.56 GB" entry in earlier notes is obsolete once 2c is in.
- **Zymo:** grammar tags make the array larger (687.4 to 750.5 MB, ratio 1.09). `results/zymo/README.md`.
- **Cliffy's figures:** read from the PDFs' vector graphics, not from result files. Draft corrections: 87.5–90.2 becomes 87.4–90.2; 93.3–94.8 becomes 93.2–94.7. The draft's 199–204 µs assumes each plotted time covers both mates, which the public files don't settle.
- **Movi 2** (`results/silva/README.md`, the Movi 2 section; this replaces the earlier table, whose times included 7 µs of the prototype's setup): on a shared BWT, all 8,532,597 MEMs, intervals and tag answers are identical to ours. BML uses its backward steps directly: its bidirectional `--mem` is about 10× slower on SILVA's wide intervals. Full tags: 44.73 µs and 1,475 MB, against the RLCSA's 46.70 and 1,402. Built-in grammar: 53.35 and 1,088, against 53.98 and 1,015. Sampled tags favour the RLCSA: s = 4 is 88.98 against 60.55. The select structure (B and E) gives identical answers but is larger (308 MB against 266) and slower (49.99 at full tags, 149.40 at s = 4). It found a `#`/`%` bug and a missing alphabet check in Movi 2, to raise with Mohsen.
- **Item 3:** `find --classify bayesem` matches the scorer call for call on real SILVA (`a4ad83f`).

**For the chat session (2026-11-05).** The desktop now owns `tool/` (`HANDOVER.md`).
- `results/SILVA_DIGEST.md` now has the 5-repeat timing and what each index size counts. Its table row for the default (52.3 µs, 145,500 reads/s) is superseded by the bullets under it.
- `find --classify bayes` now matches the scorer call for call on real SILVA: both mates (`--mate2`) and scored genera only (`--exclude-docs`).
- MicrobeMixer and ART: `results/silva/README.md`, last section. The mismatch is ART's profile, MiSeq v3 in Cliffy's Figure 2 caption against `-ss MSv1` in the code. MSv3 doesn't explain the gap to Cliffy's published LCA accuracy: range with Cliffy's scoring is 87.33 (MSv1) and 86.86 (MSv3), mean of 3 repeats each, against 87.5 to 90.2 published.
- **Better estimates than the single-repeat numbers:** the existing reads are the best of three MSv1 repeats, so range 86.68, proportional 93.59, normalised 92.66 and naive Bayes 94.57 (strict; 3 repeats) replace 86.95, 93.77, 92.76 and 94.70 where a mean is wanted. The paper's reproduced values remain as measured on those reads.
- **New:** EM, 20 iterations, gets 96.32 on the aquatic reads, the best of any scheme there.
- **Naive Bayes with an EM prior** (`--mode bayesem`; `results/silva/README.md`, last section): 95.94 on the aquatic reads against plain EM's 96.32, 92.57 / 93.94 / 94.91 on the dense dial against 92.31 / 90.20 / 90.22, 96.22 / 96.51 on the sparse dial against 96.80 / 95.40. It has the highest mean over these sets; plain EM is best on the aquatic reads.
- **Cliffy's Figure 2**: `results/silva/README.md`, "Cliffy's published times". Time per read depends on whether a plotted time covers both mates, which isn't public; the draft's 199–204 µs assumes it does.
- **Sparse dial** (`results/silva/README.md`, "The sparse dial"): the aquatic abundances on random genera. EM gets 96.32 / 96.80 / 95.40 (aquatic / uniform draws / smallest genera), proportional 93.59 / 88.98 / 83.78, normalised 92.66 / 92.10 / 91.22, naive Bayes 94.57 / 94.69 / 94.54, range 86.68 / 84.72 / 84.40. So EM's lead is sparsity, and proportional's advantage is that the aquatic genera are large.
- **Cliffy's Figure 3** read exactly from the PDF, per region and digestion: `results/silva/README.md`, "Cliffy's published accuracies", with `cliffy_fig3_values.txt`. Undigested LCA is 87.42, 90.16, 88.31, 89.76 and listing 94.06, 94.72, 93.18, 94.17 (V1–V2, V3–V4, V4, V4–V5), so the draft's 87.5 should be 87.4 and its 93.3–94.8 should be 93.2–94.7.

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

### For the desktop's SILVA queue (2026-11-04): what is prepared, so that the desktop uses it instead of writing its own

The queue is: one timing protocol, credit accuracies on all four regions, B_λ, grammar tags (built-in RePair against bigrepair), the digested index, and a uniform-abundance MicrobeMixer simulation. Three of these have something ready here, tested on synthetic or toy data:

- **The timing protocol (`find`, every run).** After its `reads=` line, `find` now prints:
  ```
  timing: N reads, T threads; mean per-read time A us; total time / reads B us (total W s)
  ```
  - **The mean per-read time** times each read alone in its thread: its search and its output formatting, then averaged over the reads.
  - **The total time / reads** is the wall time of the whole search (reading, searching, writing; loading the index excluded) divided by the number of reads.
  - With one thread the two differ only by reading and writing. With T threads the second is about the first divided by T.
  - So no post-processing is needed. Report both, as the protocol asks, and take the thread count from the same line.
- **The uniform-abundance simulation: `experiments/uniform_microbemixer.sh`.**
  ```
  KATKA=~/tool/silva/bac EXP=~/MEM-finding/katka/cliffy-experiments/src TRAV_TO_LENGTH=$S/work/ref/trav_to_length.txt \
    tool/experiments/uniform_microbemixer.sh ~/rz/MicrobeMixer $FASTA $S/exp1_data/tax_slv_ssu_138.1.txt \
    $S/work/ref/rna/doc_to_traversal.txt ~/tool/silva/uniform 200000
  ```
  (`S=~/MEM-finding/katka/silva`; `$FASTA` is the SILVA reference `run_silva.sh` simulates from.)
  - It lists every genus of the reference (`doc_to_traversal.txt`) in an abundance table, then runs MicrobeMixer per region (V1_V2, V3_V4, V4_V4, V4_V5), with the same primers and command line as `katka/silva/run_silva.sh`.
  - It checks that every genus with reads has the same number of them, to within one.
  - With `KATKA` set, it runs `katka find --preset silva` on both mates and scores them with `scripts/score_silva.py` (proportional credit, `MODE=freq`; other schemes with `MODE`). It writes `summary.txt`.
  - **Every genus gets a share of 0, on purpose.** MicrobeMixer gives each genus its share of the reads *still left*, then spreads the rest evenly. With equal nonzero shares the first genus gets about e times the reads of the last: on a 12-genus toy, 270 down to 147 of 2,400, and the factor approaches e ≈ 2.7 with SILVA's thousands of genera. With shares of 0 every read goes through the even spread: 200 each on the toy. Genera with no region between the primers get no reads, as in any MicrobeMixer simulation; `counts.txt` per region shows how many.
  - **Genera that ART can't give enough reads are dropped (2026-11-04).** On the real SILVA reference, MicrobeMixer stopped at V1_V2 with `not enough reads in file`, with 29 of 3,313 genera short. MicrobeMixer extracts regions longer than 200 bases, but ART makes no 250-base reads from a region shorter than 250. It asks for ceil(needed / sequences) × 2 reads per sequence and asserts strictly more than needed. So a genus whose regions are mostly 201–249 bases long gets too few reads: none, or exactly as many as it needs.
    - The script now drops such genera from that region's table and reruns MicrobeMixer, round by round, at most 10 rounds. It doesn't change MicrobeMixer.
    - A genus is dropped when its `genus_<id>_reads_1.fq` has no more reads than MicrobeMixer asks of it: N // g, plus one for the first N mod g of the g genera with extracted regions.
    - Each region's drops are listed in `<region>/dropped.txt`, with the reads made and needed, and counted in `OUT/dropped.txt`. The runner copies them to `uniform/<region>.dropped.txt`, and the tables report them.
    - The uniformity check runs after the drops.
    - Tested on the toy with one genus's regions all cut to 230 bases and two of another's: round 1 failed with MicrobeMixer's assertion; the two genera were dropped (0 and 134 reads, 200 needed); round 2 passed, with 240 read pairs for each of the 10 genera left.
    - Each round reruns the whole simulation for that region.
    - A partial run needs no clearing: a region without its `done` file is started afresh.
  - MicrobeMixer needs `art_illumina` (ART; on Ubuntu, `apt install art-nextgen-simulation-tools`) and the Python modules `regex`, `aiohttp` and `pandas`.
  - Tested end to end on a toy SILVA: 12 real SILVA 138.1 genus paths, 3 synthetic 16S-like sequences each with V4 primer sites, 2,400 read pairs. All 12 genera got 200 read pairs, and genus accuracy was 100%, as expected, since the toy genera share nothing.
- **Grammar tags:** `index --tags --grammar` (built-in RePair) and `--bigrepair P` exist already (see "SILVA" below). The comparison script and its desktop commands are in the next section.

### SILVA queue: runner and tables (2026-11-04)

Everything the desktop's SILVA queue needs, as one runner, `experiments/silva_queue.sh`, and one script that prints the paper's tables from its outputs, `scripts/silva_tables.py`. Please use these rather than new scripts, so that every number in the tables comes from the same protocol.

**Setup on the desktop (once):**

- bigrepair, for the grammar step:
  ```
  cd ~ && git clone https://gitlab.com/manzai/bigrepair && make -C bigrepair && pip install psutil
  ```
  (bigrepair's driver is a Python script that needs psutil.)
- The rest is in place already: the undigested index `~/tool/silva/bac` (`index --preset silva`, see "SILVA" below; it has the tag-run cache `bac.tagruns`), `~/tool/silva/list.tsv`, Big-BWT in `katka/deps`, MicrobeMixer in `katka/MicrobeMixer`, cliffy-experiments in `katka/cliffy-experiments`, and the reads and truth under `katka/silva/work`.

**The queue, in its order** (each step skips work whose output exists, so a step can be rerun after a failure):

```
cd ~/MEM-finding && git pull && make -C tool
BIGREPAIR=~/bigrepair/bigrepair tool/experiments/silva_queue.sh timing credit blam grammar digests uniform
tool/experiments/silva_queue.sh collect                  # the small outputs to tool/results/silva_queue, to commit
python3 tool/scripts/silva_tables.py tool/results/silva_queue --latex
```

Outputs go to `Q=~/tool/silva/queue`. Every timed run is capped at 11 GB (`MEMMAX`) and logs its wall time and peak memory.

- **timing.**
  - The timing protocol, used by every timed step: `katka find` with one thread pinned to CPU 9 (`PIN`), then with 12 threads (`THREADS`), on both mates of every region.
  - It checks that the 12-thread output is byte-identical to the one-thread output.
  - The tables take the time per read from the `timing:` lines: the mean per-read time with one thread, with the total time over the number of reads beside it. Reads per second is all the reads over the summed wall time of the 12-thread runs.
  - The configuration is `find --preset silva` (L = 30, counts, explicit RLCSA, `-F 10`) on `bac`.
- **credit.**
  - Scores the timing step's outputs, the same MEMs at L = 30, under six credit schemes: `lca` = range, `one`, `cliffs` = sample, `list` = even, `norm` and `freq` = proportional. It uses `scripts/score_silva.py`, one region at a time.
  - `norm` needs the number of sequences per genus, which the step counts from `$S/work/ref/dna/doc_<i>_seq.fa`.
  - These six give the SILVA column of `tab:credit`.
- **blam.**
  - Runs `index --tags --tag-lcp 0,10,30,100,1000` on `bac`. This recomputes the tag runs and B_λ with the per-genus tags, and prints ρ, r, σ, d and the lemma's bound for each λ.
  - The tables recompute the bound from these numbers as a check.
  - It peaked at 6.5 GB when `bac`'s tags were built (`results/silva_profile`), and B_λ's records are streamed (`--blam-mem`, default 6 GB), so it should fit in 11 GB. If it does not, add `--blam-external`.
- **grammar.** Runs `experiments/grammar_tags.sh` on `bac` (below). Then, with each of the two grammar tag arrays:
  - It runs the timing protocol with `--preset silva -T <tag array>`.
  - It checks that every output is byte-identical to the default's.
  - It scores the outputs under proportional credit.
- **digests.**
  - Builds the digested index `~/tool/silva/dg` if it is not there:
    - `prep --digest 4,11 list.tsv dg`;
    - Big-BWT on `dg.S`;
    - `index --csa --efcsa --tags dg`;
    - `index --tags --sample 4 --tag-out dg.s4.tag dg`.
  - For L = 15, 20, 30, 40 and 50 it runs the timing protocol with `-L $L -c -F 2` and each locate structure: `--csa` (explicit RLCSA), `--efcsa` (Elias–Fano), and, for L ≥ 30, no RLCSA with `-T dg.s4.tag` (RLBWT with tags sampled at s = 4).
  - It checks that the Elias–Fano and sampled outputs are byte-identical to the explicit RLCSA's, and records this in `digests/identical.txt`.
  - It scores the explicit RLCSA's outputs under proportional credit.
  - It then runs the undigested default (`--preset silva -L $L`) at L = 15, 20, 30 and 40.
  - Together these fill every cell of `tab:digest`, and the "minimizer digests" row of `tab:others` (L = 30).
- **uniform.**
  - Runs `experiments/uniform_microbemixer.sh` with the defaults of `katka/silva/run_silva.sh` (MicrobeMixer, `SILVA_138.1_SSURef_NR99_tax_silva.fasta`, 200,000 pairs per region; `NREADS` to change it).
  - It then copies the per-region read counts and scores to `uniform/`.
- **collect.** Copies everything except the finds' outputs (`*.out`, large) and MicrobeMixer's reads to `OUT` (default `tool/results/silva_queue`). That is everything `silva_tables.py` reads, a few MB.

**Space** in the tables is the total size in GB (10⁹ bytes) of the files a configuration loads, recorded in `Q/sizes.txt`:
- default: `bac.rlbwt`, `bac.csa`, `bac.tag` and `bac.tag.rmq`;
- grammar: `bac.rlbwt`, `bac.csa` and the grammar tag array;
- digests: `dg.rlbwt`, `dg.map` and tags, with `dg.csa`, `dg.efcsa` or `dg.s4.tag` as the row says.

**Expected cost** (from the published per-read times, 53 to 220 µs):
- Each one-thread run of 200,000 reads takes 10 to 45 s, and there are about 200 runs.
- Each scoring takes a minute or two with 12 threads, and there are about 70.
- So the queue should take a few hours, plus bigrepair (below) and MicrobeMixer's simulation.

**`silva_tables.py Q [--latex]`** prints, in order:
1. `tab:credit`, the SILVA column: strict genus accuracy, the mean over the regions and each region.
2. `tab:others`, our rows:
   - default, grammar (built-in RePair and bigrepair) and minimizer digests;
   - for each: space, time per read (and total over reads), reads/s with 12 threads, and genus accuracy (strict, and with ties counted correct).
3. `tab:digest`: space, the time per read for L = 15 … 50, and the strict accuracy rows.
4. B_λ for each λ, with the lemma's bound, ρ/bound, and whether the bound holds.
5. The grammar sizes: the plain, built-in RePair and bigrepair tag arrays; the build times and peaks; whether the answers are byte-identical.
6. The uniform simulation: reads per genus and genus accuracy per region.

`--latex` adds each table's rows in the paper's layout. **A missing input is never a blank:**
- A missing file prints `MISSING: <file> (step <step>)` in its cell.
- A file that is there but unusable prints `BAD: <file> (<why>; step <step>)`. The reasons are a missing `timing:` line, a run that did not exit 0, a bound that does not add up, or a configuration absent from `sizes.txt`.
- The script ends with the list of everything missing and exits 1. `\missing{}` stands in for these cells in the LaTeX rows.

**Tested on a toy SILVA** (12 genera, 2,400 read pairs; the same reads under all four region names):
- Every step ran; the tables were complete and the script exited 0.
- The Elias–Fano and sampled outputs were byte-identical to the explicit RLCSA's in all 64 runs.
- The grammar outputs were byte-identical to the default's.
- The lemma's bound held for every λ.
- Every accuracy was 100%, as expected, since the toy genera share nothing, so the toy checks the plumbing, not the schemes.
- The failure cases were shown correctly: with a scoring deleted, a run's log marked as failed, a configuration removed from `sizes.txt`, and the B_λ log and the uniform step's outputs deleted, the script named each of the 13 inputs and exited 1.

#### Grammar tags: built-in RePair against bigrepair (`experiments/grammar_tags.sh`)

```
BIGREPAIR=~/bigrepair/bigrepair tool/experiments/grammar_tags.sh ~/tool/silva/bac \
    katka/silva/work/compare/reads_every50/V1_V2_1.fq ~/tool/silva/queue/grammar -L 30 --csa -F 10
```

(The `grammar` step of the runner does exactly this.) It writes `WORK/summary.txt` with each step's time and peak memory:
1. `index --tags --grammar --tag-out WORK/g.tag bac`: the built-in approximate RePair, from the tag-run cache `bac.tagruns`.
2. `scripts/tagruns_to_int32.py bac.tagruns WORK/tags.int32` (the runs' tags as 32-bit integers, streamed), then `bigrepair -i -r WORK/tags.int32`.
3. `index --tags --bigrepair WORK/tags.int32 --tag-out WORK/bg.tag bac`: the tag array from bigrepair's grammar (`.R`, `.C`).
4. `find -c` and `find -l` with `bac.tag`, `g.tag` and `bg.tag`. The three outputs must be byte-identical, since only the tag array's encoding differs; the script stops if they are not.

**Tests.**
- **Toy** (`selftest` data, 1,004,033 tag runs). Answers byte-identical under `-c` and `-l` with all three tag arrays. File sizes: plain 0.98 MB (+ RMQ 0.76), built-in RePair 0.91 MB, from bigrepair 1.28 MB.
- **Synthetic** (64.0 million rows, 63.7 million tag runs, one tag per genome; 3,000 reads). Answers byte-identical under `-c` and `-l` with all three tag arrays.
  - Built-in RePair: 56 s, peak 1.62 GB.
  - bigrepair (`-m 0`): 53 s, peak 0.83 GB. Building the tag array from its grammar took 8.6 s, peak 1.03 GB.
  - Tag array files: plain 53.8 MB (+ RMQ 48.6 MB), built-in RePair 42.6 MB, from bigrepair 48.8 MB.
  - These tags are nearly random: 16 tags, one run per row on average. So this measures cost, not the compression SILVA's tags will get.

**Memory for SILVA's 359,351,512 tag runs** (n = 1,173,771,557, r = 68.3 million), scaled per run from the synthetic test. That test has 63.7 million runs, about 1/5.6 of SILVA's.
- **Plain tag array from the cache:** 0.99 GB, 15.5 bytes per run, so about 5.6 GB. This agrees with the 6.5 GB measured when `bac` was built, which also built the RLCSA.
- **Built-in RePair:** 1.62 GB, 25.4 bytes per run, so **about 9.2 GB** with SILVA's 0.16 GB RLBWT. That is under the 11 GB cap, but not by much. Item 9 of the review below estimated 16 bytes per run; the measurement replaces that.
- **The tag array from bigrepair's grammar:** 1.03 GB, 16 bytes per run, so about 5.8 GB.
- **bigrepair:** its plain RePair (`-m 0`, the script's default) peaked at 0.83 GB on 64 million runs, 13 bytes per run, so about 4.7 GB for SILVA if its dictionary is as large as the tag sequence (the worst case; it is smaller when the sequence repeats).
  - **Use `-m 0`.** bigrepair's memory-bounded RePair (`largeb_repair`) is what it uses for any `-m` above 0, and by default with 95% of the machine's RAM, which is also above the 11 GB cap. It stalled on the 64 million runs: more than 30 minutes, killed. It took 30 s on either half of them, and 0.7, 8 and 14 s on the first 1, 4 and 16 million.
  - The script also gives bigrepair a time limit, `BGTIME` (default 8 h). If bigrepair does not finish in time, the rest of the step runs without its tag array, and the tables mark the bigrepair row as not finished.

### SILVA: credit schemes when the sample doesn't match the database (2026-11-05)

On SILVA reads simulated as in Cliffy's paper, proportional credit is the best scheme. The question is whether that holds only because the sample's genera follow the database's composition, and whether a guard helps when they don't. For the desktop, after its current SILVA queue; everything below was tested on a toy only.

**What the existing read sets are.** `katka/silva/run_silva.sh` simulates them with MicrobeMixer from its `data/aquatic.tsv`: the 100 most abundant genera of an MGnify aquatic sample, with MGnify's read counts as shares. Their composition comes from that sample, not from SILVA's N_g (the number of SILVA sequences of genus g). Whether it happens to track N_g, and which α it corresponds to, is what `scripts/silva_alpha_fit.py` measures (below). The fit could not be done here: it needs the reference's N_g, which needs the SILVA FASTA, and arb-silva.de is blocked from the cloud container.

**The read sets with a dial: `ALPHA` in `experiments/uniform_microbemixer.sh`.**
- Genus g gets exactly k_g read pairs, k_g ∝ N_g^α, for α = 1, 0 and −1, with the same 200,000 pairs per region and the same four regions as the existing experiment.
- In the `dial` step every α, 0 included, goes through `scripts/simulate_dial.py` (`DIAL=1`), so the three columns use the same simulator. The uniform MicrobeMixer run (`ALPHA=0` without `DIAL=1`, the `uniform` step) is unchanged and stays a separate check.
- `scripts/simulate_dial.py`: It uses MicrobeMixer's own code to read SILVA, extract each genus's regions (2% mutated) and run ART, with MicrobeMixer's command line. It sets the counts itself, since MicrobeMixer's shares can't give exact counts:
  - k_g is the largest-remainder rounding of N · N_g^α / Σ N_h^α, in exact rationals, over the genera with at least one region of 250 or more bases (ART makes no reads from a shorter one). The counts sum to N exactly.
  - ART gets ceil(k_g / long regions) + 1 pairs per amplicon. A genus that falls short is rerun with twice the fold, and if it still falls short it is dropped and the counts are reallocated.
  - Exactly k_g pairs are then sampled per genus.
  - The truth file is MicrobeMixer's format, plus `<region>_counts.tsv` (N_g, regions, long regions, k_g).
- N_g is the same table the scorer uses: the number of sequences of each document of the reference (`Q/nseq.txt`, which the runner writes).

**The new schemes in `scripts/score_silva.py`.** All work from the existing `{doc:count}` listings and the N_g table (`--nseq`); N = Σ N_g and G is the number of genera.
- **The existing "normalised" scheme (`norm`) is exactly c_g / N_g.** It computes c_g / (2 N_g): the text holds both strands, and so does the count. Each match's credits are normalised to sum to its length, so the factor 2 cancels.
- **`--mode prior`:** credit ∝ (c_g / N_g) · ((1 − β) N_g / N + β / G), for each value in `--beta` (default 0, 0.1, 0.25, 0.5, 0.75, 1). β = 0 is proportional credit and β = 1 normalised credit, bit for bit.
- **`--mode mix`:** (1 − β) · proportional + β · even, per match, for the same β. β = 0 and β = 1 are proportional and even credit, bit for bit.
- **`--mode em`:** credit ∝ (c_g / N_g) · π_g. π starts as the database prior N_g / N and is re-estimated as each genus's share of the total credit over all reads and both mates. Results are reported after each number of iterations in `--em-iters` (default 5 and 20); 0 iterations is proportional credit, bit for bit.
- Each β or iteration count is its own output line, labelled `prior0.25`, `mix0.5`, `em20`, and so on.
- `--levels` limits the levels scored (default all six).
- Every line now ends with `strict_correct=`, the exact count behind `strict=`. Otherwise the existing schemes are unchanged: their outputs on the toy were identical, line for line, to those of the previous version.
- **Tests: `tests/test_score_schemes.py EXP TAX D2T T2L`** (42 checks, 10 s; the naive-Bayes checks are listed under "Naive Bayes with practical priors" below):
  - the endpoints, at every level, on 400 random reads full of exact ties;
  - every scheme against an independent implementation of its formula (exact rationals; floats for EM) on 400 random reads without ties;
  - hand-worked cases:
    - one match {A:1, B:3}, with N_A = 1, N_B = 9 and 12 genera. Proportional credit picks B; normalised credit picks A; the prior mixture picks A exactly when β > 9/14, so at 0.75 and 1; the linear mixture never picks A alone.
    - a read {C, D}, true genus C, among 10 reads that match only C. Proportional credit ties; EM learns that C is common and picks it.

**The runner step and the table.** `silva_queue.sh dial` does three things for α in `ALPHAS` (default "1 0 -1"):
1. It simulates the reads, every α with `simulate_dial.py`.
2. It runs `katka find --preset silva` on each mate, untimed, with `THREADS` threads.
3. It scores every scheme at genus level: lca, one, cliffs, list, norm, freq, prior (6 β), mix (6 β), em (5 and 20 iterations), and naive Bayes over the MEMs with three priors:
   - `bayes`, the sample's own prior (`--prior sample --alpha α`): the reference, since a real classifier can't know α;
   - `bayesflat`, a flat prior, for ε = 10⁻², 10⁻³, 10⁻⁴, 10⁻⁶ and 10⁻⁹ (labels `bayesflat0.01` … `bayesflat1e-09`, with `bayesflat` for 10⁻⁶);
   - `bayesdb`, the database prior N_g / N.

Outputs go to `Q/dial/a<α>/{credit,out}`, with the truth files beside them; `collect` copies everything but the reads and the listings. `scripts/dial_tables.py Q [--latex]` then prints:
- one summary table: each scheme's strict genus accuracy, the mean over the regions, for each α;
- one table per α, with each region.

Accuracy is defined as in `silva_tables.py` (strict, the mean of the regions' accuracies), from the exact counts. A missing input is printed as `MISSING:` or `BAD:` with its file, never as a blank.

**Toy test: `tests/dial_toy.sh MICROBEMIXER EXP D2T WORK`** (47 s), on `tests/make_dial_toy.py`'s SILVA:
- **The toy:** 12 real genus paths, in three groups of four. Sizes range from N_g = 1 to 50. Genera of a group differ in 0.6% of their bases. All four regions' primer sites are present.
- **The run:** the reference is built and the `dial` step runs with 2,400 pairs per region.
- **Checks:**
  - the counts are exactly those asked for, and the α fit recovers 1.002, 0.000 and −0.999;
  - the endpoints are equal;
  - at α = −1, normalised credit and EM beat proportional credit.
- **Results** (strict genus accuracy, mean of the four regions):

  | scheme | α = 1 | α = 0 | α = −1 |
  |---|---|---|---|
  | proportional | 98.1% | 87.8% | 79.4% |
  | normalised | 99.1% | 99.6% | 92.0% |
  | prior mixture, β = 0.5 | 99.1% | 96.0% | 87.2% |
  | EM, 20 iterations | 98.1% | 98.8% | 90.1% |

  This is a toy; on it, normalised credit wins even at α = 1. The real data decides.

**Why normalised credit beats proportional credit on the toy at α = 1 (2026-11-05).**
- **Not how sequences are chosen within a genus.** MicrobeMixer, both in the existing read sets and in the uniform run, extracts a region from every sequence of the genus that has exactly one match of each primer. It gives every extracted region the same fold coverage in ART, and draws the genus's reads uniformly from the pooled output. `simulate_dial.py` does the same, with equal fold per region of at least 250 bases. In both, a genus's reads are spread uniformly over its sequences that have a long enough region, not taken from one or a few.
- **It is summing credit over a read's MEMs.** The Bayes reference (`--mode bayes`) multiplies the matches' likelihoods instead, and uses the sample's own prior. On the toy's reads at α = 1 it gets 99.65%, against 99.21% for normalised and 98.16% for proportional credit. With the wrong prior (α = 0 or −1 on the same reads) it still gets 99.36% and 99.31%, so the prior hardly matters and the summing does.
  - A read from a small genus A shares some of its MEMs with the large genus B of its group, and B gets most of each shared MEM's credit. Summing rewards B for every MEM it shares and doesn't penalise it for the MEMs it lacks. The product does penalise it, through ε for each missing match.
  - A test case: matches {A:1, B:5} of length 200 and {A:1} of length 100, with N_A = 1 and N_B = 10. Proportional credit gives B 167 against A's 133; the Bayes reference picks A with any of the three priors.
- **A toy caveat.** MicrobeMixer mutates each extracted region once, at extraction, with 2% of its bases replaced (1.5% changed). On this toy, genera of a group differ in only 0.6% of their bases. So a genus with a single sequence can come out looking like a neighbour, and its accuracy depends on the random draw.
  - With seed 1, normalised credit at α = 0 got 93.1% on V1_V2, because 151 of genus 4's 200 reads went wrong; with seeds 2, 3 and 4 it got 99.5%, 99.1% and 98.1%.
  - That variance is the toy's; real genera differ by several percent in 16S.

**Naive Bayes with practical priors, and its cost (2026-11-05).**
- **On the toy** (`tests/dial_toy.sh`; strict genus accuracy, mean of the four regions):

  | scheme | α = 1 | α = 0 | α = −1 |
  |---|---|---|---|
  | proportional | 98.32% | 87.26% | 79.24% |
  | normalised | 99.22% | 95.97% | 91.85% |
  | naive Bayes, the sample's prior (reference) | 99.58% | 98.34% | 97.46% |
  | naive Bayes, flat prior, ε = 10⁻² | 95.79% | 98.61% | 99.75% |
  | naive Bayes, flat prior, ε = 10⁻³ | 99.05% | 99.36% | 99.35% |
  | naive Bayes, flat prior, ε = 10⁻⁴ | 99.23% | 98.45% | 97.48% |
  | naive Bayes, flat prior, ε = 10⁻⁶ | 99.31% | 98.34% | 97.28% |
  | naive Bayes, flat prior, ε = 10⁻⁹ | 99.30% | 98.34% | 97.28% |
  | naive Bayes, database prior | 99.58% | 98.33% | 97.20% |

  - The flat prior equals the reference at α = 0, and the database prior equals it at α = 1, as they must.
  - Below 10⁻⁴, ε hardly matters. A large ε acts as a prior toward small genera: a small genus has p_g = c_g / (2 N_g) near 1 on the MEMs it has, and pays little for those it lacks. That helps at α = −1 and hurts at α = 1.
  - On the toy, ε = 10⁻³ is the most even over α. SILVA will decide.
- **Tests** (`tests/test_score_schemes.py`):
  - the flat and database priors against the sample prior with α = 0 and α = 1, line for line;
  - each ε, and the database prior, against an independent implementation;
  - a hand case where ε decides: truth A (N_A = 1) has p = 1 on two MEMs and lacks a third, and B (N_B = 10) has p = 1/20 on all three, so A wins only for ε > 1.25 · 10⁻⁴;
  - the read where summing picks B and the product picks A.
- **What it needs.** Only what the counting listings (`find -c`) already give, plus N_g. The listings are never truncated: each MEM comes with every genus containing it and its count c_g, and a genus missing from a MEM's list has c_g = 0. No new search output is needed.
- **Cost per read.** Both schemes are linear in the listing's (genus, count) entries. Naive Bayes takes one logarithm per entry where proportional credit takes one division, plus one exponential per candidate genus.
  - Measured in Python on the toy's listings (6.0 MEMs and 11.5 entries per read pair): parsing a read pair's listings takes 4.7 µs, proportional credit 3.1 µs more, naive Bayes 4.2 µs more.
  - That is about 1.35 times proportional credit's own work, or 1.14 times with parsing.
  - Both are small next to the search itself (53 µs per read for the default SILVA configuration).

**The experiment, for the methods section** (exactly as the code does it; `scripts/simulate_dial.py`, `scripts/score_silva.py`, `experiments/silva_queue.sh dial`, `scripts/dial_tables.py`).

*Reference and notation.* The reference is the tool's SILVA index: one document per genus, built from SILVA 138.1 by `katka/silva/build_katka.sh`, with G documents (genera). N_g is the number of sequences in genus g's document, and N = Σ_g N_g over all G genera. For a MEM M of a read, of length |M|, c_g(M) is its number of occurrences in genus g's sequences and their reverse complements, as `katka find --preset silva` reports it: L = 30, counts, explicit RLCSA, a lookup table of 10-mers. D(M) = {g : c_g(M) > 0}, and C(M) = Σ_{g∈D(M)} c_g(M).

*Reads.* For each region (V1_V2, V3_V4, V4, V4_V5, with MicrobeMixer's primer files) and each α ∈ {1, 0, −1}, R = 200,000 read pairs are simulated:
1. Every SILVA 138.1 sequence whose lineage reaches the genus rank is assigned to that genus, by MicrobeMixer's code.
2. For each sequence of each genus of the reference, MicrobeMixer's in-silico PCR runs:
   - the forward primer and the reverse complement of the reverse primer are each matched with at most 3 substitutions;
   - a sequence with exactly one match of each gives the region between them, primers excluded, if it is longer than 200 bases;
   - each base of the region is replaced, with probability 0.02, by a base drawn uniformly from {A, C, G, T}, which may be the same base.
3. A genus is eligible if it has at least one region of at least 250 bases, the read length (ART makes no reads from a shorter amplicon). The eligible genera are E.
4. Genus g ∈ E gets k_g read pairs. k_g is N_g^α R / Σ_{h∈E} N_h^α, in exact rationals, rounded down; the remaining pairs go one each to the genera with the largest fractional parts, ties in document order. So Σ_g k_g = R exactly.
5. For each g with k_g > 0, ART runs on all of g's regions with MicrobeMixer's command, `art_illumina -ss MSv1 -amp -p -na -l 250 -f F` (amplicon mode, paired, MiSeq v1 profile, F read pairs per amplicon).
   - F = ⌈k_g / ℓ_g⌉ + 1, where ℓ_g is g's number of regions of at least 250 bases.
   - If ART gives fewer than k_g pairs, it is rerun with F doubled, up to 3 more times. A genus still short is removed from E and step 4 is redone; this did not happen on the toy.
6. Exactly k_g pairs are drawn uniformly without replacement from g's pairs. Since every long region gets F pairs, this is uniform over g's sequences with a long enough region.
7. Python's random generator is seeded with 1 (for the mutations and the draws). ART's generator is not seeded, as in MicrobeMixer.

*Credit schemes.* Each MEM M of either mate gives weights w_g(M) to genera. A read's weight for genus g is W_g = Σ_M w_g(M) over the MEMs of both mates, and the read is assigned the genus of largest W_g.
- Weights for genera labelled "uncultured" are discarded.
- Where a scheme normalises by a sum that is 0, the MEM gives no weight.

| scheme | w_g(M), for g ∈ D(M) |
|---|---|
| proportional | \|M\| c_g / C |
| even | \|M\| / \|D\| |
| normalised | \|M\| x_g / Σ_{h∈D} x_h, x_g = c_g / N_g |
| prior mixture (β) | \|M\| x_g / Σ_{h∈D} x_h, x_g = (c_g / N_g) · ((1 − β) N_g / N + β / G) |
| linear mixture (β) | \|M\| ((1 − β) c_g / C + β / \|D\|) |
| EM (t iterations) | \|M\| x_g / Σ_{h∈D} x_h, x_g = (c_g / N_g) · π_g^(t) |

with β ∈ {0, 0.1, 0.25, 0.5, 0.75, 1}. The endpoints are:
- prior mixture: β = 0 is proportional and β = 1 normalised;
- linear mixture: β = 0 is proportional and β = 1 even;
- EM with t = 0 is proportional.

The code computes each in a form that gives these endpoints bit for bit. The existing schemes are as before:
- **range:** |M| / (g_r − g_ℓ + 1) to every genus from g_ℓ to g_r, the leftmost and rightmost genera of D(M) in document (taxonomic) order.
- **one label:** all of |M| to one genus drawn with probability c_g / C.
- **sample:** independent uniform values for the genera of D(M) in document order; |M| / |D′| to each genus of D′, the left-to-right and right-to-left maxima.
- The random draws of the last two are seeded per block of reads.

*EM.* π^(0)_g = N_g / N, the database prior.
- **Each iteration:** from π^(t), compute w_g(M) as in the table for every MEM of every read (both mates). Let T_g be the sum of g's weights over all of them, and set π^(t+1)_g = T_g / Σ_h T_h.
- **Stopping:** after a fixed number of iterations, with no convergence test. Reads are assigned with π^(5) and with π^(20), reported as two schemes.
- **Weight per genus:** T_g is computed per genus, including genera labelled "uncultured".
- **Genera without credit:** a genus that gets no weight in some iteration has π_g = 0 from then on.

*Naive Bayes over the MEMs.* This treats a read's MEMs as independent observations, as the RDP classifier treats a read's 8-mers (Wang et al., 2007).
- For each genus g of ∪_M D(M), the score is S_g = log π_g + Σ_M λ_g(M) over the MEMs of both mates.
- λ_g(M) = log(c_g(M) / (2 N_g)) if c_g(M) > 0, the 2 because c_g counts both strands; otherwise λ_g(M) = log ε.
- The read is assigned the genus of largest S_g. The scores are turned into posteriors exp(S_g − max S), so ties and the accuracy rule are as for the other schemes.
- The code computes S_g in one pass over the listing: S_g = log π_g + Σ_{M ∋ g} log p_g(M) + (m − m_g) log ε, where m is the read's number of MEMs and m_g the number containing g.
- Three priors:
  - **reference (`bayes`):** π_g ∝ N_g^α, the sample's own, ε = 10⁻⁶;
  - **flat (`bayesflat`):** π_g = 1/G, for ε ∈ {10⁻², 10⁻³, 10⁻⁴, 10⁻⁶, 10⁻⁹};
  - **database (`bayesdb`):** π_g = N_g / N, ε = 10⁻⁶.
- The flat prior is the reference at α = 0 and the database prior is the reference at α = 1, bit for bit.

*The ε → 0 limit and the EM prior.* `--epsilon 0` is the limit: the genera missing the fewest MEMs, then the largest product of c_g / N_g. `--mode bayesem` starts from a flat prior and re-estimates it from the reads' posteriors, with the naive-Bayes likelihood at the first `--epsilon`, or the mixture likelihood below at `--mixture-eps`; a read's posterior is over the genera of its MEMs, and the prior over all genera. `bayes` and `bayesem` call only among scored genera, leaving out "uncultured".

*ε fitted to the reads (`bayesauto`).* The ε of the schemes above is not a probability, so it can't be fitted by likelihood: a larger ε raises every read's score. `bayesauto` uses the proper mixture model instead.
- **Model:** each MEM M of a read comes from the read's genus g with probability 1 − ε, and is noise otherwise. Its likelihood is (1 − ε) c_g(M) / (2 N_g) + ε C(M) / (2N), where C(M) / (2N) is the fraction of all database sequences (both strands) containing it.
- **The fit:** ε maximises Σ_reads log Σ_g (1/G) Π_M [(1 − ε) c_g / (2 N_g) + ε C / (2N)] over the grid 10⁻⁸, 3·10⁻⁸, …, 0.3, 0.5 (`--eps-grid`). The genera containing none of a read's MEMs enter as one lumped term.
- **Classification:** reads are then classified by S_g = Σ_M log[(1 − ε) c_g / (2 N_g) + ε C / (2N)], with a flat prior, over the genera of the read's MEMs.
- **Log:** the fitted ε and the log-likelihood at each grid point are written to the scoring log; `dial_tables.py` prints the ε per α and region.
- **On the toy:** it fits ε = 0.03 at α = 1, mostly 0.1 at α = 0, and 0.1–0.2 at α = −1, and gets 99.18%, 99.51% and 99.81%. That is against 99.19%, 99.46% and 99.55% for the flat prior with ε = 10⁻³ (the ε chosen on α = 0), and 99.35%, 99.46% and 99.77% for each α's own best ε. So it is as robust as ε = 10⁻³, or better.
- The prior is not fitted (flat); fitting it as well, by EM, is a possible extension.

*Naive Bayes in the tool.* `find --classify bayes --nseq FILE [--epsilon E]` writes, per read, the name and the top documents (0-based tags, comma-separated if tied; `-` without MEMs), with a flat prior.
- It accumulates in the same order and with the same arithmetic as `score_silva.py --mode bayes --prior flat`, so the calls are the same read by read: 28,800 calls on the toy (three values of α, two regions, two values of ε), checked by `tests/dial_toy.sh`.
- **Pairs:** `--mate2 FILE` classifies a pair over both mates' MEMs, mate 1's first, as the scorer does.
- **Scored genera:** `--exclude-docs FILE` keeps the listed documents from being called; a MEM still counts as missing for a genus that lacks it. `score_silva.py --exclude-out` writes the list, 373 documents on SILVA, so "scored" has one definition.
- **Checked on real SILVA** (2026-11-05): call for call the same as the scorer's `--calls` on all four aquatic regions at ε = 10⁻⁶ and 10⁻³, and on the dial's V4_V4 at α = 1 and −1, 200,000 pairs each. The exclusion changed 2,239 to 17,709 calls per set, so the check exercises it.
- `--epsilon 0`, the scorer's limit, isn't in `find`; it refuses a non-positive ε.

*Accuracy.* Strict genus accuracy of a scheme on a region is A/R, where A is the number of read pairs whose true genus alone has the largest W_g. In detail:
- a tie for the largest weight counts as wrong;
- so does a read pair with no weight at all (no MEM in either mate);
- a read whose true genus name is missing from SILVA's genus rank, or occurs at more than one node of it, counts as correct, as in Cliffy's scoring (`classify_silva_readset.py`). This applies only to reads from genera whose names are ambiguous in SILVA.

A table entry is the mean of the four regions' strict accuracies, each computed from its exact count (`strict_correct` / R), as in `silva_tables.py`.

**For the desktop, after its SILVA queue:**

```
cd ~/MEM-finding && git pull origin tool && make -C tool
S=~/MEM-finding/katka/silva
python3 tool/scripts/silva_alpha_fit.py $S/work/ref/rna/doc_to_traversal.txt ~/tool/silva/queue/nseq.txt \
    $S/work/reads/aquatic/*/*_seqtax.txt
tool/experiments/silva_queue.sh dial
tool/experiments/silva_queue.sh collect
python3 tool/scripts/dial_tables.py tool/results/silva_queue
```

(The every-50th set is a uniform sample of these reads, so the aquatic truth files give its composition too.) The α fit answers what the existing read sets correspond to. For each region it prints:
- the least-squares slope of log k_g on log N_g over the genera with reads, and Spearman's correlation;
- the fraction of SILVA's sequences in those genera;
- their ranks by N_g among all the reference's genera, genera without reads included (1 = largest, ties sharing the best rank): median, range, and how many are in the top 1% and 10%. Its output line per region goes into `tool/results/silva/README.md` with the tables.

**The dial's results (desktop, 2026-11-05; for the SILVA digest).** `results/silva/README.md` ("The α dial"), `results/silva/dial_tables.txt` and `.tex`. Strict genus accuracy, mean of the four regions, at α = 1 / 0 / −1:

| scheme | α = 1 | α = 0 | α = −1 |
|---|---:|---:|---:|
| range | 85.58 | 84.50 | 83.90 |
| sample (simulated cliffs) | 88.04 | 87.51 | 87.05 |
| normalised | 91.76 | 91.85 | 91.83 |
| proportional | 92.68 | 88.19 | 85.01 |
| naive Bayes, flat, ε = 10⁻⁶ | 92.66 | 94.23 | 94.78 |
| naive Bayes, flat, ε = 10⁻⁴ (held-out choice) | 84.49 | 95.61 | 96.18 |
| naive Bayes, sample's prior (reference) | 92.77 | 94.23 | 95.89 |
| naive Bayes, mixture, ε at the likelihood's peak | 91.42 | 94.26 | 95.19 |

- **Three things for the cloud session:**
  - **`bayesauto`'s fit fails on SILVA.** ε = 0.5, the top of the grid, in all 12 runs, with the log-likelihood still rising there; the toy fitted 0.03 to 0.2. The accuracy above is for ε = 0.5, so don't describe it as fitted.
  - **`scripts/simulate_dial.py` has no `__main__` guard.** Python 3.14's default start method, forkserver, re-imports the script, reruns `art(jobs)` at line 114, and the `Pool` at line 96 raises multiprocessing's bootstrapping `RuntimeError`. The first dial stopped there at α = 1, V1_V2, with nothing simulated. The rerun set the start method to `fork` from outside (`PYTHONPATH` with a `sitecustomize.py`), without changing the script. The script needs a guard or `multiprocessing.get_context("fork")`.
  - **`silva_queue.sh collect` copies 3.5 GB** (grammar binaries, `out.*`, tags) into `results/silva_queue`, and those can't be committed. The desktop deletes files over 1 MB after each collect. That also removes `dial/a*/<region>.seqtax.txt`, so `dial_tables.py tool/results/silva_queue`, as in the steps above, reports them missing and exits 1. It was run on `~/tool/silva/queue` instead.
- **The cloud session's answers (2026-11-05, fb7c06c and before):**
  - **`simulate_dial.py`:** its `Pool` now comes from `multiprocessing.get_context("fork")`, as `score_silva.py`'s does. Checked with the start method set to forkserver.
  - **`bayesauto`:** the fit isn't broken. The likelihood really is still rising at 0.5.
    - Near ε = 0, the log-likelihood falls by m · ln 10 per decade per read, where m is how many of a read's MEMs its most likely genus lacks. The SILVA logs give m = 4.6 to 6.5 per read pair, at every α and region; the toy gives 0.2 to 0.8.
    - So a large share of a read's MEMs are absent from its best genus. This probably comes from MicrobeMixer's 2% mutation of each region and ART's errors: a segment around a changed base can match another genus exactly over 30 or more bases, while the true genus lacks it.
    - The maximum-likelihood ε is therefore probably above 0.5. The grid now goes to 0.99, and the scoring log says when the fit is at its edge.
  - **`collect`:** it now copies only `.csv`, `.log`, `.txt` and `.tsv` files under 5 MB, and not the grammar step's `out.*`. The seqtax files (1.1 MB each) are kept, so `dial_tables.py tool/results/silva_queue` works.
  - **For the desktop, scoring only:**

    ```
    cd ~/MEM-finding && git pull origin tool
    rm -f ~/tool/silva/queue/dial/a*/credit/bayesauto.*
    tool/experiments/silva_queue.sh dial collect
    python3 tool/scripts/dial_tables.py tool/results/silva_queue
    ```

    The simulations and searches are done, so `dial` skips them and only rescores `bayesauto`: about 24 grid points × 12 read sets of likelihood passes, in the order of an hour with 12 threads. The ε fitted per region and α, and whether it is at the edge, are in `dial_tables.py`'s output.
- **The α fit doesn't place the real reads on the dial.** Its slope is 0.056 with Spearman about −0.04 (`results/silva/alpha_fit.txt`): abundance is unrelated to N_g. That is not the dial's α = 0, where every genus gets the same reads: the real sets are skewed, and proportional gets 93.77 on them against 88.19 at α = 0.

**The rescoring (desktop, 2026-11-05; for the SILVA digest).** `results/silva/README.md` ("Rescoring on the existing listings"), `rescore_tables.txt`, `timing_repeats.txt`.
- **`score_silva.py --mode bayes --epsilon 0`** is new: the ε → 0 limit, the genera missing the fewest MEMs, ties by the product of c_g / N_g. It now also chooses among scored clades only, since an unscored "uncultured" top could underflow the rest at tiny ε. ε = 10⁻⁹ reproduced the dial's row, and ε = 10⁻³⁰⁰ matched the limit to two decimals.
- **Flat naive Bayes falls steadily as ε shrinks below 10⁻⁶:** at α = 1 / 0 / −1, 92.66 / 94.23 / 94.78 at 10⁻⁶ and 89.68 / 90.80 / 91.48 at the limit. On the real aquatic reads it is 94.70, 94.36, 94.19 and 94.08 at 10⁻⁶, 10⁻⁹, 10⁻¹² and the limit, against proportional 93.77 and normalised 92.76.
- **`bayesauto`:** the likelihood peaks at ε = 0.5 in 9 runs and 0.6 in 3, inside the grid, as the cloud session predicted. Accuracy 91.42, 94.26, 95.19.
- **Excluding reads from genera with ambiguous names:** the aquatic table doesn't change, as it has none. The simulated sets have up to 9% such reads; both versions are in `rescore_tables.txt`.
- **Timing, 5 repeats:** 45.57 µs (45.46–45.87), 167,224 reads/s (164,575–167,504). These should replace 52.3 µs and 145,500 reads/s in the digest.
- **Timing for the digest:** one thread 45.57 µs (45.46–45.87), 46.88 µs total/reads, 21,332 reads/s (21,210–21,410); 12 threads 167,224 reads/s. The digest's 43.0 µs is V1_V2 mate 1 alone, not noise, so its timing-noise bullet needs rewriting: 52.3 µs is the outlier, and the read files differ by up to 22%.
- **Naive Bayes, added columns:** uniform 94.33 (ε = 10⁻⁶) and 94.32 (`bayesauto`); aquatic `bayesauto` 94.12, ε = 0.3 in V1_V2 and 0.5 elsewhere.
- **Index size, for the digest and the tool:** the digest's explanation of 1.44 GB is wrong. 1.56 GB is `.rlbwt` + `.csa` + `.tag`; 1.44 GB is RLCSA + tags + an 18 MB lookup table, with no RLBWT. `find --csa` appears not to use the RLBWT in the default configuration, so a `find` that didn't load it would reach about 1.42 GB with the lookup table. `sizes.txt` also misses the lookup table, which is built at load time. See `results/silva/README.md`, "Index size".
- **`simulate_dial.py`:** the desktop added no guard. The fork context in the script already fixes it: `tests/dial_toy.sh` passed here without the outside workaround.
- **For the cloud session:** `tests/test_score_schemes.py` gives 25 "reference" failures on the desktop with the toy's files and with the real reference, on fb7c06c's scorer as well as with the change above, while `dial_toy.sh` passes. So the change didn't cause them, but it hasn't been tested by that suite either. The failures start at "reference: freq".

**Time and memory on the desktop** (estimates; the toy runs in under a minute):
- **The simulation,** for α = 1, 0 and −1, each over 4 regions: each region costs about one MicrobeMixer run, mostly reading SILVA and extracting regions, with ART much faster. The uniform step's runs are the measure: see the timestamps of `~/tool/silva/queue/uniform/work/<region>/stdout.log`. Memory is MicrobeMixer's, within the 4.5 GB the uniform step ran in.
- **The searches:** 24 runs of 200,000 reads with 12 threads, a few seconds each.
- **The scoring** (genus level only):
  - Six single schemes at about 15 s each.
  - prior and mix, 6 β each, about 30 s each.
  - EM, 20 passes of about 3 s, so about 70 s.
  - So about 4 minutes per region and α, or about 45 minutes in all.
- **In all:** about 12 MicrobeMixer-sized simulations, plus under an hour (the three naive-Bayes files add three scorings per region and α, the flat one with five values of ε).

### Block length as genomes are added (2026-11-05)

The question: with the smallest successor, every new genome can only add successors to a phrase, so the smallest one is eventually a rare variant and block steps should fail at almost every phrase. Any rule that picks one successor per phrase works if each block symbol carries a bit ω (*More RLBWT Tricks* v2, the lemma on other continuations), and the most frequent successor should keep blocks long.

**The most frequent successor was already in the tool** (`index --twolevel --continuation frequent`, `src/twolevel_build.hpp`):
- Ties go to the smallest phrase id.
- ω = 1 when the next block's first phrase is larger than the successor of the block's last phrase; under the smallest successor it is always 1, and isn't stored.
- Symbols are ordered by δ, then ω = 0 by ascending j, then ω = 1 by descending j, as the lemma needs.
- `selftest.sh` checks it against plain search, with brute-force MEMs and intervals, on the chr19-like collection, as for the smallest successor. It passed again on 2026-11-05.

**Synthetic pangenomes: `scripts/pangenome_sim.py OUT --haplotypes N --queries Q --length L --seed S`.**
- **Model:** a coalescent with recombination (msprime 1.4.5, `pip install msprime`), with these parameters, all options of the script:
  - N_e = 20,000;
  - exponential growth to 10⁶ over the last 400 generations (recent human growth, which makes most variants rare and gives each new haplotype private ones);
  - μ = 1.25·10⁻⁸ and ρ = 10⁻⁸ per base per generation;
  - 10% of variant sites are indels of 1–10 bases (geometric, mean 2), the rest SNPs, on a random ancestral sequence.
- **Output:** `collection.txt`, the reference (sample 0) and then N haplotypes, so the first n lines are every smaller collection; `queries.txt`, Q held-out haplotypes; `stats.txt`.
- **With N = 1024 and L = 300 kb** (the run below):
  - 2,013 variant sites;
  - diversity 4.6·10⁻⁴ per base (humans: 6·10⁻⁴ to 10⁻³, with SVs);
  - singletons are 44% of the segregating sites;
  - each new haplotype adds about 24 new variants at N ≈ 16, 2.9 at 64, 2.1 at 256 and 0.9 at 1,024, and has about one private variant.
- It takes 2 s.

**The measurement: `experiments/block_scaling.sh DATA WORK "n1 n2 ..." [L]`.**
- **Input:** DATA has `collection.txt` and `queries.txt`, the generator's or chr19's. For each n, the collection is the first n sequences.
- **Indexes:** it builds the index (prep, BWT, `index --csa --twolevel`), or links to an existing one with `EXISTING=prefix%d`. It then rebuilds only `PREFIX.lvl` with `--continuation frequent`, under links to the same text, BWT and RLCSA.
- **Searches:** `find -L L --format mems --csa --twolevel --stats`, with and without block steps, single-threaded and pinned. Every output must be byte-identical, and with `CHECK_PLAIN=1`, equal to plain search.
- **Table** (also `WORK/blocks.tsv`), by n and rule:
  - the parse's phrases and distinct phrases;
  - blocks and phrases per block;
  - the fraction of phrase positions that are continuations;
  - the `.lvl` size;
  - character, parse and block steps per pattern character;
  - phrases crossed per block step, the mean block length as the patterns see it: parse steps without blocks minus parse steps with them, over block steps;
  - search times with and without block steps.

**Results on the synthetic pangenome** (`results/block_scaling/`). These are 300 kb haplotypes; the 20 held-out haplotypes are the patterns (6.0 million characters), with MEMs of at least 1000, k = 51 and s = 12. Every output was identical under both rules, with and without block steps, and equal to plain search.

| N | rule | phrases/block | continuations | phrases per block step | parse steps/char | block steps/char | time (s) |
|---|---|---|---|---|---|---|---|
| 16 | smallest | 62.9 | 98.4% | 51.5 | 0.0230 | 0.0014 | 0.17 |
| 16 | frequent | 184.9 | 99.5% | 116.7 | 0.0554 | 0.0004 | 0.17 |
| 64 | smallest | 43.0 | 97.7% | 39.7 | 0.0078 | 0.0018 | 0.14 |
| 64 | frequent | 179.5 | 99.4% | 135.9 | 0.0286 | 0.0004 | 0.15 |
| 256 | smallest | 27.1 | 96.3% | 26.3 | 0.0028 | 0.0028 | 0.13 |
| 256 | frequent | 171.8 | 99.4% | 147.9 | 0.0153 | 0.0004 | 0.14 |
| 1024 | smallest | 16.3 | 93.8% | 16.1 | 0.0010 | 0.0042 | 0.10 |
| 1024 | frequent | 169.9 | 99.4% | 157.3 | 0.0093 | 0.0004 | 0.13 |

- **The expectation holds for block length.**
  - Under the smallest successor, the mean block length falls steadily, from 63 phrases at N = 16 to 16 at N = 1024, by a factor of 0.6 to 0.7 with every fourfold increase in N. As the patterns see it, it falls from 52 to 16 phrases per block step.
  - Under the most frequent successor it stays flat, 185 to 170 phrases per block, and as the patterns see it, it even rises, from 117 to 157.
- **It doesn't make the search faster on these data.**
  - The most frequent successor takes 2.4 to 9 times as many parse steps per pattern character, so its parse and block steps together exceed the smallest successor's at every N: 0.0558 against 0.0244 at 16, and 0.0097 against 0.0052 at 1024.
  - These are probably the extra parse steps at the ends of matches and when backing up that the paper's §4 mentions for the earlier implementation; the counters don't separate them, so this is not checked.
  - The search times are the same within noise (0.10–0.18 s), because character steps dominate and are the same under both rules (0.0029–0.0245 per character).
  - So block length as genomes are added is not what limits two-level search here. Under the smallest successor, parse plus block steps per character even fall as N grows (0.0244 at 16, 0.0052 at 1024).
- **Caveat:** the synthetic haplotypes are short (300 kb), have no repeats or structural variants, and their patterns are held-out haplotypes from the same population. Real chr19 decides.

**Where the parse steps come from, and capped blocks (2026-11-05).**
- **Counters** (`find --stats`, a second line). Parse steps are split into:
  - **entry:** before the search's first block attempt, at the right end of what it searches;
  - **inside:** re-crossing a block phrase by phrase, after an attempt that failed or backed up, when another attempt followed, so the match went on past the segment;
  - **end:** the same, when the match ended in the segment, by a failed parse step or the left end of the search.

  Block attempts are reported with:
  - the attempts whose symbol was absent (not counted as steps);
  - failed block steps (an empty interval);
  - back-ups, and exact exits (Lemma "backup", first case).

  Without a cap, the steps and outputs are as before.
- **Uncapped, the most frequent successor's extra parse steps are all at match ends.** Per pattern character:

  | N | rule | parse | entry | inside | end |
  |---|---|---|---|---|---|
  | 16 | smallest | 0.0230 | 0.0142 | 0.00003 | 0.0087 |
  | 16 | frequent | 0.0554 | 0.0366 | 0.00004 | 0.0187 |
  | 1024 | smallest | 0.0010 | 0.0004 | 0.00003 | 0.0005 |
  | 1024 | frequent | 0.0093 | 0.0048 | 0.00002 | 0.0045 |

  Block attempts start only at a pattern's real block boundaries. With long blocks, the first boundary in from each end of a match is far away, and the search crosses up to it phrase by phrase.
- **`index --twolevel [--continuation ...] --max-block B`:** a block also ends after B phrases, under either rule.
  - **Ordering:** a block cut by the cap is followed by its own successor. Among the symbols with the same first phrase, its suffix therefore sorts after every ω = 0 symbol and before every ω = 1 symbol. It gets the key (B, +), with j = B larger than any real block's j ≤ B − 1. So the ordering lemma holds, with that one addition.
  - **Lemma 6 does not hold for cuts.** A pattern's block starts are text block starts in every occurrence that extends one phrase to the left only at real boundaries, where the next phrase is not the successor. A cut's position in the text depends on where its run of continuations started, and an occurrence that doesn't reach back to that start can be cut elsewhere.
  - **So the search keeps exactness differently.** Block steps never give false occurrences, only possible misses. After an inexact exit, the search backs up only to an exact interval: the attempt's own parse interval, or the block interval at a real boundary. The exact exit (the current phrase never continues) is unchanged. Without a cap every pattern block start is a real boundary, and this is the old rule.
  - **Format:** the cap is stored in the `.lvl` file under a new format version, KTLVL002, written only with a cap; KTLVL001 files load as before.
  - **Locating:** two-level locating (`--twolevel-locate`) is refused with a cap. A cut breaks the toehold's entry into the blocks: the first block row of a parse interval need no longer be its first row or a run start.
  - **Checks:** `katka-test mems PREFIX queries --twolevel --random 20` on the chr19-like test collection compares MEMs and intervals with plain search. It passes with both rules and B = 1, 2, 3, 4, 8 and 16, as do the uncapped ones; in the scaling runs every variant gave the same MEMs.
- **The scaling series with caps** (`results/block_scaling/table.txt`; most frequent successor with B = 8, 16, 32, 64). At N = 1024, per pattern character:

  | rule | phrases/block | parse | entry | inside | end | block | parse + block |
  |---|---|---|---|---|---|---|---|
  | smallest | 16.3 | 0.0010 | 0.0004 | 0.00003 | 0.0005 | 0.0042 | 0.0052 |
  | frequent | 169.9 | 0.0093 | 0.0048 | 0.00002 | 0.0045 | 0.0004 | 0.0097 |
  | frequent, B = 8 | 7.8 | 0.0059 | 0.0002 | 0.0055 | 0.0002 | 0.0163 | 0.0222 |
  | frequent, B = 16 | 15.3 | 0.0060 | 0.0003 | 0.0052 | 0.0005 | 0.0057 | 0.0117 |
  | frequent, B = 32 | 29.3 | 0.0062 | 0.0006 | 0.0047 | 0.0010 | 0.0023 | 0.0085 |
  | frequent, B = 64 | 53.2 | 0.0065 | 0.0012 | 0.0042 | 0.0011 | 0.0012 | 0.0077 |

  - A cap removes the entry and end steps, as intended, but brings as many re-crossing steps elsewhere ("inside"). The total of parse plus block steps stays above the smallest successor's at every N and every B; the same holds at N = 16 (0.0244 against 0.0413 for B = 64).
  - The likely mechanism is the one above. Near a match's left end, many occurrences sit in text runs that started earlier, with different cuts. Attempts at cuts then fail, or end at a cut that can't be trusted, and the search backs up and re-crosses the run phrase by phrase, attempting again at each cut.
  - The counter calls these steps "inside" because another attempt follows, though the match's end causes them. A counter that separates them was not added.
  - Search times are the same within noise (0.11–0.23 s), as before.

**Longer haplotypes (2026-11-05).** The run used:
- 2 Mb haplotypes from the same population model (`pangenome_sim.py --haplotypes 1024 --queries 20 --length 2000000`; diversity 5.2·10⁻⁴), at N = 64, 256 and 1024;
- 40 million pattern characters (20 held-out haplotypes);
- Big-BWT, run single-threaded: its multithreaded parser failed an internal assertion on these texts;
- `CLEAN=1` (`results/block_scaling/long2mb/`).

Every variant gave the same MEMs.

| N | rule | phrases/block | parse/char | block/char | entry | inside | end | time (s) | without block steps (s) |
|---|---|---|---|---|---|---|---|---|---|
| 64 | smallest | 40.4 | 0.0067 | 0.0024 | 0.0038 | 0.00007 | 0.0028 | 1.11 | 1.57 |
| 64 | frequent | 153.6 | 0.0281 | 0.0006 | 0.0162 | 0.00009 | 0.0118 | 1.17 | 1.37 |
| 64 | frequent, B = 64 | 52.0 | 0.0181 | 0.0017 | 0.0036 | 0.0117 | 0.0028 | 1.09 | 1.37 |
| 256 | smallest | 26.8 | 0.0024 | 0.0036 | 0.0013 | 0.00005 | 0.0011 | 0.98 | 1.63 |
| 256 | frequent | 152.1 | 0.0164 | 0.0006 | 0.0090 | 0.00005 | 0.0073 | 1.08 | 1.61 |
| 256 | frequent, B = 64 | 51.9 | 0.0102 | 0.0018 | 0.0017 | 0.0071 | 0.0014 | 0.92 | 1.61 |
| 1024 | smallest | 16.2 | 0.0011 | 0.0056 | 0.0006 | 0.00003 | 0.0005 | 0.99 | 1.66 |
| 1024 | frequent | 151.4 | 0.0107 | 0.0006 | 0.0058 | 0.00008 | 0.0048 | 0.98 | 1.81 |
| 1024 | frequent, B = 64 | 51.8 | 0.0065 | 0.0017 | 0.0011 | 0.0045 | 0.0009 | 1.09 | 1.81 |

- **The same pattern as with 300 kb haplotypes.** Blocks shrink under the smallest successor (40 → 16 phrases) and stay flat under the most frequent (154 → 151). The most frequent successor's extra parse steps are at match ends, and the cap trades them for re-crossing.
- **Block steps now matter in time:** turning them off costs 0.4 to 0.8 s of about 1 s.
- **The rules still don't differ in time.**
  - At N = 1024, turning block steps off adds about 3.35 million steps and 0.67 s, so a parse step costs about 0.2 µs.
  - The most frequent successor's extra parse and block steps over the smallest successor's come to about 184,000, or 0.04 s, which is within the noise of one run.
  - Steps of every kind account for well under half of each search's time. The rest is presumably parsing the patterns into phrases (closed syncmers, hashing, dictionary lookups over 40 million characters), which was not timed separately.
  - So on these data, the choice of continuation rule moves search time by a few percent at most, and the smallest successor takes the fewest steps.

**Entry against match-end parse steps, and a range start (2026-11-05; analysed, not implemented).**
- **The split** (2 Mb data, per pattern character; inside steps are negligible):
  - N = 64, smallest successor: 0.0038 entry against 0.0028 end; most frequent successor: 0.0162 against 0.0118;
  - N = 1024, smallest: 0.0006 against 0.0005; most frequent: 0.0058 against 0.0048.

  Entry steps are 55–58% of the steps at match ends, under both rules.
- **The proposal.** Start the search in Δ at the first block boundary i reached from the right. The pattern's phrases from i on are δ followed by x continuations, and their block-start occurrences are the Δ rows of the symbols (δ, j) with j ≥ x: under the ω order, the tail of the ω = 0 group plus the head of the ω = 1 group.
- **Why it is not exact as stated.** A search is for X[begin..end) with any end, and it starts with character steps for the partial last piece X[S[idx]..end). That piece runs past the last complete phrase by a tail X[S[idx]+K..end) of up to one phrase step. The parse interval at i is therefore for "δ, x continuations, tail", and the Δ range is a superset of it.
  - Blocks with at least x + 2 phrases continue with succ(last phrase). Whether the tail matches that phrase is one string comparison, so these are exact.
  - Blocks with exactly x + 1 phrases are followed by some other phrase y, and only those whose y starts with the syncmer and the tail match. Phrase ids are in lexicographic order, so within each symbol (δ, x, ω) these rows are probably one contiguous sub-range, ordered by y. Finding it needs y for a Δ row: a forward step (Ψ) on the block BWT, which is not stored, or a binary search over each symbol's rows.
- **So an exact range start needs more than a C-array lookup.** Either store enough to find that sub-range, or use the range only when no block with exactly x + 1 phrases can match: for example, when the tail is empty, or when the dictionary shows that no phrase other than succ starts with the syncmer and the tail. Not implemented: an inexact start would break the rule that every configuration gives exactly plain search's MEMs and intervals.

**Exact one-step match ends for blocks: the design, and what it needs (2026-11-05; not implemented).**
- **Search start (the right end).**
  - **One contiguous range:** under the ω order, the symbols of first phrase δ are (δ, b, 0), then every symbol with more than b continuations (both sides, and the capped key), then (δ, b, 1).
  - **The longer blocks:** these all continue with χ^{b+1}(δ), so they are kept or dropped together by comparing that phrase with the tail. That needs the phrase's text. The `.lvl` dictionary stores only hashes and ids, so it needs the parse oracle's 2-bit dictionary (`.ptxt`) loaded too.
  - **The exactly-b blocks:**
    - Within (δ, b, 0), and within (δ, b, 1), Δ's rows are sorted by the next block's suffix, hence by its first phrase y.
    - The ids of the phrases that start with the syncmer and the tail form one contiguous range, since ids are in lexicographic order. Finding that range is a binary search over phrase texts, so it also needs `.ptxt`.
    - The matching rows are one sub-range of the symbol's rows, found by a binary search that reads y through Ψ of Δ.
    - Ψ is a select on Δ's BWT. The block BWT (`EFRank`) supports rank but not select, so this needs a select over its per-symbol runs.
- **Match end (the left end).**
  - The rows whose preceding block ends with δ′ and b′ continuations are contiguous in the co-lex order of blocks. Counting them needs a wavelet tree over Δ's BWT relabelled by co-lex rank: a new structure, about 23 M blocks × 21 bits, roughly 60 MB at 16 chr19 individuals.
  - Turning that count into the exact parse interval also needs the occurrences where the run of continuations goes on past b′, so that the match end falls inside a text block. I have not worked out how to count those without stepping through the parse; this is the open part.
- **Status:** not started in code. The three pieces it needs are a phrase-text range lookup (`.ptxt`), Ψ on Δ (select on the block BWT), and the co-lex wavelet tree, plus the counting argument for runs that continue past b′. The rule stands: exact intervals or nothing, so no approximate version was tried.
- **What it could save:** the entry and end steps measured above, 55–58% and 42–45% of the steps at match ends. On the 2 Mb data that is 0.0011 (smallest successor) to 0.0106 (most frequent successor) parse steps per character at N = 1024. At about 0.2 µs per step that is under 0.1 s of a 1 s search, so a few percent on these data. It could be more where matches are shorter and match ends more frequent.

**For the desktop: chr19 with 16, 32 and 64 individuals.** This reuses the existing indexes (`~/tool/chr19/c33`, `c65`, `c129`, built with `index --preset chr19`) and the patterns from `~/shortcuts/chr19_scaling/queries.txt`:

```
cd ~/MEM-finding && git pull origin tool && make -C tool
EXISTING=$HOME/tool/chr19/c%d PIN=9 tool/experiments/block_scaling.sh ~/shortcuts/chr19_scaling ~/tool/chr19/blocks "33 65 129" 1000
cp ~/tool/chr19/blocks/blocks.tsv tool/results/block_scaling/chr19_blocks.tsv     # and commit, with the printed table
```

- **What it builds:** for each n, the most frequent successor's `.lvl` (under `~/tool/chr19/blocks/c<n>f`, linked to the existing text, BWT and RLCSA), and the smallest successor's `.lvl` once more, under a probe prefix, only for its statistics; the probe is deleted.
- **Cost:** each two-level build is the two-level step of the original build: about 200 s and 4.6 GB at 33, 360 s and 5.2 GB at 65, and 700 s and 8.3 GB at 129. Then four searches per n of 5–10 s each.
- **In all:** about 45 minutes, within the 11 GB cap.
- **If `c129` was deleted for disk space:** put a link `collection.txt -> hprc_chr19` in `~/shortcuts/chr19_scaling` and run without `EXISTING` for 129, with `BIGBWT=~/rz/MEM-finding/katka/deps/Big-BWT/bigbwt`. It needs about 40 GB of disk, as before.
- **Optional:** `L = 100` instead of 1000 is the paper's other threshold.
- **Capped blocks as well:** add `CAPS="16 64"`, which builds the most frequent successor with blocks of at most 16 and 64 phrases. Each cap costs one more two-level build per n.

**Results on chr19 (desktop, 2026-11-05; `results/block_scaling/chr19_table.txt`, `chr19_blocks.tsv`).** 16, 32 and 64 individuals, L = 1000, 122.4 million pattern characters, k = 51, s = 12. Every output was byte-identical across the two rules, with and without block steps. Plain search was not rerun here (`CHECK_PLAIN` unset); `results/chr19/` has that check for the same indexes.

| n | rule | phrases/block | continuations | phrases per block step | parse/char | block/char | parse + block | time (s) | no block steps (s) |
|---|---|---|---|---|---|---|---|---|---|
| 33 | smallest | 9.52 | 89.5% | 6.93 | 0.0208 | 0.0148 | 0.0356 | 5.64 | 8.64 |
| 33 | frequent | 30.59 | 96.7% | 21.54 | 0.0435 | 0.0037 | 0.0472 | 5.92 | 8.41 |
| 65 | smallest | 8.32 | 88.0% | 6.95 | 0.0128 | 0.0174 | 0.0302 | 5.52 | 8.98 |
| 65 | frequent | 30.10 | 96.7% | 24.29 | 0.0364 | 0.0040 | 0.0404 | 5.63 | 8.81 |
| 129 | smallest | 7.15 | 86.0% | 6.37 | 0.0090 | 0.0208 | 0.0298 | 5.71 | 9.50 |
| 129 | frequent | 29.67 | 96.6% | 25.51 | 0.0333 | 0.0042 | 0.0375 | 5.71 | 9.60 |

- **As on the synthetic pangenomes, the expectation holds for block length but search doesn't get faster.** Under the smallest successor, phrases per block fall from 9.52 to 7.15 as the individuals quadruple. Under the most frequent they stay at about 30, and as the patterns see them they rise from 21.5 to 25.5.
- **The fall on chr19 is mild.** As the patterns see it, the smallest successor's block length barely moves (6.93, 6.95, 6.37).
- **The most frequent successor takes more steps at every n,** in parse plus block steps, and its extra parse steps are again at the ends of searches (entry and end, in the table's second part). Times are the same within noise, and block steps save 30–40% of the time under both rules.
- **The c129 run:** `c129.S` and `c129.bwt` had been deleted for disk space, so they were rebuilt on the external drive and linked (`~/ext/chr19_c129`; prep and Big-BWT took 33 minutes, the `.tbl` and `.info` match). `c129.rlbwt`, which `index` rewrites from the BWT, kept the same SHA-256, so the rebuilt BWT gives the same RLBWT as the original. The `.lvl` builds at 129 took about 30 minutes each against the estimated 12, reading from that drive. DATA was `~/shortcuts/chr19_scaling/hprc_chr19`, so no link was needed.

### Cliffy (2026-11-04)

Cliffy, github.com/oma219/cliffy at `3763529`, built here with `katka/silva/build_cliffy.sh`'s steps (`cliffy_rz.patch`, `gxx-compat`) under GCC 13 and CMake 3.28. Indexed: 20 synthetic genomes as 20 documents, with `--taxcomp --num-col 7 --two-pass`, with and without `--revcomp`.

**Why our sizes of its SILVA indexes are about twice those in Cliffy's paper (9.01 and 4.27 GB).**
- **The paper counts 10 files** (cliffy-experiments, `workflow/rules/exp_1.smk`, rule `combine_all_index_sizes_into_csv_file_exp1`): `.bwt.cliffy`, `.F.cliffy`, `.runcnt`, `.ftab`, and the six profile files `.taxcomp.{sdap,edap,of.sdap,of.edap,ofptr.sdap,ofptr.edap}`. It reports the sum in GiB (2³⁰ bytes).
- **We counted every file the build writes** (`du -cb output*` in `katka/silva/compare_lib.sh`), in GB. That adds the text (`.fna`), the parse and dictionary (`.parse`, `.parse_old`, `.dict`), LCP and SA samples (`.lcp`, `.esa`, `.ssa`), the BWT's run heads and lengths (`.bwt.heads`, `.bwt.len`), `.occ`, `.last`, and a text dump of the ftab (`.ftab.debug.txt`).
  - On the small index these extras make the total 1.82 times the paper's count with both strands, and 1.94 times without.
  - At SILVA's scale the extras grow with the text and the number of BWT runs r. On the small index `.lcp` took 19 bytes per run and `.esa` and `.ssa` 10 each. With r of about 70 million and the text of both strands, they and `.fna` alone add roughly 5 GB to the paper's 9.67 GB (9.01 GiB). The other extras make up the rest of our 19.58 GB; the script below lists them on the server's index. GiB against GB adds another 7%.
- **Both strands are not the difference:** Cliffy's experiments build with `--revcomp` too (`exp_1.smk`), as we did.
- **On disk against in memory:**
  - The query (`cliffy run`, traced with strace) opens the paper's 10 files plus `.fna`, `.bwt.heads` and `.bwt.len`.
  - On the small index the files it opens take 67 MB, while its peak memory was 242 MB. Most of the latter is the fixed 10-mer ftab and the run's buffers, so neither number scales simply to SILVA.
  - The fair on-disk number is the paper's 10 files. For memory, measure the query's peak.
- **`scripts/cliffy_index_size.sh PREFIX [CLIFFY READS.fq [RUN OPTIONS]]`** prints all three counts for any Cliffy index: the paper's (GiB and GB), all files (with the extras listed, largest first), and the files `cliffy run` opens. **For Gonzalo:** run it on the server's two SILVA indexes. The first line should give about 9.01 and 4.27 GiB; our 19.58 GB is the second line.
  ```
  tool/scripts/cliffy_index_size.sh katka/silva/work/index/full_text/output
  tool/scripts/cliffy_index_size.sh katka/silva/work/index/minimizers/output
  ```

**Scoring Cliffy's approximate listing with one command.**
- **The query:** `cliffy run --taxcomp --ftab` lists, for each match, `[start,end] {d1,…}`: the documents of its approximate listing. Cliffy's paper scores that listing with "method2" of `classify_silva_readset.py`:
  - each match adds its full length to the clade of every listed document;
  - a read is correct if its clade is among those of maximum weight (ties count as correct, and so does a read with no match unless its clade is "uncultured").
- **The scorer:** `scripts/cliffy_listing_score.py` is a fast equivalent of method2, as `katka/silva/rz_score2.py` is of method0 (the LCA query). On a toy data set (24 clades with an "uncultured" one, 20 documents, 3,000 read pairs with ties, unmatched reads and repeated documents) its `theirs=` accuracy equals their script's method2 at all six levels.
- **Correction:** `rz_score2.py`'s comment calls its `--mode list` "the method2 weighting", but that mode gives length/k to each of k listed documents, which is our even credit, not method2. Use `cliffy_listing_score.py` for Cliffy's listing.
- **For Gonzalo**, on the server where `katka/silva/run_cliffy.sh` ran:
  ```
  EXP=~/MEM-finding/katka/cliffy-experiments/src tool/scripts/score_cliffy_listing.sh
  ```
  - It finds the newest `katka/silva/work/compare/cliffy-*`. For both indexes and all four regions it scores `runs/<index>.<region>_{1,2}.p1.all` into `scores/<index>.<region>.listing.csv`, and prints genus-level accuracy per region and the mean.
  - Each csv has all six levels, with `theirs=` (Cliffy's rule), the accuracy counting unmatched reads as wrong, and `strict` (ties counted as wrong).
  - It was tested on the toy data laid out as on the server.

### Reproducing More RLBWT Tricks v2 (2026-11-04)

`CHR19.md` is the guide for others (Eddie Ferro's group): the commit, the build, the data, the commands for Table 1, Table 2 and the locating paragraph, their expected outputs, and using two-level search with blocks from code.
- **The driver:** `experiments/chr19_paper.sh` (`build`, `table1`, `table2`, `locating`, `all`). It checks every output against plain search and writes the tables in the paper's layout.
- **The data:** `experiments/hprc_chr19.py`. It reuses `chr19/twolevel/download_hprc_chr19.py`'s download code for the 16 individuals, or takes a manifest. **For the desktop:** write `experiments/hprc_chr19_64.tsv` from the sources of `~/shortcuts/chr19_scaling/hprc_chr19`: role, sample, haplotype, sequence name, length, number of N, for CHM13, the 128 haplotypes in collection order, and the two queries. Then `CHR19.md`'s data section works for 32 and 64 individuals.
- **The API example:** `experiments/twolevel_example.cpp`. Its MEMs and intervals are byte-identical to `katka find`'s.
- **The synthetic check:** `experiments/chr19_synthetic.py`.

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
- `selftest.sh` checks `find --twolevel --locate` against plain search: the same MEMs and intervals, and occurrence counts equal to ep − sp. It also checks that its `--locate-out` positions equal STPD over the parse's, MEM by MEM.
- `find --locate-out FILE` (with any `--locate`) writes, per MEM, the name, start and end, then its sorted positions (those of its reverse complement in S, for BML and FB). On the synthetic collection, two-level locating, STPD over the parse and character-level STPD write byte-identical files: 4,410 positions at L = 1000 and 29,460 at L = 100. Every one of them holds the MEM's reverse complement in the text.

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

`index --twolevel-locate` needs `c33.bwt`, `c33.lvl` and `c33.ptxt`; the last is rebuilt if it is missing or of the old format. Check each output against plain search (`plain$L.33.tsv` from `run_n.sh`): columns 1–5 must be equal, and column 6 (occurrences) must equal ep − sp.

That checks only the number of occurrences. Also compare the positions themselves: add `--locate-out tlpos$L.33.txt` to the two-level run, and run STPD over the parse with `--csa --stpd-parse --locate --locate-out sppos$L.33.txt`. The two files must be byte-identical: the same MEMs in the same order, each with its sorted positions. Report:
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
9. **Grammar.** The built-in RePair is approximate: each round replaces every pair whose frequency is at least half the largest. The grammar differs from bigrepair's, so the size reported in the paper (0.23 GB for SILVA) should be remeasured, or the tool given bigrepair's output with `--bigrepair`. Answers are identical either way (tested). The built-in one was estimated here at 16 bytes per run. Measured since ("SILVA queue: runner and tables", Grammar tags), it takes 25 bytes per run, about 9.2 GB for SILVA's 359 million runs.
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
