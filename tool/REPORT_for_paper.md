# Report for the paper task: KATKA as a tool

From the session that built `tool/` (branch `tool`), 2026-11-01, for the task that edits *Three Examples of MEM-Finding* for the CPM second round (deadline about 13 February 2027). It says what the tool is, what changed compared with the code behind the current draft, what has been verified, and which numbers and sentences in the draft will need revisiting. It proposes no text for the paper: the paper is edited elsewhere, and every number below comes from small synthetic data unless stated otherwise.

The November submission uses the old code and its numbers. Nothing under `katka/`, `katka-readable/`, `chr19/` or `paper/` was changed.

## 1. What the tool is

One program, `katka` (`tool/src/`, about 3,000 lines of C++17). It needs only zlib; libsais is bundled, and sdsl-lite is no longer needed. It has one index of one text and one search engine, with an option for each technique in the paper and a preset for each use case:

| preset | use case | settings |
|---|---|---|
| `zymo` | noisy long reads, few taxa (§3) | homopolymer compression; θ = 14; adaptive L from 20; KeBaB pseudo-MEMs with k = 20, longest first; species tags; proportional credit; explicit RLCSA |
| `silva` | short reads, many taxa (§4) | L = 30; explicit RLCSA; 10-mer lookup table; run-length tags; proportional credit |
| `chr19` | long accurate matches (§5) | BML with L = 1000 (or 100); two-level indexing with blocks; no tags |
| `novel` | a few long MEMs in mostly novel sequence | fixed L = 1000; the inline KeBaB window check (new, §3.2 below) |

The commands are `katka prep` (the text of both strands, optionally homopolymer-compressed or digested), `katka bwt` (libsais, or external Big-BWT/pfp-merge files), `katka index` and `katka find`. The output is KATKA's listing format (Cliffy's, extended with counts), so the existing scoring scripts work unchanged. `tool/README.md` lists every option and maps paper sections to source files.

**Techniques in the tool, by section of the current draft:**
- **§2:**
  - the RLBWT; explicit and Elias–Fano RLCSAs; the k-mer lookup table;
  - both strands in one text;
  - BML with fixed L, adaptive L, or the t longest MEMs; forward–backward search;
  - KeBaB pseudo-MEMs.
- **§3:** homopolymer compression of text and reads, with spans mapped back to the read; the low-complexity filter θ; KeBaB longest first.
- **§4:**
  - the complete run-length tag array, built without a suffix array by an LF pass and φ streaming;
  - counting; Muthukrishnan–Sadakane listing without the previous-run array; range/LCA answers;
  - sampled tags via LF; grammar-compressed tag runs; threads.
- **§5:** two-level indexing at closed syncmers, and blocks of continuations with the smallest or the most frequent successor (the order bit ω).
- **App. B:** symmetric minimizer digests and BML on digests.

**Left out, as decided:** the rz-index, shortcuts, the frequency-rank encoding of the parse, RLZ verification of digested matches, whole-phrase matching, and the combined LCA-and-listing rule. If the resubmission presents KATKA as a tool, these parts of the paper describe things the tool does not do.

## 2. Changes that affect what the paper says

### 2.1 Section 5 now uses both strands

Decision 2 of the brief: the chromosome 19 experiment now indexes one text, S, containing every sequence and its reverse complement, as in §§3–4, instead of T and the reversed collection T′. Sentences that will no longer describe the tool:

- `prelim.tex:30–34` ("Both strands and the reversed collection"): the chr19 case with T and T′.
- `prelim.tex:50`: "(or in T′, or as a reverse complement)".
- `chr19.tex:20–21, 93, 225, 234–235, 302, 340, 421, 431`: T, T′ and T\*, "reverses" of phrases, the dictionary "shared by T and T′", and sizes "of T and T′".
- `chr19.tex:573`: "index of T and T′, which finds only forward-strand MEMs". The tool now answers the same question as ropebwt3 and Movi 2 (MEMs on either strand), so the strand-check paragraph comparing outputs (`chr19.tex` around lines 560–575) changes character. The comparisons become like-for-like, but the reverse-strand differences it explains no longer arise from our side.

Consequences to state or re-measure:
- MEMs are with respect to both strands. A reverse-strand match can lengthen a MEM or add one.
- A MEM's BWT interval is in the BWT of S, and its size counts occurrences on both strands.
- n doubles: 4.07 × 10⁹ for 16 individuals.

On synthetic data without long reverse-strand matches, the tool reports the same MEMs, with the same occurrence counts, as the old `query_twolevel` for BML and forward–backward search with L = 100 and 1000 (`tool/tests/regress_chr19.sh`).

### 2.2 The syncmer order is canonical under reverse complement

The draft (`chr19.tex` around lines 200–204) defines g(γ) as a hash of the smaller of the codes of γ and of its reverse, so that g(γ) = g(γ^rev). The tool uses the smaller of the codes of γ and of its reverse complement (a splitmix64 hash; s ≤ 32), so that g(γ) = g(rc(γ)). Both of the brief's checks hold:

- **Prefix-freeness and Lemma 6 (order) are unchanged.** Their proofs use only that whether a k-mer is a closed syncmer does not depend on its context, which holds for any order g on s-mers. The lemmas "window" and "bounce" and Theorem "twolevel" also make no use of the symmetry of g.
- **One dictionary serves both strands.** The s-mers of rc(κ) are the reverse complements of those of κ, in reverse order, with the same g values. So κ is a closed syncmer if and only if rc(κ) is, with first and last swapped. The syncmer starts of rc(P) are |P| − k − q for the starts q of P, and the phrases of rc(P) are the reverse complements of P's, in reverse order. This replaces the sentence "By the symmetry of g, the boundaries of T′ mirror those of T, and the phrases of T′ are the reverses of those of T, except at the ends" (`chr19.tex:234–236`).
- **The dictionary** stores, for each phrase of S over A/C/G/T, its id and the id of its reverse complement. One lookup per phrase therefore serves backward search on P and forward search (backward search on rc(P)). At construction, every phrase over A/C/G/T has its reverse complement in S, except possibly the first phrase of S, which starts at position 0 rather than at a syncmer and so equals no pattern phrase. Phrases containing a separator or $ are left out of the dictionary.

The blocks lemmas (blockorder, blockstart, backup, other continuations) concern the parse alone and are unaffected.

### 2.3 One BML variant everywhere

The tool uses KATKA's BML for every section: after a window matches, extend left to a, then forward to e, and continue at x = e − L + 1. The chr19 code used Gagie's variant: forward from x, backward from e, and continue at the start of the longest suffix of P[0..e+1). Both report exactly the MEMs of length at least L (verified by brute force), but they take different numbers of steps.

On a synthetic collection of 16 mutated copies of 2 Mbp (64 × 10⁶ characters with both strands), BML with L = 1000 took 7.3 million backward steps in the tool, against 11.5 million in the old code. So the "3.7 backward steps per pattern character" and the steps-per-character figures in `chr19.tex` ("Results") will change.

With two-level search, the left extension after a window repeats a search from the window's end (`segment_search::extend_left`), which costs extra parse steps.

### 2.4 Implementation details the draft states

- **RLBWT.** The draft (`app-katka-details.tex:58`, `prelim.tex:25`, `chr19.tex:545`) says the run heads are in a Huffman-shaped wavelet tree (KATKA), or a wavelet tree with rank9-style bitvectors (chr19). The tool stores them in a wavelet matrix over the compacted alphabet, with rank9-style counts, and does one traversal per rank. The run starts and per-character cumulative lengths are Elias–Fano coded, as before.
- **RMQ for listing.** KATKA used sdsl's `rmq_succinct_sct`. The tool uses its own structure, which also does not store the previous-run array: the balanced parentheses of the tree in which each element's parent is the last earlier element with a smaller value, and a range-min-max search over the excess. It takes about 6 bits per run: 6.1 measured on a tag array of 10⁶ runs, of which 2 are the parentheses (`tool/src/rmq.hpp`). Listing answers are identical to KATKA's. Any size quoted for the RMQ (e.g. `app-katka-digests.tex:50`: 0.03 GB) is for the old structure.
- **Grammar.** The draft (`silva.tex:243–244`) uses bigrepair's RePair, giving 0.23 GB for SILVA's runs' taxa. The tool's built-in RePair is approximate: in rounds, it replaces every pair at least half as frequent as the most frequent. It can also read bigrepair output (`--bigrepair`). Answers are identical either way, but the size with the built-in grammar has not been measured on SILVA. The encoding of the grammar is unchanged.
- **Tags** may be any 32-bit integers; KATKA stored 14 bits.
- **Digests.** Reads now collapse equal neighbouring codes, as the reference text does; KATKA collapsed equal neighbouring bytes in reads. The two differ only when two different codes that share the merged byte are adjacent, which with k = 4 means two of the 6 rarest codes. In that case KATKA's read digest could miss a match. The digested text, byte map and listings were identical to KATKA's original code on test data that includes merged codes. The draft's description (App. A "Digests", App. B) remains accurate if "equal consecutive symbols" is read as equal consecutive minimizers.

### 2.5 Numbers that come from the old code

None of the following has been re-measured with the tool. The full runs are planned on the desktop (`tool/NOTES.md`, "Commands for the desktop").

- **§3 (Zymo):** all times and sizes.
- **§4 (SILVA):** all times and sizes (Table "Time and space"), the tag-array construction time and memory (`silva.tex:216`), and the grammar sizes and times.
- **§5:** Table `tab:main` (MiB, BML and FB times), the parse-index sizes (102.5 and 139.6 MiB), the step counts, the build times and memories, and the dictionary size (46 MiB).
- **App. B:** Table `tab:digest`.

Accuracy should not change for §§3–4. The tool's listings are byte-identical to KATKA's for the configurations tested: counting, listing and range answers; sampled, per-species and grammar tags; `-X`, `-K`, `-D` and `-U`; a per-dataset BWT; and threads. Accuracy figures therefore carry over, provided the same reads and references are used.

## 3. New in the tool and not in the draft

If the resubmission presents these, here is what they do and the correctness arguments, which are checked by brute force.

### 3.1 BML for the t longest MEMs

Keep the t longest MEMs found so far. Once there are t, set L to one more than the shortest of them. Raising L only stops the search for MEMs no longer than the shortest kept one, so at the end the kept MEMs are the t longest. When several MEMs tie for t-th place, the first one found is kept. Adaptive L (keep L equal to the longest so far) is the case that reports all ties of the longest.

### 3.2 The inline KeBaB window check

This is for a fixed, large L on queries mostly absent from the text. Before taking backward steps for the window P[x..x+L), check its k-mers against the Bloom filter. If P[j..j+k) fails, with j the rightmost failing k-mer start in the window, jump to x = j + 1.

Correctness: the filter has no false negatives, so P[j..j+k) does not occur, and no MEM contains it. Any MEM of length at least L starting at some s in [x, j] covers [s, s+L), which contains [j, j+k) because j + k ≤ x + L ≤ s + L.

Implementation (as the paper task asked, 2026-11-02): the k-mers are tested lazily, from the window's right end leftward. The rolling forward and reverse-complement codes are restarted at the window's last k-mer, which takes k steps, and are then updated one character to the left per test. Testing stops at the first failure, which is the rightmost failing k-mer. Nothing is precomputed over the pattern. `--stats` reports the numbers of k-mer tests and restarts. On 50 random 5,000-base queries with L = 200 and k = 20, there were 1,350 restarts and 1,475 k-mer tests, with no backward steps.

On the tool's tests, the output is exactly the MEMs of length at least L, as for plain BML.

### 3.3 Tags with two-level search

Two-level search returns the same intervals in the BWT of S as plain backward search, so tag answers carry over. Listings with tags were byte-identical with and without two-level search, with and without block jumps, and with grammar-compressed tags. This combination had never been run before.

### 3.4 STPD search

This is a backend from *Compressing Suffix Trees by Path Decompositions* (Becker et al., ICALP 2026; `src/stpd.hpp`), added 2026-11-02. It is the mirror image of their Algorithm 1: it matches right to left with π = ISA, i.e. their st-colex⁻ of the reversed text.

- **Sample:** the sample positions are BWT run starts. A run start x is needed iff lcp(SA[y], SA[x]) < lcp(SA[x−1], SA[x]), where y is the previous row with the same BWT character.
- **Construction:** the sample is built from the r-index's LF pass plus text comparisons, without a suffix tree.
- **Size:** it is exactly 2 smaller than the paper's |st-colex⁻| computed by definition on every test collection, the 2 being boundary entries. It holds about 0.66 r positions there.
- **Against the paper's code:** it is, position for position, the st-colex⁻ sampling computed by the paper's code (`stpd-sampling -c` of github.com/regindex/STPD-index) on the reversed text, except for their position 0. On their yeast text (12.2 million characters) both have 5,479,329 positions (0.654 r), apart from that position 0. Ours took 4.7 s to build, from the BWT runs, against 28 s for theirs, which builds the suffix tree. Both timings are from the cloud container.
- **Correctness:** it is checked to give the same MEMs and intervals as plain search, with every search ending at the lexicographically smallest occurrence.
- **Speed:** on the synthetic chromosome-like collection, BML with L = 1000 took 0.20 s without intervals and 1.2 s with them, against 2.1 s for plain search and 0.40 s for two-level search with blocks.

### 3.5 Items for the draft's TODOs (2026-11-02)

- **Lines 333–334 (intervals of STPD MEMs by two-level search):** implemented as `find --stpd --twolevel`, with identical output. On the synthetic chromosome-like collection (L = 1000), it takes 0.40 s, against 0.98 s with plain backward search for the intervals and 0.34 s for two-level search alone.
- **Line 390 (B_λ):** `index --tags --tag-lcp λ1,…` reports B_λ, exact against brute force, with Lemma 11's bound. It still needs running on the Zymo and SILVA references.
- **Lines 192–193 (KeBaB experiment):** `tests/exp_kebab.sh` runs KeBaB then BML, plain BML, and KeBaB inside BML. On small synthetic data, inside BML tests 15 times fewer k-mers than the separate pass but is not faster than plain BML (0.110 s against 0.100 s), because a 64-Mb index makes failed windows cheap. The experiment needs a large index on the desktop.

### 3.6 BWTs

`katka bwt` builds BWT(S$) with libsais. `katka bwt -d` builds a per-dataset BWT like pfp-merge's, for small data, from the suffix array of the datasets with distinct terminators. That makes the whole pipeline runnable without Big-BWT or pfp-merge on small data.

## 4. What has been verified (and may be stated as checked)

All of this is on small data (synthetic: 20 genomes of 4 species, 300 reads with errors and reverse complements, and a chromosome-like collection of 8 haplotypes). `tool/selftest.sh` runs it all in about 2.5 minutes.

- **Brute force against a suffix array of S$:**
  - BWT characters and ranks;
  - all MEMs of length ≥ L (L = 1, 12, 25, 40) by BML and by forward–backward search, with their BWT intervals;
  - adaptive L (every reported MEM is a MEM, and all the longest are reported);
  - the t longest (t = 1, 3, 10);
  - KeBaB pseudo-MEMs and the inline check;
  - counting and listing with full, sampled (s = 4), grammar-compressed and per-species tag arrays;
  - homopolymer compression;
  - BML on digests: every digested MEM with its span, for L = 11, 15, 20 and 30.

  These cover about 240,000 MEMs on reads plus random patterns from both strands.
- **Every search configuration gives identical output:** RLBWT; explicit and Elias–Fano RLCSA; lookup table; two-level without and with block jumps, with either continuation rule. Two-level search is also checked query by query against plain backward search, on both strands, at k = 51, s = 12 and at k = 21, s = 8.
- **Regressions against the old code:**
  - KATKA (`katka-readable`): byte-identical text, table, KeBaB filter and listings for 12 option combinations on a plain, a homopolymer-compressed and a per-dataset text.
  - KATKA's original digest code (`katka/rz-index`): byte-identical digested text, byte map and listings, for L = 11 to 30, with and without a lookup table.
  - `chr19/twolevel`: identical MEMs and occurrence counts.
- **STPD:** the same MEMs and intervals as plain search, and the sample's size against the paper's definition (§3.4).
- **Sabotage tests:** deliberately broken jumps in BML, in block back-ups and in the digest failure jump were each detected by the tests.

Not verified: anything at the scale of the paper's data. KATKA's 40-genome `test-small` collection is not in the repository, so the KATKA regression ran only on synthetic data.

## 5. Indicative measurements here (not for the paper)

These come from a cloud container, not from the desktop the paper uses.

- **Speed against `katka-classify`:** on 12,000 short reads, with L = 30 and counting, the tool's RLBWT path takes 3.8 s against 4.3 s. Its explicit-RLCSA path, the same algorithm, is equal to up to 15% slower, within the noise here.
- **Chromosome-like collection** (16 copies of 2 Mbp, both strands; two 2-Mbp queries, L = 1000):

  | method | tool | old code (T and T′) |
  |---|---|---|
  | plain search | 1.78 s | 6.87 s |
  | two-level | 0.39 s | 0.77 s |
  | two-level with blocks | 0.38 s | 0.72 s |

  Output was identical across methods.
- **Index sizes there:**
  - the RLBWT of S takes 8.7 MB, against 16.0 MB for the old RLBWTs of T and T′. The encodings differ, so this is not a like-for-like comparison;
  - the parse index with blocks takes 9.6 MB, against 6.9 MB for the old two levels.

## 6. Decisions still open (Travis's)

1. Whether to add `test-small` to the repository, so the KATKA regression can run on real genomes.
2. The `novel` preset's values: L = 1000 and k = 20 were chosen here; the brief says only "fixed large L".
3. The `zymo` preset includes θ = 14 for the longest-MEM runs, as the brief's table says. The `-X` command line in `katka-readable/README.md` does not use it.
4. Whether the paper keeps the T/T′ description of §5 for the November version only, and switches to both strands for February.

## 7. Where to look

- `tool/README.md`: commands, options, presets, file formats, and the map from paper sections to source files.
- `tool/NOTES.md`: status, all the differences from the old code (in more detail than §2 here), and the desktop commands to rebuild the SILVA, Zymo and chromosome 19 indexes and rerun the experiments.
- `tool/selftest.sh` and `tool/tests/`: the tests, the brute-force checker (`test.cpp`) and the three regression scripts.
