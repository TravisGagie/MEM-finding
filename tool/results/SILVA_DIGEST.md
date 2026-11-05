# SILVA paper: every number we have, with its source (as of 2026-11-05, after the desktop's 6a6bd22)

Paths are relative to `tool/`. The `file:line` is the raw result each number comes from; where a number is derived
(a mean, a ratio), the derivation is given. Protocols:
- **Strict genus accuracy:** `scripts/score_silva.py`, genus level. A read is correct if its true genus alone has
  the largest weight; ties count as wrong. A read with no weight at all counts as wrong. A read whose true genus
  name is missing from SILVA's genus rank, or occurs at more than one node of it, counts as correct, as in Cliffy's
  scoring (`ambiguous_truth=` in each line). Means are over the four regions, unweighted (each region has 200,000
  read pairs).
- **Reads:** every 50th read pair of the aquatic simulation, 200,000 pairs per region, both mates
  (`results/silva/README.md:6–12`); truth from MicrobeMixer.
- **Search:** `katka find --preset silva` (L = 30, counts, explicit RLCSA, `-F 10`) on the default index (`index
  --preset silva`).
- **Timing:** the queue's protocol (`experiments/silva_queue.sh`), on an otherwise idle desktop (i7-8700), each timed run
  in an 11 GB cgroup:
  - each region and mate with one thread pinned to CPU 9, then with 12 threads;
  - time per read is the mean per-read time with one thread (the `timing:` line), over every read of every region and
    mate;
  - reads/s is all reads over the summed wall time of the 12-thread runs.

  Run 2026-11-05 with the tool at 08e8f4b (`results/silva/README.md:89–94`). All 12-thread outputs are byte-identical
  to the one-thread ones.

## 1. Credit schemes on the aquatic reads (`tab:credit`, SILVA column)

From `results/silva_queue/credit/<scheme>.<region>.csv`, line 1 (the genus line), field `strict=`; 200,000 reads
each, `ambiguous_truth=0`. Run: `silva_queue.sh timing credit`, 2026-11-04, at 3647b68 (`results/silva/README.md:27–45`).

| credit scheme | V1_V2 | V3_V4 | V4_V4 | V4_V5 | mean | paper |
|---|---|---|---|---|---|---|
| range (Cliffy's LCA query): `lca.*.csv:1` | 89.04% | 87.57% | 84.66% | 86.51% | 86.945% | 87.0% |
| one label per MEM: `one.*.csv:1` | 93.95% | 92.33% | 91.20% | 91.20% | 92.170% | 92.18% |
| sample (simulated cliffs): `cliffs.*.csv:1` | 92.25% | 91.12% | 88.41% | 90.07% | 90.4625% | 90.47% |
| even: `list.*.csv:1` | 91.90% | 91.32% | 87.79% | 90.10% | 90.2775% | 90.28% |
| proportional, normalised by taxon size: `norm.*.csv:1` | 94.36% | 93.36% | 91.15% | 92.17% | 92.760% | 92.76% |
| proportional: `freq.*.csv:1` | 95.29% | 94.20% | 92.67% | 92.92% | 93.770% | 93.77% |

- **Agreement with the paper:** three of the six means match the paper exactly (even, normalised, proportional).
  - Range: the mean of the rounded per-region values, 86.945%, sits exactly at the rounding boundary of the paper's
    one-decimal 87.0%; the exact counts would decide.
  - "One label" (92.17% against 92.18%) and "sample" (90.46% against 90.47%) differ by 0.01 point. Both are random
    draws (seeded per block of reads), and the per-region values are rounded to 4 decimals before averaging, so the
    paper's values may come from another seed or from exact counts.
  - The scorer now writes `strict_correct=` for exact means, but these CSVs predate it.
- **Rounding:** `results/silva/README.md:73` gives range's aquatic mean as 86.94% and line 35 as 86.95%; the mean of
  the four rounded values is 86.945%.

## 2. Credit schemes on uniform-abundance reads

The `uniform` step: MicrobeMixer with every genus at share 0, 200,000 pairs per region, the same primers; genera that
ART cannot give enough reads were dropped (`results/silva/README.md:47–87`). From
`results/silva_queue/uniform/credit/<scheme>.<region>.csv:1`, `strict=`.

| region | genera with reads | pairs per genus | dropped | `ambiguous_truth` (counted correct) | source |
|---|---|---|---|---|---|
| V1_V2 | 3,283 of 9,088 | 60–61 | 30 | 16,198 (8.1%) | `results/silva_queue/uniform/summary.txt:1–2`; `V1_V2.counts.txt`, `V1_V2.dropped.txt` (30 lines) |
| V3_V4 | 4,135 of 9,085 | 48–49 | 33 | 15,529 (7.8%) | `summary.txt:3–4`; `V3_V4.dropped.txt` (33 lines) |
| V4_V4 | 5,418 of 9,043 | 36–37 | 75 | 12,672 (6.3%) | `summary.txt:5–6`; `V4_V4.dropped.txt` (75 lines) |
| V4_V5 | 8,758 of 8,858 | 22–23 | 260 | 8,921 (4.5%) | `summary.txt:7–8`; `V4_V5.dropped.txt` (260 lines) |

| credit scheme | V1_V2 | V3_V4 | V4_V4 | V4_V5 | mean, uniform | mean, aquatic (§1) | change |
|---|---|---|---|---|---|---|---|
| normalised by taxon size | 95.68% | 92.35% | 90.24% | 88.74% | 91.7525% | 92.760% | −1.01 |
| proportional | 93.00% | 88.13% | 86.44% | 84.79% | 88.090% | 93.770% | −5.68 |
| sample (simulated cliffs) | 91.44% | 87.45% | 85.62% | 85.17% | 87.420% | 90.4625% | −3.04 |
| even | 91.34% | 87.62% | 85.33% | 85.15% | 87.360% | 90.2775% | −2.92 |
| one label per MEM | 91.70% | 86.24% | 84.13% | 81.21% | 85.820% | 92.170% | −6.35 |
| range (Cliffy's LCA query) | 88.56% | 84.49% | 82.63% | 81.94% | 84.405% | 86.945% | −2.54 |

- **Caution for the paper:** under uniform abundances, 4.5–8.1% of the reads come from genera whose name is
  ambiguous in SILVA, and these count as correct for every scheme. The aquatic reads have none. Uniform accuracies
  are therefore inflated relative to the aquatic ones by up to that share. Comparisons between schemes within a
  column are unaffected.
- **Proportional credit's advantage over even credit** is 3.49 points under aquatic abundances and 0.73 under uniform
  (from the means above).
- **Robustness:** normalised credit loses 1.01 points between the two, proportional 5.68 (from the means above).

## 3. The default SILVA index (desktop build)

`results/silva_profile/index_build.log`: `index --preset silva` on `bac`.
- RLBWT: n = 1,173,771,557, r = 68,324,430, 155.265 MB (line 3).
- Explicit RLCSA: 565.231 MB (line 4).
- Tag array: 359,351,512 tag runs, 1,150.312 MB (line 5).
- Build: 348.1 s wall, peak 6,535,916 KB = 6.69 GB (lines 6–7).

These are build-time figures. The space column of `tab:others` (the files a configuration loads) comes from the
queue's `sizes.txt` (pending).

## 4. Cliffy: why our sizes of its indexes are about twice the paper's (9.01 and 4.27 GB)

`NOTES.md:813–823` (the Cliffy section, 2026-11-04; Cliffy at 3763529).
- **What the paper counts:** 10 files, in GiB (`exp_1.smk`, rule `combine_all_index_sizes_into_csv_file_exp1`);
  we counted all 22 files the build writes, in GB (`katka/silva/compare_lib.sh`).
- **On the small test index:** the extra files make the total 1.82 times the paper's count with both strands, and
  1.94 times without (line 816).
- **At SILVA's scale:** `.lcp` took 19 bytes per run, and `.esa` and `.ssa` 10 each (line 817). With r ≈ 70 M and
  both strands, these and `.fna` add roughly 5 GB to the paper's 9.67 GB (9.01 GiB). The other extras and GiB
  against GB make up the rest of our 19.58 GB (line 817).
- **Both strands are not the difference:** Cliffy's experiments use `--revcomp` too (line 818).
- **Queries:** a query opens the paper's 10 files plus `.fna`, `.bwt.heads` and `.bwt.len`. On the small index these
  take 67 MB, against a 242 MB peak (line 821).
- `scripts/cliffy_index_size.sh` prints all three counts for any Cliffy index (line 823). It has not yet been run on
  the server's SILVA indexes (pending).

## 5. B_λ on Zymo (Lemma "runs")

`results/blam/zymo_blam.err`. Per-species tags, homopolymer-compressed text, 2026-11-04/05, 3ab8589
(`results/blam/README.md:6–39`).
- ρ = 448,723,822 tag runs; r = 338,097,187, σ = 6, d = 7 (line 8).

| λ | B_λ | bound (r + σ + 2d)(λ + 1) + B_λ + 1 | ρ / bound | source |
|---|---|---|---|---|
| 0 | 448,723,813 | 786,821,021 | 0.570 | line 9 |
| 10 | 448,430,416 | 4,167,499,694 | 0.108 | line 10 |
| 30 | 107,938,468 | 10,588,951,886 | 0.042 | line 11 |
| 100 | 97,171,998 | 34,244,989,906 | 0.013 | line 12 |
| 1000 | 65,409,566 | 338,500,713,774 | 0.0013 | line 13 |

- **Tag array:** n = 39,277,433,153 rows, 687.420 MB (line 14).
- **Cost:** 45,570.9 s wall, peak 10,771,576 KB = 11.03 GB (lines 15–16).
- **Check:** the tag runs recomputed from the BWT are byte-identical to the converted KATKA file
  (`results/blam/README.md:37–39`).

## 6. Time, space and accuracy (`tab:others`, our rows)

`results/silva/tables.txt` (`silva_tables.py`'s output from `results/silva_queue/`), lines 17–20. Space is the total
size of the files each configuration loads (`results/silva_queue/sizes.txt`, lines 1–4, in bytes); time per read and
reads/s are from the `timing:` lines of `results/silva_queue/timing/*.log` and `digests/csa.L30.*.log`. Accuracy: §1's
scoring of the same outputs (the grammar and digest outputs are scored separately, `credit/freq.grammar.*`,
`freq.bigrepair.*`, `digests/freq.L30.*`).

| classifier | space (GB) | µs per read, 1 thread | total / reads (µs) | reads/s, 12 threads | strict | ties correct | source |
|---|---|---|---|---|---|---|---|
| long MEMs, proportional credit | 1.87 | 52.3 | 53.9 | 145,500 | 93.77% | 93.79% | `tables.txt:17`, `sizes.txt:1` |
| grammar tags, built-in RePair | 1.17 | 53.2 | 54.6 | 142,800 | 93.77% | 93.79% | `tables.txt:18`, `sizes.txt:2` |
| grammar tags, bigrepair | 1.23 | 56.0 | 57.1 | 137,200 | 93.77% | 93.79% | `tables.txt:19`, `sizes.txt:3` |
| minimizer digests, L = 30 | 0.59 | 31.5 | 33.1 | 217,600 | 92.48% | 92.51% | `tables.txt:20`, `sizes.txt:4` |

- **Space** (`results/silva/README.md`, "Index size"): 1.87 GB is `.rlbwt` 0.155 + `.csa` 0.565 + `.tag` 0.836 + `.tag.rmq`
  0.314. Without the RMQ, which counting doesn't use, it is 1.56 GB. The paper's 1.44 GB (`silva.tex:249`) is the RLCSA 0.57 +
  tags 0.86 + the lookup table 0.018, so it never counted the RMQ, and it leaves out the RLBWT.
- **Timing, superseding the row above** (`results/silva/timing_repeats.txt`, 5 repeats over all the reads): one thread
  45.57 µs per read (45.46–45.87), 46.88 µs total / reads, 21,332 reads/s; 12 threads 167,224 reads/s (164,575–167,504).
  The earlier 52.3 µs was a slow run, 45.9 µs agrees with the repeats, and 43.0 µs was V1_V2 mate 1 alone: read files
  run from 42.48 to 51.94 µs.
- **Newer results, in `results/silva/README.md`:** smaller ε and the ε → 0 limit, accuracy excluding ambiguous genera,
  `bayesauto`'s peak, and naive Bayes on the aquatic and uniform reads.

## 7. The digested index (`tab:digest`)

`results/silva/tables.txt`, lines 26–31 (strict accuracy; mean µs per read with one thread; space from `sizes.txt:4–6` and
`:1`). The Elias–Fano and sampled rows' outputs were byte-identical to the explicit RLCSA's in all 64 runs
(`tables.txt:23`).

| | space (GB) | L = 15 | L = 20 | L = 30 | L = 40 | L = 50 | source |
|---|---|---|---|---|---|---|---|
| digests: strict accuracy | | 92.47% | 92.47% | 92.48% | 92.36% | 92.05% | `tables.txt:26` |
| explicit RLCSA, tags | 0.59 | 54.8 | 46.8 | 31.5 | 22.7 | 18.7 | `tables.txt:27` |
| Elias–Fano RLCSA, tags | 0.38 | 88.5 | 82.0 | 65.3 | 53.3 | 46.0 | `tables.txt:28` |
| RLBWT, tags sampled with s = 4 | 0.20 | | | 232.7 | 141.2 | 103.8 | `tables.txt:29` |
| undigested: strict accuracy | | 93.74% | 93.77% | 93.77% | 93.48% | | `tables.txt:30` |
| default configuration | 1.87 | 65.6 | 56.9 | 45.9 | 37.7 | | `tables.txt:31` |

**Build of the digested index** (`results/silva/README.md:126–127`): `prep` 37 s; Big-BWT 128 s (2.11 GB); `index --csa
--efcsa --tags` 81 s (1.43 GB); sampled tags 84 s (2.26 GB).

## 8. Grammar-compressed tags

`results/silva_queue/grammar/summary.txt`: 359,351,512 tag runs over 1,173,771,557 rows (line 2).

| tag array | file (MB) | build | peak | source |
|---|---|---|---|---|
| plain (+ RMQ) | 836.4 (+ 313.9) | | | line 6 |
| built-in RePair | 449.7 | 1,139.4 s | 9.33 GB | lines 3, 6 |
| bigrepair (plain RePair, `-m 0`) | 509.5 | 205.9 s + 69.5 s = 275.4 s | 6.74 GB (bigrepair itself 4.08 GB) | lines 4–6 |

- **Equivalence:** `find -c` and `find -l` outputs are byte-identical with the three tag arrays (lines 7, 11).
- **bigrepair's grammar:** `.R` 456.7 MB, `.C` 9.2 MB (line 4).

## 9. B_λ on SILVA (Lemma "runs")

`results/silva_queue/blam/blam.log`: per-genus tags of the default index.
- ρ = 359,351,512; r = 68,324,430; σ = 6; d = 9,118 (line 3).

| λ | B_λ | bound | ρ / bound | source |
|---|---|---|---|---|
| 0 | 359,351,506 | 427,694,179 | 0.840 | line 4 |
| 10 | 352,215,328 | 1,103,984,721 | 0.326 | line 5 |
| 30 | 243,877,805 | 2,362,500,638 | 0.152 | line 6 |
| 100 | 107,860,112 | 7,010,469,985 | 0.051 | line 7 |
| 1000 | 372,456 | 68,411,387,129 | 0.0053 | line 8 |

- **Cost:** 404.1 s, peak 5.35 GB (line 11).
- **The bound holds at every λ** (ρ / bound computed from the lines above).

## 10. The aquatic reads' composition

`results/silva/alpha_fit.txt` (`scripts/silva_alpha_fit.py` on the aquatic truth files, 10 million pairs per region; the
every-50th set is a uniform sample of them).
- **Read counts don't follow N_g:** 98 to 100 genera have reads (lines 1, 5, 9, 13). The fitted α of k_g ∝ N_g^α is
  0.056, with Spearman −0.037 to −0.045 (lines 2, 6, 10, 14).
- **The genera themselves are well represented:**
  - they hold 59,103 to 59,321 of SILVA's 402,378 sequences (14.7%; lines 3, 7, 11, 15);
  - their median rank by N_g among all 9,118 genera is 620.5, the range 1 to 3,881;
  - 12 of them are in the top 1%, and 63 to 64 in the top 10% (lines 4, 8, 12, 16).

## 11. Credit schemes against the sample's composition (the α dial)

`results/silva/dial_tables.txt` (`dial_tables.py` on the desktop's queue directory; CSVs and logs in
`results/silva_queue/dial/`). 200,000 pairs per region simulated with k_g ∝ N_g^α by `simulate_dial.py`; strict genus
accuracy, mean of the four regions.

| scheme | α = 1 | α = 0 | α = −1 | source |
|---|---|---|---|---|
| range (Cliffy's LCA query) | 85.58% | 84.50% | 83.90% | line 10 |
| one label per MEM | 91.63% | 85.94% | 81.87% | line 11 |
| sample (simulated cliffs) | 88.04% | 87.51% | 87.05% | line 12 |
| even | 88.20% | 87.48% | 87.01% | line 13 |
| normalised | 91.76% | 91.85% | 91.83% | line 14 |
| proportional | 92.68% | 88.19% | 85.01% | line 15 |
| prior mixture, β = 0.5 | 92.67% | 90.66% | 90.42% | line 19 |
| EM, 20 iterations | 92.31% | 90.20% | 90.22% | line 29 |
| naive Bayes, the sample's prior (reference), ε = 10⁻⁶ | 92.77% | 94.23% | 95.89% | line 30 |
| naive Bayes, flat prior, ε = 10⁻² | 41.58% | 86.28% | 95.99% | line 31 |
| naive Bayes, flat prior, ε = 10⁻³ | 65.28% | 94.54% | 96.57% | line 32 |
| naive Bayes, flat prior, ε = 10⁻⁴ (the held-out choice) | 84.49% | 95.61% | 96.18% | lines 33, 39–44 |
| naive Bayes, flat prior, ε = 10⁻⁶ | 92.66% | 94.23% | 94.78% | line 34 |
| naive Bayes, flat prior, ε = 10⁻⁹ | 91.11% | 92.44% | 93.10% | line 35 |
| naive Bayes, database prior, ε = 10⁻⁶ | 92.77% | 93.20% | 93.00% | line 36 |
| naive Bayes, mixture, ε = 0.5 (not a fit: see below) | 91.42% | 94.22% | 95.11% | line 37 |

- **Check against the uniform run:** α = 0 repeats the uniform design with the other simulator. Normalised gets 91.85
  against 91.75, proportional 88.19 against 88.09, range 84.50 against 84.41, and one label 85.94 against 85.82
  (`results/silva/README.md:169–170`).
- **ε for the flat prior:** the held-out choice (the best at α = 0) is 10⁻⁴, which gives 84.49% at α = 1
  (`dial_tables.txt:39–44`). ε = 10⁻⁶ is better at α = 1 but was not chosen held out.
- **The mixture's ε is not a fit.** It is 0.5, the top of the grid, in all 12 runs (`dial_tables.txt:49–51`), with the
  log-likelihood still rising there (`dial/a*/credit/bayesauto.*.csv.log`).
  - The logs explain why. Near ε = 0 their slope is 1.0–1.3 million per decade over 200,000 pairs: the most likely
    genus lacks 4.6 to 6.5 of a read pair's MEMs, at every α and region. On the toy it was 0.2 to 0.8.
  - The probable cause is MicrobeMixer's 2% mutation of each region and ART's errors. Around a changed base, a read
    segment can match some other genus exactly over 30 or more bases, while the true genus lacks it.
  - So the maximum-likelihood ε is probably above 0.5. The grid now goes to 0.99 (fb7c06c); the row needs rescoring
    (see Pending).

## 12. Pending

| item | what it gives | where it will be | status |
|---|---|---|---|
| `bayesauto` rescored with the grid to 0.99 | the mixture's real fitted ε on SILVA, and its accuracy | `silva_queue.sh dial` after deleting `dial/a*/credit/bayesauto.*` (scoring only) | not yet run |
| Cliffy's approximate listing on our reads | its accuracy under the paper's scoring | `scripts/score_cliffy_listing.sh` | needs Gonzalo's run of Cliffy |
| which reads Cliffy's `tab:others` numbers are on | 5,000 or 200,000 pairs per region | `results/silva/README.md:14–25` | to be settled with Gonzalo |
| Cliffy index sizes on the server | the three counts for its SILVA indexes | `scripts/cliffy_index_size.sh` | not yet run |
