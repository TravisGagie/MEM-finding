# SILVA paper: every number we have, with its source (as of 2026-11-05, commit after edb28a0)

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
- **Timing:** the queue's protocol (one thread pinned, then 12 threads; `experiments/silva_queue.sh`). **No timing
  under this protocol has been run on an idle machine yet**; see "Pending".

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

## 6. Pending (not yet run, or run but not yet in the repository)

| item | what it gives | where it will be | status |
|---|---|---|---|
| timing protocol on an idle machine | `tab:others` time per read (1 thread), reads/s (12 threads), space; `tab:digest` times | `silva_queue.sh timing` → `results/silva_queue/timing/`, `sizes.txt` | the 2026-11-04 timings ran beside B_λ and were not committed; to be rerun idle |
| B_λ on SILVA | B_λ and the bound for SILVA's per-genus tags | `silva_queue.sh blam` → `blam/blam.log` | not yet run |
| grammar tags | built-in RePair against bigrepair: sizes, times, peaks, identical answers | `silva_queue.sh grammar` → `grammar/summary.txt` | not yet run |
| digested index | `tab:digest` and the digests row of `tab:others` | `silva_queue.sh digests` → `digests/` | not yet run |
| credit schemes against the sample's composition | α fit of the aquatic reads; every scheme at α = 1, 0, −1, including naive Bayes (flat, database and sample priors, ε sweep, `bayesauto`) | `silva_alpha_fit.py`; `silva_queue.sh dial` → `dial/` | not yet run (toy results in `NOTES.md`, the dial section) |
| Cliffy's approximate listing on our reads | its accuracy under the paper's scoring | `scripts/score_cliffy_listing.sh` | needs Gonzalo's run of Cliffy |
| which reads Cliffy's `tab:others` numbers are on | 5,000 or 200,000 pairs per region | `results/silva/README.md:14–25` | to be settled with Gonzalo |
| Cliffy index sizes on the server | the three counts for its SILVA indexes | `scripts/cliffy_index_size.sh` | not yet run |

Once the queue's results are collected, `python3 scripts/silva_tables.py results/silva_queue` and `python3
scripts/dial_tables.py results/silva_queue` print every table above that is still pending, and name any input that
is missing.
