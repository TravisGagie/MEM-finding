# SILVA (Section 4, for the Cliffy paper): desktop runs with the new tool

## Reads

All accuracy runs here score the same reads:
- **The set:** every 50th read pair of MicrobeMixer's simulation (10 million pairs of 250-base reads per region,
  aquatic abundances), that is pairs 0, 50, 100, … (reads 1, 51, 101, … of `aquatic/<region>/<region>_mate_{1,2}.fq`).
- **Size:** 200,000 pairs per region, both mates, for V1_V2, V3_V4, V4_V4 and V4_V5.
- **Files:** `katka/silva/work/compare/reads_every50/<region>_{1,2}.fq`. On 2026-11-04, `cmp` showed all 8 files
  byte-identical to the release sample `work/reads/sample50/<region>_mate_{1,2}.fq.gz` (`make_samples.sh`:
  `awk 'int((NR-1)/4) % 50 == 0'`).
- **Truth:** `work/reads/aquatic/<region>/<region>_seqtax.txt`.

**Discrepancy, for Travis to raise with Gonzalo.** The paper describes two read sets:
1. **Accuracy (§4, the SILVA column of Table `tab:credit`):** every 50th pair, 200,000 pairs per region
   (`app-katka-details.tex` lines 35 and 123). The runs here use this set.
2. **The Cliffy comparison (Table `tab:others`):** every 2,000th pair, 5,000 per region (`work/reads/sample2000/`).
   `app-katka-details.tex` line 95 says Gonzalo ran Cliffy "on our 5,000 read pairs per region", and the
   `silva.tex` line 262 caption says "Cliffy's numbers are on our reads (5,000 pairs per region)". Our own desktop
   run of the old comparison (`work/compare/katka-travis-gagie-desktop/summary.txt`) also says "every 2000th pair".

According to Travis, the header of Gonzalo's Cliffy summary says "every 50th pair". That summary is not on the
desktop, so this was not checked here. If so, Gonzalo's Cliffy numbers are on the 200,000-pair set, and the
paper's description of the Cliffy runs (5,000 pairs) is wrong, or refers to a different run. The runs here use
the every-50th set, so that they are on the same reads as Gonzalo's summary.

## Accuracy of the credit schemes (item 2; 2026-11-04, 3647b68)

`experiments/silva_queue.sh timing credit` (`find --preset silva`: L = 30, counts, explicit RLCSA, `-F 10`, on the
reads above, the tool's default SILVA index), scored per region under each scheme with `scripts/score_silva.py`.
Per-region CSVs are in `../silva_queue/credit/`. The numbers are strict genus accuracy:

| credit scheme | mean | V1_V2 | V3_V4 | V4_V4 | V4_V5 |
|---|---|---|---|---|---|
| range (Cliffy's LCA query) | 86.95% | 89.04% | 87.57% | 84.66% | 86.51% |
| one label per MEM (drawn by counts) | 92.17% | 93.95% | 92.33% | 91.20% | 91.20% |
| sample (simulated cliffs) | 90.46% | 92.25% | 91.12% | 88.41% | 90.07% |
| even | 90.28% | 91.90% | 91.32% | 87.79% | 90.10% |
| proportional, normalised by taxon size | 92.76% | 94.36% | 93.36% | 91.15% | 92.17% |
| proportional | **93.77%** | 95.29% | 94.20% | 92.67% | 92.92% |

Proportional and even credit reproduce the paper's 93.77% and 90.28% at L = 30.

These ran while B_lambda on Zymo was in its φ pass. Accuracy doesn't depend on that, and the times don't count: the
timing step is rerun on an idle machine. Its outputs are deterministic, so the scores stand.

## Is proportional credit's advantage the prior? A uniform-abundance simulation (item 6; 2026-11-04, 7378117)

`experiments/uniform_microbemixer.sh` (through `silva_queue.sh uniform`): MicrobeMixer with every genus of the
reference listed at a share of 0, so reads are spread evenly over the genera it can simulate, 200,000 read pairs per
region, same primers and command line as the aquatic simulation. `find --preset silva` on both mates, then every
scheme scored as in the queue's `credit` step (`uniform_schemes.sh`). Outputs are in `../silva_queue/uniform/`.

| region | genera with reads | reads per genus | dropped for too few reads |
|---|---|---|---|
| V1_V2 | 3,283 | 60–61 | 30 |
| V3_V4 | 4,135 | 48–49 | 33 |
| V4_V4 | 5,418 | 36–37 | 75 |
| V4_V5 | 8,758 | 22–23 | 260 |

MicrobeMixer needs more reads from a genus than ART's fold coverage gives when its sequences are short, so those
genera are removed from the table and the region is rerun (`dropped.txt`).

Strict genus accuracy:

| credit scheme | V1_V2 | V3_V4 | V4_V4 | V4_V5 | **mean, uniform** | mean, aquatic | change |
|---|---|---|---|---|---|---|---|
| normalised by taxon size | 95.68% | 92.35% | 90.24% | 88.74% | **91.75%** | 92.76% | −1.01 |
| proportional | 93.00% | 88.13% | 86.44% | 84.79% | **88.09%** | 93.77% | −5.68 |
| sample (simulated cliffs) | 91.44% | 87.45% | 85.62% | 85.17% | **87.42%** | 90.46% | −3.04 |
| even | 91.34% | 87.62% | 85.33% | 85.15% | **87.36%** | 90.28% | −2.92 |
| one label per MEM | 91.70% | 86.24% | 84.13% | 81.21% | **85.82%** | 92.17% | −6.35 |
| range (Cliffy's LCA query) | 88.56% | 84.49% | 82.63% | 81.94% | **84.41%** | 86.94% | −2.53 |

**Reading:**
- **Most of proportional credit's advantage comes from the prior.** Under aquatic abundances it beats even credit by
  3.49 points; under uniform abundances, by 0.73, and on V4_V5 it falls just below even (84.79% against 85.15%).
  Counts favour genera with many sequences. That helps when, as in the aquatic community, abundant genera are the
  well-represented ones, and hurts when every genus is equally likely.
- **Normalised credit, which divides each count by the genus's number of sequences, is the most robust.** It is
  1.01 points below proportional under aquatic abundances, but 3.66 above it under uniform, and it loses only 1.01
  points between the two priors, against 5.68.
- **Against the range:** every scheme beats Cliffy's LCA query under both priors. Normalised credit leads it by 7.34
  points under uniform, proportional by 3.68.

The uniform reads have fewer reads per genus and many more rare genera than the aquatic ones, so the two columns are
not the same difficulty. All the schemes lose accuracy; the comparison that matters is between schemes within a column.

## The SILVA queue (items 1, 3, 4, 5; 2026-11-05, tool built at 08e8f4b)

`experiments/silva_queue.sh timing blam grammar digests`, on an otherwise idle machine, each timed run in an 11 GB
cgroup, one thread pinned to CPU 9, then 12 threads. `tables.txt` and `tables.tex` are `silva_tables.py`'s output, from
`../silva_queue/`. All outputs with 12 threads are byte-identical to those with one, and the grammar tags' and the
Elias-Fano and sampled digests' answers are byte-identical to the default's.

**Rescoring item 2.** `score_silva.py`'s "exact strict count" (1edd8ee) adds a column, `strict_correct`, and doesn't
change the definition. All 144 rows of item 2's scores are identical to before, apart from that column. No number changed.

**Time, space and accuracy** (genus, 200,000 pairs per region; time is the mean per read with one thread):

| classifier | space (GB) | µs per read | total / reads (µs) | reads/s, 12 threads | strict | ties correct |
|---|---|---|---|---|---|---|
| long MEMs, proportional credit | 1.87 | 52.3 | 53.9 | 145,500 | 93.77% | 93.79% |
| … grammar tags, built-in RePair | 1.17 | 53.2 | 54.6 | 142,800 | 93.77% | 93.79% |
| … grammar tags, bigrepair | 1.23 | 56.0 | 57.1 | 137,200 | 93.77% | 93.79% |
| minimizer digests, L = 30 | 0.59 | 31.5 | 33.1 | 217,600 | 92.48% | 92.51% |

The default's 1.87 GB includes the tag array's RMQ (0.31 GB), which counting does not use; the paper's 1.44 GB did not
have it as a separate file.

**Run-to-run noise.** The timing step (00:59) gave 52.3 µs for the default, while the same configuration in the digests
step (02:29, `undigested.L30`) gave 45.9 µs. The two runs give identical outputs, so a 14% difference between idle-machine
runs is a measure of the noise. Compare times here only with care; repeats would be needed for a fine comparison.

**Grammar-compressed tags (item 4):**

| tag array | file (MB) | build | peak |
|---|---|---|---|
| plain (+ RMQ) | 836.4 (+ 313.9) | | |
| built-in RePair | 449.7 | 1,139 s | 9.33 GB |
| bigrepair | 509.5 | 275 s | 6.74 GB |

The built-in RePair's tag array is smaller, at about four times the build time; bigrepair's grammar is .R 456.7 MB and
.C 9.2 MB.

**The digested index (item 5):** build: `prep` 37 s, Big-BWT 128 s (2.11 GB), `index --csa --efcsa --tags` 81 s (1.43 GB),
sampled tags 84 s (2.26 GB). Strict accuracy, then the mean µs per read with one thread:

| | space | L = 15 | L = 20 | L = 30 | L = 40 | L = 50 |
|---|---|---|---|---|---|---|
| digests: strict accuracy | | 92.47% | 92.47% | 92.48% | 92.36% | 92.05% |
| explicit RLCSA, tags | 0.59 GB | 54.8 | 46.8 | 31.5 | 22.7 | 18.7 |
| Elias-Fano RLCSA, tags | 0.38 GB | 88.5 | 82.0 | 65.3 | 53.3 | 46.0 |
| RLBWT, tags sampled with s = 4 | 0.20 GB | | | 232.7 | 141.2 | 103.8 |
| undigested: strict accuracy | | 93.74% | 93.77% | 93.77% | 93.48% | |
| default configuration | 1.87 GB | 65.6 | 56.9 | 45.9 | 37.7 | |

**The aquatic reads' composition (`silva_alpha_fit.py`, `alpha_fit.txt`).** Over the about 100 genera with reads, the
fitted α of k_g ∝ N_g^α is 0.056, with Spearman's correlation about −0.04: within those genera, read counts do not follow
SILVA's sequence counts. But the genera themselves are well represented: they hold 14.7% of SILVA's sequences, with a
median rank by N_g of 620 of 9,118, and 12 in the top 1%. This qualifies the reading of the uniform simulation above:
the aquatic reads favour well-represented genera through which genera are present, not through how many reads each gets.

## The α dial (2026-11-05, 07b5ffe)

`MM=~/rz/MicrobeMixer tool/experiments/silva_queue.sh dial`: 200,000 simulated pairs per region with k_g ∝ N_g^α,
for α = 1, 0 and −1, on the tool's default SILVA index. Full tables: `dial_tables.txt`, `dial_tables.tex`, from
`dial_tables.py ~/tool/silva/queue` (run on the queue directory: it also reads `dial/a*/<region>.seqtax.txt`, 1.1 MB
each, which were not committed). The CSVs and logs are in `../silva_queue/dial/`.

Strict genus accuracy, mean of the four regions:

| scheme | α = 1 | α = 0 | α = −1 |
|---|---:|---:|---:|
| range (Cliffy's LCA query) | 85.58 | 84.50 | 83.90 |
| one label per MEM | 91.63 | 85.94 | 81.87 |
| sample (simulated cliffs) | 88.04 | 87.51 | 87.05 |
| even | 88.20 | 87.48 | 87.01 |
| normalised | 91.76 | 91.85 | 91.83 |
| proportional | 92.68 | 88.19 | 85.01 |
| prior mixture, β = 0.5 | 92.67 | 90.66 | 90.42 |
| EM, 20 iterations | 92.31 | 90.20 | 90.22 |
| naive Bayes, the sample's prior (reference), ε = 10⁻⁶ | 92.77 | 94.23 | 95.89 |
| naive Bayes, flat prior, ε = 10⁻⁴ (best at α = 0) | 84.49 | 95.61 | 96.18 |
| naive Bayes, flat prior, ε = 10⁻⁶ | 92.66 | 94.23 | 94.78 |
| naive Bayes, database prior, ε = 10⁻⁶ | 92.77 | 93.20 | 93.00 |
| naive Bayes, mixture, ε at the likelihood's peak (0.5 or 0.6) | 91.42 | 94.26 | 95.19 |

- **Check:** α = 0 is the same design as the uniform simulation above, simulated again by `simulate_dial.py`: normalised 91.85 against
  91.75, proportional 88.19 against 88.09, range 84.50 against 84.41, one 85.94 against 85.82.
- **ε = 10⁻⁶ was not chosen held out.** The held-out choice, the best ε at α = 0, is 10⁻⁴, and it gives 84.49 at
  α = 1. The other values are in `dial_tables.txt`.
- **The mixture's ε:** the first run's grid stopped at 0.5, the peak's edge. On the grid to 0.99 the likelihood peaks
  inside it, at 0.5 or 0.6 (see below).
- **Run note:** `simulate_dial.py` fails under this machine's Python 3.14 (see `NOTES.md`), so the run set the
  multiprocessing start method to `fork` from outside the script (`PYTHONPATH` with a `sitecustomize.py`). The script
  was not changed.

## Rescoring on the existing listings (2026-11-05)

No new simulations or searches, except that the uniform listings were rescored to get exact strict counts, which
reproduced the earlier uniform column. Full tables: `rescore_tables.txt`, from `rescore_tables.py`; scoring CSVs in
`../silva_queue/`. Naive Bayes here has a flat prior.

**Smaller ε and the limit** (`score_silva.py --epsilon 0`: fewest missing MEMs, ties by the sum of log c_g / (2 N_g)).
Strict genus accuracy, mean of the four regions, all reads and excluding reads from genera with ambiguous names:

| ε | α = 1 | excl. | α = 0 | excl. | α = −1 | excl. | aquatic |
|---|---:|---:|---:|---:|---:|---:|---:|
| 10⁻⁶ | 92.66 | 91.94 | 94.23 | 93.84 | 94.78 | 94.65 | 94.70 |
| 10⁻⁹ | 91.11 | 90.23 | 92.44 | 91.92 | 93.10 | 92.90 | 94.36 |
| 10⁻¹² | 90.35 | 89.40 | 91.59 | 91.00 | 92.28 | 92.06 | 94.19 |
| limit | 89.68 | 88.66 | 90.80 | 90.16 | 91.48 | 91.23 | 94.08 |

- **Checks:** ε = 10⁻⁹ repeats the dial's earlier row exactly, and ε = 10⁻³⁰⁰ gives the limit's accuracy to two decimals.
- **Scored clades only:** the scorer now picks among scored clades, because an "uncultured" top could underflow the
  other genera's posteriors to 0 at very small ε. This doesn't change the larger ε values above.
- **Real aquatic reads:** proportional 93.77 and normalised 92.76, so naive Bayes beats both at every ε here, the
  limit included. ε = 10⁻⁶ wasn't chosen held out.

**Excluding reads from genera with ambiguous names** (the scorer counts those reads as correct). Reads per region,
of 200,000: α = 1, 17,609–18,267; α = 0, 9,000–16,215; α = −1, 2,452–9,712; uniform, 8,921–16,198; aquatic, 0.
So the aquatic table is unchanged, and the simulated tables fall by 0.1 to 1.5 points for most schemes:

| scheme | uniform | excl. | aquatic (both) |
|---|---:|---:|---:|
| range | 84.41 | 83.32 | 86.95 |
| one label per MEM | 85.83 | 84.88 | 92.17 |
| sample (simulated cliffs) | 87.42 | 86.55 | 90.46 |
| even | 87.36 | 86.49 | 90.28 |
| normalised | 91.75 | 91.20 | 92.76 |
| proportional | 88.09 | 87.28 | 93.77 |

The whole dial, both ways, is in `rescore_tables.txt`.

**`bayesauto`'s likelihood peak**, on the grid to 0.99: ε = 0.5 in 9 of the 12 runs and 0.6 in 3 (α = 0, V4_V5;
α = −1, V3_V4 and V4_V5). None is at the grid's edge. The grid is coarse there (0.3, 0.5, 0.6, 0.7). Accuracy: 91.42, 94.26, 95.19.

**The aquatic genera** (`alpha_fit.txt`): 98 genera in V1_V2 and 100 in the other regions, of 9,118. They hold 14.7% of
SILVA's sequences (59,103 of 402,378 in V1_V2, 59,321 elsewhere). Their median rank by N_g is 620.5, range 1 to 3881,
with 12 in the top 1% and 63–64 in the top 10%. N_g runs from 5 to 15,700, median 111.5.

**Timing, 5 repeats** (`timing_repeats.txt`): mean per-read time 45.57 µs (median; range 45.46–45.87), 167,224
reads/s with 12 threads (164,575–167,504). This replaces the queue's single 52.3 µs and 145,500 reads/s, which ran with
other load on the machine; the 45.9 µs seen in another run agrees with these repeats.
