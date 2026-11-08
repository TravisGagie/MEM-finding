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

The default's 1.87 GB includes the tag array's RMQ (0.31 GB), which counting does not use. See "Index size" below for the paper's 1.44 GB.

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

**On Zymo the grammar doesn't pay:** bigrepair makes the tag array larger, 687.4 to 750.5 MB, a ratio of 1.09 against 0.61 here. It has 1.33 tag runs per BWT run against SILVA's 5.3, and with 7 tags the sequence is only 3 bits per run, 24% of the array (`../zymo/README.md`, "Grammar-compressed tags").

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

**Timing, 5 repeats** (`timing_repeats.txt`), all 1.6 million reads:

| | median | range |
|---|---:|---:|
| 1 thread, mean per read (µs) | 45.57 | 45.46–45.87 |
| 1 thread, total / reads (µs) | 46.88 | 46.71–47.15 |
| 1 thread, reads/s | 21,332 | 21,210–21,410 |
| 12 threads, reads/s | 167,224 | 164,575–167,504 |

The digest's three figures for "the same configuration" are not three samples of one quantity:
- **43.0 µs** (`grammar/summary.txt:8`) is V1_V2 mate 1 alone, a single file. That file's median over the repeats is 42.48 µs; the files run from 42.48 to 51.94.
- **45.9 µs** (digests step) is over all the reads, and agrees with the repeats.
- **52.3 µs** (timing step) is the outlier: every one of its 8 files was slower than its range in the repeats. Its cause was not established. It should be replaced.

**Uniform and aquatic, naive Bayes** (flat prior):

| scheme | uniform | excl. | aquatic |
|---|---:|---:|---:|
| normalised | 91.75 | 91.20 | 92.76 |
| proportional | 88.09 | 87.28 | 93.77 |
| ε = 10⁻⁶ | 94.33 | 93.94 | 94.70 |
| `bayesauto` | 94.32 | 93.93 | 94.12 |

`bayesauto`'s ε: uniform 0.5, 0.5, 0.5, 0.6; aquatic 0.3, 0.5, 0.5, 0.5 (V1_V2, V3_V4, V4_V4, V4_V5). None is at the grid's edge.

**Index size** (bytes, `~/tool/silva/bac.*`, and `silva.tex` line 249):

| part | 1.87 GB | 1.56 GB | paper's 1.44 GB |
|---|---:|---:|---:|
| RLCSA (`.csa`) | 565.2 MB | 565.2 MB | 0.57 GB |
| tag array (`.tag`) | 836.4 MB | 836.4 MB | 0.86 GB (19.1 bits per run) |
| RLBWT (`.rlbwt`) | 155.3 MB | 155.3 MB | not counted |
| RMQ (`.tag.rmq`) | 313.9 MB | not counted | not counted |
| lookup table | not counted | not counted | 18 MB |

- **From 1.56 to 1.44:** take off the RLBWT (155.3 MB), add the lookup table (18 MB), and the tag array is 20.6 MB larger in the paper's count. The RMQ is not the cause, contrary to the note above and the digest: 1.44 never counted it, and removing it leaves 1.56.
- **The RLBWT:** `find` always loads it (`main.cpp` line 872), but with `--csa` the search steps use the RLCSA (`search.hpp` line 27), and the tags use it only when sampled. So it appears unused in the default configuration. This comes from reading the code, not from a run without it.
- **The lookup table:** the new tool builds it at load time (`main.cpp` line 962), so it's in memory but in no file. `sizes.txt` misses it. Peak memory is 1.58 GB.
- **The tag array:** both have 359,351,512 runs and 14-bit tags (628.9 MB). The paper's 0.86 is `PLAN.md`'s 0.228 + 0.629 GB. The old file on the desktop, `bac.s1.tag` (2026-10-05), is 970,965,362 bytes, which with the RLCSA gives the 1.54 GB in `compare/katka-travis-gagie-desktop/summary.txt`. That file doesn't match 0.86 GB, and this wasn't explained.

## MicrobeMixer and ART: what Cliffy's paper says and what its code runs (2026-11-05)

**MicrobeMixer's version.** No source pins one.
- **Papers:** the Genome Research paper (Ahmed, Boucher and Langmead, 2025, section "Simulating realistic 16S rRNA data sets") and the bioRxiv preprint (10.1101/2024.05.25.595899, posted 2024-05-30) name the tool and link the repository, with no version or commit.
- **`cliffy-experiments`:** `config/config.yaml` runs a private path, `/home/oahmed6/Sim_16S_Reads/src/sim_16s_reads.py`, with no commit. A public repository of that name doesn't exist.
- **Public history** (`oma219/MicrobeMixer`, 10 commits, no tags or releases): the behaviour changed once after the first days, in the last commit, 6409f5d (2024-04-25). It adds `--num-reads` and `--output-name`, and sets ART's fold coverage per genus, instead of `-f 2000`.
- **What Cliffy must have run:** `cliffy-experiments` (3b7df02, 2024-05-01) passes `--num-reads 10000000` and `--output-name`, which only 6409f5d accepts, unless the private copy differed. Our reads use 6409f5d too, so on the public record the commit is the same.
- **One loose end:** the preprint's Appendix A.4 shows `simulate` without `--num-reads`, which 6409f5d requires, but `-f {coverage}`, which only 6409f5d has. So the appendix is abbreviated, or describes neither commit exactly.

**ART's error profile. These differ, and this is the mismatch.**
- **Paper:** Figure 2's caption, "simulated with art_illumina using MiSeq v3 error model", in both the preprint and the journal.
- **Code:** every commit of MicrobeMixer runs `art_illumina -ss MSv1`, and the preprint's own Appendix A.4 prints that command. `MSv1` is ART's MiSeq v1 profile; MiSeq v3 is `MSv3`. Both are 250 bp in ART 2.5.8.
- **Not available here:** the journal's Supplemental Section 4, which holds the journal version's commands.

So the reads were most likely made with MSv1, as ours were, and the caption is wrong; this can't be settled from the paper alone.

**The two profiles compared** (`experiments/aquatic_art_profiles.sh`; `art_profile_tables.txt`; CSVs in `art_profiles/`).
Aquatic abundances, the same four regions, 10 million pairs per region, every 50th pair kept, `find --preset silva`, and every scheme.
The repeat MSv1.1 is the reads used everywhere above; the other five are new, with 3 repeats per profile. The sampled reads are on the
external drive (`moved_2026-11-05/art_profiles_samples`), since neither ART nor MicrobeMixer is seeded. Genus accuracy, mean of the four
regions, mean of the repeats (range):

| scheme | MSv1, strict | MSv3, strict | MSv1, Cliffy's scoring | MSv3, Cliffy's scoring |
|---|---:|---:|---:|---:|
| range | 86.68 (86.45–86.95) | 86.01 (85.85–86.20) | 87.33 (87.10–87.60) | 86.86 (86.73–87.00) |
| one label per MEM | 92.01 (91.81–92.17) | 91.82 (91.75–91.91) | 92.09 | 91.93 |
| sample (simulated cliffs) | 90.34 (90.22–90.46) | 89.61 (89.51–89.67) | 90.91 | 90.49 |
| even | 90.20 (90.06–90.28) | 89.38 (89.25–89.46) | 91.75 | 91.40 |
| normalised | 92.66 (92.52–92.76) | 92.06 (91.94–92.16) | 92.66 | 92.07 |
| proportional | 93.59 (93.44–93.77) | 93.28 (93.23–93.37) | 93.62 | 93.34 |
| EM, 20 iterations | 96.32 (96.13–96.43) | 96.03 (95.94–96.14) | 96.32 | 96.03 |
| naive Bayes, flat, ε = 10⁻⁶ | 94.57 (94.41–94.70) | 93.98 (93.91–94.02) | 94.57 | 93.98 |
| naive Bayes, mixture | 93.95 (93.74–94.12) | 93.22 (93.13–93.29) | 93.95 | 93.22 |

Range by region, with Cliffy's scoring, mean of the repeats:

| region | MSv1 | MSv3 |
|---|---:|---:|
| V1_V2 | 89.34 | 89.03 |
| V3_V4 | 87.70 | 87.42 |
| V4_V4 | 85.55 | 85.89 |
| V4_V5 | 86.74 | 85.11 |

- **The profile doesn't explain the gap to Cliffy's published LCA accuracy.** `silva.tex` line 256 gives 87.5 to 90.2% for Cliffy's own reads. Those values aren't in the paper's text as extracted here, so they probably come from its figures or supplement, and weren't re-traced. Range with Cliffy's scoring is 85.55 to 89.34 under MSv1 and 85.11 to 89.03 under MSv3. So MSv3 is slightly further away on average, not closer.
- **The profile's effect is small:** range is lower under MSv3 by 0.67 points strict and 0.47 with Cliffy's scoring. Every scheme is lower under MSv3, by 0.2 to 0.8 points, and the strict means' ranges don't overlap for most schemes. By region it isn't uniform: V4_V5 falls by 1.6, while V4_V4 rises by 0.3.
- **Repeats vary by about ±0.15 to 0.25 points.** Repeat MSv1.1 is the best of the three under nearly every scheme, so the numbers used so far (93.77, 94.70, 86.95) are at the top of their ranges; 93.59, 94.57 and 86.68 are the better estimates.
- **Range here isn't Cliffy itself:** it is our MEMs at L = 30 with an LCA rule, not Cliffy's own query.
- **EM is the best scheme on the aquatic reads,** 96.32 against naive Bayes's 94.57, by about 1.8 points on every repeat. It wasn't in the earlier aquatic table. On the dial it was behind naive Bayes (92.31, 90.20, 90.22). The aquatic reads come from about 100 genera, and EM learns that from the reads themselves, so this doesn't carry over to a sample with many genera.

## Cliffy's published accuracies: Figure 3 (2026-11-05)

**No result files are public.**
- `cliffy-experiments` (566f647, all 10 commits' file lists checked) holds 16 files of code and no results. `analysis/exp3.r` line 45 reads `accuracy.csv` from `/Users/omarahmed/Downloads/work_dir/cliffy_paper/exp_3/`, which was never committed; line 279 of `src/analyze_classifications.py` writes per-sample CSVs, also not committed.
- `oma219/cliffy` has `analysis/analyze_cliffy.py` and `analyze_doclisting.py`, and example data only.
- The journal's Supplemental Section 4 and Code weren't available here.

**So the values were read from the figure.** Figure 3 is vector graphics, so its points are exact coordinates: page 9 of the Genome Research PDF and page 13 of the bioRxiv PDF. `cliffy_fig3_values.py` reads them, and `cliffy_fig3_values.txt` has all 24 series at all 6 ranks. The two versions give the same values to two decimals. The axis calibration uses the 6 tick marks, labelled 0.75 to 1.00 in the preprint; ticks are 19 units per 5 points of accuracy, so reading precision is about 0.01 points. Series are identified as in `exp3.r`: `y = sensitivity` (line 69), the legend's "Cliffy (TP/n)"; line type gives the query; marker shape gives digestion (from line 88).

Genus-level accuracy on the aquatic reads:

| query | digestion | V1–V2 | V3–V4 | V4 | V4–V5 | range | mean |
|---|---|---:|---:|---:|---:|---:|---:|
| LCA | none | 87.42 | 90.16 | 88.31 | 89.76 | 87.4–90.2 | 88.91 |
| LCA | minimizers | 86.27 | 88.31 | 87.53 | 88.49 | 86.3–88.5 | 87.65 |
| LCA | DNA minimizers | 86.86 | 88.12 | 87.84 | 88.58 | 86.9–88.6 | 87.85 |
| approximate listing | none | 94.06 | 94.72 | 93.18 | 94.17 | 93.2–94.7 | 94.03 |
| approximate listing | minimizers | 93.38 | 93.35 | 92.88 | 93.30 | 92.9–93.4 | 93.23 |
| approximate listing | DNA minimizers | 93.95 | 93.15 | 92.71 | 93.22 | 92.7–94.0 | 93.26 |

- **The draft's ranges are close, with small errors.** LCA undigested: draft 87.5–90.2, figure 87.4–90.2. LCA minimizers: draft 86.3–88.5, which matches. Listing: draft 93.3–94.8, figure 93.2–94.7 undigested.
- **The script and the figure disagree on read count.** `exp3.r` line 63 onward titles the panels "n=200k", while the published panels say "n=10mil", as does the caption. The script was edited after, or the titles were; this doesn't settle how many reads were scored.
- **The gap to our range depends on the region.** With Cliffy's scoring, MSv1, mean of 3 repeats:

| region | Cliffy LCA, undigested, its reads | range, our reads | ours − Cliffy |
|---|---:|---:|---:|
| V1–V2 | 87.42 | 89.34 | +1.92 |
| V3–V4 | 90.16 | 87.70 | −2.46 |
| V4 | 88.31 | 85.55 | −2.76 |
| V4–V5 | 89.76 | 86.74 | −3.02 |
| mean | 88.91 | 87.33 | −1.58 |

  The difference changes sign by region, so a single offset, such as "our reads are harder", doesn't describe it. Cliffy's query and our range rule also differ, so this is not a test of the reads alone.

## Cliffy's published times: Figure 2 (2026-11-05)

Figure 2 plots the time to classify 10 million reads with approximate document listing, per biome and region, beside Kraken 2. It is vector graphics too: page 9 of the Genome Research PDF and page 12 of the preprint. `cliffy_fig2_values.py` reads it, and `cliffy_fig2_values.txt` has all three biomes. The y axis has 9 ticks, 500 to 4500 s, labelled in the preprint; the two versions agree to within 2 s. The series come from the figure's own legend, because `exp3.r`'s line types don't match it: `exp3.r` lines 11–35 give no digestion "longdash", while the figure has dashed for DNA minimizers, dotted for minimizers and dot-dash for no digestion. `time.csv`, which it plots, isn't public.

Aquatic reads, seconds as plotted:

| method | V1–V2 | V3–V4 | V4 | V4–V5 | mean |
|---|---:|---:|---:|---:|---:|
| Cliffy, no digestion | 3975 | 3996 | 4078 | 4068 | 4029 |
| Cliffy, DNA minimizers | 3493 | 3726 | 3850 | 3835 | 3726 |
| Cliffy, minimizers | 1524 | 1640 | 1595 | 1620 | 1595 |
| Kraken 2 | 372 | 370 | 386 | 386 | 378 |

These agree with the paper's text (minimizers "on average 4.5× slower" than Kraken 2, the other modes "closer to 10×–11× slower").

**Time per read depends on what was timed, which the public files don't settle.** Cliffy runs each mate separately under GNU time, and each has two trials (`exp_3.smk` lines 46–77); `time.csv` may hold the sum of the mates or one of them. Kraken 2 takes both mates in one run (`--paired`).

| method | µs per read if the time covers 20 million reads (both mates) | µs per read if it covers 10 million reads |
|---|---:|---:|
| Cliffy, no digestion | 199–204 | 398–408 |
| Cliffy, DNA minimizers | 175–193 | 349–385 |
| Cliffy, minimizers | 76–82 | 152–164 |
| Kraken 2 | 19 | 37–39 |

- **The draft's 199–204 and 77–82 µs** are the first column, so they assume both mates are included. Its table says these times are for approximate listing, which matches the caption. (Rounding: 1524 s gives 76, against the draft's 77.)
- **This could matter for the comparison.** `silva.tex` line 254 gives 637 µs per read for Cliffy's LCA query on the server, about 3 times the first column, and that is a different query on a different machine, so it doesn't decide the question.
- **Machine:** an Intel Xeon Gold 6248R, 24 cores, 3 GHz, with 1500 GB of RAM, timed with GNU time.

**Threads.**
- **Kraken 2:** `--threads 1` (`workflow/rules/exp_3.smk` lines 271 and 289; preprint Appendix A.3).
- **Cliffy:** no thread count anywhere. `exp_3.smk` runs `pfp_doc64 run` with `--ref`, `--pattern`, `--output`, `--taxcomp`, `--ftab` and `--num-col`. Its option table (`src/pfp_doc.cpp` lines 535–546, at 3763529) has no thread option, and there are no `std::thread`, `pthread` or OpenMP calls in `src/` or `include/`. So the query is single-threaded, at least in this version of the code.
- **Bracken's time:** `exp_3.smk` also writes Kraken 2's time plus Bracken's to `output.total.time`, so Kraken 2's plotted time may include Bracken. That isn't settled either.

## The sparse dial (2026-11-05)

Is EM's lead on the aquatic reads, and proportional credit's, due to the sample having few genera, or to which genera they are?

**Design** (`scripts/sparse_genera.py`, `sparse_draws.py`, `experiments/aquatic_art_profiles.sh` with `AB=`; `results/silva/sparse/`; sampled reads on the external drive, `moved_2026-11-05/sparse_samples`):
- The aquatic sample's 100 abundances, in order, go to 100 random genera. Then MicrobeMixer at 6409f5d, MSv1, the same regions, 10 million pairs per region, every 50th pair kept, `find --preset silva`, every scheme.
- 3 draws, each with its own simulation, of two sets. **uniform:** any usable genus. **small:** usable genera with N_g ≤ 2, which are SILVA's smallest 3,974 genera (the third, with ties).
- **Deviation from "drawn from all 9,118 genera":** each region gets its own draw from the genera that can supply reads in that region. Of 8,721 uniquely named genera, the usable ones are 3,002 in V1–V2, 3,781 in V3–V4, 5,009 in V4 and 8,135 in V4–V5, and a genus outside these gets no reads. MicrobeMixer would spread such a genus's reads evenly over the others, destroying the table. Requiring all four regions would leave about 3,000 genera, with median N_g 21 against 4 overall. So the genera differ by region. Usable means over 55% of its amplicons are at least 250 bases and only A, C, G and T: ART silently makes no reads from an amplicon with an N.
- **Median N_g:** aquatic 111.5; the uniform draws range from 4 to 24 by region and draw, the smallest being V4–V5; small 1 to 2.
- **Limit of the small set:** only 151 such genera are usable in V1–V2, so its draws share about 70 of their 100 genera there.

Strict genus accuracy, mean of the four regions, mean of the draws (range). No read is from a genus with an ambiguous name, so that exclusion changes nothing here:

| scheme | aquatic | uniform | small |
|---|---:|---:|---:|
| range | 86.68 (86.45–86.95) | 84.72 (83.36–86.23) | 84.40 (82.98–86.77) |
| one label per MEM | 92.01 (91.81–92.17) | 86.99 (85.84–88.40) | 80.46 (79.36–82.09) |
| sample (simulated cliffs) | 90.34 (90.22–90.46) | 88.07 (86.64–89.07) | 86.73 (85.08–89.45) |
| even | 90.20 (90.06–90.28) | 88.05 (86.50–89.20) | 86.71 (84.80–89.38) |
| normalised | 92.66 (92.52–92.76) | 92.10 (91.23–92.84) | 91.22 (89.65–93.59) |
| proportional | 93.59 (93.44–93.77) | 88.98 (88.05–90.23) | 83.78 (82.02–86.64) |
| prior mixture, β = 0.5 | 93.47 (93.28–93.63) | 91.22 (90.18–91.84) | 89.84 (88.32–92.70) |
| EM, 20 iterations | 96.32 (96.13–96.43) | 96.80 (96.47–97.15) | 95.40 (94.42–96.64) |
| naive Bayes, flat, ε = 10⁻⁶ | 94.57 (94.41–94.70) | 94.69 (93.75–95.67) | 94.54 (92.70–95.68) |
| `bayesauto` | 93.95 (93.74–94.12) | 94.45 (93.67–95.26) | 94.94 (93.40–96.05) |

The aquatic column is the 3 MSv1 repeats, so its range is simulation noise alone. The other two also vary by which genera were drawn, so their ranges are wider. Other schemes and each region are in `sparse_tables.txt`.

- **EM's lead comes from sparsity, not from the aquatic genera:** it stays at 95.4 to 96.8 on genera of any size. On the dial, with thousands of genera, it was 90.2 to 92.3.
- **Proportional credit's advantage on the aquatic reads comes from those genera being large:** it falls from 93.59 to 88.98 to 83.78 as N_g falls, with the same abundances and number of genera.
- **Normalised credit and naive Bayes are stable:** normalised falls by 1.4 points across the sets, and flat naive Bayes stays at 94.5 to 94.7.
- **Range falls by 2 to 2.3 points** from the aquatic genera to random ones.
- **`bayesauto`:** ε is 0.5 or 0.6 in all 24 runs, none at the edge.
- **Regions vary a lot between draws,** especially V4 and V4–V5 (range in V4, small: 74.89 to 87.45), because a few genera get most of the reads and each draw's top genera differ.

## Naive Bayes with an EM-estimated prior (2026-11-05)

`score_silva.py --mode bayesem`. The prior starts flat. Each iteration computes every read's posterior from the likelihood and sets each genus's prior to its share of the posteriors. No new simulation; listings were regenerated with `find` from the sampled reads where they had been deleted, and all 44 reproduced their stored proportional scores exactly.

- **Likelihood:** naive Bayes with ε = 10⁻⁶, and also (`--mixture-eps`) `bayesauto`'s mixture likelihood at each read set's fitted ε. Its ε belongs to that likelihood; 0.5 in the plain formula would make a missing MEM more likely than a present one.
- **Details:** a read's posterior is over the genera of its MEMs, as in `bayes`. The prior is over all genera, "uncultured" included; calls are among the scored genera.
- **Checks:** 0 iterations gives `bayesflat` and `bayesauto`; `test_score_schemes.py` compares 5 and 20 iterations with an independent implementation (full vectors over the genera), and a hand case shows the prior deciding a tie. `dial_toy.sh` compares calls with another independent implementation on the toy's real listings: 2,400 pairs at each α, for both likelihoods, 0 disagreements.

Strict genus accuracy, mean of the four regions, 20 iterations, mean of the repeats or draws (full table, with 5 iterations: `bayesem_tables.txt`):

| read set | proportional | normalised | EM | NB, ε = 10⁻⁶ | `bayesauto` | NB + EM prior | mixture + EM prior |
|---|---:|---:|---:|---:|---:|---:|---:|
| aquatic, MSv1 | 93.59 | 92.66 | **96.32** | 94.57 | 93.95 | 95.94 | 95.59 |
| aquatic, MSv3 | 93.28 | 92.06 | **96.03** | 93.98 | 93.22 | 95.57 | 95.18 |
| dense dial, α = 1 | **92.68** | 91.76 | 92.31 | 92.66 | 91.42 | 92.57 | 91.55 |
| dense dial, α = 0 | 88.19 | 91.85 | 90.20 | 94.23 | **94.26** | 93.94 | 93.93 |
| dense dial, α = −1 | 85.01 | 91.84 | 90.22 | 94.78 | 95.19 | 94.91 | **95.27** |
| sparse, uniform | 88.98 | 92.10 | **96.80** | 94.69 | 94.45 | 96.22 | 96.17 |
| sparse, smallest | 83.78 | 91.22 | 95.40 | 94.54 | 94.94 | 96.51 | **96.71** |

Excluding reads from genera with ambiguous names, on the dense dial (the other sets have none): NB + EM prior 91.84, 93.53, 94.77 at α = 1, 0, −1, and mixture + EM prior 90.71, 93.52, 95.14.

- **The EM prior helps on sparse samples and does little on dense ones.** Against naive Bayes with a flat prior it gains 1.4 to 2.0 points on the aquatic and sparse sets, and changes the dense dial by −0.3 to +0.1.
- **It has the highest mean over these seven rows,** 95.09, against 94.91 for the mixture with an EM prior, 94.21 for flat naive Bayes and 93.90 for plain EM. Its lowest entry, 92.57, is just below flat naive Bayes's 92.66, and well above plain EM's 90.20. This is a mean over sets chosen here, not over any real distribution of samples.
- **Plain EM stays best on the aquatic reads** by 0.4 points, and on the uniform draws by 0.6, but loses 3.7 to 4.7 points on the dense dial at α = 0 and −1. On the smallest genera it is 1.1 points behind.
- **Convergence:** 5 iterations are within 0.05 points of 20.
- **The prior learned doesn't go to zero outside the sample:** on V1–V2 about 3,700 genera keep a prior above 0, against about 100 in the sample. Whether a sparser estimate would help wasn't tested.
- This prior is learned from the reads being classified, like plain EM's, so it needs the whole read set before calling any read, unlike `find --classify bayes`.

## RLBWT with sampled tags (2026-11-06)

The default undigested index at L = 30, `find -L 30 -c -F 10`, with backward steps on the run-length BWT and no RLCSA, against the explicit RLCSA with full tags (`experiments/tag_sampling.sh`, `tag_sampling.txt`). The protocol is that of the timing repeats: 5 repeats, one thread on CPU 9, all 1.6 million reads, with the configurations taking turns. All 160 outputs are byte-identical to the RLCSA's.

| configuration | tags (MB) | whole index, files loaded (MB) | mean per read (µs), median (range) | total / reads (µs) | peak (GB) |
|---|---:|---:|---:|---:|---:|
| RLCSA, full tags (default) | 836.4 | 1,556.9 | 45.78 (45.60–46.27) | 47.20 | 1.58 |
| RLBWT, full tags | 836.4 | 991.7 | 178.58 (178.22–179.03) | 179.89 | 1.02 |
| RLBWT, s = 2 | 588.2 | 743.5 | 186.50 (186.09–187.09) | 187.84 | 0.77 |
| RLBWT, s = 4 | 430.3 | 585.5 | 202.16 (202.03–202.35) | 203.49 | 0.61 |
| RLBWT, s = 8 | 351.2 | 506.5 | 235.26 (234.65–235.88) | 236.62 | 0.53 |

- **Files loaded:** the RLCSA row counts `.rlbwt`, `.csa` and `.tag`, as `find` loaded all three in this build; the others `.rlbwt` and the tags.
- **The RLCSA is 3.9 times faster** than the RLBWT with the same tags, at 1.6 times the space.
- **Sampling is cheap at first:** s = 2 costs 4% over full tags and saves 25% of the index; s = 4, 13% and 41%; s = 8, 32% and 49%. The tags' size flattens because the run boundaries stay.
- **Repeats vary by under 1%,** and the RLCSA row agrees with the 45.57 µs of the earlier repeats.
- **Code path:** with s = 1, `tag_array::tag()` reads the label directly; with s > 1 it follows LF through the RLBWT until it reaches a sampled run. The rest of the search is the same.
- `bac.rlbwt`, rewritten by `index` when s = 2 and 8 were built, kept its SHA-256.

## The lookup table at every step (2026-11-06)

`find --absent` (a), `--narrow` (b) and `--compact`, on the default undigested index at L = 30 (`experiments/narrow_timing.sh`, `narrow_tables.txt`). 5 repeats, one thread on CPU 9, all 1.6 million reads, configurations taking turns. All 920 outputs are byte-identical to `find --preset silva`, and the self-test passes on every mode.

- **(a)** A t-mer that doesn't occur gives the length of its longest occurring suffix, so BML with k < L gets its failure point from one lookup. It costs no space in the plain table, whose entry is unused.
- **(b)** Our binary searches run over c's runs, not its rows, so (b) needs each k-mer's run bounds as well as its interval.
- **`--compact`:** a bitvector over the k-mers with rank; Elias-Fano lists of starts and ends; 4 bits per absent k-mer; and, for (b), two Elias-Fano lists of the run bounds, which increase in k-mer order once each first base is offset by where its runs begin.

Mean per read (µs), median of 5; ranges are in `narrow_tables.txt` and are under 2%.

| | k = 8 | k = 10 | k = 11 | k = 12 |
|---|---:|---:|---:|---:|
| k-mers occurring | 65,536 (100%) | 1,048,575 (100%) | 4,175,610 (99.6%) | 13,872,668 (82.7%) |
| plain table (MB), without (b) | 1.11 | 17.83 | 71.30 | 285.21 |
| plain table, with (b) | | 26.21 | | |
| compact table (MB), without (b) | 0.29 | 3.57 | 12.16 | 36.48 |
| compact table, with (b) | 0.50 | 5.87 | 19.24 | 54.06 |
| **current search, plain table** | **52.13** | **47.13** | **45.16** | **43.98** |
| plain, (a) | | 47.10 | | |
| plain, (b) | | 52.16 | | |
| plain, (a) and (b) | | 52.16 | | |
| compact table | 53.24 | 49.20 | 48.12 | 47.68 |
| compact, (a) | 53.20 | 49.48 | 48.28 | 47.09 |
| compact, (b) | 77.26 | 82.11 | 83.81 | 87.88 |
| compact, (a) and (b) | 77.47 | 82.04 | 83.95 | 87.37 |

- **(b) makes search slower:** by 11% with the plain table at k = 10, and by 48% to 100% with the compact one. The jump table already leaves about 8 runs per search.
- **(a) changes little:** nearly every k-mer occurs up to k = 11. At k = 12, with 17% absent, it saves 0.6 µs.
- **The compact table is 4 to 8 times smaller,** and costs 1 to 4 µs per read on the first lookup alone.
- **A larger k helps the current search,** from 52.13 to 43.98 µs, at the plain table's cost in space; compact at k = 12 is 36 MB and 47.68 µs.
- **Peak memory** with the compact table counts the plain table, which is built first and then converted. That is a build cost of this implementation.

**Probes per backward step, explicit RLCSA** (a build with `-DKATKA_COUNT_PROBES`; steps in millions, over all 1.6 million reads):

| k | current | (a) | (b) | both | steps, current | steps, both |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 8.64 | 8.64 | 8.59 | 8.59 | 814.6 | 814.6 |
| 10 | 8.46 | 8.46 | 7.95 | 7.95 | 753.1 | 753.1 |
| 11 | 8.37 | 8.37 | 7.18 | 7.17 | 722.4 | 722.3 |
| 12 | 8.32 | 8.30 | 6.45 | 6.42 | 696.7 | 690.6 |

This is why (b) doesn't pay: the jump table already leaves about 4 probes per binary search, two searches per step, and narrowing saves at most 1.9 of the 8.3, against an added table access on every step.

**Elias-Fano RLCSA** (`index --efcsa`; one run each, µs per read): current search 220.1, 202.3, 194.7, 188.8 at k = 8, 10, 11, 12, and compact with (a) and (b) 376.1, 296.5, 276.2, 265.1. All 1,120 outputs across the timing, Elias-Fano and probe runs are byte-identical. Here (b) is 40 to 71% slower: narrowing replaces Elias-Fano's own bucket search with binary search by `access`.

**Baseline: this code change slowed the current search by 3.5%.** `find --preset silva` with old (0a66fa9) and new binaries in turn, 5 repeats: 45.54 µs against 47.13 µs (`regress.txt`). The cause was refactoring `run_lists::rank` through two helper functions; undoing just that gave 47.01 against the old binary's 46.92 on one file, and forcing the helpers inline made it worse (48.16), while keeping the RLBWT loaded made no difference. With `rank` written out again the default is 45.69 µs (next section), and the self-test, `dial_toy.sh`, `test_score_schemes.py` and the block checks all pass. **The table above was measured before the fix,** so its times are about 3.5% high; comparisons within the table still hold, since every configuration ran on that binary.

## RLCSA with sampled tags, no RLBWT (2026-11-06)

`find --preset silva -T ...` (`experiments/csa_sampled_tags.sh`, `narrow_tables.txt`), 5 repeats, one thread, pinned. Unsampled tags are recovered by LF from the RLCSA alone: the character whose runs hold the row, tried by decreasing frequency, then C[c] + rank_c. `find` no longer loads the RLBWT with an RLCSA, except for two-level search, STPD and locating. All 232 outputs are byte-identical to the default's. Index size is `.csa` + tags.

| tags | boundaries (MB) | labels (MB) | tags (MB) | whole index (MB) | mean per read (µs), median (range) | peak (GB) |
|---|---:|---:|---:|---:|---:|---:|
| full (default) | 207.6 | 628.9 | 836.4 | **1,401.6** | 45.69 (45.35–46.08) | 1.43 |
| s = 2 | 207.6 | 380.7 | 588.2 | 1,153.4 | 50.79 (50.60–51.07) | 1.18 |
| s = 4 | 207.6 | 222.7 | 430.3 | 995.5 | 59.33 (59.16–59.40) | 1.02 |
| s = 8 | 207.6 | 143.7 | 351.2 | 916.5 | 78.13 (77.97–78.39) | 0.94 |
| grammar, built-in RePair | 207.6 | 242.2 | 449.7 | 1,015.0 | 53.03 (52.89–53.46) | 1.04 |
| grammar, bigrepair | 207.6 | 302.0 | 509.5 | 1,074.8 | 55.83 (55.56–56.13) | 1.10 |

- **The default's footprint:** 1.40 GB of files and a 1.43 GB peak, down from 1.56 and 1.58 GB.
- **Grammar tags beat sampling here:** the built-in RePair's array is about as small as s = 4 (449.7 against 430.3 MB) and 11% faster than it, at 53.03 µs against 59.33.
- **An index with these tags takes** 1.40 GB at 46 µs, 1.02 GB at 53 µs, and 0.92 GB at 78 µs.
- These times are on the fixed build (below), rerun after the slowdown was found.
- **Boundaries are the same** in every row; only the labels change.

## `find --classify bayesem` (2026-11-06)

Naive Bayes with an EM prior as a classification mode of `find`, ε = 10⁻⁶ by default and 20 iterations, with `--mate2` and `--exclude-docs`. It keeps each read's log-likelihoods, runs EM in the scorer's chunks and order, and writes the calls at the end.

**Calls are identical to `score_silva.py --mode bayesem --calls`** on the four aquatic regions and on the dial's V4_V4 at α = 1 and −1, 200,000 pairs each, and on the toy at α = 1, 0, −1 (`tests/dial_toy.sh`). Each real set takes 3.4 to 5.4 s with 12 threads, in 1.5 to 1.8 GB. The match on exact ties assumes the same `-j` as the scorer's `--threads`, because the sum order depends on it.

## Movi 2 for MEM-finding, and sampled tags with its LF (2026-11-06)

Movi 2 (v2.0.0 + 12 commits, 9bc0acc), `regular` mode, with the prototype in `movi2/` (`prototype.patch`, scripts, `sweep.txt`). Nothing here changes `find`.

**A shared BWT.** Our text with `X` replaced by `%`, Movi 2's separator, which must sort first; only the rows change. Big-BWT's BWT, with byte 1 turned into Movi 2's terminator, byte 0, goes to both tools. Tag runs: 359,322,005 against 359,351,512. Movi 2 builds its index from it in 22 s: 68,329,725 rows, 546.6 MB. Movi 2's own text builder isn't used: it turns non-ACGT into `A`, and orders strands per record.

**Identical answers.** Over all 8 read files, 200,000 reads each, all 8,532,597 MEMs, their intervals as rows, and their tag answers with full tags are identical to `find`'s; so are the tag answers at every s and both grammars (first repeat of each).

**What it took:**
1. `complement()` (`src/utils.cpp`) tests for `#`, while the separator is `%` (`SEPARATOR`, `commons.hpp`), so `%` gets complement `A`, and `extend_bidirectional` misplaces the reverse complement's interval when separator rows fall in the interval. MEMs then end early: one read gave (0, 38) where the text check confirms (0, 49). A one-line fix: `movi2/separator.patch`. This affects indexes built with `--separators`.
2. `query_mem_bml` has no alphabet check before `initialize_backward_search` at the end of a MEM, so a read with an ambiguity code at that point indexes `first_runs` with `alphamap` = 256 and crashes. 1% of these reads have such codes, copied from SILVA. `query_backward_search` has the check. Worked around here by splitting reads at non-ACGT, as `find` does; not patched.
3. At query time `-l` doesn't set the MEM length; `--min-mem-length` does.

**Bidirectional extension costs too much on SILVA.** Unmodified `--mem`, with the fix in (1), takes about 70 s per file. `extend_left` walks 466 runs per step on V3_V4_1 and 1,204 on V4_V5_2, 14.5 and 44.3 billion loop steps; `count()` walks only 1.8 and 2.3 million. So BML here uses Movi 2's backward steps alone: the window on the forward strand, then a fresh search of the reverse complement only when the window matches, as `find` does. That takes about 6.7 s per file, with the same MEMs. How this compares on narrow intervals was not measured. 

**Search, 5 repeats, µs per read** (total processing time over the reads, one thread, pinned, without loading):

| | median (range) |
|---|---:|
| Movi 2, MEMs with counts | 33.45 (32.92–34.50) |
| Movi 2, intervals as rows, one line per MEM | 39.82 (39.62–40.25) |
| ours, RLCSA, one line per MEM | 41.14 (41.08–41.20) |

So with the same output the two are within 3%. The first row doesn't convert to rows.

**With tags, 5 repeats** (`movi2/select.txt`). The first run of this table was wrong: the prototype's setup (loading the tags, building the row list) fell inside Movi 2's clock, about 1.4 s per file or 7 µs per read. These times subtract it, and replace that table. The RLCSA, rows at full tags and the select version took turns; the other row-list configurations ran directly afterwards. µs per read, median (range), one thread, pinned. All 80 first-repeat comparisons in the rotation are identical to `find`'s; the row-list tag arrays were compared in the earlier run.

| tags | RLCSA | Movi 2, row list | Movi 2, select (below) |
|---|---:|---:|---:|
| full | 46.70 (46.41–46.94) | **44.73** (44.63–46.00) | 49.99 (49.67–50.28) |
| s = 2 | 51.64 (51.58–52.41) | 59.48 (59.25–60.86) | 82.87 (82.37–83.33) |
| s = 4 | 60.55 (60.31–60.86) | 88.98 (88.83–90.54) | 149.40 (149.06–150.67) |
| s = 8 | 80.02 (79.94–80.19) | 156.44 (153.75–160.32) | 298.10 (295.11–304.97) |
| s = 16 | | 297.93 (297.60–309.14) | |
| built-in grammar | 53.98 (53.88–54.37) | **53.35** (53.10–55.22) | 58.12 (58.02–59.39) |
| bigrepair | | 55.86 (55.71–56.99) | |

Whole index (MB), files loaded. Boundaries are 207.5 MB in the first two columns. RLCSA: `.csa` and tags, its table built at load. Row list: Movi 2's index 546.6, ftab 33.6, row starts 58.5, tags. Select: Movi 2's index, ftab, B and E 308.0, labels.

| tags | labels | RLCSA | Movi 2, row list | Movi 2, select |
|---|---:|---:|---:|---:|
| full | 628.8 | 1,401.5 | 1,475.1 | 1,517.0 |
| s = 2 | 380.6 | 1,153.4 | 1,226.9 | 1,268.8 |
| s = 4 | 222.7 | 995.4 | 1,068.9 | 1,110.9 |
| s = 8 | 143.7 | 916.4 | 989.9 | 1,031.9 |
| s = 16 | 104.2 | | 950.5 | |
| built-in grammar | 242.0 | 1,014.7 | 1,088.3 | 1,130.2 |
| bigrepair | 302.0 | | 1,148.2 | |

- **Full tags:** Movi 2 with the row list is 4% faster than the RLCSA, at 5% more space. With grammar tags the two are level (1%), at 7% more space.
- **Sampled tags:** the RLCSA wins clearly, by 15% at s = 2 and 2× at s = 8. Each LF step on Movi 2 adds a rank and an access on the row list.
- **Peak memory** for the row-list prototype includes a temporary 547 MB array in its loader; it isn't the index's footprint.

**Tag lookup by select, with no row numbers** (`KATKA_SELECT` in `prototype.patch`). B has one 0 per Movi run start and one 1 per tag run start, in row order, BWT starts first on ties; E is the Elias-Fano list of each entry's row. A point query takes two select₀ and scans run i's tag starts; an interval query adds the tag runs between. The tie bit T (about 8.5 MB) wasn't needed: the scan finds a tied tag run at offset 0. It's not counted.

- **Answers are identical:** all 8 files, first repeat, at every tag array.
- **Size:** B is 80.2 MB and E 227.8 MB, 308.0 MB in all, against 266.1 MB for row starts and boundaries. So it's 42 MB larger, not smaller: merging the lists adds a bit per entry, and E must hold every entry.
- **Scan:** 4.99 tag starts per point query for lookups, as predicted, and 7.2 to 7.6 in LF steps, which land more often in long runs.
- **Time:** slower than the row list everywhere: by 12% at full tags, 9% with grammar tags, and 39 to 91% with sampling. A point query costs two selects plus a select in E per entry scanned, where the row list costs a rank and an access.
- **Sequential scan, tried afterwards** (`movi2/seqscan.txt`, in rotation with the RLCSA and the row list, 5 repeats, 48 comparisons identical): one select on B and one on E, then each next entry from E's high bits, up to the next 0 of B, with binary search after 16 entries. Sizes are unchanged. µs per read, median (range), against 46.68 for the RLCSA and 44.76 for the row list in the same rotation:

| tags | select, sequential scan | select, one select per entry |
|---|---:|---:|
| full | **42.83** (42.65–44.01) | 49.99 |
| built-in grammar | **50.76** (50.47–52.42) | 58.12 |
| s = 2 | 74.60 (74.17–75.84) | 82.87 |
| s = 4 | 136.54 (135.55–138.46) | 149.40 |
| s = 8 | 267.75 (266.10–274.75) | 298.10 |

  This makes it the fastest configuration at full tags and with grammar tags, 4% ahead of the row list at 42 MB more. Sampled tags stay well behind: 73% of LF-step lookups on one file (s = 4) fell in runs of over 16 tag starts, where the scan gives way to binary search.

## Cliffy's tie rule (2026-11-07)

**What the code does.** `cliffy-experiments/src/classify_silva_readset.py` counts a read as correct if the true clade is among the clades of maximum weight (lines 336–340 in `classify_method0`, 430–434 in `classify_method1`, and the same in `classify_method2`). `scripts/cliffy_listing_score.py` follows that as `theirs=`, and also reports `strict=`, with ties wrong; `score_cliffy_listing.sh` prints `theirs` and leaves `strict` in the CSVs. Reads with no match also count as correct under their rule, unless the genus is "uncultured".

**No Cliffy listings are on the desktop or the external drive,** so they couldn't be rescored. They are on the server, from Gonzalo's run (`score_cliffy_listing.sh RESULTS_DIR` gives `theirs` and `strict` per region). From the draft, for the LCA query on our reads: 85.75 strict and 86.45 under their rule, undigested; 84.15 and 85.53 digested.

**A discrepancy in the draft.** `classify_method2`, Cliffy's approximate-listing score, adds the match's **whole length** to every listed genus (`weight_dict[taxa_hit][0] += length`). The draft (`silva.tex` lines 38 and 42) describes |M|/|D′| and |M|/|D|, and our `cliffs` and `list` schemes split the credit that way. So those schemes aren't Cliffy's listing score. `score_silva.py --mode cliffsfull` and `listfull` now apply method2's rule.

Aquatic reads, MSv1.1, from `find --preset silva` at L = 30 (`tie_rule.txt`, with each region). Means of the four regions, as %; ties are reads whose top weight is shared:

| scheme | Cliffy's rule | strict | gain | reads tied |
|---|---:|---:|---:|---:|
| range | 87.60 | 86.95 | 0.65 | 0.77 |
| sample, |M|/|D′| | 91.03 | 90.46 | 0.57 | 0.73 |
| even, |M|/|D| | 91.83 | 90.28 | 1.55 | 1.60 |
| sample, full length (method2's credit) | 91.41 | 90.51 | 0.90 | 1.76 |
| every genus, full length (method2's credit) | 95.97 | 92.50 | 3.47 | 3.60 |
| normalised | 92.76 | 92.76 | 0.00 | 0.01 |
| proportional | 93.79 | 93.77 | 0.02 | 0.04 |
| EM, 20 iterations | 96.39 | 96.39 | 0.00 | 0.00 |
| naive Bayes, ε = 10⁻⁶ | 94.70 | 94.70 | 0.00 | 0.01 |
| naive Bayes, EM prior | 96.06 | 96.07 | 0.00 | 0.00 |

- **The rule matters most where credit doesn't discriminate,** and full-length credit ties most often. Ties reach 5.7% on V4 with exact listings.
- **Our best schemes don't depend on the rule.** They hardly ever tie, so they can be quoted under either.
- **These are MEMs at L = 30 with simulated lists,** not Cliffy's ZMPs and sampled profiles. They show the rule's effect, not Cliffy's own numbers. The small negative gap on the last row is rounding: `theirs` has 4 decimals.

## Naive Bayes with an EM prior learned from a sample (2026-11-07)

`score_silva.py --mode bayesem --em-sample N --em-smooth ... --em-seeds 1,2,3` (`em_sample.sh`, `em_sample.txt`). The prior is learned by EM, 20 iterations, from N random read pairs, then every read is classified with it. Smoothing: `p<c>` adds c pseudo-reads to every genus at each iteration; `m<w>` mixes the final prior with a flat one at weight w. ε = 10⁻⁶. Sets: the aquatic reads (MSv1.1, 4 regions), and the sparse dial's random and smallest genera (3 draws each). The range is over the 3 seeds.

Strict genus accuracy, mean of the four regions, and the change against the full-set prior (aquatic 96.07, random genera 96.22, smallest genera 96.51):

| reads | smoothing | aquatic | random genera | smallest genera |
|---:|---|---:|---:|---:|
| 1,000 | none | 96.87 (+0.80) | 96.99 (+0.78) | 96.43 (−0.08) |
| 1,000 | 0.5 pseudo-reads | 95.63 (−0.44) | 95.86 (−0.35) | 96.07 (−0.44) |
| 1,000 | 1 pseudo-read | 95.49 (−0.58) | 95.72 (−0.50) | 96.03 (−0.48) |
| 1,000 | flat, 1% | 96.43 (+0.36) | 96.57 (+0.36) | 96.32 (−0.19) |
| 1,000 | flat, 10% | 96.21 (+0.15) | 96.36 (+0.15) | 96.24 (−0.28) |
| 10,000 | none | 96.33 (+0.27) | 96.37 (+0.16) | 96.44 (−0.07) |
| 10,000 | 1 pseudo-read | 95.88 (−0.19) | 96.05 (−0.17) | 96.31 (−0.20) |
| 10,000 | flat, 1% | 96.16 (+0.09) | 96.26 (+0.05) | 96.39 (−0.12) |
| 100,000 | none | 96.09 (+0.02) | 96.20 (−0.01) | 96.49 (−0.02) |
| 100,000 | flat, 1% | 96.06 (−0.01) | 96.19 (−0.03) | 96.48 (−0.03) |

Every combination, with seed ranges, is in `em_sample.txt`. Seeds differ by at most 0.25 points at 1,000 reads and 0.1 at 10,000.

- **A small sample is enough, and smoothing doesn't help.** Unsmoothed, 1,000 pairs gain about 0.8 points on the aquatic and random-genera sets and lose 0.08 on the smallest genera. 10,000 pairs are within 0.3 everywhere, and 100,000 within 0.02. Every smoothing setting is worse than none at the same size.
- **Why fewer reads can do better:** the full-set prior keeps 3,709 genera above 0 on aquatic V1_V2, while 1,000 reads give 1,475 to 2,031. Spurious genera that the full set props up get no weight. Smoothing adds them back. This is one region; the others weren't counted.
- **As the streaming default for `find --classify bayesem`:** yes. A prior from about 10,000 read pairs, unsmoothed, then every read classified in one pass, is within 0.3 points of the full-set prior on all three kinds of sample, and needs memory for 10,000 pairs instead of the whole set.
- **The sample must be random.** These read files are grouped by genus (MicrobeMixer writes each genus's reads together), so the first 10,000 reads of a file would be a few genera only. A reservoir sample, or a first pass, is needed. That was not tested here.

## Facts for the paper (2026-11-07)

- **Genera:** `work/ref/done` says `9118 genera`, and `ref/rna/doc_to_traversal.txt` has 9,118 lines, as in Cliffy's paper.
- **Reference (rechecked 2026-11-07):** `write_silva_genera.py` logs only the number of genera, so the counts come from its input (`SILVA_138.1_SSURef_NR99_tax_silva.fasta`) and its output (`work/ref/rna/doc_*_seq.fa`). 510,508 sequences in, **402,378 included** (586,483,400 bases), **108,130 excluded** (160,907,699 bases), as expected. Every included header is distinct. Excluded are the sequences whose lineage has no genus-rank taxon: 68,572 Bacteria, 29,144 Eukaryota, 10,414 Archaea. Included are 363,003 Bacteria, 29,400 Eukaryota, 9,975 Archaea. The index's text has 1,173,771,556 characters = 2 × (586,483,400 + 402,378): both strands, one separator per sequence. 171,854 of the included bases are not A, C, G or T. `reference_counts.txt`.
- **Cliffy's minimizer-digested indexes:** not on the desktop, and `cliffy_index_size.sh` has not been run on them here (they were built on the server). What we have is Cliffy's Table 1 (preprint and journal): 4.266 for minimizers and 7.858 for DNA minimizers, against 9.014 undigested. By `cliffy_index_size.sh`'s description of cliffy-experiments' rule, those are GiB (2³⁰ bytes), about 4.58, 8.44 and 9.68 GB. The draft's 7.63 GB for Cliffy's digested index (`silva.tex` 255 and 276) is from the server run and couldn't be traced here; it isn't one of Table 1's figures. Gonzalo's run of `cliffy_index_size.sh` on the server would give both counts.
- **KATKA's build, Movi 2 with the row list** (`build/`, one machine, each step in an 11 GB cgroup, timed alone; i7-8700):

| step | full tags | grammar tags | peak (GB) |
|---|---:|---:|---:|
| `katka prep` (text, both strands, `X` separators) | 4.9 s | 4.9 s | 0.05 |
| separator `X` to `%` (prototype only) | 1.2 s | 1.2 s | 0.01 |
| Big-BWT | 324.1 s | 324.1 s | 5.33 |
| tag runs and tag array (`index --tags --no-rmq`) | 329.0 s | | 5.35 |
| … with the built-in RePair (`--grammar`) | | 1,454.6 s | 9.03 |
| terminator byte (prototype only) | 2.2 s | 2.2 s | 0.01 |
| Movi 2 index and ftab (`movi-regular build --ftab-k 10`) | 22.1 s | 22.1 s | 2.42 |
| row list, built at load | 0.6 s | 0.6 s | |
| **total** | **684 s (11.4 min)** | **1,810 s (30.2 min)** | **5.35 / 9.03** |

  - **Tag runs:** the grammar step computes them again, 309.5 s of its 1,454.6. From a cached `.tagruns` it would take about 1,145 s, for a total of about 25 min.
  - **Movi 2's exit status:** `movi-regular build` stores the index, then exits with signal 11 while generating null statistics, because the prototype passes an empty `ref.fa`. The index files are complete and were used for all the results here.
  - **The two prototype-only steps** go away if `prep` writes `%` and Big-BWT's terminator is written as byte 0.
  - **Querying,** the whole index takes 1,475 MB (full) or 1,088 MB (grammar) of files, from the Movi 2 section above.

## Held out: reads from sequences not in the index (2026-11-07)

`heldout/` (`heldout.py`, `heldout2.sh`, `split.tsv`, `dropped.txt`, the CSVs). For each of the 100 aquatic genera, a random 10% of its SILVA sequences, at least 1, were held out: 5,977 sequences, against 53,344 left in the index. The other 9,018 genera are unchanged. The held-out index was built as the default (`index --preset silva`, 334 s). Reads were simulated as for the aquatic reads (MicrobeMixer at 6409f5d, MSv1, its 2% mutation, 10 million pairs per region, every 50th kept), but from the held-out sequences only: MicrobeMixer's `--silva-ref` is given a FASTA of just those, so MicrobeMixer itself is unchanged.

**Genera without reads.** In V1–V2, 12 of the 100 genera were left out of the abundance table, and MicrobeMixer spread their reads over the other 88. For 11 of them the held-out set had no V1–V2 region after primer extraction; for one, too few of its amplicons were usable by ART (250 bases or more, A/C/G/T only). The other three regions kept all 100 (`dropped.txt`). The first attempt, without that step, failed on MicrobeMixer's "not enough reads".

Strict genus accuracy, mean of the four regions. In-database is the aquatic reads (MSv1, mean of 3 repeats); held-out is these reads against the held-out index; the last column is the same reads against the full index, which contains their sequences:

| scheme | in-database | held-out index | full index, same reads | change, held-out |
|---|---:|---:|---:|---:|
| range | 86.68 | 76.85 | 87.41 | −9.83 |
| one label per MEM | 92.01 | 83.59 | 92.12 | −8.42 |
| sample | 90.34 | 81.07 | 90.72 | −9.27 |
| even | 90.20 | 80.91 | 90.58 | −9.29 |
| normalised | 92.66 | 84.57 | 93.09 | −8.09 |
| proportional | 93.59 | 85.40 | 93.61 | −8.19 |
| prior mixture, β = 0.5 | 93.47 | 85.03 | 93.41 | −8.44 |
| proportional/even mix, β = 0.5 | 93.07 | 84.23 | 92.95 | −8.84 |
| EM, 20 iterations | 96.32 | 88.04 | 96.01 | −8.28 |
| naive Bayes, ε = 10⁻⁴ | | 83.68 | 90.02 | |
| naive Bayes, ε = 10⁻⁵ | | 87.77 | 95.18 | |
| naive Bayes, ε = 10⁻⁶ | 94.57 | 87.66 | 95.18 | −6.91 |
| naive Bayes, ε = 10⁻⁷ | | 87.36 | 95.00 | |
| naive Bayes, ε = 10⁻⁸ | | 87.09 | 94.81 | |
| `bayesauto` | 93.95 | 86.87 | 94.51 | −7.08 |
| naive Bayes, EM prior | 95.94 | **89.03** | 96.08 | −6.91 |

- **Holding out costs 7 to 10 points for every scheme.** The same reads against the full index score about as the in-database reads do (within ±0.7), so the reads aren't harder: the loss comes from their sequences being missing.
- **The ranking holds.** Naive Bayes with an EM prior is best held out (89.03), then EM (88.04) and naive Bayes (87.66 to 87.77). Range loses most (−9.83).
- **The naive Bayes schemes lose least,** about 7 points against 8 to 10 for the credit schemes.
- **ε:** 10⁻⁵ and 10⁻⁶ are best held out, within 0.11 of each other. ε = 10⁻⁴ costs 4 points here, so the choice of 10⁻⁶ holds without its sequences in the index.
- **One repeat, one held-out split.** In-database repeats varied by about ±0.2, so differences under a point here aren't meaningful.

## Time and space, nine configurations in one rotation (2026-11-07)

The undigested index on the shared BWT (`%` separators), L = 30, the 8 aquatic files (200,000 reads each), one thread pinned to CPU 9, 5 repeats, all nine configurations taking turns within each repeat (`movi2/movi12.sh`, `movi2/ef.txt`). RLCSA: `find -L 30 -c --csa` or `--efcsa -F 10 -T TAG`, which no longer loads the RLBWT. Movi 2: the prototype's row list, with its setup subtracted. All 72 first-repeat outputs are identical to `find --preset silva`'s. µs per read is the whole processing time over the reads, without loading.

| configuration | index without tags (MB) | boundaries (MB) | labels (MB) | whole index (MB) | µs per read, median (range) |
|---|---:|---:|---:|---:|---:|
| explicit RLCSA, full tags | 565.2 | 207.5 | 628.8 | 1,401.5 | 46.79 (46.49–46.95) |
| explicit RLCSA, grammar tags | 565.2 | 207.5 | 242.0 | 1,014.7 | 54.16 (53.90–54.70) |
| explicit RLCSA, s = 4 | 565.2 | 207.5 | 222.7 | 995.4 | 60.56 (60.12–60.61) |
| Elias-Fano RLCSA, full tags | 185.3 | 207.5 | 628.8 | 1,021.7 | 199.87 (199.49–200.58) |
| Elias-Fano RLCSA, grammar tags | 185.3 | 207.5 | 242.0 | **634.9** | 208.73 (208.48–209.34) |
| Elias-Fano RLCSA, s = 4 | 185.3 | 207.5 | 222.7 | 615.6 | 262.02 (261.82–263.84) |
| Movi 2 with the row list, full tags | 638.7 | 207.5 | 628.8 | 1,475.1 | **44.85** (44.79–45.88) |
| Movi 2 with the row list, grammar tags | 638.7 | 207.5 | 242.0 | 1,088.3 | 53.01 (52.94–53.29) |
| Movi 2 with the row list, s = 4 | 638.7 | 207.5 | 222.7 | 1,068.9 | 91.87 (89.56–93.35) |

- **Index without tags:** Movi 2 is its index (546.6), ftab (33.6) and the row list (58.5). The RLCSAs' lookup table, 17.8 MB at k = 10, is built at load and isn't in a file.
- **Labels at s = 4** include the bitvector of sampled runs.
- **The Elias-Fano RLCSA is 3.05 times smaller than the explicit one** (185.3 against 565.2 MB) and 4.3 times slower with full tags. With grammar tags it makes the smallest index here, 634.9 MB, at 208.73 µs.
- **Movi 2 is 4% faster than the explicit RLCSA** with full tags, and 2% with grammar tags, at 73.6 MB more.
- **With s = 4,** Movi 2 is 52% slower than the explicit RLCSA, because each LF step on Movi 2 adds a rank and an access on the row list.
- **Elias-Fano costs much less with sampling, relatively:** s = 4 adds 31% to its time, against 29% on the explicit RLCSA and 105% on Movi 2.

## Where the Elias-Fano RLCSA's time goes (2026-11-07)

`experiments/rlcsa_profile.cpp`, `rlcsa_profile.txt`. `find` built with `-DKATKA_RECORD_STEPS` wrote every backward step (c, sp, ep) on the first 20,000 reads of V3_V4_1 (L = 30, full tags, `--efcsa -F 10`): 11.2 million steps. They were replayed on both RLCSAs, each part timed alone, best of 3, pinned to CPU 9.

**How our Elias-Fano RLCSA differs from the paper's.** In Brown et al.'s RLCSA (Theorem 7, "Faster run-length compressed suffix arrays"), interval ends are kept as an L-run and an offset; a step searches Ψ′ in run numbers, and a rank on the sparse bitvector B_L turns the new F position back into a run and an offset. B_FL replaces that rank. Ours keeps rows. `EFRank` stores, for each character c, the run starts shifted by c·n (`KS`) and their positions in F (`FP`), and a step is:
- **a predecessor search** in `KS` for sp and for ep: `count_less`, a rank on a sparse bitvector over a universe of σn = 3·10¹¹;
- **the first run of c**, `count_less(c n)`, the same for every step by c;
- **the rank from the run found**, `rank_from`: one `KS.access` and one `FP.access2` per end, and `FP.access(b)`, all selects on the Elias-Fano lists.

| part | ns per step | share |
|---|---:|---:|
| predecessor search: 2 `KS.count_less` | 180.3 | 46% |
| rank from the runs: 2 × (`KS.access`, `FP.access2`, `FP.access(b)`) | 188.4 | 49% |
| first run of c: `KS.count_less(c n)` | 19.7 | 5% |
| **whole step** (parts sum to 388.3) | **394.7** | |
| explicit RLCSA, whole step | 71.0 | |
| … of which the jump table and binary search | 69.3 | 98% |

- **Per search:** 80.7 entries in the high bucket (KS has l = 12 low bits over 68.3 million runs), and 8.6 low-bit reads.
- **Steps per read:** this file has 559 steps per read, so 395 ns per step is about 221 µs of its 275 µs per read: 80% is backward steps, and tags, the table and the rest come to about 54 µs.
- **What Theorem 7 removes,** the rank on B_L, corresponds to our predecessor search: 46%. But it applies only if interval ends are kept as runs and offsets. In our design, rows are what make the predecessor search necessary.
- **Half the time is selects,** six per step, of which two (`count_less(c n)` and `FP.access(b)`) depend only on c and could be stored per character.

## The Elias-Fano RLCSA, made faster (2026-11-07)

Two changes (`src/bits.hpp` `EFRank`, `src/rlcsa.hpp`), each profiled on the same 11.2 million recorded steps (`rlcsa_profile_cached.txt`, `rlcsa_profile_remapped.txt`). Answers are unchanged: `find --efcsa` gives the default's output on V3_V4_1 with the new file, and with a file in the old format, which still loads. The self-test passes.

1. **Per-character values.** c's first run in `KS` and its position in F depend only on c, so they're stored per character (5 KB) instead of a `count_less` and an `FP.access` on every step.
2. **`KS` over the characters that occur.** `KS` keyed run starts by c·n with c up to 256, so its universe was 256n. With 6 characters, that gave 12-bit low parts, and each high-bit bucket held about 81 runs. Keys are now slot(c)·n, slot(c) being c's rank among the characters that occur, so the universe is 7n. That gives 6-bit low parts and buckets of 2.6 runs, which is a finer sampling of the high bits, here at no cost: `KS` shrinks from 126.8 to 85.0 MB. New files carry the mapping (format 3 of `KTCSA001`).

| version | ns per step | predecessor search | rank from the runs | Elias-Fano RLCSA (MB) |
|---|---:|---:|---:|---:|
| before | 394.7 | 180.3 (46%) | 188.4 (49%) | 185.3 |
| 1. per-character values | 346.9 | 183.2 (52%) | 166.1 (47%) | 185.3 |
| 2. and `KS` over the characters that occur | **267.1** | **90.0** (34%) | 174.6 (65%) | **143.6** |
| explicit RLCSA | 69.7 | | | 565.2 |

- **Per search:** 2.6 entries in the bucket and 2.3 low-bit reads, against 80.7 and 8.6 before.
- **What's left is selects:** two `KS.access` and two `FP.access2` per step, 65% of it. The gap to the explicit RLCSA is 3.8×, and the search no longer dominates.

## Interpolative coding, bits per value (2026-11-07)

`experiments/interp_bits.cpp`, `interp_bits.txt`: code lengths only, on the explicit RLCSA's lists, per character. Midpoint codes recurse over the same nodes as interpolative coding: the middle value's deviation from the linear interpolation of the node's two known neighbours (their average when the middle is the centre), zig-zagged. Rice's parameter is 2^round(log₂(R / (2√(2k + 3)))), for a node of 2k + 1 values whose neighbours differ by R.

| list (68.3 M values each) | Elias-Fano | interpolative, minimal binary | interpolative, centred | midpoint, Elias γ | midpoint, Rice |
|---|---:|---:|---:|---:|---:|
| run starts (universe n) | 8.14 | 6.21 | **6.16** | 9.59 | 6.76 |
| F positions (universe of c's rows) | 6.05 | **3.26** | 3.31 | 5.47 | 4.23 |

- **Interpolative coding beats Elias-Fano by about 2 and 3 bits per value.** On both lists together that's 80.9 MB against 121.1 for these per-character lists. These have no select support; `KS` and `FP` as built are 143.5 MB.
- **Centred codes help slightly on run starts** (−0.05 bits) and hurt slightly on F positions (+0.05).
- **Midpoint deviation is worse than plain interpolative coding** with either code. With γ it's worse even than Elias-Fano on run starts. Values don't sit close enough to the linear interpolation: the deviations are about as large as the range allows, so a code that favours small deviations pays for it.
- **Random access:** interpolative coding gives none by itself, so it would need blocks with sampled values. These counts leave that cost out.

## Time and space, eleven configurations, the faster Elias-Fano RLCSA (2026-11-07)

`movi2/ef_table2.py`, `movi2/ef_table2.txt`; script `~/tool/silva/movi13.sh` (as `movi12.sh`, with more configurations). One rotation, 5 repeats, the 8 files of the aquatic reads (1.6 million reads), one thread, pinned to CPU 9, setup excluded. The Elias-Fano RLCSA is the 143.6 MB file of "The Elias-Fano RLCSA, made faster", before the select changes below. **Index sizes now include each index's lookup table:** the RLCSAs' 17.8 MB table (`-F 10`, built at load) and Movi 2's ftab (k = 10) and 58.5 MB row list. Every configuration has the same 207.5 MB of tag boundaries; the labels differ. Every output of the first repeat (88 files) is identical to `find --preset silva`'s.

| index | tags | whole (MB) | µs per read, median (range) | Pareto |
|---|---|---:|---:|:-:|
| explicit RLCSA (583.0) | full | 1,419.3 | 47.13 (46.90–47.58) | * |
| | s2 | 1,171.2 | 52.00 (51.81–52.46) | * |
| | grammar | 1,032.5 | 54.42 (54.31–54.71) | * |
| | s4 | 1,013.2 | 60.84 (60.65–61.07) | * |
| | s8 | 934.2 | 79.69 (79.55–80.38) | * |
| Elias-Fano RLCSA (161.4) | full | 997.7 | 137.64 (137.33–138.42) | |
| | grammar | 610.9 | 146.25 (145.75–147.29) | * |
| | s4 | 591.6 | 183.27 (183.22–184.98) | * |
| | s8 | 512.6 | 246.77 (246.50–248.25) | * |
| Movi 2, row list (638.7) | full | 1,475.1 | 44.80 (44.75–44.88) | * |
| | grammar | 1,088.3 | 53.29 (53.08–53.73) | * |

- **Against the old Elias-Fano RLCSA** (movi12, 185.3 MB): grammar tags 208.73 → 146.25 µs (−30%), full 199.87 → 137.64, s4 262.02 → 183.27.
- **Elias-Fano with full tags is off the frontier:** the explicit RLCSA with s8 tags is smaller (934 MB) and faster (80 µs).
- **The gap is still over 2×:** Elias-Fano with grammar tags is 2.7 times slower per read than the explicit RLCSA with grammar tags, for 41% less space.

## Elias-Fano selects: shared and finer samples (2026-11-07)

`src/bits.hpp` (`EliasFano`, `EFRank::rank_pair`); `ef_selects/*.txt`, from `experiments/rlcsa_profile.cpp` on the same 11.2 million recorded steps, best of 3, CPU 9, each variant compiled separately (`-DKATKA_EF_PAIR=0/1 -DKATKA_EF_SAMPLE1=0/64/128/256/512`). In every variant, every step's two ranks equal the explicit RLCSA's (the profile now checks them all).

- **Shared selects (`KATKA_EF_PAIR`):** a step ranks sp and ep. In 42.2% of the steps both lie in the same run of c, and in 11.6% in neighbouring runs. `rank2` then reads both runs' starts and F positions from one select in KS and one in FP (`access3`: the next values come from scanning to the next one-bits), instead of two each. The two predecessor searches share their select-zero when both keys are in the same high bucket (43.3%) or neighbouring ones (8.8%) (`count_less2`). `FP.access2` already read F's next value from the first select.
- **Finer samples (`KATKA_EF_SAMPLE1`):** the position of every k-th one of the high bits, 32 bits each, built at load (the file format is unchanged; not used if the high bits pass 2³²). Select-one starts from the sample's word and scans, with no binary search over the 512-bit blocks' counts. The old 256-one samples (64 bits each, followed by that binary search) are still used for select-zero.

| variant | ns per step | extra space (KS + FP) | Elias-Fano RLCSA |
|---|---:|---:|---:|
| before | 265.6 | | 143.6 MB |
| shared selects | 232.2 | 0 | 143.6 MB |
| samples every 512 | 259.6 | 1.1 MB | 144.6 MB |
| samples every 256 | 244.4 | 2.1 MB | 145.7 MB |
| shared + every 512 | 222.8 | 1.1 MB | 144.6 MB |
| shared + every 256 | 211.3 | 2.1 MB | 145.7 MB |
| shared + every 128 | 200.6 | 4.3 MB | 147.8 MB |
| **shared + every 64 (the default now)** | **193.6** | **8.5 MB** | **152.1 MB** |

The explicit RLCSA takes 70–73 ns per step in the same runs, so the gap per step is now 2.7×, down from 3.7×. `find --efcsa` gives `find --preset silva`'s output on V3_V4_1 and V1_V2_2 (`-T g.tag`), and the self-test passes. These have not yet been timed per read: that waits for the next rotation.

## Blocked interpolative coding, bits per value (2026-11-07)

`experiments/blocked_interp.cpp`, `blocked_interp.txt`: code lengths only, on the same lists as "Interpolative coding, bits per value". Each character's list is cut into blocks of B values. Each block's first value (its head) is in an Elias-Fano list over the list's universe. The other B − 1 values are coded by BIC between the block's head and the next head, in preorder (middle, left half, right half). "Pointers" are an Elias-Fano list of each block's first bit, needed to jump to a block. Per value:

| list | B | heads | BIC (centred) | pointers | total | MB | values decoded, access | values decoded, search |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| run starts | 64 | 0.221 | 6.000 (5.952) | 0.164 | 6.385 | 54.5 | 31.50 | 31.50 |
| F positions | 64 | 0.188 | 3.094 (3.142) | 0.149 | 3.432 | 29.3 | 31.50 | 34.51 |
| run starts | 128 | 0.118 | 6.096 (6.046) | 0.090 | 6.305 | 53.8 | 63.50 | 61.17 |
| F positions | 128 | 0.102 | 3.168 (3.215) | 0.083 | 3.353 | 28.6 | 63.50 | 67.06 |

- **Space:** blocks cost almost nothing. 6.39 and 3.43 bits per value at B = 64, against 6.21 and 3.26 for BIC on whole lists and 8.14 and 6.05 for Elias-Fano. Both lists together are 83.8 MB (B = 64) or 82.4 MB (B = 128), against 143.5 MB in Elias-Fano (152.1 MB with the new samples): 42% less.
- **Values decoded per lookup:** access(j), for j uniform over the values (a head costs nothing), decodes every value before j's in preorder. Skipping a left half means decoding it, since its code lengths depend on its values. That is (B − 1)/2 on average: 31.5 at B = 64, 63.5 at B = 128. A search for x (count_less, x uniform over the universe, inside the block the heads select) decodes until x's gap is known: 31.5–34.5 at B = 64, 61–67 at B = 128. A step needs two searches in run starts and two accesses in F positions, so about 130 decoded values per step at B = 64, where the Elias-Fano RLCSA now does 3 to 6 selects. Storing each left half's code length would cut this to about log₂ B per lookup, for more bits per value; that was not estimated.
