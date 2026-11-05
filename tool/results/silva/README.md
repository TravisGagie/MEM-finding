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
