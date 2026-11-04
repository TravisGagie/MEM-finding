# KeBaB revision: Section 5.1 setting with ropebwt3 (Zymo positive reads)

Desktop (i7-8700, 14 GiB), 2026-11-03/04. Everything ran single-threaded, as the paper measured (appendix A), each
step in an 11 GB cgroup with no swap.

**Setting:**
- **Reads:** the 581,802 positive (microbial) Zymo reads (`~/zymo/all.fq`), 11,456,712,680 bases, with no
  homopolymer compression and no low-complexity filter.
- **KeBaB:** its own tool (https://github.com/drnatebrown/kebab at 4bed76f), with `kebab build -k 20 -e 0.1 -f 1`
  (one hash function) and `kebab scan -l 40 -s [--top-t N]`.
- **ropebwt3:** 3.10-r285, unmodified, `ropebwt3 mem -l 40`.
- **Scripts and logs:** `build.sh`, `build2.sh`, `run.sh`, `run2.sh`, `logs/`.

**Same setting as the paper:**
- Our pseudo-MEMs at L = 40 average 70.97 bases, exactly the paper's figure for positive reads.
- KeBaB's filter is 0.2684 GB, as in the paper.
- The paper gives the positive reads an average length of 25,378 and the yeast reads 19,693. Our positive reads
  average 19,692, so the paper seems to have swapped the two.

## Indexes

| | time | peak | size | paper (authors' server) |
|---|---|---|---|---|
| ropebwt3 index, both strands (`ropebwt3 build -t 12 -m 1g -do`, the 7 FASTAs in `list.tsv` order) | 8,707 s (2.42 h), 12 threads | 10.06 GB | 0.770 GB (`zymo.fmd`) | 162.9 min (2.71 h), 0.80 GB |
| KeBaB filter (`kebab build -k 20 -e 0.1 -f 1 -t 1`) | 2,451 s, one thread | 0.27 GB | 0.2684 GB | 4.02 min (16 threads), 0.2684 GB |

About the filter:
- **Estimates:** 308,721,725 distinct 20-mers estimated; observed false-positive rate 0.130.
- **Why one thread:** with more than one, `kebab`'s `process_sequences` leaks every sequence it reads (the copies
  are freed only once, after the loop). The 12-thread build was killed by the cap at 11.5 GB, 43% of the way in.
  The same happens in `kebab scan`. `kebab_leak_note.md` is a note for KeBaB's authors. With one thread there is
  no leak, and the filter is the same.
- **Uncompressed input:** `kebab build` reads only plain FASTA, so it ran on an uncompressed copy of the reference
  (27.2 GB) on the external drive.

## Runs

- **(a)** ropebwt3 on the t longest pseudo-MEMs only (`kebab scan --top-t t`).
- **(b)** ropebwt3 on the T = 2t longest (`--top-t T+1`; `split_frags.py` keeps the T longest and notes the
  (T+1)st's length λ). For a read with more than T pseudo-MEMs and fewer than t MEMs of length at least λ, rerun
  ropebwt3 on the whole read (`longest.py fallback`). A read with at most T pseudo-MEMs needs no check: every MEM
  of length ≥ L lies in a pseudo-MEM, and all of them were searched.
- **(c)** ropebwt3 on whole reads, the exact reference: every MEM of length ≥ 40.

For each run, the t longest MEMs of each read (ties at the t-th place included) are compared with (c)'s
(`longest.py compare`). MEMs from pseudo-MEMs are mapped back to read coordinates; the headers are
`read:START-END`, 1-based inclusive. On a test index, the MEMs of all pseudo-MEMs, mapped back, are exactly those
of the whole reads.

| run | t | time (scan + ropebwt3 [+ fallback]) | speed-up over (c) | read bases searched | fallbacks | reads whose t longest lengths differ from (c) | reads whose spans (with ties) differ |
|---|---|---|---|---|---|---|---|
| (c) whole reads | — | 7,124 s (622 ns per base) | 1 | 100% | — | reference | reference |
| (a) t longest pseudo-MEMs | 1 | 221 + 56 = 277 s | 25.7× | 1.03% | — | **30,235 (5.20%)** | 39,675 (6.82%) |
| (b) T = 2 | 1 | 221 + 102 + 39 = 362 s (+ 4 s split and check) | 19.7× | 1.90% | 2,413 (0.42%) | **0** | 1,038 (0.18%) |
| (a) | 5 | 224 + 222 = 446 s | 16.0× | 4.12% | — | **89,424 (15.37%)** | 129,141 (22.20%) |
| (b) T = 10 | 5 | 226 + 397 + 1 = 624 s (+ 9 s) | 11.4× | 7.24% | 20 (0.003%) | **0** | 16 (0.003%) |
| (a) | 10 | 225 + 397 = 621 s | 11.5× | 7.24% | — | **125,200 (21.52%)** | 190,291 (32.71%) |
| (b) T = 20 | 10 | 227 + 688 + 1 = 916 s (+ 15 s) | 7.8× | 12.40% | 3 (0.001%) | **0** | 0 |

How the table was read:
- **Peak memory:** ropebwt3 1.0–1.3 GB; kebab scan 0.27 GB.
- **Filtration:** all pseudo-MEMs at L = 40 cover 53.66% of the read bases (86,626,626 fragments, 6.15 Gbases).
  The table's percentages are those actually searched.
- **Reads with no MEM ≥ 40:** 6,585 in (c), of which 6,410 have no pseudo-MEM at all. (a) at t = 1 gets no MEM
  for 85 more: their longest pseudo-MEM contains no MEM ≥ 40. That is one more at t = 5, and none at t = 10.
- **Fallback reason:** every fallback in (b) is "fewer than t MEMs reach λ". By construction, a read with at most
  T pseudo-MEMs is accepted without a check: 3,747, 10,632 and 18,236 reads for t = 1, 5 and 10.
- **(b)'s spans:** the remaining span differences are ties. Its acceptance test is "at least as long as λ", so a
  MEM tied with the t-th longest at length exactly λ can lie in an unsearched pseudo-MEM. The lengths of the t
  longest are always right, and at t = 10 even the spans all agree.

## Reading

- **KeBaB's heuristic misses one of the t longest MEMs in 5.2% of reads at t = 1, 15.4% at t = 5 and 21.5% at
  t = 10.** On the same reads with the zymo preset (homopolymer compression, θ = 14, L = 20; below), it missed in
  12.4%, 26.6% and 35.3% (t = 1, 3, 5). Even without compression, at L = 40, the miss rate is far from negligible
  and grows with t.
- **The check rule (b) gets every read's t longest lengths right.** Its cost: it searches about twice as many
  pseudo-MEMs, 1.3–1.5 times the time of (a), with fallbacks in 0.42% of reads at t = 1 and almost none beyond.
  It is still 7.8–19.7 times faster than ropebwt3 on whole reads.

# The check rule under the zymo preset (from `../zymo_kebab/`)

Homopolymer-compressed reads, `-U -D 14 -c --csa -j 12`, KeBaB k = 20 with katka's own filter, L = 20. The exact
reference is `--longest t -L 1`. Full details are in `../zymo_kebab/README.md`.

| rule | t | wall (12 threads) | pseudo-MEMs searched / read | fallbacks | reads whose t longest lengths differ from exact | accuracy (longest MEMs vote, prop) |
|---|---|---|---|---|---|---|
| heuristic `--kebab-top t` | 1 | 89.3 s | 1.00 | — | 72,280 (12.42%) | 98.72% |
| check `--kebab-check 2 1 --strict` | 1 | 97.7 s | 1.15 | 17,213 (2.96%) | 0 | 99.61% |
| batch `--kebab-batch 1` | 1 | 91.1 s | 1.29 | 0 | 0 | 99.61% |
| heuristic | 3 | 100.5 s | 2.98 | — | 154,636 (26.58%) | 99.26% |
| check `--kebab-check 6 3 --strict` | 3 | 112.3 s | 3.41 | 10,729 (1.84%) | 0 | 99.61% |
| batch | 3 | 107.9 s | 3.72 | 0 | 0 | 99.61% |
| heuristic | 5 | 109.6 s | 4.97 | — | 205,491 (35.32%) | 99.40% |
| check `--kebab-check 10 5 --strict` | 5 | 125.4 s | 5.62 | 10,127 (1.74%) | 0 | 99.61% |
| batch | 5 | 112.9 s | 6.07 | 0 | 0 | 99.61% |
| exact, no filter (`--longest t -L 1`) | 1 / 3 / 5 | 126.8 / 186.3 / 227.4 s | — | — | reference | 99.61% |

- **Non-strict check:** without `--strict`, the check rule takes 103.3, 110.3 and 122.2 s, with 13,840, 9,026
  and 8,946 fallbacks. Its t longest lengths are always right; only ties beyond the t-th place differ, in 0.24%,
  0.10% and 0.08% of reads.
- **Fallback reasons:** every fallback was for "fewer than t long MEMs".
- **KeBaB's own pipeline at the journal's fixed L = 25:** `--kebab-top t --fixed-L`, every MEM voting, gives 98.66,
  99.23 and 99.40% for t = 1, 3 and 5. Every MEM ≥ 25 gives 99.58%, but takes 4.3–5.4 times longer.

# The check rule with smaller T (Section 5.1 setting, 2026-11-04)

Same setting and pipeline as (b) above, single-threaded (`run3.sh`, `run3b.sh`, `run3c.sh`; logs
`logs/run3*.log`). Time is scan + ropebwt3 on the T longest pseudo-MEMs + ropebwt3 on the fallback reads. Splitting
and checking add 2–9 s, not included. (a) is the heuristic: ropebwt3 on the t longest pseudo-MEMs only.

| t | T | time | × (a) | fallbacks | reads whose t longest lengths differ from (c) | reads whose spans (with ties) differ |
|---|---|---|---|---|---|---|
| 1 | (a) | 221 + 56 = 277 s | 1 | — | 30,235 (5.20%) | 39,675 (6.82%) |
| 1 | **2** (= t + 1 = ⌈1.25t⌉ = 2t) | 221 + 102 + 39 = 362 s | **1.31** | 2,413 (0.42%) | 0 | 1,038 (0.18%) |
| 1 | 3 (t + 2) | 227 + 149 + 4 = 380 s | 1.37 | 225 (0.04%) | 0 | 109 (0.02%) |
| 5 | (a) | 224 + 222 = 446 s | 1 | — | 89,424 (15.37%) | 129,141 (22.20%) |
| 5 | 6 (t + 1) | 226 + 273 + 236 = 735 s | 1.65 | 14,207 (2.47%) | 0 | 10,142 (1.74%) |
| 5 | **7** (t + 2 = ⌈1.25t⌉) | 226 + 302 + 43 = 571 s | **1.28** | 2,261 (0.39%) | 0 | 2,065 (0.36%) |
| 5 | 10 (2t) | 226 + 397 + 1 = 624 s | 1.40 | 20 (0.003%) | 0 | 16 (0.003%) |
| 10 | (a) | 225 + 397 = 621 s | 1 | — | 125,200 (21.52%) | 190,291 (32.71%) |
| 10 | 11 (t + 1) | 225 + 437 + 459 = 1,121 s | 1.81 | 27,160 (4.72%) | 0 | 25,655 (4.41%) |
| 10 | 12 (t + 2) | 228 + 467 + 105 = 800 s | 1.29 | 5,458 (0.95%) | 0 | 7,506 (1.29%) |
| 10 | **13** (⌈1.25t⌉) | 227 + 497 + 25 = 749 s | **1.21** | 1,177 (0.20%) | 0 | 1,870 (0.32%) |
| 10 | 15 (1.5t) | 227 + 558 + 2 = 787 s | 1.27 | 76 (0.01%) | 0 | 102 (0.02%) |
| 10 | 20 (2t) | 227 + 688 + 1 = 916 s | 1.48 | 3 (0.001%) | 0 | 0 |

How the table was read:
- **Fallback reasons:** all of them are "fewer than t MEMs reach the (T+1)st pseudo-MEM's length". Reads with at
  most T pseudo-MEMs need no check.
- **Exactness:** every check-rule run gets every read's t longest lengths right. The remaining span differences are
  MEMs tied at exactly λ, the length of the (T+1)st pseudo-MEM. A `--strict`-style test ("longer than λ") would
  catch them at the cost of more fallbacks.

**Reading:**
- **The cost curve is U-shaped in T.** With T = t + 1 the fallbacks dominate: 2.5–4.7% of reads go to whole-read
  ropebwt3, and each costs a full read. With T = 2t the extra pseudo-MEMs dominate.
- **T = ⌈1.25t⌉ is the best of those tried for every t:** 1.31, 1.28 and 1.21 times the heuristic, with 0.2–0.4%
  fallbacks.
- **The fixed cost:** the `kebab scan` step (about 225 s, single-threaded) is the same for every rule. Against
  ropebwt3's share alone, the check rule at T = ⌈1.25t⌉ costs 2.5×, 1.6× and 1.3× the heuristic's (141 / 56,
  345 / 222 and 522 / 397 s).

# KeBaB's pipeline at the journal's fixed L = 25, and the check rule, under the zymo preset

Homopolymer compression, θ = 14, KeBaB k = 20, RLCSA, 12 threads (`fixedL.sh`, `logs/fixedL.log`).
- **KeBaB's pipeline:** `-L 25 -K 20 --kebab-top t --fixed-L` reports every MEM of length ≥ 25 in the t longest
  pseudo-MEMs, and every MEM votes (`credit_schemes.py`).
- **The check rule:** `-L 20 -K 20 --kebab-check 2 1`, runs from 2026-11-03 (`../zymo_kebab/`), scored both ways.
  L counts compressed characters.

| run | wall | pseudo-MEMs searched / read | backward steps / read | MEMs | unclassified | accuracy, every MEM votes (prop / even) | accuracy, longest MEMs vote (prop / even) |
|---|---|---|---|---|---|---|---|
| KeBaB pipeline, t = 1 | 87.5 s | 1.00 | 322 | 758,336 | 1.05% | 98.66% / 98.37% | |
| KeBaB pipeline, t = 3 | 98.0 s | 3.00 | 836 | 2,195,477 | 0.46% | 99.23% / 99.19% | |
| KeBaB pipeline, t = 5 | 109.7 s | 5.00 | 1,281 | 3,587,812 | 0.28% | 99.40% / 99.40% | |
| KeBaB pipeline, t = 10 | 137.9 s | 10.00 | 2,248 | 6,975,234 | 0.09% | 99.53% / 99.59% | |
| KeBaB pipeline, t = 20 | 169.5 s | 20.00 | 3,864 | 13,529,233 | 0.02% | 99.60% / 99.68% | |
| KeBaB pipeline, t = 30 | 199.3 s | 29.97 | 5,244 | 19,896,633 | 0.01% | 99.61% / 99.70% | |
| all MEMs ≥ 25 (`-L 25 -K 20`) | 468.9 s | all | — | 108,366,876 | 0.003% | 99.58% / 99.69% | |
| **check rule, t = 1, T = 2** | 103.3 s | 1.15 | 743 | 591,459 | 0 | **99.63%** / 99.29% | 99.61% / 99.27% |
| check rule, t = 1, T = 2, `--strict` | 97.7 s | 1.15 | 794 | 592,877 | 0 | 99.63% / 99.29% | 99.61% / 99.27% |

**Reading:**
- **KeBaB's pipeline:** it improves with t and reaches the all-MEMs accuracy at t = 20–30, as the journal reports
  with Tagger. That takes 1.6–1.9 times the check rule's time.
- **The check rule with t = 1:** the longest MEMs alone, about 1 per read, give 99.63% with proportional credit.
  That is above the pipeline at t = 30 (99.61%) and above every MEM ≥ 25 (99.58%), in about half the time of t = 30
  and less than a quarter of all MEMs'.
- **Even credit:** it favours the pipeline (99.70% at t = 30, against 99.29%). With only the longest MEMs,
  proportional and even credit differ more.
