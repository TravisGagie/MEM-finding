# KeBaB stopping rules on Zymo (for the KeBaB journal revision)

Desktop (i7-8700, 14 GiB), 2026-11-03. The tool was at 7500ed0 (`tool/NOTES.md`, "KeBaB stopping rules for the
journal revision"). Data: the full Zymo reads (581,802; `~/zymo/all.fq`) on the new Zymo index (`~/zymo/hpc/zy`).
The common options were `-U -D 14 -c --csa -j 12 --stats`, the `zymo` preset settings, with KeBaB k = 20. The runs
went one at a time, each in a cgroup capped at 11 GB with no swap (`run.sh`; full log in `run.log`, stderr of each
run in `logs/`). The listings and the `--read-stats` files were compressed with zstd and are not committed.

Runs, for t = 1, 3, 5:
- (a) heuristic: `-L 20 -K 20 --kebab-top t`;
- (b) batch rule: `-L 20 -K 20 --kebab-batch t`;
- (b′) batch rule at L = 1: `-L 1 -K 20 --kebab-batch t`;
- (c) streaming rule: `-L 20 -K 20 --longest t`;
- (d) exact reference: `-L 1 --longest t`, no filter.

## Results

| run | t | wall | user CPU | peak | pseudo-MEMs searched / read | backward steps / read | accuracy (prop / even) | reads with a different t longest |
|---|---|---|---|---|---|---|---|---|
| (a) heuristic | 1 | 89.3 s | 901 s | 7.34 GB | 1.00 | 318 | 98.72% / 98.37% | **72,280 (12.42%)** |
| (b) batch | 1 | 91.1 s | 919 s | 7.35 GB | 1.29 | 388 | 99.61% / 99.27% | 0 |
| (b′) batch, L = 1 | 1 | 92.0 s | 926 s | 7.35 GB | 1.29 | 388 | 99.61% / 99.27% | 0 (against d) |
| (c) streaming | 1 | 91.1 s | 921 s | 7.34 GB | 1.29 | 388 | 99.61% / 99.27% | 0 (against d) |
| (d) exact, no filter | 1 | 126.8 s | 1,353 s | 7.20 GB | — | — | 99.61% / 99.27% | reference |
| (a) heuristic | 3 | 100.5 s | 1,040 s | 7.34 GB | 2.98 | 876 | 99.26% / 98.92% | **154,636 (26.58%)** |
| (b) batch | 3 | 107.9 s | 1,084 s | 7.35 GB | 3.72 | 1,016 | 99.61% / 99.27% | 0 |
| (b′) batch, L = 1 | 3 | 110.6 s | 1,078 s | 7.34 GB | 3.72 | 1,016 | 99.61% / 99.27% | 0 (against d) |
| (c) streaming | 3 | 102.5 s | 1,065 s | 7.34 GB | 3.72 | 1,016 | 99.61% / 99.27% | 0 (against d) |
| (d) exact, no filter | 3 | 186.3 s | 2,045 s | 7.20 GB | — | — | 99.61% / 99.27% | reference |
| (a) heuristic | 5 | 109.6 s | 1,147 s | 7.34 GB | 4.97 | 1,355 | 99.40% / 99.07% | **205,491 (35.32%)** |
| (b) batch | 5 | 112.9 s | 1,187 s | 7.34 GB | 6.07 | 1,539 | 99.61% / 99.27% | 0 |
| (b′) batch, L = 1 | 5 | 113.4 s | 1,189 s | 7.35 GB | 6.07 | 1,539 | 99.61% / 99.27% | 0 (against d) |
| (c) streaming | 5 | 114.6 s | 1,189 s | 7.34 GB | 6.07 | 1,539 | 99.61% / 99.27% | 0 (against d) |
| (d) exact, no filter | 5 | 227.4 s | 2,511 s | 7.20 GB | — | — | 99.61% / 99.27% | reference |

How the columns were measured:
- **Wall time** includes loading the index, about 9 s.
- **Accuracy** comes from `scripts/maxmem_roc_stream.py`, the streamed `maxmem_roc.py`, at the smallest threshold.
  Only the longest MEMs vote, measured in compressed characters; for every t, accuracy counts as "the longest
  MEMs voting".
- **The last column** comes from `scripts/compare_longest_t.py`. For (a) and (b) it is against (c), the exact
  rule at the same `-L 20`; for (b′) and (c) it is against (d). For (a), the comparison against (d) gives the same
  counts. In no read is the reference longest MEM shorter than 20.

Batches:
- **Second batches:** 17,213 reads (t = 1), 10,729 (t = 3) and 10,127 (t = 5) needed one. That is 2.96%, 1.84%
  and 1.74% of the reads.
- **Whole-read fallback:** it never ran, not even at L = 1 (b′). Every read's pseudo-MEMs held at least t MEMs
  of length ≥ 20.

In all, 600,809,096 pseudo-MEMs, 1,032.7 per read; the searched counts above are out of these.

## Reading

- **KeBaB's heuristic often misses a longest MEM on real reads:** in 12.4% of reads for t = 1, rising to 35.3%
  for t = 5. It also costs accuracy: 98.72% against 99.61% at t = 1, and still 99.40% at t = 5. For t = 1,
  0.81% of reads get no MEM at all (sensitivity 99.19% at the lowest threshold). For them the longest pseudo-MEM
  contains no MEM of length ≥ 20.
- **The batch rule matches the exact answer for every read and every t, at nearly the heuristic's speed:**
  2–8% more wall time, searching 1.29 pseudo-MEMs per read instead of 1 at t = 1. So the case for it is accuracy
  as well as the guarantee.
- **With L = 20 and k = 20, the batch rule and the streaming rule search exactly the same pseudo-MEMs.** Their
  counts and backward steps are identical, so their times differ only by noise. NOTES predicted this: a
  pseudo-MEM shorter than the current L is skipped, so the batch rule differs only in when it checks.
- **The exact reference without the filter (d)** takes 1.4 to 2.0 times longer than the KeBaB rules.
