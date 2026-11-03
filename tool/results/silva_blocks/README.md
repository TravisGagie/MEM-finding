# Order-k character-level blocks on SILVA: feasibility (`katka-blocks`, commit b87e722)

Desktop (i7-8700, 14 GiB), 2026-11-02. The SILVA default index is `bac`, from `../silva_profile/`. The reads are
the same 8,000 as in the profile: the first 1,000 of each mate of each region, all 250 bases. Their names were
made unique by prefixing the region and mate. The original names repeat across regions, and `katka-blocks` matches
listing records to reads by name, so part (c) of the first run mixed MEMs of different reads (137,331 MEMs instead
of 43,414). The listing is `katka find -j 1 --preset silva` (L = 30, `-c --csa -F 10`) on those reads: 43,414
MEMs, and `katka-blocks` confirms 0 reads whose MEMs differ from its own BML.
Commands are in `run2.sh`; full output is in `blocks_k8_12.txt` and `blocks_k16.txt`.

Run time 48–49 s per invocation. Peak memory 1.87 GB for k = 8 and 12, and 2.74 GB for k = 16, whose hash table
took 1.61 GB.

## a–d, per k (from `katka-blocks`)

| | k = 8 | k = 12 | k = 16 |
|---|---|---|---|
| a. distinct k-mers followed by a base | 65,536 | 13,861,591 | 53,297,379 |
| b. blocks in S, mean length | 318.1 M, 3.69 | 126.3 M, 9.29 | 85.7 M, 13.68 |
| b. characters in blocks ≥ 17 | 32.3% | 61.0% | 75.4% |
| c. steps per character to cross the MEMs (plain: 1) | 0.616 | 0.693 | 0.775 |
| c. MEMs 33–64 / 65–128 / 129–256 | 0.693 / 0.531 / 0.464 | 0.826 / 0.579 / 0.384 | 0.913 / 0.677 / 0.425 |
| d. BML steps in failed windows | 24.3% | 24.3% | 24.3% |
| d. at the ends of searches (still character steps) | 33.3% | 53.9% | 61.7% |
| d. inside long crossings → block steps | 42.4% → 948,659 | 21.7% → 266,888 | 14.0% → 129,623 |
| d. BML total under the model (no lookup table) | 1.33× fewer | 1.20× fewer | 1.13× fewer |

`katka-blocks`' BML has no lookup table, so it takes 5,417,049 steps (677 per read). `katka find --preset silva`
takes 3,889,139 (486 per read). The difference is exactly the table's 10 steps for each of the 152,791 searches.

## Backward steps per read with blocks, against 486 now (`estimate.py`, `estimate.txt`)

| k | A: from (d), with the table | B: from the profile's bounds × (c) | now |
|---|---|---|---|
| 8 | 318 (1.53× fewer) | 278 – 343 (1.42 – 1.75×) | 486 |
| 12 | 372 (1.31×) | 317 – 374 (1.30 – 1.53×) | 486 |
| 16 | 408 (1.19×) | 361 – 405 (1.20 – 1.35×) | 486 |

Method A: take 486, subtract the steps (d) puts inside long crossings, and add their block steps. This assumes the
10 steps the table saves per search lie at the search's right end, which the model charges as character steps
anyway, and not inside a crossing.

Method B: the crossing bounds of `../silva_profile/step_bounds.py`, split at 64 to match (c)'s length classes
(`step_bounds64.txt`: MEMs ≤ 64 take 212–261 steps per read, MEMs ≥ 65 take 167–271). Each class is scaled by its
steps per character under the model, character-weighted: 0.707 / 0.514 for k = 8, 0.839 / 0.529 for k = 12 and
0.922 / 0.613 for k = 16. Failed windows are left unchanged.

The two agree. A lies inside B for k = 8 and 12, and just above it for k = 16.

## Reading

- **Step counts:** k = 8 saves the most steps, about a third (486 → about 318).
  - **Larger k:** blocks are longer, but a 250-base read's MEMs (median 47, 90th percentile 100) contain fewer
    complete blocks whose contexts lie inside the match.
  - **Unchanged:** the failed windows (24% of BML's steps) and the ends of searches don't change. The ends grow
    with k, from 33% to 62%.
- **Time:** these are step counts. If a block step cost the same as a character step, k = 8 would cut the search
  time (89.5% of 55.1 µs per read, `../silva_profile/`) by about a third, to about 39 µs per read. A block step
  probably costs more, since it is a step on another BWT, plus the continuation lookup. So that figure is an upper
  bound on the gain, not an estimate.
- **Continuation table:** k = 8 has 65,536 entries (2 bits each in an array: 16 KB); k = 12 has 13.9 M of the
  16.8 M possible (4 MB as an array); k = 16 has 53.3 M (a hash table or a 1 GB array).

Note for `katka-blocks`: it matches listing records to reads by name. Reads with repeated names, as when
regions are concatenated, silently corrupt part (c). Matching by order, or a warning on duplicate names, would
prevent this.
