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
