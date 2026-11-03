# Where SILVA classification spends its time (preset `silva`, L = 30, one thread)

Desktop (i7-8700, 14 GiB), 2026-11-02. The tool was at commit b87e722, unchanged. The index is the new SILVA
default, built from a `prep` and `bwt` whose output is byte-identical to KATKA's: `index --preset silva`,
348 s, peak 6.5 GB (`index_build.log`). It has an RLBWT of 155 MB with r = 68,324,430, an explicit RLCSA of 565 MB,
and a tag array of 836 MB with 359,351,512 runs, plus a 314 MB RMQ.

Reads: the first 1,000 reads of each mate of each region (V1_V2, V3_V4, V4_V4, V4_V5) from the paper's SILVA
reads (`katka/silva/work/compare/reads/`): 8,000 reads of 250 bases. `find -j 1 --preset silva`, that is
`-L 30 -c --csa -F 10`. Each variant ran 5 times; times are the tool's own wall time, which excludes loading the
index (`profile.sh`, `profile.log`).

`perf` could not be used: `kernel.perf_event_paranoid` is 4 on the desktop, and lowering it needs root.
valgrind is not installed. So the time is split by difference between runs, and the backward steps by bounds
derived from the reported MEMs. The exact split needs counters the tool does not have (see the end).

## Time per read

| variant | options | µs/read (median of 5) |
|---|---|---|
| full | `--preset silva` | 55.1 |
| without tags | `-L 30 --csa -F 10 --no-tags` | 50.4 |
| without search | `-L 1000 -c --csa -F 10` (no read is that long) | 1.1 |

| share | µs/read | % of 55.1 |
|---|---|---|
| (4) reading input, parsing, writing record headers | 1.1 | 2% |
| (3) tag-run scanning, credit (counts) and writing the tag lists | 4.7 | 8.5% |
| (1) + (2) MEM finding: backward steps (RLCSA), lookup table, control | 49.3 | 89.5% |

## Backward steps per read

From `--stats`: 3,889,139 backward steps over 8,000 reads, that is **486 per read**, in 19.1 searches per read,
with 5.4 MEMs per read. The lookup table (t = 10) saves the first 10 steps of each search. MEMs are listed only
from L = 30; their median length is 47 and the 90th percentile 100.

Each reported MEM R[a..e) costs, in `bml_segment` (`src/mems.hpp`), its successful window search
(L − 10 steps) plus `extend_left` (x − a steps) plus `forward` (len − 10 steps), plus one failing step at each end.
The window start x is not in the output, only a ≤ x ≤ e − L. So each MEM costs between L + len − 20 and
2 len − 20 steps, and the failed windows get the rest (`step_bounds.py`, `step_bounds.txt`):

| category | per read | % of steps | % of time (if every step costs the same) |
|---|---|---|---|
| (1) crossing MEMs longer than 50 (2.40 per read, mean length 83) | 228 – 356 | 47 – 73% | 42 – 66% |
|     crossing MEMs of 30 to 50 (3.03 per read, mean length 38) | 151 – 176 | 31 – 36% | 28 – 32% |
| (2) failed windows | 0 – 107 | 0 – 22% | 0 – 20% |
| (2) ends of MEMs (the failing step at each end, included above) | ≤ 10.9 | ≤ 2.2% | ≤ 2% |

The failed-window share is bounded by subtraction, so its lower bound is 0. Equal cost per step is an
assumption: steps on wide intervals at the start of a search may cost more or less than steps deep inside a
crossing.

## Answer: could faster crossing matter on SILVA?

Yes, it could. At least 47% of the backward steps, and roughly 42% or more of the time, are spent crossing MEMs
longer than 50 characters. Another 31–36% goes on crossing MEMs of 30–50 characters. Failed windows take at most
22% of the steps. Reading input and output is about 2% of the time, and tags and credit about 8.5%. MEMs overlap
heavily: 5.4 per 250-base read, so a read is crossed about 1.6 times over.

## Counters needed for the exact split

To be added by the cloud session; not added here:
- In `bml_segment`, take the `cnt.char_steps` delta of each call and add it to one of:
  - failed window: the `back()` call returns j > x;
  - successful window: the `back()` call with j == x;
  - `extend_left`;
  - `forward`.
  Bucket the last three by the reported MEM's length (e − a: 30–50, 51–100, > 100), and count separately the
  final failing step of `extend_left` and `forward`, which are the MEM ends.
- Optionally, a cycle counter (`__rdtsc()`) around the same four calls, around `found()` (tag scanning and credit)
  and around parsing and output, summed per thread. This gives time shares without assuming equal step costs.
  The overhead is about 2 rdtsc per call, against about 25 backward steps per call.
- Report them under `--stats`, per read.

With `perf` (after `sudo sysctl kernel.perf_event_paranoid=1`), `perf record -g` would split tags, I/O and search
by function. It could not split steps within crossings from steps in failed windows, since both are the same
inlined `step()`.
