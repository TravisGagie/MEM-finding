#!/usr/bin/env python3
# step_bounds.py mems.tsv reads.fq L t total_char_steps [split=50]
# Bounds on BML's backward steps spent crossing each reported MEM R[a..e) (katka find --format mems: name, a, e
# exclusive, sp, ep), from the structure of bml_segment (tool/src/mems.hpp): the successful window search
# back(0, x, x+L) takes L - t steps after the t-mer lookup table, extend_left takes x - a steps plus one failing
# step if a > 0, and forward(a) takes e - a - t steps plus one failing step if e < n. The window start x is not in
# the output, only a <= x <= e - L, so a MEM costs between L + len - 2t and 2 len - 2t steps (plus the failing ones).
# The other steps (total from --stats minus these) are in failed windows.
import sys, collections
mems, fq, L, t, total = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
SPLIT = int(sys.argv[6]) if len(sys.argv) > 6 else 50     # long MEMs: longer than SPLIT
n = {}; R = 0                       # read names repeat across regions, so reads are counted, not names
with open(fq) as F:
    for i, line in enumerate(F):
        if i % 4 == 0: name = line[1:].split()[0]; R += 1
        elif i % 4 == 1:
            if n.get(name, len(line.strip())) != len(line.strip()): sys.exit('reads with the same name and different lengths')
            n[name] = len(line.strip())
lo = collections.Counter(); hi = collections.Counter(); cnt = collections.Counter(); lens = collections.Counter()
for line in open(mems):
    name, a, e = line.split('\t')[:3]; a, e = int(a), int(e); ln = e - a
    k = 'long (> %d)' % SPLIT if ln > SPLIT else 'short (30-%d)' % SPLIT
    fail = (a > 0) + (e < n[name])
    lo[k] += L + ln - 2 * t + fail; hi[k] += 2 * ln - 2 * t + fail; cnt[k] += 1; lens[k] += ln
print('reads %d, total backward steps %d (%.1f per read)' % (R, total, total / R))
for k in ('long (> %d)' % SPLIT, 'short (30-%d)' % SPLIT):
    print('  MEMs %-14s %6d (%.2f per read, mean length %.1f): crossing steps %d to %d (%.1f to %.1f per read, %.1f%% to %.1f%% of all)'
          % (k, cnt[k], cnt[k] / R, lens[k] / max(1, cnt[k]), lo[k], hi[k], lo[k] / R, hi[k] / R, 100 * lo[k] / total, 100 * hi[k] / total))
LO, HI = sum(lo.values()), min(total, sum(hi.values()))   # crossings cannot exceed the total
print('  failed windows (the rest): %d to %d (%.1f to %.1f per read, %.1f%% to %.1f%% of all)' % (total - HI, total - LO, (total - HI) / R, (total - LO) / R, 100 * (total - HI) / total, 100 * (total - LO) / total))
