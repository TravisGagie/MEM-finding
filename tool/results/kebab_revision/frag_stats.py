#!/usr/bin/env python3
# frag_stats.py reads.fq < fragments.fa: fragments, their total length, and the filtration ratio (fragment bases
# over read bases)
import sys
tot = 0
with open(sys.argv[1]) as F:
    for i, line in enumerate(F):
        if i % 4 == 1: tot += len(line.rstrip('\n'))
n = b = 0
for line in sys.stdin:
    if line[0] == '>': n += 1
    else: b += len(line.rstrip('\n'))
print('read bases %d; fragments %d, %d bases: filtration ratio %.4f (%.2f%% of the read bases kept)' % (tot, n, b, b / tot, 100 * b / tot))
