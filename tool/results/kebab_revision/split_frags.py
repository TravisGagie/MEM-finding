#!/usr/bin/env python3
# split_frags.py frags.fa T search.fa lambda.tsv: from kebab scan -s --top-t T+1 output (fragments of each read,
# longest first, headers read:START-END, 1-based inclusive), write the T longest of each read to search.fa and,
# per read, the number of fragments and the length of the (T+1)st (0 if there are at most T) to lambda.tsv.
import sys
src, T, out, lam = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
O = open(out, 'w'); Lm = open(lam, 'w')
cur, cnt, lt = None, 0, 0
def flush():
    if cur is not None: Lm.write('%s\t%d\t%d\n' % (cur, cnt, lt))
with open(src) as F:
    for h in F:
        s = F.readline()
        name, rng = h[1:].rstrip('\n').rsplit(':', 1)
        if name != cur: flush(); cur, cnt, lt = name, 0, 0
        cnt += 1
        if cnt <= T: O.write(h); O.write(s)
        elif cnt == T + 1: lt = len(s.rstrip('\n'))
flush()
