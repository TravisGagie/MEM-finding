#!/usr/bin/env python3
# maxmem_roc.py listings truth [tbl]
# For rz-classify -X output (adaptive BML: the longest MEMs of each read are among those reported), score each read
# by the length in bases of its longest MEM and print, for each threshold t, the sensitivity (non-yeast reads with a
# MEM of at least t bases) and the specificity (yeast reads without one), as in SPUMONI 2's Table 2.
# truth: read name and label per line ("yeast" for the negatives).  With tbl (index, tag, ...), also prints the
# species accuracy of the accepted non-yeast reads when only their longest MEMs vote, with proportional credit.
import re, sys, collections
out, truthf = sys.argv[1:3]; tbl = sys.argv[3] if len(sys.argv) > 3 else None
names = []
if tbl:
    for l in open(tbl):
        s = l.split('\t')[1]
        if s not in names: names.append(s)
truth = dict(l.split()[:2] for l in open(truthf))
pat = re.compile(r'\[(\d+),(\d+)\] \{([^}]*)\}')
score = {}; vote = {}; vote_even = {}; name = None
for line in open(out):
    if line[0] == '>': name = line[1:].strip(); score[name] = 0; continue
    mems = [(int(b) - int(a) + 1, d) for a, b, d in pat.findall(line)]
    if not mems: continue
    M = max(m for m, _ in mems); score[name] = M
    w = collections.Counter(); we = collections.Counter()
    for m, d in mems:
        if m < M: continue
        cs = [(int(x.split(':')[0]), int(x.split(':')[1]) if ':' in x else 1) for x in d.split(',') if x]
        tot = sum(c for _, c in cs)
        for t, c in cs: w[t] += m * c / tot; we[t] += m / len(cs)
    for W_, V_ in ((w, vote), (we, vote_even)):
        if W_:
            b = max(W_.values()); top = [t for t in W_ if W_[t] == b]
            V_[name] = top[0] if len(top) == 1 else None
pos = [k for k in score if truth[k] != 'yeast']; neg = [k for k in score if truth[k] == 'yeast']
print('reads %d (%d positive, %d yeast); longest MEM, median: positive %d, yeast %d bases'
      % (len(score), len(pos), len(neg), sorted(score[k] for k in pos)[len(pos) // 2] if pos else 0,
         sorted(score[k] for k in neg)[len(neg) // 2] if neg else 0))
print('    t   sens    spec' + ('   species acc (prop / even)' if names else ''))
for t in range(15, 81, 1 if len(sys.argv) > 4 else 5):
    sens = 100 * sum(score[k] >= t for k in pos) / max(1, len(pos))
    spec = 100 * sum(score[k] < t for k in neg) / max(1, len(neg))
    row = '  %3d  %6.2f  %6.2f' % (t, sens, spec)
    if names:
        for V_ in (vote, vote_even):
            ok = sum(score[k] >= t and V_.get(k) is not None and names[V_[k]] == truth[k] for k in pos)
            row += '   %6.2f' % (100 * ok / max(1, len(pos)))
    print(row)
