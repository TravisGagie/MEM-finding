#!/usr/bin/env python3
# maxmem_roc_stream.py [--plain | --info PREFIX.info] listings truth tbl reads.fq
# katka-readable/scripts/maxmem_roc.py with the reads given (-U: each MEM scored by its homopolymer-compressed
# length), streaming instead of loading every read: the listings must be in the reads' order (katka find and
# katka-classify write reads back in input order), which is checked name by name. Same output.
# --plain (or --info PREFIX.info, if its hpc line is 0): lengths in bases, end - start + 1, for an index built
# without -U. The default (no option) is the homopolymer-compressed length, for an index built with -U.
# With --plain the output is that of maxmem_roc.py without the reads (thresholds in steps of 5, not 2).
import re, sys, collections
def length_mode(argv):
    # --plain, or --info PREFIX.info with hpc 0 (an index built without -U): MEM length = end - start + 1
    plain, rest, i = False, [], 0
    while i < len(argv):
        if argv[i] == '--plain': plain = True
        elif argv[i] == '--info':
            l = open(argv[i + 1]).readline().split(); i += 1
            plain = plain or not (len(l) == 2 and l[0] == 'hpc' and int(l[1]))
        else: rest.append(argv[i])
        i += 1
    return plain, rest
plain, sys.argv[1:] = length_mode(sys.argv[1:])
out, truthf, tbl, fq = sys.argv[1:5]
names = []
for l in open(tbl):
    s = l.split('\t')[1]
    if s not in names: names.append(s)
truth = dict(l.split()[:2] for l in open(truthf))
pat = re.compile(r'\[(\d+),(\d+)\] \{([^}]*)\}')
F = open(fq)
def next_read():
    h = F.readline()
    if not h: sys.exit('reads ended before the listings')
    s = F.readline().strip(); F.readline(); F.readline()
    return h[1:].split()[0], s
def clen(r, a, b):
    if plain: return b - a + 1
    r = r[a:b + 1]
    return sum(1 for i in range(len(r)) if i == 0 or r[i] != r[i - 1])
score = {}; vote = {}; vote_even = {}
def finish(name, seq, lines):
    score[name] = 0
    mems = []
    for line in lines: mems += [(clen(seq, int(a), int(b)), d) for a, b, d in pat.findall(line)]
    if not mems: return
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
name = None; seq = None; lines = []
for line in open(out):
    if line[0] == '>':
        if name is not None: finish(name, seq, lines)
        name = line[1:].strip(); rn, seq = next_read(); lines = []
        if rn != name: sys.exit('listings and reads out of step: %s vs %s' % (name, rn))
        continue
    lines.append(line)
if name is not None: finish(name, seq, lines)
pos = [k for k in score if truth[k] != 'yeast']; neg = [k for k in score if truth[k] == 'yeast']
print('reads %d (%d positive, %d yeast); longest MEM, median: positive %d, yeast %d bases'
      % (len(score), len(pos), len(neg), sorted(score[k] for k in pos)[len(pos) // 2] if pos else 0,
         sorted(score[k] for k in neg)[len(neg) // 2] if neg else 0))
print('    t   sens    spec   species acc (prop / even)')
for t in range(10, 61, 5 if plain else 2):   # the steps of maxmem_roc.py without and with the reads
    sens = 100 * sum(score[k] >= t for k in pos) / max(1, len(pos))
    spec = 100 * sum(score[k] < t for k in neg) / max(1, len(neg))
    row = '  %3d  %6.2f  %6.2f' % (t, sens, spec)
    for V_ in (vote, vote_even):
        ok = sum(score[k] >= t and V_.get(k) is not None and names[V_[k]] == truth[k] for k in pos)
        row += '   %6.2f' % (100 * ok / max(1, len(pos)))
    print(row)
