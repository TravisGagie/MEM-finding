#!/usr/bin/env python3
# same_longest_mems.py listingsA listingsB reads.fq
# Checks that two -X runs (e.g. with and without KeBaB) find, for every read, the same longest MEMs: the same
# spans with the same {tag:count}, among the MEMs of maximum homopolymer-compressed length. Streams all three
# files in read order.
import re, sys
pat = re.compile(r'\[(\d+),(\d+)\] \{([^}]*)\}')
def records(f):
    name, lines = None, []
    for line in open(f):
        if line[0] == '>':
            if name is not None: yield name, lines
            name, lines = line[1:].strip(), []
        else: lines.append(line)
    if name is not None: yield name, lines
def reads(f):
    F = open(f)
    while True:
        h = F.readline()
        if not h: return
        s = F.readline().strip(); F.readline(); F.readline(); yield h[1:].split()[0], s
def clen(r, a, b):
    r = r[a:b + 1]; return sum(1 for i in range(len(r)) if i == 0 or r[i] != r[i - 1])
def longest(seq, lines):
    m = [(clen(seq, int(a), int(b)), int(a), int(b), d) for l in lines for a, b, d in pat.findall(l)]
    if not m: return set()
    M = max(x[0] for x in m); return {x for x in m if x[0] == M}
n = bad = 0
for (na, la), (nb, lb), (nr, s) in zip(records(sys.argv[1]), records(sys.argv[2]), reads(sys.argv[3])):
    if not (na == nb == nr): sys.exit('out of step at %s %s %s' % (na, nb, nr))
    n += 1
    if longest(s, la) != longest(s, lb):
        bad += 1
        if bad <= 5: print('  differs:', na, sorted(longest(s, la))[:3], sorted(longest(s, lb))[:3])
print('%d reads: %d with different longest MEMs' % (n, bad))
