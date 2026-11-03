#!/usr/bin/env python3
# topt_lengths.py run ref reads.fq t: reads whose t longest compressed lengths (the first t of the sorted lengths,
# without the extra ties) differ between two listings, and reads where only the tie set differs.
import re, sys
pat = re.compile(r'\[(\d+),(\d+)\]')
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
        s = F.readline().strip(); F.readline(); F.readline(); yield s
def clen(r, a, b):
    r = r[a:b + 1]; return sum(1 for i in range(len(r)) if i == 0 or r[i] != r[i - 1])
t = int(sys.argv[4]); top_diff = tie_only = n = 0
for (na, la), (nb, lb), s in zip(records(sys.argv[1]), records(sys.argv[2]), reads(sys.argv[3])):
    assert na == nb; n += 1
    A = sorted((clen(s, int(a), int(b)) for l in la for a, b in pat.findall(l)), reverse=True)
    B = sorted((clen(s, int(a), int(b)) for l in lb for a, b in pat.findall(l)), reverse=True)
    if A[:t] != B[:t]: top_diff += 1
    elif A != B: tie_only += 1
print('%d reads: the t = %d longest lengths differ in %d; only ties beyond the t-th differ in %d' % (n, t, top_diff, tie_only))
