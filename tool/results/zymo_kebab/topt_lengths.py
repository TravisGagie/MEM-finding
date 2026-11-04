#!/usr/bin/env python3
# topt_lengths.py [--plain | --info PREFIX.info] run ref reads.fq t: reads whose t longest compressed lengths
# (the first t of the sorted lengths, without the extra ties) differ between two listings, and reads where only
# the tie set differs.
# --plain (or --info PREFIX.info, if its hpc line is 0): lengths in bases, end - start + 1, for an index built
# without -U. The default (no option) is the homopolymer-compressed length, for an index built with -U.
import re, sys
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
    if plain: return b - a + 1
    r = r[a:b + 1]; return sum(1 for i in range(len(r)) if i == 0 or r[i] != r[i - 1])
t = int(sys.argv[4]); top_diff = tie_only = n = 0
for (na, la), (nb, lb), s in zip(records(sys.argv[1]), records(sys.argv[2]), reads(sys.argv[3])):
    assert na == nb; n += 1
    A = sorted((clen(s, int(a), int(b)) for l in la for a, b in pat.findall(l)), reverse=True)
    B = sorted((clen(s, int(a), int(b)) for l in lb for a, b in pat.findall(l)), reverse=True)
    if A[:t] != B[:t]: top_diff += 1
    elif A != B: tie_only += 1
print('%d reads: the t = %d longest lengths differ in %d; only ties beyond the t-th differ in %d' % (n, t, top_diff, tie_only))
