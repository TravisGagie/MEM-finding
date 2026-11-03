#!/usr/bin/env python3
# compare_longest_t.py run.txt reference.txt reads.fq t
# Compares two katka listings, read by read, made with --longest t, --kebab-top t, --kebab-batch t or -K
# --longest t (each reports every MEM at least as long as the t-th longest it found, ties included; lengths
# are homopolymer-compressed with -U, as the search sees them).  The reference is the exact run without a
# filter (--longest t -L 1).  For each read, the t longest of each run (with ties at the t-th) are compared:
#   lengths differ:  the sorted compressed lengths of the t longest differ (a longest MEM was missed, or a
#                    shorter one reported in its place);
#   sets differ:     the spans of the t longest differ (includes the above);
#   reference < k:   the reference's longest MEM is shorter than k (default 20), so no pseudo-MEM holds it.
# Streams the three files in read order.
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
def t_longest(seq, lines, t):
    m = sorted({(clen(seq, int(a), int(b)), int(a), int(b)) for l in lines for a, b, d in pat.findall(l)}, reverse=True)
    if len(m) > t: m = [x for x in m if x[0] >= m[t - 1][0]]
    return m
if len(sys.argv) < 5: sys.exit('usage: compare_longest_t.py run.txt reference.txt reads.fq t [k]')
t = int(sys.argv[4]); k = int(sys.argv[5]) if len(sys.argv) > 5 else 20
n = dl = ds = short = dl_short = 0
for (na, la), (nb, lb), (nr, s) in zip(records(sys.argv[1]), records(sys.argv[2]), reads(sys.argv[3])):
    if not (na == nb == nr): sys.exit('out of step at %s %s %s' % (na, nb, nr))
    n += 1
    A, B = t_longest(s, la, t), t_longest(s, lb, t)
    lens_differ = [x[0] for x in A] != [x[0] for x in B]
    ref_short = bool(B) and B[0][0] < k
    dl += lens_differ; ds += A != B; short += ref_short; dl_short += lens_differ and ref_short
print('reads %d, t = %d: lengths of the t longest differ in %d (%.4f%%), of which %d have a reference longest MEM shorter than %d; '
      'spans differ in %d (%.4f%%); reference longest shorter than %d: %d reads' % (n, t, dl, 100.0 * dl / max(1, n), dl_short, k, ds, 100.0 * ds / max(1, n), k, short))
