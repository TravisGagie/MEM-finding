#!/usr/bin/env python3
"""Which alpha a read set's genus composition corresponds to: k_g against N_g (tool/NOTES.md, "SILVA: credit
schemes when the sample doesn't match the database").

    silva_alpha_fit.py DOC_TO_TRAVERSAL NSEQ SEQTAX [SEQTAX ...]

DOC_TO_TRAVERSAL and NSEQ: the reference's genera and their numbers of sequences N_g (one line per document, as
experiments/silva_queue.sh writes Q/nseq.txt).  SEQTAX: MicrobeMixer truth files (genus path, id, read pairs k_g).
For each file, over the genera with reads: the least-squares slope of log k_g on log N_g (alpha if k_g ~ N_g^alpha),
Spearman's rank correlation of k_g and N_g, the share of the reads in the genera's share of the database (sum of
their N_g over N), and how many of the reference's genera have reads.  Genera of a truth file that are not in the
reference are counted and skipped.
"""
import math, sys

if len(sys.argv) < 4:
    sys.exit(__doc__)
genera = [l.split(None, 1)[1].strip() for l in open(sys.argv[1]) if l.strip()]
nseq = [int(l) for l in open(sys.argv[2]) if l.strip()]
assert len(genera) == len(nseq), "the two files differ in length"
N = dict(zip(genera, nseq)); NTOT = sum(nseq)


def ranks(v):
    o = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]: j += 1
        for t in range(i, j + 1): r[o[t]] = (i + j) / 2
        i = j + 1
    return r


def pearson(x, y):
    mx, my = sum(x) / len(x), sum(y) / len(y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y)); sx = sum((a - mx) ** 2 for a in x); sy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sx * sy) if sx and sy else float("nan")


for f in sys.argv[3:]:
    k = {}; unknown = 0
    for l in open(f):
        t = l.rstrip("\n").split("\t")
        if len(t) < 3 or not t[2].strip(): continue
        g, c = t[0].strip(), int(t[2])
        if c <= 0: continue
        if g not in N: unknown += 1; continue
        k[g] = k.get(g, 0) + c
    g = sorted(k)
    x = [math.log(N[h]) for h in g]; y = [math.log(k[h]) for h in g]
    mx, my = sum(x) / len(x), sum(y) / len(y)
    sxx = sum((a - mx) ** 2 for a in x)
    slope = sum((a - mx) * (b - my) for a, b in zip(x, y)) / sxx if sxx else float("nan")
    rho = pearson(ranks([N[h] for h in g]), ranks([k[h] for h in g]))
    share = sum(N[h] for h in g) / NTOT
    print("%s: %d read pairs over %d of %d genera (%d not in the reference); fitted alpha %.3f (log k_g on log N_g); "
          "Spearman %.3f; these genera hold %.1f%% of the database's sequences; N_g of the sampled genera %d..%d (median %d)"
          % (f, sum(k.values()), len(g), len(genera), unknown, slope, rho, 100 * share,
             min(N[h] for h in g), max(N[h] for h in g), sorted(N[h] for h in g)[len(g) // 2]))
