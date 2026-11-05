#!/usr/bin/env python3
"""Which alpha a read set's genus composition corresponds to: k_g against N_g (tool/NOTES.md, "SILVA: credit
schemes when the sample doesn't match the database").

    silva_alpha_fit.py DOC_TO_TRAVERSAL NSEQ SEQTAX [SEQTAX ...]

DOC_TO_TRAVERSAL and NSEQ: the reference's genera and their numbers of sequences N_g (one line per document, as
experiments/silva_queue.sh writes Q/nseq.txt).  SEQTAX: MicrobeMixer truth files (genus path, id, read pairs k_g).
For each file, over the genera with reads: the least-squares slope of log k_g on log N_g (alpha if k_g ~ N_g^alpha),
Spearman's rank correlation of k_g and N_g, the fraction of the database's sequences in those genera (sum of their
N_g over N), and their ranks by N_g among all the reference's genera, with or without reads (median, range, and how
many are in the top 1% and 10%).  Genera of a truth file that are not in the
reference are counted and skipped.
"""
import math, sys

if len(sys.argv) < 4:
    sys.exit(__doc__)
genera = [l.split(None, 1)[1].strip() for l in open(sys.argv[1]) if l.strip()]
nseq = [int(l) for l in open(sys.argv[2]) if l.strip()]
assert len(genera) == len(nseq), "the two files differ in length"
N = dict(zip(genera, nseq)); NTOT = sum(nseq)
# rank of each genus by N_g among all the reference's genera, largest first; ties share the best rank (1 + the number
# of genera with strictly more sequences)
_desc = sorted(nseq, reverse=True)
_first = {}
for i, n in enumerate(_desc): _first.setdefault(n, i + 1)
RANK = {h: _first[N[h]] for h in genera}


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
    rk = sorted(RANK[h] for h in g)
    med = (rk[(len(rk) - 1) // 2] + rk[len(rk) // 2]) / 2
    print("%s: %d read pairs over %d of %d genera (%d not in the reference)" % (f, sum(k.values()), len(g), len(genera), unknown))
    print("  fitted alpha %.3f (least squares, log k_g on log N_g, over the genera with reads); Spearman %.3f" % (slope, rho))
    print("  these genera hold %d of the database's %d sequences (%.1f%%); N_g %d..%d (median %g)"
          % (sum(N[h] for h in g), NTOT, 100 * share, min(N[h] for h in g), max(N[h] for h in g),
             (sorted(N[h] for h in g)[(len(g) - 1) // 2] + sorted(N[h] for h in g)[len(g) // 2]) / 2))
    print("  their ranks by N_g among all %d genera (1 = largest; ties share their best rank; genera without reads"
          " included): median %g, range %d..%d; %d of them in the top 1%%, %d in the top 10%%"
          % (len(genera), med, rk[0], rk[-1], sum(1 for x in rk if x <= len(genera) / 100), sum(1 for x in rk if x <= len(genera) / 10)))
