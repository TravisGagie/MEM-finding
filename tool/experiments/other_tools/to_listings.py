#!/usr/bin/env python3
"""Kraken 2's or Tagger's output as KATKA listings, one file per mate, for score_silva.py:
each mate of a read (pair) assigned to a genus gets the listing "[0,len-1] {doc}" (doc: KATKA's document number, from 0),
so that it votes for that genus alone; a read assigned above the genus level, or not at all, gets an empty listing.
  to_listings.py kraken2 OUT.kraken READS_1.fq OUT_1 OUT_2 TAX.txt DOC_TO_TRAVERSAL   (--paired output, one line per pair)
  to_listings.py tagger OUT_1.tsv OUT_2.tsv READS_1.fq READS_2.fq OUT_1 OUT_2          (one output per mate)"""
import sys, gzip
op = lambda f: gzip.open(f, "rt") if f.endswith(".gz") else open(f)
def reads(f):
    with op(f) as h:
        for i, l in enumerate(h):
            if i % 4 == 0: name = l[1:].split()[0]
            elif i % 4 == 1: yield name, len(l.strip())
def write(path, items, doc_of):
    with open(path, "w") as o:
        for name, n in items:
            d = doc_of(name)
            o.write(">%s\n%s\n" % (name, "[0,%d] {%d}" % (n - 1, d) if d is not None else ""))
if sys.argv[1] == "kraken2":
    kout, r1, o1, o2, tax, d2t = sys.argv[2:8]
    path = {}; rank = {}
    for l in open(tax):
        p, i, rk = l.rstrip("\n").split("\t")[:3]; path[int(i)] = p; rank[p] = rk
    doc = {l.split(None, 1)[1].strip(): int(l.split()[0]) - 1 for l in open(d2t) if l.strip()}
    def genus(t):
        p = path.get(t)
        while p:
            if rank.get(p) == "genus": return doc.get(p)
            p = p[:p[:-1].rfind(";") + 1]
        return None
    g = {}; nc = ng = 0
    for l in open(kout):
        f = l.split("\t")
        if f[0] == "C":
            nc += 1; d = genus(int(f[2]))
            if d is not None: g[f[1]] = d; ng += 1
    R = list(reads(r1))
    write(o1, R, lambda n: g.get(n))
    write(o2, [(n.replace("mate1", "mate2"), k) for n, k in R], lambda n: g.get(n.replace("mate2", "mate1")))
    print("pairs %d, classified %d, at a genus %d" % (len(R), nc, ng), file=sys.stderr)
else:
    t1, t2, r1, r2, o1, o2 = sys.argv[2:8]
    for t, r, o in ((t1, r1, o1), (t2, r2, o2)):
        tag = {}
        for l in open(t):
            x = l.split()
            if len(x) >= 2 and x[1].isdigit() and int(x[1]) < 65535: tag[x[0].lstrip("@>")] = int(x[1])   # 65535: no tag
        write(o, reads(r), lambda n: tag.get(n))
