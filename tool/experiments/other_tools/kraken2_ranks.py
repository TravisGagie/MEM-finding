#!/usr/bin/env python3
"""Kraken 2's accuracy by how far up the true lineage its call lies, from its --paired output (columns C/U, name, taxon).
  kraken2_ranks.py TAX.txt TRUTH_seqtax.txt OUT.kraken[.gz]
A pair is 'genus' if Kraken 2 assigns the true genus, 'below the domain' if it assigns the true genus or one of its
ancestors other than the domain (phylum to family, or any intermediate SILVA rank), 'any ancestor' if the domain counts too.
A taxon below the genus counts as its genus.  Prints the three fractions and the counts of unclassified and wrong pairs."""
import sys, gzip
tax, truth, kout = sys.argv[1:4]
path = {}; rank = {}
for l in open(tax):
    p, i, rk = l.rstrip("\n").split("\t")[:3]; path[int(i)] = p; rank[p] = rk
gen = {}
for l in open(truth):
    t = l.split(); k = max(i for i in range(len(t)) if ";" in t[i]); gen[int(t[k + 1])] = " ".join(t[:k + 1])
g = b = a = u = w = n = 0
for l in (gzip.open(kout, "rt") if kout.endswith(".gz") else open(kout)):
    f = l.rstrip("\n").split("\t"); n += 1
    tg = gen[int(f[1].split("_")[1])]
    if f[0] != "C" or f[2] == "0": u += 1; continue
    p = path.get(int(f[2]), "")
    if not p: w += 1; continue
    if tg.startswith(p) or p.startswith(tg):           # an ancestor of the true genus, or the genus or below it
        depth = p.count(";")
        if p.startswith(tg): g += 1; b += 1; a += 1
        elif depth >= 2: b += 1; a += 1
        else: a += 1
    else: w += 1
print("pairs %d genus %.4f below_domain %.4f any_ancestor %.4f unclassified %d wrong %d" % (n, g / n, b / n, a / n, u, w))
