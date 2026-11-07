#!/usr/bin/env python3
"""The held-out split: for each aquatic genus, a random 10% of its sequences, at least 1, go to the simulation and the rest to the index."""
import math, os, random
H = os.path.expanduser("~/tool/silva/heldout"); os.makedirs(H + "/docs", exist_ok=True)
S = os.path.expanduser("~/rz/MEM-finding/katka/silva"); rnd = random.Random(20261107)
d2t = [l.split(None, 1)[1].strip() for l in open(S + "/work/ref/rna/doc_to_traversal.txt") if l.strip()]
aq = [l.split("\t")[0] for l in list(open(S + "/work/abundance/aquatic.tsv"))[1:]]
docs = {g: i + 1 for i, g in enumerate(d2t)}
held = {}; nseq = [int(l) for l in open(os.path.expanduser("~/tool/silva/queue/nseq.txt"))]
lst = [l.rstrip("\n").split("\t") for l in open(os.path.expanduser("~/tool/silva/list.tsv"))]
def records(f):
    h = None; b = []
    for l in open(f):
        if l[0] == ">":
            if h: yield h, "".join(b)
            h = l; b = []
        else: b.append(l)
    if h: yield h, "".join(b)
kept = 0
with open(H + "/split.tsv", "w") as t:
    for g in aq:
        d = docs[g]; f = "%s/work/ref/dna/doc_%d_seq.fa" % (S, d)
        R = list(records(f)); rnd.shuffle(R); k = max(1, math.ceil(len(R) * 0.1))
        with open("%s/docs/doc_%d_seq.fa" % (H, d), "w") as o:
            for h, b in R[k:]: o.write(h + b)
        for h, b in R[:k]: held[h[1:].split()[0]] = g
        nseq[d - 1] = len(R) - k; lst[d - 1][1] = "%s/docs/doc_%d_seq.fa" % (H, d); kept += len(R) - k
        t.write("%s\t%d\t%d\t%d\n" % (g, d, len(R), k))
        if len(R) == k: print("genus left empty in the index:", g)
open(H + "/list.tsv", "w").write("".join("\t".join(x) + "\n" for x in lst))
open(H + "/nseq.txt", "w").write("".join("%d\n" % n for n in nseq))
n = w = 0; keep = False
with open(H + "/heldout.fasta", "w") as o:
    for l in open(S + "/exp1_data/SILVA_138.1_SSURef_NR99_tax_silva.fasta"):
        if l[0] == ">":
            keep = l[1:].split()[0] in held; n += 1; w += keep
        if keep: o.write(l)
print("aquatic genera %d: %d sequences in the index, %d held out, %d of them found in the FASTA of %d" % (len(aq), kept, len(held), w, n))
