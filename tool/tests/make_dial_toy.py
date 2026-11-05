#!/usr/bin/env python3
"""A toy SILVA for the credit-scheme experiment: genera of very different sizes N_g that share sequence.

    make_dial_toy.py DOC_TO_TRAVERSAL OUT [--seed S]

DOC_TO_TRAVERSAL: real SILVA 138.1 genus paths, one per line ("<doc> <path>", as katka/silva builds it), at least
12 of them.  OUT gets, for the first 12 genera:
  silva.fasta            SILVA-style records (RNA, ">ACC.1.1500 <path><species>"), N_g of them for genus g;
  ref/doc_<i>_seq.fa     the reference's documents (DNA), list.tsv for katka prep;
  doc_to_traversal.txt, trav_to_length.txt, nseq.txt (N_g per document).
The 12 genera come in 3 groups of 4 that descend from one ancestor each (0.6% of the bases differ between genera of
a group, 0.3% between sequences of a genus; --between, --within), so long matches hit several genera of a group, mostly in proportion to
their sizes.  Each group has one large genus and three small ones, in a different order per group.  Every record
has the V1_V2, V3_V4 and V4 primer sites of MicrobeMixer's data/*_primers.txt (one IUPAC choice each), around
regions of 300 to 330 bases.
"""
import argparse, os, random

ap = argparse.ArgumentParser(); ap.add_argument("d2t"); ap.add_argument("out"); ap.add_argument("--seed", type=int, default=3)
ap.add_argument("--between", type=float, default=0.006, help="bases that differ between genera of a group")
ap.add_argument("--within", type=float, default=0.003, help="bases that differ between sequences of a genus")
a = ap.parse_args()
rnd = random.Random(a.seed)
paths = [l.split(None, 1)[1].strip() for l in open(a.d2t) if l.strip()][:12]
assert len(paths) == 12
SIZES = [40, 2, 3, 1,   1, 30, 2, 4,   3, 1, 50, 2]     # N_g: one large genus per group of four
B = "ACGU"
def mutate(s, p):
    return "".join(rnd.choice(B.replace(c, "")) if rnd.random() < p else c for c in s)


# MicrobeMixer's primer patterns as they appear in an RNA record (its process_primer_file: the forward primer, and
# the reverse complement of the reverse one), with the first base of each IUPAC class
def site(p):
    out, i = "", 0
    while i < len(p):
        if p[i] == "[": j = p.index("]", i); out += p[i + 1]; i = j + 1
        else: out += p[i]; i += 1
    return out
F12, R12 = site("AG[AC]GUU[CU]GAU[CU][AC]UGGCUCAG"), site("ACUCCUACGGGAGGCAGC")
F34, R34 = site("CCUACGGG[ACGU]GGC[AU]GCAG"), site("GGAUUAGAUACCC[CGU][AGU]GUAGUC")
F4, R4 = site("GUGCCAGC[AC]GCCGCGGUAA"), site("GGAUUAGAUACCC[CGU][AGU]GUAGUC")
R45 = site("ACU[CU]AAA[GU]GAAUUGACGGG")
assert R12 == "ACU" + F34[:15] and R34 == R4


def ancestor():
    """V1_V2: F12 [310] R12, where R12 overlaps F34 (as 338R and 341F do); V3_V4: F34 [120] F4 [300] R4 (439 bases);
    V4: F4 [300] R4; V4_V5 (whose forward primer is V4's): F4 [300] R4 [80] R45 (401 bases)"""
    r = lambda n: "".join(rnd.choice(B) for _ in range(n))
    return r(40) + F12 + r(310) + "ACU" + F34 + r(120) + F4 + r(300) + R4 + r(80) + R45 + r(60)


os.makedirs(os.path.join(a.out, "ref"), exist_ok=True)
anc = [ancestor() for _ in range(3)]
recs = []; acc = 0
for gi, (path, n) in enumerate(zip(paths, SIZES)):
    base = mutate(anc[gi // 4], a.between)
    seqs = []
    for j in range(n):
        acc += 1
        s = mutate(base, a.within)
        seqs.append(s); recs.append((">AB%06d.1.1500 %s%s sp.\n" % (acc, path, path.rstrip(";").split(";")[-1]), s))
    with open(os.path.join(a.out, "ref", "doc_%d_seq.fa" % (gi + 1)), "w") as o:
        for j, s in enumerate(seqs):
            o.write(">g%d_s%d\n%s\n" % (gi + 1, j, s.replace("U", "T")))
with open(os.path.join(a.out, "silva.fasta"), "w") as o:
    for h, s in recs:
        o.write(h + s + "\n")
with open(os.path.join(a.out, "list.tsv"), "w") as o:
    for gi in range(12):
        o.write("silva\t%s\n" % os.path.abspath(os.path.join(a.out, "ref", "doc_%d_seq.fa" % (gi + 1))))
with open(os.path.join(a.out, "doc_to_traversal.txt"), "w") as o:
    for gi, p in enumerate(paths):
        o.write("%d %s\n" % (gi + 1, p))
with open(os.path.join(a.out, "trav_to_length.txt"), "w") as o:
    for gi, p in enumerate(paths):
        o.write("%s %d\n" % (p, len(anc[gi // 4]) * SIZES[gi]))
with open(os.path.join(a.out, "nseq.txt"), "w") as o:
    o.write("".join("%d\n" % n for n in SIZES))
print("%d records, %d genera (sizes %s)" % (len(recs), 12, SIZES))
