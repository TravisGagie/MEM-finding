#!/usr/bin/env python3
"""Which genera can MicrobeMixer simulate?  One line per genus of the reference, for the sparse dial's draws.

    sparse_genera.py --mm MM --fasta F --tax T --doc-to-traversal D2T --nseq NSEQ --primers-dir DIR --out OUT.tsv

For each region, with MicrobeMixer's own reading and extraction: how many sequences give the region (nreg) and how
many of those are at least 250 bases long (nlong), the ones ART can read.  Columns: document (0-based), N_g, path,
then nreg and nlong per region.
"""
import argparse, os, sys, random
ap = argparse.ArgumentParser()
for x in ("--mm", "--fasta", "--tax", "--doc-to-traversal", "--nseq", "--primers-dir", "--out"): ap.add_argument(x, required=True)
a = ap.parse_args()
sys.path.insert(0, os.path.join(os.path.abspath(a.mm), "src"))
from stats_requests import parse_taxonomy_file
from simulate import process_reference_file, process_primer_file, extract_variable_regions
random.seed(1)
genera = [l.split(None, 1)[1].strip() for l in open(a.doc_to_traversal) if l.strip()]
nseq = [int(l) for l in open(a.nseq) if l.strip()]
_, _, lineage = parse_taxonomy_file(a.tax, ignore_eukaryotes=False)
to_seqs = process_reference_file(a.fasta, {g: [] for g in genera}, lineage)
cols = []; REG = [("V1_V2", "V1_V2"), ("V3_V4", "V3_V4"), ("V4_V4", "V4"), ("V4_V5", "V4_V5")]
for r, p in REG:
    fwd, rev = process_primer_file(os.path.join(a.primers_dir, p + "_primers.txt"))
    regs, _ = extract_variable_regions({g: [] for g in genera}, to_seqs, fwd, rev)
    cols.append(([len(regs[g]) for g in genera], [sum(1 for s in regs[g] if len(s) >= 250) for g in genera]))
    print(r, sum(1 for x in cols[-1][0] if x), sum(1 for x in cols[-1][1] if x), file=sys.stderr, flush=True)
with open(a.out, "w") as o:
    o.write("doc\tN_g\tpath\t" + "\t".join("nreg_%s\tnlong_%s" % (r, r) for r, _ in REG) + "\n")
    for i, g in enumerate(genera):
        o.write("%d\t%d\t%s\t%s\n" % (i, nseq[i], g, "\t".join("%d\t%d" % (c[0][i], c[1][i]) for c in cols)))
