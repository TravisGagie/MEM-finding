#!/usr/bin/env python3
"""16S read pairs with exactly k_g pairs from each genus g, k_g proportional to N_g^alpha, simulated as MicrobeMixer
does (tool/NOTES.md, "SILVA: credit schemes when the sample doesn't match the database").

    simulate_dial.py --mm MICROBEMIXER --fasta SILVA_FASTA --tax SILVA_TAX --primers PRIMERS \\
        --doc-to-traversal D2T --nseq NSEQ --alpha ALPHA --num-reads N --name REGION --out-dir DIR [--seed S] [--threads T]

N_g is the number of sequences of the reference's document for genus g (NSEQ, one line per document of D2T, as
experiments/silva_queue.sh writes it), the same table the scorer's normalised and prior schemes use.

MicrobeMixer's own code does the simulation: it reads SILVA (process_reference_file), finds the region between the
primers and mutates 2% of its bases (extract_variable_regions), and runs ART on each genus's regions with its own
command line (run_art_command: art_illumina -ss MSv1 -amp -p -na -l 250 -f FOLD).  What changes is the number of
read pairs per genus, which MicrobeMixer cannot set exactly from shares:
  1. a genus is eligible if it has at least one region of at least 250 bases, the read length (ART makes no reads
     from a shorter amplicon);
  2. k_g = N * w_g / sum(w) over the eligible genera, w_g = N_g^alpha (exact rationals), rounded by largest
     remainder (ties in document order), so the k_g sum to N exactly; genera with k_g = 0 get no reads;
  3. ART runs on each genus with k_g > 0, with FOLD = ceil(k_g / long regions) + 1 read pairs per amplicon; a genus
     with too few pairs is rerun with twice the fold (up to 3 times); one that still falls short is made
     ineligible and step 2 is redone, rerunning ART only where a genus now needs more pairs than it has;
  4. exactly k_g pairs are sampled from each genus's ART output and written in document order, as MicrobeMixer
     names them (@genus_<id>_read_<n>_mate1/2, id = document number), with the truth file <name>_seqtax.txt
     (genus path, id, k_g for every genus of D2T) and <name>_counts.tsv (id, N_g, regions, long regions, k_g).
"""
import argparse, math, os, random, shutil, subprocess, sys
from fractions import Fraction
import multiprocessing
# fork, so the workers share the regions without pickling and the script is not re-run in them (Python 3.14
# defaults to forkserver on Linux, which re-imports the main module; this script has no __main__ guard)
Pool = multiprocessing.get_context("fork").Pool

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in ("--mm", "--fasta", "--tax", "--primers", "--doc-to-traversal", "--nseq", "--name", "--out-dir"):
    ap.add_argument(x, required=True)
ap.add_argument("--alpha", type=float, required=True)
ap.add_argument("--num-reads", type=int, required=True)
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--threads", type=int, default=os.cpu_count())
ap.add_argument("--read-length", type=int, default=250, help="ART's read length in MicrobeMixer's command (250)")
a = ap.parse_args()
sys.path.insert(0, os.path.join(os.path.abspath(a.mm), "src"))
from stats_requests import parse_taxonomy_file          # MicrobeMixer's
from simulate import process_reference_file, process_primer_file, extract_variable_regions, run_art_command

random.seed(a.seed)
D = os.path.abspath(a.out_dir) + "/"
os.makedirs(D, exist_ok=True)
log = lambda *x: print(*x, file=sys.stderr, flush=True)

genera = [l.split(None, 1)[1].strip() for l in open(a.doc_to_traversal) if l.strip()]
nseq = [int(l) for l in open(a.nseq) if l.strip()]
if len(nseq) != len(genera):
    sys.exit("%s has %d lines, %s %d genera" % (a.nseq, len(nseq), a.doc_to_traversal, len(genera)))
G = len(genera)

# MicrobeMixer's reading and extraction, for every genus of the reference
_, _, lineage = parse_taxonomy_file(a.tax, ignore_eukaryotes=False)
to_seqs = process_reference_file(a.fasta, {g: [] for g in genera}, lineage)
fwd, rev = process_primer_file(a.primers)
regions, _ = extract_variable_regions({g: [] for g in genera}, to_seqs, fwd, rev)
del to_seqs
nreg = [len(regions[g]) for g in genera]
nlong = [sum(1 for s in regions[g] if len(s) >= a.read_length) for g in genera]
eligible = [i for i in range(G) if nlong[i] > 0]
log("%s: %d genera; %d with a region, %d with one of at least %d bases" % (a.name, G, sum(1 for x in nreg if x), len(eligible), a.read_length))


def allocate(elig):
    """k_g for g in elig: N * w_g / sum(w), w_g = N_g^alpha, by largest remainder (ties in document order)"""
    al = a.alpha
    if al == int(al):
        w = {i: Fraction(nseq[i]) ** int(al) for i in elig}
    else:
        w = {i: Fraction(nseq[i] ** al) for i in elig}
    tot = sum(w.values())
    q = {i: a.num_reads * w[i] / tot for i in elig}
    k = {i: math.floor(q[i]) for i in elig}
    left = a.num_reads - sum(k.values())
    for i in sorted(elig, key=lambda i: (-(q[i] - k[i]), i))[:left]:
        k[i] += 1
    return k


def pairs(i):
    f = "%sgenus_%d_reads_1.fq" % (D, i + 1)
    if not os.path.exists(f):
        return 0
    with open(f) as h:
        return sum(1 for _ in h) // 4


def art(jobs):
    """jobs: (genus index, fold); MicrobeMixer's ART command for each"""
    for i, _ in jobs:
        with open("%sgenus_%d_seqs.fna" % (D, i + 1), "w") as o:     # as write_out_extracted_sequences
            for j, s in enumerate(regions[genera[i]]):
                o.write(">genus_%d_seq_%d\n%s\n" % (i + 1, j, s))
    with Pool(a.threads) as p:
        p.map(run_art_command, [(i + 1, D, fold) for i, fold in jobs])


have = {}                                       # genus index -> pairs ART has made for it
fold = {}
while True:
    k = allocate(eligible)
    need = [i for i in eligible if k[i] > have.get(i, 0)]
    for attempt in range(4):
        if not need:
            break
        jobs = []
        for i in need:
            f = math.ceil(k[i] / nlong[i]) + 1
            fold[i] = max(f, 2 * fold[i]) if i in fold else f
            jobs.append((i, fold[i]))
        log("%s: ART on %d genera (attempt %d)" % (a.name, len(jobs), attempt + 1))
        art(jobs)
        for i in need:
            have[i] = pairs(i)
        need = [i for i in need if k[i] > have[i]]
    if not need:
        break
    log("%s: %d genera still short of pairs after 4 attempts; made ineligible, reallocating" % (a.name, len(need)))
    eligible = [i for i in eligible if i not in set(need)]

# exactly k_g pairs per genus, in document order
with open(D + a.name + "_mate_1.fq", "w") as m1, open(D + a.name + "_mate_2.fq", "w") as m2:
    for i in eligible:
        if k[i] == 0:
            continue
        with open("%sgenus_%d_reads_1.fq" % (D, i + 1)) as h1, open("%sgenus_%d_reads_2.fq" % (D, i + 1)) as h2:
            r1 = [x.rstrip("\n") for x in h1]; r2 = [x.rstrip("\n") for x in h2]
        assert len(r1) == len(r2) and len(r1) % 4 == 0 and len(r1) // 4 >= k[i]
        for n, j in enumerate(random.sample(range(len(r1) // 4), k[i]), 1):
            m1.write("@genus_%d_read_%d_mate1\n%s\n%s\n%s\n" % (i + 1, n, r1[4 * j + 1], r1[4 * j + 2], r1[4 * j + 3]))
            m2.write("@genus_%d_read_%d_mate2\n%s\n%s\n%s\n" % (i + 1, n, r2[4 * j + 1], r2[4 * j + 2], r2[4 * j + 3]))
with open(D + a.name + "_seqtax.txt", "w") as o:
    for i, g in enumerate(genera):
        o.write("%s\t%d\t%d\n" % (g, i + 1, k.get(i, 0)))
with open(D + a.name + "_counts.tsv", "w") as o:
    o.write("id\tN_g\tregions\tlong_regions\tk_g\n")
    for i in range(G):
        o.write("%d\t%d\t%d\t%d\t%d\n" % (i + 1, nseq[i], nreg[i], nlong[i], k.get(i, 0)))
for f in os.listdir(D):                          # the per-genus files (os.listdir, not a glob: thousands of genera)
    if f.startswith("genus_"):
        os.remove(D + f)
total = sum(k.values())
assert total == a.num_reads, (total, a.num_reads)
nz = [k[i] for i in eligible if k[i] > 0]
log("%s: alpha = %g: %d read pairs over %d genera (%d to %d each)" % (a.name, a.alpha, total, len(nz), min(nz), max(nz)))
