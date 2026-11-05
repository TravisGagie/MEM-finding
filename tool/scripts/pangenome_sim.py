#!/usr/bin/env python3
"""Synthetic pangenomes with a dial for the number of haplotypes, from a population model (tool/NOTES.md, "Block
length as genomes are added").

    pangenome_sim.py OUT [--haplotypes N] [--queries Q] [--length L] [--seed S] [other parameters below]

Haplotypes are a sample of N + Q + 1 from a coalescent with recombination (msprime), with human-like parameters:
  --ne 20000            effective population size before the recent growth (with --mu and the growth, diversity about
                          4e-4 to 5e-4 per base in SNPs and short indels, varying with the seed; human: 6e-4 to 1e-3)
  --growth-ne 1000000   present population size; the population grows exponentially from --ne to this over the
  --growth-gens 400       last --growth-gens generations (recent human growth, which makes most variants rare and
                          gives each new haplotype its own private ones)
  --mu 1.25e-8          mutations per base per generation
  --rho 1e-8            recombinations per base per generation (about 1 cM per Mb)
  --indels 0.1          fraction of variant sites that are indels rather than SNPs; an indel is a deletion or
  --indel-max 10          an insertion (equally likely) of 1..indel-max bases, length geometric with mean 2
  --length 300000       bases per haplotype (before indels)
Sites hit more than once keep the first mutation's allele model (msprime's binary model: a site is 0 or 1 per
haplotype).  The ancestral sequence is random (A/C/G/T uniform); every variant is applied to it, so haplotype h
carries the derived allele at the sites where msprime gives it 1.  Indels that would overlap an earlier variant
are dropped (counted).

Output (OUT): collection.txt, one haplotype per line: the reference (sample 0) first, then samples 1..N, in
order, so the first n + 1 lines are the collection with n haplotypes for every n <= N, as chr19's
collection.txt (CHM13 first); queries.txt, the Q held-out haplotypes (samples N + 1 .. N + Q), one per line;
stats.txt: the parameters, segregating sites, SNPs and indels, the site frequency spectrum's low end, nucleotide
diversity, and the variants that are private to each haplotype in the collection (counted among the N + 1).
"""
import argparse, math, os, random, sys

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("out")
ap.add_argument("--haplotypes", type=int, default=1024)
ap.add_argument("--queries", type=int, default=2)
ap.add_argument("--length", type=int, default=300000)
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--ne", type=float, default=2e4)
ap.add_argument("--growth-ne", type=float, default=1e6)
ap.add_argument("--growth-gens", type=float, default=400)
ap.add_argument("--mu", type=float, default=1.25e-8)
ap.add_argument("--rho", type=float, default=1e-8)
ap.add_argument("--indels", type=float, default=0.1)
ap.add_argument("--indel-max", type=int, default=10)
a = ap.parse_args()
import msprime

os.makedirs(a.out, exist_ok=True)
n = a.haplotypes + 1 + a.queries
dem = msprime.Demography()
g = math.log(a.growth_ne / a.ne) / a.growth_gens if a.growth_gens > 0 and a.growth_ne > a.ne else 0.0
dem.add_population(name="pop", initial_size=a.growth_ne if g else a.ne, growth_rate=g)
if g:
    dem.add_population_parameters_change(time=a.growth_gens, initial_size=a.ne, growth_rate=0, population="pop")
ts = msprime.sim_ancestry(samples={"pop": n}, ploidy=1, demography=dem, sequence_length=a.length,
                          recombination_rate=a.rho, random_seed=a.seed)
ts = msprime.sim_mutations(ts, rate=a.mu, model=msprime.BinaryMutationModel(), discrete_genome=True, random_seed=a.seed)
rnd = random.Random(a.seed)
B = "ACGT"
anc = "".join(rnd.choice(B) for _ in range(a.length))

# variants: (position, kind, payload); carriers per variant
variants = []; dropped = 0; last_end = -1
for v in ts.variants():
    pos = int(v.site.position)
    if pos <= last_end:
        dropped += 1; continue
    carriers = [i for i, x in enumerate(v.genotypes) if x == 1]
    if not carriers or len(carriers) == n:
        continue
    if rnd.random() < a.indels:
        k = 1
        while k < a.indel_max and rnd.random() < 0.5: k += 1
        if rnd.random() < 0.5 and pos + k < a.length:
            variants.append((pos, "del", k, carriers)); last_end = pos + k - 1
        else:
            variants.append((pos, "ins", "".join(rnd.choice(B) for _ in range(k)), carriers)); last_end = pos
    else:
        variants.append((pos, "snp", rnd.choice(B.replace(anc[pos], "")), carriers)); last_end = pos

per = [[] for _ in range(n)]
for idx, (pos, kind, x, carriers) in enumerate(variants):
    for c in carriers: per[c].append(idx)


def haplotype(h):
    out = []; at = 0
    for idx in per[h]:
        pos, kind, x, _ = variants[idx]
        out.append(anc[at:pos])
        if kind == "snp": out.append(x); at = pos + 1
        elif kind == "del": at = pos + x
        else: out.append(anc[pos] + x); at = pos + 1
    out.append(anc[at:])
    return "".join(out)


with open(os.path.join(a.out, "collection.txt"), "w") as o:
    for h in range(a.haplotypes + 1):
        o.write(haplotype(h) + "\n")
with open(os.path.join(a.out, "queries.txt"), "w") as o:
    for h in range(a.haplotypes + 1, n):
        o.write(haplotype(h) + "\n")

# statistics, over the collection (samples 0..N)
C = a.haplotypes + 1
freq = [sum(1 for c in v[3] if c < C) for v in variants]
seg = [f for f in freq if 0 < f < C]
sfs = [sum(1 for f in seg if f == i) for i in range(1, 6)]
pi = sum(f * (C - f) for f in seg) / (C * (C - 1) / 2) / a.length
priv = [0] * C
for v, f in zip(variants, freq):
    if f == 1: priv[next(c for c in v[3] if c < C)] += 1
with open(os.path.join(a.out, "stats.txt"), "w") as o:
    o.write("pangenome_sim.py %s\n" % " ".join(sys.argv[1:]))
    o.write("parameters: haplotypes %d (+1 reference, +%d held out), length %d, seed %d, Ne %g, growth to %g over %g "
            "generations, mu %g, rho %g, indels %g (up to %d bases)\n" % (a.haplotypes, a.queries, a.length, a.seed, a.ne,
            a.growth_ne, a.growth_gens, a.mu, a.rho, a.indels, a.indel_max))
    o.write("variants: %d sites (%d SNPs, %d indels; %d dropped as overlapping); %d segregating in the collection\n" % (
        len(variants), sum(1 for v in variants if v[1] == "snp"), sum(1 for v in variants if v[1] != "snp"), dropped, len(seg)))
    o.write("site frequency spectrum in the collection, counts 1..5: %s; nucleotide diversity %.2e per base\n" % (sfs, pi))
    for k in (16, 64, 256, 1024):
        if k < C:
            o.write("private variants (among the first %d + 1) per haplotype, mean over haplotypes %d..%d: %.1f\n" % (
                k, max(1, k - 15), k, sum(priv[max(1, k - 15):k + 1]) / len(priv[max(1, k - 15):k + 1])))
    # new variants added by haplotype n, among the collection of the first n
    seen = set(); new = []
    for h in range(C):
        s = set(per[h]); new.append(len(s - seen)); seen |= s
    for k in (16, 64, 256, 1024):
        if k < C:
            o.write("variants first seen in haplotypes %d..%d: %.1f per haplotype\n" % (max(1, k - 15), k, sum(new[max(1, k - 15):k + 1]) / len(new[max(1, k - 15):k + 1])))
print(open(os.path.join(a.out, "stats.txt")).read(), end="")
