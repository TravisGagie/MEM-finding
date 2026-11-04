#!/usr/bin/env python3
"""A small synthetic stand-in for the HPRC chromosome 19 data, for checking tool/CHR19.md end to end.

    python3 chr19_synthetic.py OUTDIR [--length 300000] [--individuals 8] [--seed 1]

Writes OUTDIR/collection.txt (one sequence per line: a random "reference", then both haplotypes of each
individual) and OUTDIR/queries.txt (haplotype 1 of two more individuals), in the format of the real data
set's collection.txt and queries.txt (tool/CHR19.md, "Data").  Each haplotype is the reference with SNPs
(1 in 500), short indels (1 in 5000), a few runs of N, and variants shared within a population of
founder haplotypes, so that matches are long and the collection is repetitive, as on chromosome 19.
"""
import argparse, os, random


def mutate(seq, rng, snp=1 / 500, indel=1 / 5000):
    out, i, n = [], 0, len(seq)
    while i < n:
        r = rng.random()
        if r < snp:
            out.append(rng.choice("ACGT".replace(seq[i], ""))); i += 1
        elif r < snp + indel:
            if rng.random() < 0.5: i += rng.randint(1, 10)                       # deletion
            else: out.append(seq[i]); out.extend(rng.choice("ACGT") for _ in range(rng.randint(1, 10))); i += 1
        else:
            out.append(seq[i]); i += 1
    return "".join(out)


def with_ns(seq, rng):
    s = list(seq)
    for _ in range(2):                                                          # two runs of N
        a = rng.randrange(len(s) - 200); s[a:a + rng.randint(20, 200)] = "N" * rng.randint(20, 200)
    return "".join(s)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--length", type=int, default=300000)
    ap.add_argument("--individuals", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    os.makedirs(a.outdir, exist_ok=True)
    ref = "".join(rng.choice("ACGT") for _ in range(a.length))
    founders = [mutate(ref, rng) for _ in range(4)]                              # shared variation
    def haplotype():
        return with_ns(mutate(rng.choice(founders), rng, snp=1 / 2000, indel=1 / 20000), rng)
    with open(os.path.join(a.outdir, "collection.txt"), "w") as f:
        f.write(ref + "\n")
        for _ in range(2 * a.individuals):
            f.write(haplotype() + "\n")
    with open(os.path.join(a.outdir, "queries.txt"), "w") as f:
        for _ in range(2):
            f.write(haplotype() + "\n")
    print("wrote %s/collection.txt (%d sequences) and queries.txt (2)" % (a.outdir, 2 * a.individuals + 1))


if __name__ == "__main__":
    main()
