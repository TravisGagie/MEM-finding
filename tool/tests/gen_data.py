#!/usr/bin/env python3
# Small synthetic test data: a few "species", each a random base genome with several mutated copies
# (substitutions, indels, runs of N), one FASTA file per genome; reads sampled from mutated genomes,
# half of them reverse-complemented, with errors; and a chromosome-19-like collection.
#   gen_data.py OUTDIR [seed]
import os, random, sys
out = sys.argv[1]; seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
random.seed(seed)
os.makedirs(out, exist_ok=True)
def rnd(n): return ''.join(random.choice('ACGT') for _ in range(n))
def mut(s, rate):
    o = []; i = 0
    while i < len(s):
        r = random.random()
        if r < rate: o.append(random.choice('ACGT')); i += 1
        elif r < rate * 1.2: i += random.randint(1, 10)
        elif r < rate * 1.4: o.append(rnd(random.randint(1, 10)))
        elif r < rate * 1.402: o.append('N' * random.randint(1, 30))
        else: o.append(s[i]); i += 1
    return ''.join(o)
def rc(s): return s[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))
def fasta(path, recs):
    with open(path, 'w') as f:
        for name, s in recs:
            f.write('>' + name + '\n')
            for i in range(0, len(s), 70): f.write(s[i:i + 70] + '\n')
species = []
with open(os.path.join(out, 'list.tsv'), 'w') as L:
    for sp in range(4):
        base = rnd(random.randint(20000, 30000))
        # a shared region between species, so that tags overlap
        if sp > 0: base = base[:5000] + species[0][3000:4500] + base[5000:]
        species.append(base)
        for g in range(5):
            gs = mut(base, 0.01)
            recs = [('g%d_%d.1 chr' % (sp, g), gs[:len(gs) // 2]), ('g%d_%d.2 plasmid' % (sp, g), gs[len(gs) // 2:])]
            if sp == 3 and g == 4:   # short A/C/G/T runs between Ns: every 4-mer becomes a minimizer, so digests share a byte
                recs.append(('g3_4.3 scaffold', 'N'.join(rnd(random.randint(4, 7)) for _ in range(3000))))
            p = os.path.join(out, 'sp%d_g%d.fa' % (sp, g))
            fasta(p, recs)
            L.write('species%d\t%s\n' % (sp, p))
reads = []
for i in range(300):
    sp = random.randrange(4); g = mut(species[sp], 0.02)
    ln = random.choice([150, 250, 1000, 3000]); a = random.randrange(max(1, len(g) - ln))
    r = mut(g[a:a + ln], 0.01 if ln < 1000 else 0.05)
    if random.random() < 0.5: r = rc(r)
    if random.random() < 0.05: r = rnd(ln)              # novel
    if random.random() < 0.1: r = r[:ln // 2] + 'AAAAAAAAAAAAAAAAAAAA' + r[ln // 2:]
    if random.random() < 0.1: r = r.lower()
    reads.append(('read%d sp%d' % (i, sp), r))
with open(os.path.join(out, 'reads.fq'), 'w') as f:
    for n, r in reads: f.write('@%s\n%s\n+\n%s\n' % (n, r, 'I' * len(r)))
fasta(os.path.join(out, 'reads.fa'), reads)
# chromosome-19-like: one sequence, 8 mutated copies (one FASTA each), 2 query haplotypes
os.makedirs(os.path.join(out, 'chr'), exist_ok=True)
ref = rnd(60000)
with open(os.path.join(out, 'chr', 'list.tsv'), 'w') as L:
    for k in range(8):
        p = os.path.join(out, 'chr', 'h%d.fa' % k)
        fasta(p, [('h%d' % k, mut(ref, 0.002))])
        L.write('-\t%s\n' % p)
fasta(os.path.join(out, 'chr', 'queries.fa'), [('q%d' % k, mut(ref, 0.004)) for k in range(2)] + [('q_rc', rc(mut(ref, 0.004)))])
