#!/usr/bin/env python3
# Synthetic queries mostly absent from a collection: random sequence with a few segments copied (with a
# few substitutions) from the collection's FASTA files.
#   gen_novel.py list.tsv out.fa [queries=20] [length=100000] [segments=3] [segment length=3000] [seed=1]
import random, sys
lst, out = sys.argv[1], sys.argv[2]
nq, ln, ns, sl, seed = [int(x) for x in (sys.argv[3:] + ['20', '100000', '3', '3000', '1'][len(sys.argv) - 3:])][:5]
random.seed(seed)
seqs = []
for line in open(lst):
    path = line.rstrip('\n').split('\t')[-1]
    s = ''.join(l.strip().upper() for l in open(path) if not l.startswith('>'))
    seqs.append(''.join(c for c in s if c in 'ACGT'))
with open(out, 'w') as f:
    for q in range(nq):
        s = [random.choice('ACGT') for _ in range(ln)]
        for _ in range(ns):
            g = random.choice(seqs); a = random.randrange(len(g) - sl); seg = list(g[a:a + sl])
            for i in range(len(seg)):
                if random.random() < 0.001: seg[i] = random.choice('ACGT')
            p = random.randrange(ln - sl); s[p:p + sl] = seg
        f.write('>q%d\n%s\n' % (q, ''.join(s)))
