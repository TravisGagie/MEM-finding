#!/usr/bin/env python3
# norm_credit.py listings truth tbl
# Like credit_schemes.py, but compares proportional and even credit with proportional credit normalized by the
# number of reference sequences of each tag (norm_rec) or by their total length (norm_len); also prints both per tag.
# tbl columns: index, tag, sequence name, length, start.
import re, sys, collections
out, truthf, tbl = sys.argv[1:4]
names = []; nrec = collections.Counter(); blen = collections.Counter()
for l in open(tbl):
    f = l.rstrip('\n').split('\t'); s = f[1]
    if s not in names: names.append(s)
    nrec[s] += 1; blen[s] += int(f[3])
truth = dict(l.split()[:2] for l in open(truthf))
print('tag  records  bases')
for s in names: print(' ', s, nrec[s], blen[s])
pat = re.compile(r'\[(\d+),(\d+)\] \{([^}]*)\}')
S = ('prop', 'even', 'norm_rec', 'norm_len')
ok = collections.Counter(); n = 0; name = None
def best(w):
    if not w: return None
    b = max(w.values()); top = [t for t in w if w[t] == b]
    return top[0] if len(top) == 1 else None
for line in open(out):
    line = line.rstrip('\n')
    if line.startswith('>'): name = line[1:]; continue
    n += 1; W = {s: collections.Counter() for s in S}
    for a, b, d in pat.findall(line):
        ln = int(b) - int(a) + 1
        cs = [(int(x.split(':')[0]), int(x.split(':')[1]) if ':' in x else 1) for x in d.split(',') if x]
        if not cs: continue
        tot = sum(c for _, c in cs)
        fr = [c / nrec[names[t]] for t, c in cs]; fl = [c / blen[names[t]] for t, c in cs]
        for (t, c), x, y in zip(cs, fr, fl):
            W['prop'][t] += ln * c / tot; W['even'][t] += ln / len(cs)
            W['norm_rec'][t] += ln * x / sum(fr); W['norm_len'][t] += ln * y / sum(fl)
    for s in S:
        t = best(W[s]); ok[s] += t is not None and names[t] == truth[name]
for s in S: print('  %-9s %7.3f%%' % (s, 100 * ok[s] / n))
