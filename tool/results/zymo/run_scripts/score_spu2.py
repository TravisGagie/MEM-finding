# usage: python3 score_spu2.py listings truth tbl
import re, sys, collections
out, truthf, tbl = sys.argv[1:4]
sp = []
for l in open(tbl):
    s = l.split('\t')[1]
    if s not in sp: sp.append(s)
truth = dict(l.split() for l in open(truthf))
pat = re.compile(r'\[(\d+),(\d+)\] \{([^}]*)\}')
ok = n = unc = 0; conf = collections.Counter(); per = collections.Counter(); name = None
for line in open(out):
    line = line.rstrip('\n')
    if line.startswith('>'): name = line[1:]; continue
    n += 1; w = collections.Counter(); per[truth[name]] += 1
    for a, b, d in pat.findall(line):
        ln = int(b) - int(a) + 1; cs = [x.split(':') for x in d.split(',') if x]
        cs = [(int(x[0]), int(x[1]) if len(x) > 1 else 1) for x in cs]; tot = sum(c for _, c in cs)
        for t, c in cs: w[t] += ln * c / tot
    if not w: unc += 1; conf[(truth[name], 'none')] += 1; continue
    best = max(w.values()); top = [t for t in w if w[t] == best]
    pred = sp[top[0]] if len(top) == 1 else 'tie'
    ok += pred == truth[name]; conf[(truth[name], pred)] += 1
print('reads %d accuracy %.3f%% unclassified %d (%.3f%%)' % (n, 100 * ok / n, unc, 100 * unc / n))
for s in sorted(per): print('  %-24s %6d reads  %.3f%%' % (s, per[s], 100 * conf[(s, s)] / per[s]))
for (t, p), c in sorted(conf.items(), key=lambda x: -x[1]):
    if t != p: print('   %-24s -> %-24s %d' % (t, p, c))
